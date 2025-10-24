# ======================================================
# COMBINED 2x4 AE CORRELATION PLOTS + NATURE-STYLE CSV
# WITH UNIFORM VERTICAL SPACING
# ======================================================

# ======================================================
# Load libraries
# ======================================================
library(dplyr)
library(ggplot2)
library(readr)
library(lubridate)
library(broom)
library(patchwork)
library(devEMF)
library(grid)
library(tidyr)

# ======================================================
# File paths
# ======================================================
microbiome_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
ae_path <- "C:/Users/oofordile/Desktop/Adverse Events.csv"
pdf_path <- "C:/Users/oofordile/Desktop/Combined_AE_Freq_Dur_CloserSlope_Wide.pdf"
emf_path <- "C:/Users/oofordile/Desktop/Combined_AE_Freq_Dur_CloserSlope_Wide.emf"
csv_path <- "C:/Users/oofordile/Desktop/AE_Freq_Dur_NatureStatistics.csv"

# ======================================================
# Load data
# ======================================================
microbiome_data <- read_csv(microbiome_path)
ae_data <- read_csv(ae_path)

# ======================================================
# Compute total infectious diagnoses
# ======================================================
total_diarrhoeal <- sum(ae_data$`AE Diarrhoea` == 1, na.rm = TRUE)
total_ari      <- sum(ae_data$`AE ARI` == 1, na.rm = TRUE)
total_fever    <- sum(ae_data$`AE Fever` == 1, na.rm = TRUE)
total_infectious <- total_diarrhoeal + total_ari + total_fever

cat("Total AE Diarrhoea diagnoses: ", total_diarrhoeal, "\n")
cat("Total AE ARI diagnoses: ", total_ari, "\n")
cat("Total AE Fever diagnoses: ", total_fever, "\n")
cat("Total infectious diagnoses (all types combined, counting multiple per AE incident): ", 
    total_infectious, "\n")

# ======================================================
# Filter D85 and log-transform P. stercorea
# ======================================================
d85_data <- microbiome_data %>%
  filter(timepoints == "D85") %>%
  mutate(log_P = log1p(Prevotella_stercorea_7.05))

# ======================================================
# Function to prepare AE plots with uniform spacing
# ======================================================
prepare_ae_plot <- function(ae_column, ae_label, reduce_freq_spacing = FALSE){
  
  ae_filtered <- ae_data %>%
    mutate(StartDate = dmy(`Start date`),
           EndDate = dmy(`End date`)) %>%
    filter(.data[[ae_column]] == 1)
  
  ae_freq <- ae_filtered %>%
    group_by(`Randomisation No`) %>%
    summarise(AE_NO = n(), .groups="drop")
  
  ae_dur <- ae_filtered %>%
    filter(!is.na(StartDate) & !is.na(EndDate)) %>%
    mutate(duration = as.numeric(EndDate - StartDate)) %>%
    group_by(`Randomisation No`) %>%
    summarise(AE_duration = sum(duration, na.rm=TRUE), .groups="drop")
  
  ae_summary <- full_join(ae_freq, ae_dur, by="Randomisation No") %>%
    replace_na(list(AE_NO=0, AE_duration=0))
  
  ill_d85 <- d85_data %>%
    filter(Ill_Status != "Not-Ill") %>%
    dplyr::select(`Randomisation No`, log_P) %>%
    inner_join(ae_summary, by="Randomisation No")
  
  not_ill <- d85_data %>%
    filter(Ill_Status == "Not-Ill") %>%
    dplyr::select(`Randomisation No`, log_P) %>%
    mutate(AE_NO=0, AE_duration=0)
  
  final_data <- bind_rows(ill_d85, not_ill)
  
  # -------------------
  # Frequency plot with uniform spacing
  # -------------------
  freq_summary <- final_data %>% filter(AE_NO <= 5) %>%
    group_by(AE_NO) %>%
    summarise(mean_log_P = mean(log_P, na.rm=TRUE), .groups="drop")
  
  freq_lm_ind <- lm(log_P ~ AE_NO, data = final_data %>% filter(AE_NO <= 5))
  freq_ci_ind <- tidy(freq_lm_ind, conf.int = TRUE)
  r2_ind_freq <- glance(freq_lm_ind)$r.squared
  
  freq_lm_group <- lm(mean_log_P ~ AE_NO, data = freq_summary)
  freq_ci_group <- tidy(freq_lm_group, conf.int = TRUE)
  r2_group_freq <- glance(freq_lm_group)$r.squared
  
  # Uniform spacing: 1 unit within group/individual, 1.5 units between them
  # Doubly reduced for Infectious and Fever frequency
  y_max <- max(freq_summary$mean_log_P)
  if (reduce_freq_spacing) {
    freq_group_y_top <- y_max + 4.0
    freq_group_y_slope <- freq_group_y_top - 0.5      # 0.5 paragraph space within
    freq_ind_y_top <- freq_group_y_slope - 0.75       # 0.75 paragraph spaces between
    freq_ind_y_slope <- freq_ind_y_top - 0.5          # 0.5 paragraph space within
  } else {
    freq_group_y_top <- y_max + 4.0
    freq_group_y_slope <- freq_group_y_top - 1.0      # 1 paragraph space within
    freq_ind_y_top <- freq_group_y_slope - 1.5        # 1.5 paragraph spaces between
    freq_ind_y_slope <- freq_ind_y_top - 1.0          # 1 paragraph space within
  }
  
  p_freq <- ggplot() +
    geom_point(data=freq_summary, aes(AE_NO, mean_log_P), size=4, color="black") +
    geom_smooth(data=freq_summary, aes(AE_NO, mean_log_P), method="lm", se=TRUE, color="black", size=1.3) +
    geom_smooth(data=final_data %>% filter(AE_NO <= 5), aes(AE_NO, log_P), method="lm", se=TRUE, color="red", size=1, linetype="dashed") +
    scale_x_continuous(breaks=0:5, limits=c(0,5)) +
    labs(x=paste("Number of", ae_label, "AEs"),
         y=expression("log("*italic("P. stercorea")*") at D85"),
         title=paste(ae_label,"AE Frequency")) +
    annotate("text", x=2.5, y=freq_group_y_top,
             label=paste0("Group: p=", signif(freq_ci_group$p.value[2],3), ", R^2=", round(r2_group_freq,3)),
             color="black", hjust=0.5, size=5.5) +
    annotate("text", x=2.5, y=freq_group_y_slope,
             label=paste0("slope=", round(freq_ci_group$estimate[2],3),
                          " (95% CI ", round(freq_ci_group$conf.low[2],3), "-", round(freq_ci_group$conf.high[2],3), ")"),
             color="black", hjust=0.5, size=5.5) +
    annotate("text", x=2.5, y=freq_ind_y_top,
             label=paste0("Individual: p=", signif(freq_ci_ind$p.value[2],3), ", R^2=", round(r2_ind_freq,3)),
             color="red", hjust=0.5, size=5.5) +
    annotate("text", x=2.5, y=freq_ind_y_slope,
             label=paste0("slope=", round(freq_ci_ind$estimate[2],3),
                          " (95% CI ", round(freq_ci_ind$conf.low[2],3), "-", round(freq_ci_ind$conf.high[2],3), ")"),
             color="red", hjust=0.5, size=5.5) +
    theme_minimal(base_size=16) +
    theme(panel.grid=element_blank(),
          panel.background=element_rect(fill="white", color=NA),
          plot.background=element_rect(fill="white", color=NA),
          axis.line=element_line(color="black", size=0.5),
          axis.ticks=element_line(color="black", size=0.5),
          axis.ticks.length=unit(0.2,"cm"),
          plot.title=element_text(size=18, face="bold", hjust=0.5),
          axis.title=element_text(size=17),
          axis.text=element_text(size=15))
  
  # -------------------
  # Duration plot with uniform spacing
  # -------------------
  dur_summary <- final_data %>% filter(AE_duration <= 30) %>%
    group_by(AE_duration) %>%
    summarise(mean_log_P = mean(log_P, na.rm=TRUE), .groups="drop")
  
  dur_lm_ind <- lm(log_P ~ AE_duration, data = final_data %>% filter(AE_duration <= 30))
  dur_ci_ind <- tidy(dur_lm_ind, conf.int = TRUE)
  r2_ind_dur <- glance(dur_lm_ind)$r.squared
  
  dur_lm_group <- lm(mean_log_P ~ AE_duration, data = dur_summary)
  dur_ci_group <- tidy(dur_lm_group, conf.int = TRUE)
  r2_group_dur <- glance(dur_lm_group)$r.squared
  
  # Uniform spacing: 1 unit within group/individual, 1.5 units between them
  x_mid <- max(dur_summary$AE_duration)/2
  y_max_dur <- max(dur_summary$mean_log_P)
  dur_group_y_top <- y_max_dur + 4.0
  dur_group_y_slope <- dur_group_y_top - 1.0        # 1 paragraph space within
  dur_ind_y_top <- dur_group_y_slope - 1.5          # 1.5 paragraph spaces between
  dur_ind_y_slope <- dur_ind_y_top - 1.0            # 1 paragraph space within
  
  p_dur <- ggplot() +
    geom_point(data=dur_summary, aes(AE_duration, mean_log_P), size=4, color="black") +
    geom_smooth(data=dur_summary, aes(AE_duration, mean_log_P), method="lm", se=TRUE, color="black", size=1.3) +
    geom_smooth(data=final_data %>% filter(AE_duration <= 30), aes(AE_duration, log_P), method="lm", se=TRUE, color="red", size=1, linetype="dashed") +
    labs(x=paste("Total Duration of", ae_label, "AEs (days)"),
         y=expression("log("*italic("P. stercorea")*") at D85"),
         title=paste(ae_label,"AE Duration")) +
    annotate("text", x=x_mid, y=dur_group_y_top,
             label=paste0("Group: p=", signif(dur_ci_group$p.value[2],3), ", R^2=", round(r2_group_dur,3)),
             color="black", hjust=0.5, size=5.5) +
    annotate("text", x=x_mid, y=dur_group_y_slope,
             label=paste0("slope=", round(dur_ci_group$estimate[2],3),
                          " (95% CI ", round(dur_ci_group$conf.low[2],3), "-", round(dur_ci_group$conf.high[2],3), ")"),
             color="black", hjust=0.5, size=5.5) +
    annotate("text", x=x_mid, y=dur_ind_y_top,
             label=paste0("Individual: p=", signif(dur_ci_ind$p.value[2],3), ", R^2=", round(r2_ind_dur,3)),
             color="red", hjust=0.5, size=5.5) +
    annotate("text", x=x_mid, y=dur_ind_y_slope,
             label=paste0("slope=", round(dur_ci_ind$estimate[2],3),
                          " (95% CI ", round(dur_ci_ind$conf.low[2],3), "-", round(dur_ci_ind$conf.high[2],3), ")"),
             color="red", hjust=0.5, size=5.5) +
    theme_minimal(base_size=16) +
    theme(panel.grid=element_blank(),
          panel.background=element_rect(fill="white", color=NA),
          plot.background=element_rect(fill="white", color=NA),
          axis.line=element_line(color="black", size=0.5),
          axis.ticks=element_line(color="black", size=0.5),
          axis.ticks.length=unit(0.2,"cm"),
          plot.title=element_text(size=18, face="bold", hjust=0.5),
          axis.title=element_text(size=17),
          axis.text=element_text(size=15))
  
  return(list(freq=p_freq, dur=p_dur))
}

# ======================================================
# Generate plots for all AE types
# ======================================================
plots_infect <- prepare_ae_plot("AE Infection", "Infectious", reduce_freq_spacing = TRUE)
plots_diarr <- prepare_ae_plot("AE Diarrhoea", "Diarrhoea")
plots_ari    <- prepare_ae_plot("AE ARI", "ARI")
plots_fever  <- prepare_ae_plot("AE Fever", "Fever", reduce_freq_spacing = TRUE)

# ======================================================
# Combine 2x4 panels
# ======================================================
row1 <- wrap_plots(plots_infect$freq, plots_diarr$freq, plots_ari$freq, plots_fever$freq, ncol=4)
row2 <- wrap_plots(plots_infect$dur, plots_diarr$dur, plots_ari$dur, plots_fever$dur, ncol=4)
combined_plot <- row1 / row2 + plot_layout(heights=c(1,1), guides="collect") & theme(plot.margin=margin(5,5,5,5))

# ======================================================
# Export plots
# ======================================================
pdf(pdf_path, width=22, height=12)
print(combined_plot)
dev.off()

emf(emf_path, width=22, height=12)
print(combined_plot)
dev.off()

cat("Combined 2x4 AE plots exported successfully:\n")
cat("- PDF:", pdf_path, "\n")
cat("- EMF:", emf_path, "\n")

# ======================================================
# Compute Nature-style CSV
# ======================================================
compute_ae_stats <- function(ae_column, ae_label){
  
  ae_filtered <- ae_data %>%
    mutate(StartDate = dmy(`Start date`),
           EndDate = dmy(`End date`)) %>%
    filter(.data[[ae_column]] == 1)
  
  ae_freq <- ae_filtered %>%
    group_by(`Randomisation No`) %>%
    summarise(AE_NO = n(), .groups="drop")
  
  ae_dur <- ae_filtered %>%
    filter(!is.na(StartDate) & !is.na(EndDate)) %>%
    mutate(duration = as.numeric(EndDate - StartDate)) %>%
    group_by(`Randomisation No`) %>%
    summarise(AE_duration = sum(duration, na.rm=TRUE), .groups="drop")
  
  ae_summary <- full_join(ae_freq, ae_dur, by="Randomisation No") %>%
    replace_na(list(AE_NO=0, AE_duration=0))
  
  ill_d85 <- d85_data %>%
    filter(Ill_Status != "Not-Ill") %>%
    dplyr::select(`Randomisation No`, log_P) %>%
    inner_join(ae_summary, by="Randomisation No")
  
  not_ill <- d85_data %>%
    filter(Ill_Status == "Not-Ill") %>%
    dplyr::select(`Randomisation No`, log_P) %>%
    mutate(AE_NO=0, AE_duration=0)
  
  final_data <- bind_rows(ill_d85, not_ill)
  
  # Frequency regressions
  freq_summary <- final_data %>% filter(AE_NO <= 5) %>%
    group_by(AE_NO) %>%
    summarise(mean_log_P = mean(log_P, na.rm=TRUE), .groups="drop")
  
  freq_lm_ind <- lm(log_P ~ AE_NO, data = final_data %>% filter(AE_NO <= 5))
  freq_ci_ind <- tidy(freq_lm_ind, conf.int = TRUE)
  r2_ind_freq <- glance(freq_lm_ind)$r.squared
  
  freq_lm_group <- lm(mean_log_P ~ AE_NO, data = freq_summary)
  freq_ci_group <- tidy(freq_lm_group, conf.int = TRUE)
  r2_group_freq <- glance(freq_lm_group)$r.squared
  
  freq_stats <- data.frame(
    AE_Type = ae_label,
    Metric = "Frequency",
    Level = c("Group","Individual"),
    Slope = c(freq_ci_group$estimate[2], freq_ci_ind$estimate[2]),
    CI_low = c(freq_ci_group$conf.low[2], freq_ci_ind$conf.low[2]),
    CI_high = c(freq_ci_group$conf.high[2], freq_ci_ind$conf.high[2]),
    p_value = c(freq_ci_group$p.value[2], freq_ci_ind$p.value[2]),
    R2 = c(r2_group_freq, r2_ind_freq)
  ) %>%
    mutate(Nature_Report = paste0(Level, ": slope=", round(Slope,3),
                                  " (95% CI ", round(CI_low,3), "-", round(CI_high,3), ")",
                                  ", p=", signif(p_value,3),
                                  ", R^2=", round(R2,3)))
  
  # Duration regressions
  dur_summary <- final_data %>% filter(AE_duration <= 30) %>%
    group_by(AE_duration) %>%
    summarise(mean_log_P = mean(log_P, na.rm=TRUE), .groups="drop")
  
  dur_lm_ind <- lm(log_P ~ AE_duration, data = final_data %>% filter(AE_duration <= 30))
  dur_ci_ind <- tidy(dur_lm_ind, conf.int = TRUE)
  r2_ind_dur <- glance(dur_lm_ind)$r.squared
  
  dur_lm_group <- lm(mean_log_P ~ AE_duration, data = dur_summary)
  dur_ci_group <- tidy(dur_lm_group, conf.int = TRUE)
  r2_group_dur <- glance(dur_lm_group)$r.squared
  
  dur_stats <- data.frame(
    AE_Type = ae_label,
    Metric = "Duration",
    Level = c("Group","Individual"),
    Slope = c(dur_ci_group$estimate[2], dur_ci_ind$estimate[2]),
    CI_low = c(dur_ci_group$conf.low[2], dur_ci_ind$conf.low[2]),
    CI_high = c(dur_ci_group$conf.high[2], dur_ci_ind$conf.high[2]),
    p_value = c(dur_ci_group$p.value[2], dur_ci_ind$p.value[2]),
    R2 = c(r2_group_dur, r2_ind_dur)
  ) %>%
    mutate(Nature_Report = paste0(Level, ": slope=", round(Slope,3),
                                  " (95% CI ", round(CI_low,3), "-", round(CI_high,3), ")",
                                  ", p=", signif(p_value,3),
                                  ", R^2=", round(R2,3)))
  
  return(bind_rows(freq_stats, dur_stats))
}

# ======================================================
# Compute statistics for all AE types
# ======================================================
stats_infect <- compute_ae_stats("AE Infection", "Infectious")
stats_diarr  <- compute_ae_stats("AE Diarrhoea", "Diarrhoea")
stats_ari    <- compute_ae_stats("AE ARI", "ARI")
stats_fever  <- compute_ae_stats("AE Fever", "Fever")

all_stats <- bind_rows(stats_infect, stats_diarr, stats_ari, stats_fever)
write_csv(all_stats, csv_path)

cat("Nature-style AE statistics exported successfully:\n")
cat("- CSV:", csv_path, "\n")