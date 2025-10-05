# ======================================================
# Load required libraries
# ======================================================
library(dplyr)
library(ggplot2)
library(readr)
library(MASS)
library(tidyr)
library(lubridate)
library(broom)
library(patchwork)
library(devEMF)   # For EMF export

# ======================================================
# ---- FILE PATHS ----
# ======================================================
microbiome_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
ae_path <- "C:/Users/oofordile/Desktop/Adverse Events.csv"
pdf_path <- "C:/Users/oofordile/Desktop/Infectious_AE_Freq_Dur.pdf"
emf_path <- "C:/Users/oofordile/Desktop/Infectious_AE_Freq_Dur.emf"

# ======================================================
# ---- READ DATA ----
# ======================================================
microbiome_data <- read_csv(microbiome_path)
ae_data <- read_csv(ae_path)

# ======================================================
# ---- FILTER FOR D85 ----
# ======================================================
d85_data <- microbiome_data %>%
  filter(timepoints == "D85") %>%
  mutate(log_P = log1p(Prevotella_stercorea_7.05))

# ======================================================
# ---- NOT-ILL DATA ----
# ======================================================
not_ill_data <- d85_data %>%
  filter(Ill_Status == "Not-Ill") %>%
  dplyr::select(`Randomisation No`, log_P) %>%
  mutate(AE_NO = 0, AE_duration = 0)

# ======================================================
# ---- AE DATA (INFECTIOUS ONLY) ----
# ======================================================
ae_data <- ae_data %>%
  filter(`AE Infection` == 1) %>%
  mutate(
    StartDate = lubridate::dmy(`Start date`),
    EndDate = lubridate::dmy(`End date`)
  )

# ======================================================
# AE FREQUENCY — COUNT number of infectious episodes per child
# ======================================================
ae_freq <- ae_data %>%
  group_by(`Randomisation No`) %>%
  summarise(AE_NO = n(), .groups = "drop")

# ======================================================
# AE DURATION — SUM total infectious days per child
# ======================================================
ae_dur <- ae_data %>%
  filter(!is.na(StartDate) & !is.na(EndDate)) %>%
  mutate(duration = as.numeric(EndDate - StartDate)) %>%
  group_by(`Randomisation No`) %>%
  summarise(AE_duration = sum(duration, na.rm = TRUE), .groups = "drop")

# ======================================================
# MERGE AE SUMMARY
# ======================================================
ae_summary <- full_join(ae_freq, ae_dur, by = "Randomisation No") %>%
  replace_na(list(AE_NO = 0, AE_duration = 0))

# ======================================================
# ---- D85 ILL DATA ----
# ======================================================
ill_d85_data <- d85_data %>%
  filter(Ill_Status != "Not-Ill") %>%
  dplyr::select(`Randomisation No`, log_P) %>%
  inner_join(ae_summary, by = "Randomisation No")

# ======================================================
# ---- FINAL DATA ----
# ======================================================
final_data <- bind_rows(ill_d85_data, not_ill_data)

# ======================================================
# ---- SUMMARY DATA FOR PLOTTING ----
# ======================================================
freq_summary <- final_data %>% filter(AE_NO <= 5) %>%
  group_by(AE_NO) %>%
  summarise(mean_log_P = mean(log_P, na.rm = TRUE), .groups = "drop")

dur_summary <- final_data %>% filter(AE_duration <= 30) %>%
  group_by(AE_duration) %>%
  summarise(mean_log_P = mean(log_P, na.rm = TRUE), .groups = "drop")

# ======================================================
# ---- LINEAR MODELS ----
# ======================================================
# Frequency
freq_lm_ind <- lm(log_P ~ AE_NO, data = final_data %>% filter(AE_NO <= 5))
freq_ci_ind <- tidy(freq_lm_ind, conf.int = TRUE)
r2_ind_freq <- broom::glance(freq_lm_ind)$r.squared

freq_lm_group <- lm(mean_log_P ~ AE_NO, data = freq_summary)
freq_ci_group <- tidy(freq_lm_group, conf.int = TRUE)
r2_group_freq <- broom::glance(freq_lm_group)$r.squared

# Duration
dur_lm_ind <- lm(log_P ~ AE_duration, data = final_data %>% filter(AE_duration <= 30))
dur_ci_ind <- tidy(dur_lm_ind, conf.int = TRUE)
r2_ind_dur <- broom::glance(dur_lm_ind)$r.squared

dur_lm_group <- lm(mean_log_P ~ AE_duration, data = dur_summary)
dur_ci_group <- tidy(dur_lm_group, conf.int = TRUE)
r2_group_dur <- broom::glance(dur_lm_group)$r.squared

mid_x <- max(dur_summary$AE_duration)/2
annot_y_upper <- max(dur_summary$mean_log_P) + 1.0  # Increased from 0.60
annot_y_lower <- max(dur_summary$mean_log_P) + 0.55  # Adjusted to maintain spacing

# ======================================================
# ---- PLOTS ----
# ======================================================
# Frequency plot
p_freq <- ggplot() +
  geom_point(data=freq_summary, aes(AE_NO, mean_log_P), size=3, color="black") +
  geom_smooth(data=freq_summary, aes(AE_NO, mean_log_P), method="lm", se=TRUE, color="black", size=1.2) +
  geom_smooth(data=final_data %>% filter(AE_NO<=5), aes(AE_NO, log_P),
              method="lm", se=TRUE, color="red", size=0.8, linetype="dashed") +
  scale_x_continuous(breaks=0:5, limits=c(0,5)) +
  labs(x="Number of Infectious AEs",
       y=expression("log("*italic("P. stercorea")*") at D85"),
       title="Infectious AE Frequency") +
  annotate("text", x=2.5, y=annot_y_upper,
           label=paste0("Group: p=", signif(freq_ci_group$p.value[2],3),
                        ", slope=", round(freq_ci_group$estimate[2],3),
                        " (95% CI ", round(freq_ci_group$conf.low[2],3),"–",round(freq_ci_group$conf.high[2],3),")",
                        ", R²=", round(r2_group_freq,3)), color="black", hjust=0.5) +
  annotate("text", x=2.5, y=annot_y_lower,
           label=paste0("Individual: p=", signif(freq_ci_ind$p.value[2],3),
                        ", slope=", round(freq_ci_ind$estimate[2],3),
                        " (95% CI ", round(freq_ci_ind$conf.low[2],3),"–",round(freq_ci_ind$conf.high[2],3),")",
                        ", R²=", round(r2_ind_freq,3)), color="red", hjust=0.5) +
  theme_minimal(base_size=11) +
  theme(
    panel.grid=element_blank(),
    panel.background=element_rect(fill="white", color=NA),
    plot.background=element_rect(fill="white", color=NA),
    # Add thin axis lines with inward-pointing ticks
    axis.line = element_line(size=0.3, colour="black"),
    axis.ticks = element_line(size=0.3, colour="black"),
    axis.ticks.length = unit(-2, "mm")  # Negative value makes ticks point inward
  ) +
  expand_limits(y = annot_y_upper + 0.2)

# Duration plot
p_dur <- ggplot() +
  geom_point(data=dur_summary, aes(AE_duration, mean_log_P), size=3, color="black") +
  geom_smooth(data=dur_summary, aes(AE_duration, mean_log_P), method="lm", se=TRUE, color="black", size=1.2) +
  geom_smooth(data=final_data %>% filter(AE_duration<=30), aes(AE_duration, log_P),
              method="lm", se=TRUE, color="red", size=0.8, linetype="dashed") +
  labs(x="Total Duration of Infectious AEs (days)",
       y=expression("log("*italic("P. stercorea")*") at D85"),
       title="Infectious AE Duration") +
  annotate("text", x=mid_x, y=annot_y_upper,
           label=paste0("Group: p=", signif(dur_ci_group$p.value[2],3),
                        ", slope=", round(dur_ci_group$estimate[2],3),
                        " (95% CI ", round(dur_ci_group$conf.low[2],3),"–",round(dur_ci_group$conf.high[2],3),")",
                        ", R²=", round(r2_group_dur,3)), color="black", hjust=0.5) +
  annotate("text", x=mid_x, y=annot_y_lower,
           label=paste0("Individual: p=", signif(dur_ci_ind$p.value[2],3),
                        ", slope=", round(dur_ci_ind$estimate[2],3),
                        " (95% CI ", round(dur_ci_ind$conf.low[2],3),"–",round(dur_ci_ind$conf.high[2],3),")",
                        ", R²=", round(r2_ind_dur,3)), color="red", hjust=0.5) +
  theme_minimal(base_size=11) +
  theme(
    panel.grid=element_blank(),
    panel.background=element_rect(fill="white", color=NA),
    plot.background=element_rect(fill="white", color=NA),
    # Add thin axis lines with inward-pointing ticks
    axis.line = element_line(size=0.3, colour="black"),
    axis.ticks = element_line(size=0.3, colour="black"),
    axis.ticks.length = unit(-2, "mm")  # Negative value makes ticks point inward
  ) +
  # Expand y-axis to accommodate moved annotations
  expand_limits(y = annot_y_upper + 0.2)

# ======================================================
# ---- COMBINE PLOTS ----
# ======================================================
combined_plot <- p_freq + p_dur + plot_layout(ncol=2)

# ======================================================
# ---- EXPORT PDF ----
# ======================================================
pdf(pdf_path, width=14, height=6)
print(combined_plot)
dev.off()

# ======================================================
# ---- EXPORT EMF ----
# ======================================================
emf(emf_path, width=14, height=6)
print(combined_plot)
dev.off()

cat("Plots exported successfully:\n")
cat("- PDF:", pdf_path, "\n")
cat("- EMF:", emf_path, "\n")
