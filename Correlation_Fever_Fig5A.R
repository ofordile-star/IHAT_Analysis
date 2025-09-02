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
library(devEMF)  # For EMF export

# ======================================================
# ---- FILE PATHS ----
# ======================================================
microbiome_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
ae_path <- "C:/Users/oofordile/Desktop/Adverse Events.csv"

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
# ---- AE DATA (FEVER ONLY) ----
# ======================================================
ae_data <- ae_data %>%
  mutate(
    StartDate = dmy(`Start date`),
    EndDate = dmy(`End date`)
  ) %>%
  filter(`AE Fever` == 1)

# AE FREQUENCY
ae_freq <- ae_data %>%
  group_by(`Randomisation No`) %>%
  summarise(AE_NO = sum(`AE Fever`, na.rm = TRUE)) %>%
  ungroup()

# AE DURATION
ae_dur <- ae_data %>%
  filter(!is.na(StartDate) & !is.na(EndDate)) %>%
  mutate(duration = as.numeric(EndDate - StartDate)) %>%
  group_by(`Randomisation No`) %>%
  summarise(AE_duration = sum(duration, na.rm = TRUE)) %>%
  ungroup()

# Merge AE summaries
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
# ---- AE FREQUENCY PLOT ----
# ======================================================
freq_summary <- final_data %>%
  filter(AE_NO <= 5) %>%
  group_by(AE_NO) %>%
  summarise(mean_log_P = mean(log_P, na.rm = TRUE), .groups = "drop")

freq_lm_ind <- lm(log_P ~ AE_NO, data = final_data %>% filter(AE_NO <= 5))
freq_ci_ind <- tidy(freq_lm_ind, conf.int = TRUE)
r2_ind_freq <- glance(freq_lm_ind)$r.squared

freq_lm_group <- lm(mean_log_P ~ AE_NO, data = freq_summary)
freq_ci_group <- tidy(freq_lm_group, conf.int = TRUE)
r2_group_freq <- glance(freq_lm_group)$r.squared

# Annotation positions - red moved down by 0.5 unit
freq_ind_y   <- max(freq_summary$mean_log_P) + 4      # black annotation
freq_group_y <- max(freq_summary$mean_log_P) + 3.5    # red annotation moved down by 0.5

p_freq <- ggplot() +
  annotate("rect", xmin = -0.5, xmax = 0.5,
           ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "grey70") +
  geom_point(data = freq_summary, aes(x = AE_NO, y = mean_log_P),
             size = 3, color = "black") +
  geom_smooth(data = freq_summary, aes(x = AE_NO, y = mean_log_P),
              method = "lm", se = TRUE, color = "black", size = 1.2) +
  geom_smooth(data = final_data %>% filter(AE_NO <= 5),
              aes(x = AE_NO, y = log_P),
              method = "lm", se = TRUE,
              color = "red", size = 0.8, linetype = "dashed") +
  scale_x_continuous(breaks = 0:5, limits = c(0, 5)) +
  labs(
    x = "Number of Fever AEs",
    y = expression("log("*italic("P. stercorea")*") at D85"),
    title = "Fever AE Frequency"
  ) +
  annotate("text", x = 2.5, y = freq_group_y,
           label = paste0(
             "Group: p=", signif(freq_ci_group$p.value[2], 3),
             ", slope=", round(freq_ci_group$estimate[2], 3),
             " (95% CI ", round(freq_ci_group$conf.low[2], 3), "–", round(freq_ci_group$conf.high[2], 3), ")",
             ", R²=", round(r2_group_freq, 3)
           ), color = "red") +
  annotate("text", x = 2.5, y = freq_ind_y,
           label = paste0(
             "Individual: p=", signif(freq_ci_ind$p.value[2], 3),
             ", slope=", round(freq_ci_ind$estimate[2], 3),
             " (95% CI ", round(freq_ci_ind$conf.low[2], 3), "–", round(freq_ci_ind$conf.high[2], 3), ")",
             ", R²=", round(r2_ind_freq, 3)
           ), color = "black") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black", size = 0.3),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(-0.1, "cm")  # inward ticks
  )

# ======================================================
# ---- AE DURATION PLOT ----
# ======================================================
dur_summary <- final_data %>%
  filter(AE_duration <= 30) %>%
  group_by(AE_duration) %>%
  summarise(mean_log_P = mean(log_P, na.rm = TRUE), .groups = "drop")

dur_lm_ind <- lm(log_P ~ AE_duration, data = final_data %>% filter(AE_duration <= 30))
dur_ci_ind <- tidy(dur_lm_ind, conf.int = TRUE)
r2_ind_dur <- glance(dur_lm_ind)$r.squared

dur_lm_group <- lm(mean_log_P ~ AE_duration, data = dur_summary)
dur_ci_group <- tidy(dur_lm_group, conf.int = TRUE)
r2_group_dur <- glance(dur_lm_group)$r.squared

# Annotation positions - unchanged
dur_group_y <- max(dur_summary$mean_log_P) + 4
dur_ind_y   <- max(dur_summary$mean_log_P) + 3

p_dur <- ggplot() +
  annotate("rect", xmin = -0.5, xmax = 0.5,
           ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "grey70") +
  geom_point(data = dur_summary, aes(x = AE_duration, y = mean_log_P),
             size = 3, color = "black") +
  geom_smooth(data = dur_summary, aes(x = AE_duration, y = mean_log_P),
              method = "lm", se = TRUE, color = "black", size = 1.2) +
  geom_smooth(data = final_data %>% filter(AE_duration <= 30),
              aes(x = AE_duration, y = log_P),
              method = "lm", se = TRUE,
              color = "red", size = 0.8, linetype = "dashed") +
  labs(
    x = "Total Duration of Fever AEs (days)",
    y = expression("log("*italic("P. stercorea")*") at D85"),
    title = "Fever AE Duration"
  ) +
  annotate("text", x = 12, y = dur_group_y,
           label = paste0(
             "Group: p=", signif(dur_ci_group$p.value[2], 3),
             ", slope=", round(dur_ci_group$estimate[2], 3),
             " (95% CI ", round(dur_ci_group$conf.low[2], 3), "–", round(dur_ci_group$conf.high[2], 3), ")",
             ", R²=", round(r2_group_dur, 3)
           ), color = "black") +
  annotate("text", x = 12, y = dur_ind_y,
           label = paste0(
             "Individual: p=", signif(dur_ci_ind$p.value[2], 3),
             ", slope=", round(dur_ci_ind$estimate[2], 3),
             " (95% CI ", round(dur_ci_ind$conf.low[2], 3), "–", round(dur_ci_ind$conf.high[2], 3), ")",
             ", R²=", round(r2_ind_dur, 3)
           ), color = "red") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black", size = 0.3),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(-0.1, "cm")  # inward ticks
  )

# ======================================================
# ---- COMBINE & EXPORT ----
# ======================================================
combined_plot <- p_freq + p_dur + plot_layout(ncol = 2)

# PDF export
pdf("C:/Users/oofordile/Desktop/Fever_AE_Freq_Dur_BlackRed_AbovePoints_swapped_moved_half.pdf", width = 14, height = 6)
print(combined_plot)
dev.off()

# EMF export (Windows only)
emf("C:/Users/oofordile/Desktop/Fever_AE_Freq_Dur_BlackRed_AbovePoints_swapped_moved_half.emf", width = 14, height = 6)
print(combined_plot)
dev.off()
