# ======================================================
# Load libraries
# ======================================================
library(tidyverse)
library(ggpubr)
library(patchwork)
library(broom)
library(effectsize)
if(!requireNamespace("devEMF", quietly = TRUE)) install.packages("devEMF")
if(!requireNamespace("effectsize", quietly = TRUE)) install.packages("effectsize")

# ======================================================
# File Paths
# ======================================================
file_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
csv_path <- "C:/Users/oofordile/Desktop/ANOVA_Nature_Combined.csv"
pdf_path <- "C:/Users/oofordile/Desktop/Ill_vs_NotIll_Taxa_Boxplots.pdf"
emf_path <- "C:/Users/oofordile/Desktop/Ill_vs_NotIll_Taxa_Boxplots.emf"

# ======================================================
# Read and clean data
# ======================================================
data <- read_csv(file_path) %>%
  filter(!is.na(Ill_Status) & Ill_Status != "")

# Prepare analysis data
analysis_data <- data %>%
  dplyr::select(`Randomisation No`, timepoints, Ill_Status, `3agegroups based on age at sampling`, 13:22) %>%
  dplyr::rename(age_group = `3agegroups based on age at sampling`) %>%
  dplyr::mutate(
    timepoints = factor(timepoints, levels = c("D1", "D15", "D85")),
    Ill_Status = factor(Ill_Status, levels = c("Not-Ill", "Ill")),
    age_group = factor(age_group, levels = c("7to12 mths", "1to2years", "plus2years"))
  )

# Normalize taxa to relative abundance
taxa_names <- colnames(analysis_data)[5:14]

analysis_data_rel <- analysis_data %>%
  rowwise() %>%
  mutate(across(all_of(taxa_names), ~ .x / sum(c_across(all_of(taxa_names)), na.rm = TRUE))) %>%
  ungroup()

# ======================================================
# ANOVA and regression loop
# ======================================================
anova_list <- list()
reg_list <- list()

for(taxon in taxa_names) {
  formula <- as.formula(paste(taxon, "~ Ill_Status * timepoints + age_group"))
  model <- lm(formula, data = analysis_data_rel)
  
  anova_res <- anova(model) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    select(term, Df, `Sum Sq`, `F value`, `Pr(>F)`) %>%
    rename(df = Df, SS = `Sum Sq`, F_value = `F value`, p_value = `Pr(>F)`) %>%
    mutate(taxon = taxon)
  
  anova_list[[taxon]] <- anova_res
  
  # Age-adjusted regression for effect sizes
  lm_group <- lm(as.numeric(analysis_data_rel[[taxon]]) ~ Ill_Status + age_group, data = analysis_data_rel)
  lm_summary <- tidy(lm_group, conf.int = TRUE)
  
  ss_total <- sum(anova_res$SS, na.rm = TRUE)
  ss_ill <- anova_res$SS[anova_res$term == "Ill_Status"]
  ss_resid <- anova_res$SS[anova_res$term == "Residuals"]
  eta2_partial <- ss_ill / (ss_ill + ss_resid)
  
  f_stat <- anova_res$F_value[anova_res$term == "Ill_Status"]
  df1 <- anova_res$df[anova_res$term == "Ill_Status"]
  df2 <- anova_res$df[anova_res$term == "Residuals"]
  
  tryCatch({
    eta2_ci <- effectsize::F_to_eta2(f_stat, df1, df2, ci = 0.95)
    eta2_low <- eta2_ci$CI_low
    eta2_high <- eta2_ci$CI_high
  }, error = function(e) {
    eta2_low <- max(0, eta2_partial - 0.05)
    eta2_high <- min(1, eta2_partial + 0.05)
  })
  
  slope <- lm_summary$estimate[2]
  ci_low_slope <- lm_summary$conf.low[2]
  ci_high_slope <- lm_summary$conf.high[2]
  r2_group <- glance(lm_group)$r.squared
  
  reg_stats <- data.frame(
    Taxon = taxon,
    Effect = "Ill_Status (age-adjusted)",
    Slope = slope,
    Slope_CI_low = ci_low_slope,
    Slope_CI_high = ci_high_slope,
    R2 = r2_group,
    Eta2_partial = eta2_partial,
    Eta2_CI_low = eta2_low,
    Eta2_CI_high = eta2_high,
    stringsAsFactors = FALSE
  )
  
  reg_list[[taxon]] <- reg_stats
}

# ======================================================
# Combine ANOVA results
# ======================================================
anova_results <- bind_rows(anova_list)
terms_of_interest <- c("Ill_Status", "timepoints", "Ill_Status:timepoints", "age_group")

anova_results_filtered <- anova_results %>%
  filter(term %in% terms_of_interest) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

anova_combined <- anova_results_filtered %>%
  pivot_wider(
    id_cols = taxon,
    names_from = term,
    values_from = c(F_value, df, p_value, p_adj),
    names_sep = "_"
  )

combined_results <- anova_combined %>%
  arrange(match(taxon, taxa_names)) %>%
  transmute(
    Taxon = taxon,
    `Time Effect: F(df), adj. p` = paste0(
      round(F_value_timepoints, 3),
      "(", df_timepoints, "), ",
      signif(p_adj_timepoints, 4)
    ),
    `Group Effect: F(df), adj. p` = paste0(
      round(F_value_Ill_Status, 3),
      "(", df_Ill_Status, "), ",
      signif(p_adj_Ill_Status, 4)
    ),
    `Interaction Effect: F(df), adj. p` = paste0(
      round(`F_value_Ill_Status:timepoints`, 3),
      "(", `df_Ill_Status:timepoints`, "), ",
      signif(`p_adj_Ill_Status:timepoints`, 4)
    ),
    `Age Effect: F(df), adj. p` = paste0(
      round(F_value_age_group, 3),
      "(", df_age_group, "), ",
      signif(p_adj_age_group, 4)
    )
  )

# ======================================================
# Add Nature report for significant p_time
# ======================================================
reg_combined <- bind_rows(reg_list)

# Extract Ill_Status and Timepoints adjusted stats
ill_status_adj <- anova_results_filtered %>%
  filter(term == "Ill_Status") %>%
  select(taxon, p_adj, F_value, df) %>%
  rename(Taxon = taxon, Ill_p_adj = p_adj, Ill_F = F_value, Ill_df = df)

time_adj <- anova_results_filtered %>%
  filter(term == "timepoints") %>%
  select(taxon, p_adj, F_value, df) %>%
  rename(Taxon = taxon, Time_p_adj = p_adj, Time_F = F_value, Time_df = df)

# Merge
reg_combined <- reg_combined %>%
  left_join(ill_status_adj, by = "Taxon") %>%
  left_join(time_adj, by = "Taxon")

# Nature report for Ill_Status (existing)
reg_combined <- reg_combined %>%
  mutate(
    Nature_Report = paste0(
      "Ill_Status (age-adj): F(", Ill_df, ") = ", round(Ill_F, 2),
      ", adj p = ", signif(Ill_p_adj, 3),
      ", eta2_p = ", round(Eta2_partial, 3),
      ", 95% CI [", round(Eta2_CI_low, 3), ", ", round(Eta2_CI_high, 3), "]",
      "; slope = ", round(Slope, 3),
      " (95% CI [", round(Slope_CI_low, 3), ", ", round(Slope_CI_high, 3), "])"
    )
  )

# Nature report for significant p_time (< 0.05)
reg_combined <- reg_combined %>%
  mutate(
    Nature_Report_Time = ifelse(
      Time_p_adj < 0.05,
      paste0(
        "Time Effect: F(", Time_df, ") = ", round(Time_F, 2),
        ", adj p = ", signif(Time_p_adj, 3),
        ", eta2_p = ", round(Eta2_partial, 3),
        ", 95% CI [", round(Eta2_CI_low, 3), ", ", round(Eta2_CI_high, 3), "]",
        "; slope = ", round(Slope, 3),
        " (95% CI [", round(Slope_CI_low, 3), ", ", round(Slope_CI_high, 3), "])"
      ),
      NA
    )
  )

# ======================================================
# Merge and save CSV
# ======================================================
final_csv <- combined_results %>%
  left_join(
    reg_combined %>%
      select(Taxon, Slope, Slope_CI_low, Slope_CI_high, R2, Eta2_partial, Eta2_CI_low, Eta2_CI_high, Nature_Report, Nature_Report_Time),
    by = "Taxon"
  )

write_csv(final_csv, csv_path)
cat("CSV saved to Desktop:", csv_path, "\n")

# ======================================================
# Plot Section (unchanged)
# ======================================================
plot_taxa <- c("Prevotella_stercorea_7.05", "Escherichia_coli_3.51")

pvals_for_plot <- anova_results_filtered %>%
  filter(taxon %in% plot_taxa & term == "Ill_Status") %>%
  select(taxon, p_adj) %>%
  mutate(p.label = paste0("adj p = ", signif(p_adj, 3)))

long_data <- analysis_data %>%
  rowwise() %>%
  mutate(
    Prevotella_rel = Prevotella_stercorea_7.05 / sum(c_across(all_of(taxa_names)), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(`Randomisation No`, timepoints, Ill_Status, Prevotella_rel, Escherichia_coli_3.51) %>%
  rename(Prevotella_stercorea_7.05 = Prevotella_rel) %>%
  pivot_longer(cols = all_of(plot_taxa), names_to = "Taxa", values_to = "Abundance") %>%
  mutate(
    Value = ifelse(Taxa == "Escherichia_coli_3.51", log10(Abundance + 1e-6), Abundance),
    TaxaLabel = recode(Taxa,
                       "Prevotella_stercorea_7.05" = "italic('Prevotella stercorea')",
                       "Escherichia_coli_3.51" = "italic('Escherichia coli')"),
    Ylabel = ifelse(Taxa == "Escherichia_coli_3.51", "Log Abundance", "Relative Abundance")
  )

long_data$Taxa <- factor(long_data$Taxa, levels = plot_taxa)
long_data$TaxaLabel <- factor(long_data$TaxaLabel, levels = c("italic('Prevotella stercorea')", "italic('Escherichia coli')"))

make_taxa_plot <- function(df, ylab, force_ylim = NULL) {
  y_max <- max(df$Value, na.rm = TRUE)
  y_min <- min(df$Value, na.rm = TRUE)
  y_range <- y_max - y_min
  taxa_name <- unique(df$Taxa)
  pval <- pvals_for_plot$p.label[pvals_for_plot$taxon == taxa_name]
  annot_y <- y_max + y_range * 0.05
  
  p <- ggplot(df, aes(x = timepoints, y = Value, fill = Ill_Status)) +
    geom_boxplot(position = position_dodge(0.7), width = 0.6, alpha = 0.8, color = "black", outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), alpha = 0.4, size = 1) +
    stat_summary(fun = median, geom = "line", aes(group = Ill_Status, color = Ill_Status), size = 1, position = position_dodge(0.7)) +
    geom_text(data = df %>% distinct(TaxaLabel),
              aes(x = 2, y = annot_y, label = pval), inherit.aes = FALSE, size = 5) +
    scale_fill_manual(values = c("Ill" = "#FF6666", "Not-Ill" = "#6699FF")) +
    scale_color_manual(values = c("Ill" = "red", "Not-Ill" = "blue")) +
    labs(x = "Study Day", y = ylab, fill = "Status", color = "Status") +
    facet_wrap(~TaxaLabel, labeller = label_parsed) +
    theme_pubr(border = TRUE, base_size = 14) +
    theme(
      strip.background = element_rect(fill = "gray95", color = "black", size = 0.5),
      strip.text = element_text(size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.spacing = unit(0.5, "lines"),
      axis.ticks = element_line(size = 0.3, color = "black"),
      axis.ticks.length = unit(-2, "mm")
    )
  
  if (!is.null(force_ylim)) p <- p + coord_cartesian(ylim = force_ylim)
  else p <- p + expand_limits(y = y_max + y_range*0.1)
  
  return(p)
}

prev_plot <- long_data %>% filter(Taxa == "Prevotella_stercorea_7.05") %>% make_taxa_plot("Relative Abundance")
ecoli_plot <- long_data %>% filter(Taxa == "Escherichia_coli_3.51") %>% make_taxa_plot("Log Abundance", force_ylim = c(0, 5))

taxa_plot <- prev_plot + ecoli_plot + plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom", legend.justification = "center", legend.margin = margin(t = 10, b = 10))

taxa_plot <- taxa_plot + plot_annotation(
  title = "Taxa with significant Ill vs Not-Ill gap (FDR and Age-Adjusted p < 0.05)",
  theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
)

ggsave(pdf_path, taxa_plot, width = 10, height = 6, dpi = 300)
devEMF::emf(emf_path, width = 10, height = 6)
print(taxa_plot)
dev.off()
