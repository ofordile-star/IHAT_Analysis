# --- Load Libraries ---
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(vcd)
library(gridExtra)
library(forcats)
library(cowplot)
library(knitr)

# --- Load Data ---
data_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
df <- read_csv(data_path)

# --- Define correct age group levels ---
age_levels <- c("7to12 mths", "1to2years", "plus2years")

# --- Assign factor with exact levels ---
df$age_group <- factor(df$`3agegroups based on age at sampling`,
                       levels = age_levels)

# --- Define fill colors ---
age_colors <- c(
  "7to12 mths" = "#A6D854",   # light green
  "1to2years"  = "#FDB863",   # light orange
  "plus2years" = "#80B1D3"    # light blue
)

# ================================
# COMPREHENSIVE DEMOGRAPHIC ANALYSIS
# ================================

cat("========================================\n")
cat("DEMOGRAPHIC ANALYSIS FOR REVIEWER\n")
cat("========================================\n\n")

# Check available columns first
cat("Available columns in dataset:\n")
print(names(df))
cat("\n")

# Look for gender-related columns
gender_cols <- names(df)[grepl("gender|sex|female|male", names(df), ignore.case = TRUE)]
cat("Potential gender columns found:", paste(gender_cols, collapse = ", "), "\n\n")

# Filter for complete demographic data
df_demo <- df %>%
  filter(!is.na(Ill_Status), !is.na(age_group), 
         !is.na(`age at sampling`), !is.na(gender)) %>%
  mutate(
    Ill_Status = factor(Ill_Status, levels = c("Not-Ill", "Ill")),
    Gender = factor(gender, levels = c("Female", "Male"))  # Using lowercase 'gender' column
  )

# 1. AGE DISTRIBUTION ANALYSIS
cat("1. AGE DISTRIBUTION BY ILLNESS STATUS\n")
cat("=====================================\n")

# Summary statistics for age by illness status
age_summary <- df_demo %>%
  group_by(Ill_Status) %>%
  summarise(
    N = n(),
    Mean_Age = round(mean(`age at sampling`, na.rm = TRUE), 2),
    Median_Age = round(median(`age at sampling`, na.rm = TRUE), 2),
    SD_Age = round(sd(`age at sampling`, na.rm = TRUE), 2),
    Min_Age = min(`age at sampling`, na.rm = TRUE),
    Max_Age = max(`age at sampling`, na.rm = TRUE),
    .groups = "drop"
  )

print(kable(age_summary, caption = "Age Summary by Illness Status"))

# Statistical test for age differences
age_test <- wilcox.test(`age at sampling` ~ Ill_Status, data = df_demo)
cat("\nMann-Whitney U test for age differences between groups:\n")
cat(sprintf("W = %.2f, p-value = %.4f\n", age_test$statistic, age_test$p.value))

if(age_test$p.value < 0.05) {
  cat("*** Significant age difference between Ill and Not-Ill groups ***\n")
} else {
  cat("No significant age difference between groups (potential confounder ruled out)\n")
}

# 2. AGE GROUP DISTRIBUTION ANALYSIS
cat("\n2. AGE GROUP DISTRIBUTION BY ILLNESS STATUS\n")
cat("==========================================\n")

# Cross-tabulation
age_group_table <- table(df_demo$age_group, df_demo$Ill_Status)
print(age_group_table)

# Proportions within each illness status
prop_table <- prop.table(age_group_table, margin = 2) * 100
cat("\nProportions within each illness status (%):\n")
print(round(prop_table, 1))

# Chi-square test for age group distribution
chisq_age_group <- chisq.test(age_group_table)
cat(sprintf("\nChi-square test: X²(%d) = %.3f, p-value = %.4f\n", 
            chisq_age_group$parameter, chisq_age_group$statistic, chisq_age_group$p.value))

# 3. GENDER DISTRIBUTION ANALYSIS
cat("\n3. GENDER DISTRIBUTION BY ILLNESS STATUS\n")
cat("========================================\n")

# Gender cross-tabulation
gender_table <- table(df_demo$Gender, df_demo$Ill_Status)
print(gender_table)

# Gender proportions
gender_prop <- prop.table(gender_table, margin = 2) * 100
cat("\nGender proportions within each illness status (%):\n")
print(round(gender_prop, 1))

# Chi-square test for gender distribution
chisq_gender <- chisq.test(gender_table)
cat(sprintf("\nChi-square test for gender: X²(%d) = %.3f, p-value = %.4f\n", 
            chisq_gender$parameter, chisq_gender$statistic, chisq_gender$p.value))

# 4. COMPREHENSIVE DEMOGRAPHIC TABLE
cat("\n4. COMPREHENSIVE DEMOGRAPHIC COMPARISON TABLE\n")
cat("============================================\n")

demo_comparison <- df_demo %>%
  group_by(Ill_Status) %>%
  summarise(
    N_total = n(),
    N_unique_subjects = n_distinct(`Randomisation No`),
    Female_n = sum(Gender == "Female", na.rm = TRUE),
    Female_percent = round(100 * Female_n / N_total, 1),
    Age_7to12_n = sum(age_group == "7to12 mths", na.rm = TRUE),
    Age_7to12_percent = round(100 * Age_7to12_n / N_total, 1),
    Age_1to2_n = sum(age_group == "1to2years", na.rm = TRUE),
    Age_1to2_percent = round(100 * Age_1to2_n / N_total, 1),
    Age_plus2_n = sum(age_group == "plus2years", na.rm = TRUE),
    Age_plus2_percent = round(100 * Age_plus2_n / N_total, 1),
    Mean_Age_months = round(mean(`age at sampling`, na.rm = TRUE), 1),
    SD_Age_months = round(sd(`age at sampling`, na.rm = TRUE), 1),
    .groups = "drop"
  )

print(demo_comparison)

# 5. STATISTICAL SUMMARY FOR REVIEWER
cat("\n5. SUMMARY FOR REVIEWER RESPONSE\n")
cat("================================\n")

cat("Key findings addressing reviewer concerns:\n\n")

cat(sprintf("• Total sample: %d samples from %d unique subjects\n", 
            nrow(df_demo), n_distinct(df_demo$`Randomisation No`)))

cat(sprintf("• Age range: %.1f to %.1f months (mean ± SD: %.1f ± %.1f months)\n",
            min(df_demo$`age at sampling`, na.rm = TRUE),
            max(df_demo$`age at sampling`, na.rm = TRUE),
            mean(df_demo$`age at sampling`, na.rm = TRUE),
            sd(df_demo$`age at sampling`, na.rm = TRUE)))

if(age_test$p.value >= 0.05) {
  cat("• No significant age differences between Ill and Not-Ill groups (p = ", 
      sprintf("%.3f", age_test$p.value), "), ruling out age as a major confounder\n")
} else {
  cat("• WARNING: Significant age differences detected between groups (p = ", 
      sprintf("%.3f", age_test$p.value), "), age may be a confounder\n")
}

if(chisq_age_group$p.value >= 0.05) {
  cat("• Age group distributions are similar between Ill and Not-Ill groups (p = ", 
      sprintf("%.3f", chisq_age_group$p.value), ")\n")
} else {
  cat("• Age group distributions differ between Ill and Not-Ill groups (p = ", 
      sprintf("%.3f", chisq_age_group$p.value), ")\n")
}

if(chisq_gender$p.value >= 0.05) {
  cat("• Gender distributions are balanced between groups (p = ", 
      sprintf("%.3f", chisq_gender$p.value), ")\n")
} else {
  cat("• Gender distributions differ between groups (p = ", 
      sprintf("%.3f", chisq_gender$p.value), ")\n")
}

# ================================
# VISUALIZATION PLOTS
# ================================

# PLOT 1: Age Distribution by Illness Status (NEW)
plot_age_dist <- ggplot(df_demo, aes(x = `age at sampling`, fill = Ill_Status)) +
  geom_histogram(alpha = 0.7, bins = 15, position = "identity") +
  facet_wrap(~Ill_Status, ncol = 1) +
  scale_fill_manual(values = c("Not-Ill" = "#4575B4", "Ill" = "#D73027")) +
  labs(x = "Age at Sampling (months)", 
       y = "Number of Samples",
       title = "Age Distribution by Illness Status") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
           label = sprintf("Mann-Whitney U: p = %.3f", age_test$p.value),
           size = 4) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        panel.grid = element_blank())

# PLOT 2: Prevotella stercorea Boxplot (EXISTING)
df_d <- df %>%
  filter(!is.na(`Prevotella_stercorea_7.05`)) %>%
  mutate(logP = log(`Prevotella_stercorea_7.05` + 1))

df_d$age_group <- factor(df_d$age_group, levels = age_levels)

anova_res <- aov(logP ~ age_group, data = df_d)
anova_summary <- summary(anova_res)[[1]]
anova_F <- anova_summary[["F value"]][1]

anova_label <- paste0("ANOVA F(2) = ", round(anova_F, 2), ", P < 2e-16")

plot_d_inner <- ggplot(df_d, aes(x = age_group, y = logP, fill = age_group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 23, fill = "white", size = 2) +
  scale_fill_manual(values = age_colors, drop = FALSE) +
  labs(
    x = "Age Group",
    y = expression("Log Abundance of " * italic("Prevotella stercorea"))
  ) +
  annotate("text", x = 2, y = max(df_d$logP, na.rm = TRUE) * 1.05,
           label = anova_label, size = 5, hjust = 0.5, color = "black") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )

plot_d <- ggdraw() +
  draw_label(expression(italic("Prevotella stercorea") * " Abundance by Age Group"),
             x = 0.5, y = 0.98, hjust = 0.5, vjust = 1,
             size = 18, fontface = "bold", color = "black") +
  draw_plot(plot_d_inner, y = 0, height = 0.95)

# PLOT 3: Illness Bar Plot (EXISTING, Enhanced)
df_e <- df %>%
  filter(!is.na(Ill_Status), !is.na(age_group))

df_e$Ill_Status <- factor(df_e$Ill_Status, levels = c("Not-Ill", "Ill"))
df_e$age_group <- factor(df_e$age_group, levels = age_levels)

ill_df <- df_e %>%
  group_by(age_group, Ill_Status) %>%
  summarise(unique_ids = n_distinct(`Randomisation No`), .groups = "drop") %>%
  rename(count = unique_ids)

ill_tab <- table(df_e$age_group, df_e$Ill_Status)
if (all(dim(ill_tab) == c(length(age_levels), 2)) && all(rowSums(ill_tab) > 0)) {
  chisq_res <- chisq.test(ill_tab)
  test_label <- paste0("Chi-sq(", chisq_res$parameter, ") = ", 
                       round(chisq_res$statistic, 2), 
                       ", P = ", formatC(chisq_res$p.value, format = "e", digits = 2))
} else {
  fisher_res <- fisher.test(ill_tab)
  test_label <- paste0("Fisher's Exact test P = ", formatC(fisher_res$p.value, format = "e", digits = 2))
}

plot_e_inner <- ggplot(ill_df, aes(x = age_group, y = count, fill = age_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  facet_wrap(~Ill_Status) +
  scale_fill_manual(values = age_colors, drop = FALSE) +
  labs(
    x = "Age Group",
    y = "Number of Children"
  ) +
  coord_cartesian(ylim = c(0, max(ill_df$count) * 1.15), clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.margin = margin(t = 30, r = 10, b = 10, l = 10),
    panel.grid = element_blank()
  )

plot_e <- ggdraw() +
  draw_label("Ill Status Distribution by Age Group",
             x = 0.5, y = 0.98, hjust = 0.5, vjust = 1,
             size = 18, fontface = "bold", color = "black") +
  draw_label(test_label,
             x = 0.5, y = 0.92, hjust = 0.5, vjust = 1,
             size = 14, fontface = "plain", color = "black") +
  draw_plot(plot_e_inner, y = 0, height = 0.9)

# ================================
# Export all plots
# ================================

# 3-panel plot including new age distribution
pdf("C:/Users/oofordile/Desktop/Comprehensive_Demographics_Analysis.pdf", width = 18, height = 6)
grid.arrange(plot_age_dist, plot_d, plot_e, nrow = 1)
dev.off()

png("C:/Users/oofordile/Desktop/Comprehensive_Demographics_Analysis.png",
    width = 18, height = 6, units = "in", res = 300)
grid.arrange(plot_age_dist, plot_d, plot_e, nrow = 1)
dev.off()

# Original 2-panel plot
pdf("C:/Users/oofordile/Desktop/Figure_AgeGroup_Analysis_Aligned.pdf", width = 12, height = 6)
grid.arrange(plot_d, plot_e, nrow = 1)
dev.off()

png("C:/Users/oofordile/Desktop/Figure_AgeGroup_Analysis_Aligned.png",
    width = 12, height = 6, units = "in", res = 300)
grid.arrange(plot_d, plot_e, nrow = 1)
dev.off()

# Save demographic summary table
write.csv(demo_comparison, "C:/Users/oofordile/Desktop/Demographic_Comparison_Table.csv", row.names = FALSE)
write.csv(age_summary, "C:/Users/oofordile/Desktop/Age_Summary_by_Illness.csv", row.names = FALSE)

cat("\n\nAnalysis complete! Check console output and saved files.\n")
