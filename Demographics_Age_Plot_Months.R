# ======================================================
# Load Libraries
# ======================================================
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(knitr)
library(devEMF)  # for EMF export

# ======================================================
# Load Data
# ======================================================
data_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
df <- read_csv(data_path)

# Define correct age group levels
age_levels <- c("7to12 mths", "1to2years", "plus2years")

# Assign factor with exact levels
df$age_group <- factor(df$`3agegroups based on age at sampling`,
                       levels = age_levels)

# ======================================================
# DEMOGRAPHIC ANALYSIS
# ======================================================

cat("========================================\n")
cat("DEMOGRAPHIC ANALYSIS SUMMARY\n")
cat("========================================\n\n")

# Identify gender-related columns
gender_cols <- names(df)[grepl("gender|sex|female|male", names(df), ignore.case = TRUE)]
cat("Potential gender columns found:", paste(gender_cols, collapse = ", "), "\n\n")

# Filter for complete demographic data
df_demo <- df %>%
  filter(!is.na(Ill_Status), !is.na(age_group), 
         !is.na(`age at sampling`), !is.na(gender)) %>%
  mutate(
    Ill_Status = factor(Ill_Status, levels = c("Not-Ill", "Ill")),
    Gender = factor(gender, levels = c("Female", "Male"))
  )

# ======================================================
# 1. AGE DISTRIBUTION ANALYSIS
# ======================================================
cat("1. AGE DISTRIBUTION BY ILLNESS STATUS\n")
cat("=====================================\n")

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

# Mann–Whitney U test for age differences
age_test <- wilcox.test(`age at sampling` ~ Ill_Status, data = df_demo)
cat(sprintf("\nMann-Whitney U test for age differences between groups:\nW = %.2f, p-value = %.4f\n", 
            age_test$statistic, age_test$p.value))

if(age_test$p.value < 0.05) {
  cat("*** Significant age difference between Ill and Not-Ill groups ***\n")
} else {
  cat("No significant age difference between groups (age not a major confounder)\n")
}

# ======================================================
# 2. COMPREHENSIVE DEMOGRAPHIC COMPARISON TABLE
# ======================================================
cat("\n2. COMPREHENSIVE DEMOGRAPHIC COMPARISON TABLE\n")
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

# ======================================================
# 3. PLOT 1: Age Distribution by Illness Status
# ======================================================
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

# ======================================================
# EXPORTS
# ======================================================

# --- Plot Export (PDF + EMF) ---
pdf("C:/Users/oofordile/Desktop/Age_Distribution_by_Illness_Status.pdf",
    width = 8, height = 6)
print(plot_age_dist)
dev.off()

emf("C:/Users/oofordile/Desktop/Age_Distribution_by_Illness_Status.emf",
    width = 8, height = 6)
print(plot_age_dist)
dev.off()

# --- CSV Exports ---
write.csv(demo_comparison, 
          "C:/Users/oofordile/Desktop/Demographic_Comparison_Table.csv", 
          row.names = FALSE)

write.csv(age_summary, 
          "C:/Users/oofordile/Desktop/Age_Summary_by_Illness.csv", 
          row.names = FALSE)

cat("\nAnalysis complete. Plot (PDF + EMF) and CSVs saved to Desktop.\n")
