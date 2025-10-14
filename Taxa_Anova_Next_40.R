# ======================================================
# Load Required Libraries
# ======================================================
library(tidyverse)
library(broom)
library(effectsize)

# ======================================================
# Set File Paths
# ======================================================
file_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
csv_path <- "C:/Users/oofordile/Desktop/ANOVA_Ill_Status_Significant.csv"

# ======================================================
# Read and Clean Data
# ======================================================
data <- read_csv(file_path)

# Remove rows with blank Ill_Status
data <- data %>% filter(!is.na(Ill_Status) & Ill_Status != "")

# Subset data for analysis (columns 23-62)
analysis_data <- data %>%
  select(`Randomisation No`, timepoints, Ill_Status, `3agegroups based on age at sampling`, 23:62)

# Rename age group column
colnames(analysis_data)[4] <- "age_group"

# Convert to factors with explicit levels
analysis_data <- analysis_data %>%
  mutate(
    timepoints = factor(timepoints, levels = c("D1", "D15", "D85")),
    Ill_Status = factor(Ill_Status, levels = c("Not-Ill", "Ill")),
    age_group = factor(age_group, levels = c("7to12 mths", "1to2years", "plus2years"))
  )

# Get taxa names (columns 5 onwards)
taxa_names <- colnames(analysis_data)[5:ncol(analysis_data)]

# Normalize counts to relative abundance per sample
analysis_data <- analysis_data %>%
  rowwise() %>%
  mutate(row_sum = sum(c_across(all_of(taxa_names)), na.rm = TRUE)) %>%
  mutate(across(all_of(taxa_names), ~ ifelse(row_sum == 0, 0, .x / row_sum))) %>%
  select(-row_sum) %>%
  ungroup()

# ======================================================
# Initialize lists to store results
# ======================================================
anova_results_list <- list()
reg_list <- list()

# ======================================================
# Loop over taxa: ANOVA and regression
# ======================================================
for(taxon in taxa_names) {
  tryCatch({
    # Skip taxa with no variation
    if(var(analysis_data[[taxon]], na.rm = TRUE) == 0) {
      cat("Skipping", taxon, "- no variation in data\n")
      next
    }
    
    # Full model
    formula <- as.formula(paste0("`", taxon, "` ~ Ill_Status * timepoints + age_group"))
    model <- lm(formula, data = analysis_data)
    
    # ANOVA results
    anova_res_full <- as.data.frame(anova(model)) %>% rownames_to_column("term")
    anova_df <- anova_res_full %>%
      select(term, Df, `F value`, `Pr(>F)`) %>%
      rename(
        df = Df,
        F_value = `F value`,
        p.value = `Pr(>F)`
      ) %>%
      mutate(taxon = taxon)
    anova_results_list[[taxon]] <- anova_df
    
    # Age-adjusted effect for Ill_Status
    lm_group <- lm(as.numeric(analysis_data[[taxon]]) ~ Ill_Status + age_group, data = analysis_data)
    lm_summary <- tidy(lm_group, conf.int = TRUE)
    
    # Partial eta-squared
    ss_ill <- anova_res_full[anova_res_full$term == "Ill_Status", "Sum Sq"][[1]]
    ss_residual <- anova_res_full[anova_res_full$term == "Residuals", "Sum Sq"][[1]]
    eta2_partial <- ss_ill / (ss_ill + ss_residual)
    
    f_stat <- anova_res_full[anova_res_full$term == "Ill_Status", "F value"][[1]]
    df1 <- anova_res_full[anova_res_full$term == "Ill_Status", "Df"][[1]]
    df2 <- anova_res_full[anova_res_full$term == "Residuals", "Df"][[1]]
    
    # CI for partial eta-squared
    tryCatch({
      eta2_ci <- effectsize::F_to_eta2(f_stat, df1, df2, ci = 0.95)
      eta2_low <- eta2_ci$CI_low
      eta2_high <- eta2_ci$CI_high
    }, error = function(e) {
      eta2_low <- max(0, eta2_partial - 0.05)
      eta2_high <- min(1, eta2_partial + 0.05)
    })
    
    # Slope and CI
    slope <- lm_summary$estimate[2]
    ci_low_slope <- lm_summary$conf.low[2]
    ci_high_slope <- lm_summary$conf.high[2]
    
    # R-squared
    r2_group <- glance(lm_group)$r.squared
    
    reg_stats <- data.frame(
      Taxon = taxon,
      Effect = "Ill_Status (age-adj)",
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
    
  }, error = function(e) {
    cat("Error processing", taxon, ":", e$message, "\n")
  })
}

# ======================================================
# Combine ANOVA results
# ======================================================
anova_results <- bind_rows(anova_results_list)

# --- BH adjustment applied across all taxa together (Script 1 style) ---
terms_of_interest <- c("Ill_Status", "timepoints")
anova_results_filtered <- anova_results %>%
  filter(term %in% terms_of_interest) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))  # <-- no grouping by term

# Pivot wider for reporting
anova_combined <- anova_results_filtered %>%
  pivot_wider(
    id_cols = taxon,
    names_from = term,
    values_from = c(F_value, df, p.value, p.adj),
    names_sep = "_"
  )

combined_results <- anova_combined %>%
  transmute(
    Taxon = taxon,
    `Ill_Status Effect: F(df), adj. p` = paste0(
      round(F_value_Ill_Status, 3),
      "(", df_Ill_Status, "), ",
      signif(p.adj_Ill_Status, 4)
    ),
    `Time Effect: F(df), adj. p` = paste0(
      round(F_value_timepoints, 3),
      "(", df_timepoints, "), ",
      signif(p.adj_timepoints, 4)
    )
  )

# --- Extract Ill_Status adjusted p-values for merging
ill_status_adj_p <- anova_results_filtered %>%
  filter(term == "Ill_Status") %>%
  select(taxon, p.adj, F_value, df) %>%
  rename(Taxon = taxon, p_adj_ill = p.adj, F_ill = F_value, df_ill = df)

# --- Combine regression results with adjusted p-values
reg_combined <- bind_rows(reg_list) %>%
  left_join(ill_status_adj_p, by = "Taxon")

# --- Filter significant Ill_Status taxa
sig_taxa <- ill_status_adj_p %>% filter(p_adj_ill < 0.05) %>% pull(Taxon)
cat("Number of significant Ill_Status taxa (adj p < 0.05):", length(sig_taxa), "\n")

# --- Nature-style report for significant taxa
reg_combined_sig <- reg_combined %>%
  filter(Taxon %in% sig_taxa) %>%
  mutate(
    Nature_Report = paste0(
      "Ill_Status (age-adj): F(", df_ill, ") = ", round(F_ill, 2),
      ", adj p = ", signif(p_adj_ill, 3),
      ", eta2_p = ", round(Eta2_partial, 3),
      ", 95% CI [", round(Eta2_CI_low, 3), ", ", round(Eta2_CI_high, 3), "]",
      "; slope = ", round(Slope, 3),
      " (95% CI [", round(Slope_CI_low, 3), ", ", round(Slope_CI_high, 3), "])"
    )
  )

# Merge with combined ANOVA results
final_csv <- combined_results %>%
  filter(Taxon %in% sig_taxa) %>%
  left_join(
    reg_combined_sig %>% 
      select(Taxon, Slope, Slope_CI_low, Slope_CI_high, R2, Eta2_partial, Eta2_CI_low, Eta2_CI_high, Nature_Report),
    by = "Taxon"
  )

# ======================================================
# Write CSV
# ======================================================
write_csv(final_csv, csv_path)
cat("CSV saved to:", csv_path, "\n")

# Print summary
print(final_csv)
