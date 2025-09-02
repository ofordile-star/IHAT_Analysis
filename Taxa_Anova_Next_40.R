library(tidyverse)
library(broom)

# --- Set file path ---
file_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"

# --- Read the data ---
data <- read_csv(file_path)

# --- Remove rows with blank Ill_Status ---
data <- data %>% filter(!is.na(Ill_Status) & Ill_Status != "")

# --- Subset data for analysis (columns 23-62) ---
analysis_data <- data %>%
  select(`Randomisation No`, timepoints, Ill_Status, `3agegroups based on age at sampling`, 23:62)

# --- Rename age group column ---
colnames(analysis_data)[4] <- "age_group"

# --- Convert to factors with explicit levels ---
analysis_data <- analysis_data %>%
  mutate(
    timepoints = factor(timepoints, levels = c("D1", "D15", "D85")),
    Ill_Status = factor(Ill_Status, levels = c("Not-Ill", "Ill")),
    age_group = factor(age_group, levels = c("7to12 mths", "1to2years", "plus2years"))
  )

# --- Get taxa names (columns 5 onwards after selection) ---
taxa_names <- colnames(analysis_data)[5:ncol(analysis_data)]

# --- Normalize counts to relative abundance per sample ---
analysis_data <- analysis_data %>%
  rowwise() %>%
  mutate(
    row_sum = sum(c_across(all_of(taxa_names)), na.rm = TRUE)
  ) %>%
  mutate(
    across(all_of(taxa_names), ~ ifelse(row_sum == 0, 0, .x / row_sum))
  ) %>%
  select(-row_sum) %>%
  ungroup()

# --- Initialize list to store ANOVA results ---
anova_results_list <- list()

# --- Loop over taxa and run ANOVA safely ---
for(taxon in taxa_names) {
  tryCatch({
    # Check if there's variation in the data
    if(var(analysis_data[[taxon]], na.rm = TRUE) == 0) {
      cat("Skipping", taxon, "- no variation in data\n")
      next
    }
    
    formula <- as.formula(paste0("`", taxon, "` ~ Ill_Status * timepoints + age_group"))
    
    model <- lm(formula, data = analysis_data)
    
    anova_res <- anova(model)
    
    anova_df <- as.data.frame(anova_res) %>%
      rownames_to_column("term") %>%
      select(term, Df, `F value`, `Pr(>F)`) %>%
      rename(
        df = Df,
        F_value = `F value`,
        p.value = `Pr(>F)`
      ) %>%
      mutate(taxon = taxon)
    
    anova_results_list[[taxon]] <- anova_df
    
  }, error = function(e) {
    cat("Error processing", taxon, ":", e$message, "\n")
  })
}

# --- Combine all results ---
anova_results <- bind_rows(anova_results_list)

# --- Define terms of interest and adjust p-values ---
terms_of_interest <- c("Ill_Status", "timepoints")

anova_results_filtered <- anova_results %>%
  filter(term %in% terms_of_interest) %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    significant = p.adj < 0.05,
    term_short = case_when(
      term == "Ill_Status" ~ "p_ill",
      term == "timepoints" ~ "p_time"
    )
  )

# --- Extract only significant adjusted p-values for Ill_Status and timepoints ---
adjusted_p_significant <- anova_results_filtered %>%
  filter(significant == TRUE) %>%
  select(taxon, term_short, p.adj)

# --- Print significant adjusted p-values ---
print(adjusted_p_significant)
