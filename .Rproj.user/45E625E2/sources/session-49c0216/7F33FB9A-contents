library(tidyverse)
library(broom)

# Set file path
file_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"

# Read the data
data <- read_csv(file_path)

# Identify and remove rows with blank Ill_Status
data <- data %>% filter(!is.na(Ill_Status) & Ill_Status != "")

# Subset data for analysis
analysis_data <- data %>%
  select(`Randomisation No`, timepoints, Ill_Status, `3agegroups based on age at sampling`, 13:22)

# Rename age group column to simpler name
colnames(analysis_data)[4] <- "age_group"

# Convert to factors with explicit levels
analysis_data <- analysis_data %>%
  mutate(
    timepoints = factor(timepoints, levels = c("D1", "D15", "D85")),
    Ill_Status = factor(Ill_Status, levels = c("Not-Ill", "Ill")),
    age_group = factor(age_group, levels = c("7to12 mths", "1to2years", "plus2years"))
  )

# Get taxa names (columns 5 to 14)
taxa_names <- colnames(analysis_data)[5:14]

# Normalize counts to relative abundance within each sample (row)
analysis_data <- analysis_data %>%
  rowwise() %>%
  mutate(across(all_of(taxa_names), ~ .x / sum(c_across(all_of(taxa_names)), na.rm = TRUE))) %>%
  ungroup()

# Initialize list to store ANOVA results
anova_results_list <- list()

# Loop over taxa and run ANOVA on relative abundance
for(taxon in taxa_names) {
  # Formula with Ill_Status and interaction
  formula <- as.formula(paste(taxon, "~ Ill_Status * timepoints + age_group"))
  
  model <- lm(formula, data = analysis_data)
  
  # Get ANOVA table
  anova_res <- anova(model)
  
  # Convert to data frame and retain key statistics
  anova_df <- as.data.frame(anova_res) %>%
    rownames_to_column("term") %>%
    select(term, Df, `F value`, `Pr(>F)`) %>%
    rename(
      df = Df,
      F_value = `F value`,
      p.value = `Pr(>F)`
    ) %>%
    mutate(taxon = taxon)
  
  # Store result
  anova_results_list[[taxon]] <- anova_df
}

# Combine all results into one data frame
anova_results <- bind_rows(anova_results_list)

# Define terms of interest
terms_of_interest <- c("Ill_Status", "timepoints", "Ill_Status:timepoints", "age_group")

# Filter and adjust p-values using BH method
anova_results_filtered <- anova_results %>%
  filter(term %in% terms_of_interest) %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    significant = p.adj < 0.05,
    term_short = case_when(
      term == "Ill_Status" ~ "p_ill",
      term == "timepoints" ~ "p_time",
      term == "Ill_Status:timepoints" ~ "p_int",
      TRUE ~ term
    )
  )

# Extract and print all results for key terms
all_results <- anova_results_filtered %>%
  filter(term_short %in% c("p_ill", "p_time", "p_int")) %>%
  select(taxon, term_short, df, F_value, p.value, p.adj, significant)

# Print full results
print(all_results)

# Print only significant results
significant_results <- all_results %>% filter(significant == TRUE)
print(significant_results)
