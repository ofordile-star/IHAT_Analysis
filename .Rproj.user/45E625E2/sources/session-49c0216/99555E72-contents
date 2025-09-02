library(dplyr)
library(readr)

# Load data
data <- read_csv("C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv")

# Check what columns are available
print("Available columns:")
print(colnames(data))

# Choose age group column
age_col <- "3agegroups based on age at sampling"

# ID column identified from your dataset
id_col <- "Randomisation No"

# Function to summarise females for a grouping variable, stratified by age
summarise_gender_age <- function(df, var, category_name, age_group_col, id_column){
  df %>%
    group_by(.data[[age_group_col]], {{var}}) %>%
    summarise(
      N_total = n_distinct(.data[[id_column]]),  # Count unique IDs, not samples
      N_female = n_distinct(.data[[id_column]][gender %in% c("Female","F","female","f",1)], na.rm = TRUE),
      Female_pct = round(100 * N_female / N_total, 1),
      .groups = "drop"
    ) %>%
    rename(Age_Group = .data[[age_group_col]],
           Group = {{var}}) %>%
    mutate(Category = category_name,
           Group = as.character(Group))
}

# Apply function to each category
timepoints_summary <- summarise_gender_age(data, timepoints, "Timepoints", age_col, id_col)
adverse_summary   <- summarise_gender_age(data, Adverse_Events, "Adverse_Events", age_col, id_col)
ill_summary       <- summarise_gender_age(data, Ill_Status, "Ill_Status", age_col, id_col)

# Combine
summary_table_age <- bind_rows(timepoints_summary,
                               adverse_summary,
                               ill_summary) %>%
  select(Age_Group, Category, Group, N_total, N_female, Female_pct)

# --- Add grand total per age group ---
grand_totals_age <- data %>%
  group_by(.data[[age_col]]) %>%
  summarise(
    N_total = n_distinct(.data[[id_col]]),  # Count unique IDs, not samples
    N_female = n_distinct(.data[[id_col]][gender %in% c("Female","F","female","f",1)], na.rm = TRUE),
    Female_pct = round(100 * N_female / N_total, 1),
    .groups = "drop"
  ) %>%
  mutate(Category = "Overall",
         Group = "Total") %>%
  rename(Age_Group = .data[[age_col]]) %>%
  select(Age_Group, Category, Group, N_total, N_female, Female_pct)

# Final stratified table
summary_table_age <- bind_rows(summary_table_age, grand_totals_age)

# --- Export CSV ---
write.csv(summary_table_age,
          "C:/Users/oofordile/Desktop/Summary_by_gender_age.csv",
          row.names = FALSE)

# --- Word-friendly table ---
# Format percentages as % sign for Word
summary_table_age_word <- summary_table_age %>%
  mutate(Female_pct = paste0(Female_pct, "%"))

summary_table_age_word
