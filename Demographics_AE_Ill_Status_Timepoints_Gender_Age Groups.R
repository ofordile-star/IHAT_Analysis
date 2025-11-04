# Load libraries
library(dplyr)
library(readr)
library(tidyr)

# Load data
data <- read_csv("C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv")

# Check columns
print("Available columns:")
print(colnames(data))

# Define key columns
age_col <- "3agegroups based on age at sampling"
id_col <- "Randomisation No"

# Define correct age group order (youngest first)
age_groups <- c("7to12 mths", "1to2years", "plus2years")

# Function: compute female numerator/denominator and formatted % string
calculate_female_stats <- function(df, filter_var, filter_value,
                                   age_group_value = NULL, age_group_col, id_column,
                                   show_fraction = FALSE) {
  if (!is.null(age_group_value)) {
    filtered_df <- df %>%
      filter(.data[[age_group_col]] == age_group_value,
             .data[[filter_var]] == filter_value)
  } else {
    filtered_df <- df %>%
      filter(.data[[filter_var]] == filter_value)
  }
  
  n_total <- n_distinct(filtered_df[[id_column]])
  n_female <- n_distinct(filtered_df[[id_column]][filtered_df$gender %in% c("Female","F","female","f",1)])
  
  if (n_total > 0) {
    pct <- round(100 * n_female / n_total, 1)
    if (show_fraction) {
      return(paste0(n_female, "/", n_total, " (", pct, ")"))
    } else {
      return(paste0(n_total, " (", pct, ")"))
    }
  } else {
    if (show_fraction) {
      return("0/0 (0.0)")
    } else {
      return("0 (0.0)")
    }
  }
}

# Initialize results
results <- data.frame()

# Loop through each defined age group
for (age in age_groups) {
  age_data <- data %>% filter(.data[[age_col]] == age)
  if (nrow(age_data) == 0) next
  
  n_total <- n_distinct(age_data[[id_col]])
  n_female <- n_distinct(age_data[[id_col]][age_data$gender %in% c("Female","F","female","f",1)])
  female_pct <- round(100 * n_female / n_total, 1)
  
  # Compute female fractions for Ill and Not-Ill (with numerator/denominator)
  ill_stats <- calculate_female_stats(data, "Ill_Status", "Ill", age, age_col, id_col, show_fraction = TRUE)
  notill_stats <- calculate_female_stats(data, "Ill_Status", "Not-Ill", age, age_col, id_col, show_fraction = TRUE)
  
  # Timepoints (regular % form)
  d1_stats <- calculate_female_stats(data, "timepoints", "D1", age, age_col, id_col)
  d15_stats <- calculate_female_stats(data, "timepoints", "D15", age, age_col, id_col)
  d85_stats <- calculate_female_stats(data, "timepoints", "D85", age, age_col, id_col)
  
  # Missing Ill_Status count
  na_count <- n_distinct(age_data[[id_col]][is.na(age_data$Ill_Status)])
  
  row <- data.frame(
    `Age group` = age,
    `N total` = n_total,
    `F (female)` = n_female,
    `F %` = female_pct,
    `D1 (F %)` = d1_stats,
    `D15 (F %)` = d15_stats,
    `D85 (F %)` = d85_stats,
    `Ill (F %)` = ill_stats,
    `Not-Ill (F %)` = notill_stats,
    `NA*` = na_count,
    check.names = FALSE
  )
  
  results <- bind_rows(results, row)
}

# Add "All" row (with same formatting)
all_n_total <- n_distinct(data[[id_col]])
all_n_female <- n_distinct(data[[id_col]][data$gender %in% c("Female","F","female","f",1)])
all_female_pct <- round(100 * all_n_female / all_n_total, 1)
all_na_count <- n_distinct(data[[id_col]][is.na(data$Ill_Status)])

all_row <- data.frame(
  `Age group` = "All",
  `N total` = all_n_total,
  `F (female)` = all_n_female,
  `F %` = all_female_pct,
  `D1 (F %)` = calculate_female_stats(data, "timepoints", "D1", NULL, age_col, id_col),
  `D15 (F %)` = calculate_female_stats(data, "timepoints", "D15", NULL, age_col, id_col),
  `D85 (F %)` = calculate_female_stats(data, "timepoints", "D85", NULL, age_col, id_col),
  `Ill (F %)` = calculate_female_stats(data, "Ill_Status", "Ill", NULL, age_col, id_col, show_fraction = TRUE),
  `Not-Ill (F %)` = calculate_female_stats(data, "Ill_Status", "Not-Ill", NULL, age_col, id_col, show_fraction = TRUE),
  `NA*` = all_na_count,
  check.names = FALSE
)

results <- bind_rows(results, all_row)

# Export CSV
output_path <- "C:/Users/oofordile/Desktop/Summary_Gender_Table.csv"
write_csv(results, output_path)

print(paste("✅ CSV file created successfully at:", output_path))
print(results)
