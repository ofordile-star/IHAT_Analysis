library(dplyr)
library(readr)

# Load data
data <- read_csv("C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv")

# ID column 
id_col <- "Randomisation No"

# Function to summarise females for a grouping variable
summarise_gender <- function(df, var, category_name, id_column){
  df %>%
    group_by({{var}}) %>%
    summarise(
      N_total = n_distinct(.data[[id_column]]),  # Count unique IDs, not samples
      N_female = n_distinct(.data[[id_column]][gender %in% c("Female","F","female","f",1)], na.rm = TRUE),
      Female_pct = round(100 * N_female / N_total, 1),
      .groups = "drop"
    ) %>%
    rename(Group = {{var}}) %>%
    mutate(
      Category = category_name,
      Group = as.character(Group)
    )
}

# Apply to each category
timepoints_summary <- summarise_gender(data, timepoints, "Timepoints", id_col)
adverse_summary   <- summarise_gender(data, Adverse_Events, "Adverse_Events", id_col)
ill_summary       <- summarise_gender(data, Ill_Status, "Ill_Status", id_col)

# Combine summaries
summary_table <- bind_rows(timepoints_summary,
                           adverse_summary,
                           ill_summary) %>%
  select(Category, Group, N_total, N_female, Female_pct)

# Add grand total
grand_total <- data %>%
  summarise(
    N_total = n_distinct(.data[[id_col]]),  # Count unique IDs, not samples
    N_female = n_distinct(.data[[id_col]][gender %in% c("Female","F","female","f",1)], na.rm = TRUE),
    Female_pct = round(100 * N_female / N_total, 1)
  ) %>%
  mutate(Category = "Overall",
         Group = "Total") %>%
  select(Category, Group, N_total, N_female, Female_pct)

# Final table
summary_table <- bind_rows(summary_table, grand_total)

summary_table
# Save as CSV on Desktop
write.csv(summary_table,
          "C:/Users/oofordile/Desktop/Summary_by_gender.csv",
          row.names = FALSE)
