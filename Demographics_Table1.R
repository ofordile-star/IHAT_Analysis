# ======================================================
# COHORT DEMOGRAPHICS: MAIN + SUPPLEMENTARY TABLES
# ======================================================

library(dplyr)
library(readr)
library(tidyr)
library(officer)
library(flextable)

# ======================================================
# Load Data
# ======================================================
data_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
df <- read_csv(data_path)

# Harmonize column names
names(df) <- gsub("\\.", "_", names(df))

# Identify correct age column
age_col <- grep("age.at.sampling|age.at.enrolment", names(df), ignore.case = TRUE, value = TRUE)[1]
age_group_col <- grep("3agegroups|X3agegroups", names(df), ignore.case = TRUE, value = TRUE)[1]

# ======================================================
# Clean Working Copy with AE handling
# ======================================================
df_clean <- df %>%
  rename(
    ID = `Randomisation No`,
    Gender = gender,
    Location = location,
    Ill_Status = Ill_Status,
    AE = Adverse_Events
  ) %>%
  mutate(
    AE = as.numeric(AE),
    Ill_Status = factor(Ill_Status, levels = c("Not-Ill", "Ill")),
    Gender = factor(Gender, levels = c("Female", "Male")),
    Age_months = !!sym(age_col),
    age_group = factor(!!sym(age_group_col),
                       levels = c("7to12 mths", "1to2years", "plus2years")),
    # Define Ill_Group including Non-infectious AE
    Ill_Group = case_when(
      !is.na(Ill_Status) ~ as.character(Ill_Status),
      is.na(Ill_Status) & AE == 1 ~ "Non-infectious AE",
      TRUE ~ NA_character_
    ),
    Ill_Group = factor(Ill_Group, levels = c("Not-Ill", "Ill", "Non-infectious AE"))
  ) %>%
  filter(!is.na(age_group), !is.na(Gender))

# ======================================================
# Standardize location names
# ======================================================
df_clean <- df_clean %>%
  mutate(Location = case_when(
    Location == "CHAMOI" ~ "Chamoi",
    TRUE ~ tools::toTitleCase(tolower(Location))
  ))

# ======================================================
# Diagnostic Summary
# ======================================================
cat("========================================\n")
cat("DATA OVERVIEW\n")
cat("========================================\n")
cat(sprintf("Total samples: %d\n", nrow(df_clean)))
cat(sprintf("Unique children: %d\n", n_distinct(df_clean$ID)))
cat(sprintf("Unique locations: %d\n\n", n_distinct(df_clean$Location)))

# ======================================================
# MAIN DEMOGRAPHICS TABLE (Publication)
# ======================================================
# Core summary by Age × Illness
main_table <- df_clean %>%
  group_by(age_group, Ill_Group) %>%
  summarise(
    N_children = n_distinct(ID),
    N_Female = n_distinct(ID[Gender == "Female"]),
    Mean_Age_mo = round(mean(Age_months[!duplicated(ID)], na.rm = TRUE),1),
    SD_Age_mo = round(sd(Age_months[!duplicated(ID)], na.rm = TRUE),1),
    .groups = "drop"
  ) %>%
  mutate(
    Pct_Female = round(100 * N_Female / N_children, 1),
    Female_display = sprintf("%d (%.1f%%)", N_Female, Pct_Female),
    Age_display = sprintf("%.1f ± %.1f", Mean_Age_mo, SD_Age_mo)
  ) %>%
  select(
    `Age Group` = age_group,
    `Illness Status` = Ill_Group,
    `N Children` = N_children,
    `Female (%%)` = Female_display,
    `Age (months)` = Age_display
  )

# ======================================================
# ADD GEOGRAPHIC REGION COLUMNS
# ======================================================
region_summary <- df_clean %>%
  group_by(age_group, Ill_Group, Location) %>%
  summarise(
    N_children_loc = n_distinct(ID),
    N_Female_loc = n_distinct(ID[Gender == "Female"]),
    .groups = "drop"
  ) %>%
  mutate(
    FemalePct_loc = round(100 * N_Female_loc / N_children_loc, 1),
    Display = sprintf("%d (%.1f%%)", N_Female_loc, FemalePct_loc)
  ) %>%
  select(age_group, Ill_Group, Location, Display) %>%
  pivot_wider(names_from = Location, values_from = Display, values_fill = "0 (0%)")

# Merge geographic data
main_table <- main_table %>%
  left_join(region_summary, by = c("Age Group" = "age_group", "Illness Status" = "Ill_Group"))

# Reorder columns: main + 4 geographic columns
geo_cols <- setdiff(names(main_table), c("Age Group", "Illness Status", "N Children", "Female (%%)", "Age (months)"))
main_table <- main_table %>%
  select(`Age Group`, `Illness Status`, `N Children`, `Female (%%)`, `Age (months)`, all_of(sort(geo_cols)))

# ======================================================
# SUPPLEMENTARY TABLE (unchanged)
# ======================================================
if(!"timepoints" %in% names(df_clean)) stop("No 'timepoints' column in dataset!")

supp_table_full <- df_clean %>%
  group_by(age_group, Ill_Group, Location, timepoints) %>%
  summarise(
    N_children_tp = n_distinct(ID),
    N_samples_tp = n(),
    N_Female_tp = n_distinct(ID[Gender == "Female"]),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = timepoints,
    values_from = c(N_children_tp, N_samples_tp, N_Female_tp),
    values_fill = 0
  )

all_combinations <- expand_grid(
  age_group = levels(df_clean$age_group),
  Ill_Group = levels(df_clean$Ill_Group),
  Location = unique(df_clean$Location)
)

supp_table_full <- all_combinations %>%
  left_join(supp_table_full, by = c("age_group", "Ill_Group", "Location")) %>%
  mutate(
    N_children = pmax(N_children_tp_D1, N_children_tp_D15, N_children_tp_D85, na.rm = TRUE),
    N_Female = pmax(N_Female_tp_D1, N_Female_tp_D15, N_Female_tp_D85, na.rm = TRUE),
    Pct_Female = round(100 * N_Female / N_children, 1),
    Female_display = sprintf("%d (%.1f%%)", N_Female, Pct_Female)
  ) %>%
  filter(!(Ill_Group == "Non-infectious AE" & N_children == 0)) %>%
  select(
    `Age Group` = age_group,
    `Illness Status` = Ill_Group,
    `Location`,
    `N Children` = N_children,
    `Female (%%)` = Female_display,
    `Samples D1` = N_samples_tp_D1,
    `Samples D15` = N_samples_tp_D15,
    `Samples D85` = N_samples_tp_D85
  )

# ======================================================
# EXPORTS
# ======================================================
write.csv(main_table, "C:/Users/oofordile/Desktop/Cohort_Demographics_Main.csv", row.names = FALSE)
write.csv(supp_table_full, "C:/Users/oofordile/Desktop/Cohort_Demographics_Supplementary.csv", row.names = FALSE)

make_flex <- function(tbl, title) {
  flextable(tbl) %>%
    autofit() %>%
    bold(part = "header") %>%
    set_caption(title) %>%
    theme_vanilla()
}

doc_main <- read_docx() %>%
  body_add_par("Main Demographics Table", style = "heading 1") %>%
  body_add_flextable(make_flex(main_table, "Main Cohort Demographics by Age, Illness, and Region"))

doc_supp <- read_docx() %>%
  body_add_par("Supplementary Demographics Table", style = "heading 1") %>%
  body_add_flextable(make_flex(supp_table_full, "Detailed Cohort Demographics by Age, Illness, Location, and Timepoints"))

print(doc_main, target = "C:/Users/oofordile/Desktop/Cohort_Demographics_Main.docx")
print(doc_supp, target = "C:/Users/oofordile/Desktop/Cohort_Demographics_Supplementary.docx")

cat("========================================\n")
cat("COMPLETED: Tables exported to Desktop\n")
cat(" - Cohort_Demographics_Main.csv / .docx (includes 4 geographic region columns)\n")
cat(" - Cohort_Demographics_Supplementary.csv / .docx\n")
cat("========================================\n")
