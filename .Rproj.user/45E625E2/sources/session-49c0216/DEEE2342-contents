# Load libraries
library(dplyr)
library(Maaslin2)

# Paths to your data
ae_path <- "C:/Users/oofordile/Desktop/Adverse Events.csv"
microbiome_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
output_folder <- "C:/Users/oofordile/Desktop/maaslin2_output"

# Load AE data
ae_data <- read.csv(ae_path, stringsAsFactors = FALSE)

# Load microbiome data
microbiome_data <- read.csv(microbiome_path, stringsAsFactors = FALSE)

# Filter AE data for infection events only and get earliest infection day per child
ae_infections <- ae_data %>%
  filter(AE.Infection == 1) %>%
  group_by(Randomisation.No) %>%
  summarise(earliest_infection_day = min(Study.Day, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    Group = case_when(
      earliest_infection_day <= 35 ~ "Soon-Ill",
      earliest_infection_day > 35 & earliest_infection_day <= 70 ~ "Later-Ill",
      earliest_infection_day > 70 ~ "Much-Later-Ill"
    )
  )

# Join earliest infection group to microbiome samples
meta_joined <- microbiome_data %>%
  left_join(ae_infections %>% select(Randomisation.No, Group), by = "Randomisation.No") %>%
  mutate(
    SampleGroup = case_when(
      timepoints == "D1" & Group == "Soon-Ill" ~ "Soon-Ill",
      timepoints == "D1" & Group == "Later-Ill" ~ "Later-Ill",
      timepoints == "D1" & Group == "Much-Later-Ill" ~ "Much-Later-Ill",
      timepoints == "D85" & Ill_Status == "Not-Ill" ~ "D85 Not-Ill",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(SampleGroup)) %>%
  rename(age_group = X3agegroups.based.on.age.at.sampling)

# Print group counts
print(table(meta_joined$SampleGroup))

# Define top 10 taxa
top10_taxa_names <- c(
  "Prevotella_copri_35.15", "Faecalibacterium_prausnitzii_11.38", "Prevotella_stercorea_7.05",
  "Bacteroides__4.64", "Bifidobacterium__4.12", "Escherichia_coli_3.51",
  "Succinivibrio_dextrinosolvens_2.42", "Paraprevotella_xylaniphila_2.02",
  "Sutterella_wadsworthensis_1.75", "Streptococcus_salivarius_1.14"
)

# Extract feature and metadata tables
features <- meta_joined[, c("Randomisation.No", top10_taxa_names)]
rownames(features) <- features$Randomisation.No
features$Randomisation.No <- NULL

metadata <- meta_joined %>%
  filter(Randomisation.No %in% rownames(features)) %>%
  select(Randomisation.No, SampleGroup, age_group)
rownames(metadata) <- metadata$Randomisation.No
metadata$Randomisation.No <- NULL

# Ensure matching rownames
stopifnot(all(rownames(metadata) == rownames(features)))

# Ensure factor levels for SampleGroup
metadata$SampleGroup <- factor(metadata$SampleGroup, levels = c("D85 Not-Ill", "Soon-Ill", "Later-Ill", "Much-Later-Ill"))
metadata$age_group <- factor(metadata$age_group, levels = c("7to12 mths", "1to2years", "plus2years"))
# Run MaAsLin2
fit <- Maaslin2(
  input_data = features,
  input_metadata = metadata,
  output = output_folder,
  fixed_effects = c("SampleGroup", "age_group"),
  normalization = "TSS",
  transform = "NONE",
  min_prevalence = 0.1,
  standardize = TRUE,
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)
