# ======================================================
# Load libraries
# ======================================================
library(dplyr)
library(Maaslin2)

# ======================================================
# File paths
# ======================================================
ae_path <- "C:/Users/oofordile/Desktop/Adverse Events.csv"
microbiome_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
output_folder <- "C:/Users/oofordile/Desktop/maaslin2_Residual"

# ======================================================
# Load data
# ======================================================
ae_data <- read.csv(ae_path, stringsAsFactors = FALSE)
microbiome_data <- read.csv(microbiome_path, stringsAsFactors = FALSE)

# ======================================================
# Prepare AE data: Infection events only
# ======================================================
ae_infections <- ae_data %>%
  filter(AE.Infection == 1) %>%
  group_by(Randomisation.No) %>%
  summarise(
    earliest_day = min(Study.Day, na.rm = TRUE),
    latest_day = max(Study.Day, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    Group = case_when(
      earliest_day >= 40 & earliest_day <= 75 & latest_day <= 75 ~ "Recent-Ill",
      earliest_day <= 35 & latest_day <= 35 ~ "Early-Ill",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group))

# ======================================================
# Merge with microbiome metadata (D85 samples)
# ======================================================
meta_joined <- microbiome_data %>%
  filter(timepoints == "D85") %>%
  left_join(ae_infections %>% select(Randomisation.No, Group), by = "Randomisation.No") %>%
  mutate(
    SampleGroup = case_when(
      Group == "Recent-Ill" ~ "Recent-Ill",
      Group == "Early-Ill" ~ "Early-Ill",
      is.na(Group) & Ill_Status == "Not-Ill" ~ "D85 Not-Ill",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(SampleGroup)) %>%
  rename(age_group = X3agegroups.based.on.age.at.sampling)

print(table(meta_joined$SampleGroup))

# ======================================================
# Top 10 taxa
# ======================================================
top10_taxa_names <- c(
  "Prevotella_copri_35.15", "Faecalibacterium_prausnitzii_11.38", "Prevotella_stercorea_7.05",
  "Bacteroides__4.64", "Bifidobacterium__4.12", "Escherichia_coli_3.51",
  "Succinivibrio_dextrinosolvens_2.42", "Paraprevotella_xylaniphila_2.02",
  "Sutterella_wadsworthensis_1.75", "Streptococcus_salivarius_1.14"
)

# ======================================================
# Function to run Maaslin2, compute CI, and add Nature text
# ======================================================
run_maaslin_with_NatureText <- function(group1, group2, out_name) {
  
  subset_meta <- meta_joined %>%
    filter(SampleGroup %in% c(group1, group2))
  
  features <- subset_meta[, c("Randomisation.No", top10_taxa_names)]
  rownames(features) <- features$Randomisation.No
  features$Randomisation.No <- NULL
  
  metadata <- subset_meta %>%
    select(Randomisation.No, SampleGroup, age_group)
  rownames(metadata) <- metadata$Randomisation.No
  metadata$Randomisation.No <- NULL
  
  metadata$SampleGroup <- factor(metadata$SampleGroup, levels = c(group2, group1)) # reference = group2
  metadata$age_group <- factor(metadata$age_group, levels = c("7to12 mths", "1to2years", "plus2years"))
  
  stopifnot(all(rownames(metadata) == rownames(features)))
  
  # Run Maaslin2
  Maaslin2(
    input_data = features,
    input_metadata = metadata,
    output = file.path(output_folder, out_name),
    fixed_effects = c("SampleGroup", "age_group"),
    normalization = "TSS",
    transform = "NONE",
    min_prevalence = 0.1,
    standardize = TRUE,
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  )
  
  # Load results
  results <- read.csv(file.path(output_folder, out_name, "all_results.tsv"), sep = "\t")
  
  # Compute 95% CI and Nature-ready text
  results_filtered <- results %>%
    mutate(
      CI_lower = coef - 1.96 * stderr,
      CI_upper = coef + 1.96 * stderr,
      Nature_Formatted_Stat = sprintf(
        "%s: coef = %.3f, 95%% CI [%.3f, %.3f], p = %.3f",
        feature, coef, CI_lower, CI_upper, pval
      )
    ) %>%
    filter(metadata == "SampleGroup")
  
  # Export CSV
  write.csv(
    results_filtered,
    file.path(output_folder, out_name, "SampleGroup_with_CI_and_NatureText.csv"),
    row.names = FALSE, fileEncoding = "UTF-8"
  )
  
  message("Exported CSV with Nature text for ", out_name)
}

# ======================================================
# Ensure output folder exists
# ======================================================
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# ======================================================
# Run comparisons
# ======================================================
run_maaslin_with_NatureText("Recent-Ill", "D85 Not-Ill", "D85_RecentIll_vs_NotIll")
run_maaslin_with_NatureText("Early-Ill", "D85 Not-Ill", "D85_EarlyIll_vs_NotIll")

message("Residual Maaslin2 analyses complete with Nature-ready text!")
