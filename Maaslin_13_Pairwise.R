# ======================================================
# Load libraries
# ======================================================
library(dplyr)
library(Maaslin2)

# ======================================================
# Path to microbiome data and output base folder
# ======================================================
microbiome_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
output_base <- "C:/Users/oofordile/Desktop/maaslin2_all_13_pairwise"

# ======================================================
# Load microbiome data
# ======================================================
microbiome_data <- read.csv(microbiome_path, stringsAsFactors = FALSE)

# ======================================================
# Define top 10 taxa column names
# ======================================================
top10_taxa_names <- c(
  "Prevotella_copri_35.15", "Faecalibacterium_prausnitzii_11.38", "Prevotella_stercorea_7.05",
  "Bacteroides__4.64", "Bifidobacterium__4.12", "Escherichia_coli_3.51",
  "Succinivibrio_dextrinosolvens_2.42", "Paraprevotella_xylaniphila_2.02",
  "Sutterella_wadsworthensis_1.75", "Streptococcus_salivarius_1.14"
)

# ======================================================
# Function to run Maaslin2 and save SampleGroup q<0.1 with CIs and Nature text
# ======================================================
run_maaslin2_with_CI_and_NatureText <- function(data, sample_filter, group_levels, output_folder) {
  
  filtered_data <- data %>%
    filter(SampleGroup %in% sample_filter) %>%
    distinct(Randomisation.No, .keep_all = TRUE)
  
  features <- filtered_data[, top10_taxa_names]
  rownames(features) <- filtered_data$Randomisation.No
  
  metadata <- filtered_data %>%
    select(Randomisation.No, SampleGroup, age_group) %>%
    filter(Randomisation.No %in% rownames(features))
  rownames(metadata) <- metadata$Randomisation.No
  
  stopifnot(all(rownames(metadata) == rownames(features)))
  
  metadata$SampleGroup <- factor(metadata$SampleGroup, levels = group_levels)
  metadata$age_group <- factor(metadata$age_group, levels = c("7to12 mths", "1to2years", "plus2years"))
  
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  
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
  
  # Load results and filter SampleGroup q<0.1
  results <- read.csv(file.path(output_folder, "all_results.tsv"), sep = "\t")
  
  results_filtered <- results %>%
    mutate(
      CI_lower = coef - 1.96 * stderr,
      CI_upper = coef + 1.96 * stderr
    ) %>%
    filter(metadata == "SampleGroup" & qval < 0.1) %>%
    # Add Nature-ready formatted statistic
    mutate(
      Nature_Formatted_Stat = sprintf(
        "%s (%s, Age-Adjusted): coef = %.3f, 95%% CI [%.3f, %.3f], q = %.3f",
        feature, metadata, coef, CI_lower, CI_upper, qval
      )
    )
  
  # Export CSV with UTF-8 encoding
  write.csv(results_filtered,
            file.path(output_folder, "SampleGroup_significant_q_lt_0.1_with_CI_and_NatureText.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  
  return(fit)
}

# ======================================================
# Prepare metadata with SampleGroup
# ======================================================
meta_joined <- microbiome_data %>%
  rename(age_group = X3agegroups.based.on.age.at.sampling) %>%
  mutate(
    SampleGroup = case_when(
      timepoints == "D1" & Ill_Status == "Ill" ~ "D1 Ill",
      timepoints == "D1" & Ill_Status == "Not-Ill" ~ "D1 Not-Ill",
      timepoints == "D15" & Ill_Status == "Ill" ~ "D15 Ill",
      timepoints == "D15" & Ill_Status == "Not-Ill" ~ "D15 Not-Ill",
      timepoints == "D85" & Ill_Status == "Ill" ~ "D85 Ill",
      timepoints == "D85" & Ill_Status == "Not-Ill" ~ "D85 Not-Ill",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(SampleGroup))

print(table(meta_joined$SampleGroup))

# ======================================================
# Define the 13 explicit pairwise comparisons
# ======================================================
pairwise_comparisons <- list(
  list(groups = c("D1 Not-Ill", "D15 Not-Ill"), folder = file.path(output_base, "D1_vs_D15_NotIll"), levels = c("D1 Not-Ill", "D15 Not-Ill")),
  list(groups = c("D15 Not-Ill", "D85 Not-Ill"), folder = file.path(output_base, "D15_vs_D85_NotIll"), levels = c("D15 Not-Ill", "D85 Not-Ill")),
  list(groups = c("D1 Not-Ill", "D85 Not-Ill"), folder = file.path(output_base, "D1_vs_D85_NotIll"), levels = c("D1 Not-Ill", "D85 Not-Ill")),
  list(groups = c("D1 Ill", "D15 Ill"), folder = file.path(output_base, "D1_vs_D15_Ill"), levels = c("D1 Ill", "D15 Ill")),
  list(groups = c("D15 Ill", "D85 Ill"), folder = file.path(output_base, "D15_vs_D85_Ill"), levels = c("D15 Ill", "D85 Ill")),
  list(groups = c("D1 Ill", "D85 Ill"), folder = file.path(output_base, "D1_vs_D85_Ill"), levels = c("D1 Ill", "D85 Ill")),
  list(groups = c("D1 Ill", "D1 Not-Ill"), folder = file.path(output_base, "D1_Ill_vs_NotIll"), levels = c("D1 Not-Ill", "D1 Ill")),
  list(groups = c("D15 Ill", "D15 Not-Ill"), folder = file.path(output_base, "D15_Ill_vs_NotIll"), levels = c("D15 Not-Ill", "D15 Ill")),
  list(groups = c("D85 Ill", "D85 Not-Ill"), folder = file.path(output_base, "D85_Ill_vs_NotIll"), levels = c("D85 Not-Ill", "D85 Ill")),
  list(groups = c("D1 Ill", "D85 Not-Ill"), folder = file.path(output_base, "D1_Ill_vs_D85_NotIll"), levels = c("D85 Not-Ill", "D1 Ill")),
  list(groups = c("D15 Ill", "D85 Not-Ill"), folder = file.path(output_base, "D15_Ill_vs_D85_NotIll"), levels = c("D85 Not-Ill", "D15 Ill")),
  list(groups = c("D1 Not-Ill", "D85 Ill"), folder = file.path(output_base, "D1_NotIll_vs_D85_Ill"), levels = c("D1 Not-Ill", "D85 Ill")),
  list(groups = c("D15 Not-Ill", "D85 Ill"), folder = file.path(output_base, "D15_NotIll_vs_D85_Ill"), levels = c("D15 Not-Ill", "D85 Ill"))
)

# ======================================================
# Run all 13 pairwise comparisons with Nature-formatted CSV
# ======================================================
for (comp in pairwise_comparisons) {
  message("Running Maaslin2 for: ", paste(comp$groups, collapse = " vs "))
  run_maaslin2_with_CI_and_NatureText(
    data = meta_joined,
    sample_filter = comp$groups,
    group_levels = comp$levels,
    output_folder = comp$folder
  )
}

message("All 13 pairwise Maaslin2 analyses complete with Nature-ready statistics!")
