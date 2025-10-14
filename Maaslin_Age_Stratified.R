# ======================================================
# Load libraries
# ======================================================
library(dplyr)
library(Maaslin2)

# ======================================================
# Paths
# ======================================================
microbiome_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
output_base <- "C:/Users/oofordile/Desktop/maaslin2_output_stratified"

# ======================================================
# Load microbiome data
# ======================================================
microbiome_data <- read.csv(microbiome_path, stringsAsFactors = FALSE)

# ======================================================
# Define top 10 taxa
# ======================================================
top10_taxa_names <- c(
  "Prevotella_copri_35.15", "Faecalibacterium_prausnitzii_11.38", "Prevotella_stercorea_7.05",
  "Bacteroides__4.64", "Bifidobacterium__4.12", "Escherichia_coli_3.51",
  "Succinivibrio_dextrinosolvens_2.42", "Paraprevotella_xylaniphila_2.02",
  "Sutterella_wadsworthensis_1.75", "Streptococcus_salivarius_1.14"
)

# ======================================================
# Function to run Maaslin2, compute CIs, and add Nature-ready text
# ======================================================
run_maaslin2_stratified_NatureText <- function(data, sample_filter, group_levels, output_folder, age_group_label) {
  
  filtered_data <- data %>%
    filter(SampleGroup %in% sample_filter) %>%
    distinct(Randomisation.No, .keep_all = TRUE)
  
  features <- filtered_data[, top10_taxa_names]
  rownames(features) <- filtered_data$Randomisation.No
  
  metadata <- filtered_data %>%
    select(Randomisation.No, SampleGroup) %>%
    filter(Randomisation.No %in% rownames(features))
  rownames(metadata) <- metadata$Randomisation.No
  
  stopifnot(all(rownames(metadata) == rownames(features)))
  
  metadata$SampleGroup <- factor(metadata$SampleGroup, levels = group_levels)
  
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  
  # Run Maaslin2
  fit <- Maaslin2(
    input_data = features,
    input_metadata = metadata,
    output = output_folder,
    fixed_effects = c("SampleGroup"),
    normalization = "TSS",
    transform = "NONE",
    min_prevalence = 0.1,
    standardize = TRUE,
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  )
  
  # Load results
  results <- read.csv(file.path(output_folder, "all_results.tsv"), sep = "\t")
  
  # Compute 95% CI and Nature-ready text
  results_filtered <- results %>%
    mutate(
      CI_lower = coef - 1.96 * stderr,
      CI_upper = coef + 1.96 * stderr
    ) %>%
    filter(metadata == "SampleGroup" & qval < 0.1) %>%
    mutate(
      Nature_Formatted_Stat = sprintf(
        "%s (%s, Age group: %s): coef = %.3f, 95%% CI [%.3f, %.3f], q = %.3f",
        feature, metadata, age_group_label, coef, CI_lower, CI_upper, qval
      )
    )
  
  # Export filtered CSV with Nature text
  write.csv(results_filtered,
            file.path(output_folder, "SampleGroup_significant_q_lt_0.1_with_CI_and_NatureText.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  
  return(fit)
}

# ======================================================
# Prepare metadata with updated SampleGroup
# ======================================================
meta_joined <- microbiome_data %>%
  rename(age_group = X3agegroups.based.on.age.at.sampling) %>%
  mutate(
    SampleGroup = case_when(
      timepoints == "D1"  & Ill_Status == "Ill"     ~ "D1 Ill",
      timepoints == "D15" & Ill_Status == "Ill"     ~ "D15 Ill",
      timepoints == "D85" & Ill_Status == "Not-Ill" ~ "D85 Not-Ill",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(SampleGroup))

print(table(meta_joined$SampleGroup, meta_joined$age_group))

# ======================================================
# Define pairwise comparisons of interest
# ======================================================
pairwise_comparisons <- list(
  list(groups = c("D1 Ill", "D85 Not-Ill"),  folder = file.path(output_base, "D1Ill_vs_D85NotIll"),  levels = c("D85 Not-Ill", "D1 Ill")),
  list(groups = c("D15 Ill", "D85 Not-Ill"), folder = file.path(output_base, "D15Ill_vs_D85NotIll"), levels = c("D85 Not-Ill", "D15 Ill"))
)

# ======================================================
# Run stratified analyses by age group
# ======================================================
for (comp in pairwise_comparisons) {
  for (age in c("7to12 mths", "1to2years", "plus2years")) {
    message("Running Maaslin2 for: ", paste(comp$groups, collapse = " vs "), " | Age group: ", age)
    
    subset_data <- meta_joined %>% filter(age_group == age)
    
    if (nrow(subset_data) > 0) {
      run_maaslin2_stratified_NatureText(
        data = subset_data,
        sample_filter = comp$groups,
        group_levels = comp$levels,
        output_folder = file.path(comp$folder, paste0("Age_", age)),
        age_group_label = age
      )
    } else {
      message("Skipping age group ", age, " (no samples).")
    }
  }
}

message("Stratified analyses (by age group) complete!")
