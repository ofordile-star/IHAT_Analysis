# --- Load required libraries ---
library(dplyr)
library(tibble)
library(Maaslin2)

# --- Load your dataset ---
file_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
data <- read.csv(file_path, stringsAsFactors = FALSE)

# --- Remove rows with blank or NA Ill_Status ---
data <- data %>%
  filter(!is.na(Ill_Status) & Ill_Status != "")

# --- Prepare metadata for MaAsLin2 ---
meta_for_maaslin <- data %>%
  dplyr::select(Randomisation.No, Iron, Ill_Status, X3agegroups.based.on.age.at.sampling, timepoints) %>%
  mutate(SampleID = paste0(Randomisation.No, "_", timepoints)) %>%
  column_to_rownames("SampleID")

# --- Prepare taxa abundance matrix (columns 13 to 500) ---
taxa_for_maaslin <- data[, 13:500]
rownames(taxa_for_maaslin) <- paste0(data$Randomisation.No, "_", data$timepoints)

# --- Subset to only rows present in both metadata and taxa matrices ---
taxa_for_maaslin <- taxa_for_maaslin[rownames(meta_for_maaslin), ]

# --- Select top 50 taxa by mean abundance ---
taxa_means <- colMeans(taxa_for_maaslin)
top50_taxa <- names(sort(taxa_means, decreasing = TRUE))[1:50]
taxa_for_maaslin_top50 <- taxa_for_maaslin[, top50_taxa]

# --- Set factor levels and reference categories for categorical variables ---
meta_for_maaslin$X3agegroups.based.on.age.at.sampling <- factor(
  meta_for_maaslin$X3agegroups.based.on.age.at.sampling,
  levels = c("7to12 mths", "1to2years", "plus2years")
)

meta_for_maaslin$timepoints <- factor(
  meta_for_maaslin$timepoints,
  levels = c("D85", "D1", "D15")
)

meta_for_maaslin$Iron <- factor(meta_for_maaslin$Iron)
meta_for_maaslin$Ill_Status <- factor(meta_for_maaslin$Ill_Status)

# --- Optional sanity check: ensure rownames match exactly ---
stopifnot(all(rownames(meta_for_maaslin) == rownames(taxa_for_maaslin_top50)))

# --- Run MaAsLin2 ---
fit <- Maaslin2(
  input_data = taxa_for_maaslin_top50,
  input_metadata = meta_for_maaslin,
  output = "C:/Users/oofordile/Desktop/maaslin2_iron_illness_top50_output",
  fixed_effects = c("Iron", "Ill_Status", "X3agegroups.based.on.age.at.sampling", "timepoints"),
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  standardize = TRUE,
  correction = "BH",
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

# --- Export significant Iron treatment results (q < 0.1) with confidence intervals ---
# Read the all_results.tsv file from MaAsLin2 output
results <- read.table(
  "C:/Users/oofordile/Desktop/maaslin2_iron_illness_top50_output/all_results.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Check available column names
cat("\nAvailable columns in results:\n")
print(colnames(results))

# Calculate 95% confidence intervals
# CI = coefficient ± (1.96 × standard error)
results <- results %>%
  mutate(
    CI_lower = coef - (1.96 * stderr),
    CI_upper = coef + (1.96 * stderr)
  )

# Filter for Iron effects with q-value < 0.1
iron_significant <- results %>%
  filter(metadata == "Iron" & qval < 0.1) %>%
  arrange(qval)

# Reorder columns intelligently based on what exists
# Place CI columns after coefficient
base_cols <- c("feature", "metadata", "value", "coef", "CI_lower", "CI_upper", "stderr")
other_cols <- setdiff(colnames(iron_significant), c(base_cols, "CI_lower", "CI_upper"))
col_order <- c(base_cols, other_cols)
col_order <- col_order[col_order %in% colnames(iron_significant)]

iron_significant <- iron_significant %>%
  select(all_of(col_order))

# Export to CSV
write.csv(
  iron_significant,
  "C:/Users/oofordile/Desktop/Iron_significant_results_q0.1.csv",
  row.names = FALSE
)

# Print summary
cat("\n=== Summary of Significant Iron Treatment Results (q < 0.1) ===\n")
cat("Number of significant associations:", nrow(iron_significant), "\n")
if (nrow(iron_significant) > 0) {
  cat("\nTop significant taxa:\n")
  print(iron_significant[1:min(10, nrow(iron_significant)), 
                         c("feature", "coef", "CI_lower", "CI_upper", "pval", "qval")])
} else {
  cat("No significant associations found with q < 0.1\n")
}
cat("\nResults exported to: C:/Users/oofordile/Desktop/Iron_significant_results_q0.1.csv\n")