# --- Load required libraries ---
library(dplyr)
library(tibble)
library(Maaslin2)

# --- Load your dataset ---
file_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
data <- read.csv(file_path, stringsAsFactors = FALSE)

# --- Remove rows with blank or NA Ill_Status ---
data <- data %>% filter(!is.na(Ill_Status) & Ill_Status != "")

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
