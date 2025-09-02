# ======================================================
# Load Required Libraries
# ======================================================
library(phyloseq)
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(ggpubr)
library(devEMF)  # for EMF export
library(grid)    # for unit() used in tick length

# ======================================================
# Load Data
# ======================================================
file_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
merged_data <- read.csv(file_path, stringsAsFactors = FALSE)

# --- Remove NA or blank Ill_Status ---
merged_data <- merged_data %>% filter(!is.na(Ill_Status) & Ill_Status != "")

# ======================================================
# Prepare Metadata
# ======================================================
analysis_data <- data.frame(
  Randomisation_No = merged_data$Randomisation.No,
  Ill_Status = merged_data$Ill_Status,
  AgeGroup = merged_data$X3agegroups.based.on.age.at.sampling,
  Timepoints = merged_data$timepoints,
  SampleID = paste0(merged_data$Randomisation.No, "_", merged_data$timepoints),
  row.names = paste0(merged_data$Randomisation.No, "_", merged_data$timepoints)
)

# Fix AgeGroup labels
analysis_data$AgeGroup <- recode(analysis_data$AgeGroup,
                                 "7to12 mths" = "7 to 12 months",
                                 "1to2years" = "1 to 2 years",
                                 "plus2years" = "plus 2 years")
analysis_data$AgeGroup <- factor(analysis_data$AgeGroup,
                                 levels = c("7 to 12 months", "1 to 2 years", "plus 2 years"))

# ======================================================
# Prepare Microbiome Data
# ======================================================
microbiome_counts <- merged_data[, 13:500]
rownames(microbiome_counts) <- rownames(analysis_data)

# ======================================================
# Create Phyloseq Object
# ======================================================
otu_tab <- otu_table(as.matrix(microbiome_counts), taxa_are_rows = FALSE)
sample_df <- sample_data(analysis_data)
physeq <- phyloseq(otu_tab, sample_df)

# Factor conversion
sample_data(physeq)$Ill_Status <- factor(sample_data(physeq)$Ill_Status, levels = c("Ill", "Not-Ill"))
sample_data(physeq)$AgeGroup <- factor(sample_data(physeq)$AgeGroup,
                                       levels = c("7 to 12 months", "1 to 2 years", "plus 2 years"))
sample_data(physeq)$Timepoints <- factor(sample_data(physeq)$Timepoints)
sample_data(physeq)$Randomisation_No <- factor(sample_data(physeq)$Randomisation_No)

# ======================================================
# Normalize to Relative Abundances
# ======================================================
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# Subset to D1 and D85 only
ps_d1_d85 <- subset_samples(physeq_rel, Timepoints %in% c("D1","D85"))

# Keep only subjects with both D1 and D85
subject_counts <- table(sample_data(ps_d1_d85)$Randomisation_No)
subjects_keep <- names(subject_counts[subject_counts == 2])
ps_d1_d85 <- subset_samples(ps_d1_d85, Randomisation_No %in% subjects_keep)

# ======================================================
# Compute Bray–Curtis distances
# ======================================================
otu_mat <- as(otu_table(ps_d1_d85), "matrix")
if (taxa_are_rows(ps_d1_d85)) otu_mat <- t(otu_mat)
dist_bc <- vegdist(otu_mat, method = "bray")
sim_mat <- 1 - as.matrix(dist_bc)  # similarity = 1 - distance

meta <- as(sample_data(ps_d1_d85), "data.frame")

# ======================================================
# Extract D1 vs D85 similarity per subject
# ======================================================
similarity_df <- lapply(subjects_keep, function(subj) {
  subj_samples <- rownames(meta[meta$Randomisation_No == subj, ])
  if (length(subj_samples) == 2) {
    sim_score <- sim_mat[subj_samples[1], subj_samples[2]]
    data.frame(
      SubjectID = subj,
      Similarity = sim_score,
      AgeGroup = unique(meta[subj_samples, "AgeGroup"]),
      IllStatus = unique(meta[subj_samples, "Ill_Status"])
    )
  } else {
    NULL
  }
}) %>% bind_rows()

# ======================================================
# Summary Statistics
# ======================================================
summary_stats <- similarity_df %>%
  group_by(AgeGroup, IllStatus) %>%
  summarise(N = n(),
            Mean = mean(Similarity, na.rm = TRUE),
            SD = sd(Similarity, na.rm = TRUE),
            Median = median(Similarity, na.rm = TRUE),
            .groups = "drop")

write.csv(summary_stats, "C:/Users/oofordile/Desktop/WithinIndividual_Similarity_Summary.csv", row.names = FALSE)

# ======================================================
# Wilcoxon Tests + FDR correction
# ======================================================
wilcox_results <- similarity_df %>%
  group_by(AgeGroup) %>%
  do({
    test <- tryCatch(
      wilcox.test(Similarity ~ IllStatus, data = ., exact = FALSE),
      error = function(e) NULL
    )
    if (!is.null(test)) {
      data.frame(
        AgeGroup = unique(.$AgeGroup),
        W = test$statistic,
        P_value = test$p.value
      )
    } else {
      data.frame(
        AgeGroup = unique(.$AgeGroup),
        W = NA,
        P_value = NA
      )
    }
  }) %>%
  ungroup() %>%
  mutate(FDR_q = p.adjust(P_value, method = "fdr"))

write.csv(wilcox_results, "C:/Users/oofordile/Desktop/WithinIndividual_Similarity_Wilcoxon.csv", row.names = FALSE)

# ======================================================
# Prepare labels for boxplot (show FDR q-values)
# ======================================================
label_df <- wilcox_results %>%
  dplyr::mutate(label = paste0("q = ", signif(FDR_q, 3))) %>%
  dplyr::select(AgeGroup, label)

# ======================================================
# Boxplot with FDR q-value annotations, thin axis lines, and inward ticks
# ======================================================
p <- ggplot(similarity_df, aes(x = AgeGroup, y = Similarity, fill = IllStatus)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +
  geom_text(data = label_df, aes(x = AgeGroup, y = 1.02, label = label), inherit.aes = FALSE, size = 4) +
  labs(title = "Within-individual Bray–Curtis Similarity (D1 vs D85)",
       y = "Similarity (1 - Bray–Curtis distance)",
       x = "Age Group") +
  scale_fill_manual(values = c("Ill" = "#D73027", "Not-Ill" = "#4575B4")) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(-0.2, "cm"),   # inward ticks
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right"
  )

# ======================================================
# Export: PDF and EMF
# ======================================================
# PDF
pdf("C:/Users/oofordile/Desktop/WithinIndividual_Similarity_Boxplot.pdf", width = 7, height = 5)
print(p)
dev.off()

# EMF
emf("C:/Users/oofordile/Desktop/WithinIndividual_Similarity_Boxplot.emf", width = 7, height = 5)
print(p)
dev.off()

# ======================================================
# Print outputs
# ======================================================
print(summary_stats)
print(wilcox_results)
