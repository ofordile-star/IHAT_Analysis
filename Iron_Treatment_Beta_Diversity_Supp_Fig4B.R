# --- Load Libraries ---
library(phyloseq)
library(vegan)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(readr)
library(ggpubr)
library(ape)
library(devEMF)   # for EMF export

# --- Load Data ---
file_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
merged_data <- read.csv(file_path)

# --- Remove rows with blank or NA Ill_Status ---
df <- merged_data %>% filter(!is.na(Ill_Status) & Ill_Status != "")

# --- Prepare Metadata for Analysis ---
analysis_data <- data.frame(
  Randomisation_No = df$Randomisation.No,
  Ill_Status = df$Ill_Status,
  Iron_Treatment = df$Iron,
  Age_Group = df$X3agegroups.based.on.age.at.sampling,
  Timepoints = df$timepoints,
  SampleID = paste0(df$Randomisation.No, "_", df$timepoints),
  row.names = paste0(df$Randomisation.No, "_", df$timepoints)
)

# --- Prepare Microbiome Counts ---
microbiome_counts <- df[, 13:500]
rownames(microbiome_counts) <- paste0(df$Randomisation.No, "_", df$timepoints)
stopifnot(all(rownames(microbiome_counts) == rownames(analysis_data)))

# --- Create phyloseq object ---
otu_mat <- as.matrix(microbiome_counts)
otu_tab <- otu_table(otu_mat, taxa_are_rows = FALSE)
sample_df <- sample_data(analysis_data)
physeq <- phyloseq(otu_tab, sample_df)

# --- Normalize to Relative Abundance ---
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# --- Clean factor levels ---
sample_data(physeq_rel)$Ill_Status <- factor(trimws(sample_data(physeq_rel)$Ill_Status))
sample_data(physeq_rel)$Iron_Treatment <- factor(trimws(sample_data(physeq_rel)$Iron_Treatment))
sample_data(physeq_rel)$Age_Group <- factor(trimws(sample_data(physeq_rel)$Age_Group))
sample_data(physeq_rel)$Timepoints <- factor(trimws(sample_data(physeq_rel)$Timepoints))
sample_data(physeq_rel)$Randomisation_No <- factor(sample_data(physeq_rel)$Randomisation_No)

# --- Calculate Bray–Curtis distance ---
otu_mat_rel <- otu_table(physeq_rel)
if (taxa_are_rows(physeq_rel)) otu_mat_rel <- t(otu_mat_rel)
dist_matrix <- vegdist(otu_mat_rel, method = "bray")

# --- PCoA ---
pcoa_res <- ape::pcoa(dist_matrix)
pcoa_df <- as.data.frame(pcoa_res$vectors[, 1:2])
colnames(pcoa_df) <- c("Axis.1", "Axis.2")
pcoa_df$SampleID <- rownames(pcoa_df)

# Join metadata
pcoa_df <- pcoa_df %>% inner_join(analysis_data, by = "SampleID")
pcoa_df$Ill_Status <- factor(trimws(pcoa_df$Ill_Status))
pcoa_df$Iron_Treatment <- factor(trimws(pcoa_df$Iron_Treatment))
pcoa_df$Age_Group <- factor(trimws(pcoa_df$Age_Group))

# --- PERMANOVA ---
adonis_res_terms <- adonis2(
  dist_matrix ~ Ill_Status * Iron_Treatment * Timepoints + Age_Group,
  data = as(sample_data(physeq_rel), "data.frame"),
  permutations = 999,
  by = "terms"
)

adonis_df_terms <- as.data.frame(adonis_res_terms)
adonis_df_terms$Term <- rownames(adonis_df_terms)
colnames(adonis_df_terms) <- gsub(" ", "", colnames(adonis_df_terms))
colnames(adonis_df_terms) <- sub("Pr\\(>F\\)", "p_value", colnames(adonis_df_terms))

cat("\n--- PERMANOVA Results (all terms) ---\n")
print(adonis_df_terms, digits = 4)

# --- Export PERMANOVA Results to CSV ---
csv_path <- "C:/Users/oofordile/Desktop/PCoA_BrayCurtis_PERMANOVA_Results.csv"
write.csv(adonis_df_terms, csv_path, row.names = FALSE)

# --- Extract p-value for interaction ---
p_int <- adonis_df_terms$p_value[adonis_df_terms$Term == "Ill_Status:Iron_Treatment"]
p_int_label <- paste0("Interaction p = ", signif(p_int, 3))

# --- Aesthetics ---
color_map <- c("Ill" = "#D73027", "Not-Ill" = "#1F78B4")  
shape_map <- c("1to2years" = 21, "7to12 mths" = 22, "plus2years" = 24)

pcoa_placebo <- pcoa_df %>% filter(Iron_Treatment == "placebo")
pcoa_treatment <- pcoa_df %>% filter(Iron_Treatment == "treatment")

# --- Plot function with thin axis lines and inward ticks ---
make_pcoa_plot <- function(data, group_title, pval_text) {
  x_pos <- mean(range(data$Axis.1, na.rm = TRUE))
  y_pos <- max(data$Axis.2, na.rm = TRUE) * 1.05
  
  ggplot(data, aes(x = Axis.1, y = Axis.2, color = Ill_Status, shape = Age_Group)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = shape_map) +
    labs(
      title = group_title,
      x = paste0("PCoA1 (", round(pcoa_res$values$Relative_eig[1] * 100, 1), "%)"),
      y = paste0("PCoA2 (", round(pcoa_res$values$Relative_eig[2] * 100, 1), "%)")
    ) +
    annotate(
      "text",
      x = x_pos,
      y = y_pos,
      label = pval_text,
      hjust = 0.5,
      vjust = 0,
      size = 5
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.2),      # thin axis lines
      axis.ticks = element_line(color = "black", size = 0.2),     # thin ticks
      axis.ticks.length = unit(-0.1, "cm")                        # inward ticks
    )
}

# --- Create plots ---
plot_placebo <- make_pcoa_plot(pcoa_placebo, "Placebo", p_int_label)
plot_treatment <- make_pcoa_plot(pcoa_treatment, "Treatment", p_int_label)

# --- Export side-by-side PDF without main title ---
pdf_path <- "C:/Users/oofordile/Desktop/PCoA_BrayCurtis_by_Treatment_and_IllStatus_thinAxis.pdf"
pdf(pdf_path, width = 12, height = 6)
grid.arrange(plot_placebo, plot_treatment, nrow = 1)  # Removed `top` argument
dev.off()

# --- Export EMF without main title ---
emf_path <- "C:/Users/oofordile/Desktop/PCoA_BrayCurtis_by_Treatment_and_IllStatus_thinAxis.emf"
emf(file = emf_path, width = 12, height = 6, family = "Arial")
grid.arrange(plot_placebo, plot_treatment, nrow = 1)  # Removed `top` argument
dev.off()

# --- Confirmation ---
cat("\nAnalysis complete.\nCSV saved to:", csv_path,
    "\nPDF saved to:", pdf_path,
    "\nEMF saved to:", emf_path, "\n")
