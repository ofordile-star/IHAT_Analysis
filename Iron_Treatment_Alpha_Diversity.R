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
library(Cairo)    # for PDF
library(devEMF)   # for EMF output

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

# --- Factor levels ---
sample_data(physeq)$Ill_Status <- factor(sample_data(physeq)$Ill_Status, levels = c("Ill", "Not-Ill"))
sample_data(physeq)$Iron_Treatment <- factor(sample_data(physeq)$Iron_Treatment)
sample_data(physeq)$Age_Group <- factor(sample_data(physeq)$Age_Group)
sample_data(physeq)$Randomisation_No <- factor(sample_data(physeq)$Randomisation_No)

# --- ALPHA DIVERSITY ---
alpha_div <- estimate_richness(physeq, measures = c("Observed", "Shannon", "Simpson", "Fisher"))
colnames(alpha_div)[colnames(alpha_div) == "Observed"] <- "Richness"

# --- Calculate Pielou's evenness ---
pielou <- diversity(otu_table(physeq), index = "shannon") / log(specnumber(otu_table(physeq)))
alpha_div$Pielou <- pielou

# --- Fix SampleID ---
alpha_div$SampleID <- rownames(alpha_div)
alpha_div$SampleID <- sub("^X", "", alpha_div$SampleID)

# --- Join with metadata ---
alpha_div_data <- alpha_div %>% inner_join(analysis_data, by = "SampleID")

# --- Ensure factor types ---
alpha_div_data$Ill_Status <- factor(alpha_div_data$Ill_Status, levels = c("Ill", "Not-Ill"))
alpha_div_data$Iron_Treatment <- factor(alpha_div_data$Iron_Treatment)
alpha_div_data$Age_Group <- factor(alpha_div_data$Age_Group)
alpha_div_data$Randomisation_No <- factor(alpha_div_data$Randomisation_No)

# --- Calculate sample sizes per Ill_Status × Iron_Treatment group ---
group_counts <- alpha_div_data %>%
  group_by(Ill_Status, Iron_Treatment) %>%
  summarise(n = n(), .groups = "drop")

cat("\nSample sizes per Ill_Status × Iron_Treatment group:\n")
print(group_counts)

# --- Mixed Effects Modeling ---
metrics <- c("Richness", "Shannon", "Simpson", "Fisher", "Pielou")
alpha_results <- list()

for (metric in metrics) {
  formula <- as.formula(paste0(metric, " ~ Ill_Status * Iron_Treatment + Age_Group + (1|Randomisation_No)"))
  model <- lmer(formula, data = alpha_div_data)
  tidy_res <- broom.mixed::tidy(model, effects = "fixed")
  alpha_results[[metric]] <- tidy_res
}

all_alpha_results <- bind_rows(
  lapply(names(alpha_results), function(metric) {
    alpha_results[[metric]] %>% filter(effect == "fixed") %>% mutate(Metric = metric)
  })
) %>% group_by(Metric) %>% mutate(p.value.adjusted = p.adjust(p.value, method = "fdr")) %>% ungroup()

# --- Save results ---
csv_path <- "C:/Users/oofordile/Desktop/Alpha_Diversity_MixedModel_Results.csv"
write_csv(all_alpha_results, csv_path)

# --- Extract ONLY the specific interaction term p-value ---
interaction_results <- all_alpha_results %>% filter(term == "Ill_StatusNot-Ill:Iron_Treatmenttreatment")

# --- Plot Alpha Diversity with Annotated Interaction p-value ---
color_map <- c("Ill" = "#EF2929", "Not-Ill" = "#729FCF")
plot_list <- list()

for (metric in metrics) {
  interaction_p <- interaction_results %>% filter(Metric == metric) %>% pull(p.value.adjusted)
  subtitle_text <- if (length(interaction_p) == 1) { paste0("Interaction p (adj) = ", signif(interaction_p, 3)) } else { "" }
  
  p <- ggplot(alpha_div_data, aes(x = interaction(Ill_Status, Iron_Treatment, sep = " : "), y = .data[[metric]], fill = Ill_Status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
    scale_fill_manual(values = color_map) +
    labs(
      title = paste0(metric, " by Ill Status and Iron Treatment"),
      subtitle = subtitle_text,
      x = "Group",
      y = metric
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.line = element_line(color = "black", size = 0.3),
      axis.ticks = element_line(color = "black", size = 0.3),
      axis.ticks.length = unit(-0.1, "cm")
    )
  plot_list[[metric]] <- p
}

# --- Save Outputs ---
# PDF
pdf_path <- "C:/Users/oofordile/Desktop/Alpha_Diversity_Interaction_Boxplots_withAxis.pdf"
CairoPDF(pdf_path, width = 14, height = 8)
grid.arrange(grobs = plot_list, nrow = 2, ncol = 3,
             top = "Alpha Diversity Metrics by Ill Status × Iron Treatment Interaction")
dev.off()

# EMF (vector, PowerPoint/Word compatible)
emf_path <- "C:/Users/oofordile/Desktop/Alpha_Diversity_Interaction_Boxplots_withAxis.emf"
devEMF::emf(file = emf_path, width = 14, height = 8, family = "Arial")
grid.arrange(grobs = plot_list, nrow = 2, ncol = 3,
             top = "Alpha Diversity Metrics by Ill Status × Iron Treatment Interaction")
dev.off()

cat("\nAnalysis complete.\nCSV saved to:", csv_path,
    "\nPDF saved to:", pdf_path,
    "\nEMF saved to:", emf_path, "\n")
