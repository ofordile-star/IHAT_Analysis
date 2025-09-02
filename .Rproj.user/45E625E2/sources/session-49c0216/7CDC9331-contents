# ======================================================
# Load required libraries
# ======================================================
library(tidyverse)
library(reshape2)
library(Cairo)
library(devEMF)  # For EMF export
library(grid)     # For unit()

# ======================================================
# Set file path
# ======================================================
file_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"

# ======================================================
# Load data
# ======================================================
data <- read.csv(file_path, check.names = FALSE)

# Remove rows with blank Ill_Status
data <- data %>% filter(!is.na(Ill_Status) & Ill_Status != "")

# Extract metadata and microbiome data
microbiome_data <- data[, 13:500]
meta_data <- data[, c("Randomisation No", "timepoints", "Ill_Status")]

# Combine metadata and microbiome
full_data <- cbind(meta_data, microbiome_data)

# Melt to long format
long_data <- melt(full_data, id.vars = c("Randomisation No", "timepoints", "Ill_Status"),
                  variable.name = "Taxon", value.name = "Abundance")

# Clean Taxon names: remove numeric suffixes and underscores
long_data$Taxon <- gsub("_\\d+\\.\\d+$", "", long_data$Taxon)
long_data$Taxon <- gsub("_", " ", long_data$Taxon)

# ======================================================
# Calculate top 50 most abundant taxa overall
# ======================================================
top_taxa_ordered <- long_data %>%
  group_by(Taxon) %>%
  summarise(Total_Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total_Abundance)) %>%
  slice(1:50)

top_taxa <- top_taxa_ordered$Taxon

# Filter to top taxa only
long_data_top <- long_data %>% filter(Taxon %in% top_taxa)

# Convert Taxon to factor (ordered by abundance)
long_data_top$Taxon <- factor(long_data_top$Taxon, levels = top_taxa)

# ======================================================
# Calculate mean and SEM
# ======================================================
summary_data <- long_data_top %>%
  group_by(Taxon, timepoints, Ill_Status) %>%
  summarise(
    Mean = mean(Abundance, na.rm = TRUE),
    SEM = sd(Abundance, na.rm = TRUE)/sqrt(n()),
    .groups = "drop"
  )

# Taxon labels for italicization (as expressions)
summary_data$Taxon <- factor(summary_data$Taxon, levels = top_taxa)
summary_data$TaxonLabel <- paste0("italic('", summary_data$Taxon, "')")
summary_data$TaxonLabel <- factor(summary_data$TaxonLabel,
                                  levels = paste0("italic('", top_taxa, "')"))

# ======================================================
# Export CSV of summary data
# ======================================================
write.csv(summary_data, 
          "C:/Users/oofordile/Desktop/Top50_Taxa_Longitudinal_Summary.csv",
          row.names = FALSE)

# ======================================================
# Create plot
# ======================================================
plot <- ggplot(summary_data, aes(x = timepoints, y = Mean, group = Ill_Status, color = Ill_Status)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2) +
  facet_wrap(~ TaxonLabel, scales = "free_y", ncol = 5, labeller = label_parsed) +
  scale_color_manual(values = c("Ill" = "#FF0000", "Not-Ill" = "#0000FF")) +
  labs(x = "Timepoint", y = "Mean Abundance ± SEM", color = "Ill Status") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_line(size = 0.3, color = "black"),
    axis.ticks.length = unit(-2, "mm"),  # Negative value makes ticks point inward
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

# ======================================================
# Save PDF
# ======================================================
ggsave("C:/Users/oofordile/Desktop/Top50_Taxa_Longitudinal_Sorted.pdf",
       plot = plot, device = cairo_pdf, width = 12, height = 15, units = "in")

# ======================================================
# Save EMF
# ======================================================
emf("C:/Users/oofordile/Desktop/Top50_Taxa_Longitudinal_Sorted.emf",
    width = 12, height = 15)
print(plot)
dev.off()
