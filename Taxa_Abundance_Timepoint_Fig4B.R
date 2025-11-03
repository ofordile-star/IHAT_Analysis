# ======================================================
# Load required libraries
# ======================================================
library(tidyverse)
library(reshape2)
library(Cairo)
library(devEMF)
library(grid)

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

# ======================================================
# Convert to relative abundance per sample
# ======================================================
microbiome_rel <- microbiome_data / rowSums(microbiome_data, na.rm = TRUE) * 100  # Convert to %
full_data <- cbind(meta_data, microbiome_rel)

# ======================================================
# Melt to long format
# ======================================================
long_data <- melt(full_data,
                  id.vars = c("Randomisation No", "timepoints", "Ill_Status"),
                  variable.name = "Taxon",
                  value.name = "RelAbundance")

# Clean Taxon names: remove numeric suffixes and underscores
long_data$Taxon <- gsub("_\\d+\\.\\d+$", "", long_data$Taxon)
long_data$Taxon <- gsub("_", " ", long_data$Taxon)

# ======================================================
# Calculate top 50 most relatively abundant taxa overall
# ======================================================
top_taxa_ordered <- long_data %>%
  group_by(Taxon) %>%
  summarise(Total_RelAbundance = sum(RelAbundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total_RelAbundance)) %>%
  slice(1:50)

top_taxa <- top_taxa_ordered$Taxon

# Filter to top taxa only
long_data_top <- long_data %>% filter(Taxon %in% top_taxa)

# Convert Taxon to factor (ordered by abundance)
long_data_top$Taxon <- factor(long_data_top$Taxon, levels = top_taxa)

# ======================================================
# Calculate mean and SEM (relative abundance)
# ======================================================
summary_data <- long_data_top %>%
  group_by(Taxon, timepoints, Ill_Status) %>%
  summarise(
    Mean = mean(RelAbundance, na.rm = TRUE),
    SEM = sd(RelAbundance, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Taxon labels for italicization (as expressions)
summary_data$Taxon <- factor(summary_data$Taxon, levels = top_taxa)
summary_data$TaxonLabel <- paste0("italic('", summary_data$Taxon, "')")
summary_data$TaxonLabel <- factor(summary_data$TaxonLabel,
                                  levels = paste0("italic('", top_taxa, "')"))

# ======================================================
# Export summary data as CSV
# ======================================================
write.csv(summary_data,
          "C:/Users/oofordile/Desktop/Top50_Taxa_RelativeAbundance_Summary.csv",
          row.names = FALSE)

# ======================================================
# Create plot (relative abundance)
# ======================================================
plot <- ggplot(summary_data, aes(x = timepoints, y = Mean, group = Ill_Status, color = Ill_Status)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2) +
  facet_wrap(~ TaxonLabel, scales = "free_y", ncol = 5, labeller = label_parsed) +
  scale_color_manual(values = c("Ill" = "#FF0000", "Not-Ill" = "#0000FF")) +
  labs(
    x = "Timepoint",
    y = "Mean Relative Abundance (%) ± SEM",
    color = "Ill Status"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_line(size = 0.3, color = "black"),
    axis.ticks.length = unit(-2, "mm"),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

# ======================================================
# Save as PDF
# ======================================================
ggsave("C:/Users/oofordile/Desktop/Top50_Taxa_RelativeAbundance.pdf",
       plot = plot, device = cairo_pdf, width = 12, height = 15, units = "in")

# ======================================================
# Save as EMF
# ======================================================
emf("C:/Users/oofordile/Desktop/Top50_Taxa_RelativeAbundance.emf",
    width = 12, height = 15)
print(plot)
dev.off()
