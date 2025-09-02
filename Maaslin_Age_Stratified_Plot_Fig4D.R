# ======================================================
# Load required libraries
# ======================================================
library(ggplot2)
library(dplyr)
library(readr)
library(forcats)
library(patchwork)
library(devEMF)

# ======================================================
# Load filtered MaAsLin2 summary (q < 0.1)
# ======================================================
results <- read_tsv("C:/Users/oofordile/Desktop/maaslin2_output_stratified/Stratified_Summary_qval0.1.tsv")

# ======================================================
# DEBUG: Check unique values in your data
# ======================================================
cat("Unique AgeGroup values:\n")
print(unique(results$AgeGroup))
cat("\nUnique Comparison values:\n")
print(unique(results$Comparison))

cat("\nChecking 1to2years + D15 Ill_vs_D85 Not-Ill combination:\n")
debug_subset <- results %>% 
  filter(AgeGroup == "1to2years", Comparison == "D15 Ill_vs_D85 Not-Ill")
cat("Number of rows found:", nrow(debug_subset), "\n")

cat("\nAll AgeGroup + Comparison combinations:\n")
print(results %>% count(AgeGroup, Comparison))

# ======================================================
# Function to clean taxa names
# ======================================================
clean_taxa_names <- function(taxa_name) {
  cleaned <- gsub("_\\d+\\.\\d+$", "", taxa_name)  # remove suffixes
  cleaned <- gsub("_", " ", cleaned)               # underscores → spaces
  return(cleaned)
}

# ======================================================
# Clean and prepare data
# ======================================================
results <- results %>%
  mutate(
    AgeGroup = trimws(AgeGroup),
    Comparison = trimws(Comparison),
    Taxon = clean_taxa_names(feature),
    FillColor = ifelse(coef > 0, "#D73027", "#4575B4"),
    TaxonLabel = paste0("italic('", Taxon, "')")
  ) %>%
  arrange(feature)

# ======================================================
# Plot function
# ======================================================
plot_subset <- function(age_group, comp, custom_title = NULL) {
  age_group <- trimws(age_group)
  comp <- trimws(comp)
  
  cat(paste("Attempting AgeGroup:", age_group, "Comparison:", comp, "\n"))
  
  subset_data <- results %>% 
    filter(AgeGroup == age_group, Comparison == comp)
  
  cat("Found", nrow(subset_data), "rows\n")
  if(nrow(subset_data) == 0) {
    cat("⚠️ No significant taxa for", age_group, comp, "\n")
    return(NULL)
  }
  
  if(is.null(custom_title)) custom_title <- paste(age_group, "-", comp)
  
  subset_data <- subset_data %>%
    mutate(
      Taxon = fct_reorder(Taxon, coef),
      FillColor = ifelse(coef > 0, "#D73027", "#4575B4")
    )
  
  ggplot(subset_data, aes(x = coef, y = Taxon)) +
    geom_col(aes(fill = FillColor), width = 0.6) +
    geom_text(aes(label = paste0("q=", round(qval, 3))),
              hjust = ifelse(subset_data$coef > 0, 1.1, -0.1),
              color = "black", size = 8) +             
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    scale_fill_identity() +
    scale_y_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +
    labs(
      x = "Effect Size (Coef)",
      y = "Taxa",
      title = custom_title
    ) +
    theme_bw(base_size = 16) +
    theme(
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, size = 0.3),
      axis.ticks = element_line(size = 0.3, color = "black"),
      axis.ticks.length = unit(-2, "mm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5, 20, 5, 5)  # << tiny extra right margin
    )
}

# ======================================================
# Generate plots
# ======================================================
plot_7to12_D15 <- plot_subset("7to12 mths", "D15 Ill_vs_D85 Not-Ill",
                              "7–12 months: D15 Ill vs D85 Not-Ill")

plot_1to2_D1 <- plot_subset("1to2years", "D1 Ill_vs_D85 Not-Ill",
                            "1–2 years: D1 Ill vs D85 Not-Ill")

plot_1to2_D15 <- plot_subset("1to2years", "D15 Ill_vs_D85 Not-Ill",
                             "1–2 years: D15 Ill vs D85 Not-Ill")

# fallback if needed
if(is.null(plot_1to2_D15)) {
  alt_age_groups <- c("1-2years","1to2 years","1-2 years","1to2yrs","1-2yrs","1 to 2 years")
  for(alt_age in alt_age_groups) {
    if(alt_age %in% results$AgeGroup) {
      plot_1to2_D15 <- plot_subset(alt_age, "D15 Ill_vs_D85 Not-Ill",
                                   "1–2 years: D15 Ill vs D85 Not-Ill")
      if(!is.null(plot_1to2_D15)) break
    }
  }
}

# ======================================================
# Combine plots
# ======================================================
valid_plots <- list(plot_7to12_D15, plot_1to2_D1, plot_1to2_D15)
valid_plots <- valid_plots[!sapply(valid_plots, is.null)]

if(length(valid_plots) == 1) {
  combined_plot <- valid_plots[[1]]
} else if(length(valid_plots) == 2) {
  combined_plot <- valid_plots[[1]] | valid_plots[[2]]
} else {
  combined_plot <- valid_plots[[1]] | valid_plots[[2]] | valid_plots[[3]]
}

# Add global margin to the combined plot too
combined_plot <- combined_plot & theme(plot.margin = margin(5, 20, 5, 5))

# ======================================================
# Export
# ======================================================
ggsave("C:/Users/oofordile/Desktop/Maaslin2_Bars_All_AgeGroups.pdf",
       combined_plot, width = 18, height = 6, dpi = 600)

emf("C:/Users/oofordile/Desktop/Maaslin2_Bars_All_AgeGroups.emf",
    width = 18, height = 6)
print(combined_plot)
dev.off()

cat("Plots exported successfully\n")
