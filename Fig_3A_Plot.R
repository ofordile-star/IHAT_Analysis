# ======================================================
# Load libraries
# ======================================================
library(tidyverse)
library(ggpubr)
library(patchwork)

# ======================================================
# Set file path
# ======================================================
file_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"

# ======================================================
# Read data
# ======================================================
data <- read_csv(file_path) %>%
  filter(!is.na(Ill_Status) & Ill_Status != "")

# ======================================================
# Prepare data (use dplyr::select explicitly)
# ======================================================
analysis_data <- data %>%
  dplyr::select(`Randomisation No`, timepoints, Ill_Status, `3agegroups based on age at sampling`, 13:22) %>%
  dplyr::rename(age_group = `3agegroups based on age at sampling`) %>%
  dplyr::mutate(
    timepoints = factor(timepoints, levels = c("D1", "D15", "D85")),
    Ill_Status = factor(Ill_Status, levels = c("Not-Ill", "Ill"))
  )

# ======================================================
# Target taxa
# ======================================================
ecoli_taxon <- "Escherichia_coli_3.51"
prev_taxon <- "Prevotella_stercorea_7.05"

# ======================================================
# Prepare long format for faceting
# ======================================================
long_data <- analysis_data %>%
  dplyr::select(`Randomisation No`, timepoints, Ill_Status, 
                ecoli = all_of(ecoli_taxon), prev = all_of(prev_taxon)) %>%
  tidyr::pivot_longer(cols = c(ecoli, prev), names_to = "Taxa", values_to = "Abundance") %>%
  dplyr::mutate(
    Value = ifelse(Taxa == "ecoli", 
                   pmax(0, log10(Abundance + 1e-6)), # log abundance for E. coli
                   Abundance), # raw counts for Prevotella
    TaxaLabel = recode(Taxa, 
                       "ecoli" = "italic('Escherichia coli')", 
                       "prev" = "italic('Prevotella stercorea')"), # italic via parse
    Ylabel = ifelse(Taxa == "ecoli", "Log Abundance", "Count")
  )

# Lock order: Prevotella first, E. coli second
long_data$Taxa <- factor(long_data$Taxa, levels = c("prev", "ecoli"))
long_data$TaxaLabel <- factor(long_data$TaxaLabel, levels = c("italic('Prevotella stercorea')", "italic('Escherichia coli')"))

# ======================================================
# Function to plot each taxa separately with p-value under facet label
# ======================================================
make_taxa_plot <- function(df, ylab, pval) {
  
  y_max <- max(df$Value, na.rm = TRUE)
  y_min <- min(df$Value, na.rm = TRUE)
  y_range <- y_max - y_min
  
  ggplot(df, aes(x = timepoints, y = Value, fill = Ill_Status)) +
    geom_boxplot(position = position_dodge(0.7), width = 0.6, alpha = 0.8, color = "black", outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), alpha = 0.4, size = 1) +
    stat_summary(fun = median, geom = "line", aes(group = Ill_Status, color = Ill_Status), 
                 size = 1, position = position_dodge(0.7)) +
    # Add p-value centered below facet label
    geom_text(
      data = df %>% dplyr::distinct(TaxaLabel),
      aes(x = 2, y = y_max + y_range*0.05, label = pval),
      inherit.aes = FALSE, size = 5, fontface = "plain"
    ) +
    scale_fill_manual(values = c("Ill" = "#FF6666", "Not-Ill" = "#6699FF")) +
    scale_color_manual(values = c("Ill" = "red", "Not-Ill" = "blue")) +
    labs(x = "Study Day", y = ylab, fill = "Status", color = "Status") +
    facet_wrap(~TaxaLabel, labeller = label_parsed) + # italic via parse
    theme_pubr(border = TRUE, base_size = 14) +
    theme(
      strip.background = element_rect(fill = "gray95", color = "black", size = 0.5),
      strip.text = element_text(size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.spacing = unit(0.5, "lines"),
      # Add inward-pointing ticks
      axis.ticks = element_line(size = 0.3, color = "black"),
      axis.ticks.length = unit(-2, "mm")  # Negative value makes ticks point inward
    ) +
    expand_limits(y = y_max + y_range*0.1)
}

# ======================================================
# Build plots separately
# ======================================================
prev_plot <- long_data %>%
  dplyr::filter(Taxa == "prev") %>%
  make_taxa_plot(ylab = "Count", pval = "adj p = 0.0228")

ecoli_plot <- long_data %>%
  dplyr::filter(Taxa == "ecoli") %>%
  make_taxa_plot(ylab = "Log Abundance", pval = "adj p = 0.0249")

# ======================================================
# Combine plots with shared legend (legend on bottom)
# ======================================================
taxa_plot <- prev_plot + ecoli_plot + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.margin = margin(t = 10, b = 10)
  )

# ======================================================
# Main title
# ======================================================
taxa_plot <- taxa_plot + 
  plot_annotation(
    title = "Taxa with significant Ill vs Not Ill gap (FDR and Age adjusted p < 0.05)",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
  )

# ======================================================
# Save figure (PDF and EMF)
# ======================================================
# High-res PDF
ggsave("C:/Users/oofordile/Desktop/Ill_vs_NotIll_Taxa_Boxplots.pdf", 
       taxa_plot, width = 10, height = 6, device = "pdf", dpi = 300)

# EMF (Windows vector format)
if(!requireNamespace("devEMF", quietly = TRUE)) install.packages("devEMF")
devEMF::emf("C:/Users/oofordile/Desktop/Ill_vs_NotIll_Taxa_Boxplots.emf", width = 10, height = 6)
print(taxa_plot)
dev.off()
