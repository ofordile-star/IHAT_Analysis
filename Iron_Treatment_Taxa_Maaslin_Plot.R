# --- Load libraries ---
library(ggplot2)
library(dplyr)
library(forcats)
library(cowplot)
library(readr)
library(Cairo)      # for high-res PNG
library(devEMF)     # for EMF export

# --- Load MaAsLin2 results ---
maaslin_output <- "C:/Users/oofordile/Desktop/maaslin2_iron_illness_top50_output/significant_results.tsv"
sig_iron <- read_tsv(maaslin_output)

# --- Filter for Iron effect only ---
sig_iron <- sig_iron %>%
  filter(metadata == "Iron" & qval < 0.05) %>%
  arrange(coef)

# --- Clean feature names ---
sig_iron <- sig_iron %>%
  mutate(
    feature_clean = feature %>%
      gsub("_", " ", .) %>%
      gsub("\\s*\\d+(\\.\\d+)?$", "", .) %>%
      trimws(),
    feature_clean = fct_reorder(feature_clean, coef),
    color = ifelse(coef > 0, "#4CAF50", "#A0522D")
  )

# --- Base bar plot ---
p <- ggplot(sig_iron, aes(x = coef, y = feature_clean, fill = color)) +
  geom_col(width = 0.7) +
  scale_fill_identity() +
  geom_vline(xintercept = 0, color = "gray40", linetype = "dashed") +
  geom_text(
    aes(label = paste0("q  =  ", signif(qval, 2))),
    hjust = ifelse(sig_iron$coef > 0, 1.05, -0.05),
    color = "black",
    size = 5.5   # increased font size
  ) +
  labs(x = "Effect size (Coefficient)", y = NULL) +
  theme_minimal(base_size = 16) +   # bigger overall fonts
  theme(
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(-0.2, "cm"),  # inward ticks
    plot.margin = margin(15, 40, 15, 40),
    axis.text.y = element_text(face = "italic", color = "black", size = 14),
    axis.text.x = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 16, face = "bold")
  )

# --- Add title manually with cowplot ---
title_text <- "Significant Taxa Differentially Abundant by Iron Treatment"

final_plot <- ggdraw() +
  draw_label(title_text, x = 0.5, y = 0.98, hjust = 0.5, vjust = 1,
             size = 20, fontface = "bold", color = "black") +
  draw_plot(p, y = 0, height = 0.95)

# --- Export PDF ---
pdf("C:/Users/oofordile/Desktop/sig_iron_taxa_plot_inwardsTicks.pdf", width = 8, height = 5)
print(final_plot)
dev.off()

# --- Export EMF ---
emf("C:/Users/oofordile/Desktop/sig_iron_taxa_plot_inwardsTicks.emf", width = 8, height = 5,
    emfPlus = FALSE, pointsize = 12)
print(final_plot)
dev.off()
