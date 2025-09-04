# ======================================================
# Load libraries
# ======================================================
library(tidyverse)
library(broom)
library(ggpubr)
library(patchwork)
if(!requireNamespace("devEMF", quietly = TRUE)) install.packages("devEMF")

# ======================================================
# Set file path
# ======================================================
file_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"

# ======================================================
# Read and clean data
# ======================================================
data <- read_csv(file_path) %>%
  filter(!is.na(Ill_Status) & Ill_Status != "")

# ======================================================
# Prepare analysis data
# ======================================================
analysis_data <- data %>%
  dplyr::select(`Randomisation No`, timepoints, Ill_Status, `3agegroups based on age at sampling`, 13:22) %>%
  dplyr::rename(age_group = `3agegroups based on age at sampling`) %>%
  dplyr::mutate(
    timepoints = factor(timepoints, levels = c("D1", "D15", "D85")),
    Ill_Status = factor(Ill_Status, levels = c("Not-Ill", "Ill")),
    age_group = factor(age_group, levels = c("7to12 mths", "1to2years", "plus2years"))
  )

# ======================================================
# Normalize taxa to relative abundance
# ======================================================
taxa_names <- colnames(analysis_data)[5:14]

analysis_data <- analysis_data %>%
  rowwise() %>%
  mutate(across(all_of(taxa_names), ~ .x / sum(c_across(all_of(taxa_names)), na.rm = TRUE))) %>%
  ungroup()

# ======================================================
# Run ANOVA and get adjusted p-values
# ======================================================
anova_results_list <- list()
for(taxon in taxa_names) {
  formula <- as.formula(paste(taxon, "~ Ill_Status * timepoints + age_group"))
  model <- lm(formula, data = analysis_data)
  anova_res <- anova(model)
  anova_df <- as.data.frame(anova_res) %>%
    rownames_to_column("term") %>%
    select(term, Df, `F value`, `Pr(>F)`) %>%
    rename(df = Df, F_value = `F value`, p.value = `Pr(>F)`) %>%
    mutate(taxon = taxon)
  anova_results_list[[taxon]] <- anova_df
}

anova_results <- bind_rows(anova_results_list)

terms_of_interest <- c("Ill_Status", "timepoints", "Ill_Status:timepoints", "age_group")

anova_results_filtered <- anova_results %>%
  filter(term %in% terms_of_interest) %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),  
    significant = p.adj < 0.05,
    term_short = case_when(
      term == "Ill_Status" ~ "p_ill",
      term == "timepoints" ~ "p_time",
      term == "Ill_Status:timepoints" ~ "p_int",
      TRUE ~ term
    )
  )

# ======================================================
# Prepare combined CSV with requested headings
# ======================================================
p_ill <- anova_results_filtered %>% filter(term_short == "p_ill") %>% arrange(match(taxon, taxa_names))
p_time <- anova_results_filtered %>% filter(term_short == "p_time") %>% arrange(match(taxon, taxa_names))
p_int <- anova_results_filtered %>% filter(term_short == "p_int") %>% arrange(match(taxon, taxa_names))

combined_results <- tibble(
  Taxon = taxa_names,
  `Time Effect: F(df), adj. p` = paste0(round(p_time$F_value,3),"(",p_time$df,"), ", signif(p_time$p.adj,4)),
  `Group Effect: F(df), adj. p` = paste0(round(p_ill$F_value,3),"(",p_ill$df,"), ", signif(p_ill$p.adj,4)),
  `Interaction Effect: F(df), adj. p` = paste0(round(p_int$F_value,3),"(",p_int$df,"), ", signif(p_int$p.adj,4))
)

write_csv(combined_results, "C:/Users/oofordile/Desktop/ANOVA_Combined_Results.csv")
cat("Combined ANOVA CSV saved to Desktop as 'ANOVA_Combined_Results.csv'\n")

# ======================================================
# Prepare p-values for plotting
# ======================================================
plot_taxa <- c("Prevotella_stercorea_7.05", "Escherichia_coli_3.51")
pvals_for_plot <- anova_results_filtered %>%
  filter(taxon %in% plot_taxa & term_short == "p_ill") %>%
  select(taxon, p.adj) %>%
  mutate(p.label = paste0("adj p = ", signif(p.adj, 3)))

# ======================================================
# Prepare long-format data for plotting
# ======================================================
long_data <- analysis_data %>%
  dplyr::select(`Randomisation No`, timepoints, Ill_Status, all_of(plot_taxa)) %>%
  tidyr::pivot_longer(cols = all_of(plot_taxa), names_to = "Taxa", values_to = "Abundance") %>%
  dplyr::mutate(
    Value = ifelse(Taxa == "Escherichia_coli_3.51", log10(Abundance + 1e-6), Abundance),
    TaxaLabel = recode(Taxa, 
                       "Prevotella_stercorea_7.05" = "italic('Prevotella stercorea')",
                       "Escherichia_coli_3.51" = "italic('Escherichia coli')"),
    Ylabel = ifelse(Taxa == "Escherichia_coli_3.51", "Log Abundance", "Count")
  )

long_data$Taxa <- factor(long_data$Taxa, levels = c("Prevotella_stercorea_7.05", "Escherichia_coli_3.51"))
long_data$TaxaLabel <- factor(long_data$TaxaLabel, levels = c("italic('Prevotella stercorea')", "italic('Escherichia coli')"))

# ======================================================
# Function to create taxa plot with automatic p-value
# ======================================================
make_taxa_plot <- function(df, ylab) {
  y_max <- max(df$Value, na.rm = TRUE)
  y_min <- min(df$Value, na.rm = TRUE)
  y_range <- y_max - y_min
  taxa_name <- unique(df$Taxa)
  pval <- pvals_for_plot$p.label[pvals_for_plot$taxon == taxa_name]
  
  ggplot(df, aes(x = timepoints, y = Value, fill = Ill_Status)) +
    geom_boxplot(position = position_dodge(0.7), width = 0.6, alpha = 0.8, color = "black", outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), alpha = 0.4, size = 1) +
    stat_summary(fun = median, geom = "line", aes(group = Ill_Status, color = Ill_Status), size = 1, position = position_dodge(0.7)) +
    geom_text(data = df %>% dplyr::distinct(TaxaLabel), aes(x = 2, y = y_max + y_range*0.05, label = pval), inherit.aes = FALSE, size = 5) +
    scale_fill_manual(values = c("Ill" = "#FF6666", "Not-Ill" = "#6699FF")) +
    scale_color_manual(values = c("Ill" = "red", "Not-Ill" = "blue")) +
    labs(x = "Study Day", y = ylab, fill = "Status", color = "Status") +
    facet_wrap(~TaxaLabel, labeller = label_parsed) +
    theme_pubr(border = TRUE, base_size = 14) +
    theme(
      strip.background = element_rect(fill = "gray95", color = "black", size = 0.5),
      strip.text = element_text(size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.spacing = unit(0.5, "lines"),
      axis.ticks = element_line(size = 0.3, color = "black"),
      axis.ticks.length = unit(-2, "mm")
    ) +
    expand_limits(y = y_max + y_range*0.1)
}

# ======================================================
# Create plots
# ======================================================
prev_plot <- long_data %>% filter(Taxa == "Prevotella_stercorea_7.05") %>% make_taxa_plot(ylab = "Count")
ecoli_plot <- long_data %>% filter(Taxa == "Escherichia_coli_3.51") %>% make_taxa_plot(ylab = "Log Abundance")

# Combine plots with shared legend
taxa_plot <- prev_plot + ecoli_plot + plot_layout(ncol = 2, guides = "collect") & 
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.margin = margin(t = 10, b = 10)
  )

# Add main title
taxa_plot <- taxa_plot + plot_annotation(
  title = "Taxa with significant Ill vs Not-Ill gap (FDR and Age adjusted p < 0.05)",
  theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
)

# ======================================================
# Save figure
# ======================================================
ggsave("C:/Users/oofordile/Desktop/Ill_vs_NotIll_Taxa_Boxplots.pdf", taxa_plot, width = 10, height = 6, device = "pdf", dpi = 300)
devEMF::emf("C:/Users/oofordile/Desktop/Ill_vs_NotIll_Taxa_Boxplots.emf", width = 10, height = 6)
print(taxa_plot)
dev.off()
