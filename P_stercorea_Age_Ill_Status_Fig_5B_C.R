# --- Load Libraries --- 
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(vcd)
library(gridExtra)
library(forcats)
library(cowplot)
library(devEMF)  # For EMF export

# --- Load Data --- 
data_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
df <- read_csv(data_path)

# --- Define correct age group levels --- 
age_levels <- c("7to12 mths", "1to2years", "plus2years")

# --- Assign factor with exact levels --- 
df$age_group <- factor(df$`3agegroups based on age at sampling`, levels = age_levels)

# --- Define fill colors --- 
age_colors <- c(
  "7to12 mths" = "#A6D854",   # light green
  "1to2years"  = "#FDB863",   # light orange
  "plus2years" = "#80B1D3"    # light blue
)

# ================================
# PLOT 1: Prevotella stercorea Boxplot
# ================================
df_d <- df %>%
  filter(!is.na(`Prevotella_stercorea_7.05`)) %>%
  mutate(logP = log(`Prevotella_stercorea_7.05` + 1))

df_d$age_group <- factor(df_d$age_group, levels = age_levels)

anova_res <- aov(logP ~ age_group, data = df_d)
anova_summary <- summary(anova_res)[[1]]
anova_F <- anova_summary[["F value"]][1]

anova_label <- paste0("ANOVA F(2) = ", round(anova_F, 2), ", P < 2e-16")

plot_d_inner <- ggplot(df_d, aes(x = age_group, y = logP, fill = age_group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 23, fill = "white", size = 2) +
  scale_fill_manual(values = age_colors, drop = FALSE) +
  labs(
    x = "Age Group",
    y = expression("Log Abundance of " * italic("Prevotella stercorea"))
  ) +
  annotate("text", x = 2, y = max(df_d$logP, na.rm = TRUE) * 1.05,
           label = anova_label, size = 5, hjust = 0.5, color = "black") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.3),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(-0.1, "cm")  # inward ticks
  )

plot_d <- ggdraw() +
  draw_label(expression(italic("Prevotella stercorea") * " Abundance by Age Group"),
             x = 0.5, y = 0.98, hjust = 0.5, vjust = 1,
             size = 18, fontface = "bold", color = "black") +   # stays bold
  draw_plot(plot_d_inner, y = 0, height = 0.95)

# ================================
# PLOT 2: Illness Bar Plot
# ================================
df_e <- df %>%
  filter(!is.na(Ill_Status), !is.na(age_group))

df_e$Ill_Status <- factor(df_e$Ill_Status, levels = c("Not-Ill", "Ill"))
df_e$age_group <- factor(df_e$age_group, levels = age_levels)

# Prepare barplot data 
ill_df <- df_e %>%
  group_by(age_group, Ill_Status) %>%
  summarise(unique_ids = n_distinct(`Randomisation No`), .groups = "drop") %>%
  rename(count = unique_ids)

# Chi-squared or Fisher test 
ill_tab <- table(df_e$age_group, df_e$Ill_Status)
if (all(dim(ill_tab) == c(length(age_levels), 2)) && all(rowSums(ill_tab) > 0)) {
  chisq_res <- chisq.test(ill_tab)
  test_label <- paste0("Chi-sq(", chisq_res$parameter, ") = ", 
                       round(chisq_res$statistic, 2), 
                       ", P = ", formatC(chisq_res$p.value, format = "e", digits = 2))
} else {
  fisher_res <- fisher.test(ill_tab)
  test_label <- paste0("Fisher's Exact test P = ", formatC(fisher_res$p.value, format = "e", digits = 2))
}

plot_e_inner <- ggplot(ill_df, aes(x = age_group, y = count, fill = age_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  facet_wrap(~Ill_Status) +
  scale_fill_manual(values = age_colors, drop = FALSE) +
  labs(
    x = "Age Group",
    y = "Number of Children"
  ) +
  coord_cartesian(ylim = c(0, max(ill_df$count) * 1.15), clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.margin = margin(t = 30, r = 10, b = 10, l = 10),
    panel.grid = element_blank(),
    strip.text = element_text(size = 14),
    axis.line = element_line(color = "black", size = 0.3),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(-0.1, "cm")  # inward ticks
  )

plot_e <- ggdraw() +
  draw_label("Ill Status Distribution by Age Group",
             x = 0.5, y = 0.98, hjust = 0.5, vjust = 1,
             size = 18, fontface = "plain", color = "black") +   # no longer bold
  draw_label(test_label,
             x = 0.5, y = 0.92, hjust = 0.5, vjust = 1,
             size = 14, fontface = "plain", color = "black") +
  draw_plot(plot_e_inner, y = 0, height = 0.9)

# ================================
# Export side-by-side plots: PDF and EMF
# ================================
combined_plots <- grid.arrange(plot_d, plot_e, nrow = 1)

# PDF export 
pdf("C:/Users/oofordile/Desktop/Figure_AgeGroup_Analysis_Aligned_withAxis.pdf", width = 12, height = 6)
grid.draw(combined_plots)
dev.off()

# EMF export (Windows only) 
emf("C:/Users/oofordile/Desktop/Figure_AgeGroup_Analysis_Aligned_withAxis.emf", width = 12, height = 6)
grid.draw(combined_plots)
dev.off()
