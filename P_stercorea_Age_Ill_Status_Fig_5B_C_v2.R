# ======================================================
# Load Required Libraries
# ======================================================
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(vcd)
library(gridExtra)
library(forcats)
library(cowplot)
library(devEMF)   # For EMF export
library(grid)     # For unit()
library(broom)

# ======================================================
# Load Data
# ======================================================
data_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
df <- read_csv(data_path)

# ======================================================
# Define Age Groups and Colors
# ======================================================
age_levels <- c("7to12 mths", "1to2years", "plus2years")
df$age_group <- factor(df$`3agegroups based on age at sampling`, levels = age_levels)

age_colors <- c(
  "7to12 mths" = "#A6D854",   # light green
  "1to2years"  = "#FDB863",   # light orange
  "plus2years" = "#80B1D3"    # light blue
)

# ======================================================
# PANEL A: P. stercorea by Illness Status and Age Group
# ======================================================

# ---- ANOVA dataset (with pseudocount) ----
df_anova <- df %>%
  filter(!is.na(`Prevotella_stercorea_7.05`),
         !is.na(Ill_Status),
         !is.na(age_group)) %>%
  mutate(
    logP = log(`Prevotella_stercorea_7.05` + 1) # pseudocount
  )

anova_res <- aov(logP ~ age_group, data = df_anova)
anova_summary <- summary(anova_res)[[1]]

anova_F   <- anova_summary[["F value"]][1]
anova_p   <- anova_summary[["Pr(>F)"]][1]
df1       <- anova_summary[["Df"]][1]      # numerator df
df2       <- anova_summary[["Df"]][2]      # denominator df

# Print to console
cat("ANOVA results:\n")
cat("F(", df1, ", ", df2, ") = ", round(anova_F, 2),
    ", P = ", formatC(anova_p, format = "e", digits = 2), "\n", sep = "")
cat("------------------------------------------------------------\n")

# Create ANOVA label
if (anova_p < 2e-16) {
  anova_label <- paste0("ANOVA F(", df1, ", ", df2, ") = ",
                        round(anova_F, 2), ", P < 2e-16")
} else {
  anova_label <- paste0("ANOVA F(", df1, ", ", df2, ") = ",
                        round(anova_F, 2),
                        ", P = ", formatC(anova_p, format = "e", digits = 2))
}

# ---- Logistic regression dataset (also with pseudocount) ----
df_corr <- df %>%
  filter(!is.na(`Prevotella_stercorea_7.05`),
         !is.na(Ill_Status),
         !is.na(age_group)) %>%
  mutate(
    logP = log(`Prevotella_stercorea_7.05` + 1), # pseudocount
    Ill_num = ifelse(Ill_Status == "Ill", 1, 0)
  )

glm_fit <- glm(Ill_num ~ logP * age_group, data = df_corr, family = binomial)
glm_summary <- broom::tidy(glm_fit, conf.int = TRUE, exponentiate = TRUE)

# Save logistic regression output
write_csv(glm_summary %>%
            mutate(across(where(is.numeric), ~round(., 4))),
          "C:/Users/oofordile/Desktop/AgeStratified_LogisticRegression_4dp.csv")

# Extract interaction terms
interaction_terms <- glm_summary %>%
  filter(grepl("logP:age_group", term))

write_csv(interaction_terms %>%
            mutate(across(where(is.numeric), ~round(., 4))),
          "C:/Users/oofordile/Desktop/AgeStratified_InteractionTerms_4dp.csv")

# Interaction label
if (nrow(interaction_terms) > 0) {
  interaction_p_val <- max(interaction_terms$p.value, na.rm = TRUE)
  interaction_p_label <- paste0("Interaction P = ",
                                formatC(interaction_p_val, format = "f", digits = 3))
} else {
  interaction_p_label <- "Interaction P = NA"
}

combined_label_a <- paste0(anova_label, "; ", interaction_p_label)

# ---- Plot Panel A ----
plot_a_inner <- ggplot(df_corr, aes(x = Ill_Status, y = logP, fill = age_group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA,
               position = position_dodge(width = 0.8)) +
  geom_jitter(alpha = 0.4, size = 0.7,
              position = position_jitterdodge(jitter.width = 0.15,
                                              dodge.width = 0.8)) +
  scale_fill_manual(values = age_colors) +
  facet_wrap(~age_group) +
  labs(
    x = "Illness Status",
    y = expression("Log(count + 1) of " * italic("Prevotella stercorea"))
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    strip.text = element_text(size = 12),
    axis.line = element_line(color = "black", size = 0.3),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(-0.1, "cm")
  )

plot_a <- ggdraw() +
  draw_label(expression(italic("Prevotella stercorea") * " Abundance by Illness Status and Age Group"),
             x = 0.5, y = 0.98, hjust = 0.5, vjust = 1,
             size = 16, fontface = "bold", color = "black") +
  draw_label(combined_label_a,
             x = 0.5, y = 0.92, hjust = 0.5, vjust = 1,
             size = 12, fontface = "plain", color = "black") +
  draw_plot(plot_a_inner, y = 0, height = 0.9)

# ======================================================
# PANEL B: Illness Status Distribution by Age Group
# ======================================================

df_e <- df %>%
  filter(!is.na(Ill_Status), !is.na(age_group))

df_e$Ill_Status <- factor(df_e$Ill_Status, levels = c("Not-Ill", "Ill"))
df_e$age_group <- factor(df_e$age_group, levels = age_levels)

ill_df <- df_e %>%
  group_by(age_group, Ill_Status) %>%
  summarise(unique_ids = n_distinct(`Randomisation No`), .groups = "drop") %>%
  rename(count = unique_ids)

# Chi-square / Fisher test
ill_tab <- table(df_e$age_group, df_e$Ill_Status)
if (all(dim(ill_tab) == c(length(age_levels), 2)) && all(rowSums(ill_tab) > 0)) {
  chisq_res <- chisq.test(ill_tab)
  test_label <- paste0("Chi-sq(", chisq_res$parameter, ") = ",
                       round(chisq_res$statistic, 2),
                       ", P = ", formatC(chisq_res$p.value, format = "e", digits = 2))
} else {
  fisher_res <- fisher.test(ill_tab)
  test_label <- paste0("Fisher's Exact test P = ",
                       formatC(fisher_res$p.value, format = "e", digits = 2))
}

# ---- Plot Panel B ----
plot_b_inner <- ggplot(ill_df, aes(x = age_group, y = count, fill = age_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  facet_wrap(~Ill_Status) +
  scale_fill_manual(values = age_colors, drop = FALSE) +
  labs(
    x = "Age Group",
    y = "Number of Children"
  ) +
  coord_cartesian(ylim = c(0, max(ill_df$count) * 1.15), clip = "off") +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    strip.text = element_text(size = 12),
    axis.line = element_line(color = "black", size = 0.3),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(-0.1, "cm")
  )

plot_b <- ggdraw() +
  draw_label("Illness Status Distribution by Age Group",
             x = 0.5, y = 0.98, hjust = 0.5, vjust = 1,
             size = 16, fontface = "plain", color = "black") +
  draw_label(test_label,
             x = 0.5, y = 0.92, hjust = 0.5, vjust = 1,
             size = 12, fontface = "plain", color = "black") +
  draw_plot(plot_b_inner, y = 0, height = 0.9)

# ======================================================
# Combine Plots
# ======================================================
combined_plots <- cowplot::plot_grid(plot_a, plot_b,
                                     ncol = 2, rel_widths = c(1, 1))

# ======================================================
# Export Plots (PDF + EMF)
# ======================================================
pdf("C:/Users/oofordile/Desktop/Figure_Combined_AgeGroup_Analysis.pdf",
    width = 14, height = 6)
print(combined_plots)
dev.off()

emf("C:/Users/oofordile/Desktop/Figure_Combined_AgeGroup_Analysis.emf",
    width = 14, height = 6)
print(combined_plots)
dev.off()

cat("Analysis complete! Files exported:\n")
cat("- Figure_Combined_AgeGroup_Analysis.pdf\n")
cat("- Figure_Combined_AgeGroup_Analysis.emf\n")
cat("- AgeStratified_LogisticRegression_4dp.csv\n")
cat("- AgeStratified_InteractionTerms_4dp.csv\n")