# ======================================================
# LOAD REQUIRED LIBRARIES
# ======================================================
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(cowplot)
library(broom)
library(lubridate)
library(grid)
library(rstatix)
library(devEMF)

# ======================================================
# LOAD DATA
# ======================================================
data_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
df <- read_csv(data_path)

# ======================================================
# DEFINE AGE GROUPS AND COLORS
# ======================================================
age_levels <- c("7to12 mths", "1to2years", "plus2years")
df$age_group <- factor(df$`3agegroups based on age at sampling`, levels = age_levels)

age_colors <- c(
  "7to12 mths" = "#A6D854",
  "1to2years"  = "#FDB863",
  "plus2years" = "#80B1D3"
)

illness_colors <- c("Not-Ill" = "#6BAED6", "Ill" = "#FC8D62")

# ======================================================
# PANEL A: Prevotella stercorea ABUNDANCE BY AGE GROUP
# ======================================================
df_a <- df %>%
  filter(!is.na(`Prevotella_stercorea_7.05`),
         !is.na(Ill_Status),
         !is.na(age_group)) %>%
  mutate(logP = log(`Prevotella_stercorea_7.05` + 1))

df_a$Ill_Status <- factor(df_a$Ill_Status, levels = c("Not-Ill", "Ill"))

# ---- ANOVA ----
anova_res <- aov(logP ~ age_group, data = df_a)
anova_summary <- summary(anova_res)[[1]]
anova_F <- anova_summary[["F value"]][1]
anova_p <- anova_summary[["Pr(>F)"]][1]
df1 <- anova_summary[["Df"]][1]
df2 <- anova_summary[["Df"]][2]
eta2 <- anova_summary[["Sum Sq"]][1] / sum(anova_summary[["Sum Sq"]])

anova_p_label <- if(anova_p < 2e-16) "p < 2e-16" else paste0("p = ", signif(anova_p,3))

# ---- Logistic Regression Interaction P ----
df_a <- df_a %>% mutate(Ill_num = ifelse(Ill_Status == "Ill", 1, 0))
glm_fit <- glm(Ill_num ~ logP * age_group, data = df_a, family = binomial)
glm_summary <- broom::tidy(glm_fit, conf.int = TRUE, exponentiate = TRUE)

interaction_terms <- glm_summary %>% filter(grepl("logP:age_group", term))
interaction_p_label <- if(nrow(interaction_terms) > 0) {
  paste0("Interaction P = ", signif(max(interaction_terms$p.value, na.rm = TRUE), 3))
} else {"Interaction P = NA"}

combined_label_a <- paste0(
  "F(", df1, ", ", df2, ") = ", round(anova_F,2),
  ", ", anova_p_label,
  ", eta2_p = ", round(eta2,3),
  "; ", interaction_p_label
)

# ---- Pairwise Tukey ----
tukey_res <- TukeyHSD(anova_res)
tukey_df <- as.data.frame(tukey_res$age_group)
tukey_df$comparison <- rownames(tukey_df)
tukey_df <- tukey_df %>%
  mutate(
    Nature_Report = paste0(
      "Diff(", comparison, ") = ", round(diff,2),
      ", p = ", signif(`p adj`,3),
      ", 95% CI [", round(lwr,2), ", ", round(upr,2), "]"
    )
  )
write_csv(tukey_df, "C:/Users/oofordile/Desktop/Pairwise_AgeGroup_Tukey.csv")

# ---- Logistic Regression CSV ----
glm_summary <- glm_summary %>%
  mutate(Nature_Report = paste0(
    term, ": OR = ", round(estimate,3),
    ", 95% CI [", round(conf.low,3), ", ", round(conf.high,3), "], p = ", signif(p.value,3)
  ))
write_csv(glm_summary, "C:/Users/oofordile/Desktop/AgeStratified_LogisticRegression_4dp.csv")

# ---- Combined CSV ----
anova_tukey_combined <- bind_rows(
  data.frame(Comparison="ANOVA_overall", Nature_Report=combined_label_a, stringsAsFactors=FALSE),
  tukey_df %>% select(Comparison=comparison, Nature_Report)
)
write_csv(anova_tukey_combined, "C:/Users/oofordile/Desktop/Prevotella_AgeGroup_NatureStatistics.csv")

# ======================================================
# PANEL A PLOT
# ======================================================
pairwise_data <- tukey_df %>%
  tidyr::separate(comparison, into = c("group2","group1"), sep="-") %>%
  mutate(y.position = max(df_a$logP) * c(1.05,1.15,1.25)[1:n()],
         label = paste0("p = ", signif(`p adj`,3)))

plot_a_inner <- ggplot(df_a, aes(x = age_group, y = logP, fill = age_group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA,
               aes(group = interaction(age_group, Ill_Status)),
               position = position_dodge(width = 0.8)) +
  geom_point(aes(group = interaction(age_group, Ill_Status)),
             alpha = 0.4, size = 0.8, color = "black",
             position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.8)) +
  geom_text(data = data.frame(age_group = rep(age_levels, each = 2),
                              Ill_Status = rep(c("Not-Ill", "Ill"), times = 3),
                              y = min(df_a$logP) - 0.8,
                              label = rep(c("Ill", "Not-Ill"), times = 3)),
            aes(x = age_group, y = y, label = label,
                group = interaction(age_group, Ill_Status)),
            position = position_dodge(width = 0.8),
            size = 4.5, vjust = 1, fontface = "bold") +
  stat_pvalue_manual(pairwise_data,
                     label = "label",
                     tip.length = 0.02,
                     bracket.size = 0.4,
                     size = 3.5) +
  scale_fill_manual(values = age_colors) +
  labs(x = "Age Group",
       y = expression("Log(count + 1) of " * italic("Prevotella stercorea"))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 14) +
  theme(legend.position="none",
        panel.grid=element_blank(),
        axis.line=element_line(color="black", linewidth=0.4),
        axis.ticks=element_line(color="black", linewidth=0.4),
        axis.ticks.length=unit(-0.15,"cm"),
        axis.title=element_text(size=15),
        axis.text=element_text(size=13),
        plot.margin = margin(t=5,r=10,b=20,l=10))

plot_a <- ggdraw() +
  draw_label(expression(italic("Prevotella stercorea") *
                          " Abundance Across Age Groups"),
             x=0.5, y=0.975, hjust=0.5, vjust=1,
             size=17, fontface="bold") +
  draw_label(combined_label_a,
             x=0.5, y=0.915, hjust=0.5, vjust=1,
             size=14) +
  draw_plot(plot_a_inner, y=0, height=0.88)

# ======================================================
# PANEL B: ILLNESS STATUS DISTRIBUTION BY AGE (UNCHANGED)
# ======================================================
df_b <- df %>% filter(!is.na(Ill_Status), !is.na(age_group))
df_b$Ill_Status <- factor(df_b$Ill_Status, levels=c("Not-Ill","Ill"))

ill_df <- df_b %>%
  group_by(age_group, Ill_Status) %>%
  summarise(N_children=n_distinct(`Randomisation No`), .groups="drop")

ill_tab <- table(df_b$age_group, df_b$Ill_Status)
if(all(dim(ill_tab)==c(length(age_levels),2)) & all(rowSums(ill_tab)>0)){
  chisq_res <- chisq.test(ill_tab)
  test_label <- paste0("Chi-sq(", chisq_res$parameter, ") = ",
                       round(chisq_res$statistic,2),
                       ", P = ", signif(chisq_res$p.value,3))
} else {
  fisher_res <- fisher.test(ill_tab)
  test_label <- paste0("Fisher's Exact test P = ", signif(fisher_res$p.value,3))
}

plot_b_inner <- ggplot(ill_df, aes(x=age_group, y=N_children, fill=age_group)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9)) +
  facet_wrap(~Ill_Status) +
  geom_text(data = ill_df %>% group_by(Ill_Status) %>%
              summarise(y = max(N_children)*1.05) %>%
              mutate(age_group="1to2years"),
            aes(x=age_group, y=y, label=Ill_Status),
            size=5.5, fontface="bold", vjust=0) +
  scale_fill_manual(values=age_colors) +
  labs(x="Age Group", y="Number of Children") +
  coord_cartesian(ylim=c(0,max(ill_df$N_children)*1.15), clip="off") +
  theme_minimal(base_size=14) +
  theme(legend.position="none",
        panel.grid=element_blank(),
        strip.text=element_blank(),
        axis.line=element_line(color="black", linewidth=0.4),
        axis.ticks=element_line(color="black", linewidth=0.4),
        axis.ticks.length=unit(-0.15,"cm"),
        axis.title=element_text(size=15),
        axis.text=element_text(size=13))

plot_b <- ggdraw() +
  draw_label("Illness Status Distribution by Age Group", x=0.5, y=0.975,
             hjust=0.5, vjust=1, size=17) +
  draw_label(test_label, x=0.5, y=0.915, hjust=0.5, vjust=1, size=14) +
  draw_plot(plot_b_inner, y=0.035, height=0.88)

# ======================================================
# COMBINE PANELS AND EXPORT
# ======================================================
combined_plots <- plot_grid(plot_a, plot_b, ncol=2, rel_widths=c(1,1))

pdf("C:/Users/oofordile/Desktop/Figure_Combined_AgeGroup_Analysis_v2.pdf",
    width=14, height=6)
print(combined_plots)
dev.off()

emf("C:/Users/oofordile/Desktop/Figure_Combined_AgeGroup_Analysis_v2.emf",
    width=14, height=6)
print(combined_plots)
dev.off()

cat("Plots exported successfully!\n")
cat("- PDF: Figure_Combined_AgeGroup_Analysis_v2.pdf\n")
cat("- EMF: Figure_Combined_AgeGroup_Analysis_v2.emf\n")
cat("- CSVs: Pairwise_AgeGroup_Tukey.csv, AgeStratified_LogisticRegression_4dp.csv, Prevotella_AgeGroup_NatureStatistics.csv\n")
