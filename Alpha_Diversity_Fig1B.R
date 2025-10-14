# ======================================================
# Fig 1B — D1/D15 Ill vs D85 Not-Ill: Species Richness
# Enhanced with Nature-Style Statistics (eta2_p, UTF-8)
# ======================================================

# -------------------------
# Load Required Libraries
# -------------------------
library(phyloseq)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(devEMF)
library(broom)

# -------------------------
# Define colors
# -------------------------
group_colors <- c("Ill" = "#D73027", "Not-Ill" = "#4575B4")

# -------------------------
# Load data
# -------------------------
data_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
df <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE)
df <- df %>% filter(!is.na(Ill_Status) & Ill_Status != "")
df$SampleID <- paste0(df$`Randomisation No`, "_", df$timepoints)
stopifnot(!any(duplicated(df$SampleID)))

# -------------------------
# Count data & phyloseq
# -------------------------
count_data <- df[, 13:500]
count_mat <- as.matrix(count_data)
rownames(count_mat) <- df$SampleID

sample_data_df <- data.frame(
  SampleID = df$SampleID,
  Ill_Status = df$Ill_Status,
  timepoints = df$timepoints,
  AgeGroup = df$`3agegroups based on age at sampling`,
  SubjectID = df$`Randomisation No`
)
rownames(sample_data_df) <- df$SampleID

sample_data_ps <- sample_data(sample_data_df)
otu_table_ps <- otu_table(count_mat, taxa_are_rows = FALSE)
ps <- phyloseq(otu_table_ps, sample_data_ps)

# -------------------------
# Alpha diversity
# -------------------------
alpha_div <- estimate_richness(ps, measures = c("Observed"))
alpha_div$SampleID <- df$SampleID
alpha_data <- merge(alpha_div, sample_data_df, by = "SampleID")
alpha_data <- alpha_data %>% rename(Richness = Observed)
alpha_data$timepoints <- factor(alpha_data$timepoints, levels = c("D1", "D15", "D85"))
alpha_data$Ill_Status <- factor(alpha_data$Ill_Status, levels = c("Ill", "Not-Ill"))

# -------------------------
# Age groups
# -------------------------
age_groups <- c("7to12 mths", "1to2years", "plus2years")
existing_age_groups <- intersect(age_groups, unique(alpha_data$AgeGroup))
if(length(existing_age_groups) < length(age_groups)) age_groups <- existing_age_groups

# -------------------------
# Prepare comparison dataset
# -------------------------
prepare_comparison <- function(data, age_grp = NULL) {
  if(!is.null(age_grp)) data <- data %>% filter(AgeGroup == age_grp)
  d1_ill <- data %>% filter(timepoints == "D1" & Ill_Status == "Ill")
  d15_ill <- data %>% filter(timepoints == "D15" & Ill_Status == "Ill")
  d85_ni <- data %>% filter(timepoints == "D85" & Ill_Status == "Not-Ill")
  
  comp_d1 <- rbind(d1_ill, d85_ni)
  comp_d1$Comparison <- "D1 Ill vs D85 Not-Ill"
  
  comp_d15 <- rbind(d15_ill, d85_ni)
  comp_d15$Comparison <- "D15 Ill vs D85 Not-Ill"
  
  rbind(comp_d1, comp_d15)
}

# -------------------------
# Compute detailed statistics (Nature format, ASCII eta2_p)
# -------------------------
compute_detailed_stats <- function(df, analysis_type = "Overall") {
  comparisons <- unique(df$Comparison)
  all_stats <- data.frame()
  formatted_lines <- c()
  
  for(cmp in comparisons) {
    tmp <- df %>% filter(Comparison == cmp)
    if(length(unique(tmp$Ill_Status)) == 2) {
      tmp$Ill_Status <- factor(tmp$Ill_Status, levels = c("Not-Ill", "Ill"))
      fit <- lm(Richness ~ Ill_Status, data = tmp)
      
      coef_summary <- summary(fit)$coefficients
      coef_est <- coef_summary["Ill_StatusIll", "Estimate"]
      ci <- confint(fit, "Ill_StatusIll")
      p_val <- coef_summary["Ill_StatusIll", "Pr(>|t|)"]
      
      aov_tab <- anova(fit)
      f_stat <- aov_tab["Ill_Status", "F value"]
      df1 <- aov_tab["Ill_Status", "Df"]
      df2 <- aov_tab["Residuals", "Df"]
      
      ss_effect <- aov_tab["Ill_Status", "Sum Sq"]
      ss_resid <- aov_tab["Residuals", "Sum Sq"]
      eta2_p <- ss_effect / (ss_effect + ss_resid)
      
      means <- tmp %>% 
        group_by(Ill_Status) %>% 
        summarise(n = n(),
                  Mean = mean(Richness, na.rm = TRUE),
                  SD = sd(Richness, na.rm = TRUE),
                  .groups = 'drop')
      
      timepoint <- strsplit(cmp, " ")[[1]][1]
      total_n <- sum(means$n)
      
      formatted <- sprintf(
        "Richness (%s, %s): F(%d, %d) = %.2f, p = %.3f, eta2_p = %.3f, 95%% CI [%.2f, %.2f]; n = %d",
        timepoint, analysis_type, df1, df2, f_stat, p_val, eta2_p, ci[1], ci[2], total_n
      )
      formatted_lines <- c(formatted_lines, formatted)
      
      for(i in 1:nrow(means)) {
        stats <- data.frame(
          Analysis = analysis_type,
          Comparison = cmp,
          Timepoint = timepoint,
          Ill_Status = means$Ill_Status[i],
          n = means$n[i],
          Mean = means$Mean[i],
          SD = means$SD[i],
          Estimate = coef_est,
          CI_lower = ci[1],
          CI_upper = ci[2],
          F = f_stat,
          df1 = df1,
          df2 = df2,
          eta2_p = eta2_p,
          p_value = p_val,
          Nature_Formatted_Stat = formatted,
          stringsAsFactors = FALSE
        )
        all_stats <- rbind(all_stats, stats)
      }
    }
  }
  
  cat("\n--- Nature-Style Statistics ---\n")
  cat(paste(formatted_lines, collapse = "\n"), "\n-------------------------------\n")
  
  all_stats
}

# -------------------------
# Compute stats for all analyses
# -------------------------
stats_overall <- compute_detailed_stats(prepare_comparison(alpha_data), "Age-Adjusted")

stats_age_stratified <- data.frame()
for(age_grp in age_groups) {
  stats_age <- compute_detailed_stats(prepare_comparison(alpha_data, age_grp), age_grp)
  stats_age_stratified <- rbind(stats_age_stratified, stats_age)
}

all_detailed_stats <- rbind(stats_overall, stats_age_stratified) %>%
  group_by(Analysis, Comparison) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(across(c(Mean, SD, Estimate, CI_lower, CI_upper, F, eta2_p, p_value, p_adj),
                ~round(.x, 6))) %>%
  mutate(across(c(df1, df2), ~round(.x, 0))) %>%
  select(Analysis, Comparison, Timepoint, Ill_Status, n, Mean, SD, 
         Estimate, CI_lower, CI_upper, F, df1, df2, eta2_p, p_value, p_adj, Nature_Formatted_Stat)

# -------------------------
# Save detailed statistics CSV (UTF-8)
# -------------------------
write.csv(all_detailed_stats, 
          "C:/Users/oofordile/Desktop/SpeciesRichness_D1D15_vs_D85_DetailedStats.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

# -------------------------
# Full statistics for plotting
# -------------------------
compute_full_stats <- function(df) {
  comparisons <- unique(df$Comparison)
  all_stats <- data.frame()
  for(cmp in comparisons) {
    tmp <- df %>% filter(Comparison == cmp)
    if(length(unique(tmp$Ill_Status)) == 2) {
      tmp$Ill_Status <- factor(tmp$Ill_Status, levels = c("Not-Ill", "Ill"))
      fit <- lm(Richness ~ Ill_Status, data = tmp)
      coef_est <- coef(fit)["Ill_StatusIll"]
      ci <- confint(fit, "Ill_StatusIll")
      p_val <- summary(fit)$coefficients["Ill_StatusIll", "Pr(>|t|)"]
      means <- tmp %>% group_by(Ill_Status) %>% summarise(n = n(),
                                                          Mean = mean(Richness, na.rm = TRUE),
                                                          SD = sd(Richness, na.rm = TRUE),
                                                          .groups='drop')
      stats <- data.frame(
        Comparison = cmp,
        Ill_Status = means$Ill_Status,
        n = means$n,
        Mean = means$Mean,
        SD = means$SD,
        Estimate = ifelse(means$Ill_Status[1]=="Ill", coef_est, -coef_est),
        CI_lower = ifelse(means$Ill_Status[1]=="Ill", ci[1], -ci[2]),
        CI_upper = ifelse(means$Ill_Status[1]=="Ill", ci[2], -ci[1]),
        p_value = p_val,
        p_adj = p.adjust(p_val, method = "BH")
      )
      all_stats <- rbind(all_stats, stats)
    }
  }
  all_stats
}

create_comparison_plot <- function(comp_data, age_grp = NULL, show_legend = FALSE, y_axis_text = TRUE) {
  pvals <- compute_full_stats(comp_data)
  p <- ggplot(comp_data, aes(x = Comparison, y = Richness, fill = Ill_Status, color = Ill_Status)) +
    geom_boxplot(position = position_dodge(width = 0.6), width = 0.4, outlier.shape = NA, alpha = 0.3) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6), alpha = 0.5, size = 1.8) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    labs(y = if(y_axis_text) "Species Richness" else NULL, x = NULL, title = ifelse(is.null(age_grp), "Age-Adjusted", age_grp)) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 25, hjust = 1, size = 12),
      axis.title.y = element_text(size = 14, face = "bold"),
      legend.position = if(show_legend) "bottom" else "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 13),
      axis.line = element_line(size = 0.4, colour = "black"),
      axis.ticks = element_line(size = 0.4, colour = "black"),
      axis.ticks.length = unit(-2, "mm"),
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
    )
  for(i in 1:nrow(pvals)) {
    cmp <- pvals$Comparison[i]
    y_max <- max(comp_data$Richness[comp_data$Comparison == cmp], na.rm = TRUE)
    p <- p + annotate("text", x = cmp, y = y_max + 0.08*(y_max), 
                      label = sprintf("FDR p=%.3f", pvals$p_adj[i]), size = 5, fontface = "bold")
  }
  p
}

# -------------------------
# Generate plots
# -------------------------
plots <- list()
plots[[1]] <- create_comparison_plot(prepare_comparison(alpha_data), show_legend = FALSE, y_axis_text = TRUE)
for(i in seq_along(age_groups)) {
  plots[[i+1]] <- create_comparison_plot(prepare_comparison(alpha_data, age_groups[i]),
                                         age_grp = age_groups[i],
                                         show_legend = (i==length(age_groups)),
                                         y_axis_text = FALSE)
}

main_title <- "D1/D15 Ill vs D85 Not-Ill: Species Richness"

pdf("C:/Users/oofordile/Desktop/SpeciesRichness_D1D15_vs_D85_SingleRow_LargerP.pdf", width = 20, height = 5)
grid.arrange(grobs = plots, ncol = 4, top = textGrob(main_title, gp = gpar(fontsize = 16, fontface = "bold")))
dev.off()

emf("C:/Users/oofordile/Desktop/SpeciesRichness_D1D15_vs_D85_SingleRow_LargerP.emf", width = 20, height = 5, emfPlus = TRUE)
grid.arrange(grobs = plots, ncol = 4, top = textGrob(main_title, gp = gpar(fontsize = 16, fontface = "bold")))
dev.off()

# -------------------------
# Save original format stats CSV (UTF-8)
# -------------------------
stats_full <- compute_full_stats(prepare_comparison(alpha_data))
write.csv(stats_full, "C:/Users/oofordile/Desktop/SpeciesRichness_D1D15_vs_D85_FullStats.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\n======================================================\n")
cat("Single-row D1/D15 Ill vs D85 Not-Ill plot complete!\n")
cat("Nature-style stats (eta2_p, CI, df, p) exported in detailed CSV.\n")
cat("Full statistics CSV exported: SpeciesRichness_D1D15_vs_D85_FullStats.csv\n")
cat("Detailed statistics CSV exported: SpeciesRichness_D1D15_vs_D85_DetailedStats.csv\n")
cat("======================================================\n")
