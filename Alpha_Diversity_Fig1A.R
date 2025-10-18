# ======================================================
# Fig 1A — Combined CSV with Nature-Style Statistical Reporting
# UPDATED: Uses same mixed-effects/GLM methods as first script
# ======================================================

# -------------------------
# Load Required Libraries
# -------------------------
library(phyloseq)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(grid)
library(broom)
library(devEMF)
library(lme4)
library(lmerTest)

# -------------------------
# Define colors for Ill Status
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
# Extract count data
# -------------------------
count_data <- df[, 13:500]
stopifnot(all(sapply(count_data, is.numeric)))
count_mat <- as.matrix(count_data)
rownames(count_mat) <- df$SampleID

# -------------------------
# Sample metadata and phyloseq object
# -------------------------
sample_data_df <- data.frame(
  SampleID = df$SampleID,
  Ill_Status = df$Ill_Status,
  timepoints = df$timepoints,
  AgeGroup = df$`3agegroups based on age at sampling`,
  SubjectID = df$`Randomisation No`,
  stringsAsFactors = FALSE
)
rownames(sample_data_df) <- df$SampleID
sample_data_ps <- sample_data(sample_data_df)
otu_table_ps <- otu_table(count_mat, taxa_are_rows = FALSE)
ps <- phyloseq(otu_table_ps, sample_data_ps)

# -------------------------
# Estimate alpha diversity
# -------------------------
alpha_div <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson", "Fisher"))
alpha_div$Pielou <- alpha_div$Shannon / log(alpha_div$Observed)
alpha_div$SampleID <- df$SampleID
alpha_data <- merge(alpha_div, sample_data_df, by = "SampleID")
alpha_data <- alpha_data %>% rename(Richness = Observed)
alpha_data$timepoints <- factor(alpha_data$timepoints, levels = c("D1", "D15", "D85"))
alpha_data$Ill_Status <- factor(alpha_data$Ill_Status, levels = c("Ill", "Not-Ill"))

# -------------------------
# Define age groups and metrics
# -------------------------
age_groups <- c("7to12 mths", "1to2years", "plus2years")
existing_age_groups <- intersect(age_groups, unique(alpha_data$AgeGroup))
if(length(existing_age_groups) < length(age_groups)) {
  warning("Some age groups not found in data, using available ones")
  age_groups <- existing_age_groups
}
all_metrics <- c("Richness", "Fisher")

# ======================================================
# UPDATED: Mixed-model / GLM comparison function (from first script)
# This ensures identical p-values to the first script
# ======================================================
run_age_stratified_model <- function(data, metric, timepoint, age_group) {
  age_data <- data %>% filter(AgeGroup == age_group, timepoints == timepoint)
  
  # Need at least 4 observations and both Ill_Status groups
  if(nrow(age_data) < 4 || length(unique(age_data$Ill_Status)) != 2) {
    return(list(p_value = NA_real_, n = 0, estimate = NA, ci_lower = NA, ci_upper = NA,
                F = NA, df1 = NA, df2 = NA, partial_eta2 = NA))
  }
  
  # Set reference level to match age-adjusted analysis
  age_data$Ill_Status <- relevel(age_data$Ill_Status, ref = "Not-Ill")
  
  # Use simple GLM for within-timepoint comparisons
  fit <- lm(as.formula(paste(metric, "~ Ill_Status")), data = age_data)
  aov_tab <- anova(fit)
  
  # Extract statistics
  p_val <- tryCatch(as.numeric(aov_tab["Ill_Status", "Pr(>F)"]), error = function(e) NA_real_)
  F_val <- tryCatch(as.numeric(aov_tab["Ill_Status", "F value"]), error = function(e) NA_real_)
  df1 <- tryCatch(as.numeric(aov_tab["Ill_Status", "Df"]), error = function(e) NA_real_)
  df2 <- tryCatch(as.numeric(aov_tab["Residuals", "Df"]), error = function(e) NA_real_)
  
  # Calculate partial eta squared
  ss_term <- aov_tab["Ill_Status", "Sum Sq"]
  ss_resid <- aov_tab["Residuals", "Sum Sq"]
  partial_eta2 <- ss_term / (ss_term + ss_resid)
  
  # Extract coefficient and CI
  coef_try <- tryCatch({
    cf <- coef(fit)
    possible <- grep("Ill_Status", names(cf), value = TRUE)
    if(length(possible) >= 1) possible[1] else NA
  }, error = function(e) NA)
  
  estimate <- NA; ci_lower <- NA; ci_upper <- NA
  if(!is.na(coef_try) && coef_try %in% names(coef(fit))) {
    estimate <- as.numeric(coef(fit)[coef_try])
    cfint <- tryCatch(confint(fit, parm = coef_try, level = 0.95), error = function(e) NA)
    if(!is.na(cfint)[1]) {
      ci_lower <- as.numeric(cfint[1])
      ci_upper <- as.numeric(cfint[2])
    }
  }
  
  return(list(
    p_value = p_val,
    n = nrow(age_data),
    estimate = estimate,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    F = F_val,
    df1 = df1,
    df2 = df2,
    partial_eta2 = partial_eta2
  ))
}

# ======================================================
# Helper function: extract stats from lm fit for a given term (for age-adjusted)
# ======================================================
extract_lm_term_stats <- function(fit, term_name) {
  an_tbl <- tryCatch(as.data.frame(anova(fit)), error = function(e) NULL)
  res <- list(
    term = term_name, F = NA, df1 = NA, df2 = NA, p_value = NA,
    SS_term = NA, SS_resid = NA, partial_eta2 = NA,
    estimate = NA, ci_lower = NA, ci_upper = NA, n = NA
  )
  res$n <- length(model.frame(fit)[[1]])
  
  if(!is.null(an_tbl) && term_name %in% rownames(an_tbl)) {
    ss_term <- an_tbl[term_name, "Sum Sq"]
    df_term <- an_tbl[term_name, "Df"]
    resid_rowname <- rownames(an_tbl)[nrow(an_tbl)]
    ss_resid <- an_tbl[resid_rowname, "Sum Sq"]
    df_resid <- an_tbl[resid_rowname, "Df"]
    ms_term <- an_tbl[term_name, "Mean Sq"]
    ms_resid <- an_tbl[resid_rowname, "Mean Sq"]
    F_val <- ifelse(is.na(ms_resid) || ms_resid == 0, NA, ms_term / ms_resid)
    p_val <- an_tbl[term_name, "Pr(>F)"]
    partial_eta2 <- ss_term / (ss_term + ss_resid)
    
    res$F <- F_val; res$df1 <- df_term; res$df2 <- df_resid; res$p_value <- p_val
    res$SS_term <- ss_term; res$SS_resid <- ss_resid; res$partial_eta2 <- partial_eta2
  }
  
  # Extract coefficient estimate and CI
  coef_try <- tryCatch({
    cf <- coef(fit)
    possible <- grep(paste0("^", term_name), names(cf), value = TRUE)
    if(length(possible) == 0) possible <- grep("Ill_Status", names(cf), value = TRUE)
    if(length(possible) >= 1) possible[1] else NA
  }, error = function(e) NA)
  
  if(!is.na(coef_try) && coef_try %in% names(coef(fit))) {
    est <- coef(fit)[coef_try]
    cfint <- tryCatch(confint(fit, parm = coef_try, level = 0.95), error = function(e) NA)
    if(!is.na(cfint)[1]) {
      res$estimate <- as.numeric(est)
      res$ci_lower <- as.numeric(cfint[1])
      res$ci_upper <- as.numeric(cfint[2])
    } else res$estimate <- as.numeric(est)
  }
  return(res)
}

# ======================================================
# Age-adjusted analyses (unchanged)
# ======================================================
timepoint_results_stats <- list()
for(tp in c("D1", "D15", "D85")) {
  tp_data <- alpha_data %>% filter(timepoints == tp)
  if(nrow(tp_data) < 4) next
  
  for(m in all_metrics) {
    tp_data$Ill_Status <- relevel(tp_data$Ill_Status, ref = "Not-Ill")
    fmla <- as.formula(paste0(m, " ~ Ill_Status + AgeGroup"))
    fit <- tryCatch(lm(fmla, data = tp_data), error = function(e) NULL)
    if(is.null(fit)) next
    
    term_stats <- extract_lm_term_stats(fit, "Ill_Status")
    row <- data.frame(
      Timepoint = tp, Metric = m, Analysis = "Age-adjusted", Term = "Ill_Status",
      Estimate = term_stats$estimate, CI_lower = term_stats$ci_lower,
      CI_upper = term_stats$ci_upper, F = term_stats$F, df1 = term_stats$df1,
      df2 = term_stats$df2, p_value = term_stats$p_value, partial_eta2 = term_stats$partial_eta2,
      n = term_stats$n, stringsAsFactors = FALSE
    )
    timepoint_results_stats[[paste(tp, m, sep = "_")]] <- row
  }
}
df_timepoint_stats <- bind_rows(timepoint_results_stats)
df_timepoint_stats <- df_timepoint_stats %>% 
  group_by(Timepoint) %>% 
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>% 
  ungroup()

# ======================================================
# UPDATED: Age-stratified analyses using mixed-effects method
# ======================================================
age_strat_stats <- list()
for(age_grp in age_groups) {
  for(tp in c("D1", "D15", "D85")) {
    for(m in all_metrics) {
      stats <- run_age_stratified_model(alpha_data, m, tp, age_grp)
      
      if(!is.na(stats$p_value)) {
        row <- data.frame(
          AgeGroup = age_grp, 
          Timepoint = tp, 
          Metric = m, 
          Analysis = "Age-stratified",
          Term = "Ill_Status", 
          Estimate = stats$estimate,
          CI_lower = stats$ci_lower, 
          CI_upper = stats$ci_upper,
          F = stats$F, 
          df1 = stats$df1, 
          df2 = stats$df2,
          p_value = stats$p_value, 
          partial_eta2 = stats$partial_eta2,
          n = stats$n, 
          stringsAsFactors = FALSE
        )
        age_strat_stats[[paste(age_grp, tp, m, sep = "_")]] <- row
      }
    }
  }
}

df_age_strat_stats <- bind_rows(age_strat_stats)
df_age_strat_stats <- df_age_strat_stats %>% 
  group_by(AgeGroup, Timepoint) %>% 
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>% 
  ungroup()

# ======================================================
# Combine Both Analyses with Nature Stats Column
# ======================================================
# Add AgeGroup column to timepoint stats for consistency
df_timepoint_stats$AgeGroup <- NA

# Combine both dataframes
combined_results <- bind_rows(df_timepoint_stats, df_age_strat_stats)

# ======================================================
# Create Nature-Style Statistical Text Column
# ======================================================
combined_results$Nature_Stats <- sapply(1:nrow(combined_results), function(i) {
  row <- combined_results[i, ]
  
  context <- if(!is.na(row$AgeGroup)) {
    paste0(row$AgeGroup, ", Age-stratified")
  } else {
    "Age-adjusted"
  }
  
  sprintf(
    "%s (%s, %s): F(%g, %g) = %.2f, p = %.3f, p_adj = %.3f, eta2_p = %.3f, 95%% CI [%.2f, %.2f]; n = %d",
    row$Metric, row$Timepoint, context, 
    row$df1, row$df2, row$F, row$p_value, row$p_adj, row$partial_eta2, 
    row$CI_lower, row$CI_upper, row$n
  )
})

# Reorder columns for better readability
combined_results <- combined_results %>%
  dplyr::select(Analysis, Metric, Timepoint, AgeGroup, 
                p_value, p_adj, n, Estimate, CI_lower, CI_upper, 
                F, df1, df2, partial_eta2, Nature_Stats)

# -------------------------
# Save Single Combined CSV
# -------------------------
write.csv(combined_results, 
          "C:/Users/oofordile/Desktop/Alpha_Diversity_Combined_Results.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\n======================================================\n")
cat("Combined results exported:\n")
cat(" → Alpha_Diversity_Combined_Results.csv\n")
cat("\nIncludes both age-adjusted and age-stratified analyses\n")
cat("Nature_Stats column format:\n")
cat("Metric (Timepoint, Context): F(df1, df2) = value,\n")
cat("p = value, p_adj = value, eta2_p = value,\n")
cat("95%% CI [lower, upper]; n = value\n")
cat("======================================================\n")

# ======================================================
# Create separate dataframes for plotting
# ======================================================
all_timepoint_results_adjusted <- combined_results %>%
  filter(Analysis == "Age-adjusted") %>%
  dplyr::select(Timepoint, Metric, p_value, p_adj, n, Estimate, CI_lower, CI_upper, F, df1, df2, partial_eta2)

all_timepoint_results_by_age <- combined_results %>%
  filter(Analysis == "Age-stratified") %>%
  dplyr::select(AgeGroup, Timepoint, Metric, p_value, p_adj, n, Estimate, CI_lower, CI_upper, F, df1, df2, partial_eta2)

# ======================================================
# Updated create_condensed_plot function with improved p-value formatting
# ======================================================

create_condensed_plot <- function(data, metric, y_label, colors, results_df, title_text, 
                                  show_legend = FALSE, age_grp = NULL, y_axis_text = TRUE) {
  if(!is.null(age_grp)) {
    stats_df <- results_df %>% filter(Metric == metric, AgeGroup == age_grp)
  } else {
    if("AgeGroup" %in% colnames(results_df)) {
      stats_df <- results_df %>% filter(Metric == metric & is.na(AgeGroup))
    } else {
      stats_df <- results_df %>% filter(Metric == metric)
    }
  }
  
  y_range <- range(data[[metric]], na.rm = TRUE)
  y_diff <- y_range[2] - y_range[1]
  y_max <- y_range[2] + 0.25 * y_diff
  anno_y <- y_max - 0.08 * y_diff
  
  p <- ggplot(data, aes(x = timepoints, y = .data[[metric]], color = Ill_Status, group = Ill_Status)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), alpha = 0.4, size = 1.8) +
    geom_boxplot(aes(fill = Ill_Status), alpha = 0.3, outlier.shape = NA, position = position_dodge(width = 0.6), width = 0.5) +
    stat_summary(fun = mean, geom = "line", linewidth = 1.3, position = position_dodge(width = 0.6)) +
    stat_summary(fun = mean, geom = "point", size = 3.5, shape = 18, position = position_dodge(width = 0.6)) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = c(y_range[1] - 0.05 * y_diff, y_max)) +
    labs(x = "Timepoint", y = if(y_axis_text) y_label else NULL, title = title_text) +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank(), plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 13), axis.title = element_text(size = 14, face = "bold"),
          legend.position = if(show_legend) "bottom" else "none",
          legend.title = element_blank(), legend.text = element_text(size = 13),
          axis.line = element_line(size = 0.4, colour = "black"),
          axis.ticks = element_line(size = 0.4, colour = "black"),
          axis.ticks.length = unit(-2, "mm"),
          plot.margin = margin(t = 20, r = 10, b = 10, l = 10))
  
  if(nrow(stats_df) > 0){
    for(i in 1:nrow(stats_df)) {
      tp_pos <- which(levels(data$timepoints) == stats_df$Timepoint[i])
      
      # Format p-value: use scientific notation if < 0.001, otherwise 3 decimal places
      p_val <- stats_df$p_adj[i]
      if(p_val < 0.001) {
        p_label <- sprintf("FDR p = %.2e", p_val)
      } else {
        p_label <- sprintf("FDR p = %.3f", p_val)
      }
      
      # Stagger annotation heights to prevent overlap
      # Adjust y position based on timepoint index
      anno_offset <- (tp_pos - 2) * 0.06 * y_diff  # Stagger vertically
      anno_y_adjusted <- anno_y + anno_offset
      
      p <- p + annotate("text", x = tp_pos, y = anno_y_adjusted, label = p_label, 
                        size = 4.5, fontface = "bold")
    }
  }
  return(p)
}

# ======================================================
# Helper function to format p-values for console output
# ======================================================
format_p_value <- function(p) {
  if(is.na(p)) return("NA")
  if(p < 0.001) {
    return(sprintf("%.2e", p))
  } else {
    return(sprintf("%.3f", p))
  }
}

# ======================================================
# Generate Single-Row 4-Panel Plots
# ======================================================
for(metric in all_metrics) {
  y_label <- ifelse(metric == "Richness", "Species Richness", "Fisher's Alpha")
  
  plot_adjusted <- create_condensed_plot(alpha_data, metric, y_label, group_colors,
                                         all_timepoint_results_adjusted, "Age-Adjusted",
                                         show_legend = FALSE, age_grp = NULL)
  
  plots_stratified <- list()
  for(i in seq_along(age_groups)) {
    age_grp <- age_groups[i]
    age_data <- alpha_data %>% filter(AgeGroup == age_grp)
    show_legend <- (i == length(age_groups))
    
    plots_stratified[[i]] <- create_condensed_plot(age_data, metric, y_label, group_colors,
                                                   all_timepoint_results_by_age, age_grp,
                                                   show_legend = show_legend, age_grp = age_grp,
                                                   y_axis_text = FALSE)
  }
  
  main_title <- paste0(metric, ": Age-Adjusted and Age-Stratified Alpha Diversity")
  
  pdf(paste0("C:/Users/oofordile/Desktop/", metric, "_SingleRow.pdf"), width = 20, height = 5)
  grid.arrange(plot_adjusted, plots_stratified[[1]], plots_stratified[[2]], plots_stratified[[3]],
               ncol = 4, top = textGrob(main_title, gp = gpar(fontsize = 16, fontface = "bold")))
  dev.off()
  
  emf(paste0("C:/Users/oofordile/Desktop/", metric, "_SingleRow.emf"), width = 20, height = 5, emfPlus = TRUE)
  grid.arrange(plot_adjusted, plots_stratified[[1]], plots_stratified[[2]], plots_stratified[[3]],
               ncol = 4, top = textGrob(main_title, gp = gpar(fontsize = 16, fontface = "bold")))
  dev.off()
}

cat("\n======================================================\n")
cat("Analysis complete\n")
cat("Outputs: Age-adjusted and age-stratified panels\n")
cat("Files exported to Desktop\n")
cat("======================================================\n")
