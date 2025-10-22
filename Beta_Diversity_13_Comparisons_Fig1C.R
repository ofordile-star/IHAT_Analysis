# ====================================================== 
# Load Required Libraries
# ======================================================
library(phyloseq)
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ape)
library(patchwork)
library(devEMF)
library(boot)

set.seed(12345)

# ======================================================
# Helper Function for P-value Formatting
# ======================================================
format_p <- function(p, threshold = 0.001) {
  if (is.na(p)) return("NA")
  if (p >= threshold) {
    return(sprintf("%.3f", p))
  } else {
    return(format(p, scientific = TRUE, digits = 3))
  }
}

# ======================================================
# Load and Prepare Data
# ======================================================
file_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
merged_data <- read.csv(file_path, stringsAsFactors = FALSE)

# Prepare metadata
analysis_data <- data.frame(
  Randomisation_No = merged_data$Randomisation.No,
  Ill_Status = merged_data$Ill_Status,
  AgeGroup = merged_data$X3agegroups.based.on.age.at.sampling,
  Timepoints = merged_data$timepoints,
  SampleID = paste0(merged_data$Randomisation.No, "_", merged_data$timepoints),
  row.names = paste0(merged_data$Randomisation.No, "_", merged_data$timepoints)
)

# Harmonize AgeGroup levels
analysis_data$AgeGroup <- gsub("7to12 mths", "7 to 12 months", analysis_data$AgeGroup)
analysis_data$AgeGroup <- gsub("1to2years", "1 to 2 years", analysis_data$AgeGroup)
analysis_data$AgeGroup <- gsub("plus2years", "plus 2 years", analysis_data$AgeGroup)
analysis_data$AgeGroup <- factor(analysis_data$AgeGroup, 
                                 levels = c("7 to 12 months", "1 to 2 years", "plus 2 years"))

# Create phyloseq object
microbiome_counts <- merged_data[, 13:500]
rownames(microbiome_counts) <- rownames(analysis_data)
otu_tab <- otu_table(as.matrix(microbiome_counts), taxa_are_rows = FALSE)
physeq <- phyloseq(otu_tab, sample_data(analysis_data))

# Factor conversion
sample_data(physeq)$Ill_Status <- factor(sample_data(physeq)$Ill_Status, levels = c("Ill", "Not-Ill"))
sample_data(physeq)$Timepoints <- factor(sample_data(physeq)$Timepoints)
sample_data(physeq)$Randomisation_No <- factor(sample_data(physeq)$Randomisation_No)

# Normalize to relative abundances
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# ======================================================
# Bootstrap CI Function
# ======================================================
calculate_r2_ci <- function(dist_mat, group_var, age_var, n_boot = 199, strata_var = NULL) {
  tryCatch({
    boot_r2 <- function(data, indices) {
      boot_dist <- as.dist(as.matrix(dist_mat)[indices, indices])
      boot_group <- group_var[indices]
      boot_age <- age_var[indices]
      
      if (length(unique(boot_group)) < 2) return(NA)
      
      if (!is.null(strata_var)) {
        boot_strata <- strata_var[indices]
        ad <- adonis2(boot_dist ~ boot_group + boot_age, permutations = 49, strata = boot_strata)
      } else {
        ad <- adonis2(boot_dist ~ boot_group + boot_age, permutations = 49)
      }
      return(ad$R2[1])
    }
    
    boot_result <- boot(data = 1:length(group_var), statistic = boot_r2, R = n_boot)
    ci <- boot.ci(boot_result, type = "perc")
    
    if (!is.null(ci)) c(ci$percent[4], ci$percent[5]) else c(NA, NA)
  }, error = function(e) c(NA, NA))
}

# ======================================================
# Define Comparisons
# ======================================================
comparisons <- list(
  list(name = "D1 Ill vs D1 Not-Ill", filter = expression(Timepoints == "D1" & Ill_Status %in% c("Ill","Not-Ill")), group_var = "Ill_Status", paired = FALSE),
  list(name = "D15 Ill vs D15 Not-Ill", filter = expression(Timepoints == "D15" & Ill_Status %in% c("Ill","Not-Ill")), group_var = "Ill_Status", paired = FALSE),
  list(name = "D85 Ill vs D85 Not-Ill", filter = expression(Timepoints == "D85" & Ill_Status %in% c("Ill","Not-Ill")), group_var = "Ill_Status", paired = FALSE),
  list(name = "D1 vs D15 Not-Ill", filter = expression(Timepoints %in% c("D1","D15") & Ill_Status == "Not-Ill"), group_var = "Timepoints", paired = TRUE),
  list(name = "D15 vs D85 Not-Ill", filter = expression(Timepoints %in% c("D15","D85") & Ill_Status == "Not-Ill"), group_var = "Timepoints", paired = TRUE),
  list(name = "D1 vs D85 Not-Ill", filter = expression(Timepoints %in% c("D1","D85") & Ill_Status == "Not-Ill"), group_var = "Timepoints", paired = TRUE),
  list(name = "D1 vs D15 Ill", filter = expression(Timepoints %in% c("D1","D15") & Ill_Status == "Ill"), group_var = "Timepoints", paired = TRUE),
  list(name = "D15 vs D85 Ill", filter = expression(Timepoints %in% c("D15","D85") & Ill_Status == "Ill"), group_var = "Timepoints", paired = TRUE),
  list(name = "D1 vs D85 Ill", filter = expression(Timepoints %in% c("D1","D85") & Ill_Status == "Ill"), group_var = "Timepoints", paired = TRUE),
  list(name = "D1 Ill vs D85 Not-Ill", filter = expression((Ill_Status == "Ill" & Timepoints == "D1") | (Ill_Status == "Not-Ill" & Timepoints == "D85")), group_var = "Ill_Status", paired = FALSE),
  list(name = "D15 Ill vs D85 Not-Ill", filter = expression((Ill_Status == "Ill" & Timepoints == "D15") | (Ill_Status == "Not-Ill" & Timepoints == "D85")), group_var = "Ill_Status", paired = FALSE),
  list(name = "D1 Not-Ill vs D85 Ill", filter = expression((Ill_Status == "Not-Ill" & Timepoints == "D1") | (Ill_Status == "Ill" & Timepoints == "D85")), group_var = "Ill_Status", paired = FALSE),
  list(name = "D15 Not-Ill vs D85 Ill", filter = expression((Ill_Status == "Not-Ill" & Timepoints == "D15") | (Ill_Status == "Ill" & Timepoints == "D85")), group_var = "Ill_Status", paired = FALSE)
)

# ======================================================
# Colors and Shapes
# ======================================================
ill_colors <- c("Ill" = "#E31A1C", "Not-Ill" = "#377EB8")
gray_color <- "gray50"
age_shapes <- c("7 to 12 months" = 22, "1 to 2 years" = 21, "plus 2 years" = 24)

# ======================================================
# Analysis Loop
# ======================================================
results_list <- list()

for (i in seq_along(comparisons)) {
  comp <- comparisons[[i]]
  message("Analysing ", comp$name)
  
  ps_subset <- subset_samples(physeq_rel, eval(comp$filter))
  ps_subset <- prune_samples(sample_sums(ps_subset) > 0, ps_subset)
  if (nsamples(ps_subset) == 0) next
  
  otu_mat_sub <- as(otu_table(ps_subset), "matrix")
  if (taxa_are_rows(ps_subset)) otu_mat_sub <- t(otu_mat_sub)
  dist_mat <- vegdist(otu_mat_sub, method = "bray")
  sample_df_sub <- as(sample_data(ps_subset), "data.frame")
  
  # PERMANOVA
  set.seed(12345 + i)
  if (comp$paired) {
    adonis_res <- adonis2(dist_mat ~ sample_df_sub[[comp$group_var]] + sample_df_sub$AgeGroup, 
                          permutations = 999, strata = sample_df_sub$Randomisation_No)
  } else {
    adonis_res <- adonis2(dist_mat ~ sample_df_sub[[comp$group_var]] + sample_df_sub$AgeGroup, 
                          permutations = 999)
  }
  
  # Calculate 95% CI for R²
  message("  Calculating bootstrap CI for R²...")
  ci_r2 <- if (comp$paired) {
    calculate_r2_ci(dist_mat, sample_df_sub[[comp$group_var]], sample_df_sub$AgeGroup, 
                    n_boot = 199, strata_var = sample_df_sub$Randomisation_No)
  } else {
    calculate_r2_ci(dist_mat, sample_df_sub[[comp$group_var]], sample_df_sub$AgeGroup, n_boot = 199)
  }
  
  # Age-stratified PERMANOVA
  age_perm_df <- bind_rows(lapply(levels(sample_df_sub$AgeGroup), function(ag) {
    idx <- sample_df_sub$AgeGroup == ag
    if (sum(idx) < 2 || length(unique(sample_df_sub[[comp$group_var]][idx])) < 2) {
      return(data.frame(AgeGroup = ag, F = NA, R2 = NA, P = NA, Df_num = NA, Df_denom = NA, CI_lower = NA, CI_upper = NA))
    }
    
    dist_sub <- as.dist(as.matrix(dist_mat)[idx, idx])
    grp_sub <- sample_df_sub[[comp$group_var]][idx]
    set.seed(12345 + i + which(levels(sample_df_sub$AgeGroup) == ag))
    ad_res <- adonis2(dist_sub ~ grp_sub, permutations = 999)
    
    ci_age <- tryCatch({
      boot_r2_age <- function(data, indices) {
        boot_dist <- as.dist(as.matrix(dist_sub)[indices, indices])
        boot_group <- grp_sub[indices]
        if (length(unique(boot_group)) < 2) return(NA)
        ad <- adonis2(boot_dist ~ boot_group, permutations = 49)
        return(ad$R2[1])
      }
      boot_result <- boot(data = 1:length(grp_sub), statistic = boot_r2_age, R = 199)
      ci <- boot.ci(boot_result, type = "perc")
      if (!is.null(ci)) c(ci$percent[4], ci$percent[5]) else c(NA, NA)
    }, error = function(e) c(NA, NA))
    
    data.frame(AgeGroup = ag, F = ad_res$F[1], R2 = ad_res$R2[1], P = ad_res$`Pr(>F)`[1],
               Df_num = ad_res$Df[1], Df_denom = ad_res$Df[length(ad_res$Df) - 1],
               CI_lower = ci_age[1], CI_upper = ci_age[2])
  }))
  
  # Beta dispersion
  beta_disp <- betadisper(dist_mat, sample_df_sub[[comp$group_var]])
  kruskal_disp <- kruskal.test(beta_disp$distances ~ sample_df_sub[[comp$group_var]])
  
  # PCoA
  pcoa_res <- ape::pcoa(dist_mat)
  pcoa_df <- data.frame(
    Axis.1 = pcoa_res$vectors[,1],
    Axis.2 = pcoa_res$vectors[,2],
    SampleID = rownames(sample_df_sub)
  ) %>% left_join(sample_df_sub, by = "SampleID")
  
  # Store results
  results_list[[comp$name]] <- list(
    comp_name = comp$name,
    is_canvas1 = comp$name %in% c("D1 Ill vs D1 Not-Ill","D15 Ill vs D15 Not-Ill","D85 Ill vs D85 Not-Ill"),
    f_stat = adonis_res$F[1],
    r2 = adonis_res$R2[1],
    p_val = adonis_res$`Pr(>F)`[1],
    df_num = adonis_res$Df[1],
    df_denom = adonis_res$Df[length(adonis_res$Df) - 1],
    sum_of_sqs = adonis_res$SumOfSqs[1],
    ci_r2 = ci_r2,
    kruskal_stat = kruskal_disp$statistic,
    kruskal_df = kruskal_disp$parameter,
    disp_p = kruskal_disp$p.value,
    age_perm_df = age_perm_df,
    pcoa_res = pcoa_res,
    pcoa_df = pcoa_df,
    centroid_df = data.frame(DistanceToCentroid = beta_disp$distances, 
                             Group = sample_df_sub[[comp$group_var]]),
    sample_df_sub = sample_df_sub,
    group_var = comp$group_var
  )
}

# ======================================================
# FDR CORRECTION
# ======================================================
cat("\n=== Applying FDR Correction ===\n")

# Overall PERMANOVA
all_overall_pvals <- sapply(results_list, function(x) x$p_val)
all_overall_fdr <- p.adjust(all_overall_pvals, method = "BH")
for (i in seq_along(results_list)) results_list[[i]]$p_val_fdr <- all_overall_fdr[i]

# Age-stratified
all_age_pvals <- unlist(lapply(results_list, function(x) x$age_perm_df$P))
all_age_fdr <- p.adjust(all_age_pvals, method = "BH")
idx <- 1
for (i in seq_along(results_list)) {
  n <- nrow(results_list[[i]]$age_perm_df)
  results_list[[i]]$age_perm_df$P_FDR <- all_age_fdr[idx:(idx + n - 1)]
  idx <- idx + n
}

# Dispersion
all_disp_pvals <- sapply(results_list, function(x) x$disp_p)
all_disp_fdr <- p.adjust(all_disp_pvals, method = "BH")
for (i in seq_along(results_list)) results_list[[i]]$disp_p_fdr <- all_disp_fdr[i]

canvas1_overall_p_fdr <- results_list[["D1 Ill vs D1 Not-Ill"]]$p_val_fdr

cat("Overall PERMANOVA p-values (FDR-corrected):\n")
for (nm in names(results_list)) {
  cat("  ", nm, ": raw p =", format_p(results_list[[nm]]$p_val), 
      ", FDR p =", format_p(results_list[[nm]]$p_val_fdr), "\n")
}

# ======================================================
# CREATE PLOTS
# ======================================================
create_plot <- function(res) {
  # Define colors
  present_groups <- unique(as.character(res$sample_df_sub[[res$group_var]]))
  if (res$group_var == "Timepoints") {
    tp_levels <- sort(present_groups)
    base_color <- if (grepl("Not-Ill", res$comp_name)) "#377EB8" else "#E31A1C"
    colors <- setNames(c(base_color, gray_color)[1:length(tp_levels)], tp_levels)
  } else {
    tmp <- intersect(names(ill_colors), present_groups)
    colors <- if(length(tmp)>0) ill_colors[tmp] else setNames(rep(gray_color, length(present_groups)), present_groups)
  }
  
  # Facet labels with FDR
  facet_labels <- setNames(
    paste0(res$age_perm_df$AgeGroup, " (FDR p = ", sapply(res$age_perm_df$P_FDR, format_p), ")"), 
    res$age_perm_df$AgeGroup
  )
  
  base_font <- if (res$is_canvas1) 19 else 22
  x_breaks <- pretty(res$pcoa_df$Axis.1, n = 5)
  y_breaks <- pretty(res$pcoa_df$Axis.2, n = 5)
  
  # PCoA plot
  subtitle <- if (res$is_canvas1) {
    res$comp_name
  } else {
    paste0(res$comp_name, "\nPERMANOVA FDR p = ", format_p(res$p_val_fdr))
  }
  
  p_pcoa <- ggplot(res$pcoa_df, aes(Axis.1, Axis.2, fill = .data[[res$group_var]], shape = AgeGroup)) +
    geom_point(size = if(res$is_canvas1) 3.5 else 3, 
               alpha = if(res$is_canvas1) 0.5 else 0.8, 
               color = "black", stroke = 0.3) +
    scale_fill_manual(values = colors, na.value = gray_color) +
    scale_shape_manual(values = age_shapes) +
    guides(fill = "none", shape = "none") +
    facet_wrap(~AgeGroup, labeller = labeller(AgeGroup = facet_labels)) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks) +
    labs(subtitle = subtitle,
         x = paste0("PCoA1 (", round(res$pcoa_res$values$Relative_eig[1]*100,1), "%)"),
         y = paste0("PCoA2 (", round(res$pcoa_res$values$Relative_eig[2]*100,1), "%)")) +
    theme_minimal(base_size = base_font) +
    theme(
      plot.subtitle = element_text(size = base_font + 1, face = "plain", hjust = 0, margin = margin(b = 6)),
      strip.text = element_text(size = base_font - 3, face = "plain", margin = margin(b = 4)),
      strip.background = element_rect(fill = NA, colour = NA),
      panel.grid = element_blank(),
      axis.text = element_text(angle = 0, hjust = 0.5),
      plot.margin = margin(2, 5, 2, 5),
      panel.spacing = unit(0.4, "lines"),
      axis.line = element_line(color = "black", size = 0.3),
      axis.ticks = element_line(color = "black", size = 0.3),
      axis.ticks.length = unit(-2, "mm")
    )
  
  # Centroid plot
  p_centroid <- ggplot(res$centroid_df, aes(Group, DistanceToCentroid, fill = Group)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = if(res$is_canvas1) 0.6 else 0.8, 
                color = "black", size = 1.8) +
    scale_fill_manual(values = colors, na.value = gray_color) +
    labs(y = "Distance to centroid", x = NULL,
         subtitle = paste0("Kruskal–Wallis FDR p = ", format_p(res$disp_p_fdr))) +
    theme_minimal(base_size = base_font) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      plot.subtitle = element_text(size = base_font - 2, hjust = 0.5, margin = margin(0, 0, 4, 0)),
      axis.line = element_line(color = "black", size = 0.3),
      axis.ticks = element_line(color = "black", size = 0.3),
      axis.ticks.length = unit(-2, "mm"),
      plot.margin = margin(2, 5, 2, 5)
    )
  
  p_pcoa + p_centroid + plot_layout(widths = c(3, 1))
}

canvas1_plots <- lapply(results_list[sapply(results_list, function(x) x$is_canvas1)], create_plot)
canvas2_plots <- lapply(results_list[!sapply(results_list, function(x) x$is_canvas1)], create_plot)

# ======================================================
# Export CSV
# ======================================================
results_summary <- bind_rows(lapply(names(results_list), function(nm) {
  res <- results_list[[nm]]
  ci_text <- if (!is.na(res$ci_r2[1]) && !is.na(res$ci_r2[2])) {
    paste0("[", sprintf("%.3f", res$ci_r2[1]), ", ", sprintf("%.3f", res$ci_r2[2]), "]")
  } else "[CI not available]"
  
  data.frame(
    Comparison = nm,
    PERMANOVA_Df_numerator = res$df_num,
    PERMANOVA_Df_denominator = res$df_denom,
    PERMANOVA_SumOfSqs = round(res$sum_of_sqs, 4),
    PERMANOVA_R2 = round(res$r2, 4),
    PERMANOVA_F = round(res$f_stat, 4),
    PERMANOVA_P = format_p(res$p_val),
    PERMANOVA_P_FDR = format_p(res$p_val_fdr),
    R2_95CI_lower = round(res$ci_r2[1], 4),
    R2_95CI_upper = round(res$ci_r2[2], 4),
    `Kruskal-Wallis_H` = round(res$kruskal_stat, 4),
    `Kruskal-Wallis_Df` = res$kruskal_df,
    `Kruskal-Wallis_P` = format_p(res$disp_p),
    `Kruskal-Wallis_P_FDR` = format_p(res$disp_p_fdr),
    Nature_Format_PERMANOVA = paste0("F(", res$df_num, ", ", res$df_denom, ") = ", 
                                     sprintf("%.2f", res$f_stat), ", p = ", format_p(res$p_val),
                                     ", FDR p = ", format_p(res$p_val_fdr), ", R² = ", 
                                     sprintf("%.3f", res$r2), ", 95% CI ", ci_text),
    Nature_Format_KruskalWallis = paste0("H(", res$kruskal_df, ") = ", sprintf("%.2f", res$kruskal_stat),
                                         ", p = ", format_p(res$disp_p), ", FDR p = ", 
                                         format_p(res$disp_p_fdr))
  )
}))
write.csv(results_summary, "C:/Users/oofordile/Desktop/BetaDiversity_Overall_PERMANOVA_Summary_FDR.csv", row.names = FALSE)

age_stratified_df <- bind_rows(lapply(names(results_list), function(nm) {
  res <- results_list[[nm]]
  row <- data.frame(Comparison = nm)
  for (i in 1:nrow(res$age_perm_df)) {
    ag <- res$age_perm_df$AgeGroup[i]
    if (is.na(res$age_perm_df$F[i])) {
      text <- "NA"
    } else {
      ci_text <- if (!is.na(res$age_perm_df$CI_lower[i]) && !is.na(res$age_perm_df$CI_upper[i])) {
        paste0("95% CI [", sprintf("%.3f", res$age_perm_df$CI_lower[i]), ", ", 
               sprintf("%.3f", res$age_perm_df$CI_upper[i]), "]")
      } else "95% CI [not available]"
      text <- paste0("F(", res$age_perm_df$Df_num[i], ", ", res$age_perm_df$Df_denom[i], ") = ",
                     sprintf("%.2f", res$age_perm_df$F[i]), ", p = ", format_p(res$age_perm_df$P[i]),
                     ", FDR p = ", format_p(res$age_perm_df$P_FDR[i]), ", R² = ",
                     sprintf("%.3f", res$age_perm_df$R2[i]), ", ", ci_text)
    }
    col_name <- switch(ag, "7 to 12 months" = "7–12 months", "1 to 2 years" = "1–2 years", 
                       "plus 2 years" = ">2 years")
    row[[col_name]] <- text
  }
  row
}))
write.csv(age_stratified_df[, c("Comparison", "7–12 months", "1–2 years", ">2 years")], 
          "C:/Users/oofordile/Desktop/BetaDiversity_AgeStratified_PERMANOVA_Summary_FDR.csv", row.names = FALSE)

# ======================================================
# Save Plots
# ======================================================
# Canvas 1
permanova_plot <- ggplot() +
  labs(title = paste0("All Age-Adjusted PERMANOVA FDR p = ", format_p(canvas1_overall_p_fdr))) +
  theme_void() +
  theme(plot.title = element_text(size = 21, hjust = 0.5), plot.margin = margin(10, 0, 10, 0))

final_plot <- wrap_plots(permanova_plot, wrap_plots(canvas1_plots, ncol = 1), 
                         ncol = 1, heights = c(0.06, length(canvas1_plots)))

pdf("C:/Users/oofordile/Desktop/BetaDiversity_PCoA_Centroid_SquareCanvas_FDR.pdf", width = 18, height = 12)
print(final_plot)
dev.off()

emf("C:/Users/oofordile/Desktop/BetaDiversity_PCoA_Centroid_SquareCanvas_FDR.emf", width = 18, height = 12)
print(final_plot)
dev.off()

# Canvas 2A/B
canvas2_names <- names(canvas2_plots)
split_idx <- ceiling(length(canvas2_names)/2)

pdf("C:/Users/oofordile/Desktop/BetaDiversity_PCoA_Centroid_Canvas2A_FDR.pdf", width = 22, height = 28)
print(wrap_plots(canvas2_plots[canvas2_names[1:split_idx]], ncol = 1))
dev.off()

emf("C:/Users/oofordile/Desktop/BetaDiversity_PCoA_Centroid_Canvas2A_FDR.emf", width = 22, height = 28)
print(wrap_plots(canvas2_plots[canvas2_names[1:split_idx]], ncol = 1))
dev.off()

pdf("C:/Users/oofordile/Desktop/BetaDiversity_PCoA_Centroid_Canvas2B_FDR.pdf", width = 22, height = 28)
print(wrap_plots(canvas2_plots[canvas2_names[(split_idx+1):length(canvas2_names)]], ncol = 1))
dev.off()

emf("C:/Users/oofordile/Desktop/BetaDiversity_PCoA_Centroid_Canvas2B_FDR.emf", width = 22, height = 28)
print(wrap_plots(canvas2_plots[canvas2_names[(split_idx+1):length(canvas2_names)]], ncol = 1))
dev.off()

cat("\n=== Analysis Complete ===\n")
cat("Files saved with FDR correction applied to:\n")
cat("- Overall PERMANOVA p-values (", length(results_list), " comparisons)\n")
cat("- Age-stratified PERMANOVA p-values (", length(all_age_pvals), " comparisons)\n")
cat("- Kruskal-Wallis dispersion p-values (", length(results_list), " comparisons)\n")
cat("\nCanvas 1: General FDR p at top | Canvas 2A/B: Individual FDR p in subtitles\n")
cat("\nP-values formatted as: 3 decimal places for p ≥ 0.001, scientific notation for p < 0.001\n")