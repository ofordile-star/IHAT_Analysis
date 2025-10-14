# ====================================================== 
# Load Required Libraries
# ======================================================
library(phyloseq)
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(readr)
library(ape)
library(patchwork)
library(devEMF)
library(grid)
library(cowplot)
library(boot)

# ======================================================
# Reproducibility: Fix seed
# ======================================================
set.seed(12345)

# ======================================================
# Load Data
# ======================================================
file_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
merged_data <- read.csv(file_path, stringsAsFactors = FALSE)

# ======================================================
# Prepare Metadata
# ======================================================
analysis_data <- data.frame(
  Randomisation_No = merged_data$Randomisation.No,
  Ill_Status = merged_data$Ill_Status,
  AgeGroup = merged_data$X3agegroups.based.on.age.at.sampling,
  Timepoints = merged_data$timepoints,
  SampleID = paste0(merged_data$Randomisation.No, "_", merged_data$timepoints),
  row.names = paste0(merged_data$Randomisation.No, "_", merged_data$timepoints)
)

# Fix AgeGroup levels
analysis_data$AgeGroup <- as.character(analysis_data$AgeGroup)
analysis_data$AgeGroup[analysis_data$AgeGroup == "7to12 mths"] <- "7 to 12 months"
analysis_data$AgeGroup[analysis_data$AgeGroup == "1to2years"] <- "1 to 2 years"
analysis_data$AgeGroup[analysis_data$AgeGroup == "plus2years"] <- "plus 2 years"
analysis_data$AgeGroup <- factor(analysis_data$AgeGroup,
                                 levels = c("7 to 12 months", "1 to 2 years", "plus 2 years"))

# ======================================================
# Prepare Microbiome Data
# ======================================================
microbiome_counts <- merged_data[, 13:500]
rownames(microbiome_counts) <- rownames(analysis_data)
stopifnot(all(rownames(microbiome_counts) == rownames(analysis_data)))

# ======================================================
# Create Phyloseq Object
# ======================================================
otu_tab <- otu_table(as.matrix(microbiome_counts), taxa_are_rows = FALSE)
sample_df <- sample_data(analysis_data)
physeq <- phyloseq(otu_tab, sample_df)

# Factor Conversion
sample_data(physeq)$Ill_Status <- factor(sample_data(physeq)$Ill_Status, levels = c("Ill", "Not-Ill"))
sample_data(physeq)$AgeGroup <- factor(sample_data(physeq)$AgeGroup,
                                       levels = c("7 to 12 months", "1 to 2 years", "plus 2 years"))
sample_data(physeq)$Timepoints <- factor(sample_data(physeq)$Timepoints)
sample_data(physeq)$Randomisation_No <- factor(sample_data(physeq)$Randomisation_No)

# ======================================================
# Normalize to Relative Abundances
# ======================================================
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# ======================================================
# OPTIMIZED: Function to calculate 95% CI for R² using bootstrap
# Reduced iterations from 999 to 199 for speed
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
    
    if (!is.null(ci)) {
      return(c(ci$percent[4], ci$percent[5]))
    } else {
      return(c(NA, NA))
    }
  }, error = function(e) {
    return(c(NA, NA))
  })
}

# ======================================================
# Define Comparisons
# ======================================================
comparisons <- list(
  list(name = "D1 Ill vs D1 Not-Ill", filter = expression(Timepoints == "D1" & Ill_Status %in% c("Ill","Not-Ill")), group_var = "Ill_Status", cross_timepoint = FALSE, paired = FALSE),
  list(name = "D15 Ill vs D15 Not-Ill", filter = expression(Timepoints == "D15" & Ill_Status %in% c("Ill","Not-Ill")), group_var = "Ill_Status", cross_timepoint = FALSE, paired = FALSE),
  list(name = "D85 Ill vs D85 Not-Ill", filter = expression(Timepoints == "D85" & Ill_Status %in% c("Ill","Not-Ill")), group_var = "Ill_Status", cross_timepoint = FALSE, paired = FALSE),
  list(name = "D1 vs D15 Not-Ill", filter = expression(Timepoints %in% c("D1","D15") & Ill_Status == "Not-Ill"), group_var = "Timepoints", cross_timepoint = TRUE, paired = TRUE),
  list(name = "D15 vs D85 Not-Ill", filter = expression(Timepoints %in% c("D15","D85") & Ill_Status == "Not-Ill"), group_var = "Timepoints", cross_timepoint = TRUE, paired = TRUE),
  list(name = "D1 vs D85 Not-Ill", filter = expression(Timepoints %in% c("D1","D85") & Ill_Status == "Not-Ill"), group_var = "Timepoints", cross_timepoint = TRUE, paired = TRUE),
  list(name = "D1 vs D15 Ill", filter = expression(Timepoints %in% c("D1","D15") & Ill_Status == "Ill"), group_var = "Timepoints", cross_timepoint = TRUE, paired = TRUE),
  list(name = "D15 vs D85 Ill", filter = expression(Timepoints %in% c("D15","D85") & Ill_Status == "Ill"), group_var = "Timepoints", cross_timepoint = TRUE, paired = TRUE),
  list(name = "D1 vs D85 Ill", filter = expression(Timepoints %in% c("D1","D85") & Ill_Status == "Ill"), group_var = "Timepoints", cross_timepoint = TRUE, paired = TRUE),
  list(name = "D1 Ill vs D85 Not-Ill", filter = expression((Ill_Status == "Ill" & Timepoints == "D1") | (Ill_Status == "Not-Ill" & Timepoints == "D85")), group_var = "Ill_Status", cross_timepoint = TRUE, paired = FALSE),
  list(name = "D15 Ill vs D85 Not-Ill", filter = expression((Ill_Status == "Ill" & Timepoints == "D15") | (Ill_Status == "Not-Ill" & Timepoints == "D85")), group_var = "Ill_Status", cross_timepoint = TRUE, paired = FALSE),
  list(name = "D1 Not-Ill vs D85 Ill", filter = expression((Ill_Status == "Not-Ill" & Timepoints == "D1") | (Ill_Status == "Ill" & Timepoints == "D85")), group_var = "Ill_Status", cross_timepoint = TRUE, paired = FALSE),
  list(name = "D15 Not-Ill vs D85 Ill", filter = expression((Ill_Status == "Not-Ill" & Timepoints == "D15") | (Ill_Status == "Ill" & Timepoints == "D85")), group_var = "Ill_Status", cross_timepoint = TRUE, paired = FALSE)
)

# ======================================================
# Colors and Shapes
# ======================================================
ill_colors <- c("Ill" = "#E31A1C", "Not-Ill" = "#377EB8")
gray_color <- "gray50"
age_shapes <- c("7 to 12 months" = 22, "1 to 2 years" = 21, "plus 2 years" = 24)

# ======================================================
# Initialize outputs
# ======================================================
results_summary <- list()
age_stratified_results <- list()
canvas1_plots <- list()
canvas2_plots <- list()
canvas1_overall_p <- NULL

# ======================================================
# Analysis Loop
# ======================================================
comp_index <- 0
for (comp in comparisons) {
  comp_index <- comp_index + 1
  comp_name <- comp$name
  filter_expr <- comp$filter
  group_var <- comp$group_var
  is_cross_timepoint <- comp$cross_timepoint
  is_paired <- comp$paired
  
  message("Analysing ", comp_name)
  
  ps_subset <- subset_samples(physeq_rel, eval(filter_expr))
  ps_subset <- prune_samples(sample_sums(ps_subset) > 0, ps_subset)
  if (nsamples(ps_subset) == 0) next
  
  otu_mat_sub <- as(otu_table(ps_subset), "matrix")
  if (taxa_are_rows(ps_subset)) otu_mat_sub <- t(otu_mat_sub)
  dist_mat <- vegdist(otu_mat_sub, method = "bray")
  sample_df_sub <- as(sample_data(ps_subset), "data.frame")
  
  # PERMANOVA
  set.seed(12345 + comp_index)
  if (is_cross_timepoint && is_paired) {
    adonis_res <- adonis2(as.dist(dist_mat) ~ sample_df_sub[[group_var]] + sample_df_sub$AgeGroup, 
                          permutations = 999, strata = sample_df_sub$Randomisation_No)
  } else {
    adonis_res <- adonis2(as.dist(dist_mat) ~ sample_df_sub[[group_var]] + sample_df_sub$AgeGroup, permutations = 999)
  }
  
  f_stat <- adonis_res$F[1]
  r2 <- adonis_res$R2[1]
  p_val <- adonis_res$`Pr(>F)`[1]
  df_num <- adonis_res$Df[1]
  df_denom <- adonis_res$Df[length(adonis_res$Df) - 1]
  
  # Calculate 95% CI for R² (OPTIMIZED: 199 iterations instead of 999)
  message("  Calculating bootstrap CI for R²...")
  if (is_cross_timepoint && is_paired) {
    ci_r2 <- calculate_r2_ci(dist_mat, sample_df_sub[[group_var]], sample_df_sub$AgeGroup, 
                             n_boot = 199, strata_var = sample_df_sub$Randomisation_No)
  } else {
    ci_r2 <- calculate_r2_ci(dist_mat, sample_df_sub[[group_var]], sample_df_sub$AgeGroup, n_boot = 199)
  }
  
  if (comp_name %in% c("D1 Ill vs D1 Not-Ill","D15 Ill vs D15 Not-Ill","D85 Ill vs D85 Not-Ill") && is.null(canvas1_overall_p)) {
    canvas1_overall_p <- p_val
  }
  
  # Age-stratified PERMANOVA (OPTIMIZED)
  age_permanova <- lapply(levels(sample_df_sub$AgeGroup), function(ag){
    idx <- sample_df_sub$AgeGroup == ag
    if (sum(idx) < 2) return(data.frame(AgeGroup = ag, AgeGroup_PERMANOVA_F = NA, AgeGroup_PERMANOVA_R2 = NA, 
                                        AgeGroup_PERMANOVA_P = NA, AgeGroup_Df_num = NA, AgeGroup_Df_denom = NA,
                                        AgeGroup_CI_lower = NA, AgeGroup_CI_upper = NA))
    dist_sub <- as.dist(as.matrix(dist_mat)[idx, idx])
    grp_sub <- sample_df_sub[[group_var]][idx]
    if (length(unique(grp_sub)) > 1) {
      set.seed(12345 + comp_index + which(levels(sample_df_sub$AgeGroup) == ag))
      ad_res <- adonis2(dist_sub ~ grp_sub, permutations = 999)
      
      # Bootstrap CI for age-stratified (OPTIMIZED: 199 iterations)
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
      
      data.frame(AgeGroup = ag, 
                 AgeGroup_PERMANOVA_F = ad_res$F[1], 
                 AgeGroup_PERMANOVA_R2 = ad_res$R2[1], 
                 AgeGroup_PERMANOVA_P = ad_res$`Pr(>F)`[1],
                 AgeGroup_Df_num = ad_res$Df[1],
                 AgeGroup_Df_denom = ad_res$Df[length(ad_res$Df) - 1],
                 AgeGroup_CI_lower = ci_age[1],
                 AgeGroup_CI_upper = ci_age[2])
    } else {
      data.frame(AgeGroup = ag, AgeGroup_PERMANOVA_F = NA, AgeGroup_PERMANOVA_R2 = NA, 
                 AgeGroup_PERMANOVA_P = NA, AgeGroup_Df_num = NA, AgeGroup_Df_denom = NA,
                 AgeGroup_CI_lower = NA, AgeGroup_CI_upper = NA)
    }
  })
  age_perm_df <- bind_rows(age_permanova)
  
  # Beta dispersion and Kruskal-Wallis
  beta_disp <- betadisper(dist_mat, sample_df_sub[[group_var]])
  kruskal_disp <- kruskal.test(beta_disp$distances ~ sample_df_sub[[group_var]])
  disp_p <- kruskal_disp$p.value
  kruskal_stat <- kruskal_disp$statistic
  kruskal_df <- kruskal_disp$parameter
  
  centroid_df <- data.frame(DistanceToCentroid = beta_disp$distances,
                            Group = sample_df_sub[[group_var]],
                            SampleID = rownames(sample_df_sub))
  
  # PCoA
  pcoa_res <- ape::pcoa(dist_mat)
  pcoa_df <- data.frame(Axis.1 = pcoa_res$vectors[,1],
                        Axis.2 = pcoa_res$vectors[,2],
                        SampleID = rownames(sample_df_sub))
  pcoa_df <- left_join(pcoa_df, sample_df_sub, by = "SampleID")
  set.seed(12345 + comp_index)
  pcoa_df <- pcoa_df[sample(seq_len(nrow(pcoa_df))), , drop = FALSE]
  x_breaks <- pretty(pcoa_df$Axis.1, n = 5)
  y_breaks <- pretty(pcoa_df$Axis.2, n = 5)
  
  # Define colors
  present_groups <- unique(as.character(sample_df_sub[[group_var]]))
  if (group_var == "Timepoints") {
    tp_levels <- sort(present_groups)
    if (grepl("Not-Ill", comp_name)) {
      timepoint_colors <- setNames(c("#377EB8", gray_color)[1:length(tp_levels)], tp_levels)
    } else {
      timepoint_colors <- setNames(c("#E31A1C", gray_color)[1:length(tp_levels)], tp_levels)
    }
  } else {
    tmp <- intersect(names(ill_colors), present_groups)
    timepoint_colors <- if(length(tmp)>0) ill_colors[tmp] else setNames(rep(gray_color,length(present_groups)),present_groups)
  }
  
  facet_labels <- setNames(paste0(age_perm_df$AgeGroup, " (p = ", signif(age_perm_df$AgeGroup_PERMANOVA_P, 3), ")"), age_perm_df$AgeGroup)
  
  # Determine if this is Canvas 1
  is_canvas1 <- comp_name %in% c("D1 Ill vs D1 Not-Ill","D15 Ill vs D15 Not-Ill","D85 Ill vs D85 Not-Ill")
  base_font_size <- if (is_canvas1) 19 else 22
  
  # --- PCoA Plot ---
  if (is_canvas1) {
    p_pcoa <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, fill = .data[[group_var]], shape = AgeGroup)) +
      geom_point(size = 3.5, alpha = 0.5, color = "black", stroke = 0.3) +
      scale_fill_manual(values = timepoint_colors, na.value = gray_color) +
      scale_shape_manual(values = age_shapes) +
      guides(fill = "none", shape = "none") +
      facet_wrap(~AgeGroup, labeller = labeller(AgeGroup = facet_labels), strip.position = "top") +
      scale_x_continuous(breaks = x_breaks) +
      scale_y_continuous(breaks = y_breaks) +
      labs(subtitle = comp_name,
           x = paste0("PCoA1 (", round(pcoa_res$values$Relative_eig[1]*100,1), "%)"),
           y = paste0("PCoA2 (", round(pcoa_res$values$Relative_eig[2]*100,1), "%)")) +
      theme_minimal(base_size = base_font_size) +
      theme(
        plot.subtitle = element_text(size = base_font_size + 1, face = "plain", hjust = 0, margin = margin(b = 6)),
        strip.text = element_text(size = base_font_size - 3, face = "plain", margin = margin(b = 4)),
        strip.background = element_rect(fill = NA, colour = NA),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(t = 2, r = 5, b = 2, l = 5),
        panel.spacing = unit(0.4, "lines"),
        axis.line = element_line(color = "black", size = 0.3),
        axis.ticks = element_line(color = "black", size = 0.3),
        axis.ticks.length = unit(-2, "mm")
      )
  } else {
    permanova_p_text <- paste0("PERMANOVA p = ", signif(p_val, 3))
    
    p_pcoa <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, fill = .data[[group_var]], shape = AgeGroup)) +
      geom_point(size = 3, alpha = 0.8, color = "black", stroke = 0.3) +
      scale_fill_manual(values = timepoint_colors, na.value = gray_color) +
      scale_shape_manual(values = age_shapes) +
      guides(fill = "none", shape = "none") +
      facet_wrap(~AgeGroup, labeller = labeller(AgeGroup = facet_labels), strip.position = "top") +
      scale_x_continuous(breaks = x_breaks) +
      scale_y_continuous(breaks = y_breaks) +
      labs(subtitle = paste0(comp_name, "\n", permanova_p_text),
           x = paste0("PCoA1 (", round(pcoa_res$values$Relative_eig[1]*100,1), "%)"),
           y = paste0("PCoA2 (", round(pcoa_res$values$Relative_eig[2]*100,1), "%)")) +
      theme_minimal(base_size = base_font_size) +
      theme(
        plot.subtitle = element_text(size = base_font_size + 1, face = "plain", hjust = 0, margin = margin(b = 6)),
        strip.text = element_text(size = base_font_size - 3, face = "plain", margin = margin(b = 4)),
        strip.background = element_rect(fill = NA, colour = NA),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(t = 2, r = 5, b = 2, l = 5),
        panel.spacing = unit(0.4, "lines"),
        axis.line = element_line(color = "black", size = 0.3),
        axis.ticks = element_line(color = "black", size = 0.3),
        axis.ticks.length = unit(-2, "mm")
      )
  }
  
  # --- Centroid Plot ---
  centroid_df[[group_var]] <- factor(centroid_df$Group)
  
  if (is_canvas1) {
    p_centroid <- ggplot(centroid_df, aes(x = Group, y = DistanceToCentroid, fill = Group)) +
      geom_boxplot(alpha = 0.8, outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.6, color = "black", size = 1.8) +
      scale_fill_manual(values = timepoint_colors, na.value = gray_color) +
      labs(y = "Distance to centroid", x = NULL,
           subtitle = paste0("Kruskal–Wallis p = ", signif(disp_p, 3))) +
      theme_minimal(base_size = base_font_size) +
      theme(
        legend.position = "none",
        panel.grid = element_blank(),
        plot.subtitle = element_text(size = base_font_size - 2, hjust = 0.5, margin = margin(t = 0, b = 4)),
        axis.line = element_line(color = "black", size = 0.3),
        axis.ticks = element_line(color = "black", size = 0.3),
        axis.ticks.length = unit(-2, "mm"),
        plot.margin = margin(t = 2, r = 5, b = 2, l = 5)
      )
  } else {
    p_centroid <- ggplot(centroid_df, aes(x = Group, y = DistanceToCentroid, fill = Group)) +
      geom_boxplot(alpha = 0.8, outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.8, color = "black", size = 1.8) +
      scale_fill_manual(values = timepoint_colors, na.value = gray_color) +
      labs(y = "Distance to centroid", x = NULL,
           subtitle = paste0("Kruskal–Wallis p = ", signif(disp_p, 3))) +
      theme_minimal(base_size = base_font_size) +
      theme(
        legend.position = "none",
        panel.grid = element_blank(),
        plot.subtitle = element_text(size = base_font_size - 2, hjust = 0.5, margin = margin(t = 0, b = 4)),
        axis.line = element_line(color = "black", size = 0.3),
        axis.ticks = element_line(color = "black", size = 0.3),
        axis.ticks.length = unit(-2, "mm"),
        plot.margin = margin(t = 2, r = 5, b = 2, l = 5)
      )
  }
  
  combined <- p_pcoa + p_centroid + plot_layout(widths = c(3,1))
  
  if (is_canvas1) {
    canvas1_plots[[comp_name]] <- combined
  } else {
    canvas2_plots[[comp_name]] <- combined
  }
  
  # Create Nature format string
  ci_text <- if (!is.na(ci_r2[1]) && !is.na(ci_r2[2])) {
    paste0("[", sprintf("%.3f", ci_r2[1]), ", ", sprintf("%.3f", ci_r2[2]), "]")
  } else {
    "[CI not available]"
  }
  
  nature_format <- paste0("F(", df_num, ", ", df_denom, ") = ", 
                          sprintf("%.2f", f_stat), ", p = ", 
                          sprintf("%.3f", p_val), ", R² = ", 
                          sprintf("%.3f", r2), ", 95% CI ", ci_text)
  
  kruskal_format <- paste0("H(", kruskal_df, ") = ", 
                           sprintf("%.2f", kruskal_stat), ", p = ", 
                           sprintf("%.3f", disp_p))
  
  # Store PERMANOVA results
  results_summary[[comp_name]] <- data.frame(
    Comparison = comp_name,
    PERMANOVA_Df_numerator = df_num,
    PERMANOVA_Df_denominator = df_denom,
    PERMANOVA_SumOfSqs = round(adonis_res$SumOfSqs[1], 4),
    PERMANOVA_R2 = round(r2, 4),
    PERMANOVA_F = round(f_stat, 4),
    PERMANOVA_P = round(p_val, 4),
    R2_95CI_lower = round(ci_r2[1], 4),
    R2_95CI_upper = round(ci_r2[2], 4),
    `Kruskal-Wallis_H` = round(kruskal_stat, 4),
    `Kruskal-Wallis_Df` = kruskal_df,
    `Kruskal-Wallis_P` = round(disp_p, 4),
    Nature_Format_PERMANOVA = nature_format,
    Nature_Format_KruskalWallis = kruskal_format
  )
  
  # UPDATED: Age-stratified with Nature format
  age_stratified_row <- data.frame(Comparison = comp_name, stringsAsFactors = FALSE)
  for (i in 1:nrow(age_perm_df)) {
    age_group <- age_perm_df$AgeGroup[i]
    f_val <- age_perm_df$AgeGroup_PERMANOVA_F[i]
    r2_val <- age_perm_df$AgeGroup_PERMANOVA_R2[i]
    p_val_age <- age_perm_df$AgeGroup_PERMANOVA_P[i]
    df_n <- age_perm_df$AgeGroup_Df_num[i]
    df_d <- age_perm_df$AgeGroup_Df_denom[i]
    ci_l <- age_perm_df$AgeGroup_CI_lower[i]
    ci_u <- age_perm_df$AgeGroup_CI_upper[i]
    
    if (is.na(f_val)) {
      text_val <- "NA"
    } else {
      ci_text_age <- if (!is.na(ci_l) && !is.na(ci_u)) {
        paste0("95% CI [", sprintf("%.3f", ci_l), ", ", sprintf("%.3f", ci_u), "]")
      } else {
        "95% CI [not available]"
      }
      text_val <- paste0("F(", df_n, ", ", df_d, ") = ", sprintf("%.2f", f_val), 
                         ", p = ", sprintf("%.3f", p_val_age), 
                         ", R² = ", sprintf("%.3f", r2_val), 
                         ", ", ci_text_age)
    }
    
    if (age_group == "7 to 12 months") age_stratified_row$`7–12 months` <- text_val
    if (age_group == "1 to 2 years") age_stratified_row$`1–2 years` <- text_val
    if (age_group == "plus 2 years") age_stratified_row$`>2 years` <- text_val
  }
  age_stratified_results[[comp_name]] <- age_stratified_row
}

# ======================================================
# Export results CSV
# ======================================================
results_df <- do.call(rbind, results_summary)
rownames(results_df) <- NULL
write.csv(results_df, "C:/Users/oofordile/Desktop/BetaDiversity_Overall_PERMANOVA_Summary.csv", row.names = FALSE)

age_stratified_df <- do.call(rbind, age_stratified_results)
rownames(age_stratified_df) <- NULL
age_stratified_df <- age_stratified_df[, c("Comparison", "7–12 months", "1–2 years", ">2 years")]
write.csv(age_stratified_df, "C:/Users/oofordile/Desktop/BetaDiversity_AgeStratified_PERMANOVA_Summary.csv", row.names = FALSE)

# ======================================================
# Canvas 1: Create PERMANOVA text plot and save
# ======================================================
if (length(canvas1_plots) > 0) {
  permanova_plot <- ggplot() +
    labs(title = paste0("All Age Adjusted PERMANOVA p = ", signif(canvas1_overall_p, 3))) +
    theme_void() +
    theme(
      plot.title = element_text(size = 21, face = "plain", hjust = 0.5),
      plot.margin = margin(t = 10, r = 0, b = 10, l = 0)
    )
  
  final_plot <- wrap_plots(
    permanova_plot,
    wrap_plots(canvas1_plots, ncol = 1),
    ncol = 1,
    heights = c(0.06, length(canvas1_plots))
  )
  
  pdf("C:/Users/oofordile/Desktop/BetaDiversity_PCoA_Centroid_SquareCanvas.pdf", width = 18, height = 12)
  print(final_plot)
  dev.off()
  
  emf("C:/Users/oofordile/Desktop/BetaDiversity_PCoA_Centroid_SquareCanvas.emf", width = 18, height = 12)
  print(final_plot)
  dev.off()
}

# ======================================================
# Canvas 2: Split into A/B
# ======================================================
canvas2_names <- names(canvas2_plots)
split_idx <- ceiling(length(canvas2_names)/2)
canvas2A <- canvas2_names[1:split_idx]
canvas2B <- canvas2_names[(split_idx+1):length(canvas2_names)]

# Canvas2A
pdf("C:/Users/oofordile/Desktop/BetaDiversity_PCoA_Centroid_Canvas2A.pdf", width = 22, height = 28)
print(wrap_plots(canvas2_plots[canvas2A], ncol = 1))
dev.off()
emf("C:/Users/oofordile/Desktop/BetaDiversity_PCoA_Centroid_Canvas2A.emf", width = 22, height = 28)
print(wrap_plots(canvas2_plots[canvas2A], ncol = 1))
dev.off()

# Canvas2B
pdf("C:/Users/oofordile/Desktop/BetaDiversity_PCoA_Centroid_Canvas2B.pdf", width = 22, height = 28)
print(wrap_plots(canvas2_plots[canvas2B], ncol = 1))
dev.off()
emf("C:/Users/oofordile/Desktop/BetaDiversity_PCoA_Centroid_Canvas2B.emf", width = 22, height = 28)
print(wrap_plots(canvas2_plots[canvas2B], ncol = 1))
dev.off()

cat("\nAll files saved:\n")
cat("- BetaDiversity_Overall_PERMANOVA_Summary.csv\n")
cat("- BetaDiversity_AgeStratified_PERMANOVA_Summary.csv\n")
cat("- BetaDiversity_PCoA_Centroid_SquareCanvas.pdf (Canvas 1 - compact style)\n")
cat("- BetaDiversity_PCoA_Centroid_SquareCanvas.emf (Canvas 1 - compact style)\n")
cat("- BetaDiversity_PCoA_Centroid_Canvas2A.pdf\n")
cat("- BetaDiversity_PCoA_Centroid_Canvas2A.emf\n")
cat("- BetaDiversity_PCoA_Centroid_Canvas2B.pdf\n")
cat("- BetaDiversity_PCoA_Centroid_Canvas2B.emf\n")
cat("\nCanvas 1: Compact style - comparison labels on left, transparent points\n")
cat("Canvas 2A/B: Compact style - comparison labels and PERMANOVA p on left, Kruskal-Wallis on centroid, solid points\n")
cat("\nNature format columns added to both CSV files with full statistical reporting including:\n")
cat("- F statistic with degrees of freedom\n")
cat("- p-values\n")
cat("- R² effect sizes\n")
cat("- 95% Confidence Intervals (bootstrap-estimated)\n")
cat("- Kruskal-Wallis H statistic with degrees of freedom for dispersion tests\n")
