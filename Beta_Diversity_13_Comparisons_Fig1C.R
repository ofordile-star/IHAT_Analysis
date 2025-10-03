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
library(devEMF)  # EMF export
library(grid)    # for unit() if needed

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
# Define Comparisons with Cross-Timepoint Indicator
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
age_shapes <- c("7 to 12 months" = 15, "1 to 2 years" = 16, "plus 2 years" = 17)

# ======================================================
# Initialize outputs
# ======================================================
results_summary <- list()
age_stratified_results <- list()
all_plots <- list()
canvas1_plots <- list()
canvas2_plots <- list()

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
  
  cat("Processing comparison:", comp_name, "\n")
  
  ps_subset <- subset_samples(physeq_rel, eval(filter_expr))
  ps_subset <- prune_samples(sample_sums(ps_subset) > 0, ps_subset)
  
  otu_mat_sub <- as(otu_table(ps_subset), "matrix")
  if (taxa_are_rows(ps_subset)) otu_mat_sub <- t(otu_mat_sub)
  dist_mat <- vegdist(otu_mat_sub, method = "bray")
  sample_df_sub <- as(sample_data(ps_subset), "data.frame")
  
  # Overall PERMANOVA (with random effect ONLY for paired cross-timepoint comparisons)
  set.seed(12345 + comp_index)
  use_random_effect <- (is_cross_timepoint && is_paired)
  
  if (use_random_effect) {
    # Use strata parameter to account for repeated measures (paired data only)
    adonis_res <- adonis2(as.dist(dist_mat) ~ sample_df_sub[[group_var]] + sample_df_sub$AgeGroup, 
                          permutations = 999, 
                          strata = sample_df_sub$Randomisation_No)
    cat("  - Using random effect (strata) for SubjectID in overall PERMANOVA (paired comparison)\n")
  } else {
    # No random effect for within-timepoint or unpaired comparisons
    adonis_res <- adonis2(as.dist(dist_mat) ~ sample_df_sub[[group_var]] + sample_df_sub$AgeGroup, 
                          permutations = 999)
    if (is_cross_timepoint && !is_paired) {
      cat("  - No random effect used (unpaired comparison across timepoints/groups)\n")
    }
  }
  
  # Diagnostic output for problematic comparisons
  if (!use_random_effect && is_cross_timepoint) {
    cat("  - Sample sizes by group:", table(sample_df_sub[[group_var]]), "\n")
    cat("  - Number of unique subjects:", length(unique(sample_df_sub$Randomisation_No)), "\n")
  }
  
  df <- adonis_res$Df[1]
  sum_of_sqs <- adonis_res$SumOfSqs[1]
  r2 <- adonis_res$R2[1]
  f_stat <- adonis_res$F[1]
  p_val <- adonis_res$`Pr(>F)`[1]
  
  # Age-stratified PERMANOVA 
  age_permanova <- lapply(levels(sample_df_sub$AgeGroup), function(ag){
    idx <- sample_df_sub$AgeGroup == ag
    dist_sub <- as.dist(as.matrix(dist_mat)[idx, idx])
    grp_sub <- sample_df_sub[[group_var]][idx]
    if ( length(unique(grp_sub)) > 1 ) {
      set.seed(12345 + comp_index + which(levels(sample_df_sub$AgeGroup) == ag))
      # NO strata used in age-stratified analyses
      ad_res <- adonis2(dist_sub ~ grp_sub, permutations = 999)
      data.frame(
        AgeGroup = ag,
        AgeGroup_PERMANOVA_F = ad_res$F[1],
        AgeGroup_PERMANOVA_R2 = ad_res$R2[1],
        AgeGroup_PERMANOVA_P = ad_res$`Pr(>F)`[1]
      )
    } else {
      data.frame(
        AgeGroup = ag,
        AgeGroup_PERMANOVA_F = NA,
        AgeGroup_PERMANOVA_R2 = NA,
        AgeGroup_PERMANOVA_P = NA
      )
    }
  })
  age_perm_df <- bind_rows(age_permanova)
  
  # Beta dispersion
  beta_disp <- betadisper(dist_mat, sample_df_sub[[group_var]])
  kruskal_disp <- kruskal.test(beta_disp$distances ~ sample_df_sub[[group_var]])
  disp_p <- kruskal_disp$p.value
  
  centroid_df <- data.frame(
    DistanceToCentroid = beta_disp$distances,
    Group = sample_df_sub[[group_var]],
    SampleID = rownames(sample_df_sub)
  )
  
  # PCoA
  pcoa_res <- ape::pcoa(dist_mat)
  pcoa_df <- data.frame(
    Axis.1 = pcoa_res$vectors[,1],
    Axis.2 = pcoa_res$vectors[,2],
    SampleID = rownames(sample_df_sub)
  )
  pcoa_df <- left_join(pcoa_df, sample_df_sub, by = "SampleID")
  
  x_breaks <- pretty(pcoa_df$Axis.1, n = 5)
  y_breaks <- pretty(pcoa_df$Axis.2, n = 5)
  
  # Fix colors for groups actually present
  present_groups <- unique(sample_df_sub[[group_var]])
  if (group_var == "Timepoints") {
    tp_levels <- sort(present_groups)
    current_status <- unique(sample_df_sub$Ill_Status)
    main_color <- ill_colors[current_status]
    timepoint_colors <- if (comp_name %in% c("D1 vs D15 Not-Ill", "D15 vs D85 Not-Ill", "D1 vs D85 Not-Ill")) {
      setNames(c("#377EB8", gray_color)[1:length(tp_levels)], tp_levels)
    } else {
      setNames(c(main_color, gray_color)[1:length(tp_levels)], tp_levels)
    }
  } else {
    timepoint_colors <- ill_colors[present_groups]
  }
  
  facet_labels <- setNames(
    paste0(age_perm_df$AgeGroup, " (p = ", signif(age_perm_df$AgeGroup_PERMANOVA_P, 3), ")"),
    age_perm_df$AgeGroup
  )
  
  base_font_size <- if (comp_name %in% c("D1 Ill vs D1 Not-Ill", "D15 Ill vs D15 Not-Ill", "D85 Ill vs D85 Not-Ill")) 16 else 19
  
  # --- PCoA Plot ---
  p_pcoa <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = .data[[group_var]], shape = AgeGroup)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = timepoint_colors) +
    scale_shape_manual(values = age_shapes) +
    facet_wrap(~AgeGroup, labeller = labeller(AgeGroup = facet_labels), strip.position = "top") +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks) +
    labs(title = comp_name,
         subtitle = paste0("PERMANOVA p = ", signif(p_val, 3)),
         x = paste0("PCoA1 (", round(pcoa_res$values$Relative_eig[1]*100,1), "%)"),
         y = paste0("PCoA2 (", round(pcoa_res$values$Relative_eig[2]*100,1), "%)")) +
    theme_minimal(base_size = base_font_size) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = base_font_size - 3, face = "plain", margin = margin(b = 6)),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.placement = "outside",
      plot.title = element_text(size = base_font_size + 2, face = "plain"),
      plot.subtitle = element_text(size = base_font_size + 2, face = "plain"),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.text.y = element_text(angle = 0, hjust = 0.5),
      plot.margin = margin(t = 8, r = 5, b = 5, l = 5),
      panel.spacing = unit(1, "lines"),
      axis.line.x = element_line(size = 0.3, colour = "black"),
      axis.line.y = element_line(size = 0.3, colour = "black"),
      axis.ticks = element_line(size = 0.3),
      axis.ticks.length = unit(-5, "pt")  # inward ticks
    )
  
  # --- Centroid Plot ---
  centroid_df[[group_var]] <- factor(centroid_df$Group)
  p_centroid <- ggplot(centroid_df, aes(x = Group, y = DistanceToCentroid, fill = Group)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.6, color = "black", size = 1.8) +
    scale_fill_manual(values = timepoint_colors) +
    labs(title = "Distance to Centroid",
         subtitle = paste0("Kruskal–Wallis p = ", signif(disp_p, 3)),
         y = "Distance to centroid", x = NULL) +
    theme_minimal(base_size = base_font_size) +
    theme(
      legend.position = "bottom",
      panel.grid = element_blank(),
      axis.line.x = element_line(size = 0.3, colour = "black"),
      axis.line.y = element_line(size = 0.3, colour = "black"),
      axis.ticks = element_line(size = 0.3),
      axis.ticks.length = unit(-5, "pt")  # inward ticks
    )
  
  combined <- (p_pcoa + p_centroid + patchwork::plot_layout(widths = c(3,1))) /
    patchwork::plot_spacer() + patchwork::plot_layout(heights = c(1,0.05)) & 
    theme(legend.position = "bottom")
  
  all_plots[[comp_name]] <- combined
  if (comp_name %in% c("D1 Ill vs D1 Not-Ill", "D15 Ill vs D15 Not-Ill", "D85 Ill vs D85 Not-Ill")) {
    canvas1_plots[[comp_name]] <- combined
  } else {
    canvas2_plots[[comp_name]] <- combined
  }
  
  # --- Store PERMANOVA Results ---
  results_summary[[comp_name]] <- data.frame(
    Comparison = comp_name,
    Cross_Timepoint = is_cross_timepoint,
    Paired = is_paired,
    Random_Effect_Used = (is_cross_timepoint && is_paired),
    PERMANOVA_Df = round(df, 4),
    PERMANOVA_SumOfSqs = round(sum_of_sqs, 4),
    PERMANOVA_R2 = round(r2, 4),
    PERMANOVA_F = round(f_stat, 4),
    PERMANOVA_P = round(p_val, 4),
    `Kruskal-Wallis_Dispersion_P` = round(disp_p, 4)
  )
  
  age_stratified_row <- data.frame(Comparison = comp_name, stringsAsFactors = FALSE)
  for (i in 1:nrow(age_perm_df)) {
    age_group <- age_perm_df$AgeGroup[i]
    f_val <- age_perm_df$AgeGroup_PERMANOVA_F[i]
    r2_val <- age_perm_df$AgeGroup_PERMANOVA_R2[i]
    p_val_age <- age_perm_df$AgeGroup_PERMANOVA_P[i]
    text_val <- if (is.na(f_val)) "NA, NA, NA" else paste(round(f_val,4), round(r2_val,4), round(p_val_age,4), sep = ", ")
    if (age_group == "7 to 12 months") age_stratified_row$`7–12 mo (F, R², p)` <- text_val
    if (age_group == "1 to 2 years") age_stratified_row$`1–2 yr (F, R², p)` <- text_val
    if (age_group == "plus 2 years") age_stratified_row$`>2 yr (F, R², p)` <- text_val
  }
  age_stratified_results[[comp_name]] <- age_stratified_row
}

# ======================================================
# Create CSV Files
# ======================================================
results_df <- do.call(rbind, results_summary)
rownames(results_df) <- NULL
write.csv(results_df, "C:/Users/oofordile/Desktop/BetaDiversity_Overall_PERMANOVA_Summary.csv", row.names = FALSE)

age_stratified_df <- do.call(rbind, age_stratified_results)
rownames(age_stratified_df) <- NULL
if (!"7–12 mo (F, R², p)" %in% colnames(age_stratified_df)) age_stratified_df$`7–12 mo (F, R², p)` <- "NA, NA, NA"
if (!"1–2 yr (F, R², p)" %in% colnames(age_stratified_df)) age_stratified_df$`1–2 yr (F, R², p)` <- "NA, NA, NA"
if (!">2 yr (F, R², p)" %in% colnames(age_stratified_df)) age_stratified_df$`>2 yr (F, R², p)` <- "NA, NA, NA"
age_stratified_df <- age_stratified_df[, c("Comparison", "7–12 mo (F, R², p)", "1–2 yr (F, R², p)", ">2 yr (F, R², p)")]
write.csv(age_stratified_df, "C:/Users/oofordile/Desktop/BetaDiversity_AgeStratified_PERMANOVA_Summary.csv", row.names = FALSE)

# ======================================================
# Split Canvas 2
# ======================================================
canvas2_names <- names(canvas2_plots)
first_half <- canvas2_names[1:5]
second_half <- canvas2_names[6:10]

# ======================================================
# Save PDFs
# ======================================================
pdf("C:/Users/oofordile/Desktop/BetaDiversity_PCoA_Centroid_SquareCanvas.pdf", width = 14, height = 12)
patchwork::wrap_plots(canvas1_plots, ncol = 1)
dev.off()

pdf("C:/Users/oofordile/Desktop/BetaDiversity_PCoA_Centroid_Canvas2A.pdf", width = 18, height = 27.5)
patchwork::wrap_plots(canvas2_plots[first_half], ncol = 1)
dev.off()

pdf("C:/Users/oofordile/Desktop/BetaDiversity_PCoA_Centroid_Canvas2B.pdf", width = 18, height = 27.5)
patchwork::wrap_plots(canvas2_plots[second_half], ncol = 1)
dev.off()

# ======================================================
# Save EMF Versions
# ======================================================
emf("C:/Users/oofordile/Desktop/BetaDiversity_PCoA_Centroid_SquareCanvas.emf", width = 14, height = 12)
patchwork::wrap_plots(canvas1_plots, ncol = 1)
dev.off()

emf("C:/Users/oofordile/Desktop/BetaDiversity_PCoA_Centroid_Canvas2A.emf", width = 18, height = 27.5)
patchwork::wrap_plots(canvas2_plots[first_half], ncol = 1)
dev.off()

emf("C:/Users/oofordile/Desktop/BetaDiversity_PCoA_Centroid_Canvas2B.emf", width = 18, height = 27.5)
patchwork::wrap_plots(canvas2_plots[second_half], ncol = 1)
dev.off()

# ======================================================
# Print Summary
# ======================================================
cat("Overall PERMANOVA Results:\n")
print(results_df)
cat("\nAge-Stratified PERMANOVA Results:\n")
print(age_stratified_df)
cat("\nFiles saved:\n")
cat("- BetaDiversity_Overall_PERMANOVA_Summary.csv\n")
cat("- BetaDiversity_AgeStratified_PERMANOVA_Summary.csv\n")
cat("- BetaDiversity_PCoA_Centroid_SquareCanvas.pdf/.emf\n")
cat("- BetaDiversity_PCoA_Centroid_Canvas2A.pdf/.emf\n")
cat("- BetaDiversity_PCoA_Centroid_Canvas2B.pdf/.emf\n")
cat("\nNote: Random effects (strata) applied only to paired cross-timepoint comparisons in overall PERMANOVA.\n")
cat("Paired comparisons: D1 vs D15, D15 vs D85, D1 vs D85 (within Ill or Not-Ill groups)\n")
cat("Unpaired comparisons: Different subjects at different timepoints (e.g., D1 Ill vs D85 Not-Ill)\n")
