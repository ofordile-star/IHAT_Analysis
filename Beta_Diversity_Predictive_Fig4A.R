# ======================================================
# Load Required Libraries
# ======================================================
library(dplyr)
library(phyloseq)
library(vegan)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(devEMF)
library(grid)
library(cowplot)
library(boot)

# ======================================================
# File Paths (will be overwritten)
# ======================================================
microbiome_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
ae_path <- "C:/Users/oofordile/Desktop/Adverse Events.csv"
csv_out_nonage <- "C:/Users/oofordile/Desktop/Age_Adjusted_PERMANOVA_Results.csv"
csv_out_age <- "C:/Users/oofordile/Desktop/Age_Stratified_PERMANOVA_Results.csv"
pdf_out <- "C:/Users/oofordile/Desktop/BetaDiversity_AE_PCoA_Centroid_Compact.pdf"
emf_out <- "C:/Users/oofordile/Desktop/BetaDiversity_AE_PCoA_Centroid_Compact.emf"

# ======================================================
# Reproducibility: Fix seed
# ======================================================
set.seed(12345)

# ======================================================
# Load Data
# ======================================================
merged_data <- read.csv(microbiome_path, stringsAsFactors = FALSE)
ae_data <- read.csv(ae_path, stringsAsFactors = FALSE)

# ======================================================
# Identify earliest AE with Infection=1 per child
# ======================================================
ae_infections <- ae_data %>% filter(AE.Infection == 1)
earliest_ae <- ae_infections %>%
  group_by(Randomisation.No) %>%
  summarise(EarliestDay = min(Study.Day, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(IllnessGroup = case_when(
    !is.na(EarliestDay) & EarliestDay <= 35 ~ "Soon-Ill",
    !is.na(EarliestDay) & EarliestDay >= 36 & EarliestDay <= 70 ~ "Later-Ill",
    !is.na(EarliestDay) & EarliestDay > 70 ~ "Much-Later-Ill",
    TRUE ~ NA_character_
  ))

# ======================================================
# Filter D1 and D85
# ======================================================
d1_data <- merged_data %>% filter(timepoints == "D1")
d85_data <- merged_data %>% filter(timepoints == "D85")
d85_not_ill <- d85_data %>% filter(Ill_Status == "Not-Ill")

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
# Colors and Shapes
# ======================================================
illness_groups <- c("Soon-Ill", "Later-Ill", "Much-Later-Ill")
group_colors <- c(
  "D85 Not-Ill"    = "#2166AC",
  "Soon-Ill"       = "#FF8000",
  "Later-Ill"      = "#FFD700",
  "Much-Later-Ill" = "#33A02C"
)
age_shapes <- c("7 to 12 months" = 22, "1 to 2 years" = 21, "plus 2 years" = 24)

all_plots <- list()

# Results containers
results_summary <- list()
age_stratified_results <- list()

# ======================================================
# Loop Through Illness Groups
# ======================================================
comp_index <- 0
for (ill_group in illness_groups) {
  comp_index <- comp_index + 1
  
  # Prepare sample sets
  d1_ill_samples <- earliest_ae %>% filter(IllnessGroup == ill_group) %>% pull(Randomisation.No)
  d1_ill_data <- d1_data %>% filter(Randomisation.No %in% d1_ill_samples) %>% mutate(Group = ill_group)
  d85_ref_data <- d85_not_ill %>% mutate(Group = "D85 Not-Ill")
  combined <- bind_rows(d1_ill_data, d85_ref_data)
  
  # if combined is empty, skip
  if (nrow(combined) == 0) {
    warning(paste0("No samples for ", ill_group, " vs D85 Not-Ill — skipping."))
    next
  }
  
  combined$Group <- factor(combined$Group, levels = c("D85 Not-Ill", ill_group))
  
  # Metadata
  metadata <- data.frame(
    SampleID = paste0(combined$Randomisation.No, "_", combined$timepoints),
    Randomisation_No = combined$Randomisation.No,
    Group = combined$Group,
    Ill_Status = factor(combined$Ill_Status, levels = c("Ill","Not-Ill")),
    AgeGroup = combined$X3agegroups.based.on.age.at.sampling,
    Timepoint = combined$timepoints,
    row.names = paste0(combined$Randomisation.No, "_", combined$timepoints),
    stringsAsFactors = FALSE
  )
  
  # Normalize AgeGroup labels
  metadata$AgeGroup <- as.character(metadata$AgeGroup)
  metadata$AgeGroup[metadata$AgeGroup %in% c("7to12 mths","7to12mths","7to12","7_to_12","7to12 months")] <- "7 to 12 months"
  metadata$AgeGroup[metadata$AgeGroup %in% c("1to2years","1to2yrs","1_to_2","1to2")] <- "1 to 2 years"
  metadata$AgeGroup[metadata$AgeGroup %in% c("plus2years","plus 2 years","2+ years","2plus")] <- "plus 2 years"
  metadata$AgeGroup <- factor(metadata$AgeGroup, levels=c("7 to 12 months","1 to 2 years","plus 2 years"))
  
  # Get counts
  if (ncol(combined) < 13) stop("Data has fewer than 13 columns; adjust the counts column indices.")
  counts <- combined[,13:min(500, ncol(combined))]
  counts <- as.data.frame(lapply(counts, function(col) as.numeric(as.character(col))))
  rownames(counts) <- metadata$SampleID
  if (!all(rownames(counts) == rownames(metadata))) stop("Row names of counts and metadata do not match.")
  
  # phyloseq object
  otu_tab <- otu_table(as.matrix(counts), taxa_are_rows=FALSE)
  physeq <- phyloseq(otu_tab, sample_data(metadata))
  physeq_rel <- transform_sample_counts(physeq, function(x) if(sum(x, na.rm=TRUE)==0) x else x/sum(x, na.rm=TRUE))
  
  # Bray distance
  bray_dist <- phyloseq::distance(physeq_rel, method="bray")
  meta_df <- data.frame(sample_data(physeq_rel), stringsAsFactors = FALSE)
  
  # PERMANOVA (age-adjusted)
  message("Analysing ", ill_group, " vs D85 Not-Ill")
  set.seed(12345 + comp_index)
  
  adonis_res <- adonis2(as.dist(bray_dist) ~ meta_df$Group + meta_df$AgeGroup, permutations = 999)
  
  f_stat <- adonis_res$F[1]
  r2 <- adonis_res$R2[1]
  p_val <- adonis_res$`Pr(>F)`[1]
  df_num <- adonis_res$Df[1]
  df_denom <- adonis_res$Df[length(adonis_res$Df) - 1]
  
  # Calculate 95% CI for R² (OPTIMIZED: 199 iterations)
  message("  Calculating bootstrap CI for R²...")
  ci_r2 <- calculate_r2_ci(bray_dist, meta_df$Group, meta_df$AgeGroup, n_boot = 199)
  
  # Age-stratified PERMANOVA (OPTIMIZED)
  age_permanova <- lapply(levels(meta_df$AgeGroup), function(ag){
    idx <- meta_df$AgeGroup == ag
    if (sum(idx) < 2) return(data.frame(AgeGroup = ag, AgeGroup_PERMANOVA_F = NA, AgeGroup_PERMANOVA_R2 = NA, 
                                        AgeGroup_PERMANOVA_P = NA, AgeGroup_Df_num = NA, AgeGroup_Df_denom = NA,
                                        AgeGroup_CI_lower = NA, AgeGroup_CI_upper = NA))
    dist_sub <- as.dist(as.matrix(bray_dist)[idx, idx])
    grp_sub <- meta_df$Group[idx]
    if (length(unique(grp_sub)) > 1) {
      set.seed(12345 + comp_index + which(levels(meta_df$AgeGroup) == ag))
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
  beta_disp <- betadisper(bray_dist, meta_df$Group)
  kruskal_disp <- kruskal.test(beta_disp$distances ~ meta_df$Group)
  disp_p <- kruskal_disp$p.value
  kruskal_stat <- kruskal_disp$statistic
  kruskal_df <- kruskal_disp$parameter
  
  centroid_df <- data.frame(DistanceToCentroid = beta_disp$distances,
                            Group = meta_df$Group,
                            SampleID = rownames(meta_df))
  
  # PCoA
  pcoa_res <- ape::pcoa(bray_dist)
  pcoa_df <- data.frame(Axis.1 = pcoa_res$vectors[,1],
                        Axis.2 = pcoa_res$vectors[,2],
                        SampleID = rownames(meta_df))
  pcoa_df <- left_join(pcoa_df, meta_df, by = "SampleID")
  set.seed(12345 + comp_index)
  pcoa_df <- pcoa_df[sample(seq_len(nrow(pcoa_df))), , drop = FALSE]
  x_breaks <- pretty(pcoa_df$Axis.1, n = 5)
  y_breaks <- pretty(pcoa_df$Axis.2, n = 5)
  
  # Facet labels with p-values
  facet_labels <- setNames(paste0(age_perm_df$AgeGroup, " (p = ", signif(age_perm_df$AgeGroup_PERMANOVA_P, 3), ")"), age_perm_df$AgeGroup)
  
  base_font_size <- 19
  
  # --- PCoA Plot ---
  p_pcoa <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, fill = Group, shape = AgeGroup)) +
    geom_point(size = 3.5, alpha = 0.5, color = "black", stroke = 0.3) +
    scale_fill_manual(values = group_colors, na.value = "gray50") +
    scale_shape_manual(values = age_shapes) +
    guides(fill = "none", shape = "none") +
    facet_wrap(~AgeGroup, labeller = labeller(AgeGroup = facet_labels), strip.position = "top") +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks) +
    labs(subtitle = paste0(ill_group, " vs D85 Not-Ill"),
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
  
  # --- Centroid Plot ---
  centroid_df$Group <- factor(centroid_df$Group)
  
  p_centroid <- ggplot(centroid_df, aes(x = Group, y = DistanceToCentroid, fill = Group)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.6, color = "black", size = 1.8) +
    scale_fill_manual(values = group_colors, na.value = "gray50") +
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
  
  combined <- p_pcoa + p_centroid + plot_layout(widths = c(3,1))
  all_plots[[ill_group]] <- combined
  
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
  results_summary[[ill_group]] <- data.frame(
    Comparison = paste0(ill_group, " vs D85 Not-Ill"),
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
  
  # Age-stratified with Nature format
  age_stratified_row <- data.frame(Comparison = paste0(ill_group, " vs D85 Not-Ill"), stringsAsFactors = FALSE)
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
  age_stratified_results[[ill_group]] <- age_stratified_row
}

# ======================================================
# Export results CSV
# ======================================================
results_df <- do.call(rbind, results_summary)
rownames(results_df) <- NULL
write.csv(results_df, csv_out_nonage, row.names = FALSE)

age_stratified_df <- do.call(rbind, age_stratified_results)
rownames(age_stratified_df) <- NULL
age_stratified_df <- age_stratified_df[, c("Comparison", "7–12 months", "1–2 years", ">2 years")]
write.csv(age_stratified_df, csv_out_age, row.names = FALSE)

# ======================================================
# PERMANOVA text as plot title
# ======================================================
top_p_text <- if(nrow(results_df) > 0 && !is.na(results_df$PERMANOVA_P[1])) {
  paste0("All Age Adjusted PERMANOVA p = ", signif(results_df$PERMANOVA_P[1],3))
} else {
  "All Age Adjusted PERMANOVA p = NA"
}
permanova_plot <- ggplot() +
  labs(title = top_p_text) +
  theme_void() +
  theme(
    plot.title = element_text(size = 21, face = "plain", hjust = 0.5),
    plot.margin = margin(t = 10, r = 0, b = 10, l = 0)
  )

# ======================================================
# Combine and export
# ======================================================
final_plot <- wrap_plots(
  permanova_plot,
  wrap_plots(all_plots, ncol = 1),
  ncol = 1,
  heights = c(0.06, max(1, length(all_plots)))
)

pdf(pdf_out, width=18, height=12)
print(final_plot)
dev.off()

emf(emf_out, width=18, height=12)
print(final_plot)
dev.off()

cat("\nFiles saved:\n")
cat("- ", csv_out_nonage, "\n")
cat("- ", csv_out_age, "\n")
cat("- ", pdf_out, "\n")
cat("- ", emf_out, "\n")
cat("\nNature format columns added to both CSV files with full statistical reporting including:\n")
cat("- F statistic with degrees of freedom\n")
cat("- p-values\n")
cat("- R² effect sizes\n")
cat("- 95% Confidence Intervals (bootstrap-estimated)\n")
cat("- Kruskal-Wallis H statistic with degrees of freedom for dispersion tests\n")
