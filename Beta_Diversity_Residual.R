# === Load Required Libraries ===
library(dplyr)
library(phyloseq)
library(vegan)
library(ggplot2)
library(patchwork)
library(cowplot)
library(devEMF)
library(boot)

# === Reproducibility: Fix seed ===
set.seed(12345)

# === Load Data ===
microbiome_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
ae_path <- "C:/Users/oofordile/Desktop/Adverse Events.csv"

merged_data <- read.csv(microbiome_path, stringsAsFactors = FALSE)
ae_data <- read.csv(ae_path, stringsAsFactors = FALSE)

# === Subset D85 Samples ===
d85_data <- merged_data %>% filter(timepoints == "D85")
d85_data$SampleID <- paste0(d85_data$Randomisation.No, "_", d85_data$timepoints)

# === Define Illness Groups ===
ae_infections <- ae_data %>% filter(AE.Infection == 1)

# Recent-Ill: 40–75d ONLY
recent_ids <- ae_infections %>% filter(Study.Day >= 40 & Study.Day <= 75) %>% distinct(Randomisation.No) %>% pull(Randomisation.No)
exclude_recent <- union(
  ae_infections %>% filter(Study.Day < 40) %>% pull(Randomisation.No),
  ae_infections %>% filter(Study.Day > 75) %>% pull(Randomisation.No)
)
recent_ill_ids <- setdiff(recent_ids, exclude_recent)

# Early-Ill: ≤35d ONLY
early_ids <- ae_infections %>% filter(Study.Day <= 35) %>% distinct(Randomisation.No) %>% pull(Randomisation.No)
exclude_early <- ae_infections %>% filter(Study.Day > 35) %>% pull(Randomisation.No)
early_ill_ids <- setdiff(early_ids, exclude_early)

# D85 Not-Ill Reference
d85_not_ill_ids <- d85_data %>% filter(Ill_Status == "Not-Ill") %>% pull(Randomisation.No) %>% unique()

# === Subset and Annotate Metadata ===
d85_data_filtered <- d85_data %>% filter(Randomisation.No %in% c(recent_ill_ids, early_ill_ids, d85_not_ill_ids)) %>%
  mutate(Group = case_when(
    Randomisation.No %in% recent_ill_ids ~ "Recent-Ill",
    Randomisation.No %in% early_ill_ids ~ "Early-Ill",
    Randomisation.No %in% d85_not_ill_ids ~ "D85 Not-Ill"
  ))
stopifnot(!any(is.na(d85_data_filtered$Group)))

metadata <- data.frame(
  SampleID = paste0(d85_data_filtered$Randomisation.No, "_", d85_data_filtered$timepoints),
  Randomisation_No = d85_data_filtered$Randomisation.No,
  Group = factor(d85_data_filtered$Group, levels = c("D85 Not-Ill", "Early-Ill", "Recent-Ill")),
  Ill_Status = factor(d85_data_filtered$Ill_Status, levels = c("Ill", "Not-Ill")),
  AgeGroup = d85_data_filtered$X3agegroups.based.on.age.at.sampling,
  Timepoint = d85_data_filtered$timepoints,
  row.names = paste0(d85_data_filtered$Randomisation.No, "_", d85_data_filtered$timepoints),
  stringsAsFactors = FALSE
)

# Harmonize AgeGroup levels
metadata$AgeGroup[metadata$AgeGroup == "7to12 mths"] <- "7 to 12 months"
metadata$AgeGroup[metadata$AgeGroup == "1to2years"] <- "1 to 2 years"
metadata$AgeGroup[metadata$AgeGroup == "plus2years"] <- "plus 2 years"
metadata$AgeGroup <- factor(metadata$AgeGroup, levels = c("7 to 12 months", "1 to 2 years", "plus 2 years"))

# Microbiome Counts
microbiome_counts <- d85_data_filtered[, 13:500]
rownames(microbiome_counts) <- metadata$SampleID
stopifnot(all(rownames(microbiome_counts) == rownames(metadata)))

# === Create Phyloseq Object ===
otu_tab <- otu_table(as.matrix(microbiome_counts), taxa_are_rows = FALSE)
physeq <- phyloseq(otu_tab, sample_data(metadata))
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# === Define Group Colors & Shapes ===
group_colors <- c(
  "D85 Not-Ill" = "#4575B4",   # blue
  "Early-Ill"   = "#F0E442",   # yellow
  "Recent-Ill"  = "#FFA500"    # orange
)
age_shapes <- c("7 to 12 months" = 22, "1 to 2 years" = 21, "plus 2 years" = 24)

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

# === Group Analysis Function ===
analyze_group <- function(group_name, physeq_obj, comp_index){
  subset_samples <- sample_names(physeq_obj)[sample_data(physeq_obj)$Group %in% c("D85 Not-Ill", group_name)]
  physeq_sub <- prune_samples(subset_samples, physeq_obj)
  
  bray_dist <- phyloseq::distance(physeq_sub, method="bray")
  sample_data_df <- data.frame(sample_data(physeq_sub))
  
  # Overall PERMANOVA
  set.seed(12345 + comp_index)
  perm_res <- vegan::adonis2(bray_dist ~ Group + AgeGroup, data=sample_data_df, permutations=999)
  
  f_stat <- perm_res$F[1]
  r2 <- perm_res$R2[1]
  p_val <- perm_res$`Pr(>F)`[1]
  df_num <- perm_res$Df[1]
  df_denom <- perm_res$Df[length(perm_res$Df) - 1]
  
  # Calculate 95% CI for R²
  message("  Calculating bootstrap CI for R² for ", group_name, "...")
  ci_r2 <- calculate_r2_ci(bray_dist, sample_data_df$Group, sample_data_df$AgeGroup, n_boot = 199)
  
  # Beta dispersion
  betadisp <- vegan::betadisper(bray_dist, sample_data_df$Group)
  kw_res <- kruskal.test(betadisp$distances ~ sample_data_df$Group)
  kruskal_stat <- kw_res$statistic
  kruskal_df <- kw_res$parameter
  disp_p <- kw_res$p.value
  
  # Age-stratified PERMANOVA with bootstrap CI
  age_perm <- lapply(levels(sample_data_df$AgeGroup), function(ag){
    idx <- sample_data_df$AgeGroup == ag
    if (sum(idx) < 2) return(data.frame(AgeGroup = ag, F = NA, R2 = NA, P = NA, 
                                        Df_num = NA, Df_denom = NA,
                                        CI_lower = NA, CI_upper = NA))
    dist_sub <- as.dist(as.matrix(bray_dist)[idx, idx])
    grp_sub <- sample_data_df$Group[idx]
    if(length(unique(grp_sub)) > 1){
      set.seed(12345 + comp_index + which(levels(sample_data_df$AgeGroup) == ag))
      ad_res <- adonis2(dist_sub ~ grp_sub, permutations=999)
      
      # Bootstrap CI for age-stratified
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
                 F = ad_res$F[1], 
                 R2 = ad_res$R2[1], 
                 P = ad_res$`Pr(>F)`[1],
                 Df_num = ad_res$Df[1],
                 Df_denom = ad_res$Df[length(ad_res$Df) - 1],
                 CI_lower = ci_age[1],
                 CI_upper = ci_age[2])
    } else {
      data.frame(AgeGroup = ag, F = NA, R2 = NA, P = NA,
                 Df_num = NA, Df_denom = NA,
                 CI_lower = NA, CI_upper = NA)
    }
  })
  age_perm_df <- bind_rows(age_perm)
  
  # PCoA ordination
  ordination <- ordinate(physeq_sub, method="PCoA", distance=bray_dist)
  pcoa_df <- plot_ordination(physeq_sub, ordination, type="samples", justDF=TRUE)
  pcoa_df$AgeGroup <- factor(pcoa_df$AgeGroup, levels = levels(metadata$AgeGroup))
  pcoa_df$Group <- factor(pcoa_df$Group, levels = c("D85 Not-Ill", group_name))
  
  list(
    bray_dist = bray_dist,
    sample_data_df = sample_data_df,
    perm_res = perm_res,
    f_stat = f_stat,
    r2 = r2,
    p_val = p_val,
    df_num = df_num,
    df_denom = df_denom,
    ci_r2 = ci_r2,
    kw_res = kw_res,
    kruskal_stat = kruskal_stat,
    kruskal_df = kruskal_df,
    disp_p = disp_p,
    age_perm_df = age_perm_df,
    ordination = ordination,
    pcoa_df = pcoa_df,
    betadisp = betadisp
  )
}

# === Run Analyses ===
message("Analysing Recent-Ill vs D85 Not-Ill")
recent_results <- analyze_group("Recent-Ill", physeq_rel, comp_index = 1)
message("Analysing Early-Ill vs D85 Not-Ill")
early_results <- analyze_group("Early-Ill", physeq_rel, comp_index = 2)

# ======================================================
# FDR CORRECTION
# ======================================================
# Collect all p-values for overall PERMANOVA tests
all_overall_pvals <- c(recent_results$p_val, early_results$p_val)
all_overall_fdr <- p.adjust(all_overall_pvals, method = "BH")

recent_results$p_val_fdr <- all_overall_fdr[1]
early_results$p_val_fdr <- all_overall_fdr[2]

# Collect all age-stratified p-values for FDR correction
all_age_pvals <- c(recent_results$age_perm_df$P, early_results$age_perm_df$P)
all_age_fdr <- p.adjust(all_age_pvals, method = "BH")

# Assign back to age_perm_df
recent_results$age_perm_df$P_FDR <- all_age_fdr[1:nrow(recent_results$age_perm_df)]
early_results$age_perm_df$P_FDR <- all_age_fdr[(nrow(recent_results$age_perm_df)+1):length(all_age_fdr)]

# Collect all dispersion p-values for FDR correction
all_disp_pvals <- c(recent_results$disp_p, early_results$disp_p)
all_disp_fdr <- p.adjust(all_disp_pvals, method = "BH")

recent_results$disp_p_fdr <- all_disp_fdr[1]
early_results$disp_p_fdr <- all_disp_fdr[2]

cat("\n=== FDR Correction Applied ===\n")
cat("Overall PERMANOVA p-values (FDR-corrected):\n")
cat("  Recent-Ill: raw p =", signif(recent_results$p_val, 3), ", FDR-adjusted p =", signif(recent_results$p_val_fdr, 3), "\n")
cat("  Early-Ill: raw p =", signif(early_results$p_val, 3), ", FDR-adjusted p =", signif(early_results$p_val_fdr, 3), "\n")

# ======================================================
# CREATE PLOTS WITH FDR-CORRECTED P-VALUES
# ======================================================
create_plots <- function(results, group_name) {
  base_font_size <- 19
  
  # Use FDR-corrected p-values in facet labels
  facet_labels <- setNames(
    paste0(results$age_perm_df$AgeGroup, " (FDR p = ", signif(results$age_perm_df$P_FDR, 3), ")"), 
    results$age_perm_df$AgeGroup
  )
  
  x_breaks <- pretty(results$pcoa_df$Axis.1, n = 5)
  y_breaks <- pretty(results$pcoa_df$Axis.2, n = 5)
  
  # PCoA plot
  pcoa_plot <- ggplot(results$pcoa_df, aes(Axis.1, Axis.2, fill=Group, shape=AgeGroup)) +
    geom_point(size=3.5, alpha=0.5, color="black", stroke = 0.3) +
    scale_fill_manual(values=group_colors, guide="none") +
    scale_shape_manual(values=age_shapes, guide="none") +
    facet_wrap(~AgeGroup, labeller=labeller(AgeGroup=facet_labels), strip.position = "top") +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks) +
    labs(subtitle=paste0(group_name, " vs D85 Not-Ill"),
         x=paste0("PCoA1 (", round(results$ordination$values$Relative_eig[1]*100,1), "%)"),
         y=paste0("PCoA2 (", round(results$ordination$values$Relative_eig[2]*100,1), "%)")) +
    theme_minimal(base_size = base_font_size) +
    theme(
      plot.subtitle = element_text(size = base_font_size + 1, face = "plain", hjust = 0, margin = margin(b = 6)),
      legend.position = "none",
      strip.text = element_text(size = base_font_size - 3, face = "plain", margin = margin(b = 4)),
      strip.background = element_rect(fill = NA, colour = NA),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.text.y = element_text(angle = 0, hjust = 0.5),
      axis.line = element_line(color = "black", size = 0.3),
      axis.ticks = element_line(color = "black", size = 0.3),
      axis.ticks.length = unit(-2, "mm"),
      plot.margin = margin(t = 2, r = 5, b = 2, l = 5),
      panel.spacing = unit(0.4, "lines")
    )
  
  # Distance to centroid plot with FDR-corrected p-value
  centroid_df <- data.frame(DistanceToCentroid = results$betadisp$distances, Group = results$sample_data_df$Group)
  centroid_df$Group <- factor(centroid_df$Group, levels=c("D85 Not-Ill", group_name))
  
  centroid_plot <- ggplot(centroid_df, aes(x=Group, y=DistanceToCentroid, fill=Group)) +
    geom_boxplot(alpha=0.8, outlier.shape=NA) +
    geom_jitter(width=0.15, alpha=0.6, color="black", size=1.8) +
    scale_fill_manual(values=group_colors, guide="none") +
    labs(y="Distance to centroid",
         x=NULL,
         subtitle=paste0("Kruskal–Wallis FDR p = ", signif(results$disp_p_fdr, 3))) +
    theme_minimal(base_size = base_font_size) +
    theme(
      legend.position = "none",
      plot.subtitle = element_text(size = base_font_size - 2, hjust = 0.5, margin = margin(t = 0, b = 4)),
      axis.line = element_line(color = "black", size = 0.3),
      axis.ticks = element_line(color = "black", size = 0.3),
      axis.ticks.length = unit(-2, "mm"),
      panel.grid = element_blank(),
      plot.margin = margin(t = 2, r = 5, b = 2, l = 5)
    )
  
  list(pcoa_plot = pcoa_plot, centroid_plot = centroid_plot)
}

recent_plots <- create_plots(recent_results, "Recent-Ill")
early_plots <- create_plots(early_results, "Early-Ill")

# === Combine Plots ===
recent_combined <- recent_plots$pcoa_plot + recent_plots$centroid_plot + 
  plot_layout(widths=c(3,1))
early_combined <- early_plots$pcoa_plot + early_plots$centroid_plot + 
  plot_layout(widths=c(3,1))

# === Create PERMANOVA text plot with FDR ===
permanova_plot <- ggplot() +
  labs(title = paste0("All Age Adjusted PERMANOVA FDR p = ", signif(recent_results$p_val_fdr, 3))) +
  theme_void() +
  theme(
    plot.title = element_text(size = 21, face = "plain", hjust = 0.5),
    plot.margin = margin(t = 10, r = 0, b = 10, l = 0)
  )

# === Combine all with compact spacing ===
final_plot <- wrap_plots(
  permanova_plot,
  wrap_plots(list(recent_combined, early_combined), ncol = 1),
  ncol = 1,
  heights = c(0.06, 2)
)

# === Save PDF & EMF (compact dimensions) ===
pdf("C:/Users/oofordile/Desktop/D85_BetaDiversity_IllComparisons_FDR.pdf", width = 18, height = 10)
print(final_plot)
dev.off()

emf("C:/Users/oofordile/Desktop/D85_BetaDiversity_IllComparisons_FDR.emf", width = 18, height = 10)
print(final_plot)
dev.off()

# === Helper function to format numbers appropriately ===
format_number <- function(x, digits = 4) {
  if (is.na(x)) return(NA)
  if (x == 0) return("0.0000")  # Explicit zero
  if (x < 0.0001) {
    return(formatC(x, format = "e", digits = 2))
  } else {
    return(sprintf(paste0("%.", digits, "f"), x))
  }
}

# === Prepare CSV Outputs with Proper Formatting ===
illness_groups <- list(
  list(name = "Recent-Ill", results = recent_results),
  list(name = "Early-Ill", results = early_results)
)

results_summary <- list()
age_stratified_results <- list()

for(ill_group in illness_groups){
  group_name <- ill_group$name
  analysis_results <- ill_group$results
  
  # Format confidence intervals
  ci_text <- if (!is.na(analysis_results$ci_r2[1]) && !is.na(analysis_results$ci_r2[2])) {
    paste0("[", format_number(analysis_results$ci_r2[1]), ", ", 
           format_number(analysis_results$ci_r2[2]), "]")
  } else {
    "[CI not available]"
  }
  
  # Format p-values for Nature format string
  p_formatted <- format_number(analysis_results$p_val)
  p_fdr_formatted <- format_number(analysis_results$p_val_fdr)
  r2_formatted <- format_number(analysis_results$r2)
  
  nature_format <- paste0("F(", analysis_results$df_num, ", ", analysis_results$df_denom, ") = ", 
                          format_number(analysis_results$f_stat, 2), ", p = ", p_formatted, 
                          ", FDR p = ", p_fdr_formatted, ", R² = ", r2_formatted, 
                          ", 95% CI ", ci_text)
  
  kruskal_p_formatted <- format_number(analysis_results$disp_p)
  kruskal_p_fdr_formatted <- format_number(analysis_results$disp_p_fdr)
  
  kruskal_format <- paste0("H(", analysis_results$kruskal_df, ") = ", 
                           format_number(analysis_results$kruskal_stat, 2), ", p = ", 
                           kruskal_p_formatted, ", FDR p = ", kruskal_p_fdr_formatted)
  
  # Overall results with formatted numbers
  results_summary[[group_name]] <- data.frame(
    Comparison = paste0(group_name, " vs D85 Not-Ill"),
    PERMANOVA_Df_numerator = analysis_results$df_num,
    PERMANOVA_Df_denominator = analysis_results$df_denom,
    PERMANOVA_SumOfSqs = format_number(analysis_results$perm_res$SumOfSqs[1]),
    PERMANOVA_R2 = format_number(analysis_results$r2),
    PERMANOVA_F = format_number(analysis_results$f_stat),
    PERMANOVA_P = format_number(analysis_results$p_val),
    PERMANOVA_P_FDR = format_number(analysis_results$p_val_fdr),
    R2_95CI_lower = format_number(analysis_results$ci_r2[1]),
    R2_95CI_upper = format_number(analysis_results$ci_r2[2]),
    `Kruskal-Wallis_H` = format_number(analysis_results$kruskal_stat),
    `Kruskal-Wallis_Df` = analysis_results$kruskal_df,
    `Kruskal-Wallis_P` = format_number(analysis_results$disp_p),
    `Kruskal-Wallis_P_FDR` = format_number(analysis_results$disp_p_fdr),
    Nature_Format_PERMANOVA = nature_format,
    Nature_Format_KruskalWallis = kruskal_format,
    stringsAsFactors = FALSE
  )
  
  # Age-stratified with formatted numbers
  age_df <- analysis_results$age_perm_df
  age_stratified_row <- data.frame(Comparison = paste0(group_name, " vs D85 Not-Ill"), 
                                   stringsAsFactors = FALSE)
  
  for (i in 1:nrow(age_df)) {
    age_group <- age_df$AgeGroup[i]
    f_val <- age_df$F[i]
    r2_val <- age_df$R2[i]
    p_val_age <- age_df$P[i]
    p_fdr_age <- age_df$P_FDR[i]
    df_n <- age_df$Df_num[i]
    df_d <- age_df$Df_denom[i]
    ci_l <- age_df$CI_lower[i]
    ci_u <- age_df$CI_upper[i]
    
    if (is.na(f_val)) {
      text_val <- "NA"
    } else {
      ci_text_age <- if (!is.na(ci_l) && !is.na(ci_u)) {
        paste0("95% CI [", format_number(ci_l), ", ", format_number(ci_u), "]")
      } else {
        "95% CI [not available]"
      }
      
      text_val <- paste0("F(", df_n, ", ", df_d, ") = ", format_number(f_val, 2), 
                         ", p = ", format_number(p_val_age),
                         ", FDR p = ", format_number(p_fdr_age),
                         ", R² = ", format_number(r2_val), 
                         ", ", ci_text_age)
    }
    
    if (age_group == "7 to 12 months") age_stratified_row$`7–12 months` <- text_val
    if (age_group == "1 to 2 years") age_stratified_row$`1–2 years` <- text_val
    if (age_group == "plus 2 years") age_stratified_row$`>2 years` <- text_val
  }
  age_stratified_results[[group_name]] <- age_stratified_row
}

# === Write CSVs ===
overall_csv <- do.call(rbind, results_summary)
rownames(overall_csv) <- NULL
write.csv(overall_csv, "C:/Users/oofordile/Desktop/PERMANOVA_Overall_FDR.csv", row.names=FALSE)

age_csv <- do.call(rbind, age_stratified_results)
rownames(age_csv) <- NULL
age_csv <- age_csv[, c("Comparison", "7–12 months", "1–2 years", ">2 years")]
write.csv(age_csv, "C:/Users/oofordile/Desktop/PERMANOVA_AgeStratified_FDR.csv", row.names=FALSE)

cat("\n=== CSV Formatting Rules Applied ===\n")
cat("Numbers formatted as:\n")
cat("- 4 decimal places for values >= 0.0001\n")
cat("- Scientific notation (e.g., 1.23e-05) for values < 0.0001\n")
cat("- No cells contain just '0'\n")
cat("\n=== Analysis Complete ===\n")
cat("All files saved:\n")
cat("- D85_BetaDiversity_IllComparisons_FDR.pdf\n")
cat("- D85_BetaDiversity_IllComparisons_FDR.emf\n")
cat("- PERMANOVA_Overall_FDR.csv\n")
cat("- PERMANOVA_AgeStratified_FDR.csv\n")
cat("\nFDR correction (Benjamini-Hochberg) applied to:\n")
cat("- Overall PERMANOVA p-values (2 comparisons)\n")
cat("- Age-stratified PERMANOVA p-values (6 comparisons)\n")
cat("- Kruskal-Wallis dispersion p-values (2 comparisons)\n")
cat("\nNature format columns include both raw and FDR-adjusted p-values.\n")
