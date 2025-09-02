# ======================================================
# Load Required Libraries
# ======================================================
library(dplyr)
library(phyloseq)
library(vegan)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(devEMF)  # for EMF vector export
library(grid)    # for unit()

# ======================================================
# File Paths
# ======================================================
microbiome_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
ae_path <- "C:/Users/oofordile/Desktop/Adverse Events.csv"
csv_out_nonage <- "C:/Users/oofordile/Desktop/BetaDiversity_AE_PERMANOVA_NonAge_Summary.csv"
csv_out_age <- "C:/Users/oofordile/Desktop/BetaDiversity_AE_PERMANOVA_PerAge_Summary.csv"
pdf_out <- "C:/Users/oofordile/Desktop/BetaDiversity_AE_PCoA_Centroid_Comparisons_PerAge.pdf"
emf_out <- "C:/Users/oofordile/Desktop/BetaDiversity_AE_PCoA_Centroid_Comparisons_PerAge.emf"

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
  summarise(EarliestDay = min(Study.Day)) %>%
  ungroup() %>%
  mutate(IllnessGroup = case_when(
    EarliestDay <= 35 ~ "Soon-Ill",
    EarliestDay >= 36 & EarliestDay <= 70 ~ "Later-Ill",
    EarliestDay > 70 ~ "Much-Later-Ill",
    TRUE ~ NA_character_
  ))

# ======================================================
# Filter D1 and D85
# ======================================================
d1_data <- merged_data %>% filter(timepoints == "D1")
d85_data <- merged_data %>% filter(timepoints == "D85")
d85_not_ill <- d85_data %>% filter(Ill_Status == "Not-Ill")

# ======================================================
# Colors and Shapes
# ======================================================
illness_groups <- c("Soon-Ill", "Later-Ill", "Much-Later-Ill")
group_colors <- c(
  "D85 Not-Ill"    = "#4575B4",  # blue
  "Soon-Ill"       = "#FFA500",  # orange
  "Later-Ill"      = "#F0E442",  # yellow
  "Much-Later-Ill" = "#3CA63C"   # green
)
age_shapes <- c("7 to 12 months" = 15, "1 to 2 years" = 16, "plus 2 years" = 17)

all_plots <- list()
results_summary_nonage <- list()
results_summary_age <- list()

# ======================================================
# Set seed for reproducibility
# ======================================================
set.seed(123)

# ======================================================
# Loop Through Illness Groups
# ======================================================
for (ill_group in illness_groups) {
  
  # D1 samples for illness group
  d1_ill_samples <- earliest_ae %>% filter(IllnessGroup == ill_group) %>% pull(Randomisation.No)
  d1_ill_data <- d1_data %>% filter(Randomisation.No %in% d1_ill_samples) %>% mutate(Group = ill_group)
  
  # D85 Not-Ill reference
  d85_ref_data <- d85_not_ill %>% mutate(Group = "D85 Not-Ill")
  
  # Combine
  combined <- bind_rows(d1_ill_data, d85_ref_data)
  combined$Group <- factor(combined$Group, levels = c("D85 Not-Ill", ill_group))
  
  # --- Metadata ---
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
  
  # Harmonize AgeGroup
  metadata$AgeGroup[metadata$AgeGroup=="7to12 mths"] <- "7 to 12 months"
  metadata$AgeGroup[metadata$AgeGroup=="1to2years"] <- "1 to 2 years"
  metadata$AgeGroup[metadata$AgeGroup=="plus2years"] <- "plus 2 years"
  metadata$AgeGroup <- factor(metadata$AgeGroup, levels=c("7 to 12 months","1 to 2 years","plus 2 years"))
  
  # --- OTU Table ---
  counts <- combined[,13:500]
  rownames(counts) <- metadata$SampleID
  stopifnot(all(rownames(counts) == rownames(metadata)))
  otu_tab <- otu_table(as.matrix(counts), taxa_are_rows=FALSE)
  physeq <- phyloseq(otu_tab, sample_data(metadata))
  physeq_rel <- transform_sample_counts(physeq, function(x) x/sum(x))
  
  # --- Bray-Curtis distance ---
  bray_dist <- phyloseq::distance(physeq_rel, method="bray")
  meta_df <- data.frame(sample_data(physeq_rel))
  
  # --- Overall PERMANOVA ---
  adonis_res <- adonis2(as.dist(bray_dist) ~ Group + AgeGroup, data=meta_df, permutations=999)
  df <- adonis_res$Df[1]
  sum_of_sqs <- adonis_res$SumOfSqs[1]
  r2 <- adonis_res$R2[1]
  f_stat <- adonis_res$F[1]
  p_val <- adonis_res$`Pr(>F)`[1]
  
  # --- Robust Age-stratified PERMANOVA ---
  age_perm_df <- lapply(levels(meta_df$AgeGroup), function(ag) {
    idx <- which(meta_df$AgeGroup == ag)
    if(length(unique(meta_df$Group[idx])) <= 1){
      return(data.frame(AgeGroup = ag, F = NA, R2 = NA, P = NA))
    }
    dist_sub <- as.dist(as.matrix(bray_dist)[idx, idx])
    meta_sub <- data.frame(Group = meta_df$Group[idx])
    ad_res <- adonis2(dist_sub ~ Group, data = meta_sub, permutations = 999)
    data.frame(AgeGroup = ag, F = ad_res$F[1], R2 = ad_res$R2[1], P = ad_res$`Pr(>F)`[1])
  })
  age_perm_df <- do.call(rbind, age_perm_df)
  
  # --- Beta dispersion ---
  beta_disp <- betadisper(as.dist(bray_dist), meta_df$Group)
  kruskal_disp <- kruskal.test(beta_disp$distances ~ meta_df$Group)
  disp_p <- kruskal_disp$p.value
  centroid_df <- data.frame(DistanceToCentroid=beta_disp$distances,
                            Group=meta_df$Group,
                            SampleID=rownames(meta_df))
  
  # --- PCoA ---
  pcoa_res <- cmdscale(as.dist(bray_dist), eig=TRUE, k=2)
  var_exp1 <- round(100*pcoa_res$eig[1]/sum(abs(pcoa_res$eig)),1)
  var_exp2 <- round(100*pcoa_res$eig[2]/sum(abs(pcoa_res$eig)),1)
  pcoa_df <- data.frame(SampleID=rownames(meta_df),
                        PC1=pcoa_res$points[,1],
                        PC2=pcoa_res$points[,2],
                        Group=meta_df$Group,
                        AgeGroup=meta_df$AgeGroup)
  
  x_breaks <- pretty(pcoa_df$PC1, n=5)
  y_breaks <- pretty(pcoa_df$PC2, n=5)
  facet_labels <- setNames(
    paste0(age_perm_df$AgeGroup, " (p = ", signif(age_perm_df$P, 3), ")"),
    age_perm_df$AgeGroup
  )
  
  base_font_size <- 19
  
  # --- PCoA Plot ---
  p_pcoa <- ggplot(pcoa_df, aes(x=PC1, y=PC2, color=Group, shape=AgeGroup)) +
    geom_point(size=3, alpha=0.7) +
    scale_color_manual(values=group_colors) +
    scale_shape_manual(values=age_shapes) +
    facet_wrap(~AgeGroup, labeller=labeller(AgeGroup=facet_labels), strip.position = "top") +
    scale_x_continuous(breaks=x_breaks) +
    scale_y_continuous(breaks=y_breaks) +
    labs(title=paste0(ill_group," vs D85 Not-Ill"),
         subtitle=paste0("PERMANOVA p = ", signif(p_val, 3)),
         x=paste0("PCoA1 (",var_exp1,"%)"),
         y=paste0("PCoA2 (",var_exp2,"%)")) +
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
      axis.line = element_line(color = "black", size = 0.3),
      axis.ticks = element_line(color = "black", size = 0.3),
      axis.ticks.length = unit(-2, "mm"),
      panel.border = element_blank()
    )
  
  # --- Centroid Plot ---
  centroid_df$Group <- factor(centroid_df$Group)
  p_centroid <- ggplot(centroid_df, aes(x=Group, y=DistanceToCentroid, fill=Group)) +
    geom_boxplot(alpha=0.8, outlier.shape=NA) +
    geom_jitter(width=0.15, alpha=0.6, color="black", size=1.8) +
    scale_fill_manual(values=group_colors) +
    labs(title="Distance to Centroid",
         subtitle=paste0("Kruskal–Wallis p = ", signif(disp_p, 3)),
         y="Distance to centroid", x=NULL) +
    theme_minimal(base_size = base_font_size) +
    theme(
      legend.position = "bottom",
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 0.3),
      axis.ticks = element_line(color = "black", size = 0.3),
      axis.ticks.length = unit(-2, "mm"),
      panel.border = element_blank()
    )
  
  # --- Combine PCoA + centroid ---
  combined_plot <- (p_pcoa + p_centroid + patchwork::plot_layout(widths=c(3,1))) / 
    patchwork::plot_spacer() + patchwork::plot_layout(heights=c(1,0.05)) & 
    theme(legend.position="bottom")
  
  all_plots[[ill_group]] <- combined_plot
  
  # --- Non-age-stratified results ---
  results_summary_nonage[[ill_group]] <- data.frame(
    Comparison = paste0(ill_group," vs D85 Not-Ill"),
    PERMANOVA_Df = df,
    PERMANOVA_SumOfSqs = sum_of_sqs,
    PERMANOVA_R2 = r2,
    PERMANOVA_F = f_stat,
    PERMANOVA_P = p_val,
    `Kruskal-Wallis` = disp_p,
    stringsAsFactors = FALSE
  )
  
  # --- Age-stratified results ---
  results_summary_age[[ill_group]] <- data.frame(
    Comparison = paste0(ill_group," vs D85 Not-Ill"),
    `7-12 mo (F, R², p)` = paste(round(age_perm_df$F[age_perm_df$AgeGroup=="7 to 12 months"],4),
                                 round(age_perm_df$R2[age_perm_df$AgeGroup=="7 to 12 months"],4),
                                 round(age_perm_df$P[age_perm_df$AgeGroup=="7 to 12 months"],4),
                                 sep=", "),
    `1-2 yr (F, R², p)` = paste(round(age_perm_df$F[age_perm_df$AgeGroup=="1 to 2 years"],4),
                                round(age_perm_df$R2[age_perm_df$AgeGroup=="1 to 2 years"],4),
                                round(age_perm_df$P[age_perm_df$AgeGroup=="1 to 2 years"],4),
                                sep=", "),
    `>2 yr (F, R², p)` = paste(round(age_perm_df$F[age_perm_df$AgeGroup=="plus 2 years"],4),
                               round(age_perm_df$R2[age_perm_df$AgeGroup=="plus 2 years"],4),
                               round(age_perm_df$P[age_perm_df$AgeGroup=="plus 2 years"],4),
                               sep=", "),
    stringsAsFactors = FALSE
  )
}

# ======================================================
# Export CSVs
# ======================================================
write.csv(do.call(rbind, results_summary_nonage), csv_out_nonage, row.names=FALSE)
write.csv(do.call(rbind, results_summary_age), csv_out_age, row.names=FALSE)

# ======================================================
# Export PDF
# ======================================================
pdf(pdf_out, width=18, height=16)
patchwork::wrap_plots(all_plots, ncol=1)
dev.off()

# ======================================================
# Export EMF
# ======================================================
emf(emf_out, width=18, height=16)
patchwork::wrap_plots(all_plots, ncol=1)
dev.off()

# ======================================================
# Print Summaries
# ======================================================
print(do.call(rbind, results_summary_nonage))
print(do.call(rbind, results_summary_age))

_AE_PCoA_Centroid_Comparisons_PerAge.pdf/.emf\n")
