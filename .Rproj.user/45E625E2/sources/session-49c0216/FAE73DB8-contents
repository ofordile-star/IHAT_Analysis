# === Load Required Libraries ===
library(dplyr)
library(phyloseq)
library(vegan)
library(ggplot2)
library(patchwork)
library(cowplot)
library(devEMF)

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
  "D85 Not-Ill" = "#4575B4", "Early-Ill" = "#F0E442", "Recent-Ill" = "#D73027"
)
age_shapes <- c("7 to 12 months" = 15, "1 to 2 years" = 16, "plus 2 years" = 17)

# === Group Analysis Function ===
analyze_group <- function(group_name, physeq_obj){
  subset_samples <- sample_names(physeq_obj)[sample_data(physeq_obj)$Group %in% c("D85 Not-Ill", group_name)]
  physeq_sub <- prune_samples(subset_samples, physeq_obj)
  
  bray_dist <- phyloseq::distance(physeq_sub, method="bray")
  sample_data_df <- data.frame(sample_data(physeq_sub))
  
  set.seed(42)
  perm_res <- vegan::adonis2(bray_dist ~ Group + AgeGroup, data=sample_data_df, permutations=999)
  
  betadisp <- vegan::betadisper(bray_dist, sample_data_df$Group)
  kw_res <- kruskal.test(betadisp$distances ~ sample_data_df$Group)
  
  # PERMANOVA per AgeGroup
  age_perm <- lapply(levels(sample_data_df$AgeGroup), function(ag){
    idx <- sample_data_df$AgeGroup == ag
    dist_sub <- as.dist(as.matrix(bray_dist)[idx, idx])
    grp_sub <- sample_data_df$Group[idx]
    if(length(unique(grp_sub))>1){
      ad_res <- adonis2(dist_sub ~ grp_sub, permutations=999)
      data.frame(AgeGroup=ag, F=ad_res$F[1], R2=ad_res$R2[1], P=ad_res$`Pr(>F)`[1])
    } else data.frame(AgeGroup=ag, F=NA, R2=NA, P=NA)
  })
  age_perm_df <- bind_rows(age_perm)
  
  ordination <- ordinate(physeq_sub, method="PCoA", distance=bray_dist)
  pcoa_df <- plot_ordination(physeq_sub, ordination, type="samples", justDF=TRUE)
  pcoa_df$AgeGroup <- factor(pcoa_df$AgeGroup, levels = levels(metadata$AgeGroup))
  pcoa_df$Group <- factor(pcoa_df$Group, levels = c("D85 Not-Ill", group_name))
  
  base_font_size <- 19
  facet_labels <- setNames(
    paste0(age_perm_df$AgeGroup, " (p = ", signif(age_perm_df$P, 3), ")"), 
    age_perm_df$AgeGroup
  )
  
  x_breaks <- pretty(pcoa_df$Axis.1, n = 5)
  y_breaks <- pretty(pcoa_df$Axis.2, n = 5)
  
  pcoa_plot <- ggplot(pcoa_df, aes(Axis.1, Axis.2, color=Group, shape=AgeGroup)) +
    geom_point(size=3, alpha=0.7) +
    scale_color_manual(values=group_colors) +
    scale_shape_manual(values=age_shapes) +
    facet_wrap(~AgeGroup, labeller=labeller(AgeGroup=facet_labels), strip.position = "top") +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks) +
    labs(title=paste0("D85 Not-Ill vs ", group_name),
         subtitle=paste0("PERMANOVA p = ", signif(perm_res$`Pr(>F)`[1], 3)),
         x=paste0("PCoA1 (", round(ordination$values$Relative_eig[1]*100,1), "%)"),
         y=paste0("PCoA2 (", round(ordination$values$Relative_eig[2]*100,1), "%)")) +
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
      axis.line = element_line(size = 0.3, colour = "black"),
      axis.ticks = element_line(size = 0.3),
      axis.ticks.length = unit(-5, "pt"), # << inward ticks
      plot.margin = margin(t = 8, r = 5, b = 5, l = 5),
      panel.spacing = unit(1, "lines")
    )
  
  centroid_df <- data.frame(DistanceToCentroid = betadisp$distances, Group = sample_data_df$Group)
  centroid_df$Group <- factor(centroid_df$Group, levels=c("D85 Not-Ill", group_name))
  
  centroid_plot <- ggplot(centroid_df, aes(x=Group, y=DistanceToCentroid, fill=Group)) +
    geom_boxplot(alpha=0.8, outlier.shape=NA) +
    geom_jitter(width=0.15, alpha=0.6, color="black", size=1.8) +
    scale_fill_manual(values=group_colors) +
    labs(title="Distance to Centroid",
         subtitle=paste0("Kruskal–Wallis p = ", signif(kw_res$p.value, 3)),
         y="Distance to centroid", x=NULL) +
    theme_minimal(base_size = base_font_size) +
    theme(
      legend.position = "bottom",
      panel.grid = element_blank(),
      axis.line = element_line(size = 0.3, colour = "black"),
      axis.ticks = element_line(size = 0.3),
      axis.ticks.length = unit(-5, "pt") # << inward ticks
    )
  
  list(pcoa_plot=pcoa_plot, centroid_plot=centroid_plot,
       perm_results=perm_res, kw_disp=kw_res, age_perm_df=age_perm_df)
}

# === Run Analyses ===
recent_results <- analyze_group("Recent-Ill", physeq_rel)
early_results <- analyze_group("Early-Ill", physeq_rel)

# === Combine Plots ===
recent_combined <- (recent_results$pcoa_plot + recent_results$centroid_plot + patchwork::plot_layout(widths = c(3,1))) / 
  patchwork::plot_spacer() + patchwork::plot_layout(heights = c(1,0.05)) & 
  theme(legend.position = "bottom")

early_combined <- (early_results$pcoa_plot + early_results$centroid_plot + patchwork::plot_layout(widths = c(3,1))) / 
  patchwork::plot_spacer() + patchwork::plot_layout(heights = c(1,0.05)) & 
  theme(legend.position = "bottom")

final_plot <- patchwork::wrap_plots(list(recent_combined, early_combined), ncol = 1)

# === Save PDF & EMF ===
pdf("C:/Users/oofordile/Desktop/D85_BetaDiversity_IllComparisons_3facets_RecentTop.pdf", width = 19, height = 12)
print(final_plot)
dev.off()

emf("C:/Users/oofordile/Desktop/D85_BetaDiversity_IllComparisons_3facets_RecentTop.emf", width = 19, height = 12)
print(final_plot)
dev.off()

# === Prepare CSV Outputs ===
illness_groups <- c("Recent-Ill","Early-Ill")
age_csv <- data.frame()
overall_csv <- data.frame()

for(ill_group in illness_groups){
  if(ill_group=="Recent-Ill") analysis_results <- recent_results else analysis_results <- early_results
  perm <- analysis_results$perm_results
  disp_p <- analysis_results$kw_disp$p.value
  age_df <- analysis_results$age_perm_df
  
  # Not age-stratified CSV
  overall_csv <- rbind(overall_csv, data.frame(
    Comparison = paste0(ill_group," vs D85 Not-Ill"),
    PERMANOVA_Df = perm$Df[1],
    PERMANOVA_SumOfSqs = round(perm$SumOfSqs[1],4),
    PERMANOVA_R2 = round(perm$R2[1],4),
    PERMANOVA_F = round(perm$F[1],4),
    PERMANOVA_P = round(perm$`Pr(>F)`[1],4),
    Kruskal_Wallis = round(disp_p,4)
  ))
  
  # Age-stratified CSV
  age_csv <- rbind(age_csv, data.frame(
    Comparison = paste0(ill_group," vs D85 Not-Ill"),
    `7–12 mo (F,R2,p)` = paste0(round(age_df$F[age_df$AgeGroup=="7 to 12 months"],4),",",
                                round(age_df$R2[age_df$AgeGroup=="7 to 12 months"],4),",",
                                round(age_df$P[age_df$AgeGroup=="7 to 12 months"],4)),
    `1–2 yr (F,R2,p)` = paste0(round(age_df$F[age_df$AgeGroup=="1 to 2 years"],4),",",
                               round(age_df$R2[age_df$AgeGroup=="1 to 2 years"],4),",",
                               round(age_df$P[age_df$AgeGroup=="1 to 2 years"],4)),
    `>2 yr (F,R2,p)` = paste0(round(age_df$F[age_df$AgeGroup=="plus 2 years"],4),",",
                              round(age_df$R2[age_df$AgeGroup=="plus 2 years"],4),",",
                              round(age_df$P[age_df$AgeGroup=="plus 2 years"],4))
  ))
}

# === Write CSVs ===
write.csv(overall_csv, "C:/Users/oofordile/Desktop/PERMANOVA_Overall_4decimal.csv", row.names=FALSE)
write.csv(age_csv, "C:/Users/oofordile/Desktop/PERMANOVA_AgeStratified_4decimal.csv", row.names=FALSE)

cat("CSV files exported:\n- PERMANOVA_Overall_4decimal.csv\n- PERMANOVA_AgeStratified_4decimal.csv\n")
