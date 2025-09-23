# Load Libraries
library(phyloseq)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(grid)
library(broom)
library(lme4)       # mixed-effects models
library(lmerTest)   # p-values for mixed models
library(devEMF)     # EMF export

# ======================================================
# Define colors for Ill Status
# ======================================================
group_colors <- c("Ill" = "#D73027", "Not-Ill" = "#4575B4")

# ======================================================
# Load data
# ======================================================
data_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
df <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE)
df <- df %>% filter(!is.na(Ill_Status) & Ill_Status != "")
df$SampleID <- paste0(df$`Randomisation No`, "_", df$timepoints)
stopifnot(!any(duplicated(df$SampleID)))

# ======================================================
# Extract count data (columns 13:500)
# ======================================================
count_data <- df[, 13:500]
count_mat <- as.matrix(count_data)
rownames(count_mat) <- df$SampleID

# ======================================================
# Sample metadata and phyloseq object
# ======================================================
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

# ======================================================
# Estimate alpha diversity
# ======================================================
alpha_div <- estimate_richness(ps, measures = c("Observed","Shannon","Simpson","Fisher"))
alpha_div$Pielou <- alpha_div$Shannon / log(alpha_div$Observed)
alpha_div$SampleID <- df$SampleID
alpha_data <- merge(alpha_div, sample_data_df, by="SampleID")
alpha_data <- alpha_data %>% rename(Richness = Observed)
alpha_data$timepoints <- factor(alpha_data$timepoints, levels=c("D1","D15","D85"))
alpha_data$Ill_Status <- factor(alpha_data$Ill_Status, levels=c("Ill","Not-Ill"))
alpha_data$Group <- paste0(alpha_data$timepoints,"_",alpha_data$Ill_Status)

# ======================================================
# The 5 specific comparisons
# ======================================================
comparisons <- list(
  c("D1_Ill","D1_Not-Ill"),
  c("D15_Ill","D15_Not-Ill"),
  c("D85_Ill","D85_Not-Ill"),
  c("D1_Ill","D85_Not-Ill"),
  c("D15_Ill","D85_Not-Ill")
)
plot_metrics <- c("Richness","Fisher")

# ======================================================
# Order age groups (use available ones if some missing)
# ======================================================
age_groups <- c("7to12 mths","1to2years","plus2years")
existing_age_groups <- intersect(age_groups, unique(alpha_data$AgeGroup))
if(length(existing_age_groups) < length(age_groups)) {
  warning("Some age groups not found in data, using available ones")
  age_groups <- existing_age_groups
}

# ======================================================
# Mixed-model / GLM comparison function (age-stratified)
# ======================================================
run_age_stratified_model <- function(data, metric, group_var, age_group) {
  age_data <- data %>% filter(AgeGroup==age_group)
  groups <- unique(age_data[[group_var]])
  if(length(groups)!=2 || nrow(age_data)<4) return(NA_real_)
  
  timepoints_in_groups <- sapply(strsplit(as.character(groups),"_"), function(x) x[1])
  same_timepoint <- length(unique(timepoints_in_groups))==1
  
  if(same_timepoint) {
    fit <- lm(as.formula(paste(metric,"~",group_var)), data=age_data)
    aov_tab <- anova(fit)
    p_val <- tryCatch(as.numeric(aov_tab[group_var,"Pr(>F)"]), error=function(e) NA_real_)
  } else {
    fml <- as.formula(paste0(metric," ~ ",group_var," + (1|SubjectID)"))
    p_val <- tryCatch({
      fit <- lmer(fml, data=age_data, REML=FALSE)
      an_tab <- anova(fit, type=3)
      row_idx <- which(rownames(an_tab)==group_var)
      if(length(row_idx)>0) an_tab[row_idx,"Pr(>F)"] else NA_real_
    }, error=function(e){
      fit <- lm(as.formula(paste(metric,"~",group_var)), data=age_data)
      aov_tab <- anova(fit)
      tryCatch(as.numeric(aov_tab[group_var,"Pr(>F)"]), error=function(e2) NA_real_)
    })
  }
  return(p_val)
}

# ======================================================
# Compute age-stratified p-values and export CSV
# ======================================================
age_stratified_results <- list()
for(comp in comparisons){
  comp_name <- paste(comp, collapse=" vs ")
  comp_data <- alpha_data %>% filter(Group %in% comp)
  comp_data$ComparisonGroup <- factor(comp_data$Group, levels=comp)
  if(nrow(comp_data)>0 && length(unique(comp_data$Group))==2){
    for(age_grp in age_groups){
      comp_pvals <- sapply(plot_metrics, function(m) run_age_stratified_model(comp_data, m, "ComparisonGroup", age_grp))
      comp_pvals[is.na(comp_pvals)] <- 1e-6
      comp_p_adj <- p.adjust(comp_pvals, method="BH")
      comp_results <- data.frame(
        Comparison = comp_name,
        AgeGroup = age_grp,
        Metric = plot_metrics,
        p_value = comp_pvals,
        p_adj = comp_p_adj,
        stringsAsFactors=FALSE
      )
      age_stratified_results[[paste(comp_name,age_grp,sep="_")]] <- comp_results
    }
  }
}
all_age_stratified_results <- bind_rows(age_stratified_results)
write.csv(all_age_stratified_results,
          "C:/Users/oofordile/Desktop/Alpha_Diversity_5Comparisons_AgeStratified_MixedModels.csv",
          row.names=FALSE)

# ======================================================
# Prepare plotting: compute global y ranges (for annotation placement)
# ======================================================
y_ranges <- lapply(plot_metrics, function(m){
  r <- range(alpha_data[[m]], na.rm=TRUE)
  if(any(is.infinite(r))) r <- c(0,1)
  diff <- r[2]-r[1]
  if(diff == 0) diff <- abs(r[2])*0.1 + 1e-6
  c(r[1]-0.1*diff, r[2]+0.15*diff)
})
names(y_ranges) <- plot_metrics

empty_plot <- function() {
  ggplot() + theme_void()
}

# ======================================================
# Font sizes (reduced by 1 point from previous version)
# ======================================================
base_font_size <- 10          
title_sub_size <- base_font_size
p_font_size <- title_sub_size / 2.5

# ======================================================
# Create plots for each age group, metric, and comparison
# ======================================================
plots_store <- list()
for(age_grp in age_groups){
  for(m in plot_metrics){
    for(comp in comparisons){
      comp_name <- paste(comp, collapse=" vs ")
      comp_data <- alpha_data %>% filter(Group %in% comp, AgeGroup==age_grp)
      comp_data$ComparisonGroup <- factor(comp_data$Group, levels=comp)
      if(nrow(comp_data)>0 && length(unique(comp_data$Group))==2){
        p_val <- all_age_stratified_results %>% filter(Comparison==comp_name, Metric==m, AgeGroup==age_grp) %>% pull(p_adj)
        p_label <- ifelse(length(p_val)==0 || is.na(p_val),"",sprintf("FDR p=%.3f",p_val))
        fill_colors <- group_colors[sapply(comp, function(x) strsplit(x,"_")[[1]][2])]
        names(fill_colors) <- comp
        
        p <- ggplot(comp_data, aes(x=ComparisonGroup, y=.data[[m]], fill=ComparisonGroup)) +
          geom_boxplot(alpha=0.8, outlier.shape=NA, width=0.6) +
          geom_jitter(size=1, alpha=0.7, width=0.15) +
          scale_fill_manual(values=fill_colors, guide="none") +
          labs(title = paste("Age:", age_grp),
               subtitle = comp_name,
               x = NULL, y = m) +
          theme_minimal(base_size = base_font_size) +
          theme(
            panel.grid = element_blank(),
            axis.text.x = element_text(angle=45, hjust=1, size = base_font_size - 1),
            axis.line = element_line(size=0.3, colour="black"),
            axis.ticks = element_line(size=0.3, colour="black"),
            axis.ticks.length = unit(-2, "mm"),  # Negative value makes ticks point inward
            plot.title = element_text(face="bold", size = title_sub_size + 1, vjust = 1),
            plot.subtitle = element_text(face="bold", size = title_sub_size, vjust = 0.8),
            plot.caption = element_text(face="bold", size = base_font_size - 1),
            plot.margin = margin(t = 15, r = 5, b = 5, l = 5, unit = "pt")
          ) +
          expand_limits(y = y_ranges[[m]][2] * 1.05) +
          annotate("text", x = 1.5, 
                   y = y_ranges[[m]][2] * 0.92, 
                   label = p_label, 
                   fontface = "bold", 
                   size = p_font_size)
        plots_store[[paste(age_grp,m,comp_name,sep="_")]] <- p
      } else {
        plots_store[[paste(age_grp,m,comp_name,sep="_")]] <- empty_plot()
      }
    }
  }
}

# ======================================================
# Build and export two canvases: one for Richness, one for Fisher
# ======================================================
n_cols <- length(comparisons)
n_rows_age <- length(age_groups)

collect_metric_grobs <- function(metric_name) {
  grobs <- list()
  for(age_grp in age_groups){
    for(comp in comparisons){
      comp_name <- paste(comp, collapse=" vs ")
      g <- plots_store[[paste(age_grp,metric_name,comp_name,sep="_")]]
      if(is.null(g)) g <- empty_plot()
      grobs <- c(grobs, list(g))
    }
  }
  grobs
}

top_title_rich <- paste0("Age-Stratified Richness (rows = age groups: ",
                         paste(age_groups, collapse=", "),
                         "; columns = comparisons)")
top_title_fisher <- paste0("Age-Stratified Fisher (rows = age groups: ",
                           paste(age_groups, collapse=", "),
                           "; columns = comparisons)")

# ---------- Richness canvas ----------
rich_grobs <- collect_metric_grobs("Richness")
pdf("C:/Users/oofordile/Desktop/Canvas_Richness_5Comparisons.pdf",
    width = n_cols * 3, height = n_rows_age * 2.5, useDingbats = FALSE)
grid.arrange(grobs = rich_grobs, ncol = n_cols,
             top = textGrob(top_title_rich, gp = gpar(fontsize = base_font_size + 4, fontface = "bold")))
dev.off()

emf("C:/Users/oofordile/Desktop/Canvas_Richness_5Comparisons.emf",
    width = n_cols * 3, height = n_rows_age * 2.5, emfPlus = TRUE)
grid.arrange(grobs = rich_grobs, ncol = n_cols,
             top = textGrob(top_title_rich, gp = gpar(fontsize = base_font_size + 4, fontface = "bold")))
dev.off()

# ---------- Fisher canvas ----------
fisher_grobs <- collect_metric_grobs("Fisher")
pdf("C:/Users/oofordile/Desktop/Canvas_Fisher_5Comparisons.pdf",
    width = n_cols * 3, height = n_rows_age * 2.5, useDingbats = FALSE)
grid.arrange(grobs = fisher_grobs, ncol = n_cols,
             top = textGrob(top_title_fisher, gp = gpar(fontsize = base_font_size + 4, fontface = "bold")))
dev.off()

emf("C:/Users/oofordile/Desktop/Canvas_Fisher_5Comparisons.emf",
    width = n_cols * 3, height = n_rows_age * 2.5, emfPlus = TRUE)
grid.arrange(grobs = fisher_grobs, ncol = n_cols,
             top = textGrob(top_title_fisher, gp = gpar(fontsize = base_font_size + 4, fontface = "bold")))
dev.off()

cat("\nCanvases exported:\n - C:/Users/oofordile/Desktop/Canvas_Richness_5Comparisons.pdf/.emf\n - C:/Users/oofordile/Desktop/Canvas_Fisher_5Comparisons.pdf/.emf\nCSV of age-stratified p-values saved to Desktop.\n")
