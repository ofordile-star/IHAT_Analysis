# ======================================================
# Load Required Libraries
# ======================================================
library(phyloseq)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(grid)
library(broom)
library(devEMF)  # For EMF export
library(lme4)    # For mixed-effects models
library(lmerTest) # For p-values in mixed models

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
# Extract count data
# ======================================================
count_data <- df[, 13:500]
stopifnot(all(sapply(count_data, is.numeric)))
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
  SubjectID = df$`Randomisation No`
)
rownames(sample_data_df) <- df$SampleID
sample_data_ps <- sample_data(sample_data_df)
otu_table_ps <- otu_table(count_mat, taxa_are_rows = FALSE)
ps <- phyloseq(otu_table_ps, sample_data_ps)

# ======================================================
# Estimate alpha diversity
# ======================================================
alpha_div <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson", "Fisher"))
alpha_div$Pielou <- alpha_div$Shannon / log(alpha_div$Observed)
alpha_div$SampleID <- df$SampleID
alpha_data <- merge(alpha_div, sample_data_df, by = "SampleID")
alpha_data <- alpha_data %>% rename(Richness = Observed)
alpha_data$timepoints <- factor(alpha_data$timepoints, levels = c("D1", "D15", "D85"))
alpha_data$Ill_Status <- factor(alpha_data$Ill_Status, levels = c("Ill", "Not-Ill"))
alpha_data$Group <- paste0(alpha_data$timepoints, "_", alpha_data$Ill_Status)

# ======================================================
# Define the 5 specific comparisons
# ======================================================
comparisons <- list(
  c("D1_Ill", "D1_Not-Ill"),
  c("D15_Ill", "D15_Not-Ill"),
  c("D85_Ill", "D85_Not-Ill"),
  c("D1_Ill", "D85_Not-Ill"),
  c("D15_Ill", "D85_Not-Ill")
)

plot_metrics <- c("Richness", "Fisher")
all_metrics <- c("Richness", "Shannon", "Simpson", "Pielou", "Fisher")

# ======================================================
# GLM comparison function with mixed-effects fallback
# ======================================================
run_comparison_glm <- function(data, metric, group_var) {
  groups <- unique(data[[group_var]])
  if(length(groups) != 2) return(NA)
  timepoints_in_groups <- sapply(strsplit(as.character(groups), "_"), function(x) x[1])
  same_timepoint <- length(unique(timepoints_in_groups)) == 1
  
  if(same_timepoint) {
    fit <- lm(reformulate(c("Ill_Status","AgeGroup"), metric), data = data)
    an <- broom::tidy(anova(fit))
    p_val <- an$p.value[an$term == "Ill_Status"]
  } else {
    p_val <- tryCatch({
      fit <- lmer(reformulate(c(group_var,"AgeGroup","(1|SubjectID)"), metric), data = data)
      coef_table <- summary(fit)$coefficients
      if(nrow(coef_table) >= 2) coef_table[2,"Pr(>|t|)"] else NA
    }, error=function(e) {
      fit <- lm(reformulate(c(group_var,"AgeGroup"), metric), data = data)
      an <- broom::tidy(anova(fit))
      an$p.value[an$term == group_var]
    })
  }
  if(is.null(p_val) || length(p_val)==0 || is.na(p_val)) return(NA)
  return(p_val)
}

# ======================================================
# Run comparisons
# ======================================================
overall_results <- list()
for(comp in comparisons){
  comp_name <- paste(comp, collapse = " vs ")
  comp_data <- alpha_data %>% filter(Group %in% comp)
  comp_data$ComparisonGroup <- factor(comp_data$Group, levels = comp)
  if(nrow(comp_data) > 0 && length(unique(comp_data$Group))==2){
    comp_pvals <- sapply(all_metrics, function(m) run_comparison_glm(comp_data, m, "ComparisonGroup"))
    comp_pvals[is.na(comp_pvals)] <- 1e-6
    comp_p_adj <- p.adjust(comp_pvals, method="BH")
    comp_results <- data.frame(
      Comparison = comp_name,
      Metric = all_metrics,
      p_value = comp_pvals,
      p_adj = comp_p_adj,
      Analysis_Type="Overall",
      stringsAsFactors=FALSE
    )
    overall_results[[comp_name]] <- comp_results
  }
}
all_overall_results <- bind_rows(overall_results)
write.csv(all_overall_results,
          "C:/Users/oofordile/Desktop/Alpha Diversity 5 Comparisons Results.csv",
          row.names = FALSE)

# ======================================================
# Plotting function with thin axis lines + inward ticks
# ======================================================
plot_canvas_style <- function(plot_obj, y_limits, p_label, subtitle=NULL, xaxis_font_size=NULL, increase_font=2) {
  theme_x <- if(!is.null(xaxis_font_size)) element_text(size=xaxis_font_size+increase_font, angle=45, hjust=1) else element_text(size=9+increase_font, angle=45, hjust=1)
  
  plot_obj +
    scale_y_continuous(limits=y_limits) +
    theme_minimal(base_size=10+increase_font, base_family="Helvetica") +
    theme(
      panel.grid=element_blank(),
      plot.title=element_text(size=13+increase_font, face="bold"),
      plot.subtitle=element_text(size=10+increase_font, face="bold"),
      axis.text.y=element_text(size=10+increase_font),
      axis.text.x=theme_x,
      axis.title.y=element_text(size=11+increase_font),
      plot.margin=margin(t=10,r=5,b=5,l=5,unit="pt"),
      axis.line = element_line(size=0.3, colour="black"),
      axis.ticks = element_line(size=0.3, colour="black"),
      axis.ticks.length = unit(-2, "mm")  # NEGATIVE for inward ticks
    ) +
    annotate("text", x=1.5, y=y_limits[2]*0.95,
             label=p_label, size=5+increase_font*0.5, fontface="bold") +
    labs(subtitle=subtitle)
}

# ======================================================
# Generate overall plots
# ======================================================
create_overall_plots <- function(alpha_data, all_overall_results, plot_metrics, comparisons, group_colors) {
  plot_list <- list()
  y_ranges <- lapply(plot_metrics, function(m) {
    r <- range(alpha_data[[m]], na.rm=TRUE)
    diff <- r[2]-r[1]
    c(r[1]-0.1*diff, r[2]+0.15*diff)
  })
  names(y_ranges) <- plot_metrics
  for(comp in comparisons){
    comp_name <- paste(comp, collapse=" vs ")
    comp_data <- alpha_data %>% filter(Group %in% comp)
    comp_data$ComparisonGroup <- factor(comp_data$Group, levels=comp)
    if(nrow(comp_data) > 0){
      for(m in plot_metrics){
        p_val <- all_overall_results %>% filter(Comparison==comp_name, Metric==m) %>% pull(p_adj)
        p_label <- ifelse(length(p_val)==0 || is.na(p_val), "", sprintf("FDR p = %.3f", p_val))
        fill_colors <- group_colors[sapply(comp, function(x) strsplit(x,"_")[[1]][2])]
        names(fill_colors) <- comp
        p <- ggplot(comp_data, aes(x=ComparisonGroup, y=.data[[m]], fill=ComparisonGroup)) +
          geom_boxplot(alpha=0.8, outlier.shape=NA, width=0.6) +
          geom_jitter(size=1, alpha=0.7, width=0.15) +
          scale_fill_manual(values=fill_colors, guide="none") +
          labs(title=paste(m,"-",comp_name), x=NULL, y=m)
        plot_list[[paste(m,comp_name,sep="_")]] <- plot_canvas_style(p, y_ranges[[m]], p_label, xaxis_font_size=10.5, increase_font=2)
      }
    }
  }
  plot_list
}

overall_plots <- create_overall_plots(alpha_data, all_overall_results, plot_metrics, comparisons, group_colors)

# ======================================================
# Split into two canvases: Richness and Fisher
# ======================================================
richness_plots <- lapply(comparisons, function(comp){
  overall_plots[[paste("Richness", paste(comp, collapse=" vs "), sep="_")]]
})
fisher_plots <- lapply(comparisons, function(comp){
  overall_plots[[paste("Fisher", paste(comp, collapse=" vs "), sep="_")]]
})

# ---------- Richness canvas ----------
pdf("C:/Users/oofordile/Desktop/Alpha Diversity_Richness_5Comparisons.pdf", width=20, height=4, useDingbats=FALSE)
grid.arrange(grobs=richness_plots, ncol=5, nrow=1,
             top=textGrob("Richness: Overall Analysis (GLM, Age-Adjusted, FDR)",
                          gp=gpar(fontsize=16, fontface="bold")))
dev.off()

emf("C:/Users/oofordile/Desktop/Alpha Diversity_Richness_5Comparisons.emf", width=20, height=4, emfPlus=TRUE)
grid.arrange(grobs=richness_plots, ncol=5, nrow=1,
             top=textGrob("Richness: Overall Analysis (GLM, Age-Adjusted, FDR)",
                          gp=gpar(fontsize=16, fontface="bold")))
dev.off()

# ---------- Fisher canvas ----------
pdf("C:/Users/oofordile/Desktop/Alpha Diversity_Fisher_5Comparisons.pdf", width=20, height=4, useDingbats=FALSE)
grid.arrange(grobs=fisher_plots, ncol=5, nrow=1,
             top=textGrob("Fisher: Overall Analysis (GLM, Age-Adjusted, FDR)",
                          gp=gpar(fontsize=16, fontface="bold")))
dev.off()

emf("C:/Users/oofordile/Desktop/Alpha Diversity_Fisher_5Comparisons.emf", width=20, height=4, emfPlus=TRUE)
grid.arrange(grobs=fisher_plots, ncol=5, nrow=1,
             top=textGrob("Fisher: Overall Analysis (GLM, Age-Adjusted, FDR)",
                          gp=gpar(fontsize=16, fontface="bold")))
dev.off()

cat("\nAnalysis complete! Separate Richness and Fisher PDFs + EMFs exported to Desktop.\n")
