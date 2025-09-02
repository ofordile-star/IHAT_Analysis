# ======================================================
# Load Required Libraries
# ======================================================
library(phyloseq)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(grid)      # for textGrob and gpar
library(broom)     # for tidy/anova tables
library(scales)    # for axis labels
library(lme4)      # mixed-effects models
library(lmerTest)  # p-values in mixed models
library(devEMF)    # EMF output

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
# Metrics and comparisons
# ======================================================
metrics <- c("Richness", "Shannon", "Simpson", "Pielou", "Fisher")
comparisons <- list(
  c("D1_Ill", "D1_Not-Ill"),
  c("D15_Ill", "D15_Not-Ill"),
  c("D85_Ill", "D85_Not-Ill"),
  c("D1_Ill", "D85_Not-Ill"),
  c("D15_Ill", "D85_Not-Ill"),
  c("D1_Not-Ill", "D15_Not-Ill"),
  c("D15_Not-Ill", "D85_Not-Ill"),
  c("D1_Ill", "D15_Ill"),
  c("D15_Ill", "D85_Ill"),
  c("D1_Ill", "D85_Ill"),
  c("D1_Not-Ill", "D85_Not-Ill"),
  c("D1_Not-Ill", "D85_Ill"),
  c("D15_Not-Ill", "D85_Ill")
)

# ======================================================
# GLM / Mixed Model Comparison Function
# ======================================================
run_comparison_glm <- function(data, metric, group_var) {
  groups <- unique(data[[group_var]])
  if(length(groups) != 2) return(NA_real_)
  
  timepoints_in_groups <- sapply(strsplit(as.character(groups), "_"), function(x) x[1])
  same_timepoint <- length(unique(timepoints_in_groups)) == 1
  
  if(same_timepoint) {
    fit <- lm(as.formula(paste(metric, "~ Ill_Status + AgeGroup")), data = data)
    aov_tab <- anova(fit)
    p_val <- tryCatch(as.numeric(aov_tab["Ill_Status", "Pr(>F)"]), error = function(e) NA_real_)
  } else {
    fml <- as.formula(paste0(metric, " ~ ", group_var, " + AgeGroup + (1|SubjectID)"))
    p_val <- tryCatch({
      fit <- lmer(fml, data = data, REML = FALSE)
      an_tab <- anova(fit, type = 3)
      row_idx <- which(rownames(an_tab) == group_var)
      if(length(row_idx) > 0) an_tab[row_idx, "Pr(>F)"] else NA_real_
    }, error = function(e) {
      fit <- lm(as.formula(paste(metric, "~", group_var, "+ AgeGroup")), data = data)
      aov_tab <- anova(fit)
      tryCatch(as.numeric(aov_tab[group_var, "Pr(>F)"]), error = function(e2) NA_real_)
    })
  }
  return(p_val)
}

# ======================================================
# Run all comparisons
# ======================================================
comparison_results <- list()
for(i in seq_along(comparisons)){
  comp <- comparisons[[i]]
  comp_name <- paste(comp, collapse = " vs ")
  comp_data <- alpha_data %>% filter(Group %in% comp)
  comp_data$ComparisonGroup <- factor(comp_data$Group, levels = comp)
  
  if(nrow(comp_data) > 0 && length(unique(comp_data$Group)) == 2){
    comp_pvals <- sapply(metrics, function(m) run_comparison_glm(comp_data, m, "ComparisonGroup"))
    comp_pvals[is.na(comp_pvals)] <- 1e-6
    comp_p_adj <- p.adjust(comp_pvals, method = "BH")
    
    comp_results <- data.frame(
      Comparison = comp_name,
      Metric = metrics,
      p_value = comp_pvals,
      p_adj = comp_p_adj,
      stringsAsFactors = FALSE
    )
    comparison_results[[comp_name]] <- comp_results
  }
}
all_comparison_results <- bind_rows(comparison_results)

# ======================================================
# Export CSV (long and wide)
# ======================================================
write.csv(all_comparison_results,
          "C:/Users/oofordile/Desktop/Alpha_Diversity_Comparisons_All.csv",
          row.names = FALSE)

comparison_wide <- all_comparison_results %>%
  select(Comparison, Metric, p_adj) %>%
  pivot_wider(names_from = Metric, values_from = p_adj, names_prefix = "p adj ") %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))

write.csv(comparison_wide,
          "C:/Users/oofordile/Desktop/Alpha_Diversity_Comparisons_All_Wide.csv",
          row.names = FALSE)

# ======================================================
# Plotting function (Nature-style with thin inward axis lines)
# ======================================================
create_comparison_plots_nature <- function(alpha_data, all_comparison_results, metrics, comparisons, group_colors) {
  plot_list <- list()
  y_ranges <- lapply(metrics, function(m) {
    r <- range(alpha_data[[m]], na.rm = TRUE)
    diff <- r[2] - r[1]
    c(r[1]-0.1*diff, r[2]+0.15*diff)
  })
  names(y_ranges) <- metrics
  
  for(i in seq_along(comparisons)){
    comp <- comparisons[[i]]
    comp_name <- paste(comp, collapse = " vs ")
    comp_data <- alpha_data %>% filter(Group %in% comp)
    comp_data$ComparisonGroup <- factor(comp_data$Group, levels = comp)
    
    if(nrow(comp_data) > 0){
      for(m in metrics){
        p_val <- all_comparison_results %>%
          filter(Comparison == comp_name, Metric == m) %>% pull(p_adj)
        p_label <- ifelse(length(p_val) == 0 || is.na(p_val), "", sprintf("FDR p = %.3f", p_val))
        
        fill_colors <- group_colors[sapply(comp, function(x) strsplit(x, "_")[[1]][2])]
        names(fill_colors) <- comp
        
        p <- ggplot(comp_data, aes(x = ComparisonGroup, y = .data[[m]], fill = ComparisonGroup)) +
          geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.6) +
          geom_jitter(size = 1, alpha = 0.7, width = 0.15) +
          scale_fill_manual(values = fill_colors, guide = "none") +
          scale_y_continuous(limits = y_ranges[[m]]) +
          labs(title = paste0(m, " – ", comp_name), x = NULL, y = m) +
          theme_minimal(base_size = 10, base_family = "Helvetica") +
          theme(
            panel.grid = element_blank(),
            plot.title = element_text(size = 11, face = "bold"),
            axis.text = element_text(size = 9),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.y = element_text(size = 10),
            plot.margin = margin(t=10, r=5, b=5, l=5, unit="pt"),
            axis.line = element_line(size = 0.3, colour = "black"),
            axis.ticks = element_line(size = 0.3),
            axis.ticks.length = unit(-3, "pt")  # << inward ticks
          ) +
          annotate("text", x = 1.5, y = y_ranges[[m]][2]*0.95,
                   label = p_label, size = 4, fontface = "bold")
        
        plot_list[[paste(comp_name, m, sep = "_")]] <- p
      }
    }
  }
  return(plot_list)
}

# ======================================================
# Generate plots
# ======================================================
comparison_plots <- create_comparison_plots_nature(alpha_data, all_comparison_results, metrics, comparisons, group_colors)

# ======================================================
# Create Canvas 1 and Canvas 2 PDFs + EMFs
# ======================================================
create_canvases <- function() {
  all_plots <- comparison_plots
  ncol_canvas <- 5
  nrow_canvas1 <- 7
  nrow_canvas2 <- 6
  total_plots <- length(all_plots)
  
  plots_canvas1 <- all_plots[1:(ncol_canvas*nrow_canvas1)]
  plots_canvas2 <- all_plots[(ncol_canvas*nrow_canvas1 + 1):total_plots]
  
  # --- Canvas 1 ---
  pdf("C:/Users/oofordile/Desktop/Canvas1.pdf",
      width = ncol_canvas*3, height = nrow_canvas1*2.5, useDingbats = FALSE)
  grid.arrange(grobs = plots_canvas1, ncol = ncol_canvas,
               top = textGrob("Canvas 1", gp=gpar(fontsize=14, fontface="bold")))
  dev.off()
  
  emf("C:/Users/oofordile/Desktop/Canvas1.emf",
      width = ncol_canvas*3, height = nrow_canvas1*2.5)
  grid.arrange(grobs = plots_canvas1, ncol = ncol_canvas,
               top = textGrob("Canvas 1", gp=gpar(fontsize=14, fontface="bold")))
  dev.off()
  
  # --- Canvas 2 ---
  pdf("C:/Users/oofordile/Desktop/Canvas2.pdf",
      width = ncol_canvas*3, height = nrow_canvas2*2.5, useDingbats = FALSE)
  grid.arrange(grobs = plots_canvas2, ncol = ncol_canvas,
               top = textGrob("Canvas 2", gp=gpar(fontsize=14, fontface="bold")))
  dev.off()
  
  emf("C:/Users/oofordile/Desktop/Canvas2.emf",
      width = ncol_canvas*3, height = nrow_canvas2*2.5)
  grid.arrange(grobs = plots_canvas2, ncol = ncol_canvas,
               top = textGrob("Canvas 2", gp=gpar(fontsize=14, fontface="bold")))
  dev.off()
}

# ======================================================
# Generate Canvas 1 and Canvas 2
# ======================================================
create_canvases()

print("All analyses, CSV export, and Canvas 1/Canvas 2 PDF/EMF plots complete!")
