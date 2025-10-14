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
library(effectsize) # for effect size calculations

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
# Enhanced GLM / Mixed Model Comparison Function
# ======================================================
run_comparison_detailed <- function(data, metric, group_var) {
  groups <- unique(data[[group_var]])
  if(length(groups) != 2) return(NULL)
  
  timepoints_in_groups <- sapply(strsplit(as.character(groups), "_"), function(x) x[1])
  same_timepoint <- length(unique(timepoints_in_groups)) == 1
  n <- nrow(data)
  
  if(same_timepoint) {
    # Linear model for same timepoint comparisons
    fit <- lm(as.formula(paste(metric, "~ Ill_Status + AgeGroup")), data = data)
    aov_tab <- anova(fit)
    
    p_val <- tryCatch(as.numeric(aov_tab["Ill_Status", "Pr(>F)"]), error = function(e) NA_real_)
    f_stat <- tryCatch(as.numeric(aov_tab["Ill_Status", "F value"]), error = function(e) NA_real_)
    df1 <- tryCatch(as.numeric(aov_tab["Ill_Status", "Df"]), error = function(e) NA_real_)
    df2 <- tryCatch(as.numeric(aov_tab["Residuals", "Df"]), error = function(e) NA_real_)
    
    # Estimate & CI
    coef_summary <- summary(fit)$coefficients
    if("Ill_StatusNot-Ill" %in% rownames(coef_summary)) {
      estimate <- coef_summary["Ill_StatusNot-Ill", "Estimate"]
      ci <- confint(fit)["Ill_StatusNot-Ill", ]
      ci_lower <- ci[1]
      ci_upper <- ci[2]
    } else {
      estimate <- ci_lower <- ci_upper <- NA_real_
    }
    
    # Partial eta squared
    ss_effect <- aov_tab["Ill_Status", "Sum Sq"]
    partial_eta2 <- ss_effect / (ss_effect + aov_tab["Residuals", "Sum Sq"])
    
  } else {
    # Mixed model for different timepoint comparisons
    fml <- as.formula(paste0(metric, " ~ ", group_var, " + AgeGroup + (1|SubjectID)"))
    result <- tryCatch({
      fit <- lmer(fml, data = data, REML = FALSE)
      an_tab <- anova(fit, type = 3)
      row_idx <- which(rownames(an_tab) == group_var)
      
      if(length(row_idx) > 0){
        p_val <- an_tab[row_idx, "Pr(>F)"]
        f_stat <- an_tab[row_idx, "F value"]
        df1 <- an_tab[row_idx, "NumDF"]
        df2 <- an_tab[row_idx, "DenDF"]
      } else { p_val <- f_stat <- df1 <- df2 <- NA_real_ }
      
      coef_summary <- summary(fit)$coefficients
      comparison_var_name <- paste0(group_var, levels(data[[group_var]])[2])
      if(comparison_var_name %in% rownames(coef_summary)) {
        estimate <- coef_summary[comparison_var_name, "Estimate"]
        ci <- confint(fit, parm = comparison_var_name, method = "Wald")
        ci_lower <- ci[1]
        ci_upper <- ci[2]
      } else { estimate <- ci_lower <- ci_upper <- NA_real_ }
      
      partial_eta2 <- f_stat * df1 / (f_stat * df1 + df2)
      list(p_val=p_val, f_stat=f_stat, df1=df1, df2=df2,
           estimate=estimate, ci_lower=ci_lower, ci_upper=ci_upper, partial_eta2=partial_eta2)
      
    }, error = function(e){
      fit <- lm(as.formula(paste(metric, "~", group_var, "+ AgeGroup")), data = data)
      aov_tab <- anova(fit)
      
      p_val <- tryCatch(as.numeric(aov_tab[group_var, "Pr(>F)"]), error = function(e2) NA_real_)
      f_stat <- tryCatch(as.numeric(aov_tab[group_var, "F value"]), error = function(e2) NA_real_)
      df1 <- tryCatch(as.numeric(aov_tab[group_var, "Df"]), error = function(e2) NA_real_)
      df2 <- tryCatch(as.numeric(aov_tab["Residuals", "Df"]), error = function(e2) NA_real_)
      
      coef_summary <- summary(fit)$coefficients
      comparison_var_name <- paste0(group_var, levels(data[[group_var]])[2])
      if(comparison_var_name %in% rownames(coef_summary)) {
        estimate <- coef_summary[comparison_var_name, "Estimate"]
        ci <- confint(fit)[comparison_var_name, ]
        ci_lower <- ci[1]
        ci_upper <- ci[2]
      } else { estimate <- ci_lower <- ci_upper <- NA_real_ }
      
      ss_effect <- aov_tab[group_var, "Sum Sq"]
      partial_eta2 <- ss_effect / (ss_effect + aov_tab["Residuals", "Sum Sq"])
      
      list(p_val=p_val, f_stat=f_stat, df1=df1, df2=df2,
           estimate=estimate, ci_lower=ci_lower, ci_upper=ci_upper, partial_eta2=partial_eta2)
    })
    
    p_val <- result$p_val
    f_stat <- result$f_stat
    df1 <- result$df1
    df2 <- result$df2
    estimate <- result$estimate
    ci_lower <- result$ci_lower
    ci_upper <- result$ci_upper
    partial_eta2 <- result$partial_eta2
  }
  
  return(data.frame(
    n=n, Estimate=estimate, CI_lower=ci_lower, CI_upper=ci_upper,
    F=f_stat, df1=df1, df2=df2, partial_eta2=partial_eta2, p_value=p_val
  ))
}

# ======================================================
# Run all comparisons with detailed statistics
# ======================================================
detailed_results <- list()

for(i in seq_along(comparisons)){
  comp <- comparisons[[i]]
  comp_name <- paste(comp, collapse = " vs ")
  comp_data <- alpha_data %>% filter(Group %in% comp)
  comp_data$ComparisonGroup <- factor(comp_data$Group, levels = comp)
  
  timepoint <- unique(sapply(strsplit(comp, "_"), function(x) x[1]))
  timepoint_label <- if(length(timepoint)==1) timepoint else paste(timepoint, collapse=" vs ")
  
  if(nrow(comp_data) > 0 && length(unique(comp_data$Group))==2){
    for(m in metrics){
      stats <- run_comparison_detailed(comp_data, m, "ComparisonGroup")
      if(!is.null(stats)){
        result_row <- data.frame(Comparison=comp_name, Timepoint=timepoint_label, Metric=m, stringsAsFactors=FALSE)
        result_row <- cbind(result_row, stats)
        detailed_results[[paste(comp_name, m, sep="_")]] <- result_row
      }
    }
  }
}

all_detailed_results <- bind_rows(detailed_results)

# ======================================================
# Apply FDR correction
# ======================================================
all_detailed_results <- all_detailed_results %>%
  group_by(Comparison) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  select(Comparison, Timepoint, Metric, p_value, p_adj, n, Estimate,
         CI_lower, CI_upper, F, df1, df2, partial_eta2)

# ======================================================
# Add Nature-style reporting column in requested format
# ======================================================
all_detailed_results <- all_detailed_results %>%
  mutate(Nature_Report = paste0(
    Metric, " (", Timepoint, ", Age-Adjusted): ",
    "F(", df1, ", ", df2, ") = ", round(F, 2),
    ", p = ", signif(p_adj, 3),
    ", eta2_p = ", round(partial_eta2, 3),
    ", 95% CI [", round(CI_lower, 3), ", ", round(CI_upper, 3), "]; n = ", n
  ))

# ======================================================
# Export detailed CSV
# ======================================================
write.csv(all_detailed_results,
          "C:/Users/oofordile/Desktop/Alpha_Diversity_Detailed_Statistics_Nature.csv",
          row.names = FALSE)

# ======================================================
# Simplified outputs
# ======================================================
comparison_results <- list()
for(i in seq_along(comparisons)){
  comp <- comparisons[[i]]
  comp_name <- paste(comp, collapse = " vs ")
  comp_subset <- all_detailed_results %>% filter(Comparison==comp_name)
  if(nrow(comp_subset)>0){
    comp_results <- data.frame(Comparison=comp_name, Metric=comp_subset$Metric,
                               p_value=comp_subset$p_value, p_adj=comp_subset$p_adj, stringsAsFactors=FALSE)
    comparison_results[[comp_name]] <- comp_results
  }
}
all_comparison_results <- bind_rows(comparison_results)

write.csv(all_comparison_results,
          "C:/Users/oofordile/Desktop/Alpha_Diversity_Comparisons_All.csv",
          row.names = FALSE)

comparison_wide <- all_comparison_results %>%
  select(Comparison, Metric, p_adj) %>%
  pivot_wider(names_from=Metric, values_from=p_adj, names_prefix="p adj ") %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))

write.csv(comparison_wide,
          "C:/Users/oofordile/Desktop/Alpha_Diversity_Comparisons_All_Wide.csv",
          row.names = FALSE)

# ======================================================
# Plotting function (Nature-style)
# ======================================================
create_comparison_plots_nature <- function(alpha_data, all_comparison_results, metrics, comparisons, group_colors) {
  plot_list <- list()
  y_ranges <- lapply(metrics, function(m){
    r <- range(alpha_data[[m]], na.rm=TRUE)
    diff <- r[2]-r[1]
    c(r[1]-0.1*diff, r[2]+0.15*diff)
  })
  names(y_ranges) <- metrics
  
  for(i in seq_along(comparisons)){
    comp <- comparisons[[i]]
    comp_name <- paste(comp, collapse=" vs ")
    comp_data <- alpha_data %>% filter(Group %in% comp)
    comp_data$ComparisonGroup <- factor(comp_data$Group, levels=comp)
    
    if(nrow(comp_data)>0){
      for(m in metrics){
        p_val <- all_comparison_results %>% filter(Comparison==comp_name, Metric==m) %>% pull(p_adj)
        p_label <- ifelse(length(p_val)==0 || is.na(p_val), "", sprintf("FDR p = %.3f", p_val))
        
        fill_colors <- group_colors[sapply(comp, function(x) strsplit(x,"_")[[1]][2])]
        names(fill_colors) <- comp
        
        p <- ggplot(comp_data, aes(x=ComparisonGroup, y=.data[[m]], fill=ComparisonGroup)) +
          geom_boxplot(alpha=0.8, outlier.shape=NA, width=0.6) +
          geom_jitter(size=1, alpha=0.7, width=0.15) +
          scale_fill_manual(values=fill_colors, guide="none") +
          scale_y_continuous(limits=y_ranges[[m]]) +
          labs(title=paste0(m," – ",comp_name), x=NULL, y=m) +
          theme_minimal(base_size=10, base_family="Helvetica") +
          theme(panel.grid=element_blank(),
                plot.title=element_text(size=11, face="bold"),
                axis.text=element_text(size=9),
                axis.text.x=element_text(angle=45, hjust=1),
                axis.title.y=element_text(size=10),
                plot.margin=margin(t=10,r=5,b=5,l=5,unit="pt"),
                axis.line=element_line(size=0.3, colour="black"),
                axis.ticks=element_line(size=0.3),
                axis.ticks.length=unit(-3,"pt")) +
          annotate("text", x=1.5, y=y_ranges[[m]][2]*0.95, label=p_label, size=4, fontface="bold")
        
        plot_list[[paste(comp_name,m,sep="_")]] <- p
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
  
  pdf("C:/Users/oofordile/Desktop/Canvas1.pdf", width=ncol_canvas*3, height=nrow_canvas1*2.5, useDingbats=FALSE)
  grid.arrange(grobs=plots_canvas1, ncol=ncol_canvas, top=textGrob("Canvas 1", gp=gpar(fontsize=14, fontface="bold")))
  dev.off()
  
  emf("C:/Users/oofordile/Desktop/Canvas1.emf", width=ncol_canvas*3, height=nrow_canvas1*2.5)
  grid.arrange(grobs=plots_canvas1, ncol=ncol_canvas, top=textGrob("Canvas 1", gp=gpar(fontsize=14, fontface="bold")))
  dev.off()
  
  pdf("C:/Users/oofordile/Desktop/Canvas2.pdf", width=ncol_canvas*3, height=nrow_canvas2*2.5, useDingbats=FALSE)
  grid.arrange(grobs=plots_canvas2, ncol=ncol_canvas, top=textGrob("Canvas 2", gp=gpar(fontsize=14, fontface="bold")))
  dev.off()
  
  emf("C:/Users/oofordile/Desktop/Canvas2.emf", width=ncol_canvas*3, height=nrow_canvas2*2.5)
  grid.arrange(grobs=plots_canvas2, ncol=ncol_canvas, top=textGrob("Canvas 2", gp=gpar(fontsize=14, fontface="bold")))
  dev.off()
}

# ======================================================
# Generate Canvas 1 and Canvas 2
# ======================================================
create_canvases()

print("All analyses complete!")
print(paste("Detailed statistics exported to: Alpha_Diversity_Detailed_Statistics_Nature.csv"))
print(paste("Total comparisons analyzed:", nrow(all_detailed_results)/length(metrics)))
