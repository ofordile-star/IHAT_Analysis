# ======================================================
# Load required libraries
# ======================================================
library(readr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(cluster)
library(factoextra)
library(grid)
library(devEMF)  # for EMF export

# ---- Load the dataset ----
df <- read.csv("C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv", check.names = FALSE)

# ---- Subset microbiome data (columns 13:500) ----
microbiome_data <- df[, 13:500]

# Coerce to numeric matrix safely
microbiome_matrix <- microbiome_data %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) %>%
  replace(is.na(.), 0) %>%
  as.matrix()

# Relative abundances (protect against zero row sums)
rs <- rowSums(microbiome_matrix)
rs[!is.finite(rs) | rs == 0] <- 1e-8
rel_abund <- sweep(microbiome_matrix, 1, rs, "/")

# ---- Pick top taxa robustly ----
mean_abundances <- colMeans(rel_abund) * 100
candidates <- names(mean_abundances[mean_abundances >= 0.2])
if (length(candidates) == 0) candidates <- names(sort(mean_abundances, decreasing = TRUE))
top_taxa <- head(candidates[order(mean_abundances[candidates], decreasing = TRUE)], 50)
top_taxa <- top_taxa[!is.na(top_taxa)]
if (length(top_taxa) < 2) stop("Not enough taxa after filtering to build the heatmap.")

# ---- Subset and transform for correlation ----
top_matrix <- microbiome_matrix[, top_taxa, drop = FALSE]
log_top_matrix <- log1p(top_matrix)

# ---- Correlation matrix ----
cor_matrix <- cor(log_top_matrix, method = "spearman")

# ---- Clustering ----
set.seed(123)
dist_mat <- as.dist(1 - cor_matrix)
gap_stat <- clusGap(as.matrix(dist_mat), FUN = pam, K.max = 10, B = 100)
optimal_k <- which.max(gap_stat$Tab[, "gap"])
pam_fit <- pam(dist_mat, k = optimal_k)
cluster_assignments <- pam_fit$clustering
cluster_colors <- structure(rainbow(optimal_k), names = as.character(1:optimal_k))

# ---- Ill and Not-Ill group masks ----
ill_samples    <- !is.na(df$Ill_Status) & df$Ill_Status == "Ill"
not_ill_samples <- !is.na(df$Ill_Status) & df$Ill_Status %in% c("Not-Ill", "NotIll", "Not Ill")
n_ill <- sum(ill_samples)
n_not_ill <- sum(not_ill_samples)

# ---- Mean relative abundances (percent) ----
mean_abund_ill <- if (n_ill > 0) colMeans(rel_abund[ill_samples, top_taxa, drop = FALSE]) * 100 else rep(NA_real_, length(top_taxa))
mean_abund_not_ill <- if (n_not_ill > 0) colMeans(rel_abund[not_ill_samples, top_taxa, drop = FALSE]) * 100 else rep(NA_real_, length(top_taxa))

# ---- Convert to proportions per taxon and cap ----
total_mean_abund <- mean_abund_ill + mean_abund_not_ill
prop_abund_ill    <- mean_abund_ill    / total_mean_abund
prop_abund_not_ill <- mean_abund_not_ill / total_mean_abund
prop_abund_ill[!is.finite(prop_abund_ill)] <- 0
prop_abund_not_ill[!is.finite(prop_abund_not_ill)] <- 0
prop_abund_ill    <- pmin(prop_abund_ill, 0.6)
prop_abund_not_ill <- pmin(prop_abund_not_ill, 0.6)

# ---- Cluster colors for each taxon ----
taxa_clusters <- cluster_assignments[match(top_taxa, names(cluster_assignments))]
taxa_cluster_colors <- cluster_colors[as.character(taxa_clusters)]

# ---- Shorten and mark taxon names ----
shorten_name <- function(x) {
  x <- gsub("_", " ", x)
  parts <- unlist(strsplit(x, " "))
  if (length(parts) >= 2) paste(parts[1], parts[2], sep = " ") else x
}
legend_taxa_short <- vapply(top_taxa, shorten_name, character(1))

# italicize taxa for annotations only
italicize <- function(x) as.expression(bquote(italic(.(x))))

taxa_with_asterisk <- c("Prevotella copri", "Faecalibacterium prausnitzii", 
                        "Prevotella stercorea", "Bacteroides", 
                        "Bifidobacterium", "Escherichia coli")

legend_taxa_with_star <- vapply(legend_taxa_short, function(taxon) {
  if (any(grepl(paste0("^", taxa_with_asterisk, collapse = "|"), taxon))) {
    paste0(taxon, " *")
  } else taxon
}, character(1))

formatted_taxa <- legend_taxa_with_star
formatted_taxa_italic <- lapply(formatted_taxa, italicize)

cluster_taxon_names <- tapply(formatted_taxa, taxa_clusters, function(x) paste(x[1], collapse = ""))

# ---- Top bar annotations ----
col_ha <- HeatmapAnnotation(
  Ill = anno_barplot(prop_abund_ill, gp = gpar(fill = taxa_cluster_colors), height = unit(1, "cm"), ylim = c(0, 0.6)),
  `Not-Ill` = anno_barplot(prop_abund_not_ill, gp = gpar(fill = taxa_cluster_colors), height = unit(1, "cm"), ylim = c(0, 0.6),
                           axis_param = list(at = c(0, 0.2, 0.4), labels = c("0", "0.2", "0.4"))),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)

# ---- Color function ----
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# ---- E. coli and P. stercorea colors ----
ecoli_index <- which(top_taxa == "Escherichia coli" | grepl("Escherichia.*coli", top_taxa))
if (length(ecoli_index) == 0) stop("E. coli taxon not found in top_taxa")
e_coli_cluster <- taxa_clusters[ecoli_index]
e_coli_color <- cluster_colors[as.character(e_coli_cluster)]
p_stercorea_color <- "#9ACD32"

# ---- Heatmap ----
ht <- Heatmap(
  cor_matrix,
  name = "Spearman Correlation",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_labels = formatted_taxa,
  column_labels = formatted_taxa,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_title = "Top 50 Taxa",
  row_title = "Top 50 Taxa",
  top_annotation = col_ha,
  right_annotation = rowAnnotation(
    Cluster = factor(taxa_clusters),
    col = list(Cluster = cluster_colors),
    annotation_name_gp = gpar(fontsize = 10),
    annotation_legend_param = list(
      title = "Cluster",
      at = names(cluster_taxon_names),
      labels = paste0(names(cluster_taxon_names), ": ", unname(cluster_taxon_names)),
      labels_gp = gpar(fontsize = 8, fontface = "italic")
    ),
    show_legend = TRUE
  ),
  width = unit(12, "cm"),
  height = unit(12, "cm")
)

# ---- Draw PDF ----
pdf("C:/Users/oofordile/Desktop/Microbial_Network_Heatmap.pdf", width = 12, height = 12)
ht_drawn <- draw(ht, merge_legend = TRUE)

# ---- Annotation boxes ----
row_order_vec <- top_taxa[row_order(ht_drawn)]
col_order_vec <- top_taxa[column_order(ht_drawn)]
total_rows <- length(row_order_vec)
total_cols <- length(col_order_vec)

# ---- Adjust P. stercorea: top-left x0,y1 but keep bottom-right the same ----
x_new0 <- 0       # move top-left x to 0
y_new0 <- 1       # move top-left y to 1
x_new1 <- 30 / total_cols  # bottom-right remains
y_new1 <- (total_rows - 31 + 1) / total_rows

x_ecoli0 <- 30 / total_cols
y_ecoli1 <- (total_rows - 31 + 1) / total_rows
x_ecoli1 <- 50 / total_cols
y_ecoli0 <- (total_rows - 51 + 1) / total_rows

decorate_heatmap_body("Spearman Correlation", {
  # P. stercorea
  grid.rect(x = unit((x_new0 + x_new1) / 2, "npc"), y = unit((y_new0 + y_new1) / 2, "npc"),
            width = unit(x_new1 - x_new0, "npc"), height = unit(y_new1 - y_new0, "npc"),
            gp = gpar(col = p_stercorea_color, lwd = 2, fill = NA))
  grid.text(label = bquote(italic("P. stercorea") ~ "Group"), 
            x = unit((x_new0 + x_new1) / 2, "npc"), 
            y = unit((y_new0 + y_new1) / 2, "npc"), 
            just = "center", gp = gpar(col = "black", fontsize = 10))
  # E. coli
  grid.rect(x = unit((x_ecoli0 + x_ecoli1) / 2, "npc"), y = unit((y_ecoli0 + y_ecoli1) / 2, "npc"),
            width = unit(x_ecoli1 - x_ecoli0, "npc"), height = unit(y_ecoli1 - y_ecoli0, "npc"),
            gp = gpar(col = e_coli_color, lwd = 2, fill = NA))
  grid.text(label = bquote(italic("E. coli") ~ "Group"), 
            x = unit((x_ecoli0 + x_ecoli1) / 2, "npc"), 
            y = unit((y_ecoli0 + y_ecoli1) / 2, "npc"),
            just = "center", gp = gpar(col = "black", fontsize = 10))
})
dev.off()

# ---- Draw EMF ----
emf("C:/Users/oofordile/Desktop/Microbial_Network_Heatmap.emf", width = 12, height = 12)
ht_drawn <- draw(ht, merge_legend = TRUE)
decorate_heatmap_body("Spearman Correlation", {
  # P. stercorea
  grid.rect(x = unit((x_new0 + x_new1) / 2, "npc"), y = unit((y_new0 + y_new1) / 2, "npc"),
            width = unit(x_new1 - x_new0, "npc"), height = unit(y_new1 - y_new0, "npc"),
            gp = gpar(col = p_stercorea_color, lwd = 2, fill = NA))
  grid.text(label = bquote(italic("P. stercorea") ~ "Group"), 
            x = unit((x_new0 + x_new1) / 2, "npc"), 
            y = unit((y_new0 + y_new1) / 2, "npc"), 
            just = "center", gp = gpar(col = "black", fontsize = 10))
  # E. coli
  grid.rect(x = unit((x_ecoli0 + x_ecoli1) / 2, "npc"), y = unit((y_ecoli0 + y_ecoli1) / 2, "npc"),
            width = unit(x_ecoli1 - x_ecoli0, "npc"), height = unit(y_ecoli1 - y_ecoli0, "npc"),
            gp = gpar(col = e_coli_color, lwd = 2, fill = NA))
  grid.text(label = bquote(italic("E. coli") ~ "Group"), 
            x = unit((x_ecoli0 + x_ecoli1) / 2, "npc"), 
            y = unit((y_ecoli0 + y_ecoli1) / 2, "npc"),
            just = "center", gp = gpar(col = "black", fontsize = 10))
})
dev.off()

# ---- Console diagnostics ----
cat("Top taxa selected:", length(top_taxa), "\n")
cat("Ill samples:", n_ill, " | Not-Ill samples:", n_not_ill, "\n")
cat("Any bars non-zero? Ill:", any(prop_abund_ill > 0), " Not-Ill:", any(prop_abund_not_ill > 0), "\n")
