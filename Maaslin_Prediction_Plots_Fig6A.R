# --- Load Required Libraries ---
library(grid)
library(ggplot2) 
library(gridExtra)
library(devEMF)

# --- Set base path ---
base_path <- "C:/Users/oofordile/Desktop/Maaslin2_Predictive"

# --- Folder paths ---
folders <- list(
  illness_D1 = file.path(base_path, "maaslin2_outputD1Predictive"),
  illness_D15 = file.path(base_path, "maaslin2_outputD15Predictive"),
  total_D1 = file.path(base_path, "D1_Ill_vs_D85_NotIll"),
  total_D15 = file.path(base_path, "D15_Ill_vs_D85_NotIll")
)

# --- Group colors ---
group_colors <- c(
  "D85 Not-Ill"    = "#4575B4",  # blue
  "Soon-Ill"       = "#FFA500",  # orange
  "Later-Ill"      = "#F0E442",  # yellow
  "Much-Later-Ill" = "#3CA63C"   # green
)

# --- Load data ---
illness_d1 <- read.delim(file.path(folders$illness_D1, "significant_results.tsv"), stringsAsFactors = FALSE)
illness_d15 <- read.delim(file.path(folders$illness_D15, "significant_results.tsv"), stringsAsFactors = FALSE)
total_d1 <- read.delim(file.path(folders$total_D1, "significant_results.tsv"), stringsAsFactors = FALSE)
total_d15 <- read.delim(file.path(folders$total_D15, "significant_results.tsv"), stringsAsFactors = FALSE)

# --- Groups of interest ---
groups <- c("Not-Ill", "Soon-Ill", "Later-Ill", "Much-Later-Ill")

# --- Filter function ---
filter_results <- function(df, group, timepoint) {
  df_sub <- df[df$value == group & df$metadata == "SampleGroup", ]
  if (timepoint == "D15" && group %in% c("Later-Ill", "Much-Later-Ill")) {
    df_sub <- df_sub[df_sub$p < 0.1, ]
  } else {
    df_sub <- df_sub[df_sub$q < 0.1, ]
  }
  return(df_sub)
}

# --- Apply filters ---
filtered_d1 <- lapply(groups, function(g) filter_results(illness_d1, g, "D1"))
names(filtered_d1) <- groups
filtered_d15 <- lapply(groups, function(g) filter_results(illness_d15, g, "D15"))
names(filtered_d15) <- groups

# --- Filter total ---
filtered_total_d1 <- total_d1[total_d1$q < 0.1 & total_d1$metadata == "SampleGroup", ]
filtered_total_d15 <- total_d15[total_d15$q < 0.1 & total_d15$metadata == "SampleGroup", ]
filtered_d1[["Total"]] <- filtered_total_d1
filtered_d15[["Total"]] <- filtered_total_d15

# --- Function to clean taxa names ---
clean_taxa_names <- function(taxa_name) {
  cleaned <- gsub("_\\d+\\.\\d+$", "", taxa_name)
  cleaned <- gsub("_", " ", cleaned)
  return(cleaned)
}

# --- Get group color ---
get_group_color <- function(value) {
  if (value %in% names(group_colors)) {
    return(group_colors[[value]])
  }
  return("#000000")
}

# --- Bar plot function with italic taxa labels, borders, and inward ticks ---
plot_bar <- function(df, group_label, timepoint) {
  if (nrow(df) == 0) return(NULL)
  
  df$feature_clean <- clean_taxa_names(df$feature)
  df <- df[order(df$coef), ]
  df$feature_clean <- factor(df$feature_clean, levels = df$feature_clean)
  
  # Italic labels for y-axis
  df$feature_label <- paste0("italic('", df$feature_clean, "')")
  
  df$bar_color <- sapply(seq_len(nrow(df)), function(i) {
    coef_i <- df$coef[i]
    if (coef_i < 0) {
      return("#4575B4")  # blue for negative
    } else {
      if (group_label == "Total") {
        return("#D73027")  # red for positive in Total
      } else {
        return(get_group_color(df$value[i]))
      }
    }
  })
  
  df$label_val <- ifelse(
    timepoint == "D15" & df$value %in% c("Later-Ill", "Much-Later-Ill"),
    paste0("p=", signif(df$p, 3)),
    paste0("q=", signif(df$q, 3))
  )
  
  bar_range <- max(abs(df$coef)) * 1.2
  df$x_text <- ifelse(df$coef > 0,
                      df$coef - 0.02 * bar_range,
                      df$coef + 0.02 * bar_range)
  df$hjust_val <- ifelse(df$coef > 0, 1, 0)
  
  ggplot(df, aes(x = coef, y = feature_clean, fill = bar_color)) +
    geom_col(width = 0.7) +
    geom_text(aes(x = x_text, label = label_val, hjust = hjust_val), size = 6) +
    scale_y_discrete(labels = parse(text = df$feature_label)) +   # italic taxa
    scale_fill_identity() +
    theme_bw() +
    labs(title = paste0(timepoint, ": ", group_label),
         x = "Coefficient (Effect Size)", y = NULL) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 13, face = "bold"),
      axis.text.y = element_text(size = 11),   # italics now handled by parse()
      axis.text.x = element_text(size = 11),
      axis.title.x = element_text(size = 14, face = "bold"),
      # Add thin border around plot
      panel.border = element_rect(color = "black", fill = NA, size = 0.3),
      # Add inward-pointing ticks
      axis.ticks = element_line(size = 0.3, color = "black"),
      axis.ticks.length = unit(-2, "mm"),  # Negative value makes ticks point inward
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_cartesian(clip = "off")
}

# --- Plot order ---
plot_order <- c("Soon-Ill", "Later-Ill", "Much-Later-Ill", "Total")

# --- Generate plots ---
plots_d1 <- Filter(Negate(is.null), lapply(plot_order, function(g) plot_bar(filtered_d1[[g]], g, "D1")))
plots_d15 <- Filter(Negate(is.null), lapply(plot_order, function(g) plot_bar(filtered_d15[[g]], g, "D15")))

# --- Combined layout ---
combined_plot <- arrangeGrob(
  arrangeGrob(grobs = plots_d1, nrow = 2, ncol = 2),
  arrangeGrob(grobs = plots_d15, nrow = 2, ncol = 2),
  ncol = 2
)

# --- Export PDF ---
pdf(file.path(base_path, "Maaslin2_Illness_BarPlots_FinalColored.pdf"), width = 20, height = 8)
grid.draw(combined_plot)
dev.off()

# --- Export EMF ---
emf(file.path(base_path, "Maaslin2_Illness_BarPlots_FinalColored.emf"), width = 20, height = 8)
grid.draw(combined_plot)
dev.off()

cat("Plots exported successfully!\n- PDF and EMF saved in", base_path, "\n")
