# --- Load Required Libraries ---
library(dplyr)
library(ggplot2)
library(devEMF)   # for EMF export
library(grid)     # for unit() used in tick length

# --- Load Dataset ---
data_path <- "C:/Users/oofordile/Desktop/Merged_Illness_Cohorts.csv"
data <- read.csv(data_path)

# --- Remove rows with blank or NA Ill_Status ---
data <- data %>% filter(!is.na(Ill_Status) & Ill_Status != "")

# --- Recode Variables ---
data <- data %>%
  mutate(
    Ill_Status = factor(Ill_Status, levels = c("Not-Ill", "Ill")),
    Iron = factor(Iron, levels = c("placebo", "treatment")),
    Age_Sampling_Group = factor(X3agegroups.based.on.age.at.sampling)
  )

# --- Fit GLM Model Adjusted for Age at Sampling ---
model <- glm(Ill_Status ~ Iron + Age_Sampling_Group,
             data = data, family = binomial(link = "logit"))

# --- Extract ORs, CIs, and P-value for Iron ---
OR <- exp(coef(model))
CI <- exp(confint(model))
p_val <- summary(model)$coefficients["Irontreatment", "Pr(>|z|)"]

# --- Format Output Table (all terms) ---
results <- data.frame(
  Term = names(OR),
  Odds_Ratio = round(OR, 3),
  CI_Lower = round(CI[, 1], 3),
  CI_Upper = round(CI[, 2], 3)
)

# --- Keep Only Iron Treatment Row for Reporting/Export ---
iron_results <- results %>% filter(Term == "Irontreatment")

# --- Print Results ---
print(iron_results)
cat("\nP-value for Iron treatment effect:", p_val, "\n")

# --- Save Iron Results to CSV ---
csv_path <- "C:/Users/oofordile/Desktop/Iron_Treatment_Illness_OR_Adjusted_SamplingAge.csv"
write.csv(iron_results, csv_path, row.names = FALSE)

# --- Format P-value for Subtitle ---
p_text <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", signif(p_val, 2)))
subtitle_text <- paste("No significant association observed (", p_text, ")", sep = "")

# --- Define Output Paths ---
pdf_path <- "C:/Users/oofordile/Desktop/Iron_Treatment_Illness_OR_Adjusted_SamplingAge.pdf"
emf_path <- "C:/Users/oofordile/Desktop/Iron_Treatment_Illness_OR_Adjusted_SamplingAge.emf"

# --- Create Plot with thin axis lines and inward ticks ---
or_plot <- ggplot(iron_results, aes(x = Odds_Ratio, y = Term)) +
  geom_point(size = 4, color = "#1f77b4") +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  xlab("Odds Ratio (95% CI)") +
  ylab("Iron Treatment vs Placebo") +
  ggtitle("Effect of Iron Supplementation on Illness Risk",
          subtitle = subtitle_text) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(-0.2, "cm"),   # inward ticks
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

# --- Export PDF ---
pdf(pdf_path, width = 7, height = 5)
print(or_plot)
dev.off()

# --- Export EMF ---
emf(file = emf_path, width = 7, height = 5, family = "Arial")
print(or_plot)
dev.off()

# --- Confirmation ---
cat("Adjusted outputs saved to:\nCSV:", csv_path,
    "\nPDF:", pdf_path,
    "\nEMF:", emf_path, "\n")
