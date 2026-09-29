# tEPOR gene editing study: donor random-intercept analysis
# Model, response transformation, comparisons, and adjustment settings retained.
# Set MODEL_INPUT and MODEL_OUTPUT, or edit the settings below.
output_dir <- Sys.getenv("MODEL_OUTPUT", ".")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Mixed model for unbalanced multi-donor data (lme4 + emmeans)
# Compare selected conditions vs unedited mock
#
# Model:
#   log2(value) ~ condition + (1 | donor)
#
# Outputs:
#   1) Selected comparisons vs unedited mock
#   2) Multiplicity-adjusted p values (Dunnett)
# ============================================================

# ---------------- USER SETTINGS ------------------------------
input_csv <- Sys.getenv("MODEL_INPUT", "day3_postediting.csv")   # long-format file
value_col <- "value_day3"             # numeric response column
donor_col <- "donor"
condition_col <- "condition"

control_condition <- "unedited mock"

# Conditions to compare against unedited mock
target_conditions <- c(
  "1 x BE",
  "2 x BE",
  "3 x BE",
  "Casgevy",
  "CE",
  "H",
  "L",
  "BE Control",
  "Cas9 Control"
)

# Output file
out_vs_mock_csv <- file.path(output_dir, "mixedmodel_vs_unedited_mock.csv")
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

# ---------------- Load data ---------------------------------
df <- read.csv(input_csv, stringsAsFactors = FALSE)

required <- c(donor_col, condition_col, value_col)
missing_cols <- setdiff(required, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Coerce types
df[[donor_col]]     <- factor(df[[donor_col]])
df[[condition_col]] <- as.character(df[[condition_col]])
df[[value_col]]     <- as.numeric(df[[value_col]])

# Remove rows with missing condition/donor/value
df <- df[!is.na(df[[donor_col]]) & !is.na(df[[condition_col]]) & !is.na(df[[value_col]]), , drop = FALSE]

# Remove non-positive values (log2 requires > 0)
bad <- which(!is.finite(df[[value_col]]) | df[[value_col]] <= 0)
if (length(bad) > 0) {
  warning("Removing ", length(bad), " rows with non-finite or <=0 values in ", value_col)
  df <- df[-bad, , drop = FALSE]
}

# Check control condition
all_conditions <- unique(df[[condition_col]])
if (!(control_condition %in% all_conditions)) {
  stop("Control condition '", control_condition, "' not found.\nAvailable conditions:\n  ",
       paste(sort(all_conditions), collapse = ", "))
}

# Check target conditions
missing_targets <- setdiff(target_conditions, all_conditions)
if (length(missing_targets) > 0) {
  warning("These target conditions were not found and will be skipped:\n  ",
          paste(missing_targets, collapse = ", "))
}

target_conditions_found <- intersect(target_conditions, all_conditions)
if (length(target_conditions_found) == 0) {
  stop("None of the target conditions were found in the dataset.")
}

# Keep only control + requested targets
keep_conditions <- c(control_condition, target_conditions_found)
df_sub <- df[df[[condition_col]] %in% keep_conditions, , drop = FALSE]

# Re-factor condition and set reference
df_sub[[condition_col]] <- factor(df_sub[[condition_col]], levels = keep_conditions)
df_sub[[condition_col]] <- relevel(df_sub[[condition_col]], ref = control_condition)

# ---------------- Transform & fit model ----------------------
df_sub$log2_value <- log2(df_sub[[value_col]])

form_fit <- as.formula(
  paste0("log2_value ~ ", condition_col, " + (1|", donor_col, ")")
)

fit <- lmer(form_fit, data = df_sub)

cat("\n=== Model summary ===\n")
print(summary(fit))

# ---------------- Estimated marginal means -------------------
emm_spec <- as.formula(paste0("~ ", condition_col))
emm <- emmeans(fit, emm_spec)

# Dunnett-adjusted treatment-vs-control comparisons
res_vs_mock <- contrast(
  emm,
  method = "trt.vs.ctrl",
  ref = which(levels(df_sub[[condition_col]]) == control_condition),
  adjust = "dunnett"
)

res_vs_mock_df <- as.data.frame(res_vs_mock)

# Keep only requested targets (should already be true, but extra safe)
wanted_labels <- paste(target_conditions_found, "-", control_condition)
res_vs_mock_df <- res_vs_mock_df[res_vs_mock_df$contrast %in% wanted_labels, , drop = FALSE]

# Add fold change on original scale
res_vs_mock_df$fold_change_vs_mock <- 2^(res_vs_mock_df$estimate)

# Optional: nicer column naming/order
res_vs_mock_df <- res_vs_mock_df[, c(
  "contrast", "estimate", "SE", "df", "t.ratio", "p.value", "fold_change_vs_mock"
)]

write.csv(res_vs_mock_df, out_vs_mock_csv, row.names = FALSE)
cat("\nSaved vs-mock results to: ", out_vs_mock_csv, "\n", sep = "")

# ---------------- Print results ------------------------------
cat("\n=== Selected conditions vs unedited mock (Dunnett-adjusted) ===\n")
print(res_vs_mock_df)

cat("\nDone.\n")
