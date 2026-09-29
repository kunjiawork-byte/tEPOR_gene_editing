# tEPOR gene editing study: donor random-intercept analysis
# Model, response transformation, comparisons, and adjustment settings retained.
# Set MODEL_INPUT and MODEL_OUTPUT, or edit the settings below.
output_dir <- Sys.getenv("MODEL_OUTPUT", ".")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Mixed model for unbalanced multi-donor data (lme4 + emmeans)
# - Model: log2(value) ~ condition + (1 | donor)
# - Outputs:
#   1) Dunnett-adjusted comparisons vs a control condition
#   2) Custom "+tEPOR vs -tEPOR" contrasts (Holm-adjusted)
# ============================================================

# --------- USER SETTINGS (edit these) -----------------------
input_csv <- Sys.getenv("MODEL_INPUT", "day14_long_format.csv")   # your long-format file
value_col <- "value_day14"             # column holding numeric response
donor_col <- "donor"
condition_col <- "condition"

control_condition <- "CCR5_CBE"        # reference level for control-vs-all

# Define paired contrasts (pos = +tEPOR, neg = -tEPOR)
# IMPORTANT: names must match exactly levels(df$condition)
paired_contrasts <- list(
  "B_E vs B"                   = c(pos = "B_E",           neg = "B"),
  "H_E vs HBG_CBE"             = c(pos = "H_E",           neg = "HBG_CBE"),
  "H_B_E vs H_B"               = c(pos = "H_B_E",         neg = "H_B"),
  "Casgevy_tEPOR vs Casgevy"   = c(pos = "Casgevy_tEPOR", neg = "Casgevy")
)

# Output file names
out_vs_control_csv <- file.path(output_dir, "mixedmodel_vs_control.csv")
out_paired_csv     <- file.path(output_dir, "mixedmodel_paired_contrasts.csv")
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

# --------- Load & basic checks ------------------------------
df <- read.csv(input_csv, stringsAsFactors = FALSE)

required <- c(donor_col, condition_col, value_col)
missing_cols <- setdiff(required, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Coerce types
df[[donor_col]]     <- factor(df[[donor_col]])
df[[condition_col]] <- factor(df[[condition_col]])
df[[value_col]]     <- as.numeric(df[[value_col]])

# Remove non-positive values (log2 needs >0)
bad <- which(!is.finite(df[[value_col]]) | df[[value_col]] <= 0)
if (length(bad) > 0) {
  warning("Removing ", length(bad), " rows with non-finite or <=0 values in ", value_col)
  df <- df[-bad, , drop = FALSE]
}

# Confirm control exists
levs_raw <- levels(df[[condition_col]])
if (!(control_condition %in% levs_raw)) {
  stop("Control condition '", control_condition, "' not found. Available levels:\n  ",
       paste(levs_raw, collapse = ", "))
}

# --------- Transform & fit model -----------------------------
df$log2_value <- log2(df[[value_col]])
df[[condition_col]] <- relevel(df[[condition_col]], ref = control_condition)

# IMPORTANT: build formula as text (NO get() inside formula)
form_fit <- as.formula(
  paste0("log2_value ~ ", condition_col, " + (1|", donor_col, ")")
)

fit <- lmer(form_fit, data = df)

cat("\n=== Model summary ===\n")
print(summary(fit))

# --------- (1) Comparisons vs control (Dunnett) --------------
# IMPORTANT: emmeans spec must refer to factor name, not get()
emm_spec <- as.formula(paste0("~ ", condition_col))
emm <- emmeans(fit, emm_spec)

# For trt.vs.ctrl, ref should be index or level name that exists in emm
res_vs_ctrl <- contrast(
  emm,
  method = "trt.vs.ctrl",
  ref = which(levels(df[[condition_col]]) == control_condition)
)

res_vs_ctrl_df <- as.data.frame(res_vs_ctrl)
res_vs_ctrl_df$fold_change_vs_control <- 2^(res_vs_ctrl_df$estimate)

write.csv(res_vs_ctrl_df, out_vs_control_csv, row.names = FALSE)
cat("\nSaved vs-control results to: ", out_vs_control_csv, "\n", sep = "")

# --------- (2) Custom paired contrasts (+tEPOR vs -tEPOR) ----
mk_contrast <- function(pos, neg, levs) {
  v <- setNames(rep(0, length(levs)), levs)
  if (!(pos %in% levs)) stop("Contrast pos level not found: ", pos)
  if (!(neg %in% levs)) stop("Contrast neg level not found: ", neg)
  v[pos] <-  1
  v[neg] <- -1
  v
}

# IMPORTANT: use the levels as emmeans sees them
levs <- levels(df[[condition_col]])

ct_list <- list()
for (nm in names(paired_contrasts)) {
  pos <- paired_contrasts[[nm]][["pos"]]
  neg <- paired_contrasts[[nm]][["neg"]]
  ct_list[[nm]] <- mk_contrast(pos, neg, levs)
}

ct <- contrast(emm, method = ct_list, adjust = "holm")
ct_df <- as.data.frame(ct)
ct_df$fold_change_plus_vs_minus <- 2^(ct_df$estimate)

write.csv(ct_df, out_paired_csv, row.names = FALSE)
cat("Saved paired-contrast results to: ", out_paired_csv, "\n", sep = "")

# --------- Print key tables to console -----------------------
cat("\n=== Vs-control (Dunnett-adjusted) ===\n")
print(res_vs_ctrl_df)

cat("\n=== Paired contrasts (Holm-adjusted) ===\n")
print(ct_df)

cat("\nDone.\n")
