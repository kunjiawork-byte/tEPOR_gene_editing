# tEPOR gene editing study: donor random-intercept analysis
# Model, response transformation, comparisons, and adjustment settings retained.
# Set MODEL_INPUT and MODEL_OUTPUT, or edit the settings below.
output_dir <- Sys.getenv("MODEL_OUTPUT", ".")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# run_mixed_model_spike.R
# HbF mixed model (donor random intercept) + emmeans + contrasts + fold-change
# - vs Ctrl: explicit dunnettx adjustment (retained from the original)
# - within-group (+tEPOR vs -tEPOR): Holm adjustment across 4 tests

suppressPackageStartupMessages({
  library(lme4)
  library(emmeans)
})

# ---------------------------
# USER SETTINGS (edit if needed)
# ---------------------------
input_csv <- Sys.getenv("MODEL_INPUT", "HbF_spike.csv")

value_col     <- "HbF"
condition_col <- "condition"
donor_col     <- "donor"

control_condition <- "Ctrl"

# Define within-group paired contrasts (A = +tEPOR, B = -tEPOR)
# IMPORTANT: names must match EXACTLY levels(df$condition)
paired_contrasts <- list(
  "H_E vs H"                 = c(pos = "H_E",       neg = "H"),
  "B_E vs B"                 = c(pos = "B_E",       neg = "B"),
  "H_B_E vs H_B"             = c(pos = "H_B_E",     neg = "H_B"),
  "Casgevy_E vs Casgevy"     = c(pos = "Casgevy_E", neg = "Casgevy")
)

# ---------------------------
# LOAD DATA
# ---------------------------
df <- read.csv(input_csv, check.names = FALSE)

needed_cols <- c(value_col, condition_col, donor_col)
missing_cols <- setdiff(needed_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing columns in CSV: ", paste(missing_cols, collapse = ", "))
}

df[[condition_col]] <- factor(df[[condition_col]])
if (!(control_condition %in% levels(df[[condition_col]]))) {
  stop(
    "Control condition '", control_condition, "' not found in condition levels: ",
    paste(levels(df[[condition_col]]), collapse = ", ")
  )
}
df[[condition_col]] <- relevel(df[[condition_col]], ref = control_condition)
df[[donor_col]] <- factor(df[[donor_col]])

cat("\n=== Data check ===\n")
cat("Rows:", nrow(df), "\n")
cat("Condition levels:\n  ", paste(levels(df[[condition_col]]), collapse = " | "), "\n")
cat("Donor n:", length(levels(df[[donor_col]])), "\n")

# ---------------------------
# TRANSFORM
# ---------------------------
if (any(is.na(df[[value_col]]))) stop(value_col, " column contains NA values.")
if (any(df[[value_col]] <= 0)) stop(value_col, " contains non-positive values; cannot take log2.")
df$log2_value <- log2(df[[value_col]])

# ---------------------------
# FIT MODEL
# log2(HbF) ~ condition + (1|donor)
# ---------------------------
form <- as.formula(paste0("log2_value ~ ", condition_col, " + (1|", donor_col, ")"))
fit <- lmer(form, data = df)

cat("\n=== Model summary ===\n")
print(summary(fit))

# ---------------------------
# EMMEANS
# ---------------------------
emm <- emmeans(fit, as.formula(paste("~", condition_col)))

cat("\n=== EMMEANS (log2 scale) ===\n")
print(emm)

# ---------------------------
# CONTRASTS vs Ctrl (Dunnett) + fold-change
# ---------------------------
cat("\n=== Contrasts: each level vs Ctrl (Dunnett-adjusted) ===\n")

# Use index for ref to avoid any ambiguity
ref_idx <- which(levels(df[[condition_col]]) == control_condition)

ctr_vs_ctrl <- contrast(
  emm,
  method = "trt.vs.ctrl",
  ref = ref_idx,
  adjust = "dunnettx"   # explicit
)

print(ctr_vs_ctrl)

ctr_vs_ctrl_df <- as.data.frame(ctr_vs_ctrl)
ctr_vs_ctrl_df$fold_change_vs_ctrl <- 2^(ctr_vs_ctrl_df$estimate)

cat("\n=== Contrasts vs Ctrl with fold-change (2^estimate) ===\n")
print(ctr_vs_ctrl_df)

# ---------------------------
# WITHIN-GROUP CONTRASTS (Holm-adjusted across 4 tests)
# ---------------------------
cat("\n=== Within-group contrasts (Holm-adjusted) ===\n")

# Levels as emmeans sees them
levs <- levels(df[[condition_col]])

mk_contrast <- function(pos, neg, levs) {
  v <- setNames(rep(0, length(levs)), levs)
  if (!(pos %in% levs)) stop("Contrast pos level not found: ", pos)
  if (!(neg %in% levs)) stop("Contrast neg level not found: ", neg)
  v[pos] <-  1
  v[neg] <- -1
  v
}

# Build list of contrasts; stop early with a clear message if any level missing
ct_list <- list()
for (nm in names(paired_contrasts)) {
  pos <- paired_contrasts[[nm]][["pos"]]
  neg <- paired_contrasts[[nm]][["neg"]]
  ct_list[[nm]] <- mk_contrast(pos, neg, levs)
}

within_ct <- contrast(emm, method = ct_list, adjust = "holm")
within_df <- as.data.frame(within_ct)
within_df$fold_change_plus_vs_minus <- 2^(within_df$estimate)

cat("\n--- Within-group contrasts (log2 scale) + fold-change (2^estimate) ---\n")
print(within_df[, c("contrast","estimate","SE","df","t.ratio","p.value","fold_change_plus_vs_minus")])

cat("\nDONE.\n")

# Save the same result tables printed above; no additional statistical tests.
write.csv(ctr_vs_ctrl_df, file.path(output_dir, "HbF_spike_vs_Ctrl.csv"), row.names = FALSE)
write.csv(within_df, file.path(output_dir, "HbF_spike_within_group_contrasts.csv"), row.names = FALSE)
