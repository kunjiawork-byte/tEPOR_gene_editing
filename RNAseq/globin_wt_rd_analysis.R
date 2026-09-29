# tEPOR gene editing study: D14 globin donor adjustment
# Input: an existing d14 data frame, or a CSV specified by GLOBIN_INPUT.
# Required columns: Donor, Condition, Gamma_to_Alpha; input must already be D14.
# Original averaging, log2 transformation, fixed-effects model, and plots are retained.

input_file <- Sys.getenv("GLOBIN_INPUT", "")
output_dir <- Sys.getenv("GLOBIN_OUTPUT", ".")
if (nzchar(input_file)) {
  d14 <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
}
if (!exists("d14")) {
  stop("Provide the D14 data frame d14 or set GLOBIN_INPUT to a D14 CSV file.")
}
if (!all(c("Donor", "Condition", "Gamma_to_Alpha") %in% names(d14))) {
  stop("Input requires Donor, Condition, and Gamma_to_Alpha columns.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# D14 donor-effect correction
# Use ONLY donors 22, 24, 25
# and conditions shared by all three donors
# ============================================================

library(tidyverse)


# ------------------------------------------------------------
# 1. Define donors
# ------------------------------------------------------------

keep_donors <- c(
  "22",
  "24",
  "25"
)


# ------------------------------------------------------------
# 2. Define common conditions
# ------------------------------------------------------------

keep_conditions <- c(
  "CCR5 CBE",
  "HBG CBE",
  "tEPOR CBE",
  "HBG+tEPOR CBE",
  "BCL11A CBE",
  "BCL11A+tEPOR CBE"
)


# ------------------------------------------------------------
# 3. Select balanced D14 dataset
# ------------------------------------------------------------

d14_balanced <- d14 %>%

  filter(
    Donor %in% keep_donors,
    Condition %in% keep_conditions
  ) %>%

  group_by(
    Condition,
    Donor
  ) %>%

  summarise(

    Gamma_to_Alpha =
      mean(
        Gamma_to_Alpha,
        na.rm = TRUE
      ),

    n_rep = n(),

    .groups = "drop"
  )


# Set order
d14_balanced$Condition <- factor(
  d14_balanced$Condition,
  levels = keep_conditions
)

d14_balanced$Donor <- factor(
  d14_balanced$Donor,
  levels = keep_donors
)


# Check design
print(
  table(
    d14_balanced$Condition,
    d14_balanced$Donor
  )
)



# ============================================================
# 4. log2 transform
# ============================================================

d14_balanced <- d14_balanced %>%

  mutate(

    log2_Gamma_to_Alpha =
      log2(Gamma_to_Alpha)

  )


# ============================================================
# 5. Fit model
#
# Condition + Donor
# ============================================================

fit <- lm(
  log2_Gamma_to_Alpha ~ Condition + Donor,
  data = d14_balanced
)

summary(fit)


# ============================================================
# 6. Model matrix
# ============================================================

X <- model.matrix(fit)

beta <- coef(fit)



# ============================================================
# 7. Calculate donor contribution
# ============================================================

donor_columns <- grep(
  "^Donor",
  colnames(X)
)


donor_effect <- X[
  ,
  donor_columns,
  drop = FALSE
] %*%
  beta[donor_columns]


# ============================================================
# 8. Center donor effects
#
# IMPORTANT:
# We don't want donor 22 arbitrarily to define the absolute
# level of the adjusted data.
#
# Center donor contribution around its overall mean.
# ============================================================

donor_effect_centered <-
  as.numeric(donor_effect) -
  mean(as.numeric(donor_effect))


# ============================================================
# 9. Remove donor effect
# ============================================================

d14_balanced <- d14_balanced %>%

  mutate(

    log2_Gamma_to_Alpha_adjusted =
      log2_Gamma_to_Alpha -
      donor_effect_centered,

    Gamma_to_Alpha_adjusted =
      2^log2_Gamma_to_Alpha_adjusted

  )


# ============================================================
# 10. Check result
# ============================================================

comparison <- d14_balanced %>%

  select(

    Condition,
    Donor,

    Gamma_to_Alpha,

    Gamma_to_Alpha_adjusted

  )


print(
  comparison,
  n = Inf
)


# ============================================================
# 11. Save
# ============================================================

write.csv(
  comparison,
  file.path(output_dir, "D14_gamma_alpha_donor22_24_25_adjusted.csv"),
  row.names = FALSE
)


# ============================================================
# 12. RAW paired plot
# ============================================================

p_raw <- ggplot(
  d14_balanced,
  aes(
    x = Condition,
    y = Gamma_to_Alpha,
    group = Donor
  )
) +

  geom_line(
    linewidth = 0.5,
    alpha = 0.5
  ) +

  geom_point(
    aes(shape = Donor),
    size = 3
  ) +

  labs(

    x = NULL,

    y = expression(
      gamma / alpha ~
        "globin transcript ratio"
    ),

    shape = "Donor"

  ) +

  theme_classic(
    base_size = 13
  ) +

  theme(

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),

    axis.text.y = element_text(
      color = "black"
    ),

    legend.position = "right"

  )


print(p_raw)


# ============================================================
# 13. DONOR-ADJUSTED plot
# ============================================================

p_adjusted <- ggplot(
  d14_balanced,
  aes(
    x = Condition,
    y = Gamma_to_Alpha_adjusted,
    group = Donor
  )
) +

  geom_line(
    linewidth = 0.5,
    alpha = 0.5
  ) +

  geom_point(
    aes(shape = Donor),
    size = 3
  ) +

  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "crossbar",
    width = 0.45,
    linewidth = 0.5
  ) +

  labs(

    x = NULL,

    y = expression(
      "Donor-adjusted " *
        gamma / alpha ~
        "globin transcript ratio"
    ),

    shape = "Donor"

  ) +

  theme_classic(
    base_size = 13
  ) +

  theme(

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),

    axis.text.y = element_text(
      color = "black"
    ),

    legend.position = "right"

  )


print(p_adjusted)


# ============================================================
# 14. Save figures
# ============================================================

ggsave(
  file.path(output_dir, "D14_gamma_alpha_donor22_24_25_RAW.pdf"),
  p_raw,
  width = 9,
  height = 5
)

ggsave(
  file.path(output_dir, "D14_gamma_alpha_donor22_24_25_adjusted.pdf"),
  p_adjusted,
  width = 9,
  height = 5
)

ggsave(
  file.path(output_dir, "D14_gamma_alpha_donor22_24_25_adjusted.png"),
  p_adjusted,
  width = 9,
  height = 5,
  dpi = 600
)
