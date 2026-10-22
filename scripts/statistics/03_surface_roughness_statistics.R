# ======================================================================== #
# Surface roughness (Sq, 3D profilometry) -- diagnostics, two-way ANOVA
# with emmeans contrasts (35-day) and ART / ART-C contrasts (60-day)
#
# Implements the manuscript methods sentence:
#   "Corrosion rate and surface roughness were analyzed according to the
#    distributional characteristics of each dataset. Corrosion-rate data and
#    surface-roughness data from the 60-day incubation were analyzed using the
#    non-parametric ART procedure, with pairwise treatment contrasts obtained
#    using ART-C contrasts. Surface-roughness data from the 35-day incubation
#    were analyzed using a parametric two-way fixed-effects ANOVA. Pairwise
#    contrasts for the ANOVA were estimated using estimated marginal means, and
#    contrast P values were adjusted for multiple comparisons using the
#    Benjamini-Hochberg procedure (alpha = 0.05)."
#
# Design: 2 factors x 2 levels, in two incubation durations.
#   Sterilization  C (control) / DW (drinking water)
#   Nutrient       -nutri / +nutri
#   35 days = experiments E1 and E3, 60 days = experiment E2.
#
# Three datasets, three shapes (all in one input CSV, keyed by `Dataset`):
#   35d_polished    polished coupons, measured Before and After incubation
#                   (third factor TimeToIncubation)
#   35d_unpolished  unpolished coupons and the unpolished side of the polished
#                   coupons (third factor Polishing)
#   60d             no third factor
# Each 35-day dataset is split on its third factor and analysed as a genuine
# TWO-way ANOVA within each level -- see section 0 for why pooling is not
# defensible. The 60-day dataset fails normality on both scales and goes to
# ART instead.
#
# INPUT
#   data/physicochemical/surface_roughness.csv
#       156 measurements (77 + 41 + 38), one row per measurement.
#
# OUTPUT (into results/statistics/)
#   ANOVA_tables_surface_roughness.xlsx
#   Contrast_results_surface_roughness.xlsx
#   model_diagnostics/*_diagnostics.pdf   (4-panel base diagnostic plots)
#
# Package versions used: car 3.1.5, emmeans 1.11.1, ARTool 0.11.2.
# writexl is called namespace-qualified.
#
# Run from the repository root.
# ======================================================================== #

library(dplyr)
library(car)      # leveneTest(), Anova()
library(emmeans)  # estimated marginal means & custom contrasts
library(ARTool)   # non-parametric fallback (Aligned Rank Transform)

results_dir <- "results/statistics"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------ #
# 0a. Data
# ------------------------------------------------------------------------ #
# 156 profilometry measurements with their design factors, read from a tidy
# CSV and split into the three data frames below.
#
# NOTE on factor level order: as.factor() sorts with the locale's collation.
# In the locale this analysis was run in, sort(c("-nutri", "+nutri")) puts
# "-nutri" first. surf_rough_60_df sets its levels explicitly below; the other
# two rely on as.factor(). The contrast machinery further down is level-name
# based and therefore robust to this, EXCEPT inside ARTool -- see section 3.
surf_rough <- read.csv("data/physicochemical/surface_roughness.csv",
                       stringsAsFactors = FALSE)

sr_pol <- subset(surf_rough, Dataset == "35d_polished")
surf_rough_35_pol_df <- data.frame(
  Surf_roughness   = sr_pol$Surf_roughness,
  TimeToIncubation = as.factor(sr_pol$TimeToIncubation),
  Sterilization    = as.factor(sr_pol$Sterilization),
  Nutrient         = as.factor(sr_pol$Nutrient)
)

surf_rough_35_pol_df |>
  dplyr::group_by(TimeToIncubation, Sterilization, Nutrient) |>
  dplyr::summarise(
    n = dplyr::n(),
    median_SR = median(Surf_roughness, na.rm = TRUE),
    mean_SR = mean(Surf_roughness, na.rm = TRUE),
    range_SR = paste(sprintf("%.3f", range(Surf_roughness, na.rm = TRUE)), collapse = "-"),
    .groups = "drop"
  )

sr_unpol <- subset(surf_rough, Dataset == "35d_unpolished")
surf_rough_35_unpol_df <- data.frame(
  Surf_roughness = sr_unpol$Surf_roughness,
  Polishing      = as.factor(sr_unpol$Polishing),
  Sterilization  = as.factor(sr_unpol$Sterilization),
  Nutrient       = as.factor(sr_unpol$Nutrient)
)

sr_60 <- subset(surf_rough, Dataset == "60d")
surf_rough_60_df <- data.frame(
  Surf_roughness = sr_60$Surf_roughness,
  Sterilization  = sr_60$Sterilization,
  Nutrient       = sr_60$Nutrient
)

surf_rough_60_df |>
  dplyr::group_by(Sterilization, Nutrient) |>
  dplyr::summarise(
    n = dplyr::n(),
    median_SR = median(Surf_roughness, na.rm = TRUE),
    mean_SR = mean(Surf_roughness, na.rm = TRUE),
    range_SR = paste(sprintf("%.1f", range(Surf_roughness, na.rm = TRUE)), collapse = "-"),
    .groups = "drop"
  )

# surf_rough_60_df's Sterilization/Nutrient come in as character - make factors,
# with the same level order as in the other two data frames (C/DW, -nutri/+nutri)
surf_rough_60_df <- surf_rough_60_df |>
  dplyr::mutate(
    Sterilization = factor(Sterilization, levels = c("C", "DW")),
    Nutrient      = factor(Nutrient, levels = c("-nutri", "+nutri"))
  )

# Use sum-to-zero contrasts so Type III SS from car::Anova are interpretable.
options(contrasts = c("contr.sum", "contr.poly"))

# ------------------------------------------------------------------------ #
# 0b. Two-way or three-way ANOVA, when one factor is fixed at a single level?
# ------------------------------------------------------------------------ #
# We deliberately want the Sterilization x Nutrient contrasts *within* one
# level of the third factor (TimeToIncubation for the polished data,
# Polishing for the unpolished data) and never across it. There are two ways
# to get there:
#
#   (a) Fit ONE three-way model on the full dataset (Sterilization * Nutrient
#       * ThirdFactor) and then ask emmeans for the Sterilization*Nutrient
#       contrasts conditional on each level of ThirdFactor
#       (emmeans(m, ~ Sterilization * Nutrient | ThirdFactor)). This pools
#       the residual variance and residual df across both levels of
#       ThirdFactor into one estimate -- which is only valid if the
#       within-cell variance is roughly the same across ThirdFactor's levels.
#
#   (b) Split the data by ThirdFactor first and fit a SEPARATE two-way model
#       (Sterilization * Nutrient) within each level. This estimates a
#       separate residual variance per level of ThirdFactor, at the cost of
#       fewer residual df per model.
#
# So "is it two-way or three-way ANOVA" depends on which of these you do --
# but it is NOT a free choice here: we checked below whether pooling the
# residual variance across ThirdFactor (option a) is defensible, via
# Levene's test on Surf_roughness ~ ThirdFactor alone, ignoring Sterilization/
# Nutrient. For BOTH datasets the variances across ThirdFactor's levels
# differ far too much for pooling to be valid (see diagnostics below:
# ~460-fold variance difference between Before/After in the polished data,
# and a significant Levene test, p = 0.0012, for Polishing in the unpolished
# data). We therefore use approach (b) throughout: a genuine TWO-way,
# fixed-effects ANOVA (Sterilization * Nutrient) fit separately within each
# level of the third factor. (surf_rough_60_df never had a third factor to
# begin with, so this question doesn't arise there.)

car::leveneTest(Surf_roughness ~ TimeToIncubation, data = surf_rough_35_pol_df, center = median)
surf_rough_35_pol_df |> dplyr::group_by(TimeToIncubation) |> dplyr::summarise(var = var(Surf_roughness), sd = sd(Surf_roughness))

car::leveneTest(Surf_roughness ~ Polishing, data = surf_rough_35_unpol_df, center = median)
surf_rough_35_unpol_df |> dplyr::group_by(Polishing) |> dplyr::summarise(var = var(Surf_roughness), sd = sd(Surf_roughness))

# ------------------------------------------------------------------------ #
# Planned contrasts
# ------------------------------------------------------------------------ #
# Contrast coefficients are positional: each coefficient applies to the
# corresponding row of the reference grid. Check the grid's factor levels and
# row order rather than assuming an ordering, since this can depend on the
# factors, package, and contrast method. A correctly sized coefficient vector
# may still describe a different comparison if its positions do not match the
# intended groups, so verify contrast labels and estimates against the
# underlying group means when defining or reviewing contrasts.
#
# Here, emmeans orders this Sterilization x Nutrient grid with Sterilization
# cycling fastest:
#   (C,-nutri), (DW,-nutri), (C,+nutri), (DW,+nutri)
# This helper matches factor-level names in the returned grid, avoiding an
# assumption about its positional order. Other packages or methods may build
# their reference grids differently. In particular, the ART contrasts in
# section 3 use a separate coefficient list aligned to art.con()'s grid; keep
# those definitions specific to their respective grids.
make_wanted_contrasts <- function(grid) {
  gdf <- as.data.frame(grid)
  idx <- function(steril, nutri) which(gdf$Sterilization == steril & gdf$Nutrient == nutri)
  mk <- function(s1, n1, s2, n2) {
    v <- numeric(nrow(gdf))
    v[idx(s1, n1)] <-  1
    v[idx(s2, n2)] <- -1
    v
  }
  list(
    "(C,-nutri) - (C,+nutri)"   = mk("C", "-nutri", "C",  "+nutri"),
    "(DW,-nutri) - (DW,+nutri)" = mk("DW","-nutri", "DW", "+nutri"),
    "(C,-nutri) - (DW,-nutri)"  = mk("C", "-nutri", "DW", "-nutri"),
    "(C,+nutri) - (DW,+nutri)"  = mk("C", "+nutri", "DW", "+nutri")
  )
}

diagnostic_dir <- file.path(results_dir, "model_diagnostics")
if (!dir.exists(diagnostic_dir)) dir.create(diagnostic_dir, recursive = TRUE)

# Helper: normality + homoscedasticity diagnostics for one Sterilization *
# Nutrient subset, with a 4-panel base diagnostic plot saved to PDF.
run_diagnostics <- function(dat, label, response = "Surf_roughness") {
  dat <- droplevels(dat)
  m <- lm(reformulate("Sterilization * Nutrient", response = response), data = dat)

  sh  <- shapiro.test(residuals(m))
  lev <- car::leveneTest(
    reformulate("interaction(Sterilization, Nutrient)", response = response),
    data = dat, center = median
  )

  pdf(file.path(diagnostic_dir, paste0(label, "_diagnostics.pdf")), width = 8, height = 8)
  par(mfrow = c(2, 2))
  plot(m)
  par(mfrow = c(1, 1))
  dev.off()

  cat(sprintf(
    "%-30s n=%2d  Shapiro W=%.3f p=%.4f | Levene F=%.3f p=%.4f\n",
    label, nrow(dat), sh$statistic, sh$p.value, lev[["F value"]][1], lev[["Pr(>F)"]][1]
  ))

  list(model = m, shapiro = sh, levene = lev)
}

fit_and_contrast <- function(dat, label, response = "Surf_roughness") {
  dat <- droplevels(dat)
  m <- lm(reformulate("Sterilization * Nutrient", response = response), data = dat)
  aov_tbl <- car::Anova(m, type = 3)
  emm <- emmeans::emmeans(m, ~ Sterilization * Nutrient)
  con <- emmeans::contrast(emm, method = make_wanted_contrasts(emm), adjust = "none")
  con_df <- as.data.frame(con)
  con_df$Subset <- label
  list(model = m, anova = aov_tbl, contrasts = con_df)
}

# ------------------------------------------------------------------------ #
# 1. surf_rough_35_pol_df -- diagnostics, raw scale, split by TimeToIncubation
# ------------------------------------------------------------------------ #
cat("\n--- 35d polished: RAW scale ---\n")
diag_before <- run_diagnostics(filter(surf_rough_35_pol_df, TimeToIncubation == "Before"), "35pol_Before_raw")
diag_after  <- run_diagnostics(filter(surf_rough_35_pol_df, TimeToIncubation == "After"),  "35pol_After_raw")

# NOTE on log-transformation for this dataset: the "Before" subset contains
# negative and near-zero Sq values (range -0.028 to 0.054 um; these are tiny
# roughness readings on a freshly polished, still essentially flat, coupon,
# where measurement noise straddles zero). log(x) is undefined for x <= 0, so
# a log-transform simply cannot be applied to "Before" -- this isn't a
# modeling choice, it's a hard mathematical constraint of this subset. We
# only attempt log10 on "After" (all values > 0):
cat("\n--- 35d polished, After: LOG10 scale (Before cannot be log-transformed: contains values <= 0) ---\n")
after_log <- filter(surf_rough_35_pol_df, TimeToIncubation == "After") |>
  droplevels() |>
  dplyr::mutate(Surf_roughness = log10(Surf_roughness))
diag_after_log <- run_diagnostics(after_log, "35pol_After_log10")

# Decision, 35d polished:
#   Before: Shapiro p = 0.19, Levene p = 0.63 -> both assumptions comfortably
#            met on the raw scale (as expected: no transform possible anyway).
#   After:  Shapiro p = 0.13, Levene p = 0.055 (raw) -- borderline but not
#            significant at alpha = 0.05. log10 "fixes" Levene (p = 0.12) but
#            breaks normality (Shapiro p = 0.017) -- net worse. So we keep
#            RAW scale for "After" too, flagging the borderline Levene result
#            as a caveat rather than grounds to abandon the parametric test
#            (see 02_corrosion_rate_statistics.R for the same reasoning pattern).
# -> Two-way fixed-effects ANOVA on the RAW scale for both "Before" and
#    "After", fit separately (see section 0b above for why "separately").

res_before <- fit_and_contrast(filter(surf_rough_35_pol_df, TimeToIncubation == "Before"), "Before")
res_after  <- fit_and_contrast(filter(surf_rough_35_pol_df, TimeToIncubation == "After"),  "After")

res_before$anova
res_after$anova

# One BH family across both subsets (8 tests): "Before" and "After" pose the
# same 4 questions at two points in the experiment, so we treat them as one
# family of comparisons for FDR control (mirrors 02_corrosion_rate_statistics.R).
pol_contrasts <- bind_rows(res_before$contrasts, res_after$contrasts) |>
  mutate(q_BH = p.adjust(p.value, method = "BH"))
pol_contrasts

# ------------------------------------------------------------------------ #
# 2. surf_rough_35_unpol_df -- diagnostics, raw scale, split by Polishing
# ------------------------------------------------------------------------ #
cat("\n--- 35d unpolished: RAW scale ---\n")
diag_unpol1 <- run_diagnostics(filter(surf_rough_35_unpol_df, Polishing == "Unpolished"), "35unpol_Unpolished_raw")
diag_unpol2 <- run_diagnostics(filter(surf_rough_35_unpol_df, Polishing == "Unpolished_side_of_polished"), "35unpol_UnpolishedSideOfPolished_raw")

# Decision, 35d unpolished:
#   Unpolished:                    Shapiro p = 0.61, Levene p = 0.29 -> fine.
#   Unpolished_side_of_polished:   Shapiro p = 0.38, Levene p = 0.56 -> fine.
# No transform needed for either subset.
# -> Two-way fixed-effects ANOVA on the RAW scale, fit separately per
#    Polishing level (again see section 0b: pooling across Polishing is not
#    defensible, Levene p = 0.0012 for Surf_roughness ~ Polishing alone).

res_unpol1 <- fit_and_contrast(filter(surf_rough_35_unpol_df, Polishing == "Unpolished"), "Unpolished")
res_unpol2 <- fit_and_contrast(filter(surf_rough_35_unpol_df, Polishing == "Unpolished_side_of_polished"), "Unpolished_side_of_polished")

res_unpol1$anova
res_unpol2$anova

unpol_contrasts <- bind_rows(res_unpol1$contrasts, res_unpol2$contrasts) |>
  mutate(q_BH = p.adjust(p.value, method = "BH"))
unpol_contrasts

# ------------------------------------------------------------------------ #
# 3. surf_rough_60_df -- diagnostics (only two factors here; no subsetting)
# ------------------------------------------------------------------------ #
cat("\n--- 60d: RAW scale ---\n")
diag_60 <- run_diagnostics(surf_rough_60_df, "60d_raw")

cat("\n--- 60d: LOG10 scale ---\n")
df_60_log <- surf_rough_60_df |> dplyr::mutate(Surf_roughness = log10(Surf_roughness))
diag_60_log <- run_diagnostics(df_60_log, "60d_log10")

# Decision, 60d:
#   Raw:   Shapiro p < 0.0001 -- strongly non-normal (driven by a right-hand
#          tail: most values sit in 1.9-3.3 um, but 4 coupons are much rougher,
#          up to 5.9 um). Levene p = 0.11 (not significant).
#   Log10: Shapiro p = 0.0094 -- STILL rejects normality; log10 does not fix
#          this (the tail is not a simple multiplicative effect here).
# Neither the raw nor the log-transformed residuals satisfy normality, so per
# the plan we fall back on ARTool (non-parametric, Aligned Rank Transform),
# which makes no distributional assumption about the response. This is the
# 60-day arm of the methods sentence at the top of this file.

art_60 <- ARTool::art(Surf_roughness ~ Sterilization * Nutrient, data = surf_rough_60_df)
anova_60_art <- anova(art_60)
anova_60_art

# These ART contrasts use coefficients aligned to art.con()'s concatenated
# interaction grid, whose order here is
#     1 (C,-nutri)  2 (C,+nutri)  3 (DW,-nutri)  4 (DW,+nutri)
# The emmeans helper above targets emmeans' grid instead. Keep definitions
# specific to each method and verify that the coefficient positions match the
# intended comparisons.
wanted_contrasts <- list(
  "(C,-nutri) - (C,+nutri)"   = c( 1, -1,  0,  0),
  "(DW,-nutri) - (DW,+nutri)" = c( 0,  0,  1, -1),
  "(C,-nutri) - (DW,-nutri)"  = c( 1,  0, -1,  0),
  "(C,+nutri) - (DW,+nutri)"  = c( 0,  1,  0, -1)
)

# art.con() builds that grid with factor(paste(f1, f2, sep = ",")), and
# factor() sorts by the locale's collation. The check below makes the expected
# order explicit and stops if the actual order differs from the coefficient
# definitions.
art_c_levels <- function(d) levels(factor(paste(d$Sterilization, d$Nutrient, sep = ",")))
expected_art_c_levels <- c("C,-nutri", "C,+nutri", "DW,-nutri", "DW,+nutri")
stopifnot(identical(art_c_levels(surf_rough_60_df), expected_art_c_levels))

con_60 <- ARTool::art.con(art_60, "Sterilization:Nutrient",
                          method = wanted_contrasts, adjust = "BH")
con_60_df <- as.data.frame(con_60)
con_60_df

# ------------------------------------------------------------------------ #
# 4. Export ANOVA + contrast tables to Excel
# ------------------------------------------------------------------------ #
# NOTE: this helper is written for car::Anova,
# whose rownames ARE the term names. ARTool's anova.art object already carries
# a `Term` column and has default integer rownames, so the ART_ANOVA_60d sheet
# gets a leading column of "1", "2", "3" also headed "Term", next to the real
# one. Values are unaffected.
anova_to_df <- function(x) cbind(Term = rownames(x), as.data.frame(x), row.names = NULL)

anova_tables <- list(
  ANOVA_35pol_Before = anova_to_df(res_before$anova),
  ANOVA_35pol_After  = anova_to_df(res_after$anova),
  ANOVA_35unpol_Unpolished     = anova_to_df(res_unpol1$anova),
  ANOVA_35unpol_UnpolSideOfPol = anova_to_df(res_unpol2$anova),
  ART_ANOVA_60d = anova_to_df(anova_60_art)
)

contrast_tables <- list(
  Contrasts_35pol_Before_After               = pol_contrasts,
  Contrasts_35unpol_Polishing_levels         = unpol_contrasts,
  Contrasts_60d_ART                          = con_60_df
)

writexl::write_xlsx(anova_tables,
                    file.path(results_dir, "ANOVA_tables_surface_roughness.xlsx"))
writexl::write_xlsx(contrast_tables,
                    file.path(results_dir, "Contrast_results_surface_roughness.xlsx"))

cat("\nDone. Diagnostic PDFs in", diagnostic_dir, "\n",
    "tables written to ANOVA_tables_surface_roughness.xlsx and",
    "Contrast_results_surface_roughness.xlsx in", results_dir, "\n")
