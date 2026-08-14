# ======================================================================== #
# Corrosion rate -- diagnostics, three-way ANOVA and ART / ART-C contrasts
#
# Implements the manuscript methods sentence:
#   "Corrosion rate and surface roughness were analyzed according to the
#    distributional characteristics of each dataset. Corrosion-rate data and
#    surface-roughness data from the 60-day incubation were analyzed using the
#    non-parametric ART procedure, with pairwise treatment contrasts obtained
#    using ART-C contrasts."
#
# Design: 2 factors x 2 levels, in two incubation durations.
#   Duration       35d (experiments E1 and E3) / 60d (experiment E2)
#   Sterilization  C (control) / DW (drinking water)
#   Nutrient       -nutri / +nutri
# Cells are tiny (n = 2-4) and unbalanced; see the comments in section 1.
#
# Sections 1-4 (raw/log diagnostics, Type III ANOVA, estimated marginal means)
# are the scale-selection and parametric cross-check that justify the choice;
# section 5 (ART + ART-C) is what the manuscript reports and what is exported.
#
# INPUT
#   data/physicochemical/corrosion_rate.csv
#       21 coupons, one row per coupon.
#
# OUTPUT (into results/statistics/)
#   corrosion_rate_ART_ANOVA_and_contrasts.xlsx
#       4 sheets: ART omnibus ANOVA and ART-C contrasts, per duration.
#
# Package versions used: car 3.1.5, emmeans 1.11.1, ARTool 0.11.2.
# dplyr is called namespace-qualified.
#
# Run from the repository root.
# ======================================================================== #

library(car)       # leveneTest(), Anova()
library(emmeans)   # estimated marginal means & custom contrasts
library(ARTool)    # Aligned Rank Transform (ART / ART-C)
library(openxlsx)  # needed to write a footnote below a table on the same sheet

results_dir <- "results/statistics"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------ #
# 0. Data
# ------------------------------------------------------------------------ #
# 21 mass-loss measurements with their three design factors, read from a tidy
# CSV and assembled into CR_df below.
#
# NOTE on factor level order: as.factor() sorts with the locale's collation.
# In the locale this analysis was run in, sort(c("-nutri", "+nutri")) puts
# "-nutri" first, and that ordering is assumed by the ART contrast
# coefficients in section 5. In a C/POSIX locale "+nutri" would sort first
# and those coefficients would address the wrong cells. See
# the contrast comments below.
CR_raw <- read.csv("data/physicochemical/corrosion_rate.csv",
                   stringsAsFactors = FALSE)

CR_df <- data.frame(
  Corrosion_rate = CR_raw$Corrosion_rate,
  Duration       = as.factor(CR_raw$Duration),
  Sterilization  = as.factor(CR_raw$Sterilization),
  Nutrient       = as.factor(CR_raw$Nutrient)
)
CR_df

# Cell sizes are tiny (n = 2-4) and unbalanced -> keep this in mind when
# reading the diagnostics: per-cell normality/variance tests have almost no
# power here, so we lean mainly on the *pooled* residual diagnostics from the
# fitted model rather than per-cell tests.
CR_df |> dplyr::count(Duration, Sterilization, Nutrient)

CR_df |>
  dplyr::group_by(Duration, Sterilization, Nutrient) |>
  dplyr::summarise(
    n = dplyr::n(),
    median_CR = median(Corrosion_rate, na.rm = TRUE),
    mean_CR = mean(Corrosion_rate, na.rm = TRUE),
    range_CR = paste(sprintf("%.1f", range(Corrosion_rate, na.rm = TRUE)), collapse = "-"),
    .groups = "drop"
  )

# ------------------------------------------------------------------------ #
# 1. Diagnostics on the RAW response
# ------------------------------------------------------------------------ #
# Use sum-to-zero contrasts so that Type II/III SS from car::Anova are
# interpretable for an unbalanced factorial design.
options(contrasts = c("contr.sum", "contr.poly"))

m_raw <- lm(Corrosion_rate ~ Duration * Sterilization * Nutrient, data = CR_df)

# 4-panel base diagnostic plot: Residuals vs Fitted, Q-Q, Scale-Location, Residuals vs Leverage
# (goes to the active graphics device -- this script was written to be stepped
# through interactively)
par(mfrow = c(2, 2))
plot(m_raw)
par(mfrow = c(1, 1))

# Formal normality test on the residuals (pooled, n = 21 -> Shapiro-Wilk is fine here)
shapiro.test(residuals(m_raw))

# Formal homogeneity-of-variance test across the 8 design cells
car::leveneTest(Corrosion_rate ~ Duration * Sterilization * Nutrient, data = CR_df)

# ------------------------------------------------------------------------ #
# 2. Diagnostics on the LOG-transformed response
# ------------------------------------------------------------------------ #
CR_df <- CR_df |> dplyr::mutate(log_Corrosion_rate = log(Corrosion_rate))

m_log <- lm(log_Corrosion_rate ~ Duration * Sterilization * Nutrient, data = CR_df)

par(mfrow = c(2, 2))
plot(m_log)
par(mfrow = c(1, 1))

shapiro.test(residuals(m_log))
car::leveneTest(log_Corrosion_rate ~ Duration * Sterilization * Nutrient, data = CR_df)

# ------------------------------------------------------------------------ #
# How to read the diagnostics (recap)
# ------------------------------------------------------------------------ #
# Residuals vs Fitted:
#   - Want a random scatter around the horizontal 0 line, no funnel/fan shape,
#     no curvature. A funnel (spread of residuals growing/shrinking with the
#     fitted value) signals heteroscedasticity; curvature signals a missing
#     term or the wrong scale for the response.
# Q-Q plot:
#   - Points should fall close to the diagonal line. Systematic curvature at
#     the tails (S- or banana-shape) indicates skew or heavy tails, i.e. a
#     departure from normality.
# Scale-Location (sqrt|standardized residuals| vs fitted):
#   - Want a roughly flat, horizontal band. An upward/downward trend again
#     signals heteroscedasticity (variance changing with the mean) -- very
#     common for rate/concentration data such as corrosion rate, and usually
#     fixed by a log transform.
# Residuals vs Leverage (Cook's distance):
#   - Flags influential points. With n as small as 2-4 per cell, a single
#     observation can dominate a cell mean, so check whether any point has a
#     Cook's distance that stands out (contour lines on the plot) or unusually
#     high leverage.
# Shapiro-Wilk test:
#   - H0: residuals are normally distributed. p < 0.05 -> reject normality.
#     With n = 21 this test has decent (not huge) power; use it together with
#     the Q-Q plot rather than in isolation.
# Levene's test:
#   - H0: variances are equal across the 8 groups. p < 0.05 -> reject
#     homogeneity of variance. With cells of n = 2-4 this test itself has low
#     power, so also inspect the Scale-Location plot and the spread of raw
#     values per cell.
# Decision actually reached for this dataset:
#   Raw scale:  Shapiro-Wilk p = 0.17 (no evidence against normality);
#               Levene's p = 0.005 (variances differ across the 8 cells).
#   Log scale:  Shapiro-Wilk p = 8e-05 and Levene's p = 5.5e-08 -- BOTH worse
#               than the raw scale. A Box-Cox profile on m_raw confirms this:
#               optimal lambda ~ 0.91 with a 95% CI of roughly [0.4, 1.4],
#               which comfortably contains lambda = 1 (no transform) and
#               excludes lambda = 0 (log). So log is the wrong transform here.
# -> We proceed on the RAW scale (m_raw) below.
#
# Why the Levene flag is treated as a caveat rather than a hard stop:
#   The design has only n = 21 observations spread over 8 cells, several of
#   them with just n = 2. A "variance" computed from 2 points is extremely
#   unstable (its own sampling distribution is very wide), so a single
#   unusual pair of replicates in one cell can swing Levene's test without
#   reflecting a real, generalizable difference in variability. With so few
#   residual degrees of freedom overall (df = 13), no formal test of
#   variance homogeneity has much power to distinguish "truly unequal" from
#   "noisy estimate from 2-4 points" -- so we note the significant Levene
#   result as a limitation of the analysis rather than as grounds to abandon
#   the (otherwise reasonably normal) raw-scale ANOVA.

CR_df |>
  dplyr::group_by(Duration, Sterilization, Nutrient) |>
  dplyr::summarise(mean = mean(Corrosion_rate), sd = sd(Corrosion_rate), .groups = "drop")

# ------------------------------------------------------------------------ #
# 3. Three-way fixed-effects ANOVA (RAW scale -- chosen scale, see decision above)
# ------------------------------------------------------------------------ #
# Type III sums of squares because the design is unbalanced and we care about
# the higher-order (interaction) terms as well as main effects.
# NOTE: with only n = 21 total and 13 residual df spread across an unbalanced
# 2x2x2 design, this ANOVA has limited power, especially for the three-way
# interaction and for cells with n = 2. Treat non-significant interaction/
# main-effect terms as "not detected" rather than "absent", and lean on the
# planned contrasts in section 5 rather than over-interpreting the omnibus
# table.
anova_raw <- car::Anova(m_raw, type = 3)
anova_raw

# ------------------------------------------------------------------------ #
# 4. Estimated marginal means per Duration subset
# ------------------------------------------------------------------------ #
# Cell means for every Sterilization x Nutrient combination, separately per
# Duration, borrowing the pooled residual variance/df from the full 3-way
# model (m_raw, df = 13). With cell sizes as low as n = 2 this is preferable
# to refitting two small separate models.
emm <- emmeans::emmeans(m_raw, ~ Sterilization * Nutrient | Duration)
emm

# ------------------------------------------------------------------------ #
# 5. Non-parametric analysis: Aligned Rank Transform (ART) ANOVA + ART-C contrasts
# ------------------------------------------------------------------------ #
# This is the section the manuscript reports for corrosion rate. Given n = 9
# (35d) / n = 12 (60d) split across only 4 cells each (as small as n = 2),
# rank-based inference is itself on shaky ground too (very coarse ranks, ties
# are possible, few residual df), but it makes no distributional assumption
# about the response. As in the phys-chem ART_ANOVA workflow, fit one ART
# model per Duration subset (Sterilization * Nutrient) and use art.con() with
# the ART-C contrast procedure for the same 4 planned contrasts, adjusted with
# BH within each subset.

# Planned contrasts for the ART-C procedure. The coefficients are positional,
# written for the level order art.con() uses internally: it pastes the two
# factor columns into one character column ("C,-nutri", "DW,+nutri", ...) and
# calls factor() on it, so the levels are the sorted concatenated strings --
# Sterilization outermost, Nutrient cycling fastest (ARTool 0.11.2,
# R/art.con.internal.R):
#       1 (C,-nutri)   2 (C,+nutri)   3 (DW,-nutri)   4 (DW,+nutri)
#
# This is NOT the order emmeans uses for `~ Sterilization * Nutrient` on a
# plain lm, where the first-named factor cycles fastest:
#       1 (C,-nutri)   2 (DW,-nutri)   3 (C,+nutri)   4 (DW,+nutri)
#
# The two differ by a 2 <-> 3 transposition. Both accept length-4 vectors, so
# handing these coefficients to emmeans::contrast() would not error -- it would
# silently trade the within-Sterilization and between-Sterilization tests under
# each other's labels. Use this list with art.con() only.
#
# BH (Benjamini-Hochberg / FDR) is applied within each Duration subset. These
# four contrasts are not all pairwise comparisons of one family, so Tukey and
# Dunnett do not apply; "mvt" is a defensible alternative if exact
# multivariate-t adjustment for their correlation is wanted. Given the small
# residual df, inspect the unadjusted p-values as well rather than reading
# hard significance from values near 0.05.
wanted_contrasts <- list(
  "(C,-nutri) - (C,+nutri)"   = c( 1, -1,  0,  0),
  "(DW,-nutri) - (DW,+nutri)" = c( 0,  0,  1, -1),
  "(C,-nutri) - (DW,-nutri)"  = c( 1,  0, -1,  0),
  "(C,+nutri) - (DW,+nutri)"  = c( 0,  1,  0, -1)
)

CR_35 <- CR_df |> dplyr::filter(Duration == "35d") |> droplevels()
CR_60 <- CR_df |> dplyr::filter(Duration == "60d") |> droplevels()

art_35 <- ARTool::art(Corrosion_rate ~ Sterilization * Nutrient, data = CR_35)
art_60 <- ARTool::art(Corrosion_rate ~ Sterilization * Nutrient, data = CR_60)

# art.con() builds that grid with factor(paste(f1, f2, sep = ",")), and
# factor() sorts by the locale's collation. These coefficients assume
# "-nutri" sorts before "+nutri", which is true in the locale this analysis
# ran in but NOT in a C/POSIX locale. Fail loudly rather than silently
# mislabel (note: containers/r_inext.def sets LC_ALL=C).
art_c_levels <- function(d) levels(factor(paste(d$Sterilization, d$Nutrient, sep = ",")))
expected_art_c_levels <- c("C,-nutri", "C,+nutri", "DW,-nutri", "DW,+nutri")
stopifnot(identical(art_c_levels(CR_35), expected_art_c_levels),
          identical(art_c_levels(CR_60), expected_art_c_levels))

# ART omnibus tests (Type III, analogous to anova_raw's subset-specific terms)
anova(art_35)
anova(art_60)

# ART-C planned contrasts, same coefficients/labels as wanted_contrasts above
con_35 <- ARTool::art.con(art_35, "Sterilization:Nutrient", method = wanted_contrasts, adjust = "BH")
con_60 <- ARTool::art.con(art_60, "Sterilization:Nutrient", method = wanted_contrasts, adjust = "BH")
con_35
con_60

# ------------------------------------------------------------------------ #
# 6. Export omnibus ANOVA + ART contrast tables to Excel
# ------------------------------------------------------------------------ #
# as.data.frame() on an emmGrid/summary_emm gives the printed table (estimate,
# SE, df, t/p-values, CI) as a plain data frame; the "mesg" attribute carries
# the p-value adjustment note (e.g. "P value adjustment: BH method for 4
# tests"), which we re-attach as a footnote row below each contrast table.
# openxlsx (rather than writexl) is used here because we need to write a
# second block of text below the data on the same sheet.

# Helper: turn an "anova"/data.frame object (car::Anova or ARTool::anova)
# into a plain data frame with the term names as an explicit column, since
# rownames are otherwise lost when written to Excel.
#
# NOTE: this is written for car::Anova, whose
# rownames ARE the term names. ARTool's anova.art object already carries a
# `Term` column and has default integer rownames, so the ART sheets get a
# leading column of "1", "2", "3" also headed "Term", next to the real one.
# Values are unaffected.
anova_to_df <- function(x) {
  df <- as.data.frame(x)
  cbind(Term = rownames(x), df, row.names = NULL)
}

anova_tables <- list(
  ART_ANOVA_35d      = anova_to_df(anova(art_35)),
  ART_ANOVA_60d      = anova_to_df(anova(art_60))
)

contrast_tables <- list(
  ART_contrasts_35d = as.data.frame(con_35),
  ART_contrasts_60d = as.data.frame(con_60)
)

wb <- openxlsx::createWorkbook()

# ANOVA sheets: just the table, no adjustment footnote (no multiplicity
# correction is applied to the omnibus terms).
for (sheet_name in names(anova_tables)) {
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, anova_tables[[sheet_name]])
  openxlsx::setColWidths(
    wb, sheet_name,
    cols = seq_len(ncol(anova_tables[[sheet_name]])), widths = "auto"
  )
}

# Contrast sheets: table + BH-adjustment footnote pulled from the "mesg"
# attribute, written two rows below the table.
for (sheet_name in names(contrast_tables)) {
  tbl <- contrast_tables[[sheet_name]]
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, tbl)
  # "auto" sizes columns to their widest cell (here, the long contrast
  # labels in the first column), so nothing gets truncated on open.
  openxlsx::setColWidths(wb, sheet_name, cols = seq_len(ncol(tbl)), widths = "auto")

  note <- attr(tbl, "mesg")
  if (!is.null(note)) {
    openxlsx::writeData(
      wb, sheet_name, paste(note, collapse = "; "),
      startRow = nrow(tbl) + 3, startCol = 1
    )
  }
}

openxlsx::saveWorkbook(
  wb,
  file.path(results_dir, "corrosion_rate_ART_ANOVA_and_contrasts.xlsx"),
  overwrite = TRUE
)

cat("\nWrote corrosion_rate_ART_ANOVA_and_contrasts.xlsx into", results_dir, "\n")
