# ======================================================================== #
# Physicochemical endpoints -- ART ANOVA, batch/time reproducibility,
# Fisher-combined treatment contrasts and meta-analytic effect sizes
#
# Implements the manuscript methods paragraph:
#   "Physicochemical endpoints were checked for normality and homogeneity of
#    variance by visual inspection of quantile-quantile and residuals-vs-fitted
#    diagnostic plots; several endpoints showed departures from normality and
#    non-constant variance, therefore, treatment comparisons were assessed
#    using the aligned rank transform (ART) procedure (ARTool R package
#    v0.11.2) [...] For each of the three independent experimental batches an
#    ART model (Sterilization x Nutrient) was fit and pairwise contrasts of
#    treatment effects were extracted via ART-C contrasts (two-tailed,
#    alpha = 0.05). Per-batch p-values for each contrast were combined across
#    the three batches using Fisher's combined probability test with metap
#    (v1.13), and per-batch effect sizes (r, derived from the contrast
#    t-statistic) were pooled using a Fisher-z fixed-effect meta-analysis with
#    metafor (v4.8.0). Separately, consistency of results across the two 35-day
#    replicate batches, and across the two incubation durations, was assessed
#    with additional ART models [...] For total Fe, residual diagnostics
#    differed by incubation duration: the non-parametric ART approach used
#    throughout the panel remained appropriate for the 35-day subset, while
#    log10-transformation resolved non-normality and heteroscedasticity of the
#    60-day subset, where a fixed-effects ANOVA on log-transformed
#    concentrations was used instead implementing the car package (v3.1.5),
#    with Tukey-adjusted pairwise contrasts computed via estimated marginal
#    means (emmeans v1.11.1). Fe2+ and Fe3+ concentrations were
#    log-transformed at both durations and were analyzed with two-way fixed
#    effects ANOVA."
#
# Design note (see data/physicochemical/README.md): Batch and Time are aliased.
# Experiments #1 and #3 ran 35 days, experiment #2 ran 60 days. E1 and E3 both
# ran 35 days but used water aliquots drawn at different times -- once while
# the building was quiet (exam period), once during active semester -- which is
# what the ICC / batch-effect part of this script quantifies.
#
# INPUT
#   data/physicochemical/{redox_potential,dissolved_oxygen,total_fe,fe2,fe3,
#     total_mn,sulphate,nitrate,nitrite,total_hardness,chloride,turbidity,
#     conductivity,ph}.csv
#       one tidy CSV per endpoint (14 files); response variable in the first
#       column. See data/physicochemical/README.md.
#   data/physicochemical/ph_day0.csv
#       pH with the Day-0 -nutri DW reference measurement, for section 8b.
#
# OUTPUT (into results/statistics/physicochemical/)
#   Supplementary_Data_physchem.xlsx        individual data + summaries
#   physicochemical_endpoints_long.csv      same, flat
#   Assumption_diagnostics.xlsx             Shapiro/Levene tables per subset
#   QQ_diagnostics_*.png, Residuals_by_treatment_*.png
#   Physchem_descriptive_stats_with_units.png (+ its plotting data CSV)
#   model_diagnostics_35d/<endpoint>_diagnostics.pdf
#   <endpoint>/<endpoint>_ART_ANOVA_tables.xlsx   per-endpoint workbook
#   ART_ANOVA_random-effects_VarCorr.xlsx
#   pH_Day0_ART_ANOVA.xlsx                  pH Day-0 reference contrast
#   ART_ANOVA_Executive_Summaries_FDR.xlsx
#   Figure2_3_p_values.xlsx
# OUTPUT (into results/figures/)
#   FigS2_heatmap_treatment_contrasts.png       (+ plotting data CSV)
#   FigS3_heatmap_Batch-pair_comparisons.png    (+ plotting data CSV)
#   Heatmap_TimePairs_qBH.png                   (+ plotting data CSV)
#   Heatmap_dual_batch-pair_contrasts_q_and_r.png    (+ plotting data CSV)
#   Heatmap_dual_time-pair_contrasts_q_and_r.png     (+ plotting data CSV)
#   Heatmap_ANOVA_35d_SxNxBatch_p_and_eta2.png       (+ plotting data CSV)
#   Heatmap_ANOVA_all-batches_SxNxTime_p_and_eta2.png (+ plotting data CSV)
#
# Run from the repository root. Runs on a workstation, not on the cluster.
# ======================================================================== #

library(ARTool)      # 0.11.2 -- art(), art.con() (ART-C contrasts)
library(car)         # 3.1.5  -- Anova(type = 3), leveneTest()
library(emmeans)     # 1.11.1 -- estimated marginal means for the lm() branch
library(metap)       # 1.13   -- sumlog(), Fisher's combined probability test
library(metafor)     # 4.8.0  -- escalc()/rma(), Fisher-z meta-analysis
library(lme4)        # lmer() for the descriptive batch-variance / ICC model
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(forcats)
library(stringr)
library(ggplot2)
library(ggnewscale)
library(writexl)
# MuMIn::r.squaredGLMM() and scales::squish() are called with :: only.

# repo-relative paths -- run this script with the repo root as working directory
physchem_dir <- "data/physicochemical"

base_path <- "results/statistics/physicochemical"
fig_dir   <- "results/figures"

dir.create(base_path, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,   recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------ #
# 1. Load the measurements and rebuild the per-endpoint data frames
# ------------------------------------------------------------------------ #
# The 14 endpoint data frames are built from data/physicochemical/*.csv,
# keeping the source row order and response-column names.
#
# The order of df_list is NOT alphabetical -- it is the order the analysis
# uses, and it is load-bearing: dir_names (derived from boxplot_ylabs)
# is indexed positionally in the main loop. So the file list is given
# explicitly rather than by Sys.glob().
endpoint_files <- c(
  ORP          = "redox_potential.csv",
  DO           = "dissolved_oxygen.csv",
  Fe           = "total_fe.csv",
  Fe2          = "fe2.csv",
  Fe3          = "fe3.csv",
  Mn           = "total_mn.csv",
  SO4          = "sulphate.csv",
  NO3          = "nitrate.csv",
  NO2          = "nitrite.csv",
  TH           = "total_hardness.csv",
  Cl           = "chloride.csv",
  Turbidity    = "turbidity.csv",
  Conductivity = "conductivity.csv",
  pH           = "ph.csv"
)

stopifnot(all(file.exists(file.path(physchem_dir, endpoint_files))))

df_list <- lapply(endpoint_files, function(fname) {
  d <- read.csv(file.path(physchem_dir, fname), stringsAsFactors = FALSE,
                check.names = FALSE, fileEncoding = "UTF-8")
  # response first, then the design factors; Replicate is an index only and is
  # not part of any model
  data.frame(
    d[, 1, drop = FALSE],
    Sterilization = factor(d$Sterilization, levels = c("Ctrl", "DW")),
    Nutrient      = factor(d$Nutrient,      levels = c("Nutri-", "Nutri+")),
    Time          = factor(d$Time,          levels = c("35d", "60d")),
    Batch         = factor(d$Batch,         levels = sort(unique(d$Batch))),
    check.names   = FALSE
  )
})

stopifnot(
  vapply(df_list, nrow, integer(1)) ==
    c(ORP = 36L, DO = 36L, Fe = 36L, Fe2 = 24L, Fe3 = 24L, Mn = 36L, SO4 = 36L,
      NO3 = 36L, NO2 = 36L, TH = 36L, Cl = 36L, Turbidity = 36L,
      Conductivity = 36L, pH = 36L)
)

# Relabel the factor levels to match the rest of the manuscript figures.
df_list <- lapply(df_list, function(d) {
  d %>%
    mutate(
      Sterilization = fct_recode(Sterilization, "C" = "Ctrl"),   # DW stays "DW"
      Nutrient      = fct_recode(Nutrient,
                                 "-nutri" = "Nutri-",
                                 "+nutri" = "Nutri+")
    )
})

# pH with the Day-0 reference (section 8b). Per batch: one pre-incubation,
# buffer-free measurement of -nutri DW (Time "0d"), plus the end-of-incubation
# C +nutri, DW -nutri and DW +nutri rows, which are the same values as in
# ph.csv. C -nutri is not part of this design.
pH_d0_df <- read.csv(file.path(physchem_dir, "ph_day0.csv"), stringsAsFactors = FALSE,
                     check.names = FALSE, fileEncoding = "UTF-8") %>%
  dplyr::transmute(
    pH,
    Sterilization = factor(Sterilization, levels = c("Ctrl", "DW")),
    Nutrient      = factor(Nutrient,      levels = c("Nutri-", "Nutri+")),
    Time          = factor(Time,          levels = c("0d", "35d", "60d")),
    Batch         = factor(Batch,         levels = c("#1", "#2", "#3"))
  ) %>%
  dplyr::mutate(
    Sterilization = fct_recode(Sterilization, "C" = "Ctrl"),
    Nutrient      = fct_recode(Nutrient, "-nutri" = "Nutri-", "+nutri" = "Nutri+"),
    Stage         = factor(ifelse(Time == "0d", "Day0", "End"), levels = c("Day0", "End"))
  )

stopifnot(nrow(pH_d0_df) == 30L, sum(pH_d0_df$Stage == "Day0") == 3L)

boxplot_ylabs <- list(`Redox potential` = "Redox potential (mV)",
                      `Dissolved oxygen` = "Dissolved oxygen (mg/L)",
                      `Total Fe` = "Total Fe (mg/L)",
                      `Fe(II)` = "Fe(II) mg/L",
                      `Fe(III)` = "Fe(III) mg/L",
                      `Total Mn` = "Total Mn (mg/L)",
                      `Sulphate` = "Sulphate (mg/L)",
                      `Nitrate` = "Nitrate (mg/L)",
                      `Nitrite` = "Nitrite (mg/L)",
                      `Total_hardness` = "Total hardness (mg/L)",
                      `Chloride` = "Chloride (mg/L)",
                      `Turbidity` = "Turbidity (NTU)",
                      # NOTE: the measured values are 232-4320, i.e. uS/cm,
                      # as endpoint_units below states. This label is used only
                      # to derive output directory names; nothing numeric
                      # depends on it.
                      `Conductivity` = "Conductivity (mS/cm)",
                      `pH` = "pH")

# Endpoints that need log10 transformation (tailed non-normality). Checked for
# all of Fe, Fe2, Fe3, Mn, SO4, NO3, NO2, Cl, TH, Turbidity, Conductivity, DO;
# only these reasonably benefit:
log_transform_endpoints <- c("Fe2", "Fe3")
# For total Fe only the 60-day dataset benefits from log-transformation, so it
# is applied per batch (batch #2 is the only 60-day run) in section 8C.
log_transform_60d_only <- c("Fe")

df_list_modeling <- imap(df_list, function(d, endpoint) {
  if (endpoint %in% log_transform_endpoints) {
    resp <- names(d)[1]
    d[[resp]] <- log10(d[[resp]])
  }
  d
})

# ------------------------------------------------------------------------ #
# 2. Normality / homoscedasticity diagnostics per incubation duration
# ------------------------------------------------------------------------ #
# This is the evidence behind the methods sentence "several endpoints showed
# departures from normality and non-constant variance".

run_single_batch_diagnostics <- function(df_list, duration_label) {
  imap_dfr(df_list, function(d, endpoint) {
    response <- names(d)[1]
    dat <- d %>% filter(Time == duration_label) %>% droplevels()
    if (nrow(dat) == 0L) return(tibble(Endpoint = endpoint, n = 0L, Comment = paste0("No ", duration_label, " data")))

    model <- lm(reformulate("Sterilization * Nutrient", response = response), data = dat)
    shapiro_result <- shapiro.test(residuals(model))
    levene_result <- car::leveneTest(
      reformulate("interaction(Sterilization, Nutrient)", response = response),
      data = dat, center = median
    )
    tibble(
      Endpoint = endpoint, n = nrow(dat), residual_df = df.residual(model),
      Shapiro_W = unname(shapiro_result$statistic), Shapiro_p = shapiro_result$p.value,
      Levene_F = levene_result[["F value"]][1], Levene_p = levene_result[["Pr(>F)"]][1],
      max_abs_standardized_residual = max(abs(rstandard(model)), na.rm = TRUE),
      Comment = NA_character_
    )
  })
}

assumption_diagnostics_60d      <- run_single_batch_diagnostics(df_list, "60d")
assumption_diagnostics_35d_Fe23 <- run_single_batch_diagnostics(df_list[c("Fe2", "Fe3")], "35d")
assumption_diagnostics_35d      <- run_single_batch_diagnostics(df_list, "35d")

assumption_diagnostics_60d_log      <- run_single_batch_diagnostics(df_list_modeling, "60d")
assumption_diagnostics_35d_Fe23_log <- run_single_batch_diagnostics(df_list_modeling[c("Fe2", "Fe3")], "35d")
assumption_diagnostics_35d_log      <- run_single_batch_diagnostics(df_list_modeling, "35d")

write_xlsx(
  list(
    raw_60d        = assumption_diagnostics_60d,
    raw_35d        = assumption_diagnostics_35d,
    raw_35d_Fe2Fe3 = assumption_diagnostics_35d_Fe23,
    log_60d        = assumption_diagnostics_60d_log,
    log_35d        = assumption_diagnostics_35d_log,
    log_35d_Fe2Fe3 = assumption_diagnostics_35d_Fe23_log
  ),
  file.path(base_path, "Assumption_diagnostics.xlsx")
)

build_diag_frames <- function(df_list, duration_label, diag_table) {
  frames <- imap_dfr(df_list, function(d, endpoint) {
    response <- names(d)[1]
    dat <- d %>% filter(Time == duration_label) %>% droplevels()
    model <- lm(reformulate("Sterilization * Nutrient", response = response), data = dat)
    tibble(Endpoint = endpoint, fitted = factor(round(fitted(model), 3)), std_resid = rstandard(model))
  })
  labels <- diag_table %>% transmute(Endpoint, label = sprintf("Shapiro p=%.3f\nLevene p=%.3f", Shapiro_p, Levene_p))
  frames %>% left_join(labels, by = "Endpoint")
}

make_qq_plot <- function(diag_frames, title) {
  ggplot(diag_frames, aes(sample = std_resid)) +
    stat_qq(size = 1.6, colour = "#2C7FB8") +
    stat_qq_line(colour = "grey40") +
    geom_text(data = distinct(diag_frames, Endpoint, label),
              aes(x = -Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = -0.05, vjust = 1.2, size = 2.6) +
    facet_wrap(vars(Endpoint), scales = "free") +
    labs(x = "Theoretical quantiles", y = "Standardized residuals", title = title) +
    theme_bw() + theme(strip.text = element_text(face = "bold"))
}

make_resid_plot <- function(diag_frames, title) {
  # only 4 distinct fitted values exist here (one per treatment cell) since the
  # model is purely categorical - this is really "residuals by treatment cell",
  # the visual equivalent of what Levene's test checks numerically
  ggplot(diag_frames, aes(x = fitted, y = std_resid)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_jitter(width = 0.05, size = 1.6, colour = "#D95F02") +
    facet_wrap(vars(Endpoint), scales = "free_y") +
    labs(x = "Fitted (treatment cell mean)", y = "Standardized residuals", title = title) +
    theme_bw() + theme(strip.text = element_text(face = "bold"),
                       axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
}

diag_frames_60d <- build_diag_frames(df_list_modeling, "60d", assumption_diagnostics_60d_log)
qq_60d    <- make_qq_plot(diag_frames_60d, "Normal Q-Q – 60-day models (Sterilization × Nutrient)")
resid_60d <- make_resid_plot(diag_frames_60d, "Residuals by treatment – 60-day models")

diag_frames_35d_Fe23 <- build_diag_frames(df_list_modeling[c("Fe2","Fe3")], "35d", assumption_diagnostics_35d_Fe23_log)
qq_35d_Fe23    <- make_qq_plot(diag_frames_35d_Fe23, "Normal Q-Q – 35-day Fe(II)/Fe(III) models")
resid_35d_Fe23 <- make_resid_plot(diag_frames_35d_Fe23, "Residuals by treatment – 35-day Fe(II)/Fe(III) models")

diag_frames_35d <- build_diag_frames(df_list_modeling, "35d", assumption_diagnostics_35d_log)
qq_35d    <- make_qq_plot(diag_frames_35d, "Normal Q-Q – 35-day models (Sterilization × Nutrient)")
resid_35d <- make_resid_plot(diag_frames_35d, "Residuals by treatment – 35-day models")

ggsave(file.path(base_path, "QQ_diagnostics_60d_log.png"), qq_60d, width = 10, height = 8, dpi = 300)
ggsave(file.path(base_path, "Residuals_by_treatment_60d_log.png"), resid_60d, width = 10, height = 8, dpi = 300)
ggsave(file.path(base_path, "QQ_diagnostics_35d_Fe23_log.png"), qq_35d_Fe23, width = 6, height = 3.5, dpi = 300)
ggsave(file.path(base_path, "Residuals_by_treatment_35d_Fe23_log.png"), resid_35d_Fe23, width = 6, height = 3.5, dpi = 300)
ggsave(file.path(base_path, "QQ_diagnostics_35d_log.png"), qq_35d, width = 10, height = 8, dpi = 300)
ggsave(file.path(base_path, "Residuals_by_treatment_35d_log.png"), resid_35d, width = 10, height = 8, dpi = 300)

# NOTE: the 35-day diagnostic frames are built twice -- first from df_list (raw
# scale), then again from df_list_modeling (log scale) into the SAME object
# names -- and saved twice, as "*_35d_log.png" and as "*_35d.png". Both files
# therefore hold the LOG-scale plot; the raw-scale build is
# overwritten before it was ever used. The dead raw-scale build has been removed
# (it had no effect on any output) and the duplicate save of the same objects
# under the "_35d" name has been dropped with it. Only Fe2/Fe3 are log-scaled at
# all, so for the other 12 endpoints the two files were identical anyway.

# ------------------------------------------------------------------------ #
# 3. Long-format supplementary data table
# ------------------------------------------------------------------------ #
treatment_levels <- c("C -nutri", "C +nutri", "DW -nutri", "DW +nutri")

# machine-readable per-endpoint frames feeding the long table
df_list_export <- lapply(df_list, function(d) {
  d %>%
    mutate(
      Treatment = factor(paste(Sterilization, Nutrient), levels = treatment_levels),
      Experiment = factor(Batch, levels = c("#1", "#2", "#3")),
      Duration_days = dplyr::recode(as.character(Time), "35d" = 35L, "60d" = 60L)
    ) %>%
    select(Experiment, Duration_days, Treatment, everything(),
           -Sterilization, -Nutrient, -Time, -Batch)
})

# NOTE: these two lookups are keyed "redox_pot" while df_list names that
# endpoint "ORP", so endpoint_labels[["ORP"]] and endpoint_units[["ORP"]] are
# both NA. Every redox row in physchem_long -- and therefore in the exported
# supplementary workbook and in the facet strip of the descriptive-statistics
# figure -- carries an empty Endpoint name and Unit; the `anyNA()` check below
# returns TRUE. Renaming the key to "ORP" resolves it, at the cost of changing
# the contents of a published supplementary table.
endpoint_labels <- c(
  redox_pot    = "ORP",
  DO           = "DO",
  Fe           = "Total Fe",
  Fe2          = "Fe²⁺",
  Fe3          = "Fe³⁺",
  Mn           = "Mn",
  SO4          = "SO₄²⁻",
  NO3          = "NO₃⁻",
  NO2          = "NO₂⁻",
  TH           = "Total hardness",
  Cl           = "Cl⁻",
  Turbidity    = "Turbidity",
  Conductivity = "Conductivity",
  pH           = "pH"
)

endpoint_units <- c(
  redox_pot    = "mV",
  DO           = "mg/L",
  Fe           = "mg/L",
  Fe2          = "mg/L",
  Fe3          = "mg/L",
  Mn           = "mg/L",
  SO4          = "mg/L",
  NO3          = "mg/L",
  NO2          = "mg/L",
  TH           = "mg/L as CaCO3",
  Cl           = "mg/L",
  Turbidity    = "NTU",
  Conductivity = "μS/cm",
  pH           = "dimensionless"
)

physchem_long <- imap_dfr(
  df_list_export,
  function(d, endpoint_id) {
    response_column <- setdiff(names(d), c("Experiment", "Duration_days", "Treatment"))
    if (length(response_column) != 1L) {
      stop("Expected one response column for endpoint ", endpoint_id,
           "; found: ", paste(response_column, collapse = ", "))
    }
    d %>%
      transmute(
        Experiment, Duration_days, Treatment,
        Endpoint_ID = endpoint_id,
        Endpoint = unname(endpoint_labels[endpoint_id]),
        Value = .data[[response_column]],
        Unit = unname(endpoint_units[endpoint_id])
      )
  }
) %>%
  # observation number within each endpoint x experiment x treatment group
  group_by(Endpoint_ID, Experiment, Duration_days, Treatment) %>%
  mutate(Observation = row_number()) %>%
  ungroup() %>%
  select(Experiment, Duration_days, Treatment, Observation,
         Endpoint_ID, Endpoint, Value, Unit)

anyNA(physchem_long$Endpoint)   # TRUE -- see the note above

physchem_summary_by_experiment <- physchem_long %>%
  group_by(Endpoint_ID, Endpoint, Unit, Experiment, Duration_days, Treatment) %>%
  summarise(
    n = sum(!is.na(Value)),
    Median = median(Value, na.rm = TRUE),
    Minimum = min(Value, na.rm = TRUE),
    Maximum = max(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(range = paste0(signif(Minimum, 4), "–", signif(Maximum, 4)))

physchem_summary_by_duration <- physchem_long %>%
  group_by(Endpoint_ID, Endpoint, Unit, Duration_days, Treatment) %>%
  summarise(
    n = sum(!is.na(Value)),
    n_experiments = n_distinct(Experiment),
    Median = median(Value, na.rm = TRUE),
    Minimum = min(Value, na.rm = TRUE),
    Maximum = max(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(range = paste0(signif(Minimum, 4), "–", signif(Maximum, 4)))

supplementary_workbook <- list(
  Individual_data = physchem_long,
  Summary_by_experiment = physchem_summary_by_experiment,
  Summary_by_duration = physchem_summary_by_duration,
  Data_dictionary = tibble(
    Variable = c("Experiment", "Duration_days", "Treatment", "Observation",
                 "Endpoint_ID", "Endpoint", "Value", "Unit"),
    Description = c(
      "Independent experimental run",
      "Incubation duration in days",
      "Sterilization and nutrient treatment combination",
      "Row identifier within an endpoint–experiment–treatment group",
      "Short endpoint identifier used in analysis",
      "Full physicochemical endpoint name",
      "Observed measurement",
      "Measurement unit"
    )
  )
)

write_xlsx(supplementary_workbook, file.path(base_path, "Supplementary_Data_physchem.xlsx"))
write.csv(physchem_long, file.path(base_path, "physicochemical_endpoints_long.csv"), row.names = FALSE)

# ------------------------------------------------------------------------ #
# 4. Descriptive statistics figure (not in the manuscript, kept as a useful
#    overview): jittered raw values per treatment, coloured by experiment,
#    with the median and interquartile range overlaid.
# ------------------------------------------------------------------------ #
endpoint_unit_lookup <- physchem_long %>% distinct(Endpoint, Unit) %>% deframe()

# Custom labeller: append unit in parentheses, but skip it for dimensionless
# endpoints (e.g. pH) to avoid an awkward "pH (dimensionless)"
label_with_units <- as_labeller(function(endpoint) {
  unit <- endpoint_unit_lookup[endpoint]
  ifelse(is.na(unit) | unit %in% c("dimensionless", ""), endpoint,
         paste0(endpoint, " (", unit, ")"))
})

descript_stats_endpoints_units <- ggplot(
  data = physchem_long,
  aes(x = Treatment, y = Value, colour = Experiment)
) +
  geom_point(position = position_jitter(width = 0.15), size = 2.2, alpha = 0.6) +
  facet_wrap(vars(Endpoint), scales = "free_y", labeller = label_with_units) +
  labs(x = NULL, y = "Observed value", colour = "Experiment") +
  scale_color_manual(values = c("#7ABECC", "#D6AD87", "#666666")) +
  scale_y_continuous(labels = function(l) ifelse(l == 0, "0", l)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold", size = 8)) +
  stat_summary(
    aes(group = Treatment), colour = "black", size = 0.4,
    fun = median,
    fun.min = function(x) quantile(x, 0.25),
    fun.max = function(x) quantile(x, 0.75),
    geom = "pointrange", position = position_nudge(x = 0.25)
  )

ggsave(file.path(base_path, "Physchem_descriptive_stats_with_units.png"),
       plot = descript_stats_endpoints_units, width = 9, height = 7, dpi = 300)
write.csv(physchem_long,
          file.path(base_path, "Physchem_descriptive_stats_with_units_plotting_data.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------ #
# 5. Reproducibility diagnostics on the 35-day subset (S x N x Batch)
# ------------------------------------------------------------------------ #
# NOTE (reported): this reuses the name `assumption_diagnostics_35d` from
# section 2, where it meant the per-duration Sterilization x Nutrient model.
# Here it means the 35-day Sterilization x Nutrient x Batch model. Both are
# used before the shadowing happens, so nothing downstream is affected, but the
# name carries two meanings in one script. Left as-is.
assumption_diagnostics_35d <- imap_dfr(
  df_list_modeling,
  function(d, endpoint) {
    response <- names(d)[1]
    dat <- d %>% filter(Time == "35d") %>% droplevels()

    if (n_distinct(dat$Batch) < 2L) {
      return(tibble(Endpoint = endpoint, n = nrow(dat),
                    Shapiro_p = NA_real_, Levene_p = NA_real_,
                    Comment = "Only one 35-day experiment available"))
    }

    formula_full <- reformulate("Sterilization * Nutrient * Batch", response = response)
    model <- lm(formula_full, data = dat)

    shapiro_result <- shapiro.test(residuals(model))
    levene_result <- car::leveneTest(
      reformulate("interaction(Sterilization, Nutrient, Batch)", response = response),
      data = dat, center = median
    )

    tibble(
      Endpoint = endpoint, n = nrow(dat),
      Shapiro_W = unname(shapiro_result$statistic),
      Shapiro_p = shapiro_result$p.value,
      Levene_F = levene_result[["F value"]][1],
      Levene_p = levene_result[["Pr(>F)"]][1],
      max_abs_standardized_residual = max(abs(rstandard(model)), na.rm = TRUE),
      Comment = NA_character_
    )
  }
)

diagnostic_dir <- file.path(base_path, "model_diagnostics_35d")
dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

iwalk(
  df_list,
  function(d, endpoint) {
    response <- names(d)[1]
    dat <- d %>% filter(Time == "35d") %>% droplevels()
    if (n_distinct(dat$Batch) < 2L) return(invisible(NULL))

    model <- lm(reformulate("Sterilization * Nutrient * Batch", response = response), data = dat)

    pdf(file.path(diagnostic_dir, paste0(endpoint, "_diagnostics.pdf")), width = 8, height = 8)
    par(mfrow = c(2, 2))
    plot(model)
    dev.off()
  }
)

# summary() on an ARTool `art` object prints the alignment validity check --
# the F values of ANOVAs on the aligned responses "not of interest" must all be
# ~0, otherwise the aligned-rank transform is not appropriate for that model.
# Inspect this before trusting any of the ART results below.
art_diagnostics_35d <- imap(
  df_list,
  function(d, endpoint) {
    response <- names(d)[1]
    dat <- d %>% filter(Time == "35d") %>% droplevels()
    if (n_distinct(dat$Batch) < 2L) return(NULL)

    model <- art(reformulate("Sterilization * Nutrient * Batch", response = response), data = dat)
    list(Endpoint = endpoint, Model = model, Summary = summary(model), ANOVA = anova(model))
  }
)

# ------------------------------------------------------------------------ #
# 6. Helper functions for the main loop
# ------------------------------------------------------------------------ #
dir_names <- tolower(gsub(pattern = " ", replacement = "_", names(boxplot_ylabs)))

CreateDirIfNotYet <- function(full_path) {
  if (!dir.exists(full_path)) dir.create(full_path, recursive = TRUE)
}

# Sum-to-zero contrasts, so that car::Anova(type = 3) is interpretable for the
# unbalanced (Fe2/Fe3) designs.
options(contrasts = c("contr.sum", "contr.poly"))

sig_code <- function(p) {
  ifelse(p < 0.0001, "****",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01,  "**",
                       ifelse(p < 0.05,  "*",
                              ifelse(p < 0.10,  ".", "")))))
}

# emmeans, lmer and art print the contrast delimiters distinctly ("," and " ");
# harmonize that with a helper function
normalize_contrast_label <- function(x) {
  vapply(x, function(one) {
    groups <- strsplit(one, " - ", fixed = TRUE)[[1]]
    groups <- vapply(groups, function(g) {
      inner <- sub("^\\((.*)\\)$", "\\1", g)        # strip outer parens
      toks  <- strsplit(inner, "[ ,]+")[[1]]         # split on space(s) OR comma(s)
      paste0("(", paste(toks, collapse = ","), ")")  # always rejoin with a comma
    }, character(1))
    paste(groups, collapse = " - ")
  }, character(1), USE.NAMES = FALSE)
}

# Contrast table for either branch: ART-C contrasts for an `art` model, or
# emmeans pairwise contrasts for the log-scale lm() branch. Both Tukey-adjusted
# within the model, per the methods text.
contrastArtTable <- function(model_obj, contrast = "Time") {
  if (!is.character(contrast) || length(contrast) != 1) {
    stop("contrast must be a single character string representing a factor or interaction, e.g. 'Time' or 'Sterilization:Nutrient'")
  }

  if (inherits(model_obj, "art")) {
    contrast_art <- art.con(model_obj, contrast, adjust = "tukey")
  } else {
    # log-scale lm() branch - convert "A:B" style spec to an emmeans formula "~ A * B"
    rhs <- gsub(":", " * ", contrast)
    contrast_art <- emmeans::contrast(
      emmeans::emmeans(model_obj, as.formula(paste("~", rhs))),
      method = "pairwise", adjust = "tukey"
    )
  }

  contrast_art_df <- as.data.frame(contrast_art)
  contrast_art_df$contrast <- normalize_contrast_label(contrast_art_df$contrast)

  # effect size r from the contrast t-statistic
  contrast_art_df$r_effect <- with(contrast_art_df, t.ratio / sqrt(t.ratio^2 + df))
  contrast_art_df$abs_r <- abs(contrast_art_df$r_effect)

  contrast_art_df$magnitude <- cut(
    contrast_art_df$abs_r,
    breaks = c(-Inf, 0.1, 0.3, 0.5, Inf),
    labels = c("trivial", "small", "medium", "large"),
    right = FALSE
  )

  contrast_art_df$signif <- sig_code(contrast_art_df$p.value)
  contrast_art_df
}

# partial eta squared from a type III ART / car::Anova table
eta2_from_anova <- function(anova_df, is_art) {
  anova_df <- as.data.frame(anova_df)
  if (is_art) {
    anova_df %>% dplyr::mutate(eta2_partial = `Sum Sq` / (`Sum Sq` + `Sum Sq.res`),
                               eta2_pct = round(100 * eta2_partial, 1))
  } else {
    anova_df <- tibble::rownames_to_column(anova_df, var = "Term")
    ss_resid <- anova_df$`Sum Sq`[anova_df$Term == "Residuals"]
    anova_df %>% dplyr::filter(Term != "Residuals") %>%
      dplyr::mutate(eta2_partial = `Sum Sq` / (`Sum Sq` + ss_resid),
                    eta2_pct = round(100 * eta2_partial, 1))
  }
}

# Intraclass correlation from the Sterilization*Nutrient + (1|Batch) model:
# the share of residual variance attributable to the experimental batch.
# E1 and E3 both ran 35 days but on water aliquots drawn at different times
# (quiet building during the exam period vs. active semester), so this is the
# quantity that says how much of the spread is "which experiment" rather than
# "which treatment".
icc_from_lmer <- function(m) {
  vc <- as.data.frame(VarCorr(m))
  # Sum of all random-effect variances (here it's just the Batch intercept)
  var_rand <- sum(vc$vcov[vc$grp != "Residual"], na.rm = TRUE)
  var_res  <- vc$vcov[vc$grp == "Residual"][1]
  var_rand / (var_rand + var_res)
}

# Fisher's combined probability test across batches (metap::sumlog)
fisher_combine_metap <- function(pvals) {
  p <- pvals[is.finite(pvals)]
  if (!length(p)) return(NA_real_)
  # keep p strictly within (0,1) to avoid edge warnings
  p <- pmin(pmax(p, .Machine$double.xmin), 1 - .Machine$double.eps)
  metap::sumlog(p)$p
}

# Fisher-z fixed-effect meta-analysis of the per-batch effect sizes (metafor).
# NOTE (reported): `n` is the number of observations in the whole batch (12),
# not the number entering the individual contrast, and `r` is a partial
# correlation derived from the contrast t-statistic rather than a zero-order
# correlation. escalc(measure = "ZCOR") assumes the latter, so vi = 1/(n - 3)
# is an approximation. It is applied identically to every endpoint and
# contrast, so the ranking of effects is unaffected, but the pooled CIs are
# somewhat optimistic.
meta_r_metafor <- function(r, n, model = c("FE","REML"), conf.level = 0.95) {
  model <- match.arg(model)
  ok <- is.finite(r) & is.finite(n) & (n > 3)
  if (!any(ok)) return(c(r = NA_real_, z = NA_real_, se = NA_real_, ci95 = NA_character_))

  # 1) Fisher-z effect sizes and sampling variances
  dat <- escalc(measure = "ZCOR", ri = r[ok], ni = n[ok])  # gives yi = z_i, vi = 1/(n_i - 3)
  # 2) Fit meta-analysis model (fixed-effect or random-effects REML)
  fit <- if (model == "FE") rma(yi, vi, data = dat, method = "FE")
         else               rma(yi, vi, data = dat, method = "REML")

  zF  <- as.numeric(fit$b)                  # pooled Fisher-z
  se  <- fit$se                             # SE on z scale
  a   <- (1 - conf.level)/2
  zcrit <- qnorm(1 - a)
  zlo <- zF - zcrit*se
  zhi <- zF + zcrit*se

  # back-transform to r
  rF  <- tanh(zF); rlo <- tanh(zlo); rhi <- tanh(zhi)
  c(r = rF, z = zF, se = se, ci95 = sprintf("[%.3f, %.3f]", rlo, rhi))
}

# Batch-only pairs (same S,N; #1 vs #3) out of the S*N*Batch contrasts
pick_batch_pairs <- function(contr_df) {
  # e.g., "(DW,-nutri,#1) - (DW,-nutri,#3)"
  pattern <- "^\\(([^,]+),([^,]+),#1\\) - \\(\\1,\\2,#3\\)$"
  contr_df %>% dplyr::filter(grepl(pattern, contrast))
}

# Time-only pairs (same S,N; 35d vs 60d) out of the S*N*Time contrasts
pick_time_pairs <- function(contr_df) {
  contr_df %>%
    tidyr::extract(
      contrast,
      into = c("S1","N1","T1","S2","N2","T2"),
      regex = "^\\(([^,]+),([^,]+),([^\\)]+)\\) - \\(([^,]+),([^,]+),([^\\)]+)\\)$",
      remove = FALSE
    ) %>%
    dplyr::filter(S1 == S2, N1 == N2, T1 != T2) %>%
    dplyr::select(-S1, -N1, -T1, -S2, -N2, -T2)
}

# ------------------------------------------------------------------------ #
# 7. Accumulators
# ------------------------------------------------------------------------ #
exec35_list  <- list()   # 35d S*N*Batch ANOVA tables
execAll_list <- list()   # all-batches S*N*Time ANOVA tables

snB_list <- list()       # S*N*Batch contrasts (art.con on m_35)
snT_list <- list()       # S*N*Time  contrasts (art.con on mAll)

snB_pairs_list <- list() # Batch-only pairs at fixed S,N
snT_pairs_list <- list() # Time-only pairs at fixed S,N

syn_list      <- list()  # per-parameter synthesis (Fisher p + meta r)
varCorr_list  <- list()  # lme4::VarCorr() output
icc_rows      <- list()  # batch-effect ICC / R2
sn35_list     <- list()  # 35d S*N contrasts, for the Fig. 3 annotations
long_all_list <- list()  # per-batch contrast tables, all endpoints

# ------------------------------------------------------------------------ #
# 8. Main analysis loop
# ------------------------------------------------------------------------ #
for (i in seq_along(df_list)) {

  current_data <- df_list[[i]]
  param_name <- names(df_list)[i]
  dir_name <- dir_names[i]
  response_var_name <- names(current_data)[1]  # First column is the response variable

  output_path <- file.path(base_path, dir_name)
  CreateDirIfNotYet(output_path)

  # --- Detect availability of Batch #1 and #3 for this endpoint
  has_1_3 <- with(current_data, all(c("#1","#3") %in% levels(Batch))) &&
    all(c("#1","#3") %in% unique(as.character(current_data$Batch[current_data$Time == "35d"])))

  # ---- A) 35-day, Batch-focused ANOVA: is there a batch effect? ----
  # Fe(II) and Fe(III) were not measured in experiment #1, so they skip this.
  #
  # NOTE: `response_var_summary` is created ONLY in this branch, while
  # sections B and C below append to it unconditionally. For
  # Fe2 and Fe3 the branch is skipped, so the workbook written for those two
  # endpoints still carries the ANOVA_SNB_effects / ANOVA_SNB_contrasts /
  # ANOVA_B_contrast / Random_effect_VarCorr / Mixed_effect_R_sq /
  # Mixed_effect_ICC sheets left over from the PREVIOUS endpoint in df_list
  # (total Fe, which immediately precedes Fe2). Those six sheets in
  # fe(ii)_ART_ANOVA_tables.xlsx and fe(iii)_ART_ANOVA_tables.xlsx describe
  # total Fe, not Fe(II)/Fe(III). Only the per-endpoint workbooks are affected;
  # every table that feeds the manuscript figures is built from the
  # accumulators below, which are correctly gated on has_1_3.
  if (has_1_3) {
    n_endpoints <- 14
    data_35d <- current_data %>% dplyr::filter(Time == "35d") %>% droplevels()
    use_log_35 <- param_name %in% log_transform_endpoints
    if (use_log_35) data_35d[[response_var_name]] <- log10(data_35d[[response_var_name]])

    m_35 <- if (use_log_35) {
      lm(reformulate("Sterilization * Nutrient * Batch", response = response_var_name), data = data_35d)
    } else {
      art(reformulate("Sterilization * Nutrient * Batch", response = response_var_name), data = data_35d)
    }

    ## --- 1) ANOVA table ---
    # Type III tests by default for ART, explicit type 3 for the parametric branch
    art_aov <- if (use_log_35) car::Anova(m_35, type = 3) else anova(m_35)
    art_aov_df <- as.data.frame(art_aov)
    # add ART eta^2 (rank scale): SS_effect / (SS_effect + SS_residual)
    art_aov_df <- eta2_from_anova(anova_df = art_aov_df, is_art = !use_log_35)

    response_var_summary <- list(ANOVA_SNB_effects = art_aov_df)

    # ---- A1) Unified 35d S*N*Batch ANOVA ----
    exec35_list[[param_name]] <- art_aov_df %>%
      dplyr::transmute(Parameter = param_name, Term,
                       p_value = `Pr(>F)`, eta2_partial = eta2_partial)

    contrast_art_SNB_df <- contrastArtTable(m_35, contrast = "Sterilization:Nutrient:Batch")
    # for Fig. 3
    contrast_art_SN_35d_df <- contrastArtTable(m_35, contrast = "Sterilization:Nutrient")
    sn35_list[[param_name]] <- contrast_art_SN_35d_df %>% dplyr::mutate(Parameter = param_name)

    response_var_summary$ANOVA_SNB_contrasts <- contrast_art_SNB_df

    # ---- A2) Unified S*N*Batch contrasts (on m_35) ----
    snB_list[[param_name]] <- contrast_art_SNB_df %>%
      dplyr::mutate(Parameter = param_name) %>%
      dplyr::select(Parameter, contrast, estimate, SE, df, t.ratio, p.value, r_effect, abs_r, magnitude)

    # Batch-only pairs at fixed S,N
    snB_pairs_list[[param_name]] <- pick_batch_pairs(contrast_art_SNB_df) %>%
      dplyr::mutate(Parameter = param_name) %>%
      dplyr::select(Parameter, contrast, estimate, SE, df, t.ratio, p.value, r_effect, abs_r, magnitude)

    contrast_art_batch_df <- contrastArtTable(m_35, contrast = "Batch")
    response_var_summary$ANOVA_B_contrast <- contrast_art_batch_df

    # --- descriptive size of the batch variance (ICC, on the raw scale) ---
    # Could be "~ Sterilization * Nutrient * Time + (1|Batch)" if any single
    # batch had run more than one duration, but Batch and Time are aliased
    # (batch #2 is the 60-day run, batches #1 and #3 are 35-day), so only the
    # 35-day subset supports a Batch random effect.
    m_lmer <- lmer(
      formula = as.formula(paste(response_var_name, "~ Sterilization * Nutrient + (1|Batch)")),
      data    = data_35d)

    batch_var <- as.data.frame(VarCorr(m_lmer))
    icc_val   <- icc_from_lmer(m_lmer)
    rsq_mat   <- MuMIn::r.squaredGLMM(m_lmer)      # 1x2 matrix for LMMs: R2m, R2c
    random_var <- (unname(rsq_mat[1, "R2c"]) - unname(rsq_mat[1, "R2m"]))

    response_var_summary$Random_effect_VarCorr <- batch_var
    response_var_summary$Mixed_effect_R_sq     <- as.data.frame(rsq_mat)
    response_var_summary$Mixed_effect_ICC      <- data.frame(ICC = icc_val)

    varCorr_list[[param_name]] <- batch_var
    icc_rows[[param_name]] <- data.frame(
      Parameter = param_name,
      ICC       = icc_val,
      R2m       = unname(rsq_mat[1, "R2m"]),
      R2c       = unname(rsq_mat[1, "R2c"]),
      `Var(Batch)` = random_var,
      check.names = FALSE
    )

  } else {
    n_endpoints <- 12
    message(sprintf("Skipping Section A) for %s: missing Batch #1 at 35d", param_name))
  }

  # ---- B) All batches: Sterilization * Nutrient * Time ----
  cat("B) All batches ART: Sterilization * Nutrient * Time ...\n")
  # Time and Batch are aliased: only batch #2 ran 60 days.

  use_log_all <- param_name %in% log_transform_endpoints
  data_all <- current_data
  if (use_log_all) data_all[[response_var_name]] <- log10(data_all[[response_var_name]])

  mAll <- if (use_log_all) {
    lm(reformulate("Sterilization * Nutrient * Time", response = response_var_name), data = data_all)
  } else {
    art(reformulate("Sterilization * Nutrient * Time", response = response_var_name), data = data_all)
  }
  aovAll <- if (use_log_all) car::Anova(mAll, type = 3) else anova(mAll)
  aovAll_df <- as.data.frame(aovAll)
  aovAll_df <- eta2_from_anova(anova_df = aovAll_df, is_art = !use_log_all)

  response_var_summary$ANOVA_SNT_effects <- aovAll_df

  # ---- B1) All-batches S*N*Time ANOVA ----
  execAll_list[[param_name]] <- aovAll_df %>%
    dplyr::transmute(Parameter = param_name, Term,
                     p_value = `Pr(>F)`, eta2_partial = eta2_partial)

  contrast_art_SN_df <- contrastArtTable(mAll, contrast = "Sterilization:Nutrient")
  response_var_summary$ANOVA_SN_contrasts <- contrast_art_SN_df

  contrast_art_SNT_df <- contrastArtTable(mAll, contrast = "Sterilization:Nutrient:Time")
  response_var_summary$ANOVA_SNT_contrasts <- contrast_art_SNT_df

  # ---- B2) S*N*Time contrasts (on mAll) ----
  snT_list[[param_name]] <- contrast_art_SNT_df %>%
    dplyr::mutate(Parameter = param_name) %>%
    dplyr::select(Parameter, contrast, estimate, SE, df, t.ratio, p.value, r_effect, abs_r, magnitude)

  # Time-only pairs at fixed S,N
  snT_pairs_list[[param_name]] <- pick_time_pairs(contrast_art_SNT_df) %>%
    dplyr::mutate(Parameter = param_name) %>%
    dplyr::select(Parameter, contrast, estimate, SE, df, t.ratio, p.value, r_effect, abs_r, magnitude)

  contrast_art_time_df <- contrastArtTable(mAll, contrast = "Time")
  response_var_summary$ANOVA_T_contrast <- contrast_art_time_df

  # ---- C) Per-batch ART + ART-C, then Fisher / meta-r synthesis ----
  cat("C) Within-batch ART + ART-C + Fisher/meta-r synthesis ...\n")
  batches <- unique(current_data$Batch)
  contr_by_batch <- list()

  for (b in batches) {
    db <- subset(current_data, Batch == b)
    if (nrow(db) == 0) next

    # Total Fe: ART for the 35-day batches, log10 + fixed-effects ANOVA for the
    # 60-day batch (#2), per the methods text. Fe2/Fe3: log10 everywhere.
    use_log_b <- (param_name %in% log_transform_endpoints) ||
                 (param_name %in% log_transform_60d_only && b == "#2")
    if (use_log_b) db[[response_var_name]] <- log10(db[[response_var_name]])

    mB <- if (use_log_b) {
      lm(reformulate("Sterilization * Nutrient", response = response_var_name), data = db)
    } else {
      art(reformulate("Sterilization * Nutrient", response = response_var_name), data = db)
    }

    # n for this batch, used to weight the Fisher-z meta-analysis
    n_batch <- nrow(db)
    aovB <- if (use_log_b) car::Anova(mB, type = 3) else anova(mB)
    aovB_df <- as.data.frame(aovB)
    aovB_eta <- eta2_from_anova(anova_df = aovB_df, is_art = !use_log_b)

    response_var_summary[[paste0("ANOVA_SN_effects_withinBatch", b)]] <- aovB_eta

    # entire S*N contrast matrix within this experimental batch
    consB <- contrastArtTable(mB, "Sterilization:Nutrient")
    consB$Batch <- b
    consB$n_obs <- n_batch   # attach n to every contrast row
    contr_by_batch[[as.character(b)]] <- consB
  }

  # Synthesis: combine the same contrast across batches
  if (length(contr_by_batch) > 0) {
    long <- bind_rows(contr_by_batch) %>%
      select(Batch, contrast, estimate, t.ratio, df, p.value, r_effect, n_obs)

    # NOTE (reported): the per-batch p-values entering Fisher's method are
    # already Tukey-adjusted inside contrastArtTable(). Fisher's combination
    # assumes p-values that are uniform under the null; Tukey-adjusted ones are
    # conservative, so p_Fisher here is conservative too. This is the
    # and what the manuscript reports; flagged rather than changed.
    synth <- long %>%
      group_by(contrast) %>%
      summarise(
        k        = n(),
        p_Fisher = fisher_combine_metap(p.value),
        .meta    = list(meta_r_metafor(r_effect, n = n_obs, model = "FE")),
        .groups  = "drop") %>%
      unnest_wider(.meta, names_sep = "_") %>%  # creates .meta_r, .meta_z, .meta_se, .meta_ci95
      rename(r_meta = .meta_r, z_meta = .meta_z, se_z = .meta_se, ci95_r = .meta_ci95) %>%
      select(-any_of(c("n_obs", "k", ".meta")))

    response_var_summary$ANOVA_SN_contrasts_withinBatches <- long
    response_var_summary$ANOVA_SN_contrasts_COMBINED <- synth

    write_xlsx(response_var_summary,
               file.path(output_path, paste0(dir_name, "_ART_ANOVA_tables.xlsx")))

    syn_list[[param_name]] <- synth %>%
      dplyr::transmute(Parameter = param_name, contrast, p_Fisher, r_meta)
  }

  long_all_list[[param_name]] <- long %>% dplyr::mutate(Parameter = param_name)

  cat("Done:", param_name, "\n")

  # NOTE (reported): this line sits inside the loop and prints once per
  # endpoint, and `n_endpoints` is a hardcoded 14/12 set by the has_1_3 branch
  # above rather than a count. Cosmetic only -- no result depends on it.
  cat(paste0("\nAll tasks completed for all ", n_endpoints, " endpoints.\n"))
}

# ------------------------------------------------------------------------ #
# 8b. pH only: Day-0 referenced contrast for -nutri DW
# ------------------------------------------------------------------------ #
# pH can shift simply from the autoclave procedure itself, so the standard
# "(C,-nutri) - (DW,-nutri)" contrast (autoclaved vs. raw drinking water,
# both without nutrients) confounds a sterilization-induced pH shift with
# any biological effect. For pH only, that cell is replaced by a contrast
# that IS interpretable: the single Day-0 reference measurement of -nutri DW
# (buffer-free, pre-incubation, one measurement per batch) vs. the -nutri DW
# samples at the end of incubation (35 d or 60 d, 3 replicates per batch).
# The +nutri samples keep the ordinary "(C,+nutri) - (DW,+nutri)" comparison
# (computed in the main loop, Parameter == "pH"), because the M9 buffer keeps
# pH stable regardless of sterilization. For every other endpoint the
# sterilized samples are compared with the non-sterilized controls, as
# described in the methods.
#
# To stay inside the same statistical framework used for every other
# endpoint x contrast cell (per-batch ART -> Fisher-combined p -> Fisher-z
# meta-analytic r -> one flat BH correction across the whole family), the
# Day-0 contrast goes through that identical per-batch ART + Fisher/meta-r
# pipeline, restricted to the DW -nutri cells. Each batch contributes
# n = 4 (1 Day-0 + 3 end-of-incubation).
cat("8b) pH: Day-0 referenced contrast for -nutri DW, per-batch ART + Fisher/meta-r ...\n")

pH_day0_contrast_label <- "(DW,-nutri,Day0) - (DW,-nutri)"
pH_invalid_contrast    <- "(C,-nutri) - (DW,-nutri)"

pH_d0_dw_nutri <- pH_d0_df %>%
  dplyr::filter(Sterilization == "DW", Nutrient == "-nutri") %>%
  droplevels()

pH_d0_by_batch <- list()
pH_d0_anova_by_batch <- list()
for (b in levels(pH_d0_dw_nutri$Batch)) {
  db <- subset(pH_d0_dw_nutri, Batch == b)
  mB_pH_d0 <- art(pH ~ Stage, data = db)

  pH_d0_anova_by_batch[[b]] <- as.data.frame(anova(mB_pH_d0)) %>%
    dplyr::mutate(Batch = b, .before = 1)

  con_pH_d0 <- contrastArtTable(mB_pH_d0, contrast = "Stage")
  con_pH_d0$contrast <- pH_day0_contrast_label   # "(Day0) - (End)" -> full label
  con_pH_d0$Batch    <- b
  con_pH_d0$n_obs    <- nrow(db)
  pH_d0_by_batch[[b]] <- con_pH_d0
}

pH_d0_long <- dplyr::bind_rows(pH_d0_by_batch) %>%
  dplyr::select(Batch, contrast, estimate, SE, df, t.ratio, p.value, r_effect, n_obs)

pH_d0_synth <- pH_d0_long %>%
  dplyr::group_by(contrast) %>%
  dplyr::summarise(
    k        = dplyr::n(),
    p_Fisher = fisher_combine_metap(p.value),
    .meta    = list(meta_r_metafor(r_effect, n = n_obs, model = "FE")),
    .groups  = "drop") %>%
  tidyr::unnest_wider(.meta, names_sep = "_") %>%
  dplyr::rename(r_meta = .meta_r, z_meta = .meta_z, se_z = .meta_se, ci95_r = .meta_ci95) %>%
  dplyr::mutate(Parameter = "pH", .before = 1)

# descriptive statistics for the Day-0 reference itself (one per batch)
pH_day0_descriptives <- pH_d0_df %>%
  dplyr::filter(Stage == "Day0") %>%
  dplyr::summarise(
    n         = dplyr::n(),
    median_pH = median(pH),
    min_pH    = min(pH),
    max_pH    = max(pH)) %>%
  dplyr::mutate(range_pH = paste0(min_pH, "–", max_pH))

write_xlsx(
  list(
    Batch_ART_ANOVA   = dplyr::bind_rows(pH_d0_anova_by_batch),
    Batch_contrasts   = pH_d0_long,
    Fisher_combined   = pH_d0_synth,
    Day0_descriptives = pH_day0_descriptives
  ),
  file.path(base_path, "pH_Day0_ART_ANOVA.xlsx")
)
                     
# ------------------------------------------------------------------------ #
# 9. Binding: executive tables across all endpoints
# ------------------------------------------------------------------------ #

# ---- A) ANOVA executive tables ----
exec35  <- dplyr::bind_rows(exec35_list)
execAll <- dplyr::bind_rows(execAll_list)

# Wide matrices (rows = Parameter, columns = Terms)
an35_p_mat   <- exec35 %>% select(Parameter, Term, p_value) %>%
  tidyr::pivot_wider(names_from = Term, values_from = p_value)
an35_eta_mat <- exec35 %>% select(Parameter, Term, eta2_partial) %>%
  tidyr::pivot_wider(names_from = Term, values_from = eta2_partial)

anAll_p_mat   <- execAll %>% select(Parameter, Term, p_value) %>%
  tidyr::pivot_wider(names_from = Term, values_from = p_value)
anAll_eta_mat <- execAll %>% select(Parameter, Term, eta2_partial) %>%
  tidyr::pivot_wider(names_from = Term, values_from = eta2_partial)

# lme4::VarCorr() output per endpoint, written separately
write_xlsx(varCorr_list, file.path(base_path, "ART_ANOVA_random-effects_VarCorr.xlsx"))

# Batch-effect size: ICC plus marginal/conditional R^2 (Nakagawa & Schielzeth),
# i.e. how much of the variance is "which experiment" rather than "which
# treatment". See the note at icc_from_lmer().
icc_all <- dplyr::bind_rows(icc_rows)

# ---- B) Contrast executive tables (rows = contrast, columns = Parameter) ----
snB_all <- dplyr::bind_rows(snB_list)
snT_all <- dplyr::bind_rows(snT_list)

snB_p_mat <- snB_all %>% select(contrast, Parameter, p.value) %>%
  tidyr::pivot_wider(names_from = Parameter, values_from = p.value)
snT_p_mat <- snT_all %>% select(contrast, Parameter, p.value) %>%
  tidyr::pivot_wider(names_from = Parameter, values_from = p.value)

snB_r_mat <- snB_all %>% select(contrast, Parameter, r_effect) %>%
  tidyr::pivot_wider(names_from = Parameter, values_from = r_effect)
snT_r_mat <- snT_all %>% select(contrast, Parameter, r_effect) %>%
  tidyr::pivot_wider(names_from = Parameter, values_from = r_effect)

# ---- C) Excerpts: Batch-only and Time-only pairs at fixed S,N ----
snB_pairs_p <- dplyr::bind_rows(snB_pairs_list) %>%
  select(contrast, Parameter, p.value) %>%
  tidyr::pivot_wider(names_from = Parameter, values_from = p.value)

snT_pairs_p <- dplyr::bind_rows(snT_pairs_list) %>%
  select(contrast, Parameter, p.value) %>%
  tidyr::pivot_wider(names_from = Parameter, values_from = p.value)

snB_pairs_r <- dplyr::bind_rows(snB_pairs_list) %>%
  select(contrast, Parameter, r_effect) %>%
  tidyr::pivot_wider(names_from = Parameter, values_from = r_effect)

snT_pairs_r <- dplyr::bind_rows(snT_pairs_list) %>%
  select(contrast, Parameter, r_effect) %>%
  tidyr::pivot_wider(names_from = Parameter, values_from = r_effect)

# Datasets for the batch- and time-reproducibility heatmaps: "does batch matter,
# holding treatment fixed" / "does duration matter, holding treatment fixed".
# One flat Benjamini-Hochberg correction per family.
BatchPairs_qBH <- dplyr::bind_rows(snB_pairs_list) %>%
  dplyr::mutate(q_BH = p.adjust(p.value, method = "BH")) %>%
  select(contrast, Parameter, q_BH) %>%
  tidyr::pivot_wider(names_from = Parameter, values_from = q_BH)

TimePairs_qBH <- dplyr::bind_rows(snT_pairs_list) %>%
  dplyr::mutate(q_BH = p.adjust(p.value, method = "BH")) %>%
  select(contrast, Parameter, q_BH) %>%
  tidyr::pivot_wider(names_from = Parameter, values_from = q_BH)

# ------------------------------------------------------------------------ #
# 10. Synthesis across endpoints
# ------------------------------------------------------------------------ #
syn_all <- dplyr::bind_rows(syn_list)

# The four planned treatment contrasts, by label
wanted_contrasts <- c(
  "(C,-nutri) - (C,+nutri)",
  "(DW,-nutri) - (DW,+nutri)",
  "(C,-nutri) - (DW,-nutri)",
  "(C,+nutri) - (DW,+nutri)"
)

# pH: the Day-0 contrast (section 8b) takes the place of the C vs DW -nutri
# contrast. It is put in that contrast's column so the matrices keep their
# 14 x 4 shape; `cell_note` and the dagger in Fig. S2 mark the cell.
pH_d0_synth_swapped <- pH_d0_synth %>%
  dplyr::transmute(Parameter, contrast = pH_invalid_contrast, p_Fisher, r_meta)

syn_all <- syn_all %>%
  dplyr::filter(contrast %in% wanted_contrasts) %>%
  dplyr::filter(!(Parameter == "pH" & contrast == pH_invalid_contrast)) %>%
  dplyr::bind_rows(pH_d0_synth_swapped) %>%
  dplyr::mutate(
    # one BH adjustment across 14 endpoints x 4 contrasts, still 56 tests
    q_BH = p.adjust(p_Fisher, method = "BH"),
    cell_note = dplyr::if_else(Parameter == "pH" & contrast == pH_invalid_contrast,
                               "pH: (DW,-nutri) day 0 - end of incubation, not C vs DW",
                               NA_character_))

syn_p_mat <- syn_all %>% select(Parameter, contrast, p_Fisher) %>%
  tidyr::pivot_wider(names_from = contrast, values_from = p_Fisher)  # kept for transparency
syn_q_mat <- syn_all %>% select(Parameter, contrast, q_BH) %>%
  tidyr::pivot_wider(names_from = contrast, values_from = q_BH)
syn_r_mat <- syn_all %>% select(Parameter, contrast, r_meta) %>%
  tidyr::pivot_wider(names_from = contrast, values_from = r_meta)

# Fig. 2 annotations -- treatment comparisons
# (pH's "(C,-nutri) - (DW,-nutri)" row is the Day-0 referenced one, flagged)
fig2_stats <- syn_all %>%
  dplyr::filter(Parameter %in% c("ORP", "DO", "SO4", "pH", "Mn", "NO2", "NO3"),
                contrast %in% wanted_contrasts[3:4]) %>%
  dplyr::mutate(signif = sig_code(q_BH),
                flag   = dplyr::if_else(!is.na(cell_note), "†", ""))

# Fig. 3 annotations -- iron speciation, within incubation duration
long_all <- dplyr::bind_rows(long_all_list)

fig3_endpoints <- c("Fe", "Fe2", "Fe3")

# 60-day: always single-batch (k = 1); batch #2 is the only 60-day run
fig3_60d <- long_all %>%
  dplyr::filter(Parameter %in% fig3_endpoints, Batch == "#2", contrast %in% wanted_contrasts) %>%
  dplyr::rename(p_val = p.value) %>%
  dplyr::mutate(Duration = "60d")

# 35-day, total Fe: batches #1 and #3 both ran 35d -> use m_35's pooled S*N contrast
fig3_35d_Fe <- dplyr::bind_rows(sn35_list) %>%
  dplyr::filter(Parameter == "Fe", contrast %in% wanted_contrasts) %>%
  dplyr::rename(p_val = p.value) %>%
  dplyr::mutate(Duration = "35d")

# 35-day, Fe2/Fe3: only batch #3 ran 35d for these -> single-batch
fig3_35d_Fe23 <- long_all %>%
  dplyr::filter(Parameter %in% c("Fe2", "Fe3"), Batch == "#3", contrast %in% wanted_contrasts) %>%
  dplyr::rename(p_val = p.value) %>%
  dplyr::mutate(Duration = "35d")

fig3_all <- dplyr::bind_rows(fig3_60d, fig3_35d_Fe, fig3_35d_Fe23) %>%
  dplyr::mutate(
    # one flat correction, Fig. 3's own family (3 endpoints x 4 contrasts x 2 durations = 24 tests)
    q_BH_fig3 = p.adjust(p_val, method = "BH"),
    signif    = sig_code(q_BH_fig3),
    contrast  = factor(contrast, levels = wanted_contrasts, ordered = TRUE),
    Duration  = factor(Duration, levels = c("35d", "60d"), ordered = TRUE)
  ) %>%
  dplyr::arrange(Parameter, Duration, contrast)

write_xlsx(
  list(Fig2_p_vals = fig2_stats, Fig3_p_vals = fig3_all),
  file.path(base_path, "Figure2_3_p_values.xlsx")
)

# ------------------------------------------------------------------------ #
# 11. One executive workbook
# ------------------------------------------------------------------------ #
summary_book <- list(
  ANOVA_35d_batches_p    = an35_p_mat,
  ANOVA_35d_batches_eta2 = an35_eta_mat,
  Mixed_LMM_ICC_R2       = icc_all,
  ANOVA_All_batches_p    = anAll_p_mat,
  ANOVA_All_batches_eta2 = anAll_eta_mat,
  SNB_contrasts_p        = snB_p_mat,
  SNB_contrasts_r        = snB_r_mat,
  SNT_contrasts_p        = snT_p_mat,
  SNT_contrasts_r        = snT_r_mat,
  SNB_BatchPairs_p       = snB_pairs_p,
  SNB_BatchPairs_r       = snB_pairs_r,
  SNT_TimePairs_p        = snT_pairs_p,
  SNT_TimePairs_r        = snT_pairs_r,
  BatchPairs_qBH         = BatchPairs_qBH,
  TimePairs_qBH          = TimePairs_qBH,
  SYN_pFisher_raw        = syn_p_mat,
  SYN_qBH_global         = syn_q_mat,
  SYN_r_meta             = syn_r_mat
)

write_xlsx(summary_book, file.path(base_path, "ART_ANOVA_Executive_Summaries_FDR.xlsx"))

cat("\nExecutive summary workbook written to:\n",
    file.path(base_path, "ART_ANOVA_Executive_Summaries_FDR.xlsx"), "\n")

# ------------------------------------------------------------------------ #
# 12. Heatmap helpers
# ------------------------------------------------------------------------ #

# p-value formatter: exact value, "< .001" threshold, or base-10 scientific
# with Unicode superscripts
.superscript_map <- c("0"="⁰","1"="¹","2"="²","3"="³","4"="⁴",
                      "5"="⁵","6"="⁶","7"="⁷","8"="⁸","9"="⁹",
                      "-"="⁻")
to_superscript <- function(n) {
  s <- as.character(n)
  paste0(unlist(.superscript_map[strsplit(s, "")[[1]]]), collapse = "")
}

format_q10 <- function(x, sig_digits = 2, sci_threshold = 1e-3,
                       style = c("sci10","threshold")) {
  style <- match.arg(style)
  vapply(x, function(p) {
    if (is.na(p)) return("")
    if (p == 0)  return("0")
    if (p >= sci_threshold) {
      # exact p with 2-3 sig figs; strip trailing zeros/dot
      lab <- format(signif(p, sig_digits), trim = TRUE, scientific = FALSE)
      lab <- sub("0+$","", sub("\\.$","", lab))
      return(lab)
    } else {
      if (style == "threshold") {
        return("< .001")
      } else {
        e <- floor(log10(p))
        m <- p / (10^e)
        mlab <- format(signif(m, 1), trim = TRUE, scientific = FALSE)
        return(paste0(mlab, "×10", to_superscript(e)))
      }
    }
  }, FUN.VALUE = character(1))
}

# colour-scale legend number formatter
label_drop0 <- function(digits = NULL) {
  force(digits)
  function(x) {
    if (!is.null(digits)) x <- round(x, digits)
    prettyNum(x, big.mark = "", scientific = FALSE, drop0trailing = TRUE)
  }
}

proper_labels_unicode <- c(
  ORP          = "ORP",
  DO           = "DO",
  Fe           = "Total Fe",
  Fe2          = "Fe²⁺",
  Fe3          = "Fe³⁺",
  Mn           = "Mn",
  SO4          = "SO₄²⁻",
  NO3          = "NO₃⁻",
  NO2          = "NO₂⁻",
  TH           = "Total hardness",
  Cl           = "Cl⁻",
  Turbidity    = "Turbidity",
  Conductivity = "Conductivity",
  pH           = "pH"
)

# endpoint order for the discrete y axis
full_endpoint_order <- c(
  "ORP", "DO", "Turbidity", "TH", "Conductivity", "pH",
  "Fe", "Fe2", "Fe3", "Mn", "SO4", "Cl", "NO2", "NO3"
)

# Turn "(C,-nutri,#1) - (C,-nutri,#3)" into "C -nutri": the right-hand side of
# these contrasts repeats the left except for the factor being varied, so
# printing it again only adds noise.
abbreviate_contrast <- function(x) {
  x %>%
    stringr::str_replace("\\) - \\(.*", "") %>%
    stringr::str_replace_all("[()]", "") %>%
    stringr::str_replace_all(",",   " ") %>%
    stringr::str_replace_all("35d|60d", "") %>%
    stringr::str_replace_all("#1|#3", "") %>%
    stringr::str_replace_all("\\s-\\s", "  −  ")
}

# ---- single-block heatmap: signed r in the tile, bold where q < alpha ---- #
make_contrast_heatmap_rq <- function(q_wide, r_wide,
                                     id_col = "contrast",
                                     title  = "BatchPairs (r, qBH across endpoints)",
                                     alpha  = 0.05,
                                     abbrev = TRUE) {
  stopifnot(id_col %in% names(q_wide), id_col %in% names(r_wide))

  q_long <- q_wide %>% tidyr::pivot_longer(-tidyr::all_of(id_col), names_to = "Endpoint", values_to = "q")
  r_long <- r_wide %>% tidyr::pivot_longer(-tidyr::all_of(id_col), names_to = "Endpoint", values_to = "r")

  dat <- dplyr::left_join(q_long, r_long, by = c(id_col, "Endpoint")) %>%
    dplyr::mutate(sig = q < alpha,
                  xlab = if (abbrev) abbreviate_contrast(.data[[id_col]]) else .data[[id_col]])

  # lock x order as encountered
  x_levels <- dat %>% dplyr::distinct(xlab) %>% dplyr::pull(xlab)

  p <- ggplot() +
    # 1) grey tiles for non-significant |r|
    geom_tile(data = dat %>% dplyr::filter(!sig),
              aes(x = factor(xlab, levels = x_levels), y = Endpoint, fill = abs(r)),
              color = "white", linewidth = 0.25) +
    scale_fill_gradient(name = "|r| (ns)", low = "#F0F0F0", high = "#9E9E9E") +
    ggnewscale::new_scale_fill() +
    # 2) diverging tiles for significant, signed r
    geom_tile(data = dat %>% dplyr::filter(sig),
              aes(x = factor(xlab, levels = x_levels), y = Endpoint, fill = r),
              color = "white", linewidth = 0.25) +
    geom_text(data = dat,
              aes(x = factor(xlab, levels = x_levels), y = Endpoint,
                  label = ifelse(is.na(r), "", sprintf("%.2f", r)),
                  fontface = ifelse(sig, "bold", "plain")),
              size = 3) +
    scale_fill_gradient2(name = "r (sig)", low = "#2C7BB6", mid = "white", high = "#D7191C",
                         midpoint = 0, limits = c(-1, 1), oob = scales::squish) +
    labs(x = NULL, y = NULL, title = title,
         subtitle = sprintf("alpha = %.02f; bold labels where q < alpha", alpha)) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          panel.grid  = element_blank(),
          plot.title  = element_text(face = "bold"))

  list(plot = p, data = dat)
}

# ---- dual heatmap: q block on the left, r block on the right ---- #
make_contrast_heatmap_dual <- function(q_wide, r_wide,
                                       id_col = "contrast",
                                       title  = "Contrasts across endpoints",
                                       alpha  = 0.05,
                                       abbrev = TRUE) {
  stopifnot(id_col %in% names(q_wide), id_col %in% names(r_wide))

  q_long <- q_wide %>% pivot_longer(-all_of(id_col), names_to = "Endpoint", values_to = "q")
  r_long <- r_wide %>% pivot_longer(-all_of(id_col), names_to = "Endpoint", values_to = "r")

  dat <- q_long %>%
    left_join(r_long, by = c(id_col, "Endpoint")) %>%
    mutate(sig = q < alpha,
           xcore = if (abbrev) abbreviate_contrast(.data[[id_col]]) else .data[[id_col]])

  p_block <- dat %>% mutate(xlab = paste0("q: ", xcore),
                            q_lab = format_q10(q, sig_digits = 2, sci_threshold = 1e-3,
                                               style = "threshold"))
  r_block <- dat %>% mutate(xlab = paste0("r: ", xcore))

  # Fe2/Fe3 only appear in the Time-pair family; the Batch-pair family has 12
  # endpoints because Fe(II)/Fe(III) were not measured in experiment #1.
  endpoint_order <- full_endpoint_order[full_endpoint_order %in% dat$Endpoint]

  x_levels <- c(unique(p_block$xlab), unique(r_block$xlab))

  p <- ggplot() +
    # ----- LEFT block: q values (dark blue -> white) -----
    geom_tile(data = p_block,
              aes(x = factor(xlab, levels = x_levels), y = Endpoint, fill = q),
              color = "grey40", linewidth = 0.25) +
    geom_text(data = p_block,
              aes(x = factor(xlab, levels = x_levels), y = Endpoint, label = q_lab,
                  colour = ifelse(q < alpha, "white", "black")),
              size = 3, show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_gradient(name = "q (adj. p)", low = "#061E47", high = "white",
                        limits = c(0, 1), oob = scales::squish,
                        labels = label_drop0(digits = 2)) +
    ggnewscale::new_scale_fill() +
    # ----- RIGHT block: non-significant r as grey |r| -----
    geom_tile(data = r_block %>% filter(!sig),
              aes(x = factor(xlab, levels = x_levels), y = Endpoint, fill = abs(r)),
              color = "white", linewidth = 0.25) +
    scale_fill_gradient(name = "|r| (ns)", low = "#F0F0F0", high = "#9E9E9E",
                        labels = label_drop0(digits = 2)) +
    ggnewscale::new_scale_fill() +
    # ----- RIGHT block: significant r as diverging -----
    geom_tile(data = r_block %>% filter(sig),
              aes(x = factor(xlab, levels = x_levels), y = Endpoint, fill = r),
              color = "white", linewidth = 0.1) +
    geom_text(data = r_block,
              aes(x = factor(xlab, levels = x_levels), y = Endpoint,
                  label = ifelse(is.na(r), "", sprintf("%.2f", r)),
                  fontface = ifelse(sig, "bold", "plain")),
              size = 3) +
    scale_fill_gradient2(name = "r (sig)", low = "#2C7BB6", mid = "white", high = "#D7191C",
                         midpoint = 0, limits = c(-1, 1), oob = scales::squish,
                         labels = label_drop0(digits = 2)) +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 11) +
    scale_y_discrete(limits = endpoint_order,
                     labels = function(x) {
                       # map to pretty labels; fall back to the raw id
                       out <- proper_labels_unicode[x]
                       out[is.na(out)] <- x[is.na(out)]
                       out
                     }) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          panel.grid  = element_blank(),
          plot.title  = element_text(face = "bold"))

  list(plot = p, data = dplyr::bind_rows(p_block, r_block))
}

# ---- dual ANOVA heatmap: p on the left, partial eta^2 on the right ---- #
make_dual_anova_heatmap <- function(p_wide, eta_wide,
                                    title = "ANOVA (p + η²)", alpha = 0.05,
                                    term_order = NULL,
                                    p_label_style = c("threshold","sci10")) {
  p_label_style <- match.arg(p_label_style)
  stopifnot(identical(names(p_wide), names(eta_wide)))

  p_long   <- p_wide   %>% pivot_longer(-Parameter, names_to = "Term", values_to = "p")
  eta_long <- eta_wide %>% pivot_longer(-Parameter, names_to = "Term", values_to = "eta2")

  dat <- left_join(p_long, eta_long, by = c("Parameter","Term")) %>%
    mutate(sig = !is.na(p) & p < alpha)   # NA-safe: NA p would give fontface NA

  if (!is.null(term_order)) dat$Term <- factor(dat$Term, levels = term_order)

  p_block <- dat %>% mutate(xlab = paste0("p: ", as.character(Term)))
  e_block <- dat %>% mutate(xlab = paste0("η²: ", as.character(Term)))

  # the default order is not the wanted one, so force it: build x_levels from
  # term_order so the left and right blocks line up exactly
  if (!is.null(term_order)) {
    p_block$Term <- factor(p_block$Term, levels = term_order)
    e_block$Term <- factor(e_block$Term, levels = term_order)
    x_levels <- c(paste0("p: ", term_order), paste0("η²: ", term_order))
  } else {
    x_levels <- c(unique(p_block$xlab), unique(e_block$xlab))
  }

  p_block <- p_block %>%
    mutate(p_lab = format_q10(p, sig_digits = 2, sci_threshold = 1e-3, style = p_label_style))

  # derive the y (endpoint) order from what is actually present
  endpoint_order <- full_endpoint_order[full_endpoint_order %in% unique(p_wide$Parameter)]

  p <- ggplot() +
    # LEFT: p
    geom_tile(data = p_block,
              aes(x = factor(xlab, levels = x_levels), y = Parameter, fill = p),
              color = "grey40", linewidth = 0.25) +
    geom_text(data = p_block,
              aes(x = factor(xlab, levels = x_levels), y = Parameter, label = p_lab,
                  colour = ifelse(p < alpha, "white", "black")),
              size = 3, show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_gradient(name = "p", low = "#061E47", high = "white",
                        limits = c(0,1), oob = scales::squish,
                        labels = label_drop0(digits = 2),
                        guide = guide_colourbar(order = 1)) +
    new_scale_fill() +
    # RIGHT: non-significant eta^2 (grey)
    geom_tile(data = e_block %>% filter(!sig),
              aes(x = factor(xlab, levels = x_levels), y = Parameter, fill = eta2),
              color = "white", linewidth = 0.25) +
    scale_fill_gradient(name = expression(partial~eta^~{2}~"(ns)"),
                        low = "#F0F0F0", high = "#A0A0A0", limits = c(0, NA),
                        labels = label_drop0(digits = 2),
                        guide = guide_colourbar(order = 3)) +
    new_scale_fill() +
    # RIGHT: significant eta^2
    geom_tile(data = e_block %>% filter(sig),
              aes(x = factor(xlab, levels = x_levels), y = Parameter, fill = eta2),
              color = "white", linewidth = 0.25) +
    geom_text(data = e_block,
              aes(x = factor(xlab, levels = x_levels), y = Parameter,
                  label = sprintf("%.2f", eta2),
                  fontface = ifelse(sig, "bold","plain")),
              size = 3) +
    scale_fill_gradient(name = expression(partial~eta^~{2}~"(sig)"),
                        low = "white", high = "#DA7621",
                        labels = label_drop0(digits = 2),
                        guide = guide_colourbar(order = 2)) +
    labs(x = NULL, y = NULL, title = title,
         subtitle = sprintf("Left: p (white text where p < %.02f) | Right: partial η² (bold where p < %.02f)",
                            alpha, alpha)) +
    theme_minimal(base_size = 11) +
    scale_y_discrete(limits = endpoint_order,
                     labels = function(x) {
                       out <- proper_labels_unicode[x]
                       out[is.na(out)] <- x[is.na(out)]
                       out
                     }) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          panel.grid  = element_blank(),
          plot.title  = element_text(face = "bold"))

  list(plot = p, data = dplyr::bind_rows(p_block, e_block))
}

# ---- Fisher-combined synthesis heatmap: rows = endpoints, cols = contrasts ---- #
# Selected cells (here: pH x "(C,-nutri) - (DW,-nutri)", which holds the Day-0
# referenced contrast of section 8b) get `flag_symbol` appended to the q value
# printed in the tile; with flag_r = TRUE also to the r value. Axis and column
# labels stay unchanged.
#   flag_cells  data.frame with columns RowID and ColID, using the RAW names as
#               in q_wide, e.g. data.frame(RowID = "pH", ColID = "(C,-nutri) - (DW,-nutri)")
make_Fisher_combined_contrast_heatmap_dual <- function(q_wide, r_wide,
                                                       row_id = "Parameter",
                                                       title  = "Contrasts across parameters — Fisher qBH",
                                                       alpha  = 0.05,
                                                       abbrev = TRUE,
                                                       row_order = NULL,
                                                       row_labels = NULL,
                                                       flag_cells  = NULL,
                                                       flag_symbol = "†",
                                                       flag_r      = FALSE) {
  stopifnot(row_id %in% names(q_wide), row_id %in% names(r_wide))

  if (!is.null(flag_cells)) {
    flag_cells <- as.data.frame(flag_cells)
    stopifnot(all(c("RowID", "ColID") %in% names(flag_cells)))
    # fail loudly on typos in the row / contrast name
    miss <- flag_cells[!(flag_cells$ColID %in% names(q_wide)) |
                         !(flag_cells$RowID %in% q_wide[[row_id]]), , drop = FALSE]
    if (nrow(miss) > 0) {
      stop("flag_cells refers to cells that are not in q_wide: ",
           paste(miss$RowID, "/", miss$ColID, collapse = "; "))
    }
  }
  flag_key <- if (is.null(flag_cells)) character(0) else paste(flag_cells$RowID, flag_cells$ColID, sep = "||")

  q_long <- q_wide %>% dplyr::rename(RowID = !!row_id) %>%
    tidyr::pivot_longer(-RowID, names_to = "ColID", values_to = "q")
  r_long <- r_wide %>% dplyr::rename(RowID = !!row_id) %>%
    tidyr::pivot_longer(-RowID, names_to = "ColID", values_to = "r")

  dat <- dplyr::left_join(q_long, r_long, by = c("RowID","ColID")) %>%
    dplyr::mutate(
      # enforce numeric types for plotting
      q = suppressWarnings(as.numeric(q)),
      r = suppressWarnings(as.numeric(r)),
      flag = paste(RowID, ColID, sep = "||") %in% flag_key,
      sig = !is.na(q) & q < alpha)

  # Both sides of a treatment contrast differ here, so unlike
  # abbreviate_contrast() the right-hand side is kept -- only the parentheses
  # and commas are cleaned up.
  dat <- dat %>%
    mutate(ColID_clean = if (abbrev) {
      ColID %>%
        stringr::str_replace_all("[()]", "") %>%
        stringr::str_replace_all(",",   " ") %>%
        stringr::str_replace_all("\\s-\\s", "  −  ")
    } else ColID)

  # the marker goes only into the flagged tile's label
  p_block <- dat %>% mutate(xlab = paste0("q: ", ColID_clean),
                            q_lab = format_q10(q, sig_digits = 2, sci_threshold = 1e-3,
                                               style = "threshold"),
                            q_lab = ifelse(flag, paste0(q_lab, flag_symbol), q_lab))
  r_block <- dat %>% mutate(xlab = paste0("r: ", ColID_clean),
                            r_lab = ifelse(is.na(r), "", sprintf("%.2f", r)),
                            r_lab = ifelse(flag & flag_r, paste0(r_lab, flag_symbol), r_lab))
  x_levels <- c(unique(p_block$xlab), unique(r_block$xlab))

  if (!is.null(row_order)) {
    p_block <- p_block %>% dplyr::mutate(RowID = factor(RowID, levels = row_order))
    r_block <- r_block %>% dplyr::mutate(RowID = factor(RowID, levels = row_order))
  }

  p <- ggplot() +
    # ----- LEFT block: Fisher q -----
    geom_tile(data = p_block,
              aes(x = factor(xlab, levels = x_levels), y = RowID, fill = q),
              color = "grey40", linewidth = 0.25) +
    geom_text(data = p_block,
              aes(x = factor(xlab, levels = x_levels), y = RowID, label = q_lab,
                  colour = ifelse(sig, "white", "black"),
                  fontface = ifelse(sig, "bold", "plain")),
              size = 3, show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_gradient(name = "q (Fisher adj. p)", low = "#27335A", high = "white",
                        limits = c(0, 1), oob = scales::squish,
                        na.value = "grey90",
                        labels = label_drop0(digits = 2),
                        guide = guide_colourbar(order = 1)) +
    ggnewscale::new_scale_fill() +
    # ----- RIGHT block: non-significant meta-r as grey |r| -----
    geom_tile(data = r_block %>% filter(!sig),
              aes(x = factor(xlab, levels = x_levels), y = RowID, fill = abs(r)),
              color = "white", linewidth = 0.25) +
    scale_fill_gradient(name = "|r| (ns)", low = "#F0F0F0", high = "#9E9E9E",
                        limits = c(0,1), breaks = seq(0, 1, by = 0.25),
                        na.value = "grey90",
                        labels = label_drop0(digits = 2),
                        guide = guide_colourbar(order = 3)) +
    ggnewscale::new_scale_fill() +
    # ----- RIGHT block: significant meta-r as diverging -----
    geom_tile(data = r_block %>% filter(sig),
              aes(x = factor(xlab, levels = x_levels), y = RowID, fill = r),
              color = "white", linewidth = 0.1) +
    geom_text(data = r_block,
              aes(x = factor(xlab, levels = x_levels), y = RowID,
                  label = r_lab,
                  fontface = ifelse(sig, "bold", "plain")),
              size = 3) +
    scale_fill_gradient2(name = "r (sig)", low = "#2C7BB6", mid = "white", high = "#D7191C",
                         midpoint = 0, limits = c(-1, 1), oob = scales::squish,
                         labels = label_drop0(digits = 2),
                         guide = guide_colourbar(order = 2)) +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 11) +
    scale_y_discrete(labels = function(x) {
      if (is.null(row_labels)) return(x)
      out <- unname(row_labels[x])
      out[is.na(out)] <- x[is.na(out)]
      out
    }) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          panel.grid  = element_blank(),
          plot.title  = element_text(face = "bold"))

  list(plot = p, data = dplyr::bind_rows(p_block, r_block))
}

# ------------------------------------------------------------------------ #
# 13. Heatmap figures
# ------------------------------------------------------------------------ #

# --- Fig. S3: batch reproducibility, experiment #1 vs #3 at fixed treatment
h_batch <- make_contrast_heatmap_rq(
  q_wide = BatchPairs_qBH,        # BH-adjusted p across endpoints, per contrast
  r_wide = snB_pairs_r,           # r matrix for the same contrasts
  id_col = "contrast",
  title  = "Experiment #1 vs #3 (r with qBH across endpoints)",
  alpha  = 0.05,
  abbrev = TRUE
)
ggsave(file.path(fig_dir, "FigS3_heatmap_Batch-pair_comparisons.png"),
       plot = h_batch$plot, width = 14, height = 7, dpi = 300)
write.csv(h_batch$data,
          file.path(fig_dir, "FigS3_heatmap_Batch-pair_comparisons_plotting_data.csv"),
          row.names = FALSE)

# --- Duration reproducibility, 35 d vs 60 d at fixed treatment.
# Not in the manuscript, kept because it is the companion to Fig. S3.
h_time <- make_contrast_heatmap_rq(
  q_wide = TimePairs_qBH,
  r_wide = snT_pairs_r,
  id_col = "contrast",
  title  = "35 d vs 60 d (r with qBH across endpoints)",
  alpha  = 0.05,
  abbrev = TRUE
)
ggsave(file.path(fig_dir, "Heatmap_TimePairs_qBH.png"),
       plot = h_time$plot, width = 14, height = 7, dpi = 300)
write.csv(h_time$data,
          file.path(fig_dir, "Heatmap_TimePairs_qBH_plotting_data.csv"),
          row.names = FALSE)

# --- The same two families as q + r side by side (not in the manuscript)
h_dual_batch <- make_contrast_heatmap_dual(
  q_wide = BatchPairs_qBH,
  r_wide = snB_pairs_r,
  id_col = "contrast",
  title  = ""
)
ggsave(file.path(fig_dir, "Heatmap_dual_batch-pair_contrasts_q_and_r.png"),
       plot = h_dual_batch$plot, width = 7, height = 7, dpi = 300)
write.csv(h_dual_batch$data,
          file.path(fig_dir, "Heatmap_dual_batch-pair_contrasts_q_and_r_plotting_data.csv"),
          row.names = FALSE)

h_dual_time <- make_contrast_heatmap_dual(
  q_wide = TimePairs_qBH,
  r_wide = snT_pairs_r,
  id_col = "contrast",
  title  = "Contrasts across parameters — 35 vs 60 day duration"
)
ggsave(file.path(fig_dir, "Heatmap_dual_time-pair_contrasts_q_and_r.png"),
       plot = h_dual_time$plot, width = 7, height = 7, dpi = 300)
write.csv(h_dual_time$data,
          file.path(fig_dir, "Heatmap_dual_time-pair_contrasts_q_and_r_plotting_data.csv"),
          row.names = FALSE)

# --- Omnibus ANOVA overviews (not in the manuscript)
term_order_35 <- c(
  "Sterilization", "Nutrient", "Batch",
  "Sterilization:Nutrient", "Sterilization:Batch", "Nutrient:Batch",
  "Sterilization:Nutrient:Batch"
)

h_an35 <- make_dual_anova_heatmap(
  p_wide   = an35_p_mat,
  eta_wide = an35_eta_mat,
  title    = "ANOVA (35 d; S×N×Batch): p + η²",
  alpha    = 0.05,
  term_order = term_order_35,
  p_label_style = "threshold"
)

h_anAll <- make_dual_anova_heatmap(
  p_wide   = anAll_p_mat,
  eta_wide = anAll_eta_mat,
  title    = "ANOVA (All batches; S×N×Time): p + η²",
  alpha    = 0.05,
  p_label_style = "threshold"
)

ggsave(file.path(fig_dir, "Heatmap_ANOVA_35d_SxNxBatch_p_and_eta2.png"),
       h_an35$plot,  width = 12, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "Heatmap_ANOVA_all-batches_SxNxTime_p_and_eta2.png"),
       h_anAll$plot, width = 12, height = 8, dpi = 300)
write.csv(h_an35$data,
          file.path(fig_dir, "Heatmap_ANOVA_35d_SxNxBatch_p_and_eta2_plotting_data.csv"),
          row.names = FALSE)
write.csv(h_anAll$data,
          file.path(fig_dir, "Heatmap_ANOVA_all-batches_SxNxTime_p_and_eta2_plotting_data.csv"),
          row.names = FALSE)

# --- Fig. S2: Fisher-combined treatment contrasts + meta-analytic r
row_labels_order <- c(
  "ORP", "DO", "Turbidity", "TH", "Conductivity", "pH",
  "Fe", "Fe2", "Fe3", "Mn", "SO4", "Cl", "NO2", "NO3"
)

h_syn <- make_Fisher_combined_contrast_heatmap_dual(
  q_wide = syn_q_mat,   # rows = Parameter, cols = contrasts (qBH)
  r_wide = syn_r_mat,   # same shape (r_meta)
  row_id = "Parameter",
  title  = "",
  abbrev = TRUE,
  row_order = row_labels_order,
  row_labels = proper_labels_unicode
  # pH: the Day-0 referenced contrast sits in this cell (section 8b)
  flag_cells = data.frame(RowID = "pH", ColID = pH_invalid_contrast),
  flag_r = TRUE
)
ggsave(file.path(fig_dir, "FigS2_heatmap_treatment_contrasts.png"),
       h_syn$plot, width = 8, height = 8, dpi = 300)
write.csv(h_syn$data,
          file.path(fig_dir, "FigS2_heatmap_treatment_contrasts_plotting_data.csv"),
          row.names = FALSE)

cat("\nDone. Tables in", base_path, "- figures in", fig_dir, "\n")
