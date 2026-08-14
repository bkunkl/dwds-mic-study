# ======================================================================== #
# Fig. S4 -- Sample completeness (coverage) vs. sequencing depth
#
# iNEXT "type 2" plot: sample coverage (SC) against the number of individuals
# (reads assigned to a taxon), for the 13 metagenomes. Coverage is identical
# for q = 0, 1 and 2 at a given sample size, so only the q = 0 object is used
# for the curves; the q0/q1/q2 objects are all loaded because the optional
# panels at the end of the script use them.
#
# Panel A: full range, iNEXT layer logic -- solid line = rarefaction,
#          dashed = extrapolation, filled circle = the observed sample size.
# Panel B: the same curves zoomed to 0-10,000 individuals, where the curves
#          rise; this is where the near-complete coverage becomes visible.
#
# Facets are the five sample types; colour is the experiment.
# `pilot1`/`pilot2`/`pilot3` in the Assemblage names correspond to E1/E2/E3.
# Only E3 has the pre-incubation Base-Water control, so that facet holds a
# single curve.
#
# INPUT (produced by the HPC chain, see hpc/ and scripts/local/04_merge_rarefaction_results.R)
#   results/inext_rarefaction/inext_q[0-2]_combined_columnwise_inext_object_dense_low_sizes_result.rds
#       class `iNEXT` objects, one per Hill number, merged from the
#       per-sample (column-wise) SLURM array results.
#
# OUTPUT (into results/figures/)
#   FigS4_type2_sample_coverage.png
#   FigS4_type2_sample_coverage_panelA_plotting_data.csv
#   FigS4_type2_sample_coverage_panelB_plotting_data.csv
#   plus, optionally (section 6), the standard ggiNEXT panels, the coverage
#   threshold table and the 95%-coverage reference-line variant
#
# Run from the repository root.
# ======================================================================== #

library(iNEXT)      # only needed by the optional panels in section 6
library(ggplot2)
library(pals)

# repo-relative paths -- run this script with the repo root as working directory
inext_dir <- "results/inext_rarefaction"

fig_dir <- "results/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------ #
# 0. Axis label helpers
# ------------------------------------------------------------------------ #

# Pretty labels with plotmath for big values; plain numbers otherwise.
# - threshold: boundary where we switch to mantissa × 10^exponent
# - digits: significant digits for the mantissa
label_exp_pretty <- function(threshold = 1e5, digits = 3) {
  force(threshold); force(digits)
  function(x) {
    # Build a character vector of the same length as x
    out <- character(length(x))
    for (i in seq_along(x)) {
      v <- x[i]
      if (is.na(v)) { out[i] <- "NA"; next }
      if (v == 0)  { out[i] <- "0";  next }

      if (abs(v) < threshold) {
        # Plain numeric; no thousands separator to keep parse() happy
        # Use signif() to trim trailing .0 and excess digits
        s <- as.character(signif(v, digits = digits))
        s <- sub("\\.?0+$", "", s)
        out[i] <- s
      } else {
        pow  <- floor(log10(abs(v)))
        mant <- round(v / (10^pow), digits = digits - 1)
        mant_str <- sub("\\.?0+$", "", format(mant, trim = TRUE))
        # Plotmath: "2 %*% 10^6"
        out[i] <- paste0(mant_str, " %*% 10^", pow)
      }
    }
    # Parse to plotmath expressions; length is guaranteed to match length(x)
    scales::label_parse()(out)
  }
}

# Label function for the coverage axis: 0, 0.25, 0.5, 0.75, 1
label_coverage <- function() {
  function(x) {
    labs <- ifelse(is.na(x), NA, sprintf("%.2f", x))
    sub("\\.?0+$", "", labs)  # trim trailing ".0" / ".00"
  }
}

# ------------------------------------------------------------------------ #
# 1. Load the combined iNEXT objects
# ------------------------------------------------------------------------ #

# Locate the three RDS files for q = 0, 1, 2
inext_objs_paths <- Sys.glob(file.path(
  inext_dir, "inext_q[0-2]_combined_columnwise_inext_object_dense_low_sizes_result.rds"))

# Extract q values from filenames, order paths by q, load, and name the list accordingly
q_vals <- as.integer(sub(".*inext_q([0-2]).*", "\\1", basename(inext_objs_paths)))
ord <- order(q_vals, na.last = NA)

inext_obj_list <- setNames(
  lapply(inext_objs_paths[ord], readRDS),
  paste0("q", q_vals[ord])
)

# ------------------------------------------------------------------------ #
# 2. Panel A -- sample completeness over the full range
# ------------------------------------------------------------------------ #

# --- assemble TYPE 2 data (SC vs size) across q, then keep q = 0 only ---
size2_all <- do.call(rbind, lapply(inext_obj_list, function(obj) {
  df <- obj$iNextEst$size_based
  df$Treatment <- sub(".*: (.*)", "\\1", df$Assemblage)
  df$Batch     <- sub("^(.*?):.*$", "\\1", df$Assemblage)
  # x = m (abundance) or t (incidence), as iNEXT type=2 does
  xvar <- if (names(obj$DataInfo)[2] == "n") "m" else "t"  # abundance vs incidence
  df$x <- df[[xvar]]
  df
}))

# --- keep only q = 0 because SC is the same across q at given size ---
size2_all <- subset(size2_all, Order.q == 0)

# --- iNEXT-like layer logic: rename "Observed" -> "Rarefaction" for lines; keep Observed separately for the point ---
size2_all$Method_raw <- size2_all$Method
size2_all$Method[size2_all$Method == "Observed"] <- "Rarefaction"
size2_all$lty <- factor(size2_all$Method, levels = c("Rarefaction","Extrapolation"))

# observed-point data only
size2_obs <- subset(size2_all, Method_raw == "Observed")

# factors for aesthetics
size2_all$Treatment <- factor(size2_all$Treatment)
size2_all$Batch     <- factor(size2_all$Batch)
size2_obs$Treatment <- factor(size2_obs$Treatment, levels = levels(size2_all$Treatment))
size2_obs$Batch     <- factor(size2_obs$Batch,     levels = levels(size2_all$Batch))

# Axis label, identical to iNEXT's logic
xlab_type2 <- if (any(names(inext_obj_list[[1]]$DataInfo)[2] == "n"))
  "Number of individuals" else "Number of sampling units"

# --- plot ---
p_sc_vs_size_inext <-
  ggplot(size2_all, aes(x = x, y = SC, colour = Batch)) +
  # solid = rarefaction, dashed = extrapolation
  geom_line(aes(linetype = lty), linewidth = 0.7, alpha = 0.95) +
  # single filled circle ONLY at the observed size; map shape to merge into the 'Method' legend
  geom_point(data = size2_obs, aes(shape = "Observed"), size = 1.8, stroke = 0) +
  facet_wrap(~ Treatment, nrow = 1, scales = "fixed") +
  # Batch colours
  scale_color_manual(values = pals::stepped(n = 24)[c(5, 13, 10, 22)]) +
  # Merge linetype & shape under a single legend called "Method"
  scale_linetype_manual(
    name   = "Method",
    values = c(Rarefaction = "solid", Extrapolation = "22"),
    breaks = c("Rarefaction","Extrapolation")
  ) +
  scale_shape_manual(
    name   = "Method",
    values = c(Observed = 16),
    breaks = "Observed"
  ) +
  guides(
    linetype = guide_legend(order = 1),
    shape    = guide_legend(order = 1),
    colour   = guide_legend(order = 2, title = "Experiment") # or earlier "Batch"
  ) +
  scale_x_continuous(labels = label_exp_pretty()) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = label_coverage()
  ) +
  labs(x = xlab_type2, y = "Sample Coverage") +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.box       = "horizontal",
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    legend.key.width = grid::unit(1.2, "lines"),
    legend.key.height= grid::unit(0.8, "lines"),
    strip.text       = element_text(face = "bold"),
    axis.title.x     = element_text(margin = margin(t = 7.5, r = 0, b = 0, l = 0))
  )

# ------------------------------------------------------------------------ #
# 3. Panel B -- the same curves zoomed to 0-10,000 individuals
# ------------------------------------------------------------------------ #

# We take the q0 object, extract the size_based table, and set x=m (abundance)
# or x=t (incidence) as in iNEXT type=2.
obj_q0 <- inext_obj_list$q0

df2 <- obj_q0$iNextEst$size_based
# Determine the x-variable (sample size): m for abundance, t for incidence
xvar <- if (names(obj_q0$DataInfo)[2] == "n") "m" else "t"

# Keep only q = 0 and build Treatment / Batch from Assemblage
df2$q <- df2$Order.q
df2 <- subset(df2, q == 0)
df2$x <- df2[[xvar]]
df2$Treatment <- sub(".*: (.*)", "\\1", df2$Assemblage)      # after colon
df2$Batch     <- sub("^(.*?):.*$", "\\1", df2$Assemblage)    # before colon

# Ensure factors are consistent
df2$Treatment <- factor(df2$Treatment)
df2$Batch     <- factor(df2$Batch)

# Restrict to the zoom range [0, 10000]
df2_zoom_e4 <- subset(df2, x >= 0 & x <= 10000)

# --- Plot: one row, five columns (facet_wrap), continuous line + filled points ---
p_q0_SC_vs_size_zoom_e4 <-
  ggplot(df2_zoom_e4, aes(x = x, y = SC, group = Assemblage, colour = Batch)) +
  geom_line(linewidth = 0.7, alpha = 0.90) +                   # continuous lines
  geom_point(shape = 16, size = 1.4, alpha = 0.95) +           # filled circles at actual grid points
  facet_wrap(~ Treatment, nrow = 1, scales = "fixed") +        # 1 row × (#Treatments) columns
  scale_color_manual(values = pals::stepped(n = 24)[c(5, 13, 10, 22)]) +
  scale_x_continuous(
    limits = c(0, 10000),
    breaks = scales::breaks_extended(n = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = label_coverage()
  ) +
  labs(
    x = if (xvar == "m") "Number of individuals" else "Number of sampling units",
    y = "Sample Coverage",
    colour = "Experiment", # earlier "Batch"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.box       = "horizontal",
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    legend.key.width = grid::unit(1.2, "lines"),
    legend.key.height= grid::unit(0.8, "lines"),
    strip.text       = element_text(face = "bold"),
    plot.margin = margin(t = 5, r = 10, b = 1, l = 5),
    axis.title.x = element_text(margin = margin(t = 7.5, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))
  ) +
  guides(colour = guide_legend(title.position = "left", title.hjust = 1))

# ------------------------------------------------------------------------ #
# 4. Fig. S4 -- the two panels stacked
# ------------------------------------------------------------------------ #

p_sc_vs_size_inext_no_legend <- p_sc_vs_size_inext +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        # Necessary to align the two rows of plots,
        # because the second one has extra right margin.
        plot.margin = margin(t = 5, r = 10, b = 1, l = 5))
p_q0_SC_vs_size_zoom_e4_no_strip <- p_q0_SC_vs_size_zoom_e4 +
  ylab("")+
  theme(strip.text = element_blank())

p_type2_composit <- cowplot::plot_grid(
  p_sc_vs_size_inext_no_legend, p_q0_SC_vs_size_zoom_e4_no_strip,
  nrow = 2,
  rel_heights = c(4.35, 5.65),
  labels = c("A", "B"),
  label_size = 12,
  label_fontface = "bold",
  label_x = 0.01, label_y = 0.99, hjust = 0, vjust = 1)

out_file <- file.path(fig_dir, "FigS4_type2_sample_coverage.png")
ggsave(out_file, p_type2_composit,
       width = 30, height = 12, units = "cm", dpi = 300, bg = "white")

# ------------------------------------------------------------------------ #
# 5. Export the plotted data
# ------------------------------------------------------------------------ #

write.csv(size2_all,
          file = file.path(fig_dir, "FigS4_type2_sample_coverage_panelA_plotting_data.csv"),
          row.names = FALSE)

write.csv(df2_zoom_e4,
          file = file.path(fig_dir, "FigS4_type2_sample_coverage_panelB_plotting_data.csv"),
          row.names = FALSE)

cat("\nSaved 'FigS4_type2_sample_coverage.png' into", fig_dir, "\n")

# ======================================================================== #
# 6. Optional extras (not in the manuscript)
#
#    Everything below is independent of Fig. S4 and can be skipped. It is
#    kept because it is what the coverage claim in the manuscript rests on:
#    the standard iNEXT panels, the sequencing depth at which each sample
#    reaches 99% coverage, and a variant of panel B with the 95% coverage
#    reference line drawn in.
# ======================================================================== #

extras_dir <- file.path(fig_dir, "FigS4_extras")
dir.create(extras_dir, recursive = TRUE, showWarnings = FALSE)

# --- 6a. The three standard ggiNEXT panels for each Hill number ---------- #

# Map ggiNEXT's type to human-readable suffixes
type_labels <- c(`1` = "size_based",          # type = 1: sample-size-based
                 `2` = "sample_completeness", # type = 2: coverage vs size
                 `3` = "coverage_based")      # type = 3: coverage-based

for (nm in names(inext_obj_list)) {
  obj <- inext_obj_list[[nm]]

  for (tp in 1:3) {
    p <- ggiNEXT(obj, type = tp)

    # File name: e.g. FigS4_extras/q0_1_size_based.png
    outfile <- file.path(extras_dir, sprintf("%s_%s_%s.png", nm, tp, type_labels[as.character(tp)]))

    # Save with consistent sizing
    ggsave(filename = outfile, plot = p,
           width = 30, height = 25, units = "cm", dpi = 300, bg = "white")
    message("Saved: ", outfile)
  }
}

# --- 6b. Depth at which each sample reaches 99% coverage ----------------- #
# For interpretation: the smallest sample size whose coverage rounds to 0.99.

df2_filtered <- df2[ round(df2$SC,2) == 0.99, ]
df2_min_filtered <- df2_filtered[ df2_filtered$m == ave(df2_filtered$m, df2_filtered$Assemblage, FUN = min), ]
df2_min_filtered[,c(1,2,5,8)]

write.csv(df2_min_filtered[,c(1,2,5,8)],
          file = file.path(extras_dir, "sample_size_at_99pct_coverage.csv"),
          row.names = FALSE)

# --- 6c. Panel B zoomed to 0-1000 with the 95% coverage reference line --- #

df2_zoom <- subset(df2, x >= 0 & x <= 1000)

p_q0_SC_vs_size_zoom_y95 <-
  ggplot(df2_zoom, aes(x = x, y = SC, group = Assemblage, colour = Batch)) +
  geom_hline(yintercept = 0.95,
             linetype = "dashed",
             linewidth = 0.5,
             color = "grey35")+
  geom_line(linewidth = 0.7, alpha = 0.90) +                   # continuous lines
  geom_point(shape = 16, size = 1.4, alpha = 0.95) +           # filled circles at actual grid points
  facet_wrap(~ Treatment, nrow = 1, scales = "fixed") +        # 1 row × (#Treatments) columns
  scale_color_manual(values = pals::stepped(n = 24)[c(5, 13, 10, 22)]) +
  scale_x_continuous(
    limits = c(0, 1000),
    breaks = scales::breaks_extended(n = 6)                    # ~6 pleasant ticks across 0–1000
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(seq(0, 0.75, by = 0.25), 0.95, 1),
    minor_breaks = seq(0.125, 0.875, by = 0.25),
    labels = label_coverage()
  ) +
  labs(
    x = if (xvar == "m") "Number of individuals" else "Number of sampling units",
    y = "Sample Coverage",
    colour = "Experiment",
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.box       = "horizontal",
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    legend.key.width = grid::unit(1.2, "lines"),
    legend.key.height= grid::unit(0.8, "lines"),
    strip.text       = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))
  )

out_file <- file.path(extras_dir, "type2_q0_SC_vs_size_zoom_Treatment_1x5_y95.png")
ggsave(out_file, p_q0_SC_vs_size_zoom_y95,
       width = 30, height = 8.5, units = "cm", dpi = 300, bg = "white")
