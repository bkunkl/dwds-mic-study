# ======================================================================== #
# Fig. 7 -- Sulfur bacteria (SOB + SRB) heatmap at GENUS level
#
# Individual replicates (i.e. built from `sulf_otu_mat`, NOT merged by Label),
# with a bottom bubble-plot annotation showing each replicate's total
# sulfur-bacteria relative abundance within its whole metagenomic community.
#
# Column layout ("Sample Type" color band on top, replicated across the
# 3 experiments):
#   #1: Biofilm-no-nutri, #1: Biofilm-w/-nutri, #1: Bulk-no-nutri, #1: Bulk-w/-nutri,
#   #2: Biofilm-no-nutri, #2: Biofilm-w/-nutri, #2: Bulk-no-nutri, #2: Bulk-w/-nutri,
#   #3: Biofilm-no-nutri, #3: Biofilm-w/-nutri, #3: Bulk-no-nutri, #3: Bulk-w/-nutri,
#   Base-Water
# `pilot1`/`pilot2`/`pilot3` in the metadata correspond to experiments E1/E2/E3.
# Only E3 has the pre-incubation Base-Water control, so it is a single column.
#
# Bubble diameters are area-scaled (diameter ~ sqrt(value)).
#
# INPUT
#   results/kraken_bracken/brackenSpecies_conf01.rds
#       phyloseq object built by scripts/local/01_build_phyloseq_object.sh
#       (Kraken2/Bracken, confidence 0.1, species level)
#   data/metadata/SOB_genus_list_comprehensive.csv
#   data/metadata/SRB_genus_list_comprehensive.csv
#       curated genus lists (one column, `Genus`) defining which genera count
#       as sulfur oxidizers (SOB) and sulfate reducers (SRB).
#       These two lists are hand-curated inputs, not pipeline output -- they
#       define the row set of this figure and must accompany the repository.
#
# OUTPUT (into results/figures/)
#   Fig7_complexHeatmap_sulfur_bacteria.png
#   Fig7_complexHeatmap_sulfur_bacteria_plotting_data.csv
#
# Run from the repository root. Needs the `gridtext` package installed for
# ComplexHeatmap's gt_render() (italic genus row labels).
# ======================================================================== #

library(phyloseq)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(pals)

# repo-relative paths -- run this script with the repo root as working directory
ps_file  <- "results/kraken_bracken/brackenSpecies_conf01.rds"
sob_file <- "data/metadata/SOB_genus_list_comprehensive.csv"
srb_file <- "data/metadata/SRB_genus_list_comprehensive.csv"

fig_dir <- "results/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------ #
# 0. Shared lookup tables / colors
# ------------------------------------------------------------------------ #

SOB_genera <- read.csv(file = sob_file, stringsAsFactors = FALSE)
SOB_genera <- SOB_genera[SOB_genera$Genus != "", ]

SRB_genera <- read.csv(file = srb_file, stringsAsFactors = FALSE)

functional_colors <- c("Sulfur Oxidizer" = "#8FB3C9", "Sulfate Reducer" = "#ECD788")

stepped_cols <- pals::stepped(n = 24)[c(5, 7, 13, 15, 22)]
names(stepped_cols) <- c("Biofilm-no-nutri", "Biofilm-w/-nutri",
                          "Bulk-no-nutri", "Bulk-w/-nutri", "Base-Water")

# ------------------------------------------------------------------------ #
# 1. Load raw data, keep individual replicates (no merge_samples by Label)
# ------------------------------------------------------------------------ #

pilot123_ps <- readRDS(ps_file)
pilot123_ps <- subset_taxa(pilot123_ps, Domain == "Bacteria" | Domain == "Archaea")
tax_table(pilot123_ps)[, "Genus"] <- as.character(tax_table(pilot123_ps)[, "Genus"])

# Relative abundance BEFORE subsetting to SOB/SRB
pilot123_ps_rel <- transform_sample_counts(pilot123_ps, function(x) x / sum(x))

# ------------------------------------------------------------------------ #
# 2. SOB / SRB subsetting (individual samples)
# ------------------------------------------------------------------------ #

pilot123_ps_sob_rel <- subset_taxa(
  pilot123_ps_rel,
  !is.na(Genus) & Genus %in% SOB_genera$Genus
)

# Methylobacterium/Methylorubrum and Rhodobacter are only partly sulfur
# oxidizing genera -- keep the species that are, drop the rest.
methylo_keep <- c("Methylobacterium extorquens", "Methylorubrum extorquens",
                   "Methylobacterium fujisawaense", "Methylorubrum sp. GM97")

pilot123_ps_sob_rel <- subset_taxa(
  pilot123_ps_sob_rel,
  !(Genus %in% c("Methylobacterium", "Methylorubrum")) | Species %in% methylo_keep
)

pilot123_ps_sob_rel <- subset_taxa(
  pilot123_ps_sob_rel,
  !(Genus == "Rhodobacter" & Species != "Rhodobacter capsulatus")
)

pilot123_ps_srb_rel <- subset_taxa(
  pilot123_ps_rel,
  !is.na(Genus) & Genus %in% SRB_genera$Genus
)

# Merge SOB + SRB taxa (still one row per sample = one replicate)
pilot123_ps_sulfur <- merge_phyloseq(pilot123_ps_sob_rel, pilot123_ps_srb_rel)

# ------------------------------------------------------------------------ #
# 3. Per-replicate total sulfur-bacteria abundance (used for bottom bubbles)
#    -- this is the fraction of the WHOLE community, computed BEFORE
#       renormalizing SOB+SRB to sum to 100%.
# ------------------------------------------------------------------------ #

sulfur_totals <- sample_sums(pilot123_ps_sulfur)

sulfur_metadata <- as.data.frame(sample_data(pilot123_ps_sulfur))
stopifnot(identical(rownames(sulfur_metadata), names(sulfur_totals)))

# ------------------------------------------------------------------------ #
# 4. Derive SOB+SRB community fractions (sum to 100% per replicate)
# ------------------------------------------------------------------------ #

pilot123_ps_sulfur_derived <- transform_sample_counts(
  pilot123_ps_sulfur,
  function(x) x / sum(x)
)

# ------------------------------------------------------------------------ #
# 5. Filter low-abundance / rare taxa at SPECIES level (individual samples)
# ------------------------------------------------------------------------ #

otu_table_derived <- as.data.frame(otu_table(pilot123_ps_sulfur_derived))
taxa_abundance <- rowSums(otu_table_derived)
taxa_prevalence <- rowSums(otu_table_derived > 0)

abundance_threshold <- 0.05  # 5% of SOB+SRB community
prevalence_threshold <- 5    # present in at least 5/13 samples

taxa_to_keep <- (taxa_abundance >= abundance_threshold) &
  (taxa_prevalence >= prevalence_threshold)

pilot123_ps_sulfur_filt <- prune_taxa(taxa_to_keep, pilot123_ps_sulfur_derived)

cat("\nIndividual-sample sulfur bacteria dataset (species level):",
    ntaxa(pilot123_ps_sulfur_filt), "taxa across",
    nsamples(pilot123_ps_sulfur_filt), "samples\n")

# ------------------------------------------------------------------------ #
# 6. Aggregate to GENUS level (individual replicates, not merged by Label)
# ------------------------------------------------------------------------ #

pilot123_ps_sulfur_indiv_genus <- tax_glom(pilot123_ps_sulfur_filt, taxrank = "Genus")

sulf_otu_mat <- as.matrix(otu_table(pilot123_ps_sulfur_indiv_genus))
sulf_tax_mat <- as.data.frame(tax_table(pilot123_ps_sulfur_indiv_genus))
sulf_sample_df <- as.data.frame(sample_data(pilot123_ps_sulfur_indiv_genus))

# ------------------------------------------------------------------------ #
# 7. Determine SOB vs. SRB per genus, sort rows: SOB (A-Z), then SRB (A-Z)
# ------------------------------------------------------------------------ #

sulf_taxa_names <- sulf_tax_mat$Genus
is_SOB <- sulf_taxa_names %in% SOB_genera$Genus
is_SRB <- sulf_taxa_names %in% SRB_genera$Genus

sulf_tax_mat$Functional_Group <- ifelse(is_SOB, "Sulfur Oxidizer", "Sulfate Reducer")

sob_taxa <- rownames(sulf_tax_mat)[is_SOB]
srb_taxa <- rownames(sulf_tax_mat)[is_SRB]

sob_order <- sob_taxa[order(sulf_tax_mat[sob_taxa, "Genus"])]
srb_order <- srb_taxa[order(sulf_tax_mat[srb_taxa, "Genus"])]

sulf_taxa_order <- c(sob_order, srb_order)

sulf_otu_mat_ordered <- sulf_otu_mat[sulf_taxa_order, ]

# ------------------------------------------------------------------------ #
# 8. Reorder columns: condition-major, batch-minor, Base-Water last
# ------------------------------------------------------------------------ #

condition_order <- c("Biofilm-no-nutri", "Biofilm-w/-nutri",
                      "Bulk-no-nutri", "Bulk-w/-nutri")
batch_order <- c("pilot1", "pilot2", "pilot3")

is_water <- sulf_sample_df$Label == "Base-Water"

non_water_df <- sulf_sample_df[!is_water, ]
non_water_df <- non_water_df[order(
  factor(non_water_df$Label, levels = condition_order),
  factor(non_water_df$batch, levels = batch_order)
), ]
# If we want batch #1 with all samples then batch #2 and so on:
# non_water_df <- non_water_df[order(
#   factor(non_water_df$batch, levels = batch_order),
#   factor(non_water_df$Label, levels = condition_order)
# ), ]

water_df <- sulf_sample_df[is_water, ]

sample_order <- c(rownames(non_water_df), rownames(water_df))

sulf_otu_mat_final <- sulf_otu_mat_ordered[, sample_order]
sulf_sample_df_ordered <- sulf_sample_df[sample_order, ]

# Column labels: just the experiment number "#N" under each replicate, and
# nothing under the Base-Water control (its name is carried by the color band
# and the "Sample Type" legend).
batch_num <- gsub("pilot", "", sulf_sample_df_ordered$batch)
display_labels <- ifelse(
  sulf_sample_df_ordered$Label == "Base-Water",
  "", # Originally "Base-Water"
  # paste0("#", batch_num, ": ", sulf_sample_df_ordered$Label)
  paste0("#", batch_num)
)
names(display_labels) <- sample_order

# Row labels (italicized genus names)
row_labels <- sulf_tax_mat[sulf_taxa_order, "Genus"]

# Color band ("Sample Type") repeats the same color for replicates of the
# same Label
sample_type_vec <- sulf_sample_df_ordered$Label
stepped_cols_ordered <- stepped_cols[unique(sample_type_vec)]

# ------------------------------------------------------------------------ #
# 9. Bottom bubble data: total sulfur-bacteria abundance (%) per replicate,
#    matched to the heatmap's column order
# ------------------------------------------------------------------------ #

dot_vals <- (sulfur_totals[sample_order] * 100)
names(dot_vals) <- sample_order
dot_cols <- stepped_cols[sample_type_vec]
names(dot_cols) <- sample_order

max_diam_mm <- 11
max_limit <- max(10.1, max(dot_vals) * 1.01)  # keep >= 10.1 like the original ggplot scale

# ------------------------------------------------------------------------ #
# 10. Heatmap building blocks
# ------------------------------------------------------------------------ #

col_fun <- colorRamp2(c(0, max(sulf_otu_mat_final)), c("black", "#FEC437"))

row_ha <- rowAnnotation(
  `Metabolic Function` = sulf_tax_mat[sulf_taxa_order, "Functional_Group"],
  col = list(`Metabolic Function` = functional_colors),
  annotation_name_side = "top",
  show_annotation_name = FALSE,
  show_legend = FALSE
)

top_ha <- HeatmapAnnotation(
  `Sample Type` = sample_type_vec,
  col = list(`Sample Type` = stepped_cols_ordered),
  show_legend = FALSE,
  annotation_name_side = "right",
  show_annotation_name = FALSE
)

# Bubble annotation (bottom): area-scaled bubbles with the value above each
# bubble and the experiment number below it
bubble_center_y <- 0.5
bubble_panel_height_mm <- 22

bubble_anno <- AnnotationFunction(
  which = "column",
  n = ncol(sulf_otu_mat_final),
  height = unit(bubble_panel_height_mm, "mm"),
  fun = function(index) {
    vals <- dot_vals[index]
    cols <- dot_cols[index]
    labs <- display_labels[index]
    k <- length(index)

    pushViewport(viewport(xscale = c(0.5, k + 0.5), yscale = c(0, 1)))

    for (i in seq_len(k)) {
      # area-scaled; diameter proportional to sqrt(value)
      r_mm <- 0.5 * max_diam_mm * sqrt(vals[i] / max_limit)
      r_mm <- max(r_mm, 0)
      grid.circle(
        x = unit(i, "native"), y = unit(bubble_center_y, "native"),
        r = unit(r_mm, "mm"),
        gp = gpar(fill = cols[i], col = cols[i], lwd = 1.5)
      )
      # keep the value label from overflowing the top of the annotation strip
      label_y_native <- bubble_center_y + (r_mm + 3) / bubble_panel_height_mm
      label_y_native <- min(label_y_native, 0.98)
      grid.text(
        sprintf("%.2f", vals[i]),
        x = unit(i, "native"), y = unit(label_y_native, "native"),
        gp = gpar(fontsize = 8, fontface = "bold")
      )
    }

    # Column labels only below the dots, no labels on top
    for (i in seq_len(k)) {

      # only label rotated and adjusted to the right is the control
      # the rest are just #1, #2, #3
      # is_base_water <- labs[i] == "Base-Water"

      grid.text(
        labs[i],
        x = unit(i, "native"),
        y = unit(0.12, "native"),
        rot = 0,
        just = "center",
        # rot = ifelse(is_base_water,45, 0), # originally 45
        # just = ifelse(is_base_water, "right", "center"),
        gp = gpar(fontsize = 9, fontface = "bold")
      )
    }

    popViewport()
  }
)

bottom_ha <- HeatmapAnnotation(
  bubbles = bubble_anno,
  annotation_label = "Relative Abundance\nof Sulfur\nBacteria (%)",
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
  annotation_name_rot = 0,
  show_annotation_name = FALSE # we'll adjust the title manually
)

# ------------------------------------------------------------------------ #
# 11. The heatmap itself (genus level, individual replicates)
# ------------------------------------------------------------------------ #

ht_indiv <- Heatmap(
  sulf_otu_mat_final,
  name = "Relative Abundance\nwithin Sulfur Bacteria (%)",
  col = col_fun,

  # Rows
  row_labels = gt_render(paste0("*", row_labels, "*")),
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 9, fontface = "italic"),
  cluster_rows = FALSE,
  show_row_dend = FALSE,
  left_annotation = row_ha,

  # Columns -- no column names on the heatmap itself; labels come from
  # the bubble annotation below
  column_labels = display_labels,
  column_names_side = "bottom",
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  show_column_names = FALSE,
  top_annotation = top_ha,
  bottom_annotation = bottom_ha,
  show_heatmap_legend = FALSE,

  cell_fun = function(j, i, x, y, width, height, fill) {
    if (sulf_otu_mat_final[i, j] > 0.01) {
      grid.text(sprintf("%.1f", sulf_otu_mat_final[i, j] * 100),
                x, y, gp = gpar(fontsize = 8, col = "white", fontface = "bold"))
    }
  },

  width = unit(12, "cm"),
  height = unit(0.4 * nrow(sulf_otu_mat_final), "cm")
)

# ------------------------------------------------------------------------ #
# 12. Legends
# ------------------------------------------------------------------------ #

heatmap_legend <- Legend(
  col_fun = col_fun,
  title = "Relative Abundance\nwithin Sulfur Bacteria (%)",
  title_gp = gpar(fontsize = 9, fontface = "bold"),
  labels_gp = gpar(fontsize = 9),
  at = seq(0, max(sulf_otu_mat_final), length.out = 5),
  labels = c("0", "25", "50", "75", "100"),
  legend_height = unit(3.265, "cm"),
  direction = "vertical",
  title_gap = unit(10, "pt")
)

functional_legend <- Legend(
  labels = c("Sulfur Oxidizer", "Sulfate Reducer"),
  legend_gp = gpar(fill = functional_colors),
  title = "Metabolic Function",
  title_gp = gpar(fontsize = 9, fontface = "bold"),
  labels_gp = gpar(fontsize = 9),
  nrow = 2,
  title_gap = unit(5, "pt")
)

sample_legend <- Legend(
  labels = names(stepped_cols_ordered),
  legend_gp = gpar(fill = stepped_cols_ordered),
  title = "Sample Type",
  title_gp = gpar(fontsize = 9, fontface = "bold"),
  labels_gp = gpar(fontsize = 9),
  nrow = 5,
  title_gap = unit(5, "pt")
)

column2_legends <- packLegend(
  sample_legend,
  functional_legend,
  direction = "vertical",
  row_gap = unit(0.6, "cm")
)

legend_list <- packLegend(
  heatmap_legend,
  column2_legends,
  direction = "horizontal",
  column_gap = unit(0.9, "cm")
)

# ------------------------------------------------------------------------ #
# 13. Export the plotted data for inspection
#     One row per (genus x sample) heatmap cell. `RelAbund_within_sulfur` is
#     the cell value (fraction of the SOB+SRB community); `Total_sulfur_pct`
#     is the bubble value (% of the whole community) and is therefore repeated
#     down each sample's column.
# ------------------------------------------------------------------------ #

plot_df <- data.frame(
  Genus                 = rep(row_labels, times = ncol(sulf_otu_mat_final)),
  Functional_Group      = rep(sulf_tax_mat[sulf_taxa_order, "Functional_Group"],
                              times = ncol(sulf_otu_mat_final)),
  Sample                = rep(sample_order, each = nrow(sulf_otu_mat_final)),
  Batch                 = rep(sulf_sample_df_ordered$batch, each = nrow(sulf_otu_mat_final)),
  Sample_Type           = rep(sample_type_vec, each = nrow(sulf_otu_mat_final)),
  Column_Label          = rep(display_labels, each = nrow(sulf_otu_mat_final)),
  RelAbund_within_sulfur = as.vector(sulf_otu_mat_final),
  Total_sulfur_pct      = rep(dot_vals, each = nrow(sulf_otu_mat_final)),
  row.names = NULL,
  stringsAsFactors = FALSE
)

write.csv(plot_df,
          file = file.path(fig_dir, "Fig7_complexHeatmap_sulfur_bacteria_plotting_data.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------ #
# 14. Draw and save
# ------------------------------------------------------------------------ #

png(file.path(fig_dir, "Fig7_complexHeatmap_sulfur_bacteria.png"),
    width = 9.4, height = 2.7, units = "in", res = 300)
draw(ht_indiv,
     padding = unit(c(1, -80, 2, 2), "mm"))

decorate_annotation("bubbles", {
  grid.text("Relative Abundance\nof Sulfur Bacteria (%)",
            x = unit(0.9934, "npc") + unit(2.2, "mm"),
            y = unit(0.5, "npc"),
            just = "left",
            gp = gpar(fontsize = 9, fontface = "bold"))
})

pushViewport(
  viewport(
    x = unit(1, "npc"),
    just = "right",
    width = unit(80, "mm")
  )
)

draw(
  legend_list,
  x = unit(1.951, "npc"),
  y = unit(0.97, "npc"),
  just = c("right", "top")
)

popViewport()
dev.off()

cat("\nSaved 'Fig7_complexHeatmap_sulfur_bacteria.png' into", fig_dir, "\n")
