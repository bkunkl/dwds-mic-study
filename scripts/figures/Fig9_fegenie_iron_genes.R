# ======================================================================== #
# Fig. 8 -- FeGenie iron-related functional gene categories, all samples
#
# Implements the manuscript methods sentence:
#   "Iron-related genes were annotated with MagicLamp v1.0, integrating
#    curated HMMs from FeGenie (v1.2). Predicted proteins were aligned to the
#    NCBI nr database (release 2024-02-08) using DIAMOND BLASTp. Gene coverage
#    from depth files was normalized by cumulative sum scaling (75th
#    percentile) via metagenomeSeq (v1.50.0)."
#
# Each dot is one iron gene category in one sample; dot area is the
# depth-normalized gene coverage. `pilot1`/`pilot2`/`pilot3` correspond to
# experiments E1/E2/E3. Only E3 has the pre-incubation Base-Water control,
# which is plotted as its own single column at the right-hand end.
#
# INPUT (produced by the HPC pipeline, see hpc/ and scripts/local/03_run_fegenie.sh)
#   results/pilot?_ml100/assembly_megahit/S*_contigs.fa.depth.txt
#       per-sample contig depth tables from jgi_summarize_bam_contig_depths
#   results/pilot?_ml100/assembly_megahit/FeGenie_complete_results/FeGenie-heatmap-data.csv
#       FeGenie/MagicLamp gene-category x sample coverage table
#   data/metadata/metadata_bracken_conf01.csv
#       used only to turn "pilot1.S174" style column names into readable labels
#
# OUTPUT (into results/figures/)
#   Fig8_FeGenie_iron_related_functional_gene_categories.png
#   Fig8_FeGenie_iron_related_functional_gene_categories_plotting_data.csv
#   Fig8_FeGenie_CSS75_normalized_table.csv
#   plus, optionally (section 8), the condition-mean +/- SD companion figure
#
# NOTE: all 13 samples, including the E3 "Base-Water" pre-incubation control
# (S27), are CSS(75th percentile) normalized, so every dot area is on the same
# scale. An earlier version of this script dropped S27 from the depth-file list
# while keeping its FeGenie column, leaving it on raw rounded coverage.
#
# Run from the repository root.
# ======================================================================== #

library(metagenomeSeq)
library(stringr)
library(reshape2)
library(ggplot2)
library(paletteer)

# repo-relative paths -- run this script with the repo root as working directory
metadata_file <- "data/metadata/metadata_bracken_conf01.csv"

fig_dir <- "results/figures"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------ #
# 1. Per-sample CSS(75th percentile) normalization factors from depth files
# ------------------------------------------------------------------------ #

# Assemblies were repeated with reads filtered to >= 100 bp, hence "_ml100".
megahit_depth_files <- Sys.glob("results/pilot?_ml100/assembly_megahit/S*_contigs.fa.depth.txt")

# All 13 samples are globbed, S27 (the E3 "Base-Water" pre-incubation control)
# included, so every column of the FeGenie matrix gets its own normalization
# factor. S27 used to be filtered out here -- do not reintroduce that.

# Extract sample IDs from file paths
sample_ids <- str_extract(megahit_depth_files, "S\\d+")

# Read depth files and create MRexperiment objects to calculate
# normalization factors.
# Note: cumNormStat's auto-estimated quantile is always 0.5 here, because it
# expects multiple samples and each depth file is read on its own. So the
# 75th percentile is set explicitly instead.
megahit_depth_list <- lapply(seq_along(megahit_depth_files), function(d) {
    depth_m <- read.csv(megahit_depth_files[d], sep = "\t", row.names = 1)
    depth_m <- as.matrix(depth_m[, 2, drop = FALSE])
    colnames(depth_m) <- sample_ids[d]  # Use sample ID as column name
    m_obj <- newMRexperiment(depth_m)
    m_obj <- cumNorm(m_obj, p = 0.75)
    nf <- normFactors(m_obj)[1]  # Extract just the numeric value without column name
    return(list(
      depth_matrix = depth_m,
      norm_factor = nf
    ))
})

# Assign sample IDs as names to the list elements
# Better not, it doubles names in the norm factors named vectors:
# "S217.S217"
# names(megahit_depth_list) <- sample_ids

norm_factors_75 <- sapply(megahit_depth_list, function(x) x$norm_factor)

# ------------------------------------------------------------------------ #
# 2. Combine the three FeGenie heatmap tables into one count matrix
# ------------------------------------------------------------------------ #

megahit_files <- Sys.glob("results/pilot?_ml100/assembly_megahit/FeGenie_complete_results/FeGenie-heatmap-data.csv")
pilot_names <- paste0("pilot", 1:3)
megahit_list <- lapply(seq_along(megahit_files), function(i) {
  hmap_dat <- read.csv(megahit_files[i], row.names=1)
  hmap_dat <- hmap_dat[, -ncol(hmap_dat)]
  colnames(hmap_dat) <- gsub("_contigs.fa", "", colnames(hmap_dat))
  round(hmap_dat)
})
names(megahit_list) <- pilot_names

# Reorder columns so every experiment reads
# Biofilm-no-nutri, Biofilm-w-nutri, Bulk-no-nutri, Bulk-w-nutri
megahit_list$pilot2 <- megahit_list$pilot2[,c(1,3,2,4)]
# Let's put the BaseWater sample at the end of the table
# That way we have to modify the code the least
megahit_list$pilot3 <- megahit_list$pilot3[,c(1,4,3,5,2)]

merged_megahit <- do.call(cbind, megahit_list)
merged_megahit <- as.matrix(merged_megahit)

# Extract the sample IDs from FeGenie column names ("pilot1.S174" -> "S174")
fegenie_sample_ids <- sapply(strsplit(colnames(merged_megahit), "\\."), function(x) x[2])
print(fegenie_sample_ids)  # Should now correctly show S174, S175, etc.
new_order <- match(fegenie_sample_ids, sample_ids)

# ------------------------------------------------------------------------ #
# 3. Apply the depth-based normalization factors
# ------------------------------------------------------------------------ #

normalize_FeGenie_with_depth <- function(fegenie_matrix, norm_factors, sample_mapping) {
  # Create a copy of the original matrix
  normalized_matrix <- fegenie_matrix

  # Apply normalization factors to each column. Both `if`s are guards: a
  # column with no matching normalization factor would fall through and be
  # carried over as RAW rounded coverage, silently off-scale from the rest.
  # With section 1 globbing all 13 depth files, nothing falls through -- if
  # anything ever does, the sample is missing its depth file.
  for (col in colnames(fegenie_matrix)) {
    if (col %in% names(sample_mapping)) {
      depth_sample <- sample_mapping[[col]]
      if (depth_sample %in% names(norm_factors)) {
        normalized_matrix[, col] <- round((fegenie_matrix[, col] / norm_factors[depth_sample])*10^4)
      }
    }
  }

  return(normalized_matrix)
}

# Create mapping between FeGenie column names and depth sample IDs
sample_mapping <- setNames(
  sample_ids[new_order],
  colnames(merged_megahit)
)

# Apply normalization using the 75th percentile factors
# Most straightforward, use this
normalized_megahit_75 <- normalize_FeGenie_with_depth(
  merged_megahit,
  norm_factors_75,
  sample_mapping
)

# ------------------------------------------------------------------------ #
# 4. Readable column labels from the sample metadata
# ------------------------------------------------------------------------ #

metadata_table <- read.csv(metadata_file, stringsAsFactors = FALSE)
metadata_table$combined <- paste(metadata_table$batch, metadata_table$sample_number, sep = ".")
meta_ordered <- metadata_table[match(colnames(normalized_megahit_75), metadata_table$combined),]

# Check for any NA values (possible mismatch)
if(any(is.na(match(colnames(normalized_megahit_75), metadata_table$combined)))) {
  warning("Some sample names don't match, check closer!")
}

colnames(normalized_megahit_75) <- paste(meta_ordered$batch, meta_ordered$Label, sep = ":")
colnames(normalized_megahit_75) <- gsub(pattern = "pilot", replacement = "#",colnames(normalized_megahit_75))
# Remove batch number from Base-Water
colnames(normalized_megahit_75)[which(grepl("Base-Water", colnames(normalized_megahit_75)))] <- "Base-Water"

write.csv(normalized_megahit_75,
          file = file.path(fig_dir, "Fig8_FeGenie_CSS75_normalized_table.csv"))

# ------------------------------------------------------------------------ #
# 5. Helpers and palette
# ------------------------------------------------------------------------ #

# Capitalize firts letter function:
capitalize_first_letter <- function(x) {
  paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))
}

charizard_palette <- paletteer_d("palettetown::charizard")[c(1:2,5:14)]
charizard_palette[11] <- "#B04418FF"
# replace light grey for iron reduction ("#708238" or "#556B2F"):
charizard_palette[10] <- "#4682B4"

# ------------------------------------------------------------------------ #
# 6. Long-format data frame for plotting (all samples, no aggregation)
# ------------------------------------------------------------------------ #

generate_full_df_for_plotting <- function(genie_matrix) {

  FeCategs <- rownames(genie_matrix)
  FeCategs <- gsub("_", " ", FeCategs)
  FeCategs <- gsub("iron aquisition-", "", FeCategs)
  FeCategs <- sapply(FeCategs, capitalize_first_letter)
  # NOTE: positional -- row 8 of FeGenie-heatmap-data.csv. Check that the row
  # order of the FeGenie table has not changed before re-running.
  FeCategs[8] <- "Possible iron oxidation and reduction"

  df_all <- data.frame(IronCategory = FeCategs, genie_matrix, check.names = FALSE)

  df_all_melt <- melt(
    df_all,
    id.vars = "IronCategory",
    measure.vars = colnames(genie_matrix),
    variable.name = "Sample",
    value.name = "GeneCoverage"
  )

  # keep column order from normalized_megahit_75 as x-axis label order
  df_all_melt$Sample <- factor(df_all_melt$Sample, levels = colnames(genie_matrix))

  df_all_melt$IronCategory <- factor(
    df_all_melt$IronCategory,
    levels = unique(df_all_melt$IronCategory)
  )

  return(list(melt_df = df_all_melt, full_df = df_all))
}

FeGenie_full_CSS_75_megahit <- generate_full_df_for_plotting(normalized_megahit_75)

# Export the final data frame that goes into the plot, for inspection
write.csv(FeGenie_full_CSS_75_megahit$melt_df,
          file = file.path(fig_dir,
                           "Fig8_FeGenie_iron_related_functional_gene_categories_plotting_data.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------ #
# 7. Fig. 8 -- dot plot over all samples
# ------------------------------------------------------------------------ #

FeGenie_dotplotter_allsamples <- function(FeGenie_df_melt,
                                           size_legend_title = "Normalized Gene Coverage",
                                           outfile = file.path(fig_dir, "Fig8_FeGenie_iron_related_functional_gene_categories.png")) {

  n_samples <- length(unique(FeGenie_df_melt$Sample))
  width_unit <- n_samples * 1.2

  FG_dp_all <- ggplot(FeGenie_df_melt, aes(x = Sample, y = IronCategory, size = GeneCoverage)) +
    geom_point(aes(color = IronCategory), alpha = 0.75) +
    scale_size_area(max_size = 20) +
    scale_color_manual(values = c("Iron transport" = charizard_palette[1],
                                  "Heme transport" = charizard_palette[2],
                                  "Heme oxygenase" = charizard_palette[3],
                                  "Siderophore synthesis" = charizard_palette[4],
                                  "Siderophore transport potential" = charizard_palette[5],
                                  "Iron gene regulation" = charizard_palette[6],
                                  "Iron oxidation" = charizard_palette[7],
                                  "Possible iron oxidation and reduction" = charizard_palette[8],
                                  "Probable iron reduction" = charizard_palette[9],
                                  "Iron reduction" = charizard_palette[10],
                                  "Iron storage" = charizard_palette[11],
                                  "Magnetosome formation" = charizard_palette[12])) +
    labs(x = "", y = "Iron Category", size = size_legend_title, color = "") +
    theme(panel.background = element_rect(fill = "white", colour = "white", linewidth = 1, linetype = 0),
          panel.border = element_rect(colour="black", linewidth=1, fill=NA),
          strip.background=element_rect(fill='white', colour='white'),
          strip.text = element_text(face="bold", size=10),
          panel.grid.major = element_line(linewidth = 0),
          panel.grid.minor = element_line(linewidth = 0),
          axis.text = element_text(size=15, colour="black", face = "bold"),
          axis.title = element_text(face="bold", size=17),
          axis.text.x = element_text(vjust = 1, angle = 45, color = "black", hjust=1),
          axis.title.y = element_text(margin = margin(r = 20)),
          legend.position.inside = c(1,1),
          legend.title = element_text(size = 14),
          legend.text = element_text(size=14),
          legend.justification = c(1,1)) +
    guides(color = "none",
           size = guide_legend(order = 1))

  png(outfile,
      width = width_unit,
      height = 10,
      units = "in",
      res = 150)
  print(FG_dp_all)
  dev.off()

  return(FG_dp_all)
}

FG_dp_all_plot_object <- FeGenie_dotplotter_allsamples(
  FeGenie_df_melt = FeGenie_full_CSS_75_megahit$melt_df,
  size_legend_title = "Normalized Gene Coverage",
  outfile = file.path(fig_dir, "Fig8_FeGenie_iron_related_functional_gene_categories.png")
)

FG_dp_all_plot_object

# ======================================================================== #
# 8. Companion figure (not in the manuscript, kept because it is the
#    condition-level summary of exactly the same normalized data)
#
#    Same CSS-75 matrix, but averaged over the three experiments per
#    condition. The grey dot behind each colored dot is mean + SD, so the
#    visible grey ring is the +1 SD margin. Base-Water has a single
#    observation, so its SD dot is set to 0.
#
#    After experimenting with the first version of this analysis (without the
#    Base-Water sample) this seemed the best among the available methods.
#    ANCOMBC would be better, but it needs at least two replicates in each
#    category (including Base-Water).
# ======================================================================== #

generate_2dfs_for_plotting <- function(genie_matrix) {

  FeCategs <- rownames(genie_matrix)

  FeCategs <- gsub("_", " ", FeCategs)
  FeCategs <- gsub("iron aquisition-", "", FeCategs)
  FeCategs <- sapply(FeCategs, capitalize_first_letter)
  FeCategs[8] <- "Possible iron oxidation and reduction"

  df_means <- data.frame(
    IronCategory = FeCategs,
    Biofilm_no_nutri = apply(genie_matrix[, c(1, 5, 9)], 1, mean),
    Biofilm_w_nutri = apply(genie_matrix[, c(2, 6, 10)], 1, mean),
    Bulk_no_nutri = apply(genie_matrix[, c(3, 7, 11)], 1, mean),
    Bulk_w_nutri = apply(genie_matrix[, c(4, 8, 12)], 1, mean),
    Base_Water = genie_matrix[,13]
  )
  df_sds <- data.frame(
    IronCategory = FeCategs,
    Biofilm_no_nutri = apply(genie_matrix[, c(1, 5, 9)], 1, mean) + apply(genie_matrix[, c(1, 5, 9)], 1, sd),
    Biofilm_w_nutri = apply(genie_matrix[, c(2, 6, 10)], 1, mean) + apply(genie_matrix[, c(2, 6, 10)], 1, sd),
    Bulk_no_nutri = apply(genie_matrix[, c(3, 7, 11)], 1, mean) + apply(genie_matrix[, c(3, 7, 11)], 1, sd),
    Bulk_w_nutri = apply(genie_matrix[, c(4, 8, 12)], 1, mean) + apply(genie_matrix[, c(4, 8, 12)], 1, sd),
    Base_Water = 0
  )

  df_means_melt <- melt(
    df_means,
    id.vars = "IronCategory",
    measure.vars = c(
      "Biofilm_no_nutri", "Biofilm_w_nutri",
      "Bulk_no_nutri", "Bulk_w_nutri", "Base_Water"
    ),
    variable.name = "Sample",
    value.name = "GeneCoverage"
  )

  df_sds_melt <- melt(
    df_sds,
    id.vars = "IronCategory",
    measure.vars = c(
      "Biofilm_no_nutri", "Biofilm_w_nutri",
      "Bulk_no_nutri", "Bulk_w_nutri", "Base_Water"
    ),
    variable.name = "Sample",
    value.name = "GeneCoverage"
  )

  df_means_melt$Sample <- gsub("_", "-", df_means_melt$Sample)
  df_means_melt$Sample <- gsub("w", "w/", df_means_melt$Sample)
  df_sds_melt$Sample <- gsub("_", "-", df_sds_melt$Sample)
  df_sds_melt$Sample <- gsub("w", "w/", df_sds_melt$Sample)

  df_means_melt$Sample <- factor(
    df_means_melt$Sample,
    levels = c("Biofilm-no-nutri", "Biofilm-w/-nutri", "Bulk-no-nutri", "Bulk-w/-nutri", "Base-Water")
  )
  df_sds_melt$Sample <- factor(
    df_sds_melt$Sample,
    levels = c("Biofilm-no-nutri", "Biofilm-w/-nutri", "Bulk-no-nutri", "Bulk-w/-nutri", "Base-Water")
  )

  df_means_melt$IronCategory <- factor(
    df_means_melt$IronCategory,
    levels = unique(df_means_melt$IronCategory)
  )
  df_sds_melt$IronCategory <- factor(
    df_sds_melt$IronCategory,
    levels = unique(df_sds_melt$IronCategory)
  )

  return(list(melt_means_df = df_means_melt, melt_sds_df = df_sds_melt,
              means_df = df_means, sds_df = df_sds))
}

FeGenie_dotplotter <- function(FeGenie_df_sds_melt,
                               FeGenie_df_means_melt,
                               size_legend_title = "Normalized Gene Coverage",
                               outfile = file.path(fig_dir, "FeGenie_iron_gene_categories_condition_means.png")){

  n_samples <- length(unique(FeGenie_df_means_melt$Sample))
  width_unit <- n_samples*2.1

  FG_dp <- ggplot(FeGenie_df_sds_melt, aes(x = Sample, y = IronCategory, size = GeneCoverage)) +
    geom_point(aes(color = "Standard Deviation", size = GeneCoverage), alpha = 0.75, show.legend = TRUE) +
    scale_size_area(max_size = 20) +
    geom_point(data = FeGenie_df_means_melt, aes(color = IronCategory), show.legend = FALSE)+
    scale_color_manual(values = c("Standard Deviation" = "grey40",
                                  "Iron transport" = charizard_palette[1],
                                  "Heme transport" = charizard_palette[2],
                                  "Heme oxygenase" = charizard_palette[3],
                                  "Siderophore synthesis" = charizard_palette[4],
                                  "Siderophore transport potential" = charizard_palette[5],
                                  "Iron gene regulation" = charizard_palette[6],
                                  "Iron oxidation" = charizard_palette[7],
                                  "Possible iron oxidation and reduction" = charizard_palette[8],
                                  "Probable iron reduction" = charizard_palette[9],
                                  "Iron reduction" = charizard_palette[10],
                                  "Iron storage" = charizard_palette[11],
                                  "Magnetosome formation" = charizard_palette[12]),
                       breaks = "Standard Deviation")+
    labs(x = "", y = "Iron Category", size = size_legend_title, color = "")+
    theme(panel.background = element_rect(fill = "white", colour = "white", linewidth = 1, linetype = 0),
          panel.border = element_rect(colour="black", linewidth=1, fill=NA),
          strip.background=element_rect(fill='white', colour='white'),
          strip.text = element_text(face="bold", size=10),
          panel.grid.major = element_line(linewidth = 0),
          panel.grid.minor = element_line(linewidth = 0),
          axis.text = element_text(size=15, colour="black", face = "bold"),
          axis.title = element_text(face="bold", size=17),
          axis.text.x = element_text(vjust = 1, angle = 45, color = "black", hjust=1),
          axis.title.y = element_text(margin = margin(r = 20)),
          legend.position.inside = c(1,1),
          legend.title = element_text(size = 14),
          legend.text = element_text(size=14),
          legend.justification = c(1,1))+
    guides(size = guide_legend(order = 1),
           color = guide_legend(order = 2, override.aes = list(size = 4)))

  png(outfile,
      width = width_unit,
      height = 10,
      units = "in",
      res = 150)
  print(FG_dp)
  dev.off()

  return(FG_dp)
}

FeGenie_CSS_75_megahit <- generate_2dfs_for_plotting(normalized_megahit_75)

FG_dp_plot_object <- FeGenie_dotplotter(
  FeGenie_df_sds_melt = FeGenie_CSS_75_megahit$melt_sds_df,
  FeGenie_df_means_melt = FeGenie_CSS_75_megahit$melt_means_df,
  size_legend_title = "Normalized Gene Coverage",
  outfile = file.path(fig_dir, "FeGenie_iron_gene_categories_condition_means.png")
)

write.csv(FeGenie_CSS_75_megahit$melt_means_df,
          file = file.path(fig_dir, "FeGenie_iron_gene_categories_condition_means_plotting_data.csv"),
          row.names = FALSE)

cat("\nSaved 'Fig8_FeGenie_iron_related_functional_gene_categories.png' into", fig_dir, "\n")
