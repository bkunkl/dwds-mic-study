# ======================================================================== #
# Fig. 5 and Figs. S5-S8 -- alpha diversity of the pilot 1-3 metagenomes
#
# Implements the manuscript's alpha-diversity analysis:
#   Hill numbers (q = 0 richness, q = 1 Shannon, q = 2 Simpson) estimated
#   with iNEXT at coverage-based standardization, Pielou's evenness derived
#   from the q=0 and q=1 Hill numbers, and a diversity profile over
#   q = 0 ... 8.
#
# `pilot1`/`pilot2`/`pilot3` correspond to experiments E1/E2/E3. Two factors
# at two levels each: Nutrient (no-nutri / w-nutri) and Location
# (bulk = planktonic / biofilm = carbon steel coupon surface). E1 and E3 ran
# 35 days, E2 ran 60 days. Only E3 has the pre-incubation "Base-Water"
# control (S27), so that group has n = 1 -- which is why it is used as the
# model reference level but is dropped from the Cohen's d contrasts.
# 13 samples in total (4 + 4 + 5).
#
# ------------------------------------------------------------------------ #
# THIS SCRIPT RUNS IN TWO STAGES, WITH THE HPC JOBS IN BETWEEN.
#
#   STAGE 1  needs only the phyloseq object. It writes the OTU matrix and
#            the sample metadata table that the cluster jobs consume.
#   ...      submit and merge the iNEXT jobs (see the STOP banner below)
#   STAGE 2  needs the merged diversity profiles. It makes all five figures.
#
# The two stages are separated by a banner: run stage 1, submit the cluster
# jobs, then run stage 2 once their output has been merged.
# ------------------------------------------------------------------------ #
#
# INPUT (stage 1)
#   results/kraken_bracken/brackenSpecies_conf01.rds
#       phyloseq object built by scripts/local/01_build_phyloseq_object.sh
#
# INPUT (stage 2)
#   results/diversity_profile/unified_extended_diversity_profiles.rds
#       merged q x sample Hill-number table from
#       scripts/local/05_merge_diversity_profiles.R
#   data/kraken2/kraken2_classification_rates.csv
#       per-sample Kraken2 classification rates, used as a covariate in the
#       model comparison of section 10
#
# OUTPUT (stage 1)
#   results/inext_rarefaction/pilot123_diversity.RData
#   results/pilot123_diversity_metadata_table.rds
#
# OUTPUT (stage 2, into results/figures/)
#   Fig5_alpha_diversity_metrics.png                      (+ 2 CSVs)
#   FigS5_model_coefficients_across_alpha_diversity_metrics.png   (+ CSV)
#   FigS6_standardized_treatment_contrasts_cohens_d.png           (+ CSV)
#   FigS7_alpha_diversity_profile.png                             (+ CSV)
#   FigS8_alpha_diversity_profile_derivative.png                  (+ CSV)
#   plus, from section 10, the non-manuscript companion analyses
#
# Run from the repository root.
# ======================================================================== #

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggrepel)
library(gridExtra)
library(effectsize)
# pals::, scales:: and ggpubr:: are called with their namespace prefix only,
# so they need to be installed but not attached.
# multcomp and broom are attached inside section 10 only -- multcomp pulls in
# MASS, and nothing above section 10 needs either of them.

# repo-relative paths -- run this script with the repo root as working directory
phyloseq_file  <- "results/kraken_bracken/brackenSpecies_conf01.rds"
div_prof_file  <- "results/diversity_profile/unified_extended_diversity_profiles.rds"
class_rate_file <- "data/kraken2/kraken2_classification_rates.csv"

fig_dir <- "results/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("results/inext_rarefaction", recursive = TRUE, showWarnings = FALSE)


# ######################################################################## #
#                                                                          #
#  STAGE 1 -- RUNS BEFORE THE HPC JOBS                                     #
#                                                                          #
#  Everything from here down to the STOP banner needs only the phyloseq    #
#  object. It writes the two files that the iNEXT cluster chain and the    #
#  local merge scripts read.                                               #
#                                                                          #
# ######################################################################## #

# ------------------------------------------------------------------------ #
# 1. Phyloseq object -> OTU matrix + sample metadata
# ------------------------------------------------------------------------ #

ps_c01 <- readRDS(phyloseq_file)
ps_c01_prok <- subset_taxa(ps_c01, Domain == "Bacteria" | Domain == "Archaea")

# Sample labels and batches, in the order used by every figure below.
# Note: "Base-Water" is the E3-only pre-incubation control, so it is last.
meta_c1 <- as(sample_data(ps_c01_prok), "data.frame")
meta_c1$Label <- factor(meta_c1$Label,
                         levels = c("Biofilm-no-nutri", "Biofilm-w/-nutri",
                                    "Bulk-no-nutri", "Bulk-w/-nutri", "Base-Water"))

meta_c1$batch <- factor(meta_c1$batch,
                        levels = unique(meta_c1$batch))

# iNEXT wants species as rows.
otu_matrix <- as(otu_table(ps_c01_prok), "matrix")
if (!taxa_are_rows(ps_c01_prok)) {
  otu_matrix <- t(otu_matrix)  # Transpose so species are rows
}

# Check sequencing depths (the range drives the iNEXT endpoint choice made
# in scripts/hpc/03_inext_rarefaction_worker.R: 2x the median depth).
sample_depths <- colSums(otu_matrix)
print(range(sample_depths))
print(median(sample_depths))

# ------------------------------------------------------------------------ #
# 2. Save the inputs the HPC iNEXT chain needs  (SAVE POINT)
#
# Only the two objects the cluster jobs and the merge scripts read are written
# out, at the point where they become available.
# ------------------------------------------------------------------------ #

# read by scripts/hpc/03_inext_rarefaction_worker.R (and by the diversity-profile
# worker described in scripts/hpc/README.md), which load() this file and index into
# otu_matrix's columns with the SLURM array task ID
save(otu_matrix, meta_c1, file = "results/inext_rarefaction/pilot123_diversity.RData")

# read by scripts/local/04_merge_rarefaction_results.R and
# scripts/local/05_merge_diversity_profiles.R -- they take the sample name, batch and
# Label of column N of otu_matrix from row N of this table, so the row order
# of meta_c1 must match the column order of otu_matrix.
saveRDS(meta_c1, "results/pilot123_diversity_metadata_table.rds")


# ######################################################################## #
#                                                                          #
#  >>> STOP HERE. THE HPC JOBS RUN NOW. <<<                                #
#                                                                          #
#  Nothing below this line can run until the cluster jobs have finished    #
#  and their output has been merged. In this order:                        #
#                                                                          #
#    1. scripts/hpc/03_inext_rarefaction_array.sbatch   -- 13 samples x q=0,1,2     #
#                                                 = 39 jobs, ~1 day each    #
#    2. scripts/hpc/04_submit_diversity_profile_jobs.sh -- one array per q value,   #
#                                                 35 q values x 13 samples  #
#    3. scripts/local/04_merge_rarefaction_results.R    -- per q value (0, 1, 2)    #
#    4. scripts/local/05_merge_diversity_profiles.R     -- once, over the whole     #
#                                                 q sweep                   #
#                                                                          #
#  Then come back and run STAGE 2 below. Stage 2 does not depend on any     #
#  object created in stage 1, so it can also be run on its own in a fresh   #
#  session (run the library and path block at the top first).              #
#                                                                          #
# ######################################################################## #


# ######################################################################## #
#                                                                          #
#  STAGE 2 -- NEEDS THE MERGED HPC OUTPUT                                  #
#                                                                          #
#  Input: results/diversity_profile/unified_extended_diversity_profiles.rds #
#  written by scripts/local/05_merge_diversity_profiles.R. One row per             #
#  (q value x sample), columns SC, m, Method, Order.q, qD, qD.LCL,         #
#  qD.UCL, sample_name, batch, Label.                                      #
#                                                                          #
#  All five manuscript figures are produced from here down.                #
#                                                                          #
# ######################################################################## #

# ------------------------------------------------------------------------ #
# 3. Load the merged diversity profiles
# ------------------------------------------------------------------------ #

div_profiles <- readRDS(div_prof_file)

# Define sample type colors
stepped_cols <- pals::stepped(n=24)[c(5,7,13,15,22)]

# ------------------------------------------------------------------------ #
# 4. Per-metric tables derived from the diversity profiles
#
# Two sets of subsets are taken from the same table:
#   *_df / *_data -- sample depth, sample name, batch, Label, qD. Used by
#                    the depth/treatment models in section 6 (Fig. S5).
#   *_iNEXT       -- qD, batch, Label only. Used by the boxplots in section
#                    5 (Fig. 5) and the effect sizes in section 7 (Fig. S6).
#
# NOTE ON ORDER: the *_iNEXT block is kept here, immediately before its first
# use in the effect-size section, so the script runs top to bottom.
# ------------------------------------------------------------------------ #

# retain only q values 0,1, 2
i_q0 <- which(div_profiles$Order.q %in% 0)
richness_df <- div_profiles[i_q0, c(2, 5, 8:10)]
richness_df$Diversity_Index <- "Richness"

i_q1 <- which(div_profiles$Order.q %in% 1)
shannon_df <- div_profiles[i_q1, c(2, 5, 8:10)]
shannon_df$Diversity_Index <- "Shannon"

i_q2 <- which(div_profiles$Order.q %in% 2)
simpson_df <- div_profiles[i_q2, c(2, 5, 8:10)]
simpson_df$Diversity_Index <- "Simpson"

# change column name:
# List of data frame names
df_names <- c("richness_df", "shannon_df", "simpson_df")

# Loop through and rename
for (df_name in df_names) {
  df <- get(df_name)
  names(df)[names(df) == "m"] <- "sample_depth"
  df$sample_name <- gsub("_conf01_k2$", "", df$sample_name)
  assign(df_name, df)
}

# kraken2 classification rates -- the fraction of reads Kraken2 could assign
# to any taxon, read off the run summaries and kept alongside the other
# hand-transcribed data (see data/kraken2/README.md). Originally an .rds.
# `Sample` matches sample_name after the "_conf01_k2" strip above.
# Only the model comparison in section 10 uses the classification_rate
# column; the models behind Fig. S5 do not. The merges below are kept
# exactly as they ran, for consistency.
classification_rate_df <- read.csv(class_rate_file, stringsAsFactors = FALSE)

# Example merge (adjust column names as needed):
richness_data <- merge(richness_df, classification_rate_df, 
                       by.x = "sample_name", by.y = "Sample", 
                       all.x = TRUE)

# Use Base-Water as reference in the models:
richness_data$Label <- relevel(richness_data$Label, ref = "Base-Water")

# Interaction Models for All Diversity Indices ####
# Prepare data frames for Shannon, inverse Simpson and Pielou
# classification rate is no longer important after we saw it does not
# have significant effect, but for consistency's sake, let's keep it.
shannon_data <- merge(shannon_df, classification_rate_df, 
                       by.x = "sample_name", by.y = "Sample", 
                       all.x = TRUE)

simpson_data <- merge(simpson_df, classification_rate_df,
                      by.x = "sample_name", by.y = "Sample",
                      all.x = TRUE)

pielou_data <- richness_data
pielou_data$qD <- log(shannon_data$qD)/log(richness_data$qD)

# Use Base-Water as reference in the models:
shannon_data$Label <- relevel(shannon_data$Label, ref = "Base-Water")
simpson_data$Label <- relevel(simpson_data$Label, ref = "Base-Water")
pielou_data$Label <- relevel(pielou_data$Label, ref = "Base-Water")

# Hill-number subsets for the boxplots and the effect sizes
i_q1 <- which(div_profiles$Order.q %in% 1)
shannon_iNEXT <- div_profiles[i_q1, c(5, 9:10)]
shannon_iNEXT$Diversity_Index <- "Shannon"

i_q0 <- which(div_profiles$Order.q %in% 0)
richness_iNEXT <- div_profiles[i_q0, c(5, 9:10)]
richness_iNEXT$Diversity_Index <- "Richness"

i_q2 <- which(div_profiles$Order.q %in% 2)
simpson_iNEXT <- div_profiles[i_q2, c(5, 9:10)]
simpson_iNEXT$Diversity_Index <- "Simpson"


# Pielou's evenness
# Watch out, iNEXT gives exp(Shannon entropy) so we have to take the logarithm of that too.
pielou_evenness <- log(shannon_iNEXT$qD) / log(richness_iNEXT$qD)
pielou_iNEXT <- data.frame(
  qD = pielou_evenness,
  batch = richness_iNEXT$batch,
  Label = richness_iNEXT$Label,
  Diversity_Index = "Pielou"
)

# Pielou's evenness is not on the same scale as the other three
# They must be plotted separately
# Previously I did it altogether in the div_4_indices objects, but overwrote them.
div_3_indices <- rbind(richness_iNEXT,
                       shannon_iNEXT,
                       simpson_iNEXT)


div_3_indices$Diversity_Index <- factor(div_3_indices$Diversity_Index,
                                           level = c("Richness", "Shannon",
                                                      "Simpson"))

pielou_iNEXT$Diversity_Index <- "Pielou"

# ------------------------------------------------------------------------ #
# 5. Fig. 5 -- alpha diversity metrics by sample type
#
# Three Renyi entropies (H0 = log richness, H1 = Shannon, H2 = log Simpson)
# on a shared scale in one panel row, plus Pielou's evenness in its own
# panel because it is on a different scale. The batch label ("#1", "#2",
# "#3", "Ctrl") is repelled next to each point so the reader can see which
# experiment each replicate came from.
# ------------------------------------------------------------------------ #

# Boxplot Rényi ####
# Transform Hill numbers to Rényi entropies
div_3_indices_renyi <- div_3_indices
div_3_indices_renyi$qD_renyi <- log(div_3_indices$qD)

# Update diversity index names to reflect transformation
div_3_indices_renyi$Diversity_Index <- factor(
  div_3_indices_renyi$Diversity_Index,
  levels = c("Richness", "Shannon", "Simpson"),
  labels = c("Rényi H₀ (log Richness)", "Rényi H₁ (Shannon)", "Rényi H₂ (log Simpson)")
)


# Create boxplots with Rényi entropies
#  fixate jitter points first:
p_div_3_boxplots_renyi <- ggplot(div_3_indices_renyi, aes(x = Label, y = qD_renyi, fill = Label)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_point(aes(color = Label), 
             position = position_jitter(width = 0.2, seed = 2), 
             size = 3, 
             alpha = 0.8) +
  facet_wrap(~ Diversity_Index, scales = "fixed", nrow = 1) +
  scale_fill_manual(values = stepped_cols) +
  scale_color_manual(values = stepped_cols) +
  labs(
    # title = "Alpha Diversity Metrics (Rényi Entropies) by Sample Type",
    # x = "Sample Type",
    y = "Rényi Entropy"
  ) +
  theme_bw() +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 13), # later added
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),#12
    strip.text = element_text(face = "bold", size = 14), #12
    axis.title.y = element_text(face = "bold", size = 14), #12
    axis.title.x = element_blank(),
    # plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14) #11
  ) +
  # for combined plot, otherwise legend.title = element_blank()
  guides(fill = guide_legend(nrow = 1, title = "Sample"),
         color = guide_legend(nrow = 1, title = "Sample"))

# Export the data frame that went into the plot object above, for inspection.
# Written here, before the Base-Water relabelling below, because that is the
# state ggplot() captured (see the note on that step).
write.csv(div_3_indices_renyi,
          file = file.path(fig_dir, "Fig5_alpha_diversity_metrics_plotting_data.csv"),
          row.names = FALSE)

# add geom_label_repel to the data points in the boxplot:
#
# CAUTION -- this blanking has no effect on the figure. ggplot() has already
# copied div_3_indices_renyi into p_div_3_boxplots_renyi above, and the
# geom_label_repel layer added below inherits that stored copy, not the
# modified data frame. The Base-Water point therefore keeps its batch label
# ("Ctrl") in the published figure. Removing the blanking, or moving it above
# the ggplot() call, changes the figure.
bw_i <- which(div_3_indices_renyi$Label == "Base-Water")
if(length(bw_i) > 0) {
  div_3_indices_renyi$batch[bw_i] <- ""
}

p_div_3_boxplots_renyi_repel <- p_div_3_boxplots_renyi +
  geom_label_repel(aes(label = batch),
                   position = position_jitter(width = 0.2, seed = 2),
                   size = 4,
                   max.overlaps = 5,
                   show.legend = FALSE,
                   # nudge_x = 0.15,
                   # nudge_y = 0.1,
                   label.r = 0.2,
                   label.padding = 0.2,
                   label.size = 0.2, #0.2 originally
                   fill = "white")

# Generate Pielou's evennes plotting object separately,
# but imitate the renyi boxplots, they are meant to be on a composit plot
pielou_iNEXT$Diversity_Index <- factor(
  pielou_iNEXT$Diversity_Index,
  levels = c("Pielou"),
  labels = c("Pielou's Evenness")
)

p_pielou <- ggplot(data = pielou_iNEXT, aes(x = Label, y = qD, fill = Label))+
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  geom_point(aes(color = Label),
             position = position_jitter(width = 0.2, seed = 2),
             size = 3,
             alpha = 0.8)+
  facet_wrap(~ Diversity_Index, scales = "free_y", nrow = 1)+
  scale_fill_manual(values = stepped_cols)+
  scale_color_manual(values = stepped_cols)+
  labs(y = "Evenness Index")+
  # NOTE: hardcoded y limits, chosen to make this panel the same height as the
  # Renyi panels. On re-run with different data any point outside this window
  # would be dropped (with a warning) rather than rescaling the axis.
  scale_y_continuous(limits = c(0.045,0.788))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 13), # later added
        strip.text = element_text(face = "bold", size = 14),#12
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 14))+ #12
  guides(fill = "none", color = "none")

write.csv(pielou_iNEXT,
          file = file.path(fig_dir, "Fig5_alpha_diversity_metrics_pielou_plotting_data.csv"),
          row.names = FALSE)

# add repel labels to pielou boxplot:
# (same caveat as above -- p_pielou already holds its own copy of the data)
bw_i <- which(pielou_iNEXT$Label == "Base-Water")
if(length(bw_i) > 0) {
  pielou_iNEXT$batch[bw_i] <- ""
}

p_pielou_repel <- p_pielou +
  geom_label_repel(aes(label = batch),
                   position = position_jitter(width = 0.2, seed = 2),
                   size = 4,
                   max.overlaps = 5,
                   show.legend = FALSE,
                   # nudge_x = 0.15,
                   # nudge_y = 0.01,
                   label.r = 0.2,
                   label.padding = 0.2,
                   label.size = 0.2, #0.2
                   fill = "white")

# The legend of the combined plot makes it unproportional, so it is built
# separately with wider spacing and glued underneath the two panels.
# Instead of extracting, create a new legend with proper spacing
p_legend_source <- p_div_3_boxplots_renyi + 
  theme(
    legend.position = "bottom",
    legend.spacing.x = unit(1.8, "cm"),      # Horizontal spacing between items; originally 1.5
    legend.key.width = unit(0.55, "cm"),     # Width of legend keys
    legend.text = element_text(size = 15, margin = margin(l = 8, r = 16, unit = "pt")), # size was 11, l = 4, r = 8
    legend.title = element_text(size = 14, face = "bold", #12
                                margin = margin(r = 15, unit = "pt"))
  ) +
  guides(
    fill = guide_legend(nrow = 1, title = "Sample"),
    color = guide_legend(nrow = 1, title = "Sample"))


# NOTE: ggpubr::get_legend() reaches into the built gtable for the guide box.
# If a future ggplot2 changes the guide grob names again, this is the line
# that will break first.
div_legend_spaced <- ggpubr::get_legend(p_legend_source)

# combine plots:
p_div3_plus_pielou_repel <- grid.arrange(arrangeGrob(p_div_3_boxplots_renyi_repel+
                                                       theme(legend.position = "none",
                                                             axis.title.x = element_blank()),
                                                     p_pielou_repel,
                                                     nrow = 1, widths = c(3,1.1)),
                                         div_legend_spaced, nrow = 2, heights = c(10,1))

# ggsave() confuses the legends here (the object is a gtable, not a ggplot),
# so the png device is used instead.
png(filename = file.path(fig_dir, "Fig5_alpha_diversity_metrics.png"),
    width = 16, height = 8, res = 300, units = "in")
print(p_div3_plus_pielou_repel)
dev.off()

# ------------------------------------------------------------------------ #
# 6. Fig. S5 -- model coefficients across alpha diversity metrics
#
# Sequencing depth varies by more than an order of magnitude across the 13
# samples (see the range printed in stage 1), so every treatment effect is
# estimated alongside a depth term. Three models feed the figure:
#   richness  -- depth x treatment interaction
#   Shannon   -- depth x treatment interaction
#   Simpson   -- treatment only (it carries no depth signal)
# Base-Water is the reference level in all of them.
# ------------------------------------------------------------------------ #

# first layer, sample depth model:
richness_ancova <- lm(qD ~ sample_depth + Label, data = richness_data)

# interaction models:
# sample depth interaction model:
# Check interaction of Label and sample:depth:
richness_int <- lm(qD ~ sample_depth * Label, data = richness_data)
summary(richness_int)

# Run interaction models for all indices
shannon_int <- lm(qD ~ sample_depth * Label, data = shannon_data)
simpson_int <- lm(qD ~ sample_depth * Label, data = simpson_data)
pielou_int <- lm(qD ~ sample_depth * Label, data = pielou_data)

# Since Simpson and Pielou don't depend on sequencing depth,
# let's create depth-free models:
simpson_anova <- lm(qD ~ Label, data = simpson_data)
pielou_anova <- lm(qD ~ Label, data = pielou_data)

# Pielou ANOVA in itself is not a significant model
# only extract coefficients from significant models:

extract_coefs_selective <- function(model, diversity_name, include_interactions = TRUE) {
  coef_table <- summary(model)$coefficients
  result_list <- list()
  
  # Get the main depth coefficient
  base_depth_coef <- 0
  base_depth_se <- 0
  
  if("sample_depth" %in% rownames(coef_table)) {
    base_depth_coef <- coef_table["sample_depth", "Estimate"] * 1000000
    base_depth_se <- coef_table["sample_depth", "Std. Error"] * 1000000
    
    # Add main depth effect
    depth_df <- data.frame(
      Diversity = diversity_name,
      Term = "Sample depth (per million reads)",
      Estimate = base_depth_coef,
      Std_Error = base_depth_se,
      P_value = coef_table["sample_depth", "Pr(>|t|)"],
      Effect_Type = "Depth effect"
    )
    result_list$depth <- depth_df
  }
  
  # Treatment main effects
  treatment_rows <- grepl("^Label", rownames(coef_table)) & !grepl(":", rownames(coef_table))
  if(any(treatment_rows)) {
    treatment_coefs <- coef_table[treatment_rows, ]
    treatment_df <- data.frame(
      Diversity = diversity_name,
      Term = gsub("Label", "", rownames(treatment_coefs)),
      Estimate = treatment_coefs[, "Estimate"],
      Std_Error = treatment_coefs[, "Std. Error"],
      P_value = treatment_coefs[, "Pr(>|t|)"],
      Effect_Type = "Treatment effect"
    )
    result_list$treatment <- treatment_df
  }
  
  # Interaction effects - SHOW TOTAL SLOPES, not just interaction terms
  if(include_interactions) {
    interaction_rows <- grepl("sample_depth:", rownames(coef_table))
    if(any(interaction_rows)) {
      interaction_coefs <- coef_table[interaction_rows, ]
      
      # Handle NAs
      interaction_coefs[is.na(interaction_coefs[, "Estimate"]), "Estimate"] <- 0
      interaction_coefs[is.na(interaction_coefs[, "Pr(>|t|)"]), "Pr(>|t|)"] <- 1
      
      # Calculate TOTAL slopes (main + interaction)
      total_slopes <- (base_depth_coef/1000000) + interaction_coefs[, "Estimate"]
      
      # For standard errors of sums, use: SE(A+B) = sqrt(SE(A)² + SE(B)²) 
      # (assuming independence, which is reasonable here)
      total_se <- sqrt((base_depth_se/1000000)^2 + interaction_coefs[, "Std. Error"]^2)
      
      interaction_df <- data.frame(
        Diversity = diversity_name,
        Term = paste("Sample depth:", gsub("sample_depth:Label", "", rownames(interaction_coefs))),
        Estimate = total_slopes * 1000000,  # Scale back to per million
        Std_Error = total_se * 1000000,
        P_value = interaction_coefs[, "Pr(>|t|)"],  # Keep interaction p-value
        Effect_Type = "Depth sensitivity"
      )
      result_list$interaction <- interaction_df
    }
  }
  
  return(do.call(rbind, result_list))
}


# Extract coefficients from your models
richness_coef_df <- extract_coefs_selective(richness_int, "Species Richness", TRUE)
shannon_coef_df <- extract_coefs_selective(shannon_int, "Shannon Diversity", TRUE) 
simpson_coef_df <- extract_coefs_selective(simpson_anova, "inverse Simpson Diversity", FALSE)  # No interactions

# Combine all coefficients
all_coef_df <- rbind(richness_coef_df, shannon_coef_df, simpson_coef_df)

# Set factor levels for consistent ordering
all_coef_df$Diversity <- factor(all_coef_df$Diversity, 
                                levels = c("Species Richness", "Shannon Diversity", "inverse Simpson Diversity"))
all_coef_df$Effect_Type <- factor(all_coef_df$Effect_Type,
                                  levels = c("Depth effect", "Treatment effect", "Depth sensitivity"))
all_coef_df$Term <- factor(all_coef_df$Term,
                           levels = rev(unique(all_coef_df$Term)))

# Export the final data frame that goes into the plot, for inspection.
write.csv(all_coef_df,
          file = file.path(fig_dir,
                           "FigS5_model_coefficients_across_alpha_diversity_metrics_plotting_data.csv"),
          row.names = FALSE)

# Create horizontal coefficient plot
p_coef_horizontal <- ggplot(all_coef_df, aes(x = Estimate, y = Term, fill = Effect_Type)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_errorbar(aes(xmin = Estimate - Std_Error, xmax = Estimate + Std_Error),
                width = 0.3, alpha = 0.7) +
  
  # Facet by diversity index
  facet_wrap(~ Diversity, scales = "free_x", nrow = 1) +
  
  # Add significance stars
  geom_text(aes(label = ifelse(P_value < 0.05, "*", ""),
                x = ifelse(Estimate > 0, Estimate + Std_Error*1.2, Estimate - Std_Error*1.2)),
            hjust = ifelse(all_coef_df$Estimate > 0, -0.2, 1.2),
            size = 5) +
  
  # Colors for effect types
  scale_fill_manual(values = c(
    "Depth effect" = pals::stepped(n=24)[20],
    "Treatment effect" = pals::stepped(n=24)[7], 
    "Depth sensitivity" = pals::stepped(n=24)[14]
  )) +
  
  # Reference line at zero
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.7) +
  
  # Labels and theme
  labs(
    # title = "Model Coefficients Across Diversity Indices",
    subtitle = "* p < 0.05",
    x = "Coefficient Estimate",
    y = "",
    fill = "Effect Type"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold", size = 11)
  )

# Save the plot
ggsave(file.path(fig_dir, "FigS5_model_coefficients_across_alpha_diversity_metrics.png"),
       p_coef_horizontal, width = 14, height = 7, dpi = 300)

# ------------------------------------------------------------------------ #
# 7. Fig. S6 -- standardized treatment contrasts (Cohen's d)
#
# Four contrasts x four alpha diversity metrics. Base-Water is excluded from
# the contrasts: it has a single observation, so its SD is 0 and every
# comparison against it is inflated. Tiles whose 95% CI excludes zero carry
# a black dot in the top-right corner.
# ------------------------------------------------------------------------ #

# Function to calculate Cohen's d with proper labels
calculate_cohens_d_with_labels <- function(group1_data, group2_data, comparison_name) {
  result <- cohens_d(group1_data, group2_data)
  return(data.frame(
    Comparison = comparison_name,
    Cohens_D = result$Cohens_d,
    CI_low = result$CI_low,
    CI_high = result$CI_high,
    stringsAsFactors = FALSE
  ))
}

# Extract data for the 4 key diversity metrics
richness_data_es <- richness_iNEXT$qD
shannon_data_es <- shannon_iNEXT$qD  
simpson_data_es <- simpson_iNEXT$qD
pielou_data_es <- pielou_iNEXT$qD

# Create label vectors for grouping
labels_es <- richness_iNEXT$Label  # Should be same for all metrics

# Split data by treatment
biofilm_no_idx <- which(labels_es == "Biofilm-no-nutri")
biofilm_w_idx <- which(labels_es == "Biofilm-w/-nutri")
bulk_no_idx <- which(labels_es == "Bulk-no-nutri")
bulk_w_idx <- which(labels_es == "Bulk-w/-nutri")
# base_water_idx <- which(labels_es == "Base-Water")

# Function to calculate all effect sizes for one diversity metric
# exclude base-water from comparisons, because SD = 0, inflates metrics
calculate_all_effects <- function(metric_data, metric_name) {
  
  # Extract groups
  biofilm_no <- metric_data[biofilm_no_idx]
  biofilm_w <- metric_data[biofilm_w_idx]
  bulk_no <- metric_data[bulk_no_idx]
  bulk_w <- metric_data[bulk_w_idx]
  # base_water <- metric_data[base_water_idx]
  
  # Calculate 4 comparisons, 
  effects_list <- list(
    calculate_cohens_d_with_labels(biofilm_w, biofilm_no, "Biofilm: w/-nutri vs no-nutri"),
    calculate_cohens_d_with_labels(bulk_w, bulk_no, "Bulk: w/-nutri vs no-nutri"),
    calculate_cohens_d_with_labels(biofilm_no, bulk_no, "no-nutri: Biofilm vs Bulk"),
    calculate_cohens_d_with_labels(biofilm_w, bulk_w, "w/-nutri: Biofilm vs Bulk")
    # calculate_cohens_d_with_labels(biofilm_no, base_water, "Biofilm-no vs Base-Water"),
    # calculate_cohens_d_with_labels(biofilm_w, base_water, "Biofilm+nutri vs Base-Water"),
    # calculate_cohens_d_with_labels(bulk_no, base_water, "Bulk-no vs Base-Water"),
    # calculate_cohens_d_with_labels(bulk_w, base_water, "Bulk+nutri vs Base-Water")
  )
  
  # Combine and add metric name
  effects_df <- do.call(rbind, effects_list)
  effects_df$Diversity_Metric <- metric_name
  return(effects_df)
}

# Calculate effect sizes for all 4 metrics
richness_effects <- calculate_all_effects(richness_data_es, "Species Richness")
shannon_effects <- calculate_all_effects(shannon_data_es, "Shannon Diversity")
simpson_effects <- calculate_all_effects(simpson_data_es, "Simpson Diversity") 
pielou_effects <- calculate_all_effects(pielou_data_es, "Pielou's Evenness")

# Combine all effect sizes
all_effects_4metrics <- rbind(richness_effects, shannon_effects, 
                              simpson_effects, pielou_effects)

# Set factor levels for proper ordering
all_effects_4metrics$Diversity_Metric <- factor(
  all_effects_4metrics$Diversity_Metric,
  levels = c("Species Richness", "Shannon Diversity", "Simpson Diversity", "Pielou's Evenness")
)

all_effects_4metrics$Comparison <- factor(
  all_effects_4metrics$Comparison,
  levels = c("Biofilm: w/-nutri vs no-nutri", "Bulk: w/-nutri vs no-nutri",
             "no-nutri: Biofilm vs Bulk", "w/-nutri: Biofilm vs Bulk"))
             # "Biofilm-no vs Base-Water", "Biofilm+nutri vs Base-Water",
             # "Bulk-no vs Base-Water", "Bulk+nutri vs Base-Water")
# )

# create significance flag for data points where CI spans over 0:
all_effects_4metrics_sig <- all_effects_4metrics %>%
  mutate(
    # TRUE if entire 95% CI is on one side of 0 (no overlap with 0)
    ci_excludes_zero = (CI_low > 0 & !is.na(CI_low)) | (CI_high < 0 & !is.na(CI_high))
  )

# Export the final data frame that goes into the plot, for inspection
write.csv(all_effects_4metrics_sig,
          file = file.path(fig_dir,
                           "FigS6_standardized_treatment_contrasts_cohens_d_plotting_data.csv"),
          row.names = FALSE)

p_heatmap_4metrics_corner <-
  ggplot(all_effects_4metrics_sig,
         aes(x = Diversity_Metric, y = Comparison, fill = Cohens_D)) +
  geom_tile(color = pals::stepped(24)[8], linewidth = 0.5) +
  # -- dot nudged to top-right corner of each significant tile --
  geom_point(
    data = subset(all_effects_4metrics_sig, ci_excludes_zero),
    aes(x = Diversity_Metric, y = Comparison),
    inherit.aes = FALSE,
    position = position_nudge(x = 0.35, y = -0.35),
    shape = 21, size = 1.8, stroke = 0.15, colour = "black", fill = "black"
  ) +
  geom_text(aes(label = sprintf("%.2f", Cohens_D)),
            color = "black", fontface = "bold", size = 3) +
  scale_fill_gradient2(
    low = pals::stepped(24)[13], mid = "white", high = pals::stepped(24)[1],
    midpoint = 0, breaks = seq(-5, 1), limits = c(-6, 1.02),
    name = "Cohen's d"
  ) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    axis.text.y  = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 10),
    panel.grid   = element_blank(),
    plot.margin  = margin(10, 10, 10, 10)
  )


ggsave(file.path(fig_dir, "FigS6_standardized_treatment_contrasts_cohens_d.png"),
       p_heatmap_4metrics_corner, width = 6, height = 4, dpi = 300)

# ------------------------------------------------------------------------ #
# 8. Fig. S7 -- alpha diversity profile (observed mean +/- SE)
#
# Hill diversity against the order q, averaged over the replicates of each
# sample type, log10 y axis. The q grid the HPC sweep covered has 35 values
# (dense between 0 and 2, coarser out to 8), so the connecting line between
# observed means already traces the profile without any smoothing. This is
# deliberately the observed version: no smoothing parameter to defend.
# ------------------------------------------------------------------------ #

# Calculate means by q value and treatment
div_profiles_summary <- div_profiles %>%
  group_by(Order.q, Label) %>%
  summarise(
    qD_mean = mean(qD, na.rm = TRUE),
    qD_se = sd(qD, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
)

# Export the final data frame that goes into the plot, for inspection.
write.csv(div_profiles_summary,
          file = file.path(fig_dir, "FigS7_alpha_diversity_profile_plotting_data.csv"),
          row.names = FALSE)

# Plot with manual means
p_diversity_profile_mean <- ggplot(div_profiles_summary,
                                   aes(x = Order.q, y = qD_mean, color = Label)) +
  # Individual sample points (more transparent)
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = qD_mean - qD_se,
                    ymax = qD_mean + qD_se), width = 0.1)+
  scale_color_manual(values = stepped_cols) +
  scale_fill_manual(values = stepped_cols) +
  scale_y_log10(labels = scales::comma_format()) +
  scale_x_continuous(
    breaks = c(0:8)
  ) +
  labs(
    x = "Hill Number Parameter (q)",
    y = "Hill's Diversity"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14)
  )

# Save the plot
ggsave(file.path(fig_dir, "FigS7_alpha_diversity_profile.png"),
       p_diversity_profile_mean,
       width = 10, height = 8, dpi = 300)

# ------------------------------------------------------------------------ #
# 9. Fig. S8 -- derivative of the alpha diversity profile
#
# The derivative is taken off a LOESS fit of qD against q, per sample type,
# span = 0.2, evaluated on the grid seq(0, 8, by = 0.1). Both numbers belong
# in Fig. S8's caption -- they are the only smoothing choices in the whole
# analysis.
#
# The LOESS fit exists ONLY to give the derivative a smooth, evenly spaced
# curve to difference. It is not plotted on its own: at span 0.2 over a
# 35-point q grid it simply retraces the observed means of Fig. S7, so the
# fitted variant of that figure was dropped.
# ------------------------------------------------------------------------ #

# Calculate first derivatives ####
# Step 1: Fit LOESS models per Label (or per Label + batch for replicates)
# Group by Label (condition) for means; adjust to ~Label + batch if per-batch fits needed
loess_fits <- div_profiles %>%
  group_by(Label) %>%  # Or group_by(Label, batch) for batch-specific
  do(model = loess(qD ~ Order.q, data = ., span = 0.2))  # Match your span

# Step 2: Predict on a fine grid of q (for smooth derivatives)
q_grid <- seq(0, 8, by = 0.1)  # Dense grid from 0 to 8
predictions <- loess_fits %>%
  rowwise() %>%
  mutate(pred_df = list(
    tibble(
      Order.q = q_grid,
      qD_pred = predict(model, newdata = tibble(Order.q = q_grid))
    )
  )) %>%
  unnest(pred_df)

# Step 3: Calculate numerical first derivative (central difference for accuracy)
# For edge points, use forward/backward; central elsewhere
predictions <- predictions %>%
  group_by(Label) %>%  # Or group_by(Label, batch)
  arrange(Order.q) %>%
  mutate(
    # Forward diff for first point
    derivative = ifelse(row_number() == 1, 
                        (lead(qD_pred) - qD_pred) / (lead(Order.q) - Order.q),
                        NA),
    # Backward diff for last point
    derivative = ifelse(row_number() == n(), 
                        (qD_pred - lag(qD_pred)) / (Order.q - lag(Order.q)),
                        derivative),
    # Central diff for middle points
    derivative = ifelse(is.na(derivative), 
                        (lead(qD_pred) - lag(qD_pred)) / (lead(Order.q) - lag(Order.q)),
                        derivative)
  ) %>%
  ungroup()

# Higher absolute richness in bulk-no-nutri produces scale artifact:
# lower overall diversity means smaller absolute changes, even if the community is more uneven
# absolute derivatives are confounded by the magnitude of qD ->
# high-richness systems inherently have more "room to fall," masking true unevenness.
# we need relative derivatives (normalize derivatives by qD_pred:
predictions$rel_derivative <- predictions$derivative/predictions$qD_pred

# Export the final data frame that goes into the plot, for inspection.
# `predictions` still carries the loess model objects in a list column, which
# write.csv() cannot encode, so only the plotted columns are written out.
write.csv(as.data.frame(predictions)[, c("Label", "Order.q", "qD_pred",
                                         "derivative", "rel_derivative")],
          file = file.path(fig_dir,
                           "FigS8_alpha_diversity_profile_derivative_plotting_data.csv"),
          row.names = FALSE)

# relative derivatives
p_rel_derivative <- ggplot(predictions, aes(x = Order.q, y = rel_derivative, color = Label)) +
  geom_line(linewidth = 1) +
  # geom_point()+
  scale_color_manual(values = stepped_cols) +
  scale_x_continuous(breaks = c(0:8), limits = c(0,8)) +
  labs(
    title = "Normalized First Derivative of Diversity Profiles",
    x = "Hill Number Parameter (q)",
    y = expression(Delta ~~ "qD" / ("qD" %*% Delta ~~ "q"))) +
    # y = expression("Normalized" ~ frac(d, dq) ~ "qD"))+
    # y = expression(d ~ "log(qD)" / d ~ "q"))+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave(file.path(fig_dir, "FigS8_alpha_diversity_profile_derivative.png"),
       p_rel_derivative,
       width = 10, height = 8, dpi = 300)

rel_derivatives_summary <- predictions %>%
  filter(Order.q >= 0 & Order.q <= 2) %>%  # Example: Focus on q=0-2
  group_by(Label) %>%
  summarize(mean_deriv = mean(rel_derivative, na.rm = TRUE),
            sd_deriv = sd(rel_derivative, na.rm = TRUE))
print(rel_derivatives_summary)


# ######################################################################## #
#                                                                          #
#  10. ADDITIONAL STATISTICS -- NOT MANUSCRIPT FIGURES                     #
#                                                                          #
#  Nothing below is needed to reproduce Fig. 5 or Figs. S5-S8. These are    #
#  the supporting analyses that justify the modelling choices made above,   #
#  kept because they answer the obvious reviewer questions. They can be     #
#  skipped without affecting any figure.                                    #
#                                                                          #
# ######################################################################## #

# ------------------------------------------------------------------------ #
# 10a. Does the Kraken2 classification rate add explanatory power?
#
# It does not -- which is why the classification_rate column is carried
# through section 4 but never enters the models behind Fig. S5.
# Needs data/kraken2/kraken2_classification_rates.csv (see section 4).
# ------------------------------------------------------------------------ #

# extended ANCOVA model:
richness_ancova_full <- lm(qD ~ sample_depth + classification_rate + Label, 
                           data = richness_data)

# Test if classification rate adds explanatory power
anova(richness_ancova, richness_ancova_full)

summary(richness_ancova)
summary(richness_ancova_full)

# Decompose the variance explained by different components
richness_depth_only <- lm(qD ~ sample_depth, data = richness_data)
richness_label_only <- lm(qD ~ Label, data = richness_data)

# Compare R² values
cat("Depth only R²:", summary(richness_depth_only)$r.squared, "\n")
cat("Label only R²:", summary(richness_label_only)$r.squared, "\n") 
cat("Full model R²:", summary(richness_ancova)$r.squared, "\n")
cat("Full interaction model R²:", summary(richness_int)$r.squared, "\n")

# ------------------------------------------------------------------------ #
# 10b. How big is the sequencing depth effect, in species?
# ------------------------------------------------------------------------ #

# Look at the actual magnitude of the depth effect
depth_coef <- coef(richness_ancova)["sample_depth"]
depth_range <- range(richness_data$sample_depth)
depth_effect_range <- depth_coef * diff(depth_range)

cat("Sequencing depth coefficient:", depth_coef, "\n")
cat("Range of sequencing depths:", depth_range, "\n") 
# The number of additional species gained with across the increasing depth range:
cat("Total effect across depth range:", depth_effect_range, "species\n")

# Compare to biological differences
label_effects <- coef(richness_ancova)[grepl("Label", names(coef(richness_ancova)))]
cat("Biological effect magnitudes:", label_effects, "\n")

# Calculate unique contributions
anova(richness_ancova)  # This shows the variance explained by each term

# ------------------------------------------------------------------------ #
# 10c. Tukey HSD over the sample types
#
# With 13 samples the interaction model leaves only 4 error degrees of
# freedom, so its post-hoc test has almost no power. Both are run and
# compared; the additive model is the one worth reporting.
# ------------------------------------------------------------------------ #

# Sample comparisons - Tukey test
# since the sample_depth:Label interaction test had the best results
# we use that model:
library(multcomp)
richness_int_aov <- aov(richness_int)
richness_int_tukey <- glht(richness_int_aov, linfct = mcp("Label" = "Tukey"))
summary(richness_int_tukey)

library(broom)
tukey_tidy <- tidy(richness_int_tukey)
write.table(tukey_tidy,
            file = file.path(fig_dir, "Tukey_HSD_richness_interaction_model.tsv"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)


# With only 3-4 observations per group and an interaction model with many parameters,
# you have very little power to detect differences
# Degrees of freedom for error: only 4!

cat("Model degrees of freedom:\n")
cat("Total observations:", nrow(richness_data), "\n")
cat("Model parameters:", length(coef(richness_int)[!is.na(coef(richness_int))]), "\n")
cat("Error df:", summary(richness_int)$df[2], "\n")

# Compare power between models:
cat("Additive model error df:", summary(richness_ancova)$df[2], "\n")
cat("Interaction model error df:", summary(richness_int)$df[2], "\n")

# The interaction model has only 4 error df vs 7 for additive model
# This drastically reduces power for post-hoc tests

# Check Tukey test with additive model:
richness_aov <- aov(richness_ancova)
richness_tukey <- glht(richness_aov, linfct = mcp("Label" = "Tukey"))
summary(richness_tukey)


# ------------------------------------------------------------------------ #
# 10d. Cohen's d across the whole Hill number spectrum
#
# The same four contrasts as Fig. S6, but computed at every q value of the
# sweep instead of only at q = 0, 1, 2. It is the companion of Fig. S6 and
# Fig. S7: it shows at which order of diversity each treatment contrast is
# strongest.
# ------------------------------------------------------------------------ #

# Function to calculate effect sizes for diversity profile data
calculate_profile_effects <- function(div_profiles_data) {
  
  # Get unique q values
  q_values <- unique(div_profiles_data$Order.q)
  
  # Initialize results list
  all_profile_effects <- list()
  
  # Loop through each q value
  for (q_val in q_values) {
    
    # Filter data for this q value
    q_data <- div_profiles_data[div_profiles_data$Order.q == q_val, ]
    
    # Extract groups
    biofilm_no <- q_data$qD[q_data$Label == "Biofilm-no-nutri"]
    biofilm_w <- q_data$qD[q_data$Label == "Biofilm-w/-nutri"]
    bulk_no <- q_data$qD[q_data$Label == "Bulk-no-nutri"]
    bulk_w <- q_data$qD[q_data$Label == "Bulk-w/-nutri"]
    # base_water <- q_data$qD[q_data$Label == "Base-Water"]
    
    # Skip if any group is empty
    if (length(biofilm_no) == 0 || length(biofilm_w) == 0 || 
        length(bulk_no) == 0 || length(bulk_w) == 0 ) {
      next
    }
    
    # Calculate all 8 comparisons for this q value
    q_effects <- list(
      calculate_cohens_d_with_labels(biofilm_w, biofilm_no, "Biofilm: +Nutri vs -Nutri"),
      calculate_cohens_d_with_labels(bulk_w, bulk_no, "Bulk: +Nutri vs -Nutri"),
      calculate_cohens_d_with_labels(biofilm_no, bulk_no, "No Nutri: Biofilm vs Bulk"),
      calculate_cohens_d_with_labels(biofilm_w, bulk_w, "+Nutri: Biofilm vs Bulk")
      # calculate_cohens_d_with_labels(biofilm_no, base_water, "Biofilm-no vs Base-Water"),
      # calculate_cohens_d_with_labels(biofilm_w, base_water, "Biofilm+nutri vs Base-Water"),
      # calculate_cohens_d_with_labels(bulk_no, base_water, "Bulk-no vs Base-Water"),
      # calculate_cohens_d_with_labels(bulk_w, base_water, "Bulk+nutri vs Base-Water")
    )
    
    # Combine and add q value
    q_effects_df <- do.call(rbind, q_effects)
    q_effects_df$q_value <- q_val
    
    all_profile_effects[[as.character(q_val)]] <- q_effects_df
  }
  
  # Combine all results
  return(do.call(rbind, all_profile_effects))
}

# Calculate effect sizes for all q values
profile_effects <- calculate_profile_effects(div_profiles)

# Set factor levels for proper ordering
profile_effects$Comparison <- factor(
  profile_effects$Comparison,
  levels = c("Biofilm: +Nutri vs -Nutri", "Bulk: +Nutri vs -Nutri",
             "No Nutri: Biofilm vs Bulk", "+Nutri: Biofilm vs Bulk")
             # "Biofilm-no vs Base-Water", "Biofilm+nutri vs Base-Water",
             # "Bulk-no vs Base-Water", "Bulk+nutri vs Base-Water")
)



# heatmap 2 plotting
p_heatmap_profile <- ggplot(profile_effects, 
                            aes(x = q_value, y = Comparison, fill = Cohens_D)) +
  geom_tile(color = pals::stepped(24)[8], linewidth = 0.1) +
  
  # Add text only for specific q values to avoid overcrowding
  geom_text(data = subset(profile_effects, q_value %in% c(0, 1, 2, 4, 6, 8)),
            aes(label = sprintf("%.1f", Cohens_D)), 
            color = "black", fontface = "bold", size = 3) +
  
  # Color scale: same as first heatmap for consistency
  scale_fill_gradient2(
    low = pals::stepped(24)[13], mid = "white", high = pals::stepped(24)[1],
    midpoint = 0,
    breaks = seq(-7, 1),
    limits = c(-7.01, 1.02),
    name = "Cohen's d"
  ) +
  
  # X-axis with key q values highlighted
  scale_x_continuous(
    breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8),
    minor_breaks = seq(0, 8, by = 0.5)
  ) +
  
  # Labels and theme
  labs(
    title = "Effect Sizes Across Hill Number Spectrum (q = 0-8)",
    subtitle = "Cohen's d for pairwise comparisons across diversity orders",
    x = "Hill Number Parameter (q)",
    y = "Comparison"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_line(color = "gray95", size = 0.25),
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )

write.csv(profile_effects,
          file = file.path(fig_dir, "effect_sizes_across_hill_spectrum_plotting_data.csv"),
          row.names = FALSE)

# Save heatmap 2
ggsave(file.path(fig_dir, "effect_sizes_across_hill_spectrum.png"),
       p_heatmap_profile, width = 14, height = 8, dpi = 300)

cat("\nSaved Fig5 and FigS5-S8 into", fig_dir, "\n")
