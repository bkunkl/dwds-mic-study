# ======================================================================== #
# Fig. 7 and Fig. S9 -- taxonomic composition of the 13 metagenomes
#
# Both figures start from the same Kraken2/Bracken phyloseq object
# (confidence 0.1, species level) restricted to Bacteria + Archaea, but are
# otherwise independent analyses:
#
#   PART I  (Fig. S9) -- per sample type, the genera making up the first 80%
#           of the community, as mean +/- SE over the replicates of that type,
#           with the individual replicates overplotted; plus a panel counting
#           how many species sit in that top-80% group versus in the
#           remaining ~20% ("Other").
#
#   PART II (Fig. 7)  -- CLR PCA of the species table (panel a, with 90%
#           normal-theory ellipses per nutrient condition) paired with a
#           microViz "iris" plot of the 15 most abundant species (panel b).
#
# Experiments E1/E2/E3 (the `batch` column of the metadata) are shown as
# #1/#2/#3 in the figures. Only E3 has the pre-incubation Base-Water
# control, so that sample type has a single replicate (no error bars, no jittered points).
#
# NOTE ON THE ORIGINAL SCRIPT: it looped over three Kraken2 confidence
# levels (0.0, 0.1, 0.2) and carried every phyloseq object in a named list.
# Only confidence 0.1 is used in the manuscript, so the loop is gone and the
# helper functions that took a list of phyloseq objects now take one object.
# Their bodies are unchanged.
#
# INPUT
#   results/kraken_bracken/brackenSpecies_conf01.rds
#       phyloseq object built by scripts/local/01_build_phyloseq_object.sh
#
# OUTPUT (into results/figures/)
#   FigS9_top80_pct_genera_distribution.png
#   FigS9_top80_pct_genera_distribution_plotting_data.csv
#   FigS9_top80_pct_genera_individual_replicates.csv
#   FigS9_species_count_80_20_ratio_plotting_data.csv
#   Fig7_PCA_and_iris_paired_plot_top_species.png
#   Fig7_PCA_and_iris_paired_plot_top_species_plotting_data.csv
#   Fig7_iris_panel_bar_data.csv
#   plus the single panels of both figures, saved separately
#
# Run from the repository root.
# ======================================================================== #

library(tidyverse)
library(phyloseq)
library(cowplot)
library(microViz)
library(ggthemes)
library(pals)
library(ggrepel)
library(grid)

# repo-relative paths -- run this script with the repo root as working directory
ps_file <- "results/kraken_bracken/brackenSpecies_conf01.rds"

fig_dir <- "results/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------ #
# 0. The phyloseq object, shared by both figures
# ------------------------------------------------------------------------ #

conf01 <- readRDS(ps_file)
ps_prok <- subset_taxa(conf01, Domain == "Bacteria" | Domain == "Archaea")

# ======================================================================== #
# PART I -- Fig. S9: top-80% genera per sample type
# ======================================================================== #

# ------------------------------------------------------------------------ #
# 1. Top genera per sample type, with SE across the replicates
# ------------------------------------------------------------------------ #

# Calculate statistics WITHIN each sample type, then aggregate
get_top_taxa_with_stats_by_label <- function(ps_obj, tax_level = "Genus", 
                                             cumulative_threshold = 0.8, 
                                             min_taxa = 4) {
  
  # Normalize the data first
  ps_norm <- transform_sample_counts(ps_obj, function(z) z/sum(z))
  
  # Convert to data frame for easier manipulation
  otu_norm_df <- as.data.frame(as(otu_table(ps_norm), "matrix"))
  tax_df <- as.data.frame(tax_table(ps_norm))
  meta_df <- as(sample_data(ps_norm), "data.frame")
  
  # Add taxonomy info
  otu_norm_df$Taxa <- rownames(otu_norm_df)
  otu_norm_df <- merge(otu_norm_df, tax_df, by.x = "Taxa", by.y = "row.names")
  
  # Get sample names (excluding Taxa and taxonomy columns)
  sample_cols <- setdiff(names(otu_norm_df), c("Taxa", colnames(tax_df)))
  
  results <- list()
  
  # Process each Label (sample type)
  for(label in unique(meta_df$Label)) {
    
    # Get samples belonging to this label
    samples_for_label <- rownames(meta_df[meta_df$Label == label, ])
    samples_for_label <- intersect(samples_for_label, sample_cols)
    
    if(length(samples_for_label) == 0) next
    
    # IMPORTANT: For each genus, we have abundance values in each replicate
    # Create a list to store all abundance values per genus
    genus_abundances <- list()
    genus_species_count <- list()  # NEW: Track species count per genus
    
    # Iterate through all taxa
    for(i in seq_len(nrow(otu_norm_df))) {
      genus_name <- otu_norm_df[[tax_level]][i]
      
      # Skip NA taxa
      if(is.na(genus_name) || genus_name == "") next
      
      # Get abundances for this taxon in all replicates of this label
      abundances <- as.numeric(otu_norm_df[i, samples_for_label])
      
      # Add to the list (or sum if genus already exists)
      if(genus_name %in% names(genus_abundances)) {
        genus_abundances[[genus_name]] <- genus_abundances[[genus_name]] + abundances
        genus_species_count[[genus_name]] <- genus_species_count[[genus_name]] + 1
      } else {
        genus_abundances[[genus_name]] <- abundances
        genus_species_count[[genus_name]] <- 1  # later added: initialize species count
      }
    }
    
    # Now calculate statistics for each genus within this label
    agg_df_stats <- data.frame(
      Taxon = names(genus_abundances),
      Mean = sapply(genus_abundances, mean),
      SD = sapply(genus_abundances, sd),
      Species_count = sapply(genus_species_count, function(x) x),  # NEW: Add species count
      stringsAsFactors = FALSE
    )
    
    # Calculate SE
    n_samples <- length(samples_for_label)
    agg_df_stats$SE <- agg_df_stats$SD / sqrt(n_samples)
    agg_df_stats$N <- n_samples
    
    # Remove NA taxa
    agg_df_stats <- agg_df_stats[!is.na(agg_df_stats$Taxon) & agg_df_stats$Taxon != "", ]
    
    # Sort by mean abundance descending
    agg_df_stats <- agg_df_stats[order(-agg_df_stats$Mean), ]
    
    # Calculate cumulative abundance
    agg_df_stats$Cumulative <- cumsum(agg_df_stats$Mean)
    
    # Find how many taxa needed to reach threshold
    n_for_threshold <- which(agg_df_stats$Cumulative >= cumulative_threshold)[1]
    if(is.na(n_for_threshold)) n_for_threshold <- nrow(agg_df_stats)
    
    # Take at least min_taxa, but add more if needed for threshold
    n_taxa <- max(min_taxa, n_for_threshold)
    n_taxa <- min(n_taxa, nrow(agg_df_stats))  # Don't exceed available taxa
    
    # Get top taxa
    top_taxa <- agg_df_stats[1:n_taxa, ]
    
    # Calculate "Other" category including species count
    if(n_taxa < nrow(agg_df_stats)) {
      other_mean <- sum(agg_df_stats$Mean[(n_taxa + 1):nrow(agg_df_stats)])
      other_sd <- NA
      other_se <- NA
      other_species_count <- sum(agg_df_stats$Species_count[(n_taxa + 1):nrow(agg_df_stats)])  # NEW
    } else {
      other_mean <- 0
      other_sd <- NA
      other_se <- NA
      other_species_count <- 0  # NEW
    }
    
    results[[label]] <- list(
      top_taxa = top_taxa,
      other_mean = other_mean,
      other_sd = other_sd,
      other_se = other_se,
      other_species_count = other_species_count,  # NEW
      total_abundance_covered = sum(top_taxa$Mean),
      n_taxa_selected = n_taxa,
      n_samples = n_samples
    )
  }
  
  return(results)
}

# Get top genera for each sample type WITH statistics calculated within each type
top_genera_with_se <- get_top_taxa_with_stats_by_label(ps_prok, 
                                                       tax_level = "Genus", 
                                                       cumulative_threshold = 0.8, 
                                                       min_taxa = 4)

# Print summary for each sample type
for(label in names(top_genera_with_se)) {
  cat("\n=== Sample Type:", label, "===\n")
  result <- top_genera_with_se[[label]]
  cat("Number of samples:", result$n_samples, "\n")
  cat("Number of genera selected:", result$n_taxa_selected, "\n")
  cat("Total abundance covered:", round(result$total_abundance_covered * 100, 1), "%\n")
  cat("Other abundance:", round(result$other_mean * 100, 1), "%\n")
  
  # Show top genera with stats
  cat("\nTop genera:\n")
  for(i in 1:nrow(result$top_taxa)) {
    genus <- result$top_taxa$Taxon[i]
    mean_ab <- round(result$top_taxa$Mean[i] * 100, 2)
    se <- round(result$top_taxa$SE[i] * 100, 2)
    sd <- round(result$top_taxa$SD[i] * 100, 2)
    cat(sprintf("%2d. %s: %.2f%% ± %.2f%% (SD: %.2f%%)\n", i, genus, mean_ab, se, sd))
  }
}

# ------------------------------------------------------------------------ #
# 2. The same abundances per individual replicate (the jittered points)
# ------------------------------------------------------------------------ #

# Create individual genus-abundance dataframe
# Aggregate counts to genus level
ps_genus <- tax_glom(
  ps_prok,
  taxrank = "Genus",
  NArm = FALSE
)

# Convert each sample to relative abundance
ps_genus_rel <- transform_sample_counts(
  ps_genus,
  function(x) {
    if (sum(x) == 0) {
      x
    } else {
      x / sum(x)
    }
  }
)

# Convert phyloseq object to a long dataframe
genus_individual_raw <- psmelt(ps_genus_rel) %>%
  transmute(
    Sample_ID = as.character(Sample),
    Label = as.character(Label),
    Genus = as.character(Genus),
    Abundance = Abundance
  ) %>%
  mutate(
    # Keep unclassified genera because they belong in "Other"
    Genus = ifelse(
      is.na(Genus) | Genus == "",
      "Unclassified",
      Genus
    )
  ) %>%
  # Safety aggregation in case multiple taxa have the same genus name
  group_by(Label, Sample_ID, Genus) %>%
  summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  )


# Convert nonselected genera to sample-specific "Other"
# The category "Other" must be calculated separately for each label and sample
create_individual_genus_dataframe <- function(
    genus_abundance_data,
    top_genera_results
) {
  
  individual_list <- list()
  
  for (label_name in names(top_genera_results)) {
    
    selected_genera <- top_genera_results[[label_name]]$top_taxa$Taxon
    
    label_individual <- genus_abundance_data %>%
      filter(Label == label_name)
    
    # Individual abundances of the selected genera
    selected_data <- label_individual %>%
      filter(Genus %in% selected_genera) %>%
      select(Label, Sample_ID, Genus, Abundance)
    
    # Ensure that absent selected genera are represented by zero
    selected_data <- selected_data %>%
      tidyr::complete(
        Sample_ID = unique(label_individual$Sample_ID),
        Genus = selected_genera,
        fill = list(Abundance = 0)
      ) %>%
      mutate(Label = label_name) %>%
      select(Label, Sample_ID, Genus, Abundance)
    
    # Sum all genera outside the selected top genera within each sample
    other_data <- label_individual %>%
      filter(!Genus %in% selected_genera) %>%
      group_by(Label, Sample_ID) %>%
      summarise(
        Abundance = sum(Abundance, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(Genus = "Other") %>%
      select(Label, Sample_ID, Genus, Abundance)
    
    individual_list[[label_name]] <- bind_rows(
      selected_data,
      other_data
    )
  }
  
  bind_rows(individual_list)
}

# Run the function to create the plotting dataframe:
genus_individual_df <- create_individual_genus_dataframe(
  genus_abundance_data = genus_individual_raw,
  top_genera_results = top_genera_with_se
)

# ------------------------------------------------------------------------ #
# 3. Plotting data frame (means, SE, species counts)
# ------------------------------------------------------------------------ #

# Create plotting dataframe with SE
  create_plotting_dataframe_with_se <- function(top_genera_results) {
    plot_df_list <- list()
    
    for(label in names(top_genera_results)) {
      result <- top_genera_results[[label]]
      
      # Add top taxa with species count
      sample_df <- data.frame(
        Label = label,
        Genus = result$top_taxa$Taxon,
        Mean_Abundance = result$top_taxa$Mean,
        SE = result$top_taxa$SE,
        SD = result$top_taxa$SD,
        N = result$top_taxa$N,
        Species_count = result$top_taxa$Species_count,  # NEW
        stringsAsFactors = FALSE
      )
      
      # Add "Other" category with species count
      other_df <- data.frame(
        Label = label,
        Genus = "Other",
        Mean_Abundance = result$other_mean,
        SE = result$other_se,
        SD = result$other_sd,
        N = result$n_samples,
        Species_count = result$other_species_count,  # NEW
        stringsAsFactors = FALSE
      )
      
      # Combine this sample's data
      sample_combined <- rbind(sample_df, other_df)
      plot_df_list[[label]] <- sample_combined
    }
    
    # Combine all samples
    if(length(plot_df_list) > 0) {
      plot_df <- do.call(rbind, plot_df_list)
      rownames(plot_df) <- NULL
      return(plot_df)
    } else {
      stop("No valid sample data could be processed")
    }
  }

# Create plotting dataframe with SE
  genus_plot_df_se <- create_plotting_dataframe_with_se(top_genera_with_se)

# Add genus names with species count
# genus_plot_df_se <- genus_plot_df_se %>%
#   mutate(
#     Genus_with_count = ifelse(
#       Genus != "Other",
#       paste0(Genus, " (", Species_count, ")"),
#       Genus
#     )
#   )
  
  # Add genus names with species count
genus_plot_df_se <- genus_plot_df_se %>%
  mutate(
    Genus_with_count = case_when(
      Genus == "Other" & !is.na(Species_count) & Species_count > 0 ~ 
        paste0("Other (", Species_count, ")"),
      Genus == "Other" ~ 
        "Other",
      TRUE ~ 
        paste0(Genus, " (", Species_count, ")")
    )
  )

# ------------------------------------------------------------------------ #
# 4. Genus colours and the barplot function
# ------------------------------------------------------------------------ #

# Define color scheme for genera
genus_colors <- c(
  "Acidovorax" = "#F1CE63",
  "Aeromonas" = "#DE9D9A",
  "Alicycliphilus" = "#B07AA1",
  "Azonexus" = "#FFBE7D", 
  "Azospira" = "#A0CBE8",
  "Bradyrhizobium" = "#617336",
  "Burkholderia" = "#4E79A7",
  "Caulobacter" = "#710019",
  "Comamonas" = "#59A14F",
  "Dechloromonas" = "#B6992D",
  "Diaphorobacter" = "#86BCB6",
  "Enterobacter" = "#AB6B51",
  "Ferribacterium" = "#F28E2B",
  "Microbacterium" = "#8CD17D",
  "Mycobacterium" = "#E16A59",
  "Pseudomonas" = "#499894",
  "Sphingobium" = "#523A28",
  "Sphingomonas" = "#9D7660",
  "Sphingopyxis" = "#D7B5A6",
  "Thiobacillus" = "#E0A030",
  "Other" = "grey20"
)

# Updated barplot function with error bars
create_genus_barplot_with_se <- function(
    label_name,
    plot_data,
    individual_data,
    genus_colors,
    specific_x_limits = c(
      "Bulk-no-nutri" = 36,
      "Base-Water" = 31)
    ) {
  
  # Filter data for this label
  label_data <- plot_data %>%
    filter(Label == label_name) %>%
    arrange(desc(Mean_Abundance))
  
  # Calculate number of genera for plot height
  n_categories <- nrow(label_data)
  
  # Assign colors
  label_data$color <- ifelse(
    label_data$Genus %in% names(genus_colors),
    genus_colors[label_data$Genus],
    "#BAB0AC"
  )
  
  # Separate "Other" from regular genera for proper ordering
  other_data <- label_data[label_data$Genus == "Other", ]
  regular_data <- label_data[label_data$Genus != "Other", ]
  
  # Create factor levels: "Other" at bottom, then regular genera in descending abundance order
  if(nrow(other_data) > 0) {
    level_order <- c(other_data$Genus_with_count, rev(regular_data$Genus_with_count))
  } else {
    level_order <- rev(regular_data$Genus_with_count)
  }
  
  # Convert Genus to factor with levels ordered by abundance
  label_data$Genus_with_count <- factor(label_data$Genus_with_count, levels = level_order)
  
  # Prepare individual sample-level values for this label
  individual_label_data <- individual_data %>%
    filter(
      Label == label_name,
      Genus %in% label_data$Genus
    ) %>%
    left_join(
      label_data %>%
        select(Genus, Genus_with_count),
      by = "Genus"
    ) %>%
    mutate(
      Genus_with_count = factor(
        Genus_with_count,
        levels = level_order
      )
    )
  
  n_individual_samples <- dplyr::n_distinct(
    individual_label_data$Sample_ID
  )
  
  point_layer <- if (n_individual_samples > 2) {
    geom_point(
      data = individual_label_data,
      aes(
        x = Abundance * 100,
        y = Genus_with_count,
        fill = Genus
      ),
      position = position_jitter(
        width = 0,
        height = 0.2,
        seed = 123
      ),
      shape = 21,
      color = "black",
      stroke = 0.25,
      size = 3,
      alpha = 0.5,
      inherit.aes = FALSE,
      na.rm = TRUE
    )
  } else {
    NULL
  }
  
    # Calculate smart label position with SE handling
  label_data$SE_for_label <- ifelse(is.na(label_data$SE), 0, label_data$SE)
  label_data$label_x <- label_data$Mean_Abundance + label_data$SE_for_label + 0.01
  
  # handle limits, when there is an outlier
  global_x_upper <- max(
    71,
    genus_individual_df$Abundance * 100,
    na.rm = TRUE
  )
  
  # global_x_upper <- ceiling(global_x_upper / 5) * 5
  global_x_upper <- global_x_upper + 2
  
  # Use a label-specific x-axis limit when defined;
  # otherwise use the common global limit
  if (label_name %in% names(specific_x_limits)) {
    x_upper <- unname(specific_x_limits[label_name])
  } else {
    x_upper <- global_x_upper
  }
  
  # Create the plot with error bars
  p <- ggplot(label_data, aes(x = Mean_Abundance * 100, y = Genus_with_count)) +
    # Add error bar geom (only for taxa with SE values)
    geom_errorbar(
      aes(xmin = (Mean_Abundance - SE) * 100, 
          xmax = (Mean_Abundance + SE) * 100),
      orientation = "y",
      width = 0.3,
      color = "black",
      linewidth = 0.5,
      na.rm = TRUE
    ) +
    # Add bars
    geom_col(fill = label_data$color, 
             color = "black", 
             linewidth = 0.2,
             width = 0.8) +
    
    point_layer +
    # Add individual sample-level abundances
    # geom_point(
    #   data = individual_label_data,
    #   aes(
    #     x = Abundance * 100,
    #     y = Genus_with_count,
    #     fill = Genus
    #   ),
    #   position = position_jitter(
    #     width = 0,
    #     height = 0.2,
    #     seed = 123
    #   ),
    #   shape = 21,
    #   color = "black",
    #   stroke = 0.5,
    #   size = 3,
    #   alpha = 0.5,
    #   inherit.aes = FALSE, 
    #   na.rm = TRUE
    # ) +
    
    # Add percentage labels at the end of bars, positioned outside error bars
    geom_text(aes(x = label_x * 100,
                  label = paste0(round(Mean_Abundance * 100, 1), "%")),
              hjust = -0.2,
              size = 3.3,
              fontface = "bold",
              na.rm = TRUE) +  # <- KEY: Add na.rm = TRUE
    
    # ...rest of theme and styling remains the same...
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.8),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5,
                                margin = margin(b = 10)),
      plot.margin = margin(t = 15, r = 25, b = 10, l = 22)
    ) +
    labs(
      title = label_name,
      x = "Relative Abundance (%)",
      y = NULL
    ) +
    scale_fill_manual(
      values = genus_colors,
      na.value = "#BAB0AC",
      guide = "none"
    ) +
    scale_x_continuous(
      limits = c(0, x_upper),
      expand = c(0, 0)
    )
  
  return(list(plot = p, n_categories = n_categories))
}

# Create plots for each sample type (now including Base-Water)
sample_types <- unique(genus_plot_df_se$Label)
# Reorder so Base-Water is last (optional)
sample_types <- c(setdiff(sample_types, "Base-Water"), "Base-Water")

genus_plots_se <- list()

for(label in sample_types) {
  cat("Creating plot for:", label, "\n")
  
  plot_result <- create_genus_barplot_with_se(
    label_name = label,
    plot_data = genus_plot_df_se,
    individual_data = genus_individual_df,
    genus_colors = genus_colors)
  genus_plots_se[[label]] <- plot_result
  
  # Calculate proportional height
  base_height <- 1
  height_per_category <- 0.4
  plot_height <- base_height + (plot_result$n_categories * height_per_category)
  
  # Save individual plot
  filename <- file.path(fig_dir, paste0("genus_barplot_se_", gsub("[^A-Za-z0-9]", "_", label), ".png"))
  
  ggsave(
    filename = filename,
    plot = plot_result$plot,
    width = 8,
    height = plot_height,
    dpi = 300,
    bg = "white"
  )
  
  cat("  Saved:", filename, "- Height:", round(plot_height, 1), "inches\n")
}

# ------------------------------------------------------------------------ #
# 5. Species count in the top-80% group vs. in "Other" (top panel of Fig. S9)
# ------------------------------------------------------------------------ #

ps_prok_merge <- merge_samples(ps_prok, "Label")
if (!taxa_are_rows(ps_prok_merge)) {
  ps_prok_merge <- t(ps_prok_merge)
}
# Normalize AFTER merging (this handles the replicate issue naturally)
ps_prok_norm <- transform_sample_counts(ps_prok_merge, function(z) z/sum(z))
# Create new object to include Base-Water
ps_prok_with_base <- ps_prok_norm

get_top_taxa_by_sample <- function(ps_obj, tax_level = "Genus", cumulative_threshold = 0.8, min_taxa = 5) {
  # Convert to data frame for easier manipulation
  otu_norm_df <- as.data.frame(as(otu_table(ps_obj), "matrix"))
  tax_df <- as.data.frame(tax_table(ps_obj))
  
  # Add taxonomy info
  otu_norm_df$Taxa <- rownames(otu_norm_df)
  otu_norm_df <- merge(otu_norm_df, tax_df, by.x = "Taxa", by.y = "row.names")
  
  # Get sample names (excluding Taxa and taxonomy columns)
  sample_cols <- setdiff(names(otu_norm_df), c("Taxa", colnames(tax_df)))
  
  results <- list()
  
  for(sample in sample_cols) {
    # Aggregate by taxonomic level (sum abundances for same genus)
    agg_df <- aggregate(otu_norm_df[[sample]], 
                        by = list(Taxon = otu_norm_df[[tax_level]]), 
                        FUN = sum)
    names(agg_df)[2] <- "Abundance"
    
    # Remove NA taxa
    agg_df <- agg_df[!is.na(agg_df$Taxon) & agg_df$Taxon != "", ]
    
    # Sort by abundance descending
    agg_df <- agg_df[order(-agg_df$Abundance), ]
    
    # Calculate cumulative abundance
    agg_df$Cumulative <- cumsum(agg_df$Abundance)
    
    # Find how many taxa needed to reach threshold
    n_for_threshold <- which(agg_df$Cumulative >= cumulative_threshold)[1]
    if(is.na(n_for_threshold)) n_for_threshold <- nrow(agg_df)
    
    # Take at least min_taxa, but add more if needed for threshold
    n_taxa <- max(min_taxa, n_for_threshold)
    n_taxa <- min(n_taxa, nrow(agg_df))  # Don't exceed available taxa
    
    # Get top taxa
    top_taxa <- agg_df[1:n_taxa, ]
    
    # Calculate species count for each genus
    species_counts <- table(otu_norm_df[[tax_level]][otu_norm_df[[tax_level]] %in% top_taxa$Taxon])
    top_taxa$Species_count <- species_counts[top_taxa$Taxon]
    
    # Calculate "Other" category
    other_abundance <- sum(agg_df$Abundance[(n_taxa + 1):nrow(agg_df)])
    other_species_count <- sum(otu_norm_df[[tax_level]] %in% agg_df$Taxon[(n_taxa + 1):nrow(agg_df)])
    
    # Add cumulative percentage for the selected taxa
    top_taxa$Percentage <- top_taxa$Abundance * 100
    top_taxa$Cumulative_pct <- top_taxa$Cumulative * 100
    
    results[[sample]] <- list(
      top_taxa = top_taxa,
      other_abundance = other_abundance,
      other_species_count = other_species_count,
      total_abundance_covered = sum(top_taxa$Abundance),
      n_taxa_selected = n_taxa
    )
  }
  
  return(results)
}
# Get top genera for ALL sample types including Base-Water
top_genera_with_base <- get_top_taxa_by_sample(ps_prok_with_base, 
                                               tax_level = "Genus", 
                                               cumulative_threshold = 0.8, 
                                               min_taxa = 4)

# Calculate species count ratios for stacked barplot
species_ratio_data <- data.frame()

for(sample in names(top_genera_with_base)) {
  result <- top_genera_with_base[[sample]]
  
  # Species in top genera (~80% group)
  top_species_count <- sum(result$top_taxa$Species_count)
  
  # Species in "Other" category (~20% group)  
  other_species_count <- result$other_species_count
  
  # Total species
  total_species <- top_species_count + other_species_count
  
  # Create data for stacked bar
  sample_data <- data.frame(
    Sample = sample,
    Category = c("Top genera (~80%)", "Other genera (~20%)"),
    Species_count = c(top_species_count, other_species_count),
    Percentage = c(top_species_count/total_species * 100, 
                   other_species_count/total_species * 100),
    Abundance_covered = c(result$total_abundance_covered * 100,
                          result$other_abundance * 100)
  )
  
  species_ratio_data <- rbind(species_ratio_data, sample_data)
}

# Print summary
cat("Species count summary by sample type:\n")
species_summary <- species_ratio_data %>%
  pivot_wider(names_from = Category, values_from = c(Species_count, Percentage)) %>%
  mutate(Total_species = `Species_count_Top genera (~80%)` + `Species_count_Other genera (~20%)`) %>%
  select(Sample, Total_species, everything())

print(species_summary)


# Define colors
colors_species_ratio <- c(
  "Top genera (~80%)" = "#F8A058FF",
  "Other genera (~20%)" = "#9B9696A6"
)

# Reorder sample types (put Base-Water first if desired)
sample_order <- c("Base-Water", "Biofilm-no-nutri", "Biofilm-w/-nutri", 
                  "Bulk-no-nutri", "Bulk-w/-nutri")
species_ratio_data$Sample <- factor(species_ratio_data$Sample, levels = sample_order)

# Create the plot
p_species_ratio <- ggplot(species_ratio_data, 
                          aes(x = Sample, y = Species_count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = colors_species_ratio) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "",
    y = "Number of Species", 
    fill = "",
    title = ""
  ) +
  # Add text labels showing counts
  geom_text(aes(label = Species_count), 
            position = position_stack(vjust = 0.5),
            color = "black", fontface = "bold", size = 3.5)

# ------------------------------------------------------------------------ #
# 6. Fig. S9 -- 2 + 3 + 1 panel layout
# ------------------------------------------------------------------------ #

Bflmnn_no_ax <- genus_plots_se[["Biofilm-no-nutri"]]$plot+
theme(axis.title.x = element_blank())

Blkwn_no_ax <- genus_plots_se[["Bulk-w/-nutri"]]$plot+
theme(axis.title.x = element_blank())

# My favorite:
# 3 columns : 1 + 3 + 1 plots
middle_col <- plot_grid(
  Bflmnn_no_ax,
  Blkwn_no_ax,
  genus_plots_se[["Biofilm-w/-nutri"]]$plot,
  nrow = 3,
  align = "v",
  rel_heights = c(0.85, 0.5, 1.1)
)


# Add top top 80% : 20% species ration count barplot
# adjust species ratio plot:
p_species_ratio_cwp <- p_species_ratio+
  ylab("Number\nof Species")+
  guides(fill = guide_legend(
    title = "Relative Abundance",   # Set legend title
    title.position = "top",        # Move title to the left of the keys
    label.position = "right",
    label.hjust = -0.6,        # move text away from the key
    keywidth = unit(0.6, "cm"),
    keyheight = unit(0.6, "cm")
  )) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
    axis.title.y = element_text(angle = 0, face = "bold", 
                                vjust = 0.6, hjust = 0.6,
                                margin = margin(r = -70), size = 12.5),
    legend.direction = "horizontal",
    legend.key.spacing.x = unit(0.1, "cm"), # gap between legend keys
    legend.text = element_text(margin = margin(r = 30), size = 10), # shifts legend position
    legend.key.width = unit(0.6, "cm"),
    legend.margin = margin(b = 0, t = -10)
  )

Blknn_topSpec80 <- plot_grid(
  p_species_ratio_cwp,
  genus_plots_se[["Bulk-no-nutri"]]$plot,
  nrow = 2,
  rel_heights = c(1.6,3),
  align = "v"
)

combined_231 <- plot_grid(
  Blknn_topSpec80,
  middle_col,
  genus_plots_se[["Base-Water"]]$plot,
  ncol = 3,
  rel_widths = rep(1, 3)
)

ggsave(file.path(fig_dir, "FigS9_top80_pct_genera_distribution.png"), combined_231,
       width = 18, height = 14, dpi = 300, bg = "white")

# ------------------------------------------------------------------------ #
# 7. Export the data behind Fig. S9
# ------------------------------------------------------------------------ #

# bars and error bars
write.csv(genus_plot_df_se,
          file = file.path(fig_dir, "FigS9_top80_pct_genera_distribution_plotting_data.csv"),
          row.names = FALSE)

# the overplotted individual replicates
write.csv(genus_individual_df,
          file = file.path(fig_dir, "FigS9_top80_pct_genera_individual_replicates.csv"),
          row.names = FALSE)

# the species-count panel
write.csv(species_ratio_data,
          file = file.path(fig_dir, "FigS9_species_count_80_20_ratio_plotting_data.csv"),
          row.names = FALSE)

cat("\nSaved 'FigS9_top80_pct_genera_distribution.png' into", fig_dir, "\n")

# ======================================================================== #
# PART II -- Fig. 7: CLR PCA + iris plot of the top 15 species
# ======================================================================== #

# ------------------------------------------------------------------------ #
# 8. Remove taxa that break the taxonomic hierarchy
#
#    microViz needs a strictly nested taxonomy. These entries were reported
#    by tax_fix()/ord_plot_iris() as having the same Class/Order/Family under
#    different parents, and are dropped before ordination. They carry very
#    few reads.
#
#    They are removed by rank below, which is deterministic and needs no
#    pre-collected list of rownames.
# ------------------------------------------------------------------------ #

# There were some entries where different 'Phylum' had same 'Family' or 'Class'
# Define problematic families to remove
problematic_taxon_names <- c("Candidatus Bipolaricaulis anaerobius",
                             "Candidatus Fervidibacter sacchari",
                             "Lewinella sp. LCG006",
                             "Candidatus Sulfidibacterium hydrothermale",
                             "Candidatus Xiphinematobacter sp. Idaho Grape",
                             "Candidatus Methylacidiphilum infernorum",
                             "Candidatus Velamenicoccus archaeovorus",
                             "Tichowtungia aerotolerans")

remove_problematic_taxa <- function(ps, problematic_names, level = "rownames", verbose = TRUE) {
  # Remove problematic taxa from a phyloseq object
  # 
  # Parameters:
  # ps: phyloseq object
  # problematic_names: character vector of problematic taxon names to remove
  # level: taxonomic level to search in ("rownames", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  # verbose: logical, whether to print removal statistics
  # 
  # Returns:
  # cleaned phyloseq object
  
  if(length(problematic_names) == 0) {
    if(verbose) message("No problematic taxa provided, returning original object.")
    return(ps)
  }
  
  if(level == "rownames") {
    # Search in rownames (taxa names)
    taxa_to_keep <- !rownames(tax_table(ps)) %in% problematic_names
  } else {
    # Search in specific taxonomic level
    if(!level %in% colnames(tax_table(ps))) {
      warning(paste("Level", level, "not found in tax_table. Skipping removal."))
      return(ps)
    }
    taxa_to_keep <- !tax_table(ps)[, level] %in% problematic_names
  }
  
  # Print removal statistics
  if(verbose) {
    n_remove <- sum(!taxa_to_keep)
    level_label <- ifelse(level == "rownames", "taxa names", paste(level, "level"))
    message("Removing ", n_remove, " taxa with problematic ", level_label, " (", 
            round(n_remove/ntaxa(ps)*100, 1), "% of total)")
  }
  
  # Prune the phyloseq object - this updates ALL components
  prune_taxa(taxa_to_keep, ps)
}

# The two Class-level and four Order-level offenders, and the one Family:
problematic_taxon_names_class <- c("Minisyncoccia", "Fimbriimonadia")

# These were listed when running the iris plot function:
problematic_taxon_names_order <- c("Candidatus Micrarchaeales",
                                   "Methanotrichales",
                                   "Promethearchaeales",
                                   "Thermoplasmatales")

problematic_taxon_names_family <- "Pseudobdellovibrionaceae"

ps_prok_clean <- ps_prok %>%
  remove_problematic_taxa(problematic_taxon_names, level = "rownames") %>%
  remove_problematic_taxa(problematic_taxon_names_class, level = "Class") %>%
  remove_problematic_taxa(problematic_taxon_names_order, level = "Order") %>%
  remove_problematic_taxa(problematic_taxon_names_family, level = "Family")

# Verify remaining taxa counts
ntaxa(ps_prok)
ntaxa(ps_prok_clean)

# ------------------------------------------------------------------------ #
# 9. Sample labels used on the plots
#
#    ggrepel label preparation: E1/E2/E3 become #1/#2/#3, and the single
#    Base-Water sample gets no batch label at all.
# ------------------------------------------------------------------------ #

# Get the sample data as a data frame
# (`sd` elsewhere, where it is local to an lapply() over the
#  confidence levels -- renamed here so it cannot shadow stats::sd)
sam_df <- as(sample_data(ps_prok_clean), "data.frame")

# Modify the batch column
sam_df$batch <- sub("^E", "#", sam_df$batch)

# Set empty batch for Base-Water samples
bw_i <- which(sam_df$Label == "Base-Water")
if(length(bw_i) > 0) {
  sam_df$batch[bw_i] <- ""
}

# Assign the modified data frame back to the phyloseq object
sample_data(ps_prok_clean) <- sample_data(sam_df)

# ------------------------------------------------------------------------ #
# 10. PCA + iris plot
#
#     Confidence 0.1 at Species rank, the combination used throughout.
# ------------------------------------------------------------------------ #

tableau_with20 <- ggthemes::tableau_color_pal(palette = "Tableau 20", type = "regular", direction = 1)(20)

generate_ord_plots <- function(ps,
                                  tax_level = "Species",
                                  n_taxa = 15,
                                  lepalette = c(tableau_with20[c(1,3,2,5,6,7,4,8:15)], "grey20"),
                                  label_color = stepped(n=24)[c(22,5,7,13, 15)])
{
  set.seed(137)
  # 1) Ordination without tax_agg()
  ps_obj <- ps %>%
    tax_fix() %>%
    tax_transform("clr", rank = tax_level, keep_counts = TRUE) %>%
    ord_calc("PCA")
  
  # 1) ordination
  p_ord <- ord_plot(
    ps_obj,
    color = "Label",
    axes  = 1:2,
    size = 4.5, # originally 3
    auto_caption = NA # if we don't want to see info about ordination
  )+
    # if we had enough data points and wanted to draw 4 ellipses
    # stat_ellipse(
    #   aes(linetype = Label, colour = Label),
    #   level = 0.95,
    #   geom = "path",
    #   linewidth = 0.5
    # ) +
    # scale_linetype_manual(
    #   name = "Label",
    #   values = rep("solid", 5)
    # ) +
    # Biofilm ellipse 2 factor variables
    stat_ellipse(
      aes(linetype = "no-nutri"),
      type = "norm",
      level = 0.9,
      geom = "path",
      linewidth = 0.5,
      colour = "#806554",
      data = function(d) d[d$condition == "no-nutri", ]
    ) +
    # Bulk ellipse
    stat_ellipse(
      aes(linetype = "w/-nutri"),
      type = "norm",
      level = 0.9,
      geom = "path",
      linewidth = 0.5,
      colour = "#CCAF9B",
      data = function(d) d[d$condition == "w/-nutri", ]
    ) +
    scale_linetype_manual(
      name = "Type",
      values = c("no-nutri" = "solid", "w/-nutri" = "solid"),
      guide = "none"
    ) +
    theme(legend.title = element_text(face = "bold"))

    # Get ordination coordinates
    # geom_point(aes(color = Label), size = 3)+
    p_ord <- p_ord + 
    geom_label_repel(aes(label = sample_data(ps_obj)$batch), 
                     size = 4.5, # originally 3
                     max.overlaps = 10,
                     show.legend = FALSE,
                     nudge_x = 0.1,
                     nudge_y = 0.1,
                     label.r = 0.2,
                     label.padding = 0.2, label.size = 0.2)+
      scale_color_manual(values = label_color,
                       guide = guide_legend(title = "Sample"))+
    theme(legend.title = element_text(face = "bold"),
          axis.text = element_text(size = 11.5),
          axis.title = element_text(size = 12))
  
  p_ord_for_pair <- p_ord +
    theme(legend.position = "none")
  
  # 2) iris
  p_iris <- ord_plot_iris(
    data               = ps_obj,
    tax_level          = tax_level,
    axes               = 1:2,
    n_taxa             = n_taxa,
    palette            = lepalette,
    anno_colour        = NULL,     ## ← NO segment‐ring
    anno_binary        = NULL,     ## ← NO point‐ring
    ord_plot           = "none",
    bar_width          = 0.8,
    bar_outline_colour = "grey20"
  )
  
  # Figure out where the bars end so we can place our ring just outside
  build   <- ggplot2::ggplot_build(p_iris)
  bar_max <- max(build$data[[1]]$y, na.rm = TRUE)
  ring_y  <- bar_max * 1.15
  label_y <- bar_max * 1.35 # originally 1.35
  
  # extract the same sample order that ord_plot_iris used:
  ps_sorted <- ps_sort_ord(
    ps_obj, 
    ord     = ps_obj,
    axes    = 1:2,
    scaling = 2
  )
  sample_order <- sample_names(ps_sorted)
  
  # extract sample_order + Label into a df
  iris_df <- as(sample_data(ps_sorted), "data.frame") %>%
    tibble::rownames_to_column("Sample") %>%
    dplyr::mutate(sample_order = match(Sample, sample_order))
  
  # Add our own point‐ring, mapping BOTH colour & shape to Label
  p_iris <- p_iris +
    scale_fill_manual(
      name = tax_level,
      values = lepalette)+
    #   guide = guide_legend(
    #     # title = tax_level,
    #     # ncol = 1,
    #     order = 1,
    #     override.aes = list(shape = NA, colour = NA))
    # )+
    geom_point(
      data        = iris_df,
      aes(
        x      = sample_order,
        y      = ring_y,
        colour = Label
        # shape  = Label
      ),
      shape = 16,
      size        = 4.5, # originally 3
      inherit.aes = FALSE,
      show.legend = TRUE
    ) +
    # batch labels above the bars:
    geom_text(
      data = iris_df,
      aes(x = sample_order,
          y = label_y,
          label = batch),
      size = 4.5, # originally 3
      inherit.aes = FALSE,
      show.legend = FALSE
    )+
    # Supply the scales so you get solid circles of your colours
    scale_colour_manual(
      name   = "Label",
      values = label_color,
      guide = guide_legend(
        order = 2,
        override.aes = list(shape = 16, size = 5) # originally 4
      )
    ) +
    guides(
      fill   = guide_legend(title = tax_level, order = 1,
                            override.aes = list(shape = NA, colour = NA)),
      colour = guide_legend(title = "Sample", order = 2, override.aes = list(shape = 16, size = 4.5)), #originally 4
      #   # shape  = guide_legend(override.aes = list(shape = 16, size = 4))
    ) +
    # labs(fill = tax_level) +
    theme(legend.key = element_blank(),
          legend.text = element_text(size = 14),
          legend.box = "vertical",
          legend.title = element_text(face = "bold", size = 13),
          plot.margin = margin(5.5, 20, 5.5, 5.5, "pt"))
  # labs(fill = tax_level)
  # p_iris$labels$fill <- tax_level
  
  # paired
  p_pair <- plot_grid(p_ord_for_pair, p_iris, ncol = 2,
                        rel_widths = c(1, 2),
                        labels = c("a", "b"), # originally c("A", "B")
                        label_size = 13,
                        label_fontface = "bold",
                        label_x = 0.01, label_y = 0.99, hjust = 0, vjust = 1)
  pair_width <- 14 # originally 12

  # save
  ggsave(
    filename = file.path(fig_dir, "Fig7_panel_a_PCA_ellipses.png"),
    plot     = p_ord,
    width    = 6, height = 4, dpi = 300
  )
  ggsave(
    filename = file.path(fig_dir, "Fig7_panel_b_iris_top_species.png"),
    plot     = p_iris,
    width    = 6, height = 6, dpi = 300
  )
  ggsave(
    filename = file.path(fig_dir, "Fig7_PCA_and_iris_paired_plot_top_species.png"),
    plot     = p_pair,
    width    = pair_width, height = 6, dpi = 300
  )

  # return a list if we want to inspect them interactively
  list(ord = p_ord, iris = p_iris, paired = p_pair)
}

ord_plots_prok <- generate_ord_plots(ps_prok_clean,
                                     tax_level = "Species",
                                     n_taxa = 15)

# ------------------------------------------------------------------------ #
# 11. Export the data behind Fig. 7
# ------------------------------------------------------------------------ #

# PCA sample scores (panel a)
pca_scores <- ord_plots_prok$ord$data
head(pca_scores)

# The 90% ellipses are normal-theory ones, so check the scores for normality:
shapiro.test(pca_scores$PC1[pca_scores$type == "Biofilm"])
shapiro.test(pca_scores$PC2[pca_scores$type == "Biofilm"])

write.csv(pca_scores,
          file = file.path(fig_dir, "Fig7_PCA_and_iris_paired_plot_top_species_plotting_data.csv"),
          row.names = FALSE)

# The bars of the iris panel (panel b): one row per species x sample wedge
iris_bar_data <- ggplot2::ggplot_build(ord_plots_prok$iris)$data[[1]]

write.csv(iris_bar_data,
          file = file.path(fig_dir, "Fig7_iris_panel_bar_data.csv"),
          row.names = FALSE)

cat("\nSaved 'Fig7_PCA_and_iris_paired_plot_top_species.png' into", fig_dir, "\n")
