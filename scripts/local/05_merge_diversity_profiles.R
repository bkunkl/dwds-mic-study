# Combines the hundreds of per-(q, sample) diversity-profile .rds files
# written by 04_submit_diversity_profile_jobs.sh's array jobs into two
# unified tables: the coverage-based estimate at each sample's observed
# sample size (one row per q x sample), and every coverage-based point
# iNEXT computed along the way. Run once, after the whole q sweep has
# finished. Runs on a local workstation, not the cluster.
#
# The resulting "unified_extended_diversity_profiles.rds" is what
# stage 2 of scripts/figures/Fig6_FigS5-S8_alpha_diversity.R loads to plot the
# alpha-diversity profile (FigS7/FigS8).
library(dplyr)

# q values from the sweep launched by 04_submit_diversity_profile_jobs.sh
q_values <- c(
    seq(0, 2, by = 0.1),
    seq(2.25, 4, by = 0.25),
    4.5, 5, 5.5, 6:8
)

# Run from results/diversity_profile/ (where the array jobs wrote their
# per-(q, sample) .rds files). The metadata table is saved once by the
# alpha-diversity analysis (stage 1 of scripts/figures/Fig6_FigS5-S8_alpha_diversity.R)
# alongside otu_matrix_prok.rds, one level up from here.
meta_c1 <- readRDS("../diversity_metadata_table.rds")

# Function to combine results for a specific q value
combine_q_results <- function(q_value, n_columns = 13) {
  estimates_df <- data.frame()
  all_points_df <- data.frame()
  est_df_list <- list()
    missing_files <- c()

    for (col_idx in 1:n_columns) {
        filename <- paste0("diversity_profile_q", q_value, "_col", col_idx, "_result.rds")

        if (file.exists(filename)) {
            result <- readRDS(filename)

            # take diversity estimate at full sample size
            observed_idx <- which(result$iNextEst$coverage_based$Method=="Observed")
            if(length(observed_idx)==1){
              cb_estimate <- result$iNextEst$coverage_based[observed_idx,]
            } else {
              close_to_observed_idx <- nrow(result$iNextEst$coverage_based)
              cb_estimate <- result$iNextEst$coverage_based[close_to_observed_idx,]
            }

            # Create a proper sample identifier
            sample_name <- rownames(meta_c1)[col_idx]

            # Add metadata
            cb_estimate$sample_name <- sample_name
            cb_estimate$batch <- meta_c1$batch[col_idx]
            cb_estimate$Label <- meta_c1$Label[col_idx]

            # Clean up batch names
            if(cb_estimate$Label == "Base-Water"){
              cb_estimate$batch <- "Ctrl"
            } else {
              cb_estimate$batch <- sub(pattern = "^E", replacement = "#", cb_estimate$batch)
            }
            # Remove Assemblage column
            cb_estimate <- cb_estimate[,-1]

            # Add to estimates dataframe
            estimates_df <- rbind(estimates_df, cb_estimate)

            # Process all data points
            all_cb_estimates <- result$iNextEst$coverage_based
            all_cb_estimates$sample_name <- sample_name
            all_cb_estimates$batch <- meta_c1$batch[col_idx]
            all_cb_estimates$Label <- meta_c1$Label[col_idx]

            # Clean up batch names for all points
            if(all_cb_estimates$Label[1] == "Base-Water"){
              all_cb_estimates$batch <- "Ctrl"
            } else {
              all_cb_estimates$batch <- sub(pattern = "^E", replacement = "#", all_cb_estimates$batch)
            }

            # Remove Assemblage column
            all_cb_estimates <- all_cb_estimates[,-1]

            # Add to all points dataframe using rbind
            all_points_df <- rbind(all_points_df, all_cb_estimates)

        } else {
          missing_files <- c(missing_files, filename)
        }
    }

    if (length(missing_files) > 0) {
      cat("Missing files for q =", q_value, ":", length(missing_files), "files\n")
    }

    if (nrow(estimates_df) > 0) {
      est_df_list[[1]] <- estimates_df
      est_df_list[[2]] <- all_points_df
      return(est_df_list)
    } else {
      return(NULL)
    }
}

# Combine results for all q values
cat("Starting to combine diversity profile results...\n")

div_profile_df <- data.frame()
cb_iNextEst_df <- data.frame()

# Combine all q values into one dataframe
for (q_idx in seq_along(q_values)) {
  q_val <- q_values[q_idx]
  cat("Processing q =", q_val, "\n")

  q_result <- combine_q_results(q_val)

  if (!is.null(q_result)) {
    div_profile_df <- rbind(div_profile_df, q_result[[1]])
    div_profile_df <- div_profile_df %>%
      arrange(sample_name, Order.q, Label)
    cb_iNextEst_df <- rbind(cb_iNextEst_df, q_result[[2]])
    cat("  Successfully combined results for q =", q_val, "\n")
  } else {
    cat("  No results found for q =", q_val, "\n")
  }
}

if (nrow(div_profile_df) > 0) {
  # Clean up row names
  rownames(div_profile_df) <- NULL
  rownames(cb_iNextEst_df) <- NULL

  # Save unified dataframe
  saveRDS(div_profile_df, "unified_extended_diversity_profiles.rds")
  write.csv(div_profile_df, "unified_extended_diversity_profiles.csv", row.names = FALSE)

  saveRDS(cb_iNextEst_df, "all_q_extended_iNextEstimates.rds")
  write.csv(cb_iNextEst_df, "all_q_extended_iNextEstimates.csv", row.names = FALSE)

  cat("Unified diversity profiles saved to:\n")
  cat("- unified_extended_diversity_profiles.rds\n")
  cat("- unified_extended_diversity_profiles.csv\n")

  # Summary
  cat("\nSummary of unified results:\n")
  cat("- Total observed/extrapolated estimates:", nrow(div_profile_df), "\n")
  cat("- Total all data points:", nrow(cb_iNextEst_df), "\n")
  cat("- Samples:", length(unique(div_profile_df$sample_name)), "\n")
  cat("- q values:", length(unique(div_profile_df$Order.q)), "\n")
  cat("- q range:", min(div_profile_df$Order.q), "to", max(div_profile_df$Order.q), "\n")

} else {
  cat("No results found for any q value\n")
}
