# Combines the per-sample iNEXT rarefaction results written by the
# 03_inext_rarefaction_array.sbatch / 03_inext_rarefaction_worker.R SLURM
# array (13 samples x one q value) into a single object of class "iNEXT",
# usable directly with ggiNEXT(). Run once per q value (0, 1, 2 for this
# study). Runs on a local workstation, not the cluster.
#
# Run from results/inext_rarefaction/ (where the array job wrote its
# per-sample .rds files), e.g.:
#   for q in 0 1 2; do Rscript ../../local/04_merge_rarefaction_results.R $q; done
#
# INPUT   ./inext_q<Q>_col<1..13>_dense_low_sizes_result.rds
#         ../pilot123_diversity_metadata_table.rds
# OUTPUT  ./inext_q<Q>_combined_columnwise_inext_object_dense_low_sizes_result.rds
#         (class "iNEXT"; read by scripts/figures/FigS4_rarefaction_sample_coverage.R)
library(iNEXT)
library(dplyr)

# Metadata table saved once by stage 1 of
# scripts/figures/Fig5_FigS5-S8_alpha_diversity.R, alongside pilot123_diversity.RData. Kept at
# results/pilot123_diversity_metadata_table.rds, one level above this
# script's working directory.
metadata_table_path <- "../pilot123_diversity_metadata_table.rds"
meta_c1 <- readRDS(metadata_table_path)

# Function to create proper iNEXT object
combine_inext_results <- function(inext_list) {
  # Create empty structure for combined result
  combined <- list()

  # Combine DataInfo with proper sample names
  combined$DataInfo <- do.call(rbind, lapply(names(inext_list), function(sample_name) {
    df <- inext_list[[sample_name]]$DataInfo
    df$Assemblage <- sample_name  # Use the name directly
    return(df)
  }))

  # Combine AsyEst
  combined$AsyEst <- do.call(rbind, lapply(names(inext_list), function(sample_name) {
    df <- inext_list[[sample_name]]$AsyEst
    # Add sample name as row names if needed
    rownames(df) <- paste0(sample_name, ".", rownames(df))
    return(df)
  }))

  # Initialize the iNextEst structure
  combined$iNextEst <- list()

  # Combine iNextEst's size_based data
  size_data <- lapply(names(inext_list), function(sample_name) {
    df <- inext_list[[sample_name]]$iNextEst$size_based
    df$Assemblage <- sample_name  # Use actual sample name
    return(df)
  })
  combined$iNextEst$size_based <- do.call(rbind, size_data)

  # Combine iNextEst's coverage_based data
  coverage_data <- lapply(names(inext_list), function(sample_name) {
    df <- inext_list[[sample_name]]$iNextEst$coverage_based
    df$Assemblage <- sample_name  # Use actual sample name
    return(df)
  })
  combined$iNextEst$coverage_based <- do.call(rbind, coverage_data)

  # Set the class to match iNEXT output
  class(combined) <- "iNEXT"

  return(combined)
}

# Create a function to replace assemblage names in the iNEXT object
replace_assemblage_names <- function(inext_obj, mapping) {
  # Create a named vector for easy lookup
  name_map <- setNames(mapping$display_name, mapping$original)

  # Replace names in DataInfo
  inext_obj$DataInfo$Assemblage <- name_map[inext_obj$DataInfo$Assemblage]

  # Replace names in size_based
  inext_obj$iNextEst$size_based$Assemblage <- name_map[inext_obj$iNextEst$size_based$Assemblage]

  # Replace names in coverage_based
  inext_obj$iNextEst$coverage_based$Assemblage <- name_map[inext_obj$iNextEst$coverage_based$Assemblage]

  return(inext_obj)
}

combine_inext_columnwise_results <- function(q_value, n_columns = 13) {

    cat("Combining iNEXT results for q =", q_value, "\n")

    # Initialize list to store results
    sample_results <- list()
    missing_files <- c()

    # Read individual files
    for (col_idx in 1:n_columns) {
        filename <- paste0("inext_q", q_value, "_col", col_idx, "_dense_low_sizes_result.rds")

        if (file.exists(filename)) {
            cat("Reading", filename, "\n")
            result <- readRDS(filename)

            # Extract sample name from the result or use column index
            if (exists("meta_c1") && col_idx <= nrow(meta_c1)) {
                sample_name <- rownames(meta_c1)[col_idx]
            } else {
                sample_name <- paste0("Column_", col_idx)
            }

            sample_results[[sample_name]] <- result

        } else {
            cat("Warning: File", filename, "not found\n")
            missing_files <- c(missing_files, filename)
        }
    }

    if (length(missing_files) > 0) {
        cat("Missing files:\n")
        cat(paste(missing_files, collapse = "\n"), "\n")
    }

    # Create combined structure similar to original parallel output
    if (length(sample_results) > 0) {
        cat("Successfully read", length(sample_results), "individual sample results\n")

        # Create proper iNEXT object using your function
        cat("Creating combined iNEXT object...\n")
        combined_inext_object <- combine_inext_results(sample_results)

        # Apply metadata mapping
        if (exists("meta_c1")) {
            # Create mapping table
            sample_mapping <- data.frame(
                original = rownames(meta_c1),  # Match the names you used
                batch = gsub("pilot", "#", meta_c1$batch),
                Label = meta_c1$Label,
                stringsAsFactors = FALSE
            )

            # Set control batch label
            ctrl_idx <- which(sample_mapping$original == "pilot3_S27_conf01_k2")
            if (length(ctrl_idx) > 0) {
                sample_mapping$batch[ctrl_idx] <- "Ctrl"
            }

            # Create display names
            sample_mapping$display_name <- paste(sample_mapping$batch, sample_mapping$Label, sep=": ")

            # Apply name replacement (pilot1_S174_conf01_k2 -> #1: Biofilm-no-nutri):
            combined_inext_object <- replace_assemblage_names(combined_inext_object, sample_mapping)

            # Check if replacement worked
            cat("First few assemblage names after replacement:\n")
            print(head(unique(combined_inext_object$DataInfo$Assemblage)))
        }

        # Save combined results
        output_base <- paste0("inext_q", q_value, "_combined_columnwise")

        # Save the proper iNEXT object (compatible with ggiNEXT)
        saveRDS(combined_inext_object, paste0(output_base, "_inext_object_dense_low_sizes_result.rds"))

        cat("Combined results saved:\n")
        cat("- Main iNEXT object (for ggiNEXT):", paste0(output_base, "_inext_object_dense_low_sizes_result.rds\n"))

        # Test that the object works with ggiNEXT
        cat("Testing compatibility with ggiNEXT...\n")
        tryCatch({
            test_plot <- ggiNEXT(combined_inext_object, type = 1)
            cat("Combined object is compatible with ggiNEXT.\n")
        }, error = function(e) {
            cat("Warning: Object may not be fully compatible with ggiNEXT:", e$message, "\n")
        })

        return(list(
            inext_object = combined_inext_object,
            individual_results = sample_results
        ))

    } else {
        stop("No valid result files found")
    }
}

# Usage example - combine results for specified q value
if (length(commandArgs(trailingOnly = TRUE)) > 0) {
    q_val <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
} else {
    q_val <- 0  # default
}

cat("Starting combination process for q =", q_val, "\n")
results <- combine_inext_columnwise_results(q_val)

cat("Combination completed for q =", q_val, "\n")
cat("Main iNEXT object ready for plotting with ggiNEXT()\n")

# Example usage after loading:
# combined_data <- readRDS("inext_q0_combined_columnwise_inext_object_dense_low_sizes_result.rds")
# plot <- ggiNEXT(combined_data, type = 1)
