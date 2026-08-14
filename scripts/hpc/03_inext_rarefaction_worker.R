# iNEXT rarefaction/extrapolation for one sample (one column of otu_matrix),
# with a dense grid of sample sizes at the low end and a balanced endpoint.
# Run by 03_inext_rarefaction_array.sbatch, one SLURM array task per sample:
# args are the q value (order of diversity) and the array task ID, which
# indexes into otu_matrix's columns.
#
# *** RUNTIME WARNING: ~1 day per sample on this HPC. *** iNEXT does not
# parallelise at the version used here -- this per-sample split across a
# SLURM array is the workaround, not a speedup within a single run.
#
# INPUT   ./pilot123_diversity.RData  (otu_matrix, meta_c1) in the working
#         directory the sbatch cd'd into
# OUTPUT  ./inext_q<Q>_col<N>_dense_low_sizes_result.rds, one per array task;
#         merged afterwards by scripts/local/04_merge_rarefaction_results.R
#
# The same worker also serves the diversity-profile sweep submitted by
# 04_submit_diversity_profile_jobs.sh, which overrides the output name.

library(iNEXT)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("Error: Need two arguments - q value and column index\n")
    quit(status = 1)
}

hill_n <- as.numeric(args[1])
col_index <- as.numeric(args[2])

# Load the OTU matrix (species x 13 samples) + metadata, built from the
# phyloseq object by stage 1 of scripts/figures/Fig5_FigS5-S8_alpha_diversity.R.
# Expected in the current working directory --
# see scripts/hpc/03_inext_rarefaction_array.sbatch.
load("./pilot123_diversity.RData")

# Calculate smart endpoint (2x median sequencing depth -- balanced between
# giving every sample room to extrapolate and not extrapolating so far that
# the estimate is meaningless).
sample_depths <- colSums(otu_matrix)
median_depth <- median(sample_depths)
endpoint_size <- median_depth * 2

cat("Sample depth range:", min(sample_depths), "to", max(sample_depths), "\n")
cat("Selected endpoint:", endpoint_size, "\n")
cat("Max extrapolation factor:", round(endpoint_size/min(sample_depths), 1), "\n")

# Create well-distributed sample sizes with much denser sampling in lower ranges
# At m=1 (size) some samples may have coverage above 0.5.
target_sizes <- c(
    # Very dense sampling in the critical lower range (where curves change most)
    1, seq(5, 45, by = 5),
    seq(50, 90, by = 10),
    seq(100, 500, by = 50),
    seq(600, 1000, by = 100),
    seq(1000, 10000, by = 1000),
    seq(10000, 45000, by = 5000),
    seq(50000, min(sample_depths), by = 10000),
    seq(min(sample_depths), median_depth, by = 30000),
    seq(median_depth, endpoint_size, by = 50000)
)

# Remove duplicates and ensure we don't exceed endpoint
target_sizes <- sort(unique(target_sizes[target_sizes <= endpoint_size]))

cat("Using", length(target_sizes), "sample sizes\n")
cat("Size range:", min(target_sizes), "to", max(target_sizes), "\n")

# Get sample name and data
sname <- colnames(otu_matrix)[col_index]
sample_data <- otu_matrix[, col_index]
sample_data <- sample_data[sample_data > 0]

# Run analysis
cat("Starting iNEXT with balanced endpoint for q =", hill_n,
    "and sample =", sname, "\n")
start_time <- Sys.time()

# Process only one sample at a time:
inext_out <- iNEXT(sample_data,
                    q = hill_n,
                    datatype = "abundance",
                    size = target_sizes,
                    endpoint = endpoint_size,  # Balanced endpoint
                    knots = length(target_sizes),
                    se = TRUE,
                    conf = 0.95,
                    nboot = 100)

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")
cat("Runtime:", runtime, "minutes\n")

# Save results
output_filename <- paste0("inext_q", hill_n, "_col", col_index, "_dense_low_sizes_result.rds")
saveRDS(inext_out, output_filename)
cat("Results saved to:", output_filename, "\n")
