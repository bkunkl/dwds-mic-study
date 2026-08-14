#!/usr/bin/env bash

#SBATCH --job-name=inext_dense_low_sizes
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kunkli.balazs.tibor@uni-nke.hu
#SBATCH --qos=normal
#SBATCH --output=slurm-%A_%a_pilot123_dense_low_sizes.out
#SBATCH --array=1-13

# iNEXT rarefaction/extrapolation, one array task per sample (one column of
# the OTU matrix). iNEXT could not parallelise at the version used here, so
# instead of R-level parallelism this splits the 13-sample table across a
# SLURM array and gives each task one CPU.
#
# *** RUNTIME WARNING: budget ~1 day PER SAMPLE. *** --time=72:00:00 gives one
# task headroom; if it's still running near that limit, raise it rather than
# let it get killed. With --array=1-13 that's up to ~13 days of wall-clock
# queue time per q value if the cluster can't run all 13 tasks concurrently.
#
# "dense_low_sizes": 03_inext_rarefaction_worker.R samples the LOW end of each
# rarefaction curve densely (single-read steps up through a few thousand
# reads, coarsening from there to the extrapolation endpoint) because that's
# where the curve shape changes fastest -- coarse steps alone would give an
# unusable, blocky curve near the origin.
#
# Launch once per order of diversity q (Hill number: 0 = richness, 1 = Shannon
# entropy, 2 = Simpson diversity). Q_VALUE is
# not read from the environment automatically -- it must be exported at
# submit time:
#
#   sbatch --export=Q_VALUE=0 03_inext_rarefaction_array.sbatch
#   sbatch --export=Q_VALUE=1 03_inext_rarefaction_array.sbatch
#   sbatch --export=Q_VALUE=2 03_inext_rarefaction_array.sbatch
#
# 13 samples x 3 q values = 39 jobs total for this step alone.

module load singularity

REPO_ROOT="$(pwd)"
R_CONTAINER="${REPO_ROOT}/containers/r_inext.sif"
CONTAINER_CALL="singularity exec -B ${REPO_ROOT}:${REPO_ROOT} ${R_CONTAINER}"

# The worker loads "pilot123_diversity.RData" (an `otu_matrix` of species x 13
# samples, built from the phyloseq object -- see scripts/hpc/README.md for exactly
# where that file comes from) from its working directory, and writes its
# per-sample .rds result there too. Place/symlink that .RData into
# results/inext_rarefaction/ before submitting.
WORKDIR="${REPO_ROOT}/results/inext_rarefaction"
mkdir -p "${WORKDIR}"
cd "${WORKDIR}" || exit

echo "Running iNEXT analysis for q = ${Q_VALUE}, column = ${SLURM_ARRAY_TASK_ID}"
${CONTAINER_CALL} R --vanilla --args "${Q_VALUE}" "${SLURM_ARRAY_TASK_ID}" < "${REPO_ROOT}/scripts/hpc/03_inext_rarefaction_worker.R"
