#!/bin/bash
# Writes and submits one SLURM array job per diversity-profile q value.
#
# This is a script that writes a script -- the one deliberate exception to
# this repo's "no generated sbatch files" convention (see scripts/hpc/README.md).
# Reproducing 35 q values as 35 hand-checked-in sbatch files (or reaching for
# a shared config file, which this repo otherwise avoids on purpose) is worse
# than generating each one at submit time from a single template. Kept as-is.
#
# Same per-sample-array workaround as 03_inext_rarefaction_array.sbatch
# (iNEXT has no internal parallelism at the version used): each array task
# handles one of the 13 samples, one CPU, at the q value given on the command
# line.
#
# *** RUNTIME WARNING ***
#   integer q  (0, 1, 2, ...)   -- a few hours per task
#   fractional q (e.g. 1.7)     -- 2-3 DAYS per task
# The full sweep below is 35 q values x 13 samples = 455 tasks, producing
# hundreds of .rds files. Do not launch the whole sweep unless you have days
# of queue time to spend -- consider a handful of q values first.
#
# Usage: ./04_submit_diversity_profile_jobs.sh <q_value>
#
# Launched for the full profile exactly like this, 2s apart so as not to spam
# the scheduler:
#
#   for q in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 \
#            1.8 1.9 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.5 5 5.5 6 7 8; do
#     ./04_submit_diversity_profile_jobs.sh $q
#     sleep 2
#   done
#
# The worker this submits is 03_inext_rarefaction_worker.R -- the same script
# the rarefaction array uses. It takes the q value and the column (sample)
# index as arguments, so one worker serves both steps.
#
# (An older version of this submitter referenced a separate
# "iNEXT_div_prof_colwise.R"; that name is obsolete.)

if [ -z "$1" ]; then
  echo "Usage: $0 <q_value>"
  exit 1
fi

Q_VALUE=$1
REPO_ROOT="$(pwd)"
WORKDIR="${REPO_ROOT}/results/diversity_profile"
mkdir -p "${WORKDIR}"

# The worker reads otu_matrix_prok.rds from whatever directory it is
# started in. That file is written into results/inext_rarefaction/ by stage 1
# of scripts/figures/Fig6_FigS5-S8_alpha_diversity.R, so link it in here as well.
if [ ! -e "${WORKDIR}/otu_matrix_prok.rds" ]; then
  ln -s ../inext_rarefaction/otu_matrix_prok.rds "${WORKDIR}/otu_matrix_prok.rds"
fi

# Create temporary script with expanded variables
temp_script="${WORKDIR}/temp_diversity_profile_q${Q_VALUE}.sbatch"

cat << EOF > ${temp_script}
#!/usr/bin/env bash

#SBATCH --job-name=iNEXT_divprof
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=6-23:59
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kunkli.balazs.tibor@uni-nke.hu
#SBATCH --qos=normal
#SBATCH --output=slurm-%A_%a_diversity_profile.out
#SBATCH --array=1-13

module load singularity

R_CONTAINER="${REPO_ROOT}/containers/r_inext.sif"
CONTAINER_CALL="singularity exec -B ${REPO_ROOT}:${REPO_ROOT} \${R_CONTAINER}"

cd "${WORKDIR}" || exit

# Q_VALUE is baked in below (this file is generated fresh per q); only the
# array task index is resolved at runtime.
echo "Running diversity profile calculation for q = ${Q_VALUE}, column = \${SLURM_ARRAY_TASK_ID}"
\${CONTAINER_CALL} R --vanilla --args "${Q_VALUE}" "\${SLURM_ARRAY_TASK_ID}" < "${REPO_ROOT}/scripts/hpc/03_inext_rarefaction_worker.R"

# The shared worker names its output after the rarefaction step. Rename it to
# the diversity-profile convention, which is what
# scripts/local/05_merge_diversity_profiles.R globs for.
mv "inext_q${Q_VALUE}_col\${SLURM_ARRAY_TASK_ID}_dense_low_sizes_result.rds" \\
   "diversity_profile_q${Q_VALUE}_col\${SLURM_ARRAY_TASK_ID}_result.rds"
EOF

# Submit the job
sbatch ${temp_script}
echo "Job submitted with q_value = ${Q_VALUE}"
echo "Temporary script saved as: ${temp_script}"
