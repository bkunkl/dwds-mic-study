#!/bin/bash

# Runs MagicLamp's FeGenie genie (iron-cycling gene annotation) over the MEGAHIT
# contigs + depth files from scripts/hpc/01_preprocess_and_assemble.sbatch, keyed by the
# bamfiles_path.txt tables from scripts/local/02_make_bam_path_tables.sh.
#
# Runs on a workstation/faculty server, not the cluster -- no SLURM here.
#
# Prerequisites:
#   conda activate magiclamp
# (a Singularity alternative also exists -- containers/magiclamp_fegenie.def --
#  substitute "singularity exec containers/magiclamp_fegenie.def MagicLamp.py ..."
#  for the MagicLamp.py calls below if you'd rather not set up conda)
#
# MagicLamp's stock genies/FeGenie.py crashes on this exact invocation (-ref, no
# --gbk) with an UnboundLocalError, before it ever writes FeGenie-heatmap-data.csv.
# Get the patched script with scripts/local/fegenie/get_fegenie.sh and copy it over your
# MagicLamp install's genies/FeGenie.py before running this -- see
# scripts/local/fegenie/README.md for the full analysis.

REPO_ROOT="$(pwd)"

# ---- optional DIAMOND-vs-nr validation ----
# FeGenie's -ref flag BLASTs (via DIAMOND) every candidate hit against a reference
# protein database to validate it -- here, NCBI nr (release 2024-02-08, ~200 GB).
# Set NR_DB to a local nr path to run it; leave empty to skip. Skipping still
# produces the core annotation (FeGenie-heatmap-data.csv), everything Fig9 needs.
NR_DB=""

run_FeGenie () {
    local exp="$1"
    local w_dir="${REPO_ROOT}/results/${exp}/assembly_megahit"
    cd "${w_dir}" || exit 1

    local logFile="${exp}_megahit_FeGenie_complete.log"
    echo -e "\nStarting $(date) \nMagicLamp.py FeGenie\n\n" >> "$logFile"

    local ref_args=()
    [ -n "$NR_DB" ] && ref_args=(-ref "$NR_DB")

    MagicLamp.py FeGenie \
        -bin_dir "${w_dir}/" \
        -bin_ext fa \
        -out "${w_dir}/FeGenie_complete_results/" \
        -bams "${w_dir}/bamfiles_path.txt" \
        "${ref_args[@]}" \
        --nohup \
        --meta -t 4 >> "$logFile" 2>&1

    cd "${REPO_ROOT}"
}

for exp in E{1..3}_ml100; do
    run_FeGenie "$exp"
done
