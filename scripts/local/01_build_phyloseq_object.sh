#!/bin/bash

# Combines the per-sample Kraken2/Bracken output (confidence 0.1, produced by
# scripts/hpc/02_kraken2_bracken.sbatch) into a single phyloseq object, imported into
# R for downstream diversity/statistics scripts. Runs on a local workstation,
# not the cluster -- no SLURM here.
#
# Implements the manuscript sentence: "Taxonomic profiles from Kraken2/Bracken
# were combined using scripts from KrakenTools and Bracken and imported into R
# as a merged phyloseq objects."
#
# Uses two external script collections, neither vendored into this repo:
#   - KrakenTools (Lu et al.)   https://github.com/jenniferlu717/KrakenTools
#     -- combine_bracken_outputs.py, combine_kreports.py
#   - KBKen (Tiwari)            https://github.com/SumeetTiwari07/KBKen
#     -- extract_report.py, bracken2phyloseq.R
# Clone both and point KrakenTools/KBKen below at your local checkouts.

# ---- repo-relative paths (edit these if your layout differs) ----
REPO_ROOT="$(pwd)"
k2_report_dir="${REPO_ROOT}/results/kraken_bracken"
metadata="${REPO_ROOT}/data/metadata/metadata_bracken_conf01.csv"

# path to the external script checkouts (see citations above)
KrakenTools="${REPO_ROOT}/../KrakenTools"
KBKen="${REPO_ROOT}/../KBKen"

# python with pandas installed (KBKen's extract_report.py needs it)
pandas_python_path="python3"

# Function to verify file creation
check_file() {
    local file="$1"
    local object_type="$2"
    if [[ -f "$file" ]]; then
        echo "${object_type} created successfully."
    else
        echo "Error: ${object_type} not created."
        exit 1
    fi
}

cd "${k2_report_dir}" || exit

echo "Running 'combine_bracken_outputs.py'..."

"${KrakenTools}/combine_bracken_outputs.py" --files \
    ./*conf01_k2.bracken.report \
    -o ./bracken_estimates_conf01.txt

check_file "./bracken_estimates_conf01.txt" "Bracken estimates file"

echo "-------------------------"
echo "Running 'combine_kreports.py'..."

"${KrakenTools}/combine_kreports.py" \
    -r ./*conf01_k2.report \
    -o ./merged_conf01_kraken_reports.txt

check_file "./merged_conf01_kraken_reports.txt" "Merged Kraken reports file"

# Note: scripts/hpc/02_kraken2_bracken.sbatch already checks each per-sample report
# for the R2/D domain-rank quirk (some Kraken2 PlusPF DB builds tag domain
# level "R2" instead of "D") and fixes it before Bracken ever sees it, so the
# merged file above should already be clean. If you're merging older reports
# that predate that check, re-run the same fix on merged_conf01_kraken_reports.txt
# before the extract_report.py step below:
#   sed 's/\tR2\t/\tD\t/g' merged_conf01_kraken_reports.txt > merged_conf01_kraken_reports_fixed.txt

echo "-------------------------"
echo "Running 'extract_report.py'..."

"${pandas_python_path}" "${KBKen}/extract_report.py" \
    --input ./merged_conf01_kraken_reports.txt \
    --output ./filtered_conf01_kraken_reports.txt \
    --taxonomy-output ./taxonomy_conf01.txt

check_file "./filtered_conf01_kraken_reports.txt" "Filtered Kraken reports file"
check_file "./taxonomy_conf01.txt" "Taxonomy file"

echo "-------------------------"
echo "Running 'bracken2phyloseq.R'..."

Rscript "${KBKen}/bracken2phyloseq.R" \
    -i ./bracken_estimates_conf01.txt \
    -m "${metadata}" \
    -t ./taxonomy_conf01.txt \
    -o brackenSpecies_conf01

check_file "./brackenSpecies_conf01.rds" "Phyloseq object"
