#!/usr/bin/env bash

#SBATCH --job-name=E123_kraken2_bracken
#SBATCH --partition=bigdata
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=110G
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kunkli.balazs.tibor@uni-nke.hu
#SBATCH --qos=normal
#SBATCH --output=slurm-%j__E123_tax.out

# Taxonomic profiling of the E1-E3 samples with reads at least 140 bp long
# (scripts/hpc/01_preprocess_and_assemble.sbatch's minlength=140 regime). Kraken2 at
# confidence 0.1 (Szóstak et al. 2022; see the top-level README for why), then
# Bracken at species level. The original version of this script also swept
# confidence 0.0 and 0.2 -- dropped, this study only used 0.1.

# ---- repo-relative paths ----
REPO_ROOT="$(pwd)"
WORK140="${REPO_ROOT}/results/working_fastqs_ml140"   # produced by scripts/hpc/01, cd into it below
k2DB="${REPO_ROOT}/reference/k2_pluspf_20250402"       # prebuilt Kraken2 PlusPF, built 2025-04-02,
                                                        # download from https://benlangmead.github.io/aws-indexes/k2
                                                        # -- not shipped here, too large

out_krbr="${REPO_ROOT}/results/kraken_bracken"

# Function to create directory if it doesn't exist
create_dir_if_not_exists() {
	if [ ! -d "$1" ]; then
		mkdir -p "$1"
	fi
}

cleanup_files() {
    local file="$1"
        if [ -f "$file" ]; then
            rm -f "$file"
            echo "Removed $file"
        fi
}

create_dir_if_not_exists "$out_krbr"

cd "${WORK140}" || exit

# load Singularity
module load singularity

# avoid errors due to local language settings:
export LC_ALL=C

# Define container path
KRAKEN2_BRACKEN_CONTAINER="${REPO_ROOT}/containers/kraken2_bracken_vmtouch.sif"
K2B_CONTAINER_CALL="singularity exec -B ${REPO_ROOT}:${REPO_ROOT} ${KRAKEN2_BRACKEN_CONTAINER}"

# Some Kraken2 PlusPF DB builds tag the domain-level row with the expanded
# NCBI clade rank code "R2" instead of the classic "D" that KrakenTools/KBKen
# (see scripts/local/01_build_phyloseq_object.sh) expect. combine_kreports.py just
# concatenates per-sample reports without touching rank codes, so checking
# each sample's own report here -- right after it's written, before Bracken
# reads it -- catches this exactly as reliably as checking the merged report
# later would, without adding a merge step to this script. Detect, report,
# fix only if needed:
check_taxlevel() {
    local report="$1"
    local fixed="${report%.report}_fixed.report"
    if grep -qP '\tR2\t' "${report}"; then
        echo "NOTE: ${report} tags the domain-level rank as R2 (expanded NCBI clade rank), not D. Fixing for KrakenTools/KBKen compatibility:" >&2
        echo "  sed 's/\tR2\t/\tD\t/g' ${report} > ${fixed}" >&2
        sed 's/\tR2\t/\tD\t/g' "${report}" > "${fixed}"
        echo "${fixed}"
    else
        echo "NOTE: ${report} already uses standard D-level rank codes, no fix needed." >&2
        echo "${report}"
    fi
}

# Function to run kraken2 and bracken at confidence 0.1
run_k2_bracken() {
    conf_value="0.1"
    suffix="conf01"
    echo "${s} at confidence level ${conf_value}"

    # Check if input files exist
    echo "Checking input files:"
    ls -la "${s}"_decontam.1.fastq.gz "${s}"_decontam.2.fastq.gz
    if [ ! -f "${s}"_decontam.1.fastq.gz ] || [ ! -f "${s}"_decontam.2.fastq.gz ]; then
        echo "ERROR: Input files for ${s} not found!"
        exit 1
    fi

    # Run kraken2 with verbose output
    echo "Starting Kraken2 at $(date)"
    ${K2B_CONTAINER_CALL} k2 classify \
        --db ${k2DB} \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --output /dev/null \
        --confidence "${conf_value}" \
        --report ${out_krbr}/"${s}"_"${suffix}"_k2.report \
        --paired \
        --use-names \
        "${s}"_decontam.1.fastq "${s}"_decontam.2.fastq
    echo "Finished Kraken2 at $(date)"

    # Check if kraken2 succeeded
    if [ -f "${out_krbr}/${s}_${suffix}_k2.report" ]; then
        echo "Kraken2 report generated successfully"
        k2report=$(check_taxlevel "${out_krbr}/${s}_${suffix}_k2.report")
        ${K2B_CONTAINER_CALL} bracken -d ${k2DB} \
            -i ${k2report} \
            -o ${out_krbr}/"${s}"_"${suffix}"_k2.bracken.report \
            -w ${out_krbr}/"${s}"_"${suffix}"_k2.bracken.newreport \
            -l S -r 150
    else
        echo "ERROR: Kraken2 report not generated!"
        exit 1
    fi
}

# vmtouch is load-bearing here: this Kraken2 version can't hold the ~200 GB
# PlusPF database in page cache across the sample loop below without an
# explicit preload, and would otherwise re-hit disk for every sample.
echo "Preloading Kraken2 database into memory with vmtouch ..."
date
${K2B_CONTAINER_CALL} vmtouch -t "${k2DB}"/*.k2d
echo -e "Finished preloading at $(date)\n"

echo -e "Verifying the vmtouch preload status ... \n"
${K2B_CONTAINER_CALL} vmtouch -v "${k2DB}"/*.k2d

for s in $(ls *_decontam.[1-2].fastq.gz | cut -d '_' -f -2 | uniq); do

    echo "Decompressing ${s}_decontam.[1-2].fastq.gz files..."
    pigz -p "${SLURM_CPUS_PER_TASK}" -k -d "${s}"_decontam.[1-2].fastq.gz

    run_k2_bracken

    # Cleanup decontaminated fastq files
    cleanup_files "${s}"_decontam.1.fastq
    cleanup_files "${s}"_decontam.2.fastq
    echo "Removed ${s}_decontam.[1-2].fastq"
done
