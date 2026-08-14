#!/bin/bash

# Generate the bamfiles_path.txt tables that scripts/local/03_run_fegenie.sh reads with
# MagicLamp's -bams flag. One table per pilot, written into its assembly_megahit
# directory: two columns, tab-separated, no header -- contig FASTA, then the BAM
# it maps to. See data/metadata/E{1,2,3}_bamfiles_path.txt for worked examples.
#
# Matches the naming from scripts/hpc/01_preprocess_and_assemble.sbatch: contigs
# "<sample>_contigs.fa", BAM "<sample>.bam", both under
# results/<pilot>_ml100/assembly_megahit/.

REPO_ROOT="$(pwd)"
ASSEMBLY_ROOT="${REPO_ROOT}/results"

for pilot in pilot{1..3}_ml100; do
    assembly_dir="${ASSEMBLY_ROOT}/${pilot}/assembly_megahit"
    [ -d "$assembly_dir" ] || continue

    cd "$assembly_dir"
    for ctg in *_contigs.fa; do
        [ -e "$ctg" ] || continue
        echo -e "${ctg}\t${assembly_dir}/${ctg%%_*}.bam"
    done > bamfiles_path.txt
    cd "$REPO_ROOT"

    echo "Created: ${assembly_dir}/bamfiles_path.txt"
done
