#!/usr/bin/env bash

#SBATCH --job-name=Pilot123-assembly
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --time=03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kunkli.balazs.tibor@uni-nke.hu
#SBATCH --qos=normal
#SBATCH --array=0-12
#SBATCH --output=slurm-%A_%a__pilot123_assembly.out

# Preprocessing + assembly. One array task = one sample (13 samples, E1-E3).
#
# Pipeline: BBDuk trim (adapters + poly-G, Q20) -> Bowtie2 decontamination
# (GRCh38 + UniVec) -> MEGAHIT assembly (meta-sensitive) -> contig rename ->
# Bowtie2 map reads back to contigs -> SAMtools -> jgi_summarize_bam_contig_depths.
# FastQC brackets the trimming step; run MultiQC once at the end (see bottom).
#
# Two read-length regimes come out of trimming, both Q20:
#   - minlength=140 -> taxonomic-profiling input, consumed by
#     scripts/hpc/02_kraken2_bracken.sbatch
#   - minlength=100 -> assembly input, everything below in this script uses it
# Both regimes are trimmed and decontaminated here; only the 100 bp set is
# carried on to assembly/mapping/depth in this script.
#
# MetaBAT2 is used for jgi_summarize_bam_contig_depths, which produces the
# per-contig depth tables that FeGenie's normalisation consumes later.

# ---- repo-relative paths (edit these if your layout differs) ----
RAW_DIR="data/raw_reads"
REF_ADAPTERS="reference/adapters_polyG.fasta"
# Bowtie2 index for host + vector decontamination. Build once from GRCh38
# (RefSeq GCF_000001405.40) + UniVec (accessed 2024-09-30) -- not shipped
# here, it's too large. Concatenate both fastas and bowtie2-build them, then
# point this at the resulting index prefix.
DECON_INDEX="reference/GRCh38_UniVec/GRCh38_UniVec"

WORK140="results/working_fastqs_ml140"   # taxonomic-profiling regime
WORK100="results/working_fastqs_ml100"   # assembly regime
QC_DIR="results/qc/fastqc"
ASSEMBLY_ROOT="results"

# Function to create directory if it doesn't exist
create_dir_if_not_exists() {
	if [ ! -d "$1" ]; then
		mkdir -p "$1"
	fi
}

create_dir_if_not_exists "${WORK140}"
create_dir_if_not_exists "${WORK100}"
create_dir_if_not_exists "${QC_DIR}"

# load Singularity
module load singularity

# avoid errors due to local language settings:
export LC_ALL=C

REPO_ROOT="$(pwd)"

# Define container paths
BBTOOLS_CONTAINER="containers/bbtools.sif"
BOWTIE_CONTAINER="containers/bowtie2_samtools.sif"
MEGAHIT_CONTAINER="containers/megahit.sif"
METABAT2_CONTAINER="containers/metabat2.sif"
QC_CONTAINER="containers/fastqc_multiqc.sif"

# Singularity container calls
BB_CALL="singularity exec -B ${REPO_ROOT}:${REPO_ROOT} ${BBTOOLS_CONTAINER}"
BT_CALL="singularity exec -B ${REPO_ROOT}:${REPO_ROOT} ${BOWTIE_CONTAINER}"
MH_CALL="singularity exec -B ${REPO_ROOT}:${REPO_ROOT} ${MEGAHIT_CONTAINER}"
MB_CALL="singularity exec -B ${REPO_ROOT}:${REPO_ROOT} ${METABAT2_CONTAINER}"
QC_CALL="singularity exec -B ${REPO_ROOT}:${REPO_ROOT} ${QC_CONTAINER}"

# ---- identify this array task's sample from the raw R1 files ----
# Raw file naming: <experiment>_<sample>_<label>_<read>.fastq.gz, e.g.
# E1_S174_biofilm_no-nutri_R1.fastq.gz (see data/raw_reads/README.md).
mapfile -t R1_FILES < <(ls "${RAW_DIR}"/E*_R1.fastq.gz | sort)
r1="${R1_FILES[$SLURM_ARRAY_TASK_ID]}"
base=$(basename "${r1}" "_R1.fastq.gz")            # E1_S174_biofilm_no-nutri
r2="${RAW_DIR}/${base}_R2.fastq.gz"

exp=$(echo "${base}" | cut -d_ -f1)                 # E1
snum=$(echo "${base}" | cut -d_ -f2)                # S174
pilot="pilot${exp#E}"                               # pilot1 -- kept for continuity
                                                     # with metadata_bracken_conf01.csv
                                                     # and the pilot1/2/3 naming the
                                                     # rest of the pipeline (incl. the
                                                     # example */bamfiles_path.txt
                                                     # tables) already uses
s="${pilot}_${snum}"                                # pilot1_S174 -- long sample id
sn="${snum}"                                        # S174 -- short id, used below for
                                                     # contig/BAM names

echo "Processing sample: ${base} (s=${s}, sn=${sn})"

out_megahit="${ASSEMBLY_ROOT}/${pilot}_ml100/assembly_megahit"
create_dir_if_not_exists "${out_megahit}"

# ---- FastQC on the raw reads ----
${QC_CALL} fastqc -t "${SLURM_CPUS_PER_TASK}" -o "${QC_DIR}" "${r1}" "${r2}"

# ---- trim + decontaminate one read-length regime ----
# $1 = minlength  $2 = working directory
trim_and_decontam() {
	local minlen="$1"
	local wd="$2"
	local tag="ml${minlen}"

	${BB_CALL} bbduk.sh \
		in="${r1}" \
		in2="${r2}" \
		out="${wd}/${s}_${tag}_R1_bbduk.fastq" \
		out2="${wd}/${s}_${tag}_R2_bbduk.fastq" \
		ref="${REF_ADAPTERS}" \
		k=23 mink=11 hdist=1 hdist2=0 ktrim=r qtrim=rl \
		minlength="${minlen}" trimq=20 tpe tbo \
		-Xmx38g

	${QC_CALL} fastqc -t "${SLURM_CPUS_PER_TASK}" -o "${QC_DIR}" \
		"${wd}/${s}_${tag}_R1_bbduk.fastq" "${wd}/${s}_${tag}_R2_bbduk.fastq"

	${BT_CALL} bowtie2 \
		--quiet \
		-1 "${wd}/${s}_${tag}_R1_bbduk.fastq" \
		-2 "${wd}/${s}_${tag}_R2_bbduk.fastq" \
		-x "${DECON_INDEX}" \
		--very-sensitive-local \
		--dovetail \
		--threads "${SLURM_CPUS_PER_TASK}" \
		-S /dev/null \
		--un-conc "${wd}/${s}_decontam.fastq"

	rm -f "${wd}/${s}_${tag}_R1_bbduk.fastq" "${wd}/${s}_${tag}_R2_bbduk.fastq"
}

trim_and_decontam 140 "${WORK140}"
trim_and_decontam 100 "${WORK100}"

# scripts/hpc/02_kraken2_bracken.sbatch reads the 140 bp set gzip-compressed:
gzip -f "${WORK140}/${s}_decontam.1.fastq" "${WORK140}/${s}_decontam.2.fastq"

# ---- assemble the 100 bp set with MEGAHIT ----
${MH_CALL} megahit \
	-1 "${WORK100}/${s}_decontam.1.fastq" \
	-2 "${WORK100}/${s}_decontam.2.fastq" \
	-o "${out_megahit}/${s}" \
	--out-prefix "${s}" \
	--memory 0.95 \
	--mem-flag 2 \
	--min-contig-len 1000 \
	--presets meta-sensitive \
	-t "${SLURM_CPUS_PER_TASK}"

# Remove intermediate contigs:
rm -rf "${out_megahit}/${s}/intermediate_contigs"

# ---- rename contigs ----
# The original pipeline used anvi'o's anvi-script-reformat-fasta (run inside a
# datascience_genomics.sif container) for this step. That's the dedicated
# tool for it if you're already in an anvi'o workflow -- e.g.
#   anvi-script-reformat-fasta "$fa_in" -o "$fa_out" --simplify-names \
#     --prefix "$sn" --min-len 1000 --report-file "${sn}_contig-rename-report.txt"
# It's replaced here with a plain awk rename so this repo doesn't force an
# anvi'o dependency (and its header format) on anyone who doesn't want one.
# --min-len is dropped too since MEGAHIT's --min-contig-len 1000 above
# already enforces it.
fa_in="${out_megahit}/${s}/${s}.contigs.fa"
fa_out="${out_megahit}/${sn}_contigs.fa"
# Naming matters here: the Fig. 8 script (scripts/figures/Fig8_fegenie_iron_genes.R) strips
# "_contigs.fa" off this filename to get the sample column name, and the
# data/metadata/E{1,2,3}_bamfiles_path.txt example tables key contigs by
# this same "<sn>_contigs.fa" pattern. Keep it.
awk -v pref="${sn}" '/^>/{n++; print ">"pref"_"n; next}{print}' "${fa_in}" > "${fa_out}"

indx="${fa_out%.*}"
${BT_CALL} bowtie2-build --threads "${SLURM_CPUS_PER_TASK}" "${fa_out}" "${indx}"

# ---- map the 100 bp decontaminated reads back to the renamed contigs ----
${BT_CALL} bowtie2 --threads "${SLURM_CPUS_PER_TASK}" \
	-x "${indx}" \
	-1 "${WORK100}/${s}_decontam.1.fastq" \
	-2 "${WORK100}/${s}_decontam.2.fastq" \
	--no-unal \
	--very-sensitive-local \
	-S "${out_megahit}/${sn}.sam"

${BT_CALL} samtools view -F 4 -bS -@ "${SLURM_CPUS_PER_TASK}" "${out_megahit}/${sn}.sam" > "${out_megahit}/${sn}_RAW.bam"
${BT_CALL} samtools sort -@ "${SLURM_CPUS_PER_TASK}" "${out_megahit}/${sn}_RAW.bam" -o "${out_megahit}/${sn}.bam"
${BT_CALL} samtools index -@ "${SLURM_CPUS_PER_TASK}" "${out_megahit}/${sn}.bam"

# Clean up intermediate files
rm -f "${out_megahit}/${sn}_RAW.bam" "${out_megahit}/${sn}.sam"

# ---- per-contig depth, feeds FeGenie normalisation downstream ----
${MB_CALL} jgi_summarize_bam_contig_depths --outputDepth "${fa_out}.depth.txt" "${out_megahit}/${sn}.bam"

echo "Done: ${s} (${sn})"

# MultiQC aggregates every FastQC report under ${QC_DIR} -- run once, after
# the whole array has finished, not per-task:
#   singularity exec containers/fastqc_multiqc.sif multiqc results/qc/fastqc -o results/qc/multiqc
