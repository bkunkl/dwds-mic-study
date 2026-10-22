# HPC scripts

Scripts that run on the cluster under SLURM (`sbatch ...`). Each script defines its own
`#SBATCH` parameters inline -- there is no shared config file to source.

## `01_preprocess_and_assemble.sbatch`

Array job, one task per sample (13 samples, `--array=0-12`). BBDuk trim (adapters + poly-G,
Q20) -> Bowtie2 decontamination against GRCh38 + UniVec -> MEGAHIT assembly (meta-sensitive) ->
contig rename -> Bowtie2 map reads back to contigs -> SAMtools -> `jgi_summarize_bam_contig_depths`.
FastQC runs on the raw and trimmed reads; run MultiQC once after the array finishes (command is
in a comment at the bottom of the script).

Produces two trimmed/decontaminated read sets per sample, both Q20:

- **140 bp minimum length** -- taxonomic profiling input, consumed by `02_kraken2_bracken.sbatch`.
- **100 bp minimum length** -- assembly input, used for everything else in this script.

MetaBAT2 is used for `jgi_summarize_bam_contig_depths`, which produces the per-contig depth
tables that FeGenie's normalisation consumes downstream.

### Contig renaming

MEGAHIT's contig headers are simplified and a minimum length enforced before mapping and depth
calculation. anvi'o provides a dedicated script for this, which is the natural choice inside an
anvi'o pipeline:

```
anvi-script-reformat-fasta "$fa_in" -o "$fa_out" --simplify-names --prefix "$sn" \
  --min-len 1000 --report-file "${sn}_contig-rename-report.txt"
```

This repo replaces it with a three-line `awk` rename instead, so reproducing the pipeline doesn't
require installing anvi'o or accepting its header format if you don't otherwise need it. The
`--min-len 1000` filter is dropped too, since MEGAHIT's own `--min-contig-len 1000` already
enforces it at assembly time.

**Naming is load-bearing downstream:** contig FASTA files are named `<sample>_contigs.fa` (e.g.
`S174_contigs.fa`) and their BAMs `<sample>.bam` (e.g. `S174.bam`), where `<sample>` is the bare
sample number (no experiment prefix). The Fig. 9 script (`scripts/figures/Fig9_fegenie_iron_genes.R`) strips
`_contigs.fa` off the column-source filename to recover the sample id, and the
`data/metadata/E{1,2,3}_bamfiles_path.txt` example tables key off this same pattern. Keep it if you
change anything upstream of this step.

## `02_kraken2_bracken.sbatch`

Kraken2 (confidence 0.1 only -- Szóstak et al. 2022) + Bracken (species level, read length 150,
default threshold) on the 140 bp read set from `01_preprocess_and_assemble.sbatch`. The original
version of this script swept confidence 0.0/0.1/0.2; only 0.1 was used in the manuscript, so the
sweep is gone here.

`vmtouch` preloads the ~200 GB PlusPF database into page cache before the sample loop starts. This
is not optional: this Kraken2 build cannot hold the DB in cache across loop iterations without it,
and would otherwise re-hit disk for every one of the 13 samples.

This script also checks every per-sample Kraken2 report for a domain-rank quirk seen in some
PlusPF DB builds -- tagging the domain level `R2` (an expanded NCBI clade rank) instead of the
classic `D` that KrakenTools/KBKen expect downstream -- and fixes it before Bracken runs, only if
the quirk is actually present. See `scripts/local/01_build_phyloseq_object.sh` for where that matters.

## Databases and indices (not shipped, build/download once)

- **Decontamination index**: GRCh38 (RefSeq `GCF_000001405.40`) + UniVec (accessed 2024-09-30),
  concatenated and `bowtie2-build`, pointed at by `DECON_INDEX` in script `01`.
- **Kraken2 database**: prebuilt PlusPF, built 2025-04-02, from
  <https://benlangmead.github.io/aws-indexes/k2>, pointed at by `k2DB` in script `02`.

## iNEXT rarefaction / diversity-profile chain (`03`, `04`)

**Read this before submitting anything here -- these jobs are a multi-day commitment, not a
quick sbatch.** Sequencing depth varied a lot between the three experiments (E1's fastqs are much
larger than E2/E3's). This chain runs iNEXT rarefaction/extrapolation and a full diversity profile
to show that's not what's driving the treatment effects seen downstream.

**iNEXT does not parallelise at the version used for this study (2H2025).** Everything below
works around that by handing each sample its own CPU via a SLURM array, instead of relying on
R-level parallelism inside a single job.

### Runtime -- read before submitting

| Step | Cost | Notes |
|---|---|---|
| Rarefaction (`03`), per sample | **~1 day** | One array task = one sample = one CPU |
| Diversity profile (`04`), integer q | a few hours | e.g. q = 0, 1, 2 |
| Diversity profile (`04`), fractional q | **2-3 days** | e.g. q = 1.7 |

Job counts:

- Rarefaction: 13 samples x 3 orders of diversity (q = 0, 1, 2) = **39 jobs**.
- Diversity profile: 13 samples x 35 q values (see sweep below) = **455 tasks**, producing
  hundreds of `.rds` files. Fractional q values in that sweep are the expensive ones -- don't
  launch the whole thing unless you have days of queue time to spare. Try one or two q values
  first.

### `03_inext_rarefaction_array.sbatch` + `03_inext_rarefaction_worker.R`

Array job (`--array=1-13`), one task per sample. The worker computes a rarefaction/extrapolation
curve for that sample at a given q, sampling the low end of the curve densely (single-read steps
up through a few thousand reads, coarsening from there to a balanced extrapolation endpoint at 2x
the median sequencing depth) -- that's what "dense_low_sizes" in the filenames means, and it's
there because the curve changes fastest near the origin; coarse-only sampling gives a blocky,
unusable curve there.

Q_VALUE is exported at submit time, not read from the environment automatically -- launched once
per order of diversity:

```
sbatch --export=Q_VALUE=0 03_inext_rarefaction_array.sbatch
sbatch --export=Q_VALUE=1 03_inext_rarefaction_array.sbatch
sbatch --export=Q_VALUE=2 03_inext_rarefaction_array.sbatch
```

Each task writes `inext_q<Q>_col<N>_dense_low_sizes_result.rds` into
`results/inext_rarefaction/`. `scripts/local/04_merge_rarefaction_results.R` (workstation, not cluster)
combines the 13 per-sample files for a given q into one `iNEXT`-class object,
`inext_q<Q>_combined_columnwise_inext_object_dense_low_sizes_result.rds`, ready for `ggiNEXT()`.

### `04_submit_diversity_profile_jobs.sh`

**The one deliberate exception to "no scripts that generate scripts"** in this repo -- kept
because it's how the profile sweep was actually run, and 35 near-identical checked-in sbatch
files (or a shared config file, also avoided here on purpose) would be worse. It writes a fresh
per-q sbatch from a template and submits it.

The standard three alpha-diversity metrics (species richness, Shannon entropy, Simpson diversity)
are q = 0, 1, 2. To see the trend in between, the full profile was computed from q = 0 to 8 at
varying density -- dense near 0-2 (0.1 steps), coarsening further out. Launched like this, one
call per q value, 2 seconds apart:

```
for q in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.5 5 5.5 6 7 8; do
  ./04_submit_diversity_profile_jobs.sh $q
  sleep 2
done
```

Each generated sbatch is itself a 13-task array (one task per sample, same workaround as `03`).
Output naming: `diversity_profile_q<Q>_col<N>_result.rds`, into `results/diversity_profile/`.

The worker it runs is `03_inext_rarefaction_worker.R` -- the same script the rarefaction array
uses. It takes the q value and the column (sample) index as its two arguments, so one worker serves
both steps. (An older version of the submitter referenced a separate `iNEXT_div_prof_colwise.R`;
that name is obsolete and no such file exists.)

Two consequences of sharing the worker, both handled by the submitter:

- The worker reads `otu_matrix_prok.rds` from its own working directory, so the submitter
  symlinks it into `results/diversity_profile/` from `results/inext_rarefaction/`.
- The worker names its output after the rarefaction step, so each generated sbatch renames the
  result to `diversity_profile_q<Q>_col<N>_result.rds` afterwards -- that is what
  `scripts/local/05_merge_diversity_profiles.R` globs for.

`scripts/local/05_merge_diversity_profiles.R` (workstation, not cluster) combines all `.rds` outputs
across the full q sweep into `unified_extended_diversity_profiles.{rds,csv}` (one row per
observed q x sample estimate) and `all_q_extended_iNextEstimates.{rds,csv}` (every coverage-based
point iNEXT computed, not just the observed-size one). The figure script that plots the alpha
diversity profile (`scripts/figures/Fig6_FigS5-S8_alpha_diversity.R`) loads the first
of those two.

### Where `otu_matrix_prok.rds` and the metadata table come from

Both `03`'s worker and the two `local/` merge scripts read two small files derived from the phyloseq
object. They are not shipped, because stage 1 rebuilds them in seconds:

- `otu_matrix_prok.rds` -- an OTU matrix (prokaryotic species x 13 samples), built from the phyloseq
  object (`scripts/local/01_build_phyloseq_object.sh`'s output) by stage 1 of the alpha-diversity analysis
  (`scripts/figures/Fig6_FigS5-S8_alpha_diversity.R`). Needed in `results/inext_rarefaction/` before submitting
  `03`.
- `diversity_metadata_table.rds` -- saved by that same script. Needed at
  `results/diversity_metadata_table.rds` (one level above both `results/inext_rarefaction/`
  and `results/diversity_profile/`) before running either `scripts/local/04` or `scripts/local/05`.

Both are written by stage 1 of that script, which runs before the HPC array jobs; stage 2 resumes
once the results are merged.

### Container

Both `03` and the sbatch template inside `04` expect `containers/r_inext.sif`, built from
`containers/r_inext.def` (R 4.5.1 / Bioconductor 3.21, iNEXT 3.0.1 and the rest of the R analysis
stack). Build it before submitting either step -- see `containers/README.md`.
