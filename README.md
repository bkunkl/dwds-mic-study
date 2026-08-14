# Microbially influenced corrosion in drinking water distribution systems

Analysis code for a study of microbially influenced corrosion (MIC) in drinking water distribution
systems (DWDS), combining multiple lines of evidence (MLOE):

- **shotgun metagenomics** of bulk water and carbon-steel-coupon biofilm,
- **physicochemical parameters** of the bulk water, and
- **material degradation** -- mass-loss corrosion rate and surface roughness of the coupons.

The repository holds the scripts that produced the manuscript's figures and statistics, the
container recipes for the tools they call, the physicochemical and material data as tidy CSVs, and
the sample metadata. It does not hold the sequencing reads (see
[Data availability](#data-availability)) or the built container images.

The pipeline spans a cluster and a workstation, and the order matters: see
[Running the analysis](#running-the-analysis-end-to-end) and
[Runtime notes](#runtime-notes).

## The study

Three independent experiments (E1, E2, E3), each a 2 x 2 factorial:

| Factor | Levels |
|---|---|
| **Nutrient** | none (`no-nutri`) / added (`w-nutri`) |
| **Location** | bulk water (planktonic) / carbon steel coupon biofilm |

Carbon steel coupons were submerged in 400 mL water bottles and incubated for **35 days** in E1 and
E3 and **60 days** in E2. E3 additionally has a pre-incubation **Base-Water** control -- source
water filtered, DNA-isolated and sequenced before incubation. E1 and E2 have no such control.

| Experiment | Incubation | Metagenomes | Base-Water control |
|---|---|---|---|
| E1 | 35 d | 4 | no |
| E2 | 60 d | 4 | no |
| E3 | 35 d | 5 | yes (`E3_S27_base-water`) |

**13 metagenome samples in total.** In the metadata and in the R scripts the three experiments
appear as `pilot1` / `pilot2` / `pilot3` (E1 / E2 / E3), shown in the figures as #1 / #2 / #3.

The physicochemical and material-degradation measurements come from the same three experiments and
have their own replication: three replicate bottles per treatment cell for the physicochemical
endpoints, 21 coupons for corrosion rate, 156 profilometry measurements for surface roughness. See
`data/physicochemical/README.md`.

## Data availability

Raw sequencing reads are deposited at the NCBI Sequence Read Archive under BioProject accession
**[PRJNA1480046](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1480046)**.

**The dataset is not yet public** -- release is pending manuscript review. `data/raw_reads/`
contains 26 **0-byte placeholder files** (13 samples x R1/R2) showing the expected filenames and
layout. Once the accession is released, download the reads with `prefetch` + `fasterq-dump`
(sra-tools) or from the ENA mirror into that directory under the same names; nothing else needs
changing. See `data/raw_reads/README.md` for the filename convention.

Everything else the analysis needs, apart from the reference databases and the container images, is
in this repository. The physicochemical, corrosion-rate and surface-roughness data are included as
CSVs, so `scripts/statistics/` can be run directly after cloning.

## Repository layout

```
.
├── README.md
├── LICENSE
├── containers/          Singularity .def recipes and software versions (images not shipped)
├── reference/           adapter/poly-G FASTA; DATABASES.md -- what to download and from where
├── scripts/
│   ├── hpc/             cluster scripts (SLURM) and their R worker
│   ├── local/           workstation scripts outside the figure/statistics sets
│   │   └── fegenie/     pinned upstream FeGenie and the patch used here
│   ├── figures/         R scripts producing the manuscript figures
│   └── statistics/      R scripts producing the statistical results
├── data/
│   ├── raw_reads/       0-byte FASTQ placeholders (real reads: SRA PRJNA1480046)
│   ├── metadata/        sample sheet, Kraken2/Bracken metadata, curated SOB/SRB genus lists
│   ├── kraken2/         per-sample Kraken2 classification rates
│   └── physicochemical/ tidy CSVs: 14 water endpoints, corrosion rate, surface roughness
└── results/             all output lands here; gitignored except .gitkeep
```

Directories with details of their own carry a README, linked from here rather than repeated:
[`containers/`](containers/README.md), [`scripts/hpc/`](scripts/hpc/README.md),
[`scripts/local/fegenie/`](scripts/local/fegenie/README.md),
[`data/metadata/`](data/metadata/README.md), [`data/raw_reads/`](data/raw_reads/README.md),
[`data/kraken2/`](data/kraken2/README.md),
[`data/physicochemical/`](data/physicochemical/README.md),
[`reference/DATABASES.md`](reference/DATABASES.md).

## Which scripts run where

| Location | Scripts | Notes |
|---|---|---|
| **Cluster** (SLURM) | `scripts/hpc/*` | Each sbatch carries its own `#SBATCH` block; there is no shared config file. Adjust partition, QOS and `--mail-user` to your site. Tools run from Singularity images built from `containers/`. |
| **Workstation** | `scripts/local/*`, `scripts/figures/*`, `scripts/statistics/*` | None of these need cluster resources. |

Run everything from the repository root unless a script's header says otherwise --
`scripts/local/04` and `scripts/local/05` run from the `results/` subdirectory their inputs were
written to. Paths resolve relative to the working directory (`REPO_ROOT="$(pwd)"` in the shell
scripts, plain relative paths in the R scripts).

## Before you start

1. **Build the containers** you need -- `containers/README.md`, which also lists every tool and
   package version used. Recipes only; images are not distributed.
2. **Get the reference databases** -- `reference/DATABASES.md`. Summary:

   | Database | Version | Needed for | Size | Required? |
   |---|---|---|---|---|
   | Kraken2 **PlusPF** | prebuilt, built **2025-04-02**, <https://benlangmead.github.io/aws-indexes/k2> | `scripts/hpc/02` | ~200 GB | **yes** |
   | **GRCh38** (`GCF_000001405.40`) + **UniVec** (accessed 2024-09-30), Bowtie2 index | build once from both FASTAs | `scripts/hpc/01` decontamination | ~4 GB index | **yes** |
   | FeGenie curated iron HMMs | bundled with MagicLamp `V1.0` | `scripts/local/03` | a few MB | **yes**, no separate download |
   | NCBI **nr** for DIAMOND | release **2024-02-08** | `scripts/local/03` validation step only | ~200 GB | **no -- optional** |

   **nr is optional.** It is used only for FeGenie's `-ref` validation pass, which
   `scripts/local/03_run_fegenie.sh` makes opt-in through its `NR_DB` variable (empty = skip, the
   default). Skipping it still produces the core annotation and `FeGenie-heatmap-data.csv`, which
   is what Fig. 8 needs.
3. **Clone the two external script collections** used by
   `scripts/local/01_build_phyloseq_object.sh` --
   [KrakenTools](https://github.com/jenniferlu717/KrakenTools) and
   [KBKen](https://github.com/SumeetTiwari07/KBKen) -- as siblings of this repository, or edit the
   two path variables at the top of that script.
4. **Fetch and patch FeGenie** -- `scripts/local/fegenie/get_fegenie.sh`. The patch is required:
   stock MagicLamp raises `UnboundLocalError` on the invocation used here, before any heatmap data
   is written. See `scripts/local/fegenie/README.md`.

## Running the analysis, end to end

Four parts. Parts 1-3 are strictly ordered.

### Part 1 -- reads to phyloseq object

```
sbatch scripts/hpc/01_preprocess_and_assemble.sbatch     # cluster, array 0-12
sbatch scripts/hpc/02_kraken2_bracken.sbatch             # cluster
./scripts/local/01_build_phyloseq_object.sh              # workstation
```

`scripts/hpc/01` trims (BBDuk, adapters + poly-G, Q20), decontaminates against GRCh38 + UniVec,
assembles with MEGAHIT, renames contigs, maps reads back and computes contig depths with
`jgi_summarize_bam_contig_depths`. It produces two trimmed read sets per sample: **>= 140 bp** for
taxonomic profiling and **>= 100 bp** for assembly. `scripts/hpc/02` runs Kraken2 at confidence
**0.1** (Lu & Salzberg 2022) and Bracken (species level, read length 150) on the 140 bp set.
`scripts/local/01` combines the per-sample output into the phyloseq object every downstream R
script reads, `results/kraken_bracken/brackenSpecies_conf01.rds`.

### Part 2 -- iron gene annotation

```
./scripts/local/02_make_bam_path_tables.sh               # workstation
./scripts/local/03_run_fegenie.sh                        # workstation
```

`scripts/local/02` writes the `bamfiles_path.txt` tables that MagicLamp reads with `-bams`;
`scripts/local/03` runs the FeGenie genie over the MEGAHIT contigs and depth files from
`scripts/hpc/01`. Patch FeGenie first (`scripts/local/fegenie/`) and set `NR_DB` before launching.

### Part 3 -- alpha diversity

One script is split around the cluster jobs.
`scripts/figures/Fig5_FigS5-S8_alpha_diversity.R` marks the split with banner comments:

```
#  STAGE 1 -- RUNS BEFORE THE HPC JOBS
#  >>> STOP HERE. THE HPC JOBS RUN NOW. <<<
#  STAGE 2 -- NEEDS THE MERGED HPC OUTPUT
```

In order:

1. **Stage 1** of `scripts/figures/Fig5_FigS5-S8_alpha_diversity.R` (workstation) -- from the top
   down to the `STOP HERE` banner. Needs only the phyloseq object. Writes
   `results/inext_rarefaction/pilot123_diversity.RData` and
   `results/pilot123_diversity_metadata_table.rds`, read by the cluster jobs and the merge scripts.
2. `sbatch --export=Q_VALUE=0 scripts/hpc/03_inext_rarefaction_array.sbatch`, and the same for
   `Q_VALUE=1` and `Q_VALUE=2` (cluster) -- 13 samples x 3 q values = **39 jobs**.
3. `./scripts/hpc/04_submit_diversity_profile_jobs.sh <q>` (cluster), once per q value in the sweep
   -- 35 q values x 13 samples. The loop used is in the script's header comment.
4. `scripts/local/04_merge_rarefaction_results.R` (workstation), once per q value, run from
   `results/inext_rarefaction/`.
5. `scripts/local/05_merge_diversity_profiles.R` (workstation), once, run from
   `results/diversity_profile/`, after the sweep has finished.
6. **Stage 2** of `scripts/figures/Fig5_FigS5-S8_alpha_diversity.R` (workstation) -- everything
   below the `STAGE 2` banner. Produces Fig. 5 and Figs. S5-S8. It depends on no object created in
   stage 1 and can be run on its own in a fresh session after the library/path block at the top.
   Besides the merged profiles from step 5 it reads
   `data/kraken2/kraken2_classification_rates.csv`.

Read `scripts/hpc/README.md` and the [runtime notes](#runtime-notes) before starting step 2.

### Part 4 -- remaining figures and statistics

No ordering constraints among these; each needs only what its header lists.

| Script | Needs |
|---|---|
| `scripts/figures/Fig6_FigS9_taxonomic_composition.R` | the phyloseq object (part 1) |
| `scripts/figures/Fig7_sulfur_bacteria_heatmap.R` | the phyloseq object (part 1) + the curated SOB/SRB genus lists in `data/metadata/` |
| `scripts/figures/Fig8_fegenie_iron_genes.R` | contig depth files (part 1) + FeGenie output (part 2) |
| `scripts/figures/FigS4_rarefaction_sample_coverage.R` | the merged iNEXT objects from part 3, step 4 |
| `scripts/statistics/01_physicochemical_ART_ANOVA.R` | `data/physicochemical/` only |
| `scripts/statistics/02_corrosion_rate_statistics.R` | `data/physicochemical/corrosion_rate.csv` |
| `scripts/statistics/03_surface_roughness_statistics.R` | `data/physicochemical/surface_roughness.csv` |

Every figure script also writes the data frame it plotted from as a CSV next to the image, so the
plotted values can be inspected without re-running the analysis.

## Runtime notes

**iNEXT 3.0.1 does not parallelise** -- neither `mc.cores` nor `pbmcapply` has any effect. Part 3
therefore splits the 13-sample table column-wise across SLURM arrays, one sample per task.

| Step | Cost | Jobs |
|---|---|---|
| Rarefaction, `scripts/hpc/03` | **~1 day per sample** | 13 samples x q = 0, 1, 2 = **39 jobs** |
| Diversity profile, `scripts/hpc/04`, **integer** q | a few hours per task | |
| Diversity profile, `scripts/hpc/04`, **fractional** q (e.g. 1.7) | **2-3 days per task** | full sweep = 35 q values x 13 samples = **455 tasks** |

Run one or two q values first and measure them on your own cluster before launching the full sweep.
If a newer iNEXT release provides internal parallelism, the array-based approach can be replaced by
a single job.

`scripts/hpc/01` and `scripts/hpc/02` run in hours, but `scripts/hpc/02` needs a node that can hold
the ~200 GB Kraken2 database in page cache; `vmtouch` preloads it.

## Figure map

| Figure | Output file | Produced by |
|---|---|---|
| Fig. 5 | `Fig5_alpha_diversity_metrics.png` | `scripts/figures/Fig5_FigS5-S8_alpha_diversity.R` |
| Fig. 6 | `Fig6_PCA_and_iris_paired_plot_top_species.png` | `scripts/figures/Fig6_FigS9_taxonomic_composition.R` |
| Fig. 7 | `Fig7_complexHeatmap_sulfur_bacteria.png` | `scripts/figures/Fig7_sulfur_bacteria_heatmap.R` |
| Fig. 8 | `Fig8_FeGenie_iron_related_functional_gene_categories.png` | `scripts/figures/Fig8_fegenie_iron_genes.R` |
| Fig. S2 | `FigS2_heatmap_treatment_contrasts.png` | `scripts/statistics/01_physicochemical_ART_ANOVA.R` |
| Fig. S3 | `FigS3_heatmap_Batch-pair_comparisons.png` | `scripts/statistics/01_physicochemical_ART_ANOVA.R` |
| Fig. S4 | `FigS4_type2_sample_coverage.png` | `scripts/figures/FigS4_rarefaction_sample_coverage.R` |
| Fig. S5 | `FigS5_model_coefficients_across_alpha_diversity_metrics.png` | `scripts/figures/Fig5_FigS5-S8_alpha_diversity.R` |
| Fig. S6 | `FigS6_standardized_treatment_contrasts_cohens_d.png` | `scripts/figures/Fig5_FigS5-S8_alpha_diversity.R` |
| Fig. S7 | `FigS7_alpha_diversity_profile.png` | `scripts/figures/Fig5_FigS5-S8_alpha_diversity.R` |
| Fig. S8 | `FigS8_alpha_diversity_profile_derivative.png` | `scripts/figures/Fig5_FigS5-S8_alpha_diversity.R` |
| Fig. S9 | `FigS9_top80_pct_genera_distribution.png` | `scripts/figures/Fig6_FigS9_taxonomic_composition.R` |

All are written to `results/figures/`. The remaining manuscript figures are schematics and
photographs, not script output. Fig. 6 also writes its two panels separately, and several scripts
produce companion figures and tables beyond those listed; statistical workbooks go to
`results/statistics/`. Each script header lists its full output.

## Software versions

Every tool and R package version, and the recipes that pin them, are in
[`containers/README.md`](containers/README.md).

Read filters: minimum length **140 bp** for taxonomic profiling, **100 bp** for assembly, Q20 for
both.

## Citations

If you use parts of this pipeline, please cite the underlying tools as well as the manuscript:

- **KrakenTools** -- <https://github.com/jenniferlu717/KrakenTools> (Kraken2/Bracken report
  merging).
- **KBKen** (S. Tiwari) -- <https://github.com/SumeetTiwari07/KBKen> (Bracken -> phyloseq import).
- **FeGenie** -- <https://github.com/Arkadiy-Garber/FeGenie>.
- **MagicLamp** -- <https://github.com/Arkadiy-Garber/MagicLamp>. FeGenie is run through MagicLamp
  here; script provenance is pinned in `scripts/local/fegenie/`.
- **Lu & Salzberg (2022)** -- basis for the Kraken2 confidence threshold of 0.1.
- **Sequence data** -- NCBI SRA BioProject **PRJNA1480046**.

## License

[MIT](LICENSE) -- permissive and standard for research code, anyone can freely reuse or adapt the
analysis.
