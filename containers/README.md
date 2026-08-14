# Containers

Singularity `.def` recipes for the tools used by this pipeline. Images are not committed to the
repository (see `.gitignore`); build them before running anything in `scripts/hpc/`.

## Build

```
sudo singularity build bbtools.sif bbtools.def
sudo singularity build bowtie2_samtools.sif bowtie2_samtools.def
sudo singularity build megahit.sif megahit.def
sudo singularity build assembly_stats.sif assembly_stats.def
sudo singularity build metabat2.sif metabat2.def
sudo singularity build kraken2_bracken_vmtouch.sif kraken2_bracken_vmtouch.def
sudo singularity build fastqc_multiqc.sif fastqc_multiqc.def
sudo singularity build magiclamp_fegenie.sif magiclamp_fegenie.def
sudo singularity build r_inext.sif r_inext.def
```

Without local root on the cluster, build on a machine you control (or with
`singularity build --remote`) and copy the `.sif` files over.

## Pipeline tools

| `.def` | Tool(s) | Version(s) |
|---|---|---|
| `bbtools.def` | BBTools / BBDuk | 39.19 |
| `bowtie2_samtools.def` | Bowtie2, SAMtools | 2.5.4, 1.21 |
| `fastqc_multiqc.def` | FastQC, MultiQC | 0.12.1, 1.29 |
| `megahit.def` | MEGAHIT | 1.2.9 |
| `assembly_stats.def` | assembly_stats (Trizna) | 0.1.4 |
| `metabat2.def` | MetaBAT2 (`jgi_summarize_bam_contig_depths`) | 2.18 |
| `kraken2_bracken_vmtouch.def` | Kraken2, Bracken, vmtouch | 2.1.5, 3.1, 1.3.1 |
| `magiclamp_fegenie.def` | MagicLamp (FeGenie genie), DIAMOND | 1.0, 2.1.13 |
| `r_inext.def` | R and the analysis/plotting stack | see below |

## R packages

`r_inext.def` builds this stack.

| Package | Version | Package | Version |
|---|---|---|---|
| R | 4.5.1 | metagenomeSeq | 1.50.0 |
| Bioconductor | 3.21 | ARTool | 0.11.2 |
| iNEXT | 3.0.1 | car | 3.1.5 |
| phyloseq | 1.52.0 | emmeans | 1.11.1 |
| ComplexHeatmap | 2.28.0 | metap | 1.13 |
| microViz | 0.12.7 | metafor | 4.8.0 |
| fantaxtic | 0.2.1 | BiodiversityR | 2.17-2 |
| tidyverse | 2.0.0 | ggplot2 | 3.5.2 |
| tibble | 3.3.0 | purrr | 1.0.4 |
| cowplot | 1.1.3 | paletteer | 1.6.0 |
| ggrepel | 0.9.6 | ggnewscale | 0.5.2 |
| circlize | 0.4.18 | | |

## Notes

- Every recipe except `r_inext.def` uses a `mambaorg/micromamba` base and pins each package through
  `micromamba install` from `conda-forge` + `bioconda`. `r_inext.def` starts from
  `rocker/r-ver:4.5.1` and pins through Bioconductor and `remotes::install_version()`.
- `magiclamp_fegenie.def` builds MagicLamp from its `V1.0` git tag; there is no bioconda package.
  MagicLamp's `FeGenie` genie bundles the curated iron HMMs that `MagicLamp.py FeGenie` runs. The
  optional R plotting stack from MagicLamp's `setup.sh` (`--makeplots`) is omitted, as this
  pipeline does not use it.
- `r_inext.def` sets `LC_ALL=C`. The ART-C contrasts in `scripts/statistics/` assume a collation
  order that differs under that locale; both scripts assert the expected order before running, so a
  mismatch stops the script rather than mislabelling results.
- The FeGenie step (`scripts/local/03_run_fegenie.sh`) runs from a conda environment
  (`conda activate magiclamp`), which is the route documented there. This `.def` reproduces the
  same environment under Singularity/Apptainer.
