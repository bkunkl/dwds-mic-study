# Results

Every script writes its output under `results/`. None of that output is tracked, with one
exception: the intermediate files the figure scripts read. They are shipped here at the exact
paths the pipeline writes them to, so the figures can be redrawn right after cloning, without
the sequencing reads, the reference databases or the cluster jobs.

Running the pipeline from the raw reads regenerates every one of these files at the same path.

## Shipped intermediate files

| File | Written by | Read by |
|---|---|---|
| `kraken_bracken/brackenSpecies_conf01.rds` | `scripts/local/01_build_phyloseq_object.sh` | `scripts/figures/Fig6_FigS5-S8_alpha_diversity.R` (stage 1), `Fig7_FigS9_taxonomic_composition.R`, `Fig8_sulfur_bacteria_heatmap.R` |
| `diversity_profile/unified_extended_diversity_profiles.rds` | `scripts/local/05_merge_diversity_profiles.R` | `scripts/figures/Fig6_FigS5-S8_alpha_diversity.R` (stage 2) |
| `inext_rarefaction/inext_q{0,1,2}_combined_columnwise_inext_object_dense_low_sizes_result.rds` | `scripts/local/04_merge_rarefaction_results.R` | `scripts/figures/FigS4_rarefaction_sample_coverage.R` |
| `E{1,2,3}_ml100/assembly_megahit/S*_contigs.fa.depth.txt` (13 files) | `scripts/hpc/01_preprocess_and_assemble.sbatch` (`jgi_summarize_bam_contig_depths`) | `scripts/figures/Fig9_fegenie_iron_genes.R` |
| `E{1,2,3}_ml100/assembly_megahit/FeGenie_complete_results/FeGenie-heatmap-data.csv` (3 files) | `scripts/local/03_run_fegenie.sh` | `scripts/figures/Fig9_fegenie_iron_genes.R` |

The depth files per experiment:

| Directory | Samples |
|---|---|
| `E1_ml100/assembly_megahit/` | S174, S175, S176, S177 |
| `E2_ml100/assembly_megahit/` | S216, S217, S218, S219 |
| `E3_ml100/assembly_megahit/` | S1, S2, S3, S4, S27 (Base-Water) |

`brackenSpecies_conf01.rds` is the phyloseq object: Bracken species-level counts, the
Kraken2 taxonomy and the sample metadata of `data/metadata/metadata_bracken_conf01.csv`.
The three combined iNEXT objects are of class `iNEXT` and can be passed straight to
`iNEXT::ggiNEXT()`.

## What is not shipped, and why

- **The per-sample iNEXT results** -- `inext_q<Q>_col<N>_dense_low_sizes_result.rds` (39 files)
  and `diversity_profile_q<Q>_col<N>_result.rds` (455 files) from the cluster jobs. The merge
  scripts combine them into the shipped files above, and the figure scripts read only the
  combined files.
- **`inext_rarefaction/otu_matrix_prok.rds` and `diversity_metadata_table.rds`** -- the input
  of the cluster jobs and the merge scripts. Stage 1 of
  `scripts/figures/Fig6_FigS5-S8_alpha_diversity.R` rebuilds both from the phyloseq object in
  seconds.
- **Everything the scripts write** -- figures, plotting-data CSVs, statistics workbooks.
