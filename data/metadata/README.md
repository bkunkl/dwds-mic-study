# Metadata

## Experiment design

Three independent experiments (E1, E2, E3), each a 2x2 factorial: **Nutrient** (no-nutri /
w-nutri) x **Location** (bulk = planktonic water / biofilm = carbon steel coupon surface).
Carbon steel coupons were submerged in 400 mL water bottles.

| Experiment | Incubation | Samples | Base-Water control |
|---|---|---|---|
| E1 | 35 d | 4 (biofilm/bulk x no-nutri/w-nutri) | no |
| E2 | 60 d | 4 (biofilm/bulk x no-nutri/w-nutri) | no |
| E3 | 35 d | 5 (biofilm/bulk x no-nutri/w-nutri + pre-incubation Base-Water) | yes (`E3_S27_base-water`) |

13 samples total. The Base-Water control is filtered source water, DNA-isolated and sequenced
before incubation, unique to E3.

## Files

- **`metadata_bracken_conf01.csv`** — per-sample metadata used downstream in R (`sample`, `batch`,
  `sample_number`, `original_name`, `condition`, `type`, `Label`), keyed to the Kraken2/Bracken
  output at confidence 0.1 (`sample` column, e.g. `E1_S174_conf01_k2`).
- **`E123_metadata.ods`** — original sample sheet (single sheet, `Sheet1`), source of truth.
- **`E123_metadata.csv`** — CSV export of the above via
  `soffice --headless --convert-to csv`. The `.ods` has only one sheet, so nothing is lost in the
  export. The `batch` column in the source spreadsheet contains the string `E1`/`E2`/`E3`.
- **`E1_bamfiles_path.txt`**, **`E2_bamfiles_path.txt`**, **`E3_bamfiles_path.txt`** — two-column
  tables (contig FASTA name -> BAM path, tab-separated, no header) used as the `-bams` input table
  for FeGenie/MagicLamp, one per experiment. Shipped here as worked examples of the expected format;
  the second column must be an absolute path, so `/path/to/repo` stands in for wherever the
  assemblies actually live. `scripts/local/02_make_bam_path_tables.sh` regenerates these tables in place
  with real paths.

## Sulfur bacteria genus lists

`SOB_genus_list_comprehensive.csv` and `SRB_genus_list_comprehensive.csv` are hand-curated genus
lists defining the row set of Fig. 8 (`scripts/figures/Fig8_sulfur_bacteria_heatmap.R`). Columns:

| Column | Meaning |
|---|---|
| `Group` | Ecophysiological grouping (e.g. "Acidophilic SOB", "Thermophilic SRB") |
| `Genus` | Genus name, matched against the phyloseq object's `Genus` rank |
| `Notes` | Free-text justification / provenance for including the genus |

39 genera each. The two lists are disjoint, so every matched genus is unambiguously either a sulfur
oxidiser (SOB) or a sulfate reducer (SRB) — which is what the figure's row ordering relies on.
Only the `Genus` column is read by the script; `Group` and `Notes` document the curation.
