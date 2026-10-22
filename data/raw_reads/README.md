# Raw reads

The 26 files in this directory (13 samples x R1/R2) are **0-byte placeholders**. They exist only
to show the expected directory layout and naming convention. The real reads are deposited at NCBI
SRA under BioProject accession **PRJNA1480046**. The sequencing dataset is scheduled for public release upon publication of the associated article.

Once public, download with `prefetch` + `fasterq-dump` (sra-tools) or via the ENA mirror, and drop
the files into this directory under the same names to reproduce the pipeline unchanged.

## Filename convention

```
<experiment>_<sample>_<label>_<read>.fastq.gz
```

Example: `E1_S174_biofilm_no-nutri_R1.fastq.gz`

| Field | Meaning | Example |
|---|---|---|
| 1 | Experiment (`E1`, `E2`, `E3`) | `E1` |
| 2 | Original sequencing sample number, as assigned by the sequencing core | `S174` |
| 3+ | Experiment label used downstream in R (location/nutrient condition, `_`-separated) | `biofilm_no-nutri` |
| last | Read direction | `R1` / `R2` |

Field 3+ encodes location (`biofilm` = carbon steel coupon surface, `bulk` = planktonic water) and
nutrient treatment (`no-nutri` / `w-nutri`), except for `E3_S27_base-water`, which is the
pre-incubation base-water control unique to experiment E3 (see `data/metadata/README.md` for the
full experiment design).
