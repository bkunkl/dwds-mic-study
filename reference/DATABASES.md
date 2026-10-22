# Reference databases and indices

None of these are shipped with the repository -- they are far too large. Download or build them
once, then point the variable named below at them. Sizes are approximate.

| Database | Version / build | Needed for | Size | Required? |
|---|---|---|---|---|
| Kraken2 **PlusPF** | prebuilt, built **2025-04-02** | `scripts/hpc/02_kraken2_bracken.sbatch` | ~200 GB on disk | **yes** |
| **GRCh38** + **UniVec**, Bowtie2 index | RefSeq `GCF_000001405.40`; UniVec accessed **2024-09-30** | `scripts/hpc/01_preprocess_and_assemble.sbatch` | ~4 GB index (from ~3 GB of FASTA) | **yes** |
| FeGenie curated iron HMMs | bundled with MagicLamp `V1.0` | `scripts/local/03_run_fegenie.sh` | a few MB | **yes**, but no separate download |
| NCBI **nr** (DIAMOND) | release **2024-02-08** | `scripts/local/03_run_fegenie.sh` validation step only | ~200 GB | **no -- optional** |

## Kraken2 PlusPF

Prebuilt indices are published at <https://benlangmead.github.io/aws-indexes/k2>. This study used
the PlusPF build dated **2025-04-02** (bacteria, archaea, viruses, plasmid, UniVec_Core, plus
protozoa and fungi). Download the `.tar.gz` for that date, unpack it, and point `k2DB` in
`scripts/hpc/02_kraken2_bracken.sbatch` at the resulting directory.

The script preloads the database into page cache with `vmtouch` before the sample loop starts.
That is not an optimisation you can skip: the Kraken2 build used here cannot hold the database in
cache across loop iterations by itself, and would re-read it from disk for each of the 13 samples.
Give the job a node with enough RAM to hold it.

Bracken was run against the same database directory (species level, read length 150, default
threshold).

## Decontamination index (GRCh38 + UniVec)

Human host and vector sequence are removed with Bowtie2 before assembly. Build the index once:

```
# GRCh38 -- RefSeq assembly GCF_000001405.40
# UniVec -- https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec  (accessed 2024-09-30)
cat GCF_000001405.40_GRCh38.p14_genomic.fna UniVec > GRCh38_UniVec.fa
bowtie2-build --threads 20 GRCh38_UniVec.fa GRCh38_UniVec
```

Point `DECON_INDEX` in `scripts/hpc/01_preprocess_and_assemble.sbatch` at the index prefix (the path
without the `.1.bt2` suffix). The default there is `reference/GRCh38_UniVec/GRCh38_UniVec`.

## FeGenie HMMs

No separate download. The curated iron HMM set is bundled with MagicLamp `V1.0`, which is what
`MagicLamp.py FeGenie` runs. See `scripts/local/fegenie/README.md`, which also documents an upstream
HMM-lookup issue affecting how many of those HMMs are searched.

## NCBI nr -- optional

`scripts/local/03_run_fegenie.sh` can pass FeGenie's `-ref` flag, which DIAMOND-BLASTps every candidate
gene call against a reference protein database to validate it. The published analysis used NCBI
**nr**, release **2024-02-08**.

**This step is optional.** In `scripts/local/03_run_fegenie.sh`:

```
NR_DB=""                  # empty  -> validation skipped (default)
NR_DB="/path/to/nr.dmnd"  # set    -> validation runs
```

With `NR_DB` empty the core annotation still runs and still writes
`FeGenie-heatmap-data.csv`, which is everything `scripts/figures/Fig9_fegenie_iron_genes.R` needs. Skipping
it saves roughly 200 GB of disk and a large amount of compute; the cost is that gene calls are not
cross-checked against nr.

If you do want it, fetch `nr.gz` from
<https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/> and build a DIAMOND database with
`diamond makedb --in nr.gz --db nr`. Archived releases are needed to reproduce the exact
2024-02-08 version; a current release will give slightly different validation results.
