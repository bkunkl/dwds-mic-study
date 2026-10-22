# Kraken2 classification rates

`kraken2_classification_rates.csv` -- the fraction of reads Kraken2 could assign to any
taxon, one row per sample. Read by `scripts/figures/Fig6_FigS5-S8_alpha_diversity.R` (stage 2,
section 4).

| Column | Meaning |
|---|---|
| `Sample` | `<batch>_<sample number>`, matching the `sample` column of `../metadata/metadata_bracken_conf01.csv` with the `_conf01_k2` suffix stripped |
| `classification_rate_percent` | percent of reads classified, 0-100 |
| `classification_rate` | the same number as a fraction, 0-1 |

Both columns are carried because the figure script merges the whole table in and its models
use `classification_rate` (the fraction).

## Where the numbers come from

Kraken2 prints a run summary to stderr, e.g.

```
1234567 sequences (370.37 Mbp) processed in 42.123s
  510617 sequences classified (41.36%)
  723950 sequences unclassified (58.64%)
```

`classification_rate_percent` is that classified percentage. Classified + unclassified is
100%; the classified share is what is then split across the taxonomic ranks. If the stderr
logs are gone, the same number can be recomputed from a Kraken2 report
(`results/kraken_bracken/*_conf01_k2.report`): the `U` (unclassified) row's percentage
subtracted from 100, or equivalently the percentage on the root row.

The values here were read off the run summaries by hand and were originally kept as an
`.rds`. They are a CSV here for the same reason as the physicochemical data: 13 numbers in
a text file can be checked by eye, a serialized R object cannot.

## Caveat

These rates are **not comparable across batches as a quality metric alone**. Classification
rate depends on the database version and on how much of the community is represented in it;
low rates here reflect drinking-water biofilm communities that are poorly covered by the
standard database, not failed sequencing. The alpha-diversity models include the rate as a
covariate precisely to check that -- and it comes out non-significant (see section 10 of the
figure script).
