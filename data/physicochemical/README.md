# Physicochemical endpoint data

One tidy CSV per endpoint, long format: one row per measurement. These are the
bulk-water physicochemical measurements taken at the end of each incubation.

`scripts/statistics/01_physicochemical_ART_ANOVA.R` reads the 14 endpoint files listed
below, plus `ph_day0.csv`, the pH Day-0 reference described further down. Two further
datasets from the same experiments live here as well and are described at the end of
this file: `corrosion_rate.csv`
(`scripts/statistics/02_corrosion_rate_statistics.R`) and `surface_roughness.csv`
(`scripts/statistics/03_surface_roughness_statistics.R`).
Every value in these files was verified against the source measurement record.

## Files

| File                 | Column       | Endpoint                | Unit          | n  | Experiments |
|----------------------|--------------|-------------------------|---------------|----|-------------|
| redox_potential.csv  | Redox        | Redox potential (ORP)   | mV            | 36 | #1 #2 #3    |
| dissolved_oxygen.csv | DO           | Dissolved oxygen        | mg/L          | 36 | #1 #2 #3    |
| total_fe.csv         | Fe           | Total iron              | mg/L          | 36 | #1 #2 #3    |
| fe2.csv              | Fe2          | Fe(II)                  | mg/L          | 24 | #2 #3       |
| fe3.csv              | Fe3          | Fe(III)                 | mg/L          | 24 | #2 #3       |
| total_mn.csv         | Mn           | Total manganese         | mg/L          | 36 | #1 #2 #3    |
| sulphate.csv         | SO4          | Sulphate                | mg/L          | 36 | #1 #2 #3    |
| nitrate.csv          | NO3          | Nitrate                 | mg/L          | 36 | #1 #2 #3    |
| nitrite.csv          | NO2          | Nitrite                 | mg/L          | 36 | #1 #2 #3    |
| total_hardness.csv   | TH           | Total hardness          | mg/L as CaCO3 | 36 | #1 #2 #3    |
| chloride.csv         | Cl           | Chloride                | mg/L          | 36 | #1 #2 #3    |
| turbidity.csv        | Turbidity    | Turbidity               | NTU           | 36 | #1 #2 #3    |
| conductivity.csv     | Conductivity | Electrical conductivity | uS/cm         | 36 | #1 #2 #3    |
| ph.csv               | pH           | pH                      | dimensionless | 30 | #1 #2 #3    |

14 endpoints, 474 observations in total.

## Columns

Every file has the same six columns, in this order:

| Column          | Values                                | Meaning |
|-----------------|---------------------------------------|---------|
| *(endpoint)*    | numeric                               | the measured value; always the **first** column, named as in the table above. The analysis script takes `names(d)[1]` as the response variable. |
| `Sterilization` | `Ctrl`, `DW`                          | control vs. disinfected/sterilized water. The script relabels `Ctrl` to `C` for the manuscript figures. |
| `Nutrient`      | `Nutri-`, `Nutri+`                    | without / with nutrient amendment. The script relabels these to `-nutri` / `+nutri`. |
| `Time`          | `35d`, `60d`                          | incubation duration. |
| `Batch`         | `#1`, `#2`, `#3`                      | independent experimental run (= experiment E1, E2, E3). |
| `Replicate`     | `1`, `2`, `3`                         | replicate index within a Batch x Sterilization x Nutrient cell. |

`Time` is fully determined by `Batch` (see below) and is kept as its own column
only because the analysis fits a `Sterilization * Nutrient * Time` model
alongside the `Sterilization * Nutrient * Batch` one.

## Experimental design

Three independent experiments (batches) were run: E1 (`#1`), E2 (`#2`) and
E3 (`#3`). Each is a full 2 x 2 factorial of Sterilization (`Ctrl` / `DW`) and
Nutrient (`Nutri-` / `Nutri+`), with three replicate bottles per treatment
combination, i.e. 12 bottles per experiment. Carbon steel coupons were submerged
in 400 mL water bottles and incubated for 35 days in E1, 60 days in E2 and
35 days in E3, so incubation duration is aliased with batch: `#1` and `#3` are
the 35-day runs and `#2` is the only 60-day run. This is why the analysis treats
"does batch matter" (E1 vs E3, both 35 d) and "does duration matter" (35 d vs
60 d) as two separate questions, and why Batch and Time can never appear in the
same model. E1 and E3 ran for the same length of time but drew source water at
different points in the academic year, which is what the intraclass correlation
coefficients in the analysis script quantify.

Fe(II) and Fe(III) were only measured in E2 and E3, hence 24 rather than 36
rows for those two endpoints, and no E1-vs-E3 batch comparison for them.

## Provenance and known issues

- The unit for total hardness is given as `mg/L as CaCO3` here, matching the
  supplementary data workbook the analysis script writes. Note that the source
  script also carried a second, less specific unit lookup (`mg/L`) that was the
  one in effect for the descriptive-statistics figure's facet strips.

---

## `ph_day0.csv` — 30 rows

pH with a pre-incubation reference measurement. Read by
`scripts/statistics/01_physicochemical_ART_ANOVA.R`, section 8b, for one contrast only.

Autoclaving can itself shift the pH of a solution, so for pH the
`(C,-nutri) - (DW,-nutri)` contrast used for every other endpoint mixes a
sterilization-induced shift with any biological effect. For pH that contrast is
replaced by Day 0 vs. end of incubation within `DW` / `Nutri-`. The `Nutri+`
samples keep the ordinary C vs. DW comparison, because the M9 buffer keeps pH
stable regardless of sterilization.

Same six columns as the endpoint files, with one extra `Time` level:

| Column          | Values                     |
|-----------------|----------------------------|
| `pH`            | numeric                    |
| `Sterilization` | `Ctrl`, `DW`               |
| `Nutrient`      | `Nutri-`, `Nutri+`         |
| `Time`          | `0d`, `35d`, `60d`         |
| `Batch`         | `#1`, `#2`, `#3`           |
| `Replicate`     | `1`, `2`, `3`              |

Per batch, 10 rows: one `0d` measurement of `DW` / `Nutri-` (buffer-free, before
incubation), then three replicates each of `Ctrl` / `Nutri+`, `DW` / `Nutri-` and
`DW` / `Nutri+` at the end of incubation. The end-of-incubation rows are the same
values as in `ph.csv`; `Ctrl` / `Nutri-` is not part of this file. The analysis
uses only the `DW` / `Nutri-` rows (n = 4 per batch).

## `corrosion_rate.csv` — 21 rows

Mass-loss corrosion rate of the carbon steel coupons, um/year. Read by
`scripts/statistics/02_corrosion_rate_statistics.R`. Same convention as above: response
in the first column, design factors after it.

| Column | Values |
|---|---|
| `Corrosion_rate` | numeric, um/year |
| `Duration` | `35d`, `60d` |
| `Sterilization` | `C`, `DW` |
| `Nutrient` | `-nutri`, `+nutri` |

Note the factor levels here are already the manuscript labels (`C`, `-nutri`),
not the `Ctrl` / `Nutri-` forms used by the 14 endpoint files — that is how the
source script wrote them, and the analysis relies on it.

Cells are unbalanced, n = 2..4:

| Duration | Sterilization | Nutrient | n |
|---|---|---|---:|
| 35d | C  | -nutri | 3 |
| 35d | C  | +nutri | 2 |
| 35d | DW | -nutri | 2 |
| 35d | DW | +nutri | 2 |
| 60d | C  | -nutri | 2 |
| 60d | C  | +nutri | 2 |
| 60d | DW | -nutri | 4 |
| 60d | DW | +nutri | 4 |

One 35-day DW/+nutri coupon has a negative rate (-5.13), i.e. a small apparent
mass gain rather than a loss. It is kept as measured. It makes
`log(Corrosion_rate)` NaN for that row, so the log-scale diagnostics in
section 2 of the analysis script are computed on 20 of the 21 observations.

## `surface_roughness.csv` — 156 rows

Areal root-mean-square height (Sq) from 3D profilometry, um. Read by
`scripts/statistics/03_surface_roughness_statistics.R`. Three datasets of different
shape share the file, keyed by `Dataset`; columns that do not apply to a given
dataset are left empty.

| Column | Values |
|---|---|
| `Row` | position within that dataset's original R vector |
| `Dataset` | `35d_polished` (77), `35d_unpolished` (41), `60d` (38) |
| `Duration` | `35d`, `60d` |
| `TimeToIncubation` | `Before`, `After` — `35d_polished` only |
| `Polishing` | `Unpolished`, `Unpolished_side_of_polished` — `35d_unpolished` only |
| `Sterilization` | `C`, `DW` |
| `Nutrient` | `-nutri`, `+nutri` |
| `Surf_roughness` | numeric, um |

The `35d_polished` / `Before` subset contains values at and below zero: these
are readings on a freshly polished, essentially flat coupon, where measurement
noise straddles zero. That is why it cannot be log-transformed.

