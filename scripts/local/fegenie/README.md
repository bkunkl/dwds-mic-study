# FeGenie

Iron-related genes were annotated with **MagicLamp v1.0**, which wraps the curated HMMs of
**FeGenie v1.2**. Predicted proteins were aligned to the NCBI nr database (release 2024-02-08)
with DIAMOND BLASTp.

FeGenie itself is not redistributed here. `get_fegenie.sh` fetches the exact upstream revision
and applies the patch in this directory:

```
./get_fegenie.sh /some/dir
```

It produces one file

| File | What it is |
|---|---|
| `FeGenie_asrun.py` | Exactly what produced the published results (Fig. 9) |

Copy it over MagicLamp's `genies/FeGenie.py`.

Upstream pin: MagicLamp `genies/FeGenie.py` @ commit
`fc7b6f7e3b99897bab744e9cc0363a6f6c46c459` (2022-11-25),
sha256 `39d285dff9d84deebf1ae6f0c61c4323d649e43d977d5e763717f3de7831415b`.
That file has not changed since 2022, but the repository HEAD has, so the commit is pinned rather
than `master`.

## Patch — `FeGenie_asrun.patch` (required)

27 hunks, 63 lines. Applying it to the pinned upstream reproduces the analysis script used in this
study byte-for-byte. Without it the pipeline does not run to completion:

- **Crash fix (load-bearing).** `idxDict` is bound only inside the `if args.gbk:` branch but is
  referenced in the `else:` branch. This analysis passes `-ref` and no `--gbk`, so the `else:`
  branch is always taken, raising `UnboundLocalError`. The crash occurs *before*
  `FeGenie-heatmap-data.csv` is written, so an unpatched run produces no heatmap data at all.
- **CSV formatting.** 20 sites were missing the `","` separator between the `top_blast_hit` and
  `blast_hit_evalue` fields, corrupting the output table. Six neighbouring sites in the same file
  are correct, so this is an upstream copy-paste slip.
- Two depth-file path corrections (FeGenie looks for `.depth.txt` files in the working directory),
  four `metaDict` field changes, and one commented-out no-op `rm`.

The newer standalone FeGenie repository (https://github.com/Arkadiy-Garber/FeGenie) fixes the
crash but **still contains all 20 missing-comma sites**, and differs from the MagicLamp version by
93 hunks overall — it would not reproduce the results reported here.

