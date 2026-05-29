# Module Library

Annotated parts catalogue of DNA-Binding Domains (DBDs), Effector Domains (EDs),
and Chromatin Remodelers (CRs) for engineered transcription factor design.

Built to support the computational MEF→iPSC reprogramming project.
Every build step is logged and checksummed so the library can be reconstructed
from scratch at any future point.

---

## Structure

```
module_library/
├── config/
│   └── sources.yaml          # all external URLs, versions, DOIs
├── scripts/
│   ├── utils/
│   │   └── provenance.py     # logging, checksums, manifests
│   ├── schema.py             # SQLite schema + ModuleLibrary class
│   ├── 01_fetch_dbd.py       # AnimalTFDB + UniProt + JASPAR
│   ├── 02_seed_ed.py         # canonical EDs + screen data
│   ├── 03_fetch_cr.py        # EpiFactors + curated CRs
│   ├── 04_build_library.py   # assemble SQLite + TSV export
│   └── 05_validate.py        # QC checks + report
├── data/
│   ├── raw/                  # downloaded files (git-ignored; manifested)
│   │   └── download_manifest.json   # checksums + URLs (commit this)
│   ├── processed/            # intermediate TSVs (git-ignored)
│   └── manual/
│       ├── ed_curated.yaml   # canonical EDs (commit this)
│       ├── cr_curated.yaml   # curated CRs (commit this)
│       └── README.md         # instructions for manual downloads
├── library/
│   ├── module_library.db     # SQLite (git-ignored if large)
│   ├── module_library.tsv    # TSV snapshot (commit this)
│   └── build_manifest.json   # full build provenance (commit this)
├── logs/                     # timestamped run logs (git-ignored)
├── env.yml                   # conda environment
└── README.md
```

---

## Setup

```bash
# 1. Create and activate conda environment
mamba env create -f env.yml
conda activate module_library

# 2. Clone repo and enter project root
cd module_library
```

---

## Running the pipeline

Run scripts in order from the project root:

```bash
# Fetch DBDs from AnimalTFDB, UniProt, JASPAR
python scripts/01_fetch_dbd.py

# Seed canonical EDs (+ screen data if supplementary tables are present)
python scripts/02_seed_ed.py

# Fetch CRs from EpiFactors + curated YAML
python scripts/03_fetch_cr.py

# Assemble into SQLite + export TSV
python scripts/04_build_library.py

# QC checks
python scripts/05_validate.py
```

For a quick test run without network calls:

```bash
python scripts/01_fetch_dbd.py --dry-run --limit 10
python scripts/02_seed_ed.py   --dry-run
python scripts/03_fetch_cr.py  --dry-run
python scripts/04_build_library.py
python scripts/05_validate.py
```

To rebuild from scratch (e.g. after schema change):

```bash
python scripts/04_build_library.py --rebuild
```

---

## Manual downloads required

Some data sources do not permit automated bulk download.
Place these files in `data/manual/` before running the relevant script.

| File | Source | Script that uses it |
|---|---|---|
| `Alerasool_2022_SupTable.tsv` | Suppl. Table 2 from https://doi.org/10.1038/s41588-022-01119-9 | 02_seed_ed.py |
| `DelRosso_2023_SupTable.tsv` | Suppl. Table 2 from https://doi.org/10.1038/s41586-023-06415-8 | 02_seed_ed.py |
| `epifactors.tsv` | https://epifactors.autosome.org → "Download table" | 03_fetch_cr.py (auto-attempted) |

The scripts print clear instructions if a required file is missing.
Download dates and file checksums are recorded in `library/build_manifest.json`.

---

## What gets committed to git

| File / directory | Commit? | Reason |
|---|---|---|
| `data/manual/*.yaml` | ✓ yes | human-curated source of truth |
| `data/manual/*.tsv` | ✓ yes | manually obtained; record for reproducibility |
| `data/raw/download_manifest.json` | ✓ yes | URL + checksum record |
| `data/raw/*.txt`, `*.fasta`, `*.json` | ✗ no | large, re-downloadable |
| `data/processed/*.tsv` | ✗ no | regenerated from scripts |
| `library/module_library.tsv` | ✓ yes | human-readable snapshot; enables git diff |
| `library/module_library.db` | ✗ no | binary; regenerated from TSV |
| `library/build_manifest.json` | ✓ yes | full provenance record |
| `logs/` | ✗ no | per-run logs |

Add to `.gitignore`:
```
data/raw/*.txt
data/raw/*.fasta
data/raw/**/*.json
data/processed/
library/module_library.db
logs/
```

---

## Querying the library

```python
import sys
sys.path.insert(0, "scripts")
from schema import ModuleLibrary

with ModuleLibrary("library/module_library.db") as lib:
    # All activator EDs
    df = lib.to_dataframe("ED")
    activators = df[df["subtype"] == "activator"]

    # Screen-validated DBDs for Homo sapiens
    dbds = lib.to_dataframe("DBD")
    validated = dbds[dbds["validation_level"] == "screen-validated"]

    # Summary counts
    print(lib.counts())
```

---

## Validation levels

| Level | Meaning |
|---|---|
| `predicted` | Computational prediction only; no experimental support |
| `motif-only` | PWM/motif available (e.g. JASPAR); no direct binding evidence |
| `ChIP-validated` | Genome-wide binding data available (ChIP-seq / CUT&RUN) |
| `screen-validated` | Functionally tested in a pooled screen (e.g. Alerasool 2022) |
| `structurally-resolved` | Crystal or cryo-EM structure available for the domain |

The scoring pipeline in `05_score_pipeline.py` (Phase 2) multiplies each component's
contribution by a weight derived from its `validation_level`. Noisy inputs give noisy
rankings — hence confidence over count.

---

## Adding entries manually

Add a new block to `data/manual/ed_curated.yaml` or `cr_curated.yaml`, then re-run:

```bash
python scripts/02_seed_ed.py   # or 03_fetch_cr.py
python scripts/04_build_library.py
python scripts/05_validate.py
```

The `date_modified` field is updated automatically.

---

## Version history

See `library/build_manifest.json` for per-run provenance.
See `CHANGELOG.md` for human-readable version notes.
