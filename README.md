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
│   ├── utils.py              # logging, checksums, manifests
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
│       ├── ed_curated.yaml              # canonical EDs
│       ├── cr_curated.yaml              # curated CRs
│       ├── Alerasool_2022_SupTable.xlsx  # screen data (committed)
│       ├── DelRosso_2023_SupTable.xlsx   # screen data (committed)
│       ├── Compendium_2021_SupTable.xlsx # TF effector domain compendium (committed)
│       ├── HiTEff_SupTable.xlsx          # high-throughput effector screen (committed)
│       └── EpiGenes_main.csv             # EpiFactors v2.0 main table (committed)
├── library/
│   ├── module_library.db     # SQLite (git-ignored)
│   ├── module_library.tsv    # TSV snapshot (commit this)
│   └── build_manifest.json   # full build provenance (commit this)
├── logs/                     # timestamped run logs (git-ignored)
├── requirements.txt          # pip dependencies
├── env.yml                   # conda environment
└── README.md
```

---

## Setup

```bash
# Clone and enter project root
git clone https://github.com/JimYuhaoWu/eCR_mod_lib.git
cd eCR_mod_lib

# Install dependencies (pip)
pip install -r requirements.txt

# Or with conda
mamba env create -f env.yml
conda activate module_library
```

Requires Python 3.8+.

---

## Running the pipeline

Run scripts in order from the project root:

```bash
# Fetch DBDs from AnimalTFDB (files committed), UniProt sequences, JASPAR motif IDs
python scripts/01_fetch_dbd.py

# Seed canonical EDs + screen data (Alerasool 2022, DelRosso 2023)
python scripts/02_seed_ed.py

# Fetch CRs from EpiFactors v2.0 + curated YAML
python scripts/03_fetch_cr.py

# Assemble into SQLite + export TSV
python scripts/04_build_library.py

# QC checks
python scripts/05_validate.py
```

For a quick test run (no network calls, 10 TFs per species):

```bash
python scripts/01_fetch_dbd.py --dry-run --limit 10
python scripts/02_seed_ed.py   --dry-run
python scripts/03_fetch_cr.py  --dry-run
python scripts/04_build_library.py
python scripts/05_validate.py
```

To rebuild from scratch (e.g. after a schema change):

```bash
python scripts/04_build_library.py --rebuild
```

---

## Data sources

| Module type | Source | Version | DOI |
|---|---|---|---|
| DBD | [AnimalTFDB](https://guolab.wchscu.cn/AnimalTFDB4/#/Download) | 4.0 | 10.1093/nar/gkad625 |
| DBD sequences | [UniProt](https://www.uniprot.org) | 2024_02 | — |
| DBD motifs | [JASPAR](https://jaspar.elixir.no) | 2024 | 10.1093/nar/gkad1059 |
| ED (curated) | `data/manual/ed_curated.yaml` | manual_v1 | — |
| ED (screen) | [Alerasool et al. 2022](https://doi.org/10.1016/j.molcel.2021.12.005) — tAD-seq sheet | — | 10.1016/j.molcel.2021.12.005 |
| ED (screen) | [DelRosso et al. 2023](https://doi.org/10.1038/s41586-023-05906-y) — Activation + Repression Domains sheets | — | 10.1038/s41586-023-05906-y |
| ED (ChIP-validated) | [Compendium of human TF effector domains, 2021](https://doi.org/10.1016/j.molcel.2021.11.007) — Table S2, Activity H+M only | — | 10.1016/j.molcel.2021.11.007 |
| ED (screen) | High-throughput effector screen (mmc4) — NucAct_data + NucRepr_data sheets | — | TODO — confirm DOI |
| CR (curated) | `data/manual/cr_curated.yaml` | manual_v1 | — |
| CR (EpiFactors) | [EpiFactors](https://epifactors.autosome.org) — EpiGenes_main.csv | 2.0 | 10.1093/nar/gkab1193 |

All supplementary files are committed to `data/manual/` for reproducibility.
The AnimalTFDB `.txt` files are committed to `data/raw/animaltfdb/`.

---

## What gets committed to git

| File / directory | Commit? | Reason |
|---|---|---|
| `data/manual/*.yaml` | ✓ yes | human-curated source of truth |
| `data/manual/*.xlsx` | ✓ yes | supplementary tables; fixed version |
| `data/manual/*.csv` | ✓ yes | EpiFactors main table; fixed version |
| `data/raw/animaltfdb/*.txt` | ✓ yes | AnimalTFDB TF lists; fixed version |
| `data/raw/download_manifest.json` | ✓ yes | URL + checksum record |
| `data/raw/uniprot/*.fasta` | ✗ no | re-downloadable |
| `data/raw/jaspar/*.json` | ✗ no | re-downloadable |
| `data/processed/*.tsv` | ✗ no | regenerated from scripts |
| `library/module_library.tsv` | ✓ yes | human-readable snapshot; enables git diff |
| `library/module_library.db` | ✗ no | binary; regenerated from TSV |
| `library/build_manifest.json` | ✓ yes | full provenance record |
| `logs/` | ✗ no | per-run logs |

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
| `screen-validated` | Functionally tested in a pooled screen (Alerasool 2022, DelRosso 2023) |
| `structurally-resolved` | Crystal or cryo-EM structure available for the domain |

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
