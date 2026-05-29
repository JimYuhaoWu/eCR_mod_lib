# Module Library

Annotated parts catalogue of DNA-Binding Domains (DBDs), Effector Domains (EDs),
and Chromatin Remodelers (CRs) for engineered transcription factor design.

Built to support the computational MEF→iPSC reprogramming project.
Every build step is logged and checksummed so the library can be reconstructed
from scratch at any future point.

---

## Current library size

**8,675 records** — last built 2026-05-30

| Type | Count |
|---|---|
| DBD | ~3,250 |
| ED | ~4,900 |
| CR | ~540 |

| Validation level | Meaning |
|---|---|
| `screen-validated` | Functionally tested in a pooled screen |
| `ChIP-validated` | Experimentally characterized (Gal4-fusion assays, ChIP-seq) |
| `motif-only` | PWM/motif available (JASPAR); no direct binding evidence |
| `predicted` | Computational prediction only |
| `structurally-resolved` | Crystal or cryo-EM structure available |

---

## Structure

```
module_library/
├── config/
│   └── sources.yaml               # all external URLs, versions, DOIs
├── scripts/
│   ├── utils.py                   # logging, checksums, manifests
│   ├── schema.py                  # SQLite schema + ModuleLibrary class
│   ├── 01_fetch_dbd.py            # AnimalTFDB + UniProt + JASPAR
│   ├── 02_seed_ed.py              # canonical EDs + screen data
│   ├── 03_fetch_cr.py             # EpiFactors + curated CRs
│   ├── 04_build_library.py        # assemble SQLite + TSV export
│   └── 05_validate.py             # QC checks + report
├── literature/
│   ├── papers.yaml                # master record of all assessed papers
│   ├── search_queries.yaml        # PubMed + Semantic Scholar queries
│   ├── 01_search.py               # fetch new candidates
│   ├── 02_triage.py               # score + generate review report
│   ├── 03_record.py               # validate review decisions
│   ├── candidates/                # raw JSON search results (git-ignored)
│   └── reviews/                   # markdown review reports (committed)
├── data/
│   ├── raw/                       # downloaded files (git-ignored; manifested)
│   │   ├── animaltfdb/            # AnimalTFDB .txt files (committed)
│   │   └── download_manifest.json # checksums + URLs (committed)
│   ├── processed/                 # intermediate TSVs (git-ignored)
│   └── manual/                    # supplementary files (committed)
│       ├── ed_curated.yaml
│       ├── cr_curated.yaml
│       ├── Alerasool_2022_SupTable.xlsx
│       ├── DelRosso_2023_SupTable.xlsx
│       ├── Compendium_2021_SupTable.xlsx
│       ├── HiTEff_SupTable.xlsx
│       ├── Staller_2022_mmc2.csv
│       ├── Staller_2022_mmc4.csv
│       ├── Staller_2022_mmc5.csv
│       └── EpiGenes_main.csv
├── library/
│   ├── module_library.db          # SQLite (git-ignored)
│   ├── module_library.tsv         # TSV snapshot (committed)
│   ├── build_manifest.json        # full build provenance (committed)
│   └── qc_report.tsv              # QC check results (committed)
├── logs/                          # timestamped run logs (git-ignored)
├── requirements.txt               # pip dependencies
├── env.yml                        # conda environment
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
```

Requires Python 3.8+.

---

## Running the pipeline

Run scripts in order from the project root:

```bash
python scripts/01_fetch_dbd.py        # AnimalTFDB files committed; fetches UniProt + JASPAR
python scripts/02_seed_ed.py          # all ED sources
python scripts/03_fetch_cr.py         # EpiFactors + curated CRs
python scripts/04_build_library.py    # assemble SQLite + export TSV
python scripts/05_validate.py         # QC checks + report
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

### DBD

| Source | Version | DOI |
|---|---|---|
| [AnimalTFDB 4.0](https://guolab.wchscu.cn/AnimalTFDB4/#/Download) | 4.0 | 10.1093/nar/gkad625 |
| [UniProt](https://www.uniprot.org) (sequences) | 2024_02 | — |
| [JASPAR](https://jaspar.elixir.no) (motif IDs) | 2024 | 10.1093/nar/gkad1059 |

### ED

| Source | Description | Validation | DOI |
|---|---|---|---|
| `ed_curated.yaml` | Canonical EDs (VP64, VPR, KRAB, p65, etc.) | ChIP-validated | manual |
| [Alerasool et al. 2022](https://doi.org/10.1016/j.molcel.2021.12.005) | tAD-seq screen, hit activators | screen-validated | 10.1016/j.molcel.2021.12.005 |
| [DelRosso et al. 2023](https://doi.org/10.1038/s41586-023-05906-y) | Activation + Repression Domains sheets | screen-validated | 10.1038/s41586-023-05906-y |
| [Tycko et al. 2020](https://doi.org/10.1016/j.cell.2020.11.024) | NucAct + NucRepr Pfam-domain screens | screen-validated | 10.1016/j.cell.2020.11.024 |
| [Compendium 2021](https://doi.org/10.1016/j.molcel.2021.11.007) | Table S2, Activity H+M (404 AD + 324 RD) | ChIP-validated | 10.1016/j.molcel.2021.11.007 |
| [Staller et al. 2022](https://doi.org/10.1016/j.cels.2022.01.002) | 50 predicted novel ADs (Z≥0.5) + 5 WT domains | screen-validated | 10.1016/j.cels.2022.01.002 |

### CR

| Source | Description | DOI |
|---|---|---|
| `cr_curated.yaml` | Key CRs with full annotation | manual |
| [EpiFactors 2.0](https://epifactors.autosome.org) | EpiGenes_main.csv, enzymatic CRs only | 10.1093/nar/gkab1193 |

---

## What gets committed to git

| File / directory | Commit? | Reason |
|---|---|---|
| `data/manual/*.yaml` | ✓ yes | human-curated source of truth |
| `data/manual/*.xlsx` | ✓ yes | supplementary tables; fixed version |
| `data/manual/*.csv` | ✓ yes | supplementary tables + EpiFactors |
| `data/raw/animaltfdb/*.txt` | ✓ yes | AnimalTFDB TF lists; fixed version |
| `data/raw/download_manifest.json` | ✓ yes | URL + checksum record |
| `data/raw/uniprot/*.fasta` | ✗ no | re-downloadable |
| `data/raw/jaspar/*.json` | ✗ no | re-downloadable |
| `data/processed/*.tsv` | ✗ no | regenerated from scripts |
| `library/module_library.tsv` | ✓ yes | human-readable snapshot; enables git diff |
| `library/module_library.db` | ✗ no | binary; regenerated |
| `library/build_manifest.json` | ✓ yes | full provenance record |
| `library/qc_report.tsv` | ✓ yes | QC summary |
| `literature/papers.yaml` | ✓ yes | master record of paper decisions |
| `literature/search_queries.yaml` | ✓ yes | search configuration |
| `literature/reviews/*.md` | ✓ yes | review reports |
| `literature/candidates/*.json` | ✗ no | re-generatable search results |
| `logs/` | ✗ no | per-run logs |

---

## Literature search

Semi-automated pipeline for discovering new papers. Run every 1–3 months.

```bash
# 1. Query PubMed + Semantic Scholar; deduplicate against known papers
python literature/01_search.py

# 2. Score candidates and generate a Markdown review report
python literature/02_triage.py
# → literature/reviews/review_YYYY-MM-DD.md

# 3. Add accepted/rejected entries to literature/papers.yaml, then validate:
python literature/03_record.py
```

**Rejection patterns established so far** (auto-filtered by triage scorer):
- Inhibitor / drug / therapeutic studies
- Cancer biology without domain characterization
- Synthetic lethality studies
- Non-human/mouse organisms (plants, yeast, bacteria, Drosophila, etc.)
- Reviews and methods book chapters
- Assembly/complex stoichiometry screens (wrong question)
- RAS/signaling effector domains (not transcriptional effectors)

**After accepting a paper:** place its supplementary data in `data/manual/`,
add an entry to `SCREEN_FILES` in the relevant script, and re-run the pipeline.
Update `status: integrated` in `papers.yaml` when done.

See [`literature/papers.yaml`](literature/papers.yaml) for all decisions.
See [`literature/search_queries.yaml`](literature/search_queries.yaml) to tune search terms.

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

## Adding entries manually

Add a new block to `data/manual/ed_curated.yaml` or `cr_curated.yaml`, then re-run:

```bash
python scripts/02_seed_ed.py   # or 03_fetch_cr.py
python scripts/04_build_library.py
python scripts/05_validate.py
```

---

## Version history

See `library/build_manifest.json` for per-run provenance.
See `CHANGELOG.md` for human-readable version notes.
