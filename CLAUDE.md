# CLAUDE.md — Project Context for Claude Code

## What this project is

A curated **module library** of DNA-Binding Domains (DBDs), Effector Domains (EDs),
and Chromatin Remodelers (CRs) for engineered transcription factor design, supporting
MEF→iPSC reprogramming research. GitHub: https://github.com/JimYuhaoWu/eCR_mod_lib

## Current state (~10,700 records, last rebuilt 2026-05-30)

| Type | Count | Key sources |
|---|---|---|
| DBD | ~3,250 | AnimalTFDB v4 (human + mouse), UniProt sequences, JASPAR motifs |
| ED | ~7,600 | 10 sources — see SCREEN_FILES in scripts/02_seed_ed.py |
| CR | ~540 | EpiFactors v2.0 + curated YAML |

## Running the pipeline (Linux server)

```bash
python scripts/01_fetch_dbd.py        # slow (~30-60 min); skip if dbd_raw.tsv exists
python scripts/02_seed_ed.py          # all ED sources; fast
python scripts/03_fetch_cr.py         # EpiFactors; fast
python scripts/04_build_library.py    # or --rebuild for full rebuild
python scripts/05_validate.py         # QC; review library/qc_report.tsv
```

Python is at `C:\Users\86199\AppData\Local\Programs\Python\Python312\python.exe` on Windows.
Use `python` on the Linux server (confirmed Python 3.8+, all deps in requirements.txt).

## How to add a new ED source

Add an entry to `SCREEN_FILES` dict in `scripts/02_seed_ed.py`. Supported keys:

| Key | Purpose |
|---|---|
| `path` | Path to file in `data/manual/` |
| `sheet` | Excel sheet name (None for CSV) |
| `header_row` | Header row index if non-standard (default 0) |
| `col_gene`, `col_fragment` | Name fields (combined as `gene_fragment`) |
| `col_fragment_start/end` | Use coordinates as fragment name |
| `col_sequence` | AA sequence column |
| `col_score` | Quantitative activity score column |
| `col_score_fallback` | Fallback score when primary is null |
| `col_hit` + `col_hit_values` | Filter rows to hits (list, coerced to str) |
| `row_prefilter` | `{"col": "...", "value": "..."}` pre-filter |
| `min_score_threshold` | Numeric score filter |
| `score_map` | `{"H": 1.0, "M": 0.5}` categorical → numeric |
| `col_subtype` + `subtype_map` | Derive activator/repressor from a column |
| `subtype_override` | Fixed subtype for all rows |
| `col_uniprot` | Column holding UniProt ID |
| `deduplicate_on` | Deduplicate df on this column before creating records |
| `merge_file` + `merge_on` + `merge_extract` | Join a second file to add UniProt IDs etc. |
| `validation_level_override` | `screen-validated` or `ChIP-validated` |

Then copy the supplementary file to `data/manual/`, add to `literature/papers.yaml`
with `status: integrated`, and re-run the pipeline.

## Literature search workflow (every 1-3 months)

```bash
python literature/01_search.py        # PubMed + Semantic Scholar
python literature/02_triage.py        # generates literature/reviews/review_YYYY-MM-DD.md
# → review the report, update literature/papers.yaml
python literature/03_record.py        # validate decisions
```

Review reports are in `literature/reviews/`. All decisions (accepted/rejected/pending)
are tracked in `literature/papers.yaml`.

## Pending work

### Papers awaiting human review (from 2026-05 search cycle)

| DOI | Title | Why interesting |
|---|---|---|
| 10.1038/s41587-022-01624-4 | Universal DL model for ZF design (Ichikawa 2023, NatBiotech) | 49B ZF-DNA interactions; designed ZF sequences |
| 10.1016/j.cels.2023.07.001 | Combinations of activators + repressors (Mukund 2023, Cell Systems) | mmc4-2 + mmc5 in data/manual; individual domain baselines |
| 10.1016/j.cels.2023.05.008 | **Already integrated** (Ludwig 2023 viral EDs) | — |
| 10.15252/embj.2022111179 | KRAB-KAP1 structure + mapping (Stoll 2022, EMBO J) | Structural data, mutation-function table for ZNF93 KRAB |
| 10.1016/j.molcel.2023.02.011 | Massively parallel CRISPRa in iPSC (Wu 2023, Mol Cell) | Tests existing effectors; no new domain characterization |
| 10.1016/j.gene.2017.10.058 | CHD8short truncated remodeler (Kunkel 2018, Gene) | **Zebrafish** — likely reject |
| 10.3791/64403 | In vitro selection of engineered repressors (JoVE 2023) | Protocol paper — likely reject |
| 10.1016/j.jbc.2021.100947 | ZBTB2 recruits chromatin remodelers (Olivieri 2021, JBC) | Mechanistic — likely reject |

### Data not yet integrated (files in data/manual/)

| File | Content | Notes |
|---|---|---|
| `TENet_2024_media-5.xlsx` | 220 TENet-designed repressor variants | Scored low at test loci — defer until published |
| `Mukund_2023_mmc5.xlsx` | 4,168 Pfam nuclear domain screen | Overlaps HiTEff heavily; check for unique domains |
| `Mukund_2023_mmc4.xlsx` | 8,400 ED pair combinations | Combinatorial data, not individual domain characterization |

### TENet preprint (10.1101/2024.09.21.614253)
- Currently integrated as `status: integrated` with 51 WT RD sequences
- **Update DOI to published journal DOI** when paper is peer-reviewed
- The designed sequences (media-5) could be integrated once published

## Known issues / QC flags

- `dbd_no_identifier` (~421 DBDs): no UniProt or JASPAR hit — real data gaps in AnimalTFDB
- The `python` command on Windows resolves to Windows Store stub; use full path or fix PATH:
  `C:\Users\86199\AppData\Local\Programs\Python\Python312\python.exe`
- `library/build_manifest.json` conflict: always prefer server version (authoritative rebuild)

## Key people

- **Josh Tycko** (Stanford/Bassik lab): multiple integrated sources — HiTEff (2020),
  compact effectors (2025), TENet DMS (2024 preprint). Watch his lab for new papers.
- **Lacramioara Bintu** (Stanford): corresponding author on Tycko 2025, Ludwig 2023,
  Mukund 2023. Same lab, same HT-recruit assay platform.
