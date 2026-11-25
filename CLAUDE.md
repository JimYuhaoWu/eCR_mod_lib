# CLAUDE.md — Project Context for Claude Code

## Coding principles

- **No features beyond what was asked.** No speculative abstractions, configurability, or error handling for impossible scenarios.
- **Surgical changes.** Touch only what the task requires. Don't improve adjacent code, comments, or formatting. Match existing style.
- **Surface tradeoffs before coding.** If multiple interpretations exist, present them — don't pick silently. If something is unclear, ask.
- **If you notice unrelated dead code or issues, mention them — don't silently fix them.**

## What this project is

A curated **module library** of DNA-Binding Domains (DBDs), Effector Domains (EDs),
and Chromatin Remodelers (CRs) for engineered transcription factor design, supporting
MEF→iPSC reprogramming research. GitHub: https://github.com/JimYuhaoWu/eCR_mod_lib

## Current state (11,486 records, last rebuilt 2026-11-25)

| Type | Count | Key sources |
|---|---|---|
| DBD | ~3,250 | AnimalTFDB v4 (human + mouse), UniProt sequences + zinc-finger counts, JASPAR motifs |
| ED | 7,701 | 11 sources — see SCREEN_FILES in scripts/02_seed_ed.py. `02` writes 7,701 and all of them reach the library since F10; colliding ids carry a content digest |
| CR | ~540 | EpiFactors v2.0 + curated YAML |

## Environment & installation

Shared conda environment (`ecr`) covers both repos:

```bash
conda env create -f environment.yml   # once per machine
conda activate ecr
pip install -e /path/to/eCR_mod_lib
```

### Public API — `load()` (zero-setup)

```python
from scripts import load

df = load()                                   # all records from library/module_library.tsv
eds = load(type="ED")                         # 'DBD', 'ED', or 'CR'
human_dbds = load(type="DBD", organism="Homo sapiens")
```

`load()` reads the committed TSV snapshot — no pipeline run needed. Implemented in
`scripts/__init__.py`. Raises `FileNotFoundError` with pipeline instructions if the
TSV is missing.

### Lower-level API — `ModuleLibrary` (pipeline use)

```python
from scripts.schema import ModuleLibrary

with ModuleLibrary("library/module_library.db") as lib:
    lib.insert_or_replace(records)
    df = lib.to_dataframe("DBD")
```

Used internally by the pipeline scripts and as a dependency by
**[eCR_predictor](https://github.com/JimYuhaoWu/eCR_predictor)**.

## Running the pipeline (Linux server)

```bash
python scripts/01_fetch_dbd.py        # slow (~30-60 min); skip if dbd_raw.tsv exists
python scripts/02_seed_ed.py          # all ED sources; fast
python scripts/03_fetch_cr.py         # EpiFactors; fast
python scripts/04_build_library.py    # or --rebuild for full rebuild
python scripts/05_validate.py         # QC; review library/qc_report.tsv
```

Use `python` on both Windows (Python 3.12.4, PATH configured) and the Linux server (confirmed Python 3.8+, all deps in requirements.txt).

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
| `exclude_labels_from` | List of `{path, sheet, col}` dicts — rows whose label matches any entry are dropped (used to avoid re-importing entries already present in another source) |
| `col_label` | Column in the current file to match against `exclude_labels_from` references |

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

### 2026-05 literature cycle — complete
All 17 candidates reviewed. One paper integrated:
- **Mukund 2023** (10.1016/j.cels.2023.07.001): **55 rows per readout, all 110 in the library.** The dual-readout screen is ingested through two SCREEN_FILES entries, so each domain gets one `activator` record (Avg Act) and one `repressor` record (Avg ReprD13). Since F8 each carries `readout` (activation/repression) and `paired_module_id` pointing at its sibling, so the opposite `subtype` reads as two measurements rather than a contradiction. `exclude_labels_from` filters the 4,113 HiTEff duplicates automatically. ✓ integrated 2026-05-30.

### zinc_finger_count column — complete
`scripts/01_fetch_dbd.py` fetches a UniProt-annotated zinc-finger count
(`ZN_FING` features) for C2H2 families (`zf-C2H2`, `ZBTB`) into the
`zinc_finger_count` column. It has been populated in the committed TSV since
the 2026-09-04 server rebuild (`d7c38b4`): **1,010 of the 3,248 DBD records
carry a count** — the C2H2/ZBTB families the fetch targets. ✓ done 2026-09-04.

Since F3 (2026-11-18) there is a second column beside it,
**`dbd_zinc_finger_count`** — the fingers the *extracted span* actually kept,
populated on 2,439 records. `zinc_finger_count` counts every `ZN_FING` feature
over the whole protein, PHD and RING fingers included, so it is an upper bound:
ZFPM2 claims 8 for a 102-aa span of a 1,151-aa protein. eCR_predictor's
`af3.py` prefers `dbd_zinc_finger_count` when adding Zn ions to AF3 inputs and
falls back to `zinc_finger_count` with a warning (correct for `refine.py`,
which folds full-length sequences).

### Data not yet integrated (files in data/manual/)

| File | Content | Notes |
|---|---|---|
| `TENet_2024_media-5.xlsx` | 220 TENet-designed repressor variants | Scored low at test loci — defer until published |
| `Mukund_2023_mmc4.xlsx` | 8,400 ED pair combinations | Combinatorial data — not integrated |

### TENet preprint (10.1101/2024.09.21.614253)
- Currently integrated as `status: integrated` with 51 WT RD sequences
- **Update DOI to published journal DOI** when paper is peer-reviewed
- The designed sequences (media-5) could be integrated once published

## Known issues / QC flags

- `dbd_no_identifier` (~421 DBDs): no UniProt or JASPAR hit — real data gaps in AnimalTFDB
- `python` on Windows now resolves to Python 3.12.4 (PATH configured, App Execution Aliases disabled)
- `library/build_manifest.json` conflict: always prefer server version (authoritative rebuild)

## Key people

- **Josh Tycko** (Stanford/Bassik lab): multiple integrated sources — HiTEff (2020),
  compact effectors (2025), TENet DMS (2024 preprint). Watch his lab for new papers.
- **Lacramioara Bintu** (Stanford): corresponding author on Tycko 2025, Ludwig 2023,
  Mukund 2023. Same lab, same HT-recruit assay platform.
