"""
scripts/02_seed_ed.py
---------------------
Populate Effector Domain (ED) entries in the module library from two sources:

  1. data/manual/ed_curated.yaml
     – Manually curated canonical EDs (VP64, VPR, KRAB, p65-TAD, etc.)
     – Sequences fetched from UniProt for entries with a uniprot_id;
       sequence_manual used directly for synthetic constructs.

  2. Published activation-domain screens (Alerasool et al. 2022, DelRosso et al. 2023)
     – These require manual download of supplementary tables from the journal
       website; the script checks for them and processes them if present.
     – Instructions for manual download are printed if files are missing.

Outputs
-------
  data/raw/uniprot/{id}.fasta         cached sequences
  data/processed/ed_raw.tsv           merged ED records
  library/build_manifest.json         updated build record

Run
---
  python scripts/02_seed_ed.py [--config config/sources.yaml] [--dry-run]
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests
import yaml
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))
from utils import get_logger, write_build_manifest

ROOT = Path(__file__).resolve().parents[1]
LOGS_DIR = ROOT / "logs"
RAW_DIR = ROOT / "data" / "raw"
PROCESSED_DIR = ROOT / "data" / "processed"
MANUAL_DIR = ROOT / "data" / "manual"

# Mapping: curated subtype string → canonical group for chromatin_state_effect
_EFFECT_MAP = {
    "activator": "open / H3K27ac deposition",
    "repressor": "close / H3K9me2-3 or H3K27me3",
    "dual": "context-dependent",
}


# ── UniProt sequence fetch (shared with 01_fetch_dbd.py) ─────────────────────

def _get(url: str, session: requests.Session, timeout: int = 20) -> requests.Response | None:
    import time
    for attempt in range(3):
        try:
            r = session.get(url, timeout=timeout)
            if r.status_code == 200:
                return r
            if r.status_code == 404:
                return None
        except requests.RequestException:
            pass
        time.sleep(2 ** attempt)
    return None


def fetch_domain_sequence(uniprot_id: str, domain_residues: str | None,
                           session: requests.Session, log) -> str | None:
    """
    Fetch the full protein sequence from UniProt, then slice to domain_residues
    if provided (e.g. '411-490').  Returns the domain AA sequence or None.
    """
    if not uniprot_id:
        return None

    cache = RAW_DIR / "uniprot" / f"{uniprot_id}.fasta"
    if cache.exists():
        fasta = cache.read_text()
    else:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
        resp = _get(url, session, timeout=20)
        if resp is None:
            log.warning(f"[UniProt] could not fetch {uniprot_id}")
            return None
        cache.write_text(resp.text, encoding="utf-8")
        fasta = resp.text

    # Parse FASTA
    lines = fasta.strip().splitlines()
    full_seq = "".join(l.strip() for l in lines if not l.startswith(">"))
    if not full_seq:
        return None

    # Slice domain
    if domain_residues:
        try:
            start, end = (int(x) for x in domain_residues.split("-"))
            return full_seq[start - 1 : end]   # 1-indexed, inclusive
        except Exception:
            log.warning(f"[UniProt] could not parse domain_residues='{domain_residues}'")
    return full_seq


# ── Process curated YAML ──────────────────────────────────────────────────────

def process_curated_eds(cfg_path: Path, session: requests.Session,
                         log, dry_run: bool) -> list[dict]:
    with open(cfg_path, encoding="utf-8") as f:
        raw = yaml.safe_load(f)

    records = []
    now = datetime.now(timezone.utc).isoformat()

    for category in ("activators", "repressors"):
        entries = raw.get(category, [])
        for entry in tqdm(entries, desc=f"Curated {category}", ncols=80):
            name = entry["name"]
            subtype = entry.get("subtype", "activator")
            uniprot_id = entry.get("uniprot_id") or None
            domain_residues = entry.get("domain_residues") or None
            seq_manual = entry.get("sequence_manual") or None

            # Sequence resolution order: manual > UniProt slice
            seq = seq_manual
            if seq is None and uniprot_id and not dry_run:
                seq = fetch_domain_sequence(uniprot_id, domain_residues, session, log)
            length = len(seq) if seq else None

            module_id = f"ED_{name.upper().replace('-','_').replace(' ','_')}"

            record = {
                "module_id": module_id,
                "type": "ED",
                "subtype": subtype,
                "name": name,
                "organism": entry.get("organism", ""),
                "source_species": entry.get("organism", ""),
                "gene_symbol": None,
                "uniprot_id": uniprot_id,
                "sequence_aa": seq,
                "length_aa": length,
                "target_or_mechanism": entry.get("mechanism", ""),
                "quantitative_metric": entry.get("activity_score"),
                "quantitative_metric_label": entry.get("activity_metric"),
                "quantitative_metric_source": entry.get("activity_source"),
                "validation_level": entry.get("validation_level", "screen-validated"),
                "source": "Curated: " + entry.get("doi", ""),
                "source_doi": entry.get("doi"),
                "source_version": "manual_v1",
                "source_date": now[:10],
                "known_interactors": entry.get("known_interactors", ""),
                "linker_notes": entry.get("linker_notes", ""),
                "engineering_compatibility": entry.get("engineering_compatibility", ""),
                "chromatin_state_effect": _EFFECT_MAP.get(subtype, ""),
                "notes": str(entry.get("notes", "")).strip(),
                "date_added": now,
                "date_modified": now,
            }
            records.append(record)
            log.debug(f"[curated ED] {name}  seq={'yes' if seq else 'no'}  len={length}")

    return records


# ── Process Alerasool / DelRosso screen data ──────────────────────────────────

SCREEN_FILES = {
    "Alerasool_2022": {
        "path": MANUAL_DIR / "Alerasool_2022_SupTable.xlsx",
        "doi": "10.1016/j.molcel.2021.12.005",
        "sheet": "tAD-seq",
        "col_name": None,           # built from Gene + Fragment below
        "col_gene": "Gene",
        "col_fragment": "Fragment",
        "col_sequence": "Sequence",
        "col_score": "logFC high GFP",
        "col_hit": "Hit high GFP",
        "hit_value": "hit",
        "notes": (
            "Download mmc2.xlsx (Supplementary Table 2) from "
            "https://doi.org/10.1016/j.molcel.2021.12.005 "
            "and save as data/manual/Alerasool_2022_SupTable.xlsx"
        ),
    },
    "DelRosso_2023_activators": {
        "path": MANUAL_DIR / "DelRosso_2023_SupTable.xlsx",
        "doi": "10.1038/s41586-023-05906-y",
        "sheet": "Activation Domains",
        "col_name": None,
        "col_gene": "HGNC symbol",
        "col_fragment": "Domain",
        "col_sequence": "Sequence",
        "col_score": "Max avg act",
        "col_hit": None,
        "hit_value": None,
        "subtype_override": "activator",
        "notes": (
            "Download Supplementary Table 6 from "
            "https://doi.org/10.1038/s41586-023-05906-y "
            "and save as data/manual/DelRosso_2023_SupTable.xlsx"
        ),
    },
    "DelRosso_2023_repressors": {
        "path": MANUAL_DIR / "DelRosso_2023_SupTable.xlsx",
        "doi": "10.1038/s41586-023-05906-y",
        "sheet": "Repression Domains",
        "col_name": None,
        "col_gene": "HGNC symbol",
        "col_fragment": "Domain",
        "col_sequence": "Sequence",
        "col_score": "Max avg pEF",
        "col_score_fallback": "Max avg PGK",
        "col_hit": None,
        "hit_value": None,
        "subtype_override": "repressor",
        "notes": (
            "Download Supplementary Table 6 from "
            "https://doi.org/10.1038/s41586-023-05906-y "
            "and save as data/manual/DelRosso_2023_SupTable.xlsx"
        ),
    },
    # Compendium of human TF effector domains — experimentally characterised ADs/RDs.
    # Filter: Activity H or M (drop L, NaN, S). Validation: ChIP-validated.
    "Compendium_2021_activators": {
        "path": MANUAL_DIR / "Compendium_2021_SupTable.xlsx",
        "doi": "10.1016/j.molcel.2021.11.007",
        "sheet": "Table S2",
        "col_name": None,
        "col_gene": "TF name",
        "col_fragment": "Effector domain ID",
        "col_sequence": "Sequence",
        "col_score": "Activity (H, M or L)",
        "score_map": {"H": 1.0, "M": 0.5},
        "col_hit": "Activity (H, M or L)",
        "col_hit_values": ["H", "M"],
        "row_prefilter": {"col": "Domain type", "value": "AD"},
        "subtype_override": "activator",
        "validation_level_override": "ChIP-validated",
        "notes": (
            "Download Supplementary Table 2 (mmc8.xlsx) from "
            "https://doi.org/10.1016/j.molcel.2021.11.007 "
            "and save as data/manual/Compendium_2021_SupTable.xlsx"
        ),
    },
    "Compendium_2021_repressors": {
        "path": MANUAL_DIR / "Compendium_2021_SupTable.xlsx",
        "doi": "10.1016/j.molcel.2021.11.007",
        "sheet": "Table S2",
        "col_name": None,
        "col_gene": "TF name",
        "col_fragment": "Effector domain ID",
        "col_sequence": "Sequence",
        "col_score": "Activity (H, M or L)",
        "score_map": {"H": 1.0, "M": 0.5},
        "col_hit": "Activity (H, M or L)",
        "col_hit_values": ["H", "M"],
        "row_prefilter": {"col": "Domain type", "value": "RD"},
        "subtype_override": "repressor",
        "validation_level_override": "ChIP-validated",
        "notes": (
            "Download Supplementary Table 2 (mmc8.xlsx) from "
            "https://doi.org/10.1016/j.molcel.2021.11.007 "
            "and save as data/manual/Compendium_2021_SupTable.xlsx"
        ),
    },
    # High-throughput nuclear domain screens (Pfam-domain library).
    # DOI: TODO — confirm from paper. Validation: screen-validated.
    "HiTEff_activators": {
        "path": MANUAL_DIR / "HiTEff_SupTable.xlsx",
        "doi": "10.1016/j.cell.2020.11.024",
        "sheet": "NucAct_data",
        "col_name": None,
        "col_gene": "Gene entry name",
        "col_fragment": "Domain ID",
        "col_sequence": "Extended Domain sequence",
        "col_score": "Avg Act",
        "col_hit": "Hit",
        "col_hit_values": [True],
        "subtype_override": "activator",
        "notes": (
            "Download mmc4.xlsx from https://doi.org/10.1016/j.cell.2020.11.024 "
            "and save as data/manual/HiTEff_SupTable.xlsx"
        ),
    },
    "HiTEff_repressors": {
        "path": MANUAL_DIR / "HiTEff_SupTable.xlsx",
        "doi": "10.1016/j.cell.2020.11.024",
        "sheet": "NucRepr_data",
        "col_name": None,
        "col_gene": "Gene entry name",
        "col_fragment": "Domain ID",
        "col_sequence": "Extended Domain sequence",
        "col_score": "Avg ReprD13",
        "col_hit": "Hit",
        "col_hit_values": [True],
        "subtype_override": "repressor",
        "notes": (
            "Download mmc4.xlsx from https://doi.org/10.1016/j.cell.2020.11.024 "
            "and save as data/manual/HiTEff_SupTable.xlsx"
        ),
    },
}


def process_screen_data(log) -> list[dict]:
    """
    Parse available screen supplementary tables. Prints instructions for
    missing files; skips silently if absent.
    """
    records = []
    now = datetime.now(timezone.utc).isoformat()

    for screen_name, info in SCREEN_FILES.items():
        if not info["path"].exists():
            log.warning(
                f"[{screen_name}] file not found: {info['path']}\n"
                f"  → {info['notes']}"
            )
            continue

        log.info(f"[{screen_name}] parsing {info['path']}")
        try:
            sheet = info.get("sheet")
            if info["path"].suffix == ".xlsx":
                df = pd.read_excel(info["path"], sheet_name=sheet)
            else:
                df = pd.read_csv(info["path"], sep="\t")
        except Exception as e:
            log.error(f"[{screen_name}] parse error: {e}")
            continue

        # Apply row pre-filter (e.g. restrict to a specific domain type)
        pf = info.get("row_prefilter")
        if pf and pf["col"] in df.columns:
            df = df[df[pf["col"]].astype(str).str.strip() == str(pf["value"])]

        # Filter to hits — supports single value (str) or list of values
        col_hit = info.get("col_hit")
        hit_val = info.get("hit_value")
        hit_vals = info.get("col_hit_values")
        if col_hit and col_hit in df.columns:
            if hit_vals is not None:
                df = df[df[col_hit].isin(hit_vals)]
            elif hit_val:
                df = df[df[col_hit].astype(str).str.strip().str.lower() == hit_val.lower()]

        col_n = info["col_name"]
        col_s = info.get("col_sequence")
        col_q = info.get("col_score")
        col_gene = info.get("col_gene")
        col_frag = info.get("col_fragment")

        for _, row in df.iterrows():
            # Build name from gene+fragment if col_name is None
            if col_n is None and col_gene and col_frag:
                name = f"{str(row.get(col_gene, '')).strip()}_{str(row.get(col_frag, '')).strip()}"
            else:
                name = str(row.get(col_n, "")).strip()
            if not name or name.lower() in ("nan", "_", "nan_nan"):
                continue

            seq = str(row.get(col_s, "")).strip() if col_s and col_s in df.columns else None
            if seq in (None, "nan", ""):
                seq = None
            def _parse_score(val):
                try:
                    f = float(val)
                    return None if pd.isna(f) else f
                except (ValueError, TypeError):
                    return None

            score = None
            col_score_used = col_q
            score_map = info.get("score_map")
            if col_q and col_q in df.columns:
                raw_val = row[col_q]
                if score_map and str(raw_val).strip() in score_map:
                    score = score_map[str(raw_val).strip()]
                else:
                    score = _parse_score(raw_val)
            # Fall back to secondary score column if primary is null
            if score is None:
                col_fb = info.get("col_score_fallback")
                if col_fb and col_fb in df.columns:
                    fb_score = _parse_score(row[col_fb])
                    if fb_score is not None:
                        score = fb_score
                        col_score_used = col_fb

            if info.get("subtype_override"):
                subtype = info["subtype_override"]
            else:
                is_repressor = "repression" in screen_name.lower() or "repression" in (col_q or "")
                subtype = "repressor" if is_repressor else "activator"

            module_id = f"ED_{screen_name.upper()}_{name.upper().replace(' ','_')[:40]}"
            record = {
                "module_id": module_id,
                "type": "ED",
                "subtype": subtype,
                "name": name,
                "organism": "Homo sapiens",
                "sequence_aa": seq,
                "length_aa": len(seq) if seq else None,
                "quantitative_metric": score,
                "quantitative_metric_label": col_score_used,
                "quantitative_metric_source": info["doi"],
                "validation_level": info.get("validation_level_override", "screen-validated"),
                "source": screen_name,
                "source_doi": info["doi"],
                "source_version": "supplementary_table",
                "source_date": now[:10],
                "chromatin_state_effect": _EFFECT_MAP.get(subtype, ""),
                "notes": f"From {screen_name}; {info['doi']}",
                "date_added": now,
                "date_modified": now,
            }
            records.append(record)

        log.info(f"[{screen_name}] loaded {len([r for r in records if screen_name in r['source']]):,} records")

    return records


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="config/sources.yaml")
    parser.add_argument("--dry-run", action="store_true",
                        help="Skip network calls; parse local files only")
    args = parser.parse_args()

    log = get_logger("02_seed_ed", LOGS_DIR)
    log.info("=== 02_seed_ed.py started ===")

    session = requests.Session()
    session.headers.update({"User-Agent": "module_library/1.0 (research)"})

    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Canonical curated EDs
    curated = process_curated_eds(MANUAL_DIR / "ed_curated.yaml", session, log, args.dry_run)
    log.info(f"Curated EDs: {len(curated)} entries")

    # 2. Screen data (skipped if files missing)
    screen = process_screen_data(log)
    log.info(f"Screen EDs: {len(screen)} entries")

    all_records = curated + screen
    df = pd.DataFrame(all_records)
    out = PROCESSED_DIR / "ed_raw.tsv"
    df.to_csv(out, sep="\t", index=False)
    log.info(f"Wrote {len(df):,} ED records to {out}")

    write_build_manifest(
        ROOT / "library",
        script_name="02_seed_ed.py",
        sources_used=[
            {"name": "ed_curated.yaml", "version": "manual_v1"},
            {"name": "Alerasool_2022", "doi": SCREEN_FILES["Alerasool_2022"]["doi"]},
            {"name": "DelRosso_2023",  "doi": SCREEN_FILES["DelRosso_2023_activators"]["doi"]},
        ],
        records_added=len(all_records),
    )

    log.info("=== 02_seed_ed.py complete ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
