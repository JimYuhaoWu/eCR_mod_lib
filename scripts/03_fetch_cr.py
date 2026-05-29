"""
scripts/03_fetch_cr.py
----------------------
Populate Chromatin Remodeler (CR) entries from two sources:

  1. data/manual/cr_curated.yaml
     – Key CRs with full annotation; sequences fetched from UniProt.

  2. EpiFactors 2.0 (https://epifactors.autosome.org)
     – Comprehensive list of epigenetic regulators.
     – Entries with an enzymatic activity that remodels chromatin
       (HAT, HDAC, HMT, HDM, DNMTs, ATP-dependent remodelers) are
       extracted; others are skipped.
     – EpiFactors may require manual download (see notes in
       data/manual/README.md if the direct download URL is unavailable).

Outputs
-------
  data/raw/epifactors/epifactors.tsv    raw EpiFactors download
  data/raw/uniprot/*.fasta              cached sequences
  data/processed/cr_raw.tsv            merged CR records

Run
---
  python scripts/03_fetch_cr.py [--config config/sources.yaml] [--dry-run]
"""

from __future__ import annotations

import argparse
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests
import yaml
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))
from utils import (
    get_logger, load_manifest, save_manifest,
    record_download, write_build_manifest,
)

ROOT = Path(__file__).resolve().parents[1]
LOGS_DIR = ROOT / "logs"
RAW_DIR = ROOT / "data" / "raw"
PROCESSED_DIR = ROOT / "data" / "processed"
MANUAL_DIR = ROOT / "data" / "manual"

# EpiFactors v2.0 Function values to include as CRs
_EPIFACTORS_INCLUDE = {
    "Histone modification write",
    "Histone modification erase",
    "Chromatin remodeling",
    "DNA modification",
}

# Map EpiFactors Function → chromatin_state_effect (substring match)
_EFFECT_MAP = {
    "Histone modification write": "variable / histone mark deposition (HAT/HMT/ubiquitin ligase)",
    "Histone modification erase": "variable / histone mark removal (HDAC/HDM/deubiquitinase)",
    "Chromatin remodeling": "open or close / ATP-dependent nucleosome positioning",
    "DNA modification": "variable / DNA methylation or demethylation",
}


# ── Helper ────────────────────────────────────────────────────────────────────

def _get(url, session, timeout=20):
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


def fetch_uniprot_fasta(uniprot_id: str, session, log) -> str | None:
    """Return full sequence from UniProt; cache locally."""
    cache = RAW_DIR / "uniprot" / f"{uniprot_id}.fasta"
    if cache.exists():
        text = cache.read_text()
    else:
        resp = _get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta", session)
        if resp is None:
            log.warning(f"[UniProt] could not fetch {uniprot_id}")
            return None
        cache.write_text(resp.text, encoding="utf-8")
        text = resp.text
    lines = text.strip().splitlines()
    seq = "".join(l.strip() for l in lines if not l.startswith(">"))
    return seq or None


def slice_domain(full_seq: str | None, domain_residues: str | None) -> str | None:
    """Slice full_seq to domain_residues range ('start-end', 1-indexed, inclusive)."""
    if full_seq is None:
        return None
    if not domain_residues:
        return full_seq
    try:
        start, end = (int(x) for x in domain_residues.split("-"))
        return full_seq[start - 1: end]
    except Exception:
        return full_seq


# ── Curated YAML ──────────────────────────────────────────────────────────────

def process_curated_crs(yaml_path: Path, session, log, dry_run: bool) -> list[dict]:
    with open(yaml_path, encoding="utf-8") as f:
        raw = yaml.safe_load(f)

    records = []
    now = datetime.now(timezone.utc).isoformat()

    for entry in tqdm(raw.get("remodelers", []), desc="Curated CRs", ncols=80):
        name = entry["name"]
        uniprot_id = entry.get("uniprot_id") or None
        domain_residues = entry.get("domain_residues") or None

        full_seq = None
        if uniprot_id and not dry_run:
            full_seq = fetch_uniprot_fasta(uniprot_id, session, log)
        seq = slice_domain(full_seq, domain_residues)

        module_id = f"CR_{name.upper().replace('-','_').replace(' ','_')}"
        record = {
            "module_id": module_id,
            "type": "CR",
            "subtype": entry.get("subtype", ""),
            "name": name,
            "organism": entry.get("organism", ""),
            "source_species": entry.get("organism", ""),
            "gene_symbol": entry.get("source_protein", "").split("(")[-1].rstrip(")") if "(" in entry.get("source_protein", "") else None,
            "uniprot_id": uniprot_id,
            "sequence_aa": seq,
            "length_aa": len(seq) if seq else None,
            "target_or_mechanism": entry.get("mechanism", ""),
            "validation_level": entry.get("validation_level", "ChIP-validated"),
            "source": f"Curated: {entry.get('doi', '')}",
            "source_doi": entry.get("doi"),
            "source_version": "manual_v1",
            "source_date": now[:10],
            "known_interactors": entry.get("known_interactors", ""),
            "linker_notes": entry.get("linker_notes", ""),
            "engineering_compatibility": entry.get("engineering_compatibility", ""),
            "chromatin_state_effect": entry.get("chromatin_state_effect", ""),
            "notes": (f"complex: {entry.get('complex_membership','')}  "
                      f"target: {entry.get('target_modification','')}  "
                      f"{entry.get('notes','')}").strip(),
            "date_added": now,
            "date_modified": now,
        }
        records.append(record)
        log.debug(f"[curated CR] {name}  seq={'yes' if seq else 'no'}  len={len(seq) if seq else 0}")

    return records


# ── EpiFactors download + parse ───────────────────────────────────────────────

EPIFACTORS_URL = "https://epifactors.autosome.org/public_data/v2.0.zip"
EPIFACTORS_CACHE = RAW_DIR / "epifactors" / "EpiGenes_main.csv"
EPIFACTORS_MANUAL = ROOT / "data" / "manual" / "EpiGenes_main.csv"


def fetch_epifactors(session, log, manifest: dict, dry_run: bool) -> pd.DataFrame | None:
    """Load EpiFactors EpiGenes_main.csv — from manual/, cache, or zip download."""
    # Prefer manually placed file
    for src in (EPIFACTORS_MANUAL, EPIFACTORS_CACHE):
        if src.exists():
            log.info(f"[EpiFactors] using {src}")
            try:
                return pd.read_csv(src, dtype=str)
            except Exception as e:
                log.warning(f"[EpiFactors] parse error for {src}: {e}")

    if dry_run:
        log.info("[EpiFactors] DRY-RUN: skipping download")
        return None

    log.info(f"[EpiFactors] downloading zip from {EPIFACTORS_URL}")
    resp = _get(EPIFACTORS_URL, session, timeout=60)
    if resp is None:
        log.warning(
            "[EpiFactors] automatic download failed. "
            "Download https://epifactors.autosome.org/public_data/v2.0.zip, "
            "extract EpiGenes_main.csv and place at data/manual/EpiGenes_main.csv"
        )
        return None

    import io, zipfile
    EPIFACTORS_CACHE.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(resp.content)) as z:
        with z.open("v2.0/EpiGenes_main.csv") as f:
            content = f.read()
    EPIFACTORS_CACHE.write_bytes(content)
    record_download(manifest, "epifactors", EPIFACTORS_URL, EPIFACTORS_CACHE,
                    source_version="2.0", notes="EpiFactors 2.0 EpiGenes_main.csv")
    save_manifest(manifest, RAW_DIR)

    try:
        return pd.read_csv(EPIFACTORS_CACHE, dtype=str)
    except Exception as e:
        log.error(f"[EpiFactors] parse error after download: {e}")
        return None


def process_epifactors(df: pd.DataFrame, session, log, dry_run: bool,
                        curated_ids: set[str]) -> list[dict]:
    """
    Filter EpiFactors to enzymatic CRs; enrich with UniProt sequences.
    Skips entries already present in curated_ids (avoids duplicates).
    """
    now = datetime.now(timezone.utc).isoformat()
    records = []

    df.columns = [c.strip() for c in df.columns]

    # EpiFactors v2.0 column names
    col_gene    = "HGNC_symbol"
    col_uniprot = "UniProt_AC"
    col_func    = "Function"
    col_mod     = "Modification"

    if col_gene not in df.columns or col_uniprot not in df.columns:
        log.warning(f"[EpiFactors] unexpected columns: {list(df.columns)[:10]}. "
                    f"Skipping EpiFactors integration.")
        return []

    log.info(f"[EpiFactors] {len(df):,} total entries; filtering to enzymatic CRs")

    for _, row in tqdm(df.iterrows(), total=len(df), desc="EpiFactors", ncols=80):
        gene = str(row.get(col_gene, "")).strip()
        uniprot_id = str(row.get(col_uniprot, "")).strip()
        function = str(row.get(col_func, "")).strip() if col_func in df.columns else ""

        if not any(inc.lower() in function.lower() for inc in _EPIFACTORS_INCLUDE):
            continue
        if not gene or gene.lower() == "nan":
            continue

        module_id = f"CR_EPIFACTORS_{gene.upper()}"
        # Skip if already covered by curated entries
        if module_id in curated_ids:
            continue

        seq = None
        if uniprot_id and uniprot_id.lower() not in ("nan", "") and not dry_run:
            seq = fetch_uniprot_fasta(uniprot_id, session, log)

        # Map function to effect
        effect = next((v for k, v in _EFFECT_MAP.items() if k.lower() in function.lower()), "")

        record = {
            "module_id": module_id,
            "type": "CR",
            "subtype": function,
            "name": gene,
            "organism": "Homo sapiens",
            "source_species": "Homo sapiens",
            "gene_symbol": gene,
            "uniprot_id": uniprot_id if uniprot_id.lower() not in ("nan","") else None,
            "sequence_aa": seq,
            "length_aa": len(seq) if seq else None,
            "target_or_mechanism": function,
            "validation_level": "predicted",
            "source": "EpiFactors 2.0",
            "source_doi": "10.1093/nar/gkab1193",
            "source_version": "2.0",
            "source_date": now[:10],
            "chromatin_state_effect": effect,
            "notes": f"EpiFactors 2.0 functional annotation: {function}",
            "date_added": now,
            "date_modified": now,
        }
        records.append(record)

    return records


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="config/sources.yaml")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    log = get_logger("03_fetch_cr", LOGS_DIR)
    log.info("=== 03_fetch_cr.py started ===")

    session = requests.Session()
    session.headers.update({"User-Agent": "module_library/1.0 (research)"})

    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    for subdir in ("epifactors",):
        (RAW_DIR / subdir).mkdir(parents=True, exist_ok=True)

    manifest = load_manifest(RAW_DIR)

    # 1. Curated CRs
    curated = process_curated_crs(MANUAL_DIR / "cr_curated.yaml", session, log, args.dry_run)
    curated_ids = {r["module_id"] for r in curated}
    log.info(f"Curated CRs: {len(curated)}")

    # 2. EpiFactors
    epif_df = fetch_epifactors(session, log, manifest, args.dry_run)
    epif_records = []
    if epif_df is not None:
        epif_records = process_epifactors(epif_df, session, log, args.dry_run, curated_ids)
        log.info(f"EpiFactors CRs: {len(epif_records)}")

    all_records = curated + epif_records
    df = pd.DataFrame(all_records)
    out = PROCESSED_DIR / "cr_raw.tsv"
    df.to_csv(out, sep="\t", index=False)
    log.info(f"Wrote {len(df):,} CR records to {out}")

    write_build_manifest(
        ROOT / "library",
        script_name="03_fetch_cr.py",
        sources_used=[
            {"name": "cr_curated.yaml", "version": "manual_v1"},
            {"name": "EpiFactors", "version": "2.0", "doi": "10.1093/nar/gkab1193"},
        ],
        records_added=len(all_records),
    )

    log.info("=== 03_fetch_cr.py complete ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
