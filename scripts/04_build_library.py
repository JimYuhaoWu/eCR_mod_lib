"""
scripts/04_build_library.py
---------------------------
Assemble the three processed TSV files into the SQLite module library
and export a TSV snapshot for git tracking.

Inputs
------
  data/processed/dbd_raw.tsv
  data/processed/ed_raw.tsv
  data/processed/cr_raw.tsv

Outputs
-------
  library/module_library.db      SQLite database
  library/module_library.tsv     full TSV export (commit to git)
  library/build_manifest.json    updated build provenance

Run
---
  python scripts/04_build_library.py [--rebuild]

  --rebuild  : drop and recreate the database from scratch
               (use when changing schema or doing a clean build)
"""

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from schema import ModuleLibrary
from utils import get_logger, write_build_manifest

ROOT = Path(__file__).resolve().parents[1]
LOGS_DIR = ROOT / "logs"
PROCESSED_DIR = ROOT / "data" / "processed"
LIBRARY_DIR = ROOT / "library"

SOURCES = [
    ("DBD", PROCESSED_DIR / "dbd_raw.tsv"),
    ("ED",  PROCESSED_DIR / "ed_raw.tsv"),
    ("CR",  PROCESSED_DIR / "cr_raw.tsv"),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rebuild", action="store_true",
                        help="Drop and recreate the database before importing")
    args = parser.parse_args()

    log = get_logger("04_build_library", LOGS_DIR)
    log.info("=== 04_build_library.py started ===")

    LIBRARY_DIR.mkdir(parents=True, exist_ok=True)
    db_path = LIBRARY_DIR / "module_library.db"

    if args.rebuild and db_path.exists():
        db_path.unlink()
        log.info("[rebuild] deleted existing database")

    total_inserted = total_updated = 0
    sources_meta = []

    with ModuleLibrary(db_path) as lib:

        for module_type, tsv_path in SOURCES:
            if not tsv_path.exists():
                log.warning(f"[{module_type}] {tsv_path} not found — run script "
                            f"0{SOURCES.index((module_type,tsv_path))+1} first")
                continue

            df = pd.read_csv(tsv_path, sep="\t", dtype=str)
            df = df.where(df.notna(), None)   # convert NaN → None for SQLite NULL

            # Ensure required timestamps exist
            now = datetime.now(timezone.utc).isoformat()
            for col in ("date_added", "date_modified"):
                if col not in df.columns:
                    df[col] = now
                else:
                    df[col] = df[col].fillna(now)

            records = df.to_dict(orient="records")
            # Drop rows without module_id
            records = [r for r in records if r.get("module_id")]

            inserted, updated = lib.insert_or_replace(records)
            total_inserted += inserted
            total_updated += updated
            log.info(f"[{module_type}] {inserted} inserted  {updated} updated  "
                     f"(from {tsv_path.name})")

            sources_meta.append({
                "module_type": module_type,
                "source_file": str(tsv_path),
                "records_in_file": len(records),
                "inserted": inserted,
                "updated": updated,
            })

        # Export TSV
        tsv_out = LIBRARY_DIR / "module_library.tsv"
        lib.to_tsv(tsv_out)
        log.info(f"Exported TSV: {tsv_out}")

        # Summary
        counts = lib.counts()
        log.info(f"Library counts: {counts}")

    write_build_manifest(
        LIBRARY_DIR,
        script_name="04_build_library.py",
        sources_used=sources_meta,
        records_added=total_inserted,
        notes=f"updated={total_updated}  rebuild={args.rebuild}",
    )

    log.info(f"Total: {total_inserted} inserted  {total_updated} updated")
    log.info("=== 04_build_library.py complete ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
