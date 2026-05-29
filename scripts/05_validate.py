"""
scripts/05_validate.py
----------------------
Quality-control checks on the assembled module library.

Checks
------
  1. Required fields present and non-null for every record
  2. validation_level values match the allowed enum
  3. sequence_aa length matches length_aa (where both are set)
  4. No duplicate module_ids
  5. type values match the allowed set
  6. ED screen-validated entries have a quantitative_metric
  7. DBD entries have a jaspar_id or uniprot_id (at least one)

Outputs
-------
  library/qc_report.tsv   one row per failed check (empty = all passed)
  Prints a summary table to stdout

Run
---
  python scripts/05_validate.py [--fix]

  --fix : attempt to auto-correct minor issues (e.g. recalculate length_aa)
"""

import argparse
import sys
import io
from pathlib import Path

# Ensure Unicode output works on Windows GBK terminals
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from schema import ModuleLibrary
from utils import get_logger

ROOT = Path(__file__).resolve().parents[1]
LOGS_DIR = ROOT / "logs"
LIBRARY_DIR = ROOT / "library"

VALID_TYPES = {"DBD", "ED", "CR"}
VALID_LEVELS = {"predicted", "motif-only", "ChIP-validated",
                "screen-validated", "structurally-resolved"}
REQUIRED_COLS = ["module_id", "type", "name", "validation_level",
                 "source", "date_added", "date_modified"]


def run_checks(df: pd.DataFrame, fix: bool, log) -> pd.DataFrame:
    """Run all checks. Returns a DataFrame of failures."""
    failures = []

    def fail(module_id, check, detail):
        failures.append({"module_id": module_id, "check": check, "detail": detail})

    for _, row in df.iterrows():
        mid = row.get("module_id", "?")

        # 1. Required fields
        for col in REQUIRED_COLS:
            val = row.get(col)
            if val is None or (isinstance(val, float) and pd.isna(val)) or str(val).strip() == "":
                fail(mid, "required_field_missing", f"column={col}")

        # 2. Valid type
        t = str(row.get("type", "")).strip()
        if t not in VALID_TYPES:
            fail(mid, "invalid_type", f"type='{t}'")

        # 3. Valid validation_level
        vl = str(row.get("validation_level", "")).strip()
        if vl not in VALID_LEVELS:
            fail(mid, "invalid_validation_level", f"validation_level='{vl}'")

        # 4. sequence_aa / length_aa consistency
        seq = row.get("sequence_aa")
        length = row.get("length_aa")
        if seq and not (isinstance(seq, float) and pd.isna(seq)):
            actual_len = len(str(seq))
            if length and not (isinstance(length, float) and pd.isna(length)):
                stored_len = int(float(length))
                if abs(actual_len - stored_len) > 2:   # ±2 tolerance for stripped chars
                    if fix:
                        df.loc[df["module_id"] == mid, "length_aa"] = actual_len
                        log.info(f"[fix] {mid}: corrected length_aa {stored_len}→{actual_len}")
                    else:
                        fail(mid, "length_mismatch",
                             f"stored={stored_len}  actual={actual_len}")

        # 5. ED from a screen dataset should have a quantitative_metric.
        # Curated canonical entries (source starts with 'Curated:') are exempt —
        # they are engineering-validated but lack a single numeric score.
        if t == "ED" and vl == "screen-validated":
            source = str(row.get("source", ""))
            if not source.startswith("Curated:"):
                qm = row.get("quantitative_metric")
                if qm is None or (isinstance(qm, float) and pd.isna(qm)):
                    fail(mid, "missing_metric",
                         "screen-sourced ED has no quantitative_metric")

        # 6. DBD should have jaspar_id or uniprot_id
        if t == "DBD":
            ji = row.get("jaspar_id")
            ui = row.get("uniprot_id")
            def _empty(x):
                return x is None or (isinstance(x, float) and pd.isna(x)) or str(x).strip() in ("", "nan")
            if _empty(ji) and _empty(ui):
                fail(mid, "dbd_no_identifier",
                     "DBD entry has neither jaspar_id nor uniprot_id")

    # 7. Duplicate module_ids
    dupes = df[df.duplicated("module_id", keep=False)]["module_id"].unique()
    for mid in dupes:
        fail(mid, "duplicate_module_id", "module_id appears more than once")

    return pd.DataFrame(failures)


def print_summary(df: pd.DataFrame, failures: pd.DataFrame, log) -> None:
    print("\n" + "═" * 60)
    print("  Module Library QC Report")
    print("═" * 60)

    by_type = df.groupby("type").size().to_dict()
    by_level = df.groupby("validation_level").size().sort_values(ascending=False).to_dict()

    print(f"\n  Total entries : {len(df):,}")
    print(f"\n  By type:")
    for k, v in sorted(by_type.items()):
        print(f"    {k:<10} {v:>6,}")
    print(f"\n  By validation_level:")
    for k, v in by_level.items():
        print(f"    {k:<30} {v:>6,}")

    # Sequence coverage
    has_seq = df["sequence_aa"].notna().sum()
    print(f"\n  Sequence coverage : {has_seq:,} / {len(df):,} ({100*has_seq/len(df):.1f}%)")

    if len(failures) == 0:
        print("\n  ✓  All checks passed.\n")
    else:
        print(f"\n  ✗  {len(failures)} issue(s) found:\n")
        by_check = failures.groupby("check").size().sort_values(ascending=False)
        for check, count in by_check.items():
            print(f"    {check:<35} {count:>4} record(s)")
        print()

    print("═" * 60 + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fix", action="store_true",
                        help="Auto-correct minor issues (recalculate length_aa)")
    args = parser.parse_args()

    log = get_logger("05_validate", LOGS_DIR)
    log.info("=== 05_validate.py started ===")

    db_path = LIBRARY_DIR / "module_library.db"
    if not db_path.exists():
        log.error(f"Library not found at {db_path}. Run 04_build_library.py first.")
        sys.exit(1)

    with ModuleLibrary(db_path) as lib:
        df = lib.to_dataframe()
        log.info(f"Loaded {len(df):,} records from library")

        failures = run_checks(df, args.fix, log)

        if args.fix and len(failures) < len(run_checks(df, False, log)):
            # Re-insert fixed records
            lib.insert_or_replace(df.to_dict(orient="records"))
            lib.to_tsv(LIBRARY_DIR / "module_library.tsv")
            log.info("[fix] re-exported TSV after corrections")

    print_summary(df, failures, log)

    qc_out = LIBRARY_DIR / "qc_report.tsv"
    failures.to_csv(qc_out, sep="\t", index=False)
    log.info(f"QC report written to {qc_out}")

    if len(failures) > 0:
        log.warning(f"{len(failures)} QC issue(s) found — review {qc_out}")
        return 1

    log.info("=== 05_validate.py complete — all checks passed ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
