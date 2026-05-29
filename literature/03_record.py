"""
literature/03_record.py
-----------------------
Validate and finalise review decisions in papers.yaml.

After updating papers.yaml with accepted/rejected decisions, run this script
to verify the file is well-formed and print a summary.

Run from project root:
    python literature/03_record.py
"""

from __future__ import annotations

import io
import sys
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

import yaml

ROOT = Path(__file__).resolve().parents[1]
PAPERS_YAML = ROOT / "literature" / "papers.yaml"

VALID_STATUSES = {"integrated", "accepted", "rejected", "pending"}
VALID_CATEGORIES = {"DBD", "ED", "CR"}


def main():
    if not PAPERS_YAML.exists():
        print(f"papers.yaml not found at {PAPERS_YAML}")
        sys.exit(1)

    with open(PAPERS_YAML, encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}

    papers = data.get("papers", [])
    errors = []
    counts: dict[str, int] = {}

    for i, p in enumerate(papers):
        ref = p.get("doi") or p.get("title", f"entry #{i}")

        if not p.get("doi"):
            errors.append(f"{ref}: missing doi")
        if not p.get("title"):
            errors.append(f"{ref}: missing title")
        if not p.get("year"):
            errors.append(f"{ref}: missing year")

        status = p.get("status", "")
        if status not in VALID_STATUSES:
            errors.append(f"{ref}: invalid status '{status}' (must be one of {VALID_STATUSES})")

        cats = p.get("categories", [])
        for c in cats:
            if c not in VALID_CATEGORIES:
                errors.append(f"{ref}: unknown category '{c}'")

        counts[status] = counts.get(status, 0) + 1

    print("=" * 60)
    print("  papers.yaml validation")
    print("=" * 60)
    print(f"  Total entries : {len(papers)}")
    for status in ["integrated", "accepted", "pending", "rejected"]:
        n = counts.get(status, 0)
        if n:
            print(f"  {status:<12} : {n}")

    if errors:
        print(f"\n  {len(errors)} error(s) found:")
        for e in errors:
            print(f"    ✗  {e}")
        sys.exit(1)
    else:
        print("\n  ✓  All entries valid")

    # Remind about accepted papers not yet integrated
    accepted = [p for p in papers if p.get("status") == "accepted"]
    if accepted:
        print(f"\n  {len(accepted)} accepted paper(s) ready to integrate:")
        for p in accepted:
            print(f"    - [{p.get('year','')}] {p.get('title','')[:60]}")
            if p.get("notes"):
                print(f"        notes: {p['notes']}")
        print("\n  Add a new entry to SCREEN_FILES in scripts/02_seed_ed.py (or")
        print("  03_fetch_cr.py / 01_fetch_dbd.py) and re-run the pipeline.")
    print()


if __name__ == "__main__":
    main()
