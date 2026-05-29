"""
literature/02_triage.py
-----------------------
Score and rank candidate papers from 01_search.py.
Generates a Markdown review report in literature/reviews/.

Run from project root:
    python literature/02_triage.py [--candidates literature/candidates/YYYY-MM-DD.json]

If --candidates is omitted, uses the most recent candidates file.
"""

from __future__ import annotations

import argparse
import json
from datetime import date
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
LIT_DIR = ROOT / "literature"
CANDIDATES_DIR = LIT_DIR / "candidates"
REVIEWS_DIR = LIT_DIR / "reviews"


# ── Relevance scoring ─────────────────────────────────────────────────────────
# Keywords → score weight. Matched against title + abstract (case-insensitive).

SCORE_WEIGHTS: dict[str, float] = {
    # Strong indicators
    "high-throughput": 2.0,
    "screen": 1.5,
    "tiling": 2.0,
    "massively parallel": 2.0,
    "deep mutational": 2.0,
    "compendium": 1.5,
    # DBD
    "dna-binding domain": 2.0,
    "transcription factor": 1.0,
    "zinc finger": 1.5,
    "homeodomain": 1.5,
    "bzip": 1.5,
    "helix-turn-helix": 1.5,
    # ED
    "activation domain": 2.0,
    "effector domain": 2.0,
    "repression domain": 2.0,
    "transcriptional activator": 1.5,
    "transcriptional repressor": 1.5,
    "krab": 2.0,
    "vp16": 1.5,
    "p65": 1.0,
    # CR
    "chromatin remodel": 2.0,
    "histone acetyltransferase": 2.0,
    "histone deacetylase": 2.0,
    "histone methyltransferase": 2.0,
    "swi/snf": 1.5,
    "polycomb": 1.5,
    # Synthetic biology context
    "synthetic transcription factor": 2.0,
    "crispra": 1.5,
    "crispri": 1.5,
    "gene regulation": 1.0,
    "engineered": 1.0,
    # Negative indicators (reduce score)
    "review": -1.0,
    "commentary": -1.5,
    "erratum": -3.0,
    "correction": -3.0,
    # Loss-of-function / disease / inhibitor studies — not domain characterisation
    "complex assembly": -2.0,
    "loss of function": -1.5,
    "inhibitor": -1.5,
    "drug": -1.5,
    "therapeutic": -1.5,
    "cancer": -1.0,
    "tumor": -1.0,
    "oncogene": -1.0,
    # Non-human/mouse organisms
    "arabidopsis": -3.0,
    "yeast": -2.0,
    "drosophila": -2.0,
    "zebrafish": -2.0,
    "plant": -2.0,
    "fungal": -2.0,
    "bacterial": -2.0,
    "caenorhabditis": -2.0,
    "xenopus": -2.0,
}

CATEGORY_BONUS: dict[str, dict[str, float]] = {
    "DBD": {"dna-binding domain": 1.0, "transcription factor": 0.5, "zinc finger": 1.0},
    "ED":  {"activation domain": 1.0, "repression domain": 1.0, "effector": 0.5},
    "CR":  {"chromatin remodel": 1.0, "histone": 0.5, "epigenetic": 0.5},
}


def score_paper(paper: dict) -> float:
    text = f"{paper.get('title', '')} {paper.get('abstract', '')}".lower()
    score = 0.0
    for kw, weight in SCORE_WEIGHTS.items():
        if kw in text:
            score += weight
    # Category-specific bonus
    cat = paper.get("category", "")
    for kw, bonus in CATEGORY_BONUS.get(cat, {}).items():
        if kw in text:
            score += bonus
    return round(score, 2)


# ── Report generation ─────────────────────────────────────────────────────────

def render_report(candidates: list[dict], run_date: str, source_file: Path) -> str:
    lines = [
        f"# Literature Review — {date.today().isoformat()}",
        f"",
        f"**Search date:** {run_date[:10]}  ",
        f"**Candidates file:** `{source_file.name}`  ",
        f"**Papers to review:** {len(candidates)}",
        f"",
        "## Instructions",
        "",
        "For each paper below, add your decision to `literature/papers.yaml`:",
        "",
        "```yaml",
        "- doi: \"<DOI>\"",
        "  title: \"<title>\"",
        "  year: <year>",
        "  journal: \"<journal>\"",
        "  categories: [DBD|ED|CR]",
        "  status: accepted   # or: rejected",
        "  notes: \"supplementary file: mmc2.xlsx; relevant sheet: Table S3\"",
        "```",
        "",
        "Then run `python literature/03_record.py` to validate.",
        "",
        "---",
        "",
    ]

    # Group by category
    by_cat: dict[str, list[dict]] = {}
    for p in candidates:
        cat = p.get("category", "OTHER")
        by_cat.setdefault(cat, []).append(p)

    for cat in ["DBD", "ED", "CR", "OTHER"]:
        papers = by_cat.get(cat, [])
        if not papers:
            continue
        lines.append(f"## {cat} candidates ({len(papers)})")
        lines.append("")
        for p in papers:
            title = p.get("title") or "_(no title)_"
            doi = p.get("doi") or "_(no DOI)_"
            year = p.get("year") or "?"
            journal = p.get("journal") or "?"
            score = p.get("_score", 0)
            abstract = (p.get("abstract") or "").strip()
            abstract_preview = (abstract[:300] + "…") if len(abstract) > 300 else abstract

            lines += [
                f"### {title}",
                f"",
                f"- **DOI:** [{doi}](https://doi.org/{doi})",
                f"- **Year:** {year}  **Journal:** {journal}",
                f"- **Relevance score:** {score}  **Query:** `{p.get('query_id', '?')}`",
            ]
            if abstract_preview:
                lines += ["", f"> {abstract_preview}", ""]
            lines += [
                "**Decision:** `[ ] accepted` / `[ ] rejected`  ",
                "**Notes:** _(supplementary file, relevant sheets, why rejected, etc.)_",
                "",
                "---",
                "",
            ]

    return "\n".join(lines)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--candidates", default=None,
                        help="Path to candidates JSON. Defaults to most recent.")
    parser.add_argument("--min-score", type=float, default=2.0,
                        help="Minimum relevance score to include in report (default: 2.0)")
    args = parser.parse_args()

    # Find candidates file
    if args.candidates:
        cand_path = Path(args.candidates)
    else:
        files = sorted(CANDIDATES_DIR.glob("*.json"))
        if not files:
            print("No candidate files found. Run 01_search.py first.")
            return
        cand_path = files[-1]

    print(f"Loading candidates from {cand_path}")
    with open(cand_path, encoding="utf-8") as f:
        data = json.load(f)

    candidates = data.get("candidates", [])
    run_date = data.get("run_date", "unknown")

    # Score and filter
    for p in candidates:
        p["_score"] = score_paper(p)

    candidates = [p for p in candidates if p["_score"] >= args.min_score]
    candidates.sort(key=lambda p: (-p["_score"], p.get("year") or 0))

    print(f"After scoring (min={args.min_score}): {len(candidates)} candidates")

    if not candidates:
        print("No candidates above score threshold.")
        return

    # Generate report
    REVIEWS_DIR.mkdir(parents=True, exist_ok=True)
    report_path = REVIEWS_DIR / f"review_{date.today().isoformat()}.md"
    report = render_report(candidates, run_date, cand_path)
    report_path.write_text(report, encoding="utf-8")

    print(f"Review report written to {report_path}")
    print(f"\nTop 5 candidates:")
    for p in candidates[:5]:
        print(f"  [{p['_score']:4.1f}] {p.get('title','')[:70]}")
    print(f"\nNext step: review {report_path.name}, then update literature/papers.yaml")
    print("           and run: python literature/03_record.py")


if __name__ == "__main__":
    main()
