"""
literature/01_search.py
-----------------------
Query PubMed and Semantic Scholar for new papers on DBDs, EDs, and CRs.
Deduplicates against papers already in papers.yaml.
Saves raw candidates to literature/candidates/YYYYMMDD.json.

Run from project root:
    python literature/01_search.py [--config literature/search_queries.yaml]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import date, datetime, timezone
from pathlib import Path

import requests
import yaml

ROOT = Path(__file__).resolve().parents[1]
LIT_DIR = ROOT / "literature"
CANDIDATES_DIR = LIT_DIR / "candidates"


# ── HTTP helper ───────────────────────────────────────────────────────────────

def _get(url: str, params: dict | None = None, session: requests.Session | None = None,
         timeout: int = 20) -> dict | str | None:
    s = session or requests.Session()
    for attempt in range(3):
        try:
            r = s.get(url, params=params, timeout=timeout)
            if r.status_code == 200:
                ct = r.headers.get("Content-Type", "")
                return r.json() if "json" in ct else r.text
            if r.status_code == 429:
                wait = int(r.headers.get("Retry-After", 10))
                print(f"  [rate limit] waiting {wait}s...")
                time.sleep(wait)
        except requests.RequestException as e:
            print(f"  [warn] attempt {attempt+1} failed: {e}")
        time.sleep(2 ** attempt)
    return None


# ── PubMed ────────────────────────────────────────────────────────────────────

PUBMED_SEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_FETCH  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
PUBMED_SUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def search_pubmed(query: str, min_year: int, max_results: int) -> list[dict]:
    params = {
        "db": "pubmed", "term": query,
        "mindate": str(min_year), "maxdate": str(date.today().year),
        "datetype": "pdat", "retmax": max_results,
        "retmode": "json", "usehistory": "y",
    }
    data = _get(PUBMED_SEARCH, params=params)
    if not data or "esearchresult" not in data:
        return []

    ids = data["esearchresult"].get("idlist", [])
    if not ids:
        return []

    time.sleep(0.4)  # NCBI rate limit: max 3 req/s without API key

    summary = _get(PUBMED_SUMMARY, params={
        "db": "pubmed", "id": ",".join(ids), "retmode": "json"
    })
    if not summary:
        return []

    results = []
    uids = summary.get("result", {}).get("uids", [])
    for uid in uids:
        rec = summary["result"].get(uid, {})
        doi = next(
            (a["value"] for a in rec.get("articleids", []) if a["idtype"] == "doi"),
            None
        )
        pub_date = rec.get("pubdate", "")
        year = int(pub_date[:4]) if pub_date and pub_date[:4].isdigit() else None
        results.append({
            "source": "pubmed",
            "pmid": uid,
            "doi": doi.lower() if doi else None,
            "title": rec.get("title", "").rstrip("."),
            "journal": rec.get("fulljournalname", rec.get("source", "")),
            "year": year,
            "abstract": None,  # fetched separately only if needed
        })
    return results


# ── Semantic Scholar ──────────────────────────────────────────────────────────

SS_SEARCH = "https://api.semanticscholar.org/graph/v1/paper/search"


def search_semantic_scholar(query: str, min_year: int, max_results: int) -> list[dict]:
    params = {
        "query": query,
        "limit": min(max_results, 100),
        "fields": "title,year,journal,externalIds,abstract",
        "yearFilter": f"{min_year}-{date.today().year}",
    }
    data = _get(SS_SEARCH, params=params)
    if not data or "data" not in data:
        return []

    results = []
    for paper in data["data"]:
        ext = paper.get("externalIds") or {}
        doi = ext.get("DOI", "").lower() or None
        journal_info = paper.get("journal") or {}
        results.append({
            "source": "semantic_scholar",
            "pmid": ext.get("PubMed"),
            "doi": doi,
            "title": paper.get("title", ""),
            "journal": journal_info.get("name", ""),
            "year": paper.get("year"),
            "abstract": (paper.get("abstract") or "")[:500],
        })
    return results


# ── Deduplication ─────────────────────────────────────────────────────────────

def load_known_dois(papers_yaml: Path) -> set[str]:
    if not papers_yaml.exists():
        return set()
    with open(papers_yaml, encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    return {
        p["doi"].lower()
        for p in data.get("papers", [])
        if p.get("doi")
    }


def deduplicate(results: list[dict]) -> list[dict]:
    seen: set[str] = set()
    out = []
    for r in results:
        key = r.get("doi") or r.get("pmid") or r.get("title", "").lower()[:80]
        if key and key not in seen:
            seen.add(key)
            out.append(r)
    return out


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="literature/search_queries.yaml")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print query plan without making API calls")
    args = parser.parse_args()

    cfg_path = ROOT / args.config
    with open(cfg_path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    known_dois = load_known_dois(LIT_DIR / "papers.yaml")
    print(f"Known papers: {len(known_dois)}")

    min_year = cfg.get("min_year", 2018)
    max_results = cfg.get("max_results_per_query", 50)

    all_candidates: list[dict] = []

    for q in cfg["queries"]:
        qid = q["id"]
        category = q["category"]
        print(f"\n[{qid}] category={category}")

        if args.dry_run:
            print(f"  PubMed query: {q['pubmed'][:80]}...")
            print(f"  SS query:     {q['semantic_scholar']}")
            continue

        # PubMed
        pm_results = search_pubmed(q["pubmed"], min_year, max_results)
        print(f"  PubMed: {len(pm_results)} results")
        for r in pm_results:
            r["query_id"] = qid
            r["category"] = category
        all_candidates.extend(pm_results)
        time.sleep(0.4)

        # Semantic Scholar
        ss_results = search_semantic_scholar(q["semantic_scholar"], min_year, max_results)
        print(f"  Semantic Scholar: {len(ss_results)} results")
        for r in ss_results:
            r["query_id"] = qid
            r["category"] = category
        all_candidates.extend(ss_results)
        time.sleep(1.0)

    if args.dry_run:
        print("\n[dry-run] no API calls made")
        return

    # Deduplicate within this run, then filter out known papers
    all_candidates = deduplicate(all_candidates)
    new_candidates = [
        c for c in all_candidates
        if not (c.get("doi") and c["doi"] in known_dois)
    ]

    print(f"\nTotal fetched: {len(all_candidates)}")
    print(f"Already known: {len(all_candidates) - len(new_candidates)}")
    print(f"New candidates: {len(new_candidates)}")

    if not new_candidates:
        print("Nothing new — no candidate file written.")
        return

    CANDIDATES_DIR.mkdir(parents=True, exist_ok=True)
    out_path = CANDIDATES_DIR / f"{date.today().isoformat()}.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump({
            "run_date": datetime.now(timezone.utc).isoformat(),
            "total_fetched": len(all_candidates),
            "new_candidates": len(new_candidates),
            "candidates": new_candidates,
        }, f, indent=2, ensure_ascii=False)

    print(f"\nSaved to {out_path}")
    print("Next step: python literature/02_triage.py")


if __name__ == "__main__":
    main()
