#!/usr/bin/env python
"""Build the human review sheet for flagged / unreviewed DBD spans.

`06_extract_dbd.py` writes data/processed/dbd_review_queue.tsv -- 465 records,
worst kind first. That file says WHICH records to look at. It does not carry
the evidence needed to judge one, so reviewing from it means opening the
UniProt entry and the extraction report side by side for every row.

This writes docs/DBD_REVIEW_QUEUE.md instead: the same records, grouped, with
the evidence inline and a verdict column to fill in. It is read-only -- it
re-derives every span from the committed TSV plus the caches and touches
neither the database nor the library.

Two tiers are expanded in full, the two the external validation singled out:

  pfam    30 records. Arm 1 of docs/DBD_EXTERNAL_VALIDATION.md compares spans
          against Pfam, so for these the comparison is circular -- the span IS
          Pfam. Nothing has ever checked them. The question per record is
          whether Pfam named the right family, and the sheet shows what UniProt
          offered instead (never a DNA-binding feature, or the fallback would
          not have run).
  domain  36 flagged records whose span rests on a DOMAIN name alone, with no
          DNA_BIND or ZN_FING feature under it. Weakest kind on both arms: 93%
          recall >= 0.95 against 98-99% for the others, and 40 of 51 with full
          contact coverage against 91/97 for dna_bind. The question is whether
          the matched domain is really the DNA-binding module.

The rest of the queue is summarised by kind rather than expanded; those spans
rest on a UniProt DNA_BIND or ZN_FING line, which is the evidence arm 1 and
arm 2 already test.

Run:
    python scripts/review_sheet.py
    python scripts/review_sheet.py --tsv   # also write the flat sheet
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from pfam_dbd import dbd_domains  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TSV = ROOT / "library" / "module_library.tsv"
UNIPROT = ROOT / "data" / "raw" / "uniprot"
DOC_OUT = ROOT / "docs" / "DBD_REVIEW_QUEUE.md"
SHEET_OUT = ROOT / "data" / "processed" / "dbd_review_sheet.tsv"

# imported late so this module can be read without the extractor's constants
_ex = None


def _extractor():
    """Load 06_extract_dbd.py by path -- its name is not an identifier."""
    global _ex
    if _ex is None:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "_ex06", Path(__file__).parent / "06_extract_dbd.py")
        _ex = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = _ex
        spec.loader.exec_module(_ex)
    return _ex


def _feature_summary(feats: dict, m) -> dict:
    """What the UniProt table offered, split the way extract_one splits it."""
    kept_zn, other_zn = [], []
    for s, e, note in feats.get("ZN_FING", []):
        low = (note or "").lower()
        (kept_zn if any(t in low for t in m.ZNF_KEEP_TYPES)
         else other_zn).append(f"{note} {s}-{e}")
    keep_dom, excl_dom, other_dom = [], [], []
    for s, e, note in feats.get("DOMAIN", []):
        v = m.classify_domain(note or "")
        tag = f"{note} {s}-{e}"
        (keep_dom if v == "keep" else excl_dom if v == "exclude"
         else other_dom).append(tag)
    dna_bind = [f"{note or 'DNA_BIND'} {s}-{e}"
                for s, e, note in feats.get("DNA_BIND", [])]
    regions = []
    for s, e, note in feats.get("REGION", []):
        clean = " ".join((note or "").split()).lower().rstrip(".")
        if clean in m.KEEP_REGION_EXACT:
            regions.append(f"{note} {s}-{e}")
    return {"zn_kept": kept_zn, "zn_other": other_zn, "dom_keep": keep_dom,
            "dom_excluded": excl_dom, "dom_other": other_dom,
            "dna_bind": dna_bind, "region_keep": regions}


def collect() -> list[dict]:
    m = _extractor()
    df = pd.read_csv(TSV, sep="\t", dtype=str, low_memory=False)
    df = df[(df["type"] == "DBD") & df["uniprot_id"].notna()
            & (df["uniprot_id"].str.strip() != "") & df["sequence_aa"].notna()]
    rows = []
    for r in df.to_dict("records"):
        acc = str(r["uniprot_id"]).strip()
        p = UNIPROT / f"{acc}_feat.tsv"
        if not p.exists():
            continue
        seq = str(r["sequence_aa"])
        feats = m.parse_features(p.read_text(encoding="utf-8"))
        res = m.extract_one(seq, feats, acc)
        flagged = bool(res["dbd_review_flag"] or res["dbd_dropped_review_flag"])
        kind = res["dbd_extraction"]
        promoted = bool(res.get("_region_promoted"))
        if not (flagged or kind == "pfam" or promoted):
            continue
        rows.append({
            "gene_symbol": r["gene_symbol"], "uniprot_id": acc,
            "organism": r.get("organism", ""), "family": r.get("subtype", ""),
            "protein_aa": len(seq),
            "dbd_start": res["dbd_start"], "dbd_end": res["dbd_end"],
            "dbd_length_aa": res["dbd_length_aa"],
            "fraction": (round(res["dbd_length_aa"] / len(seq), 3)
                         if res["dbd_length_aa"] else None),
            "kind": kind,
            "review_flag": res["dbd_review_flag"],
            "dropped_review_flag": res["dbd_dropped_review_flag"],
            "flag_reasons": "; ".join(res.get("_inside", [])),
            "dropped_needs_review": "; ".join(res.get("_dropped_flagged", [])),
            "excluded_domains": res["dbd_excluded_domains"] or "",
            "pfam_families": "; ".join(res.get("_pfam_families", [])),
            "region_promoted": int(promoted),
            "zn_count_span": res["dbd_zinc_finger_count"],
            **{f"uni_{k}": "; ".join(v)
               for k, v in _feature_summary(feats, m).items()},
        })
    return rows


PRIORITY = {"pfam": 1, "domain": 2, "bipartite": 3, "znf_array": 4,
            "dna_bind": 5, "dispersed": 6, "none": 7, "no_features": 8}


def _fmt_reasons(r) -> str:
    """Group reasons that differ only by coordinates.

    The self-association detector emits one entry per run, so a record with
    several runs of the same module repeats the whole description: ARNT carried
    the same 4-part reason twice, 300 characters each, differing in two
    residues. Collapse to one description with its coordinate list.
    """
    raw = r["flag_reasons"] or r["dropped_needs_review"]
    if not raw:
        return "—"
    groups: dict = {}
    for part in raw.split("; "):
        head, sep, body = part.partition(":")
        if sep and head.replace("-", "").replace(" ", "").isdigit():
            groups.setdefault(body, []).append(head)
        else:
            groups.setdefault(part, [])
    out = []
    for body, coords in groups.items():
        out.append(f"{', '.join(coords)}: {body}" if coords else body)
    return "; ".join(out)


def _span(r) -> str:
    if r["dbd_start"] is None:
        return "— (no span)"
    return (f"`{r['dbd_start']}-{r['dbd_end']}` "
            f"({r['dbd_length_aa']} aa, {r['fraction']:.0%} of "
            f"{r['protein_aa']})")


def _pfam_block(rows, out):
    by_family = defaultdict(list)
    for r in rows:
        by_family[r["pfam_families"] or "(none)"].append(r)
    out.append("The fallback runs **only** when the UniProt table names no "
               "keep-feature, so the \"UniProt offered\" column never contains "
               "a DNA-binding feature. It is there to show what the span is "
               "standing in for.\n")
    for fam, recs in sorted(by_family.items(), key=lambda kv: -len(kv[1])):
        out.append(f"#### {fam} — {len(recs)} record(s)\n")
        out.append("| gene | span | UniProt offered | flagged for | verdict |")
        out.append("|---|---|---|---|---|")
        for r in sorted(recs, key=lambda x: x["gene_symbol"] or ""):
            offered = (r["uni_dom_other"] or r["uni_dom_excluded"]
                       or r["uni_zn_other"] or "nothing")
            out.append(f"| {r['gene_symbol']} ({r['uniprot_id']}) | "
                       f"{_span(r)} | {offered} | "
                       f"{_fmt_reasons(r)} | |")
        out.append("")


def _region_block(rows, out):
    out.append("A `REGION` whose note names it as the DNA-binding module, "
               "promoted to a keep-feature on abstention (`KEEP_REGION_EXACT`). "
               "Like the Pfam spans these have never been reviewed: the rule "
               "is newer than both validation arms. It has a lower bound "
               "(`MIN_PLAUSIBLE_SPAN_AA`) and **no upper bound**, so the "
               "question per record is whether UniProt's region is a domain or "
               "a tract — MLLT10's is 364 aa and clears every gate "
               "unflagged.\n")
    out.append("| gene | span | the REGION | also present | flagged for | verdict |")
    out.append("|---|---|---|---|---|---|")
    for r in sorted(rows, key=lambda x: -(x["dbd_length_aa"] or 0)):
        also = "; ".join(x for x in (r["uni_dom_other"], r["uni_dom_excluded"],
                                     r["uni_zn_other"]) if x) or "—"
        out.append(f"| {r['gene_symbol']} ({r['uniprot_id']}) | {_span(r)} | "
                   f"{r['uni_region_keep'] or '—'} | {also} | "
                   f"{_fmt_reasons(r)} | |")
    out.append("")


def _domain_block(rows, out):
    out.append("Each span rests on the DOMAIN name in \"matched on\", with no "
               "`DNA_BIND` or `ZN_FING` feature beneath it. \"also present\" "
               "lists what the extraction dropped, which is where an "
               "over-inclusive span usually shows up.\n")
    out.append("| gene | span | matched on | also present | flagged for | verdict |")
    out.append("|---|---|---|---|---|---|")
    for r in sorted(rows, key=lambda x: x["gene_symbol"] or ""):
        also = "; ".join(x for x in (r["uni_dom_excluded"],
                                     r["uni_dom_other"]) if x) or "—"
        out.append(f"| {r['gene_symbol']} ({r['uniprot_id']}) | {_span(r)} | "
                   f"{r['uni_dom_keep'] or '—'} | {also} | "
                   f"{_fmt_reasons(r)} | |")
    out.append("")


def render(rows) -> str:
    rows.sort(key=lambda r: (PRIORITY.get(r["kind"], 9), r["gene_symbol"] or ""))
    counts = Counter(r["kind"] for r in rows)
    pfam = [r for r in rows if r["kind"] == "pfam"]
    domain = [r for r in rows if r["kind"] == "domain"]
    region = [r for r in rows if r["region_promoted"]]

    out = ["# DBD review queue",
           "",
           "Generated by `scripts/review_sheet.py` from the committed "
           "`library/module_library.tsv` and the UniProt / InterPro caches. "
           "Read-only: it re-derives each span rather than reading the "
           "database.",
           "",
           f"**{len(rows)} records** need a look — every record carrying "
           "`dbd_review_flag` or `dbd_dropped_review_flag`, plus every span "
           "whose provenance is a source the extractor did not use before "
           "(Pfam, or a promoted REGION), flagged or not.",
           "",
           "| kind | records | why it is in this queue |",
           "|---|---|---|"]
    why = {
        "pfam": "arm 1 cannot check these — the span **is** Pfam",
        "domain": "weakest kind on both arms; rests on a DOMAIN name alone",
        "bipartite": "two modules bracketed; check the gap between them",
        "znf_array": "finger whitelist and array boundaries",
        "dna_bind": "a UniProt DNA_BIND line, plus a flag reason",
        "dispersed": "no span written; the extraction refused",
        "none": "abstained",
        "no_features": "fetch returned nothing",
    }
    for kind, n in sorted(counts.items(), key=lambda kv: PRIORITY.get(kv[0], 9)):
        out.append(f"| `{kind}` | {n} | {why.get(kind, '')} |")
    out += ["",
            "The two tiers below are expanded in full. Every other record "
            "stands on a UniProt `DNA_BIND` or `ZN_FING` line — the evidence "
            "arm 1 and arm 2 already test — so those are left to the flat "
            "sheet (`--tsv`).",
            "",
            "Fill in **verdict** per row: `ok`, `wrong`, or a note. A `wrong` "
            "verdict is a keep-list or threshold change, not a hand-edit of "
            "the library — the extraction is reproducible and a manual span "
            "would be overwritten by the next rebuild.",
            ""]

    out += ["## Tier 1 — `pfam` "
            f"({len(pfam)} records, none ever reviewed)", ""]
    _pfam_block(pfam, out)
    out += [f"## Tier 1b — promoted `REGION` ({len(region)} records, "
            "none ever reviewed)", ""]
    _region_block(region, out)
    out += [f"## Tier 2 — `domain` ({len(domain)} flagged records)", ""]
    _domain_block(domain, out)
    return "\n".join(out) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tsv", action="store_true",
                    help=f"also write the flat sheet to {SHEET_OUT}")
    args = ap.parse_args()

    rows = collect()
    DOC_OUT.write_text(render(rows), encoding="utf-8", newline="\n")
    print(f"[review] {len(rows):,} records -> {DOC_OUT}")
    if args.tsv:
        SHEET_OUT.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(SHEET_OUT, sep="\t", index=False)
        print(f"[review] flat sheet -> {SHEET_OUT}")


if __name__ == "__main__":
    main()
