"""
scripts/06_extract_dbd.py
-------------------------
Extract the DNA-binding span from the full-length sequences stored on
`type='DBD'` records.

Why
---
`01_fetch_dbd.py` writes the ENTIRE UniProt sequence into `sequence_aa`, so the
DBD type holds whole transcription factors: `DBD_HO_HNF1A` is all 631 aa of
HNF1A, whose homeodomain is 81. For motif scanning that is harmless (the PWM is
used, not the protein), but `fuse.py` concatenates `sequence_aa` to an effector
domain -- which builds a full-length-TF fusion carrying every native domain the
TF has.

That is not cosmetic. The fib->beta run (2026-09-03) selected ZNF793 and ZNF850
as OPENERS; both carry native KRAB domains, which recruit KAP1/TRIM28 -> SETDB1
-> H3K9me3. Fusing a BAF recruiter to a KRAB-bearing protein builds a construct
that opens and silences the same locus at once -- self-antagonism inside a single
chain, which no design-level rule (such as opener/closer DBD disjointness) can
see. zf-C2H2 is ~39% of the library, so this is systematic.

Method
------
Source of truth is the UniProt feature table (ft_zn_fing, ft_dna_bind, ft_domain,
ft_region), fetched in batches and cached per accession.

  KEEP  the DNA-contacting apparatus:
        - Zinc finger spans, taken as ONE CONTIGUOUS BLOCK from first start to
          last end. The inter-finger TGEKP linkers are part of the array;
          splicing fingers together individually destroys the register.
        - DNA binding spans (homeobox, HMG box, fork head, T-box, ...).
        - Partner domains that form part of the interface (KEEP_DOMAIN_HINTS).
          HNF1A binds as a bipartite POU-specific + homeodomain unit, so taking
          the homeobox alone would remove half the interface.
  DROP  everything else, explicitly including KRAB, SCAN, BTB/POZ and
        interaction/activation regions.

Then span = min(start) - FLANK .. max(end) + FLANK, clamped to the sequence.

Rules that matter
-----------------
1. Contiguous span, never concatenated fragments. Splicing across an excluded
   internal domain creates a junction that exists in no real protein -- a novel
   neoepitope surface and an unpredictable fold.
2. A dropped domain lying INSIDE the kept span is a conflict, not a silent pass:
   `dbd_review_flag` is set and the offender named in `dbd_excluded_domains`.
3. Never fabricate boundaries. No usable features -> `dbd_extraction='none'` and
   the DBD columns stay NULL. A guessed domain is worse than an honest
   full-length record.
4. Additive. `sequence_aa` is not touched; the original is also copied to
   `sequence_full_aa` so consumers can migrate deliberately.

Run
---
  python scripts/06_extract_dbd.py --survey     # domain vocabulary, no writes
  python scripts/06_extract_dbd.py --dry-run    # extract + report, no writes
  python scripts/06_extract_dbd.py              # extract and update the library
  python scripts/06_extract_dbd.py --validate   # the three acceptance gates

Cost: offline apart from the feature fetch (~30 batched requests); minutes, no GPU.
"""

from __future__ import annotations

import argparse
import re
import sys
import uuid
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import requests
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from schema import ModuleLibrary
from utils import get_logger

ROOT = Path(__file__).resolve().parents[1]
LOGS_DIR = ROOT / "logs"
RAW_DIR = ROOT / "data" / "raw"
PROCESSED_DIR = ROOT / "data" / "processed"
DB_PATH = ROOT / "library" / "module_library.db"

FLANK = 6
BATCH = 100          # accessions per UniProt search request

# Feature-table keys -> the TSV columns they arrive in.
FEATURE_COLUMNS = ["Zinc finger", "DNA binding", "Domain [FT]", "Region"]
FEATURE_FIELDS = "accession,ft_zn_fing,ft_dna_bind,ft_domain,ft_region"

# Domains that are part of the DNA-binding apparatus but are annotated as
# DOMAIN rather than DNA_BIND. Matched case-insensitively as substrings of the
# UniProt /note. Curated from the domain-description survey (--survey); a family
# whose binding domain is not on this list gets dbd_extraction='none' rather
# than a guessed span.
KEEP_DOMAIN_HINTS = (
    "pou-specific", "pou domain", "homeobox", "homeo",
    "bhlh", "bzip", "leucine zipper", "basic motif",
    "hmg box", "fork-head", "forkhead", "t-box",
    "myb-like", "myb-type", "sant", "mads-box", "ipt/tig",
    "grh/cp2", "cold-shock", "at-hook", "ap2/erf",
    "nuclear receptor", "nr c4-type", "cenp-b",
)

# Short domain names that must match the WHOLE description, not a substring of
# it: 'tea' is inside 'Protease' and 'rel' inside 'Related', so substring
# matching here would keep domains that have nothing to do with DNA. The note is
# compared with its trailing repeat index stripped ('Myb 2' -> 'myb').
KEEP_DOMAIN_EXACT = frozenset({
    "tea", "ets", "rhd", "mh1", "mbd", "cut", "pax", "sand", "csd",
    "cbf", "runt", "stat", "tbp", "gcm", "arid", "psq", "hth psq-type",
})

# Domains that must never ride along on a fusion. KRAB/SCAN/BTB recruit their own
# corepressors; carrying one into an opener is the failure this script exists for.
EXCLUDE_DOMAIN_HINTS = ("krab", "scan", "btb", "poz", "bromo", "chromo", "phd",
                        "set ", "jmjc", "sam ", "wd 1", "ank ", "ring")

# Regions with an effector function of their own. Excluded always (they are
# Regions, never kept), but flagged for review when they fall inside the kept
# span. "Interaction with DNA" is deliberately NOT here -- those ARE the contacts.
HAZARD_REGION_HINTS = ("activation", "transactivation", "repress", "silenc",
                       "corepressor", "coactivator")

FEATURE_RE = re.compile(
    r"(ZN_FING|DNA_BIND|DOMAIN|REGION)\s+[<>?]?(\d+)\.\.[<>?]?(\d+)"
    r'(?:;\s*/note="([^"]*)")?'
)


# ── UniProt feature fetch ────────────────────────────────────────────────────

def _cache_path(acc: str) -> Path:
    return RAW_DIR / "uniprot" / f"{acc}_feat.tsv"


def fetch_features(accessions: list[str], session: requests.Session,
                   api_base: str, log) -> dict[str, str]:
    """
    Return {accession: cached TSV line}. Accessions already on disk are not
    re-fetched; the rest go out in batches of BATCH via the search endpoint.
    A batch is ~1 request per 100 entries, so the whole library is ~30 requests
    -- the network cost that deferred this work does not really exist.
    """
    out: dict[str, str] = {}
    todo = []
    for acc in accessions:
        p = _cache_path(acc)
        if p.exists():
            out[acc] = p.read_text(encoding="utf-8")
        else:
            todo.append(acc)

    log.info(f"[features] {len(out):,} cached, {len(todo):,} to fetch "
             f"({(len(todo) + BATCH - 1) // BATCH} requests)")

    for i in range(0, len(todo), BATCH):
        chunk = todo[i:i + BATCH]
        query = "+OR+".join(f"accession:{a}" for a in chunk)
        url = (f"{api_base}/search?query={query}&format=tsv"
               f"&fields={FEATURE_FIELDS}&size={BATCH}")
        try:
            resp = session.get(url, timeout=60)
            resp.raise_for_status()
        except requests.RequestException as exc:
            log.warning(f"[features] batch {i // BATCH}: {exc}")
            continue

        lines = resp.text.splitlines()
        got = set()
        for line in lines[1:]:
            if not line.strip():
                continue
            acc = line.split("\t", 1)[0].strip()
            # An accession may be secondary/merged: UniProt answers under the
            # primary accession, so only trust rows we actually asked for.
            if acc not in chunk:
                continue
            _cache_path(acc).write_text(line, encoding="utf-8")
            out[acc] = line
            got.add(acc)
        missing = set(chunk) - got
        if missing:
            log.debug(f"[features] no entry returned for {len(missing)} accessions")
        log.info(f"[features] batch {i // BATCH + 1}: {len(got)}/{len(chunk)}")

    return out


def parse_features(line: str) -> dict[str, list[tuple[int, int, str]]]:
    """Parse one feature TSV row into {ZN_FING|DNA_BIND|DOMAIN|REGION: [(s,e,note)]}."""
    feats: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for kind, s, e, note in FEATURE_RE.findall(line):
        feats[kind].append((int(s), int(e), note or ""))
    return feats


# ── Extraction ───────────────────────────────────────────────────────────────

_TRAILING_INDEX = re.compile(r"\s+\d+$")


def classify_domain(note: str) -> str:
    """'keep' | 'exclude' | 'other' for a DOMAIN description.

    Exclusion wins over keeping: a description that reads as both should not be
    carried into a fusion on the strength of the weaker signal."""
    low = note.lower().strip()
    if any(h in low for h in EXCLUDE_DOMAIN_HINTS):
        return "exclude"
    if any(h in low for h in KEEP_DOMAIN_HINTS):
        return "keep"
    if _TRAILING_INDEX.sub("", low) in KEEP_DOMAIN_EXACT:
        return "keep"
    return "other"


def extract_one(seq: str, feats: dict) -> dict:
    """
    Compute the DNA-binding span for one record. Returns a dict of the DBD
    columns; `dbd_extraction='none'` when no keep-feature was found.
    """
    zn = [(s, e) for s, e, _ in feats.get("ZN_FING", [])]
    dna = [(s, e) for s, e, _ in feats.get("DNA_BIND", [])]
    keep_dom, dropped = [], []
    for s, e, note in feats.get("DOMAIN", []):
        verdict = classify_domain(note)
        if verdict == "keep":
            keep_dom.append((s, e))
        else:
            dropped.append((note or "?", s, e))

    spans = zn + dna + keep_dom
    if not spans:
        return {"dbd_extraction": "none", "dbd_start": None, "dbd_end": None,
                "dbd_sequence_aa": None, "dbd_length_aa": None,
                "dbd_excluded_domains": None, "dbd_review_flag": 0,
                "n_znf": len(zn)}

    lo = max(1, min(s for s, _ in spans) - FLANK)
    hi = min(len(seq), max(e for _, e in spans) + FLANK)
    sub = seq[lo - 1:hi]

    # Rule 2: a dropped domain inside the kept span is carried along anyway.
    inside = [n for n, s, e in dropped if s >= lo and e <= hi]
    for s, e, note in feats.get("REGION", []):
        low = (note or "").lower()
        if any(h in low for h in HAZARD_REGION_HINTS) and s >= lo and e <= hi:
            inside.append(note)

    if zn and (dna or keep_dom):
        kind = "bipartite"
    elif zn:
        kind = "znf_array"
    elif dna and keep_dom:
        kind = "bipartite"
    elif dna:
        kind = "dna_bind"
    else:
        kind = "domain"

    return {
        "dbd_extraction": kind,
        "dbd_start": lo,
        "dbd_end": hi,
        "dbd_sequence_aa": sub,
        "dbd_length_aa": len(sub),
        "dbd_excluded_domains": ";".join(n for n, _, _ in dropped) or None,
        "dbd_review_flag": 1 if inside else 0,
        "n_znf": len(zn),
        "_inside": inside,
    }


# ── Survey (no writes) ───────────────────────────────────────────────────────

def survey(rows, features, log) -> None:
    """Print the DOMAIN vocabulary, split by verdict and by family.

    The keep list above is only trustworthy if the domains it does NOT name are
    genuinely not DNA-binding, so this prints what would be dropped and what
    families would end up with no span at all."""
    by_verdict = Counter()
    unknown = Counter()
    families_none = Counter()
    families_total = Counter()

    for r in rows:
        acc = r["uniprot_id"]
        fam = r["subtype"] or "?"
        families_total[fam] += 1
        feats = parse_features(features.get(acc, ""))
        for _, _, note in feats.get("DOMAIN", []):
            v = classify_domain(note)
            by_verdict[v] += 1
            if v == "other":
                unknown[note] += 1
        if not (feats.get("ZN_FING") or feats.get("DNA_BIND")
                or any(classify_domain(n) == "keep" for _, _, n in feats.get("DOMAIN", []))):
            families_none[fam] += 1

    log.info(f"[survey] DOMAIN verdicts: {dict(by_verdict)}")
    log.info("[survey] top unclassified DOMAIN descriptions:")
    for note, n in unknown.most_common(40):
        log.info(f"    {n:5d}  {note}")
    log.info("[survey] records with NO keep-feature, by family "
             "(these would get dbd_extraction='none'):")
    for fam, n in families_none.most_common(30):
        log.info(f"    {n:5d} / {families_total[fam]:5d}  {fam}")


# ── Validation gates ─────────────────────────────────────────────────────────

# Hand-checked spans from the fib->beta run (2026-09-03), used as a regression
# fixture. If a keep-list edit moves one of these, that is a deliberate decision
# to be re-justified, not a silent drift.
ANCHORS = {
    "HNF1A":  (81, 285),
    "EGR3":   (269, 361),
    "ZNF793": (221, 395),
    "ZNF850": (191, 1090),
    "ZNF540": (181, 660),
    "ZNF444": (173, 306),
    "HINFP":  (9, 374),
}

KRAB_HINT = "krab"


def validate(lib, log) -> int:
    """Three acceptance gates. Returns the number of failures."""
    df = lib.to_dataframe("DBD")
    done = df[df["dbd_extraction"].notna() & (df["dbd_extraction"] != "none")]
    fails = 0

    # Gate 1 -- length sanity against the independently stored finger count.
    # A C2H2 array runs about 23 aa per finger plus ~5 aa TGEKP linkers, so
    # 28 aa/finger is the expected pitch; allow generous slack for degenerate
    # fingers and split arrays, and only flag spans that are too SHORT to hold
    # the fingers claimed (too long is legitimate -- split arrays span a gap).
    znf = done[done["zinc_finger_count"].notna() & (done["zinc_finger_count"] > 0)
               & (done["dbd_extraction"].isin(["znf_array", "bipartite"]))]
    short = znf[znf["dbd_length_aa"] < 20 * znf["zinc_finger_count"]]
    log.info(f"[gate 1/length] {len(znf):,} ZnF records checked, "
             f"{len(short):,} shorter than 20 aa/finger")
    if len(short):
        fails += 1
        for _, r in short.head(10).iterrows():
            log.warning(f"    {r['gene_symbol']}: {r['dbd_length_aa']} aa for "
                        f"{int(r['zinc_finger_count'])} fingers")

    # Gate 2 -- known anchors.
    bad = []
    for gene, (s, e) in ANCHORS.items():
        hit = df[(df["gene_symbol"] == gene) & (df["organism"].astype(str).str.contains("sapiens", case=False, na=False))]
        if hit.empty:
            hit = df[df["gene_symbol"] == gene]
        if hit.empty:
            bad.append(f"{gene}: not in library"); continue
        r = hit.iloc[0]
        if (r["dbd_start"], r["dbd_end"]) != (s, e):
            bad.append(f"{gene}: got {r['dbd_start']}-{r['dbd_end']}, expected {s}-{e}")
    log.info(f"[gate 2/anchors] {len(ANCHORS) - len(bad)}/{len(ANCHORS)} match")
    for b in bad:
        log.warning(f"    {b}")
    if bad:
        fails += 1

    # Gate 3 -- KRAB audit. This is the specific failure that motivated the work,
    # so it gets its own gate rather than being assumed away by gate 2.
    krab = done[done["dbd_excluded_domains"].fillna("").str.lower().str.contains(KRAB_HINT)]
    carried = krab[krab["dbd_review_flag"] == 1]
    log.info(f"[gate 3/KRAB] {len(krab):,} records dropped a KRAB domain; "
             f"{len(carried):,} still have one inside the kept span")
    if len(carried):
        fails += 1
        for _, r in carried.head(10).iterrows():
            log.warning(f"    {r['gene_symbol']} ({r['uniprot_id']}): "
                        f"{r['dbd_start']}-{r['dbd_end']}")

    log.info(f"[validate] {fails} gate(s) failed")
    return fails


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=str(ROOT / "config" / "sources.yaml"))
    ap.add_argument("--db", default=str(DB_PATH))
    ap.add_argument("--survey", action="store_true",
                    help="print the DOMAIN vocabulary and exit (no writes)")
    ap.add_argument("--dry-run", action="store_true",
                    help="extract and report, but do not write the library")
    ap.add_argument("--validate", action="store_true",
                    help="run the three acceptance gates against the library")
    args = ap.parse_args()

    log = get_logger("06_extract_dbd", LOGS_DIR)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    if args.validate:
        with ModuleLibrary(args.db) as lib:
            sys.exit(1 if validate(lib, log) else 0)

    cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
    api_base = cfg["uniprot"]["api_base"]

    with ModuleLibrary(args.db) as lib:
        df = lib.to_dataframe("DBD")

    total = len(df)
    df = df[df["uniprot_id"].notna() & (df["uniprot_id"].astype(str).str.strip() != "")
            & df["sequence_aa"].notna() & (df["sequence_aa"].astype(str).str.len() > 0)]
    rows = df.to_dict("records")
    log.info(f"[library] {total:,} DBD rows, {len(rows):,} extractable "
             f"(the rest have no accession or no sequence and are left untouched)")

    session = requests.Session()
    session.headers["User-Agent"] = "eCR_mod_lib/06_extract_dbd"
    accs = sorted({str(r["uniprot_id"]).strip() for r in rows})
    features = fetch_features(accs, session, api_base, log)

    if args.survey:
        survey(rows, features, log)
        return

    updates, report, kinds = [], [], Counter()
    by_family = defaultdict(Counter)
    for r in rows:
        acc = str(r["uniprot_id"]).strip()
        seq = str(r["sequence_aa"])
        line = features.get(acc)
        if line is None:
            res = {"dbd_extraction": "none", "dbd_start": None, "dbd_end": None,
                   "dbd_sequence_aa": None, "dbd_length_aa": None,
                   "dbd_excluded_domains": None, "dbd_review_flag": 0,
                   "n_znf": 0, "_inside": []}
        else:
            res = extract_one(seq, parse_features(line))

        kinds[res["dbd_extraction"]] += 1
        by_family[r["subtype"] or "?"][res["dbd_extraction"]] += 1

        updates.append((
            seq, res["dbd_start"], res["dbd_end"], res["dbd_sequence_aa"],
            res["dbd_length_aa"], res["dbd_extraction"],
            res["dbd_excluded_domains"], res["dbd_review_flag"],
            datetime.now(timezone.utc).isoformat(timespec="seconds"),
            r["module_id"],
        ))
        report.append({
            "module_id": r["module_id"], "gene_symbol": r["gene_symbol"],
            "uniprot_id": acc, "family": r["subtype"], "organism": r["organism"],
            "full_length_aa": len(seq), "dbd_start": res["dbd_start"],
            "dbd_end": res["dbd_end"], "dbd_length_aa": res["dbd_length_aa"],
            "fraction": (round(res["dbd_length_aa"] / len(seq), 3)
                         if res["dbd_length_aa"] else None),
            "dbd_extraction": res["dbd_extraction"], "n_znf": res["n_znf"],
            "zinc_finger_count": r["zinc_finger_count"],
            "excluded_domains": res["dbd_excluded_domains"],
            "review_flag": res["dbd_review_flag"],
            "inside_kept_span": ";".join(res.get("_inside", [])) or None,
        })

    log.info(f"[extract] by kind: {dict(kinds)}")
    log.info("[extract] coverage by family (extracted / total):")
    for fam, c in sorted(by_family.items(), key=lambda kv: -sum(kv[1].values()))[:20]:
        tot = sum(c.values())
        log.info(f"    {tot - c['none']:5d} / {tot:5d}  {fam}")
    flagged = sum(1 for x in report if x["review_flag"])
    log.info(f"[extract] {flagged:,} records need review "
             f"(an excluded domain or effector region lies inside the kept span)")

    import pandas as pd
    out = PROCESSED_DIR / "dbd_extraction_report.tsv"
    pd.DataFrame(report).to_csv(out, sep="\t", index=False)
    log.info(f"[extract] wrote {out}")

    if args.dry_run:
        log.info("[dry-run] library not modified")
        return

    with ModuleLibrary(args.db) as lib:
        with lib.transaction() as cur:
            cur.executemany(
                "UPDATE modules SET sequence_full_aa=?, dbd_start=?, dbd_end=?, "
                "dbd_sequence_aa=?, dbd_length_aa=?, dbd_extraction=?, "
                "dbd_excluded_domains=?, dbd_review_flag=?, date_modified=? "
                "WHERE module_id=?", updates)
        lib.insert_provenance({
            "run_id": str(uuid.uuid4()),
            "script": "06_extract_dbd.py",
            "run_date": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "source_name": "UniProt feature table",
            "source_url": api_base,
            "source_version": cfg["uniprot"]["version"],
            "file_path": str(out),
            "records_updated": len(updates),
            "notes": (f"DBD span extraction, FLANK={FLANK}; kinds={dict(kinds)}; "
                      f"{flagged} flagged for review"),
        })
    log.info(f"[extract] updated {len(updates):,} records in {args.db}")
    log.info("Next: python scripts/06_extract_dbd.py --validate")


if __name__ == "__main__":
    main()
