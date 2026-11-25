"""
scripts/07_validate_dbd_external.py
-----------------------------------
Validate the extracted DBD spans against evidence from OUTSIDE this pipeline.

Why
---
`06_extract_dbd.py --validate` has three acceptance gates, and none of them
reaches outside the pipeline's own inputs:

  gate 1  span length against the kept-finger count -- re-read from the same
          cached UniProt feature table the extraction used
  gate 2  seven hand-checked ANCHORS from the fib->beta run (2026-09-03), so it
          shows the script reproduces the by-hand method, not that either is
          right
  gate 3  KRAB spans, again from the same cached feature table

Gates 1 and 3 do test something real, but a reviewer can say -- correctly --
that the extraction is being checked against its own source. This script
answers that objection with two independent arms.

Arm 1: Pfam via InterPro
------------------------
Pfam domain boundaries come from HMM matches against the sequence. UniProt's
`DNA_BIND` lines on these records are largely PROSITE-ProRule annotations
(Nanog's reads `ECO:0000255|PROSITE-ProRule:PRU00108`), so Pfam is a different
method on the same sequence rather than a copy of the same curation.

Reported as a DISTRIBUTION, not a pass/fail gate, because the two disagree
legitimately by design: Pfam matches each C2H2 finger individually, while this
extraction keeps a finger array contiguous because the inter-finger TGEKP
linkers are part of the array. So a `znf_array` span SHOULD contain residues no
Pfam domain covers, and precision below 1.0 is expected for that kind and not a
defect. Recall -- how much of the Pfam DNA-binding domain the span contains --
is the number that must be high for every kind.

Arm 1 also runs the test the gates cannot: for the records where the extraction
ABSTAINED (`dbd_extraction = 'none'`), does Pfam see a DNA-binding domain
anyway? Those would be chassis the method is losing.

Arm 2: DNA contacts in solved complexes
---------------------------------------
For chassis with a protein-DNA complex in the PDB, compute which protein
residues lie within CONTACT_CUTOFF_A of any DNA atom, map them to UniProt
numbering through SIFTS, and ask what fraction fall inside the extracted span.
This is the strongest arm -- experimental coordinates, no annotation in the
loop -- and the narrowest, since it only covers genes anyone has crystallised
with DNA.

Run
---
  python scripts/07_validate_dbd_external.py --fetch      # InterPro, cached per accession
  python scripts/07_validate_dbd_external.py --survey     # Pfam name frequencies (no verdicts)
  python scripts/07_validate_dbd_external.py --report     # arm 1 distribution
  python scripts/07_validate_dbd_external.py --pdb [--limit N]   # arm 2

The per-record tables land in data/processed/ (gitignored); the committed
deliverable is docs/DBD_EXTERNAL_VALIDATION.md, written by hand from them
because the interesting part is which disagreements are legitimate.

Cost: arm 1 is ~2,803 cached GETs against the InterPro REST API, minutes. Arm 2
downloads one mmCIF per selected structure.
"""

from __future__ import annotations

import argparse
import io
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent))
from utils import get_logger

ROOT = Path(__file__).resolve().parents[1]
LOGS_DIR = ROOT / "logs"
LIBRARY_TSV = ROOT / "library" / "module_library.tsv"
INTERPRO_DIR = ROOT / "data" / "raw" / "interpro"
PDB_DIR = ROOT / "data" / "raw" / "pdb"
PROCESSED_DIR = ROOT / "data" / "processed"

INTERPRO_API = "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot"
PDBE_BEST = "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures"
PDBE_MOLECULES = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules"
PDBE_SIFTS = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot"
CIF_URL = "https://files.rcsb.org/download/{pdb}.cif"

# Which Pfam families count as the DNA-binding apparatus: scripts/pfam_dbd.py.
#
# That list used to live here. It moved when 06_extract_dbd.py started using the
# same answer for its Pfam fallback (F5 recommendation 2) -- the lists must not
# drift apart, because arm 1 below would then be scoring spans against a
# different definition of "DNA-binding" than the extraction used to build them.
#
# Recall is measured against a DNA-binding family rather than against every
# Pfam hit, because a span that correctly excludes a KRAB or a nuclear-receptor
# ligand-binding domain must not be scored as having missed it.
from pfam_dbd import PFAM_DBD_GENERIC, PFAM_DBD_HINTS  # noqa: E402,F401
from pfam_dbd import is_dbd_family as _is_dbd_family   # noqa: E402

CONTACT_CUTOFF_A = 4.0      # heavy-atom distance for a protein-DNA contact
DNA_RESIDUES = {"DA", "DC", "DG", "DT", "DU", "DI"}


# ── shared ────────────────────────────────────────────────────────────────

def _library() -> pd.DataFrame:
    df = pd.read_csv(LIBRARY_TSV, sep="\t", dtype=str, low_memory=False)
    df = df[df["type"] == "DBD"].copy()
    df["uniprot_id"] = df["uniprot_id"].fillna("").str.strip()
    return df[df["uniprot_id"] != ""]


def _session(agent: str) -> requests.Session:
    s = requests.Session()
    s.headers["User-Agent"] = agent
    s.headers["Accept"] = "application/json"
    return s


def _span(row) -> "tuple[int, int] | None":
    s, e = pd.to_numeric(row["dbd_start"], errors="coerce"), pd.to_numeric(row["dbd_end"], errors="coerce")
    if pd.isna(s) or pd.isna(e):
        return None
    return int(s), int(e)


# ── arm 1: fetch ─────────────────────────────────────────────────────────────

def fetch(log) -> None:
    INTERPRO_DIR.mkdir(parents=True, exist_ok=True)
    accs = sorted(set(_library()["uniprot_id"]))
    todo = [a for a in accs if not (INTERPRO_DIR / f"{a}.json").exists()]
    log.info(f"[interpro] {len(accs):,} accessions, {len(accs) - len(todo):,} cached, "
             f"{len(todo):,} to fetch")
    sess = _session("eCR_mod_lib/07_validate_dbd_external")
    fail = 0
    for i, acc in enumerate(todo, 1):
        try:
            r = sess.get(f"{INTERPRO_API}/{acc}/", timeout=60)
            if r.status_code == 204 or r.status_code == 404:
                payload = {"accession": acc, "domains": []}
            else:
                r.raise_for_status()
                data = r.json()
                domains = []
                for res in data.get("results", []):
                    md = res.get("metadata", {})
                    frags = []
                    for prot in res.get("proteins", []):
                        for loc in prot.get("entry_protein_locations", []):
                            for f in loc.get("fragments", []):
                                frags.append([int(f["start"]), int(f["end"])])
                    domains.append({"pfam": md.get("accession"),
                                    "name": md.get("name"),
                                    "fragments": frags})
                payload = {"accession": acc, "domains": domains}
            (INTERPRO_DIR / f"{acc}.json").write_text(json.dumps(payload),
                                                      encoding="utf-8")
        except Exception as exc:                      # noqa: BLE001
            fail += 1
            log.warning(f"[interpro] {acc}: {exc}")
        if i % 200 == 0:
            log.info(f"[interpro] {i:,}/{len(todo):,} fetched ({fail} failed)")
        time.sleep(0.1)
    log.info(f"[interpro] done, {fail} failures; cache in {INTERPRO_DIR}")


def _pfam_for(acc: str) -> "list[dict] | None":
    p = INTERPRO_DIR / f"{acc}.json"
    if not p.exists():
        return None
    return json.loads(p.read_text(encoding="utf-8")).get("domains", [])


def _residues(frags) -> set:
    out = set()
    for s, e in frags:
        out.update(range(s, e + 1))
    return out


# ── arm 1: survey ────────────────────────────────────────────────────────────

def survey(log) -> None:
    """Pfam name frequencies, split by whether they overlap an extracted span.

    This is what PFAM_DBD_HINTS was built from: a family that consistently
    overlaps spans is part of the DNA-binding apparatus as this method defines
    it, and one that never does is what the extraction drops."""
    df = _library()
    overlapping, disjoint, missing = Counter(), Counter(), 0
    for _, row in df.iterrows():
        doms = _pfam_for(row["uniprot_id"])
        if doms is None:
            missing += 1
            continue
        sp = _span(row)
        for d in doms:
            res = _residues(d["fragments"])
            if sp and res & set(range(sp[0], sp[1] + 1)):
                overlapping[d["name"]] += 1
            else:
                disjoint[d["name"]] += 1
    log.info(f"[survey] {missing:,} records with no cached InterPro entry")
    log.info("[survey] Pfam families OVERLAPPING an extracted span (top 40):")
    for name, n in overlapping.most_common(40):
        mark = "DBD " if _is_dbd_family(name) else "    "
        log.info(f"    {mark}{n:5d}  {name}")
    log.info("[survey] Pfam families never overlapping a span (top 25):")
    for name, n in disjoint.most_common(25):
        mark = "DBD?" if _is_dbd_family(name) else "    "
        log.info(f"    {mark}{n:5d}  {name}")


# ── arm 1: report ────────────────────────────────────────────────────────────

def _agreement_rows(df) -> "list[dict]":
    """One row per record.

    Recall is scored against the SINGLE Pfam DNA-binding domain the span
    overlaps most, not against the union of all of them. A protein can carry
    two DNA-binding-family domains in different places -- ZSCAN20 has a
    Myb/SANT-like domain hundreds of residues from its C2H2 array -- and
    scoring against the union would mark a correct span as having missed a
    module it was never meant to contain. The omitted ones are counted
    separately instead, which is the honest way to say "the span covers one of
    the two".
    """
    rows = []
    for _, row in df.iterrows():
        doms = _pfam_for(row["uniprot_id"])
        if doms is None:
            continue
        dbd_doms, all_res = [], set()
        for d in doms:
            res = _residues(d["fragments"])
            all_res |= res
            if _is_dbd_family(d["name"]):
                dbd_doms.append((d["name"], res))
        sp = _span(row)
        rec = {
            "module_id": row["module_id"],
            "gene_symbol": row["gene_symbol"],
            "organism": row["organism"],
            "family": row["subtype"],
            "kind": row["dbd_extraction"],
            "review_flag": row["dbd_review_flag"],
            "span_start": sp[0] if sp else None,
            "span_end": sp[1] if sp else None,
            "n_pfam_dbd_domains": len(dbd_doms),
            "pfam_dbd_domains": ";".join(sorted({n for n, _ in dbd_doms})) or None,
            "n_pfam_any_residues": len(all_res),
        }
        if not sp:
            rec.update({"span_len": None, "best_domain": None,
                        "recall_best": None, "precision_span": None,
                        "jaccard": None, "start_offset": None,
                        "end_offset": None, "n_dbd_domains_outside": len(dbd_doms),
                        "verdict": "no_span"})
            rows.append(rec)
            continue

        span_res = set(range(sp[0], sp[1] + 1))
        overlapping = [(len(span_res & res), name, res) for name, res in dbd_doms
                       if span_res & res]
        outside = len(dbd_doms) - len(overlapping)
        if not dbd_doms:
            verdict = "no_pfam_dbd_domain"
        elif not overlapping:
            verdict = "disjoint_from_every_pfam_dbd_domain"
        else:
            verdict = "overlaps"
        if overlapping:
            _n, name, res = max(overlapping)
            inter = span_res & res
            rec.update({
                "span_len": len(span_res),
                "best_domain": name,
                "recall_best": len(inter) / len(res),
                "precision_span": len(inter) / len(span_res),
                "jaccard": len(inter) / len(span_res | res),
                "start_offset": sp[0] - min(res),
                "end_offset": sp[1] - max(res),
                "n_dbd_domains_outside": outside,
                "verdict": verdict,
            })
        else:
            rec.update({"span_len": len(span_res), "best_domain": None,
                        "recall_best": None, "precision_span": None,
                        "jaccard": None, "start_offset": None,
                        "end_offset": None,
                        "n_dbd_domains_outside": outside, "verdict": verdict})
        rows.append(rec)
    return rows


def _q(series) -> str:
    s = pd.to_numeric(series, errors="coerce").dropna()
    if s.empty:
        return "n/a"
    return (f"median {s.median():.3f}  p25 {s.quantile(.25):.3f}  "
            f"p75 {s.quantile(.75):.3f}  min {s.min():.3f}")


def report(log) -> pd.DataFrame:
    df = _library()
    rows = _agreement_rows(df)
    out = pd.DataFrame(rows)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    dest = PROCESSED_DIR / "dbd_external_validation.tsv"
    out.to_csv(dest, sep="\t", index=False)
    log.info(f"[report] {len(out):,} records compared; wrote {dest}")

    log.info(f"[report] verdicts: {dict(out['verdict'].value_counts())}")
    withspan = out[out["recall_best"].notna()]
    log.info(f"[report] {len(withspan):,} records have a span overlapping a Pfam "
             f"DNA-binding domain")
    log.info(f"[report] recall   (best Pfam DBD domain inside the span): {_q(withspan['recall_best'])}")
    log.info(f"[report] precision(span residues in that domain)        : {_q(withspan['precision_span'])}")
    log.info(f"[report] jaccard                                        : {_q(withspan['jaccard'])}")
    log.info(f"[report] start offset (span - domain, aa): {_q(withspan['start_offset'])}")
    log.info(f"[report] end offset   (span - domain, aa): {_q(withspan['end_offset'])}")

    log.info("[report] by extraction kind:")
    for kind, g in withspan.groupby("kind"):
        r = pd.to_numeric(g["recall_best"], errors="coerce")
        p = pd.to_numeric(g["precision_span"], errors="coerce")
        log.info(f"    {kind:10s} n={len(g):5d}  recall median {r.median():.3f} "
                 f"(>=0.95 in {100 * (r >= .95).mean():.0f}%)  "
                 f"precision median {p.median():.3f}")

    log.info("[report] recall by family (top 15 by n):")
    for family, g in sorted(withspan.groupby("family"), key=lambda kv: -len(kv[1]))[:15]:
        r = pd.to_numeric(g["recall_best"], errors="coerce")
        log.info(f"    {str(family)[:18]:18s} n={len(g):5d}  recall median {r.median():.3f}  "
                 f"worst {r.min():.3f}")

    # Spans that overlap NO Pfam DNA-binding domain: the serious class.
    disjoint = out[out["verdict"] == "disjoint_from_every_pfam_dbd_domain"]
    log.info(f"[report] {len(disjoint):,} spans overlap NO Pfam DNA-binding domain:")
    for _, r in disjoint.head(20).iterrows():
        log.info(f"    {r['gene_symbol']:12s} {r['kind']:10s} "
                 f"span {int(r['span_start'])}-{int(r['span_end'])}  "
                 f"flag {r['review_flag']}  {r['pfam_dbd_domains']}")

    nopfam = out[out["verdict"] == "no_pfam_dbd_domain"]
    log.info(f"[report] {len(nopfam):,} spans where Pfam names no DNA-binding "
             f"family at all (Pfam coverage gap, or a span with no basis)")

    # The test the gates cannot run: abstentions Pfam disagrees with.
    abstained = out[out["verdict"] == "no_span"]
    missed = abstained[abstained["n_pfam_dbd_domains"] > 0]
    log.info(f"[report] {len(abstained):,} records without a span; Pfam sees a "
             f"DNA-binding domain in {len(missed):,} of them")
    for _, r in missed.head(20).iterrows():
        log.info(f"    {r['gene_symbol']:12s} {str(r['kind']):10s} "
                 f"{r['pfam_dbd_domains']}")

    log.info("[report] lowest recall among overlapping spans:")
    for _, r in withspan.nsmallest(15, "recall_best").iterrows():
        log.info(f"    {r['gene_symbol']:12s} {r['kind']:10s} "
                 f"recall {float(r['recall_best']):.2f}  "
                 f"span {int(r['span_start'])}-{int(r['span_end'])}  "
                 f"flag {r['review_flag']}  {r['best_domain']}")
    return out


# ── arm 2: DNA contacts in the PDB ───────────────────────────────────────────

def _cif_atoms(text: str) -> "list[tuple]":
    """(auth_asym, auth_seq, label_asym, label_seq, comp_id, x, y, z) per atom.

    Both numbering schemes are kept because SIFTS populates them
    inconsistently: 1DSZ maps RARA with `author_residue_number: null` and only
    the label numbering filled in, and requiring author numbers silently
    dropped every such structure."""
    lines = text.splitlines()
    i = 0
    atoms = []
    while i < len(lines):
        if lines[i].strip() == "loop_":
            cols, j = [], i + 1
            while j < len(lines) and lines[j].lstrip().startswith("_"):
                cols.append(lines[j].strip())
                j += 1
            if any(c.startswith("_atom_site.") for c in cols):
                idx = {c.split(".", 1)[1]: k for k, c in enumerate(cols)}
                need = ("auth_asym_id", "auth_seq_id", "label_asym_id",
                        "label_seq_id", "label_comp_id",
                        "Cartn_x", "Cartn_y", "Cartn_z")
                if all(n in idx for n in need):
                    while j < len(lines) and not lines[j].lstrip().startswith(("loop_", "_", "#")):
                        parts = lines[j].split()
                        if len(parts) >= len(cols):
                            try:
                                lseq = parts[idx["label_seq_id"]]
                                atoms.append((
                                    parts[idx["auth_asym_id"]],
                                    int(parts[idx["auth_seq_id"]]),
                                    parts[idx["label_asym_id"]],
                                    int(lseq) if lseq.isdigit() else None,
                                    parts[idx["label_comp_id"]],
                                    float(parts[idx["Cartn_x"]]),
                                    float(parts[idx["Cartn_y"]]),
                                    float(parts[idx["Cartn_z"]]),
                                ))
                            except (ValueError, KeyError):
                                pass
                        j += 1
                i = j
                continue
            i = j
            continue
        i += 1
    return atoms


def _contacts(atoms) -> set:
    """{(auth_asym, auth_seq, label_asym, label_seq)} for protein residues
    within CONTACT_CUTOFF_A of any DNA atom."""
    dna = [a for a in atoms if a[4] in DNA_RESIDUES]
    prot = [a for a in atoms if a[4] not in DNA_RESIDUES and len(a[4]) == 3]
    if not dna or not prot:
        return set()
    cut2 = CONTACT_CUTOFF_A ** 2
    # Bucket DNA atoms on a 4 A grid so this stays linear rather than N*M.
    grid = defaultdict(list)
    for _ca, _cs, _la, _ls, _n, x, y, z in dna:
        grid[(int(x // CONTACT_CUTOFF_A), int(y // CONTACT_CUTOFF_A),
              int(z // CONTACT_CUTOFF_A))].append((x, y, z))
    hits = set()
    for chain, seq, lasym, lseq, _n, x, y, z in prot:
        gx, gy, gz = (int(x // CONTACT_CUTOFF_A), int(y // CONTACT_CUTOFF_A),
                      int(z // CONTACT_CUTOFF_A))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for (px, py, pz) in grid.get((gx + dx, gy + dy, gz + dz), ()):
                        if (x - px) ** 2 + (y - py) ** 2 + (z - pz) ** 2 <= cut2:
                            hits.add((chain, seq, lasym, lseq))
                            break
                    else:
                        continue
                    break
                else:
                    continue
                break
    return hits


def _pdb_candidates(accs, sess, log) -> dict:
    """{accession: (pdb_id, resolution)} for accessions with a protein-DNA entry.

    Batched, because the naive route -- ask PDBe for every accession's best
    structures, then ask for each entry's molecule list -- is thousands of
    round trips and most accessions have no structure at all. Instead: one
    UniProt search per 100 accessions for the PDB cross-references it already
    holds, then one RCSB GraphQL query per 200 entries for the DNA chain count.
    Cached, so a re-run costs nothing."""
    cache = PDB_DIR / "_candidates.json"
    if cache.exists():
        return {k: tuple(v) for k, v in
                json.loads(cache.read_text(encoding="utf-8")).items()}

    xrefs = {}
    for i in range(0, len(accs), 100):
        batch = accs[i:i + 100]
        q = " OR ".join(f"accession:{a}" for a in batch)
        r = sess.get("https://rest.uniprot.org/uniprotkb/search",
                     params={"query": q, "fields": "accession,xref_pdb",
                             "format": "tsv", "size": 500}, timeout=120)
        r.raise_for_status()
        for line in r.text.splitlines()[1:]:
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            ids = [p for p in parts[1].split(";") if p]
            if ids:
                xrefs[parts[0]] = ids
        log.info(f"[pdb] UniProt xrefs {min(i + 100, len(accs))}/{len(accs)}")
        time.sleep(0.1)
    log.info(f"[pdb] {len(xrefs):,} accessions carry a PDB cross-reference")

    entries = sorted({e for ids in xrefs.values() for e in ids})
    info = {}
    for i in range(0, len(entries), 200):
        batch = entries[i:i + 200]
        ids = ",".join(f'"{e}"' for e in batch)
        query = ("{entries(entry_ids:[" + ids + "]){rcsb_id rcsb_entry_info"
                 "{polymer_entity_count_DNA resolution_combined}}}")
        try:
            r = sess.post("https://data.rcsb.org/graphql", json={"query": query},
                          timeout=180)
            r.raise_for_status()
            for ent in r.json().get("data", {}).get("entries", []) or []:
                if not ent:
                    continue
                ei = ent.get("rcsb_entry_info") or {}
                res = (ei.get("resolution_combined") or [None])[0]
                info[ent["rcsb_id"]] = (ei.get("polymer_entity_count_DNA") or 0,
                                        res if res is not None else 99.0)
        except Exception as exc:                      # noqa: BLE001
            log.warning(f"[pdb] graphql batch at {i}: {exc}")
        log.info(f"[pdb] entry info {min(i + 200, len(entries))}/{len(entries)}")
        time.sleep(0.1)

    out = {}
    for acc, ids in xrefs.items():
        withdna = [(info[e][1], e) for e in ids
                   if e in info and info[e][0] > 0]
        if withdna:
            res, pdb_id = sorted(withdna)[0]      # best resolution, id as tie-break
            out[acc] = (pdb_id, res)
    PDB_DIR.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps({k: list(v) for k, v in out.items()}),
                     encoding="utf-8")
    log.info(f"[pdb] {len(out):,} accessions have a protein-DNA structure")
    return out


def pdb_arm(log, limit: int) -> None:
    PDB_DIR.mkdir(parents=True, exist_ok=True)
    df = _library()
    df = df[df["dbd_sequence_aa"].notna()]
    accs = sorted(set(df["uniprot_id"]))
    sess = _session("eCR_mod_lib/07_validate_dbd_external")

    candidates = _pdb_candidates(accs, sess, log)
    # Deterministic, reproducible selection: best-resolution structures first,
    # accession as the tie-break, capped at `limit`.
    chosen = sorted(((res, acc, pdb) for acc, (pdb, res) in candidates.items()
                     if acc in set(accs)))[:limit]
    log.info(f"[pdb] folding contacts for {len(chosen)} structures "
             f"(of {len(candidates)} available)")

    results = []
    for res, acc, pdb_id in chosen:
        cif_path = PDB_DIR / f"{pdb_id}.cif"
        if not cif_path.exists():
            try:
                cif = sess.get(CIF_URL.format(pdb=pdb_id), timeout=180,
                               headers={"Accept": "text/plain"})
                cif.raise_for_status()
                cif_path.write_text(cif.text, encoding="utf-8")
            except Exception as exc:                  # noqa: BLE001
                log.warning(f"[pdb] {pdb_id}: {exc}")
                continue
        try:
            sifts = sess.get(f"{PDBE_SIFTS}/{pdb_id}", timeout=60).json()
        except Exception as exc:                      # noqa: BLE001
            log.warning(f"[pdb] {pdb_id} sifts: {exc}")
            continue
        mappings = (sifts.get(pdb_id.lower(), {}).get("UniProt", {})
                    .get(acc, {}).get("mappings", []))
        if not mappings:
            log.info(f"[pdb] {pdb_id}: no SIFTS mapping for {acc}")
            continue

        contacts = _contacts(_cif_atoms(cif_path.read_text(encoding="utf-8")))
        if not contacts:
            log.info(f"[pdb] {pdb_id}: no protein-DNA contacts found")
            continue

        # Author numbering first; label numbering where SIFTS leaves the author
        # fields null (1DSZ/RARA does, and every such structure used to drop out).
        unp_contacts = set()
        for m in mappings:
            unp_start = m.get("unp_start")
            if unp_start is None:
                continue
            a_chain = m.get("chain_id")
            a_start = (m.get("start") or {}).get("author_residue_number")
            a_end = (m.get("end") or {}).get("author_residue_number")
            if None not in (a_chain, a_start, a_end):
                for chain, seq, _la, _ls in contacts:
                    if chain == a_chain and a_start <= seq <= a_end:
                        unp_contacts.add(unp_start + (seq - a_start))
                continue
            l_chain = m.get("struct_asym_id")
            l_start = (m.get("start") or {}).get("residue_number")
            l_end = (m.get("end") or {}).get("residue_number")
            if None in (l_chain, l_start, l_end):
                continue
            for _ca, _seq, lasym, lseq in contacts:
                if lasym == l_chain and lseq is not None and l_start <= lseq <= l_end:
                    unp_contacts.add(unp_start + (lseq - l_start))
        if not unp_contacts:
            log.info(f"[pdb] {pdb_id}: contacts outside the mapped segments "
                     f"for {acc} (another chain binds the DNA, or the "
                     f"DNA-binding region is unresolved)")
            continue

        for _, row in df[df["uniprot_id"] == acc].iterrows():
            sp = _span(row)
            if not sp:
                continue
            inside = {u for u in unp_contacts if sp[0] <= u <= sp[1]}
            results.append({
                "module_id": row["module_id"], "gene_symbol": row["gene_symbol"],
                "uniprot_id": acc, "pdb_id": pdb_id, "resolution_a": res,
                "kind": row["dbd_extraction"], "span": f"{sp[0]}-{sp[1]}",
                "n_contacts": len(unp_contacts),
                "n_contacts_in_span": len(inside),
                "contact_coverage": len(inside) / len(unp_contacts),
                "missed": ",".join(str(u) for u in sorted(unp_contacts - inside)[:12]),
            })
        log.info(f"[pdb] {acc} / {pdb_id} ({res} A): {len(unp_contacts)} "
                 f"DNA-contacting residues")

    if not results:
        log.warning("[pdb] no structures resolved")
        return
    out = pd.DataFrame(results)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    dest = PROCESSED_DIR / "dbd_pdb_contact_validation.tsv"
    out.to_csv(dest, sep="\t", index=False)
    cov = pd.to_numeric(out["contact_coverage"], errors="coerce")
    log.info(f"[pdb] {len(out):,} records against {out['pdb_id'].nunique()} "
             f"structures; contact coverage {_q(cov)}")
    log.info(f"[pdb] fully covered (1.00): {int((cov == 1).sum())}/{len(cov)}")
    for _, r in out.nsmallest(10, "contact_coverage").iterrows():
        log.info(f"    {r['gene_symbol']:12s} {r['pdb_id']}  span {r['span']}  "
                 f"{r['n_contacts_in_span']}/{r['n_contacts']} contacts  "
                 f"missed {r['missed']}")
    log.info(f"[pdb] wrote {dest}")



def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fetch", action="store_true", help="fetch InterPro/Pfam locations")
    ap.add_argument("--survey", action="store_true", help="Pfam family frequencies")
    ap.add_argument("--report", action="store_true", help="arm 1 agreement distribution")
    ap.add_argument("--pdb", action="store_true", help="arm 2 DNA-contact coverage")
    ap.add_argument("--limit", type=int, default=40,
                    help="max structures for --pdb (default 40)")
    args = ap.parse_args()

    log = get_logger("07_validate_dbd_external", LOGS_DIR)
    if not (args.fetch or args.survey or args.report or args.pdb):
        ap.error("pick at least one of --fetch / --survey / --report / --pdb")
    if args.fetch:
        fetch(log)
    if args.survey:
        survey(log)
    if args.report:
        report(log)
    if args.pdb:
        pdb_arm(log, args.limit)


if __name__ == "__main__":
    main()
