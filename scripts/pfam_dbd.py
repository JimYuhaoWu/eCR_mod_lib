"""
scripts/pfam_dbd.py
-------------------
Which Pfam families are DNA-binding, and where they sit on an accession.

Two scripts need the same answer and must not drift apart:

  07_validate_dbd_external.py  scores a span against the best-overlapping Pfam
                               DNA-binding domain (arm 1 of the external
                               validation)
  06_extract_dbd.py            falls back to the Pfam envelope for records
                               where the UniProt feature table names no
                               keep-feature at all (F5 recommendation 2)

The lists were built by classifying every one of the 513 distinct family names
the InterPro fetch returns for this library -- NOT from Pfam short names. The
API returns DESCRIPTIVE names ("Zinc finger, C2H2 type", not "zf-C2H2"), and a
list keyed on short names silently matched 0 of 722 C2H2 hits.
"""
from __future__ import annotations

import json
from pathlib import Path

CACHE_DIR = Path(__file__).resolve().parents[1] / "data" / "raw" / "interpro"

PFAM_DBD_GENERIC = ("dna-binding", "dna binding", "dnabinding")
PFAM_DBD_HINTS = (
    # zinc fingers that contact DNA
    "zinc finger, c2h2", "c2h2-type zinc finger", "c2h2 zinc finger",
    "zinc-finger of c2h2", "zinc-finger double domain",
    "double treble clef zinc finger", "gata zinc finger",
    "zic proteins zinc finger",
    # helical / winged-helix folds
    "homeodomain", "homeobox", "pou domain", "forkhead domain",
    "forkhead n-terminal", "t-box", "e2f/dp family winged-helix",
    "rel homology", "runt domain", "sand domain", "cut domain",
    "thap domain", "'paired box' domain", "ap2 domain",
    # basic / leucine-zipper and bHLH
    "bzip transcription factor", "bzip maf transcription factor",
    "basic region leucine zipper", "myogenic basic domain",
    "hypoxia-inducible factor 1-alpha bhlh",
    # groove binders and others
    "hmg (high mobility group) box", "hmg-box domain", "mh1 domain",
    "methyl-cpg binding domain", "cp2 transcription factor",
    "srf-type transcription factor", "ets-domain",
    "interferon regulatory factor transcription factor",
    "tale/meis homeobox",
)
# Deliberately NOT DNA-binding, though they sit in these proteins and would be
# caught by a looser rule: KRAB box, SCAN, BTB/POZ, ELM2, LIM, SAM/Pointed,
# PAS, MH2, OAR motif, PR/SET, and the nuclear-receptor ligand-binding domain.
# Also excluded: "SOX transcription factor" (PF12336 is SOX's C-terminal
# region, not its HMG box) and "GRHL1-like, C-terminal domain".


def is_dbd_family(name: str) -> bool:
    low = (name or "").lower()
    if "grhl1-like" in low or "sox transcription factor" in low:
        return False
    return (any(h in low for h in PFAM_DBD_GENERIC)
            or any(h in low for h in PFAM_DBD_HINTS))


def dbd_domains(acc: str, cache_dir: Path | None = None
                ) -> "tuple[list[tuple[int, int]], list[str]] | None":
    """The Pfam DNA-binding domains of `acc` as [(start, end), ...] plus their
    family names, or None when the cache has no entry or no DNA-binding family.

    Fragments are returned individually rather than as one envelope, because
    the caller needs both: min(start)..max(end) for the span, and the separate
    fragments to measure how much of that span is actually domain -- the same
    distinction `_merge` draws on the UniProt side.

    Returns None rather than raising when the cache is absent. The Pfam
    fallback is an improvement where the cache exists and must not become a
    hard dependency of the extraction where it does not.
    """
    p = (cache_dir or CACHE_DIR) / f"{acc}.json"
    if not p.exists():
        return None
    try:
        doms = json.loads(p.read_text(encoding="utf-8")).get("domains", [])
    except (OSError, ValueError):
        return None
    hits = [d for d in doms if is_dbd_family(d.get("name", ""))]
    frs = [(int(fr[0]), int(fr[1])) for d in hits
           for fr in d.get("fragments", []) if len(fr) >= 2]
    if not frs:
        return None
    return frs, sorted({d["name"] for d in hits})
