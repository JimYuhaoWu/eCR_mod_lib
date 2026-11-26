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
        - On abstention only: a REGION whose note NAMES it as the
          DNA-binding module and is long enough to be one (KEEP_REGION_EXACT).
          HSF1 and TP73 file their DBD this way; reading only the first three
          feature types abstained on them while TP53 next door extracted
          cleanly from a DNA_BIND line. 9 records, all `dna_bind`.
  DROP  everything else, explicitly including KRAB, SCAN, BTB/POZ and
        interaction/activation regions.

Then span = min(start) - FLANK .. max(end) + FLANK, clamped to the sequence.

When the feature table yields NO keep-feature at all, the two fallbacks run
in order: a self-naming REGION first, then Pfam (scripts/pfam_dbd.py, cached
under data/raw/interpro/). Both are gated on abstention, so neither can
override a UniProt keep-feature, and a Pfam-derived span is marked
`dbd_extraction='pfam'`. STAT1's whole feature table is one SH2 domain, so
without this the entire STAT family abstained.

Neither fallback changes a span that already existed: measured over the 2,818
extractable records, 35 abstentions gain a span (18 Pfam, 9 REGION, 8 from
"matrin" joining ZNF_KEEP_TYPES), 18 change because 'ipt/tig' left
KEEP_DOMAIN_HINTS, and nothing else moves. FLANK is deliberately unchanged at
15. The gates do not see it: each measures the annotated feature extent, so
a larger flank changes construct boundaries without changing any verdict.

Rules that matter
-----------------
1. Contiguous span, never concatenated fragments. Splicing across an excluded
   internal domain creates a junction that exists in no real protein -- a novel
   neoepitope surface and an unpredictable fold.
2. A dropped domain lying INSIDE the kept span is a conflict, not a silent pass:
   `dbd_review_flag` is set and the offender named in `dbd_excluded_domains`.
3. Never fabricate boundaries. No usable features -> `dbd_extraction='none'` and
   the DBD columns stay NULL. A guessed domain is worse than an honest
   full-length record. Rule 7 does not weaken this: a Pfam HMM match is a
   measurement, not a guess, and it is labelled as a different source.
4. Additive. `sequence_aa` is not touched; the original is also copied to
   `sequence_full_aa` so consumers can migrate deliberately.
5. A span is a DOMAIN, not a protein. `min(start)..max(end)` over two distant
   features can sweep up almost the whole chain -- ZNF462 extracted 2,442 aa of
   2,506 -- which is not an extraction, it is a relabelling. Spans that large
   get `dbd_extraction='dispersed'` and no span columns, and implausibly short
   ones are flagged (see MAX_SPAN_FRACTION / MIN_PLAUSIBLE_SPAN_AA). Both
   gates measure the annotated feature extent, not the flanked span, so FLANK
   cannot move a verdict.
6. Say what was dropped. Every residue outside the kept span is recorded in
   `dbd_dropped_regions`, labelled from the feature table where UniProt names
   it, and `dbd_dropped_review_flag` is set when the dropped material looks
   like a self-association or repeat module. No keep rule changes: whether such
   a module belongs in a construct is a biology question, and the pipeline's
   job is only to stop dropping it silently.
7. One source per span, named. The UniProt feature table is tried first,
   then a self-naming REGION, then Pfam, each only on abstention, so a span
   never mixes two sources and `dbd_extraction` says which answered. 30
   records are 'pfam': STAT (14), EBF/COE (8), GCFC2/PAXBP1 (4), RBPJ/CSL (4).
   A promoted REGION reads as `dna_bind`, because that is the assertion
   UniProt made, on a feature line this script did not read; the
   `region_promoted` column of the report carries the provenance, and both
   sources enter data/processed/dbd_review_queue.tsv whether flagged or not,
   since neither has ever been reviewed.

Run
---
  python scripts/06_extract_dbd.py --survey     # domain vocabulary, no writes
  python scripts/06_extract_dbd.py --dry-run    # extract + report, no writes
  python scripts/06_extract_dbd.py              # extract and update the library
  python scripts/06_extract_dbd.py --validate   # the three acceptance gates

Cost: offline apart from the feature fetch (~30 batched requests); minutes, no GPU.

The Pfam fallback (rule 7) reads data/raw/interpro/, which
scripts/07_validate_dbd_external.py --fetch populates. Without that cache the
fallback is inert and those records abstain as they did before -- so a fresh
checkout should run the fetch once before the extraction.
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
from pfam_dbd import dbd_domains
from schema import ModuleLibrary
from utils import get_logger

ROOT = Path(__file__).resolve().parents[1]
LOGS_DIR = ROOT / "logs"
RAW_DIR = ROOT / "data" / "raw"
PROCESSED_DIR = ROOT / "data" / "processed"
DB_PATH = ROOT / "library" / "module_library.db"

# F5 recommendation 3. Arm 2 of the external validation showed FLANK=6 trims
# real DNA contacts: of 220 records scored against solved protein-DNA
# complexes, 30 were partial, and 93 of their 97 missed contacts lay within
# 30 aa of a span edge -- GATA3 (3DFX) misses 348-365 with a span ending at
# 347, THRB (2NLL) misses 188-198 ending at 187.
#
# Raising it was tried once and reverted, because FLANK was in the numerator
# of all three gates below: padding loosened the 40-aa floor (a 13-aa A.T hook
# cleared it on padding alone) and tightened the dispersed gate (ZNF469 lost a
# 983-aa span at 25% of its protein because 983 + 18 crossed MAX_SPAN_AA).
# None of that was a boundary decision; it was the flank being counted as
# evidence. So the gates now measure `ext`, the annotated feature extent, and
# FLANK is invisible to them -- see "Measuring the gates" below.
FLANK = 15
BATCH = 100          # accessions per UniProt search request
MIN_SPAN_COVERAGE = 0.5   # below this, the span is mostly not DNA-binding feature

# F2 -- a span that is effectively the whole protein is not an extraction.
# `min(start)..max(end)` is applied to every kind, so two distant keep-features
# bracket everything between them: ZNF462 (znf_array) came out 2,442 aa of
# 2,506, ZNF521 1,271 of 1,311 with dbd_review_flag = 0. MIN_SPAN_COVERAGE does
# not see these -- it measures feature density INSIDE the span, and a long C2H2
# array is densely packed with fingers.
#
# The test is deliberately conjunctive. Fraction alone would condemn small
# proteins that genuinely are almost all DBD and make perfectly good chassis:
# HOPX is 68 aa of 73 (93%), HMGB3 165 of 200, Zfp239 201 of 201. Length alone
# would condemn legitimate long arrays whose span is a modest share of the
# chain. So: a big span that is also most of the protein, OR a span so long it
# cannot be a domain at any fraction.
#
# Measured over the 2,556 records with a span: 117 become `dispersed` (4.6%),
# 42 of them carrying dbd_review_flag = 0 today. It catches every case this
# repo's audit named -- ZNF462, ZFHX3/4, MGA, ZNF521, ZNF423 -- and spares
# HOPX, HMGB3/4, ZNF146, ZNF664 and Zfp239.
MAX_SPAN_FRACTION = 0.80      # span / protein above this is "the whole protein"
MIN_DISPERSED_SPAN_AA = 500   # ... but only when the span is also this long
MAX_SPAN_AA = 1000            # a span this long is not a domain at any fraction

# F2b -- the other end. 27 spans are under 40 aa, 23 of them unflagged: KMT2C
# extracts 25 aa from a 4,911-aa protein, AHCTF1 25 of 2,266. A C2H2 finger is
# ~28 aa and a homeodomain ~60, so below 40 aa the span cannot hold one
# complete DNA-binding module. Flagged, not dropped: the span may still be the
# right answer for a fragment-level question, and F1 in eCR_predictor refuses
# flagged spans by default anyway.
MIN_PLAUSIBLE_SPAN_AA = 40

# Extraction verdicts that carry NO span. `none` is a rule-3 abstention,
# `no_features` is a fetch that returned nothing (F9), `dispersed` is F2's
# verdict above. Consumers must test membership here rather than `!= 'none'`.
SPANLESS_VERDICTS = ("none", "no_features", "dispersed")

# Dropped-region annotations that read as self-association. The extraction is
# meant to remove effector function, not the module through which a factor
# assembles on DNA: Nanog dimerizes through the tryptophan repeat UniProt
# annotates at 198-247, which lies outside its 72-aa kept span, and p53
# tetramerizes through 325-356, likewise dropped. Activation and interaction
# regions are NOT on this list -- removing those is the point of the
# extraction, and flagging them would flag most of the library.
SELF_ASSOCIATION_HINTS = ("dimeriz", "oligomeriz", "trimeriz", "tetrameriz",
                          "self-associ", "self associ", "coiled coil")
# A region UniProt describes as a tandem repeat counts only once it is long
# enough to be a module: 'repeat' also names compositional annotations such as
# MESP1's 4-residue "2 X 2 AA tandem repeats of G-Q", which is not a structure
# anyone would reconstitute. Nanog's tryptophan repeat is 45-50 aa.
REPEAT_REGION_HINT = "repeat"
MIN_REPEAT_REGION_AA = 20

# Repeat screen for dropped stretches UniProt does not name, so a periodic
# module is surfaced without a manual read. A ladder is MIN_REPEAT_RUNGS
# occurrences of one residue at constant spacing, tolerating one skipped rung
# (gap 2k) because human Nanog's tryptophan ladder has gaps 5,5,10,5,5,5,5.
# MAX_REPEAT_RESIDUE_FRACTION is what keeps the screen from being swamped by
# composition: without it, poly-Q and poly-G stretches at spacing 2 dominate.
# Measured against live UniProt feature tables for all 2,818 extractable DBD
# records: this screen accounts for 3 of the 97 dbd_dropped_review_flag hits --
# it is the net for stretches no annotation covers, not the main detector. The
# Nanog family is caught by its annotated repeat region instead.
MIN_REPEAT_RUNGS = 6
REPEAT_SPACINGS = range(3, 9)
MAX_REPEAT_RESIDUE_FRACTION = 0.08

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
    "myb-like", "myb-type", "sant", "mads-box",
    "grh/cp2", "cold-shock", "at-hook", "ap2/erf",
    "nuclear receptor", "nr c4-type", "cenp-b",
)
# F5 recommendation 1 -- "ipt/tig" was on the list above and is not.
#
# IPT/TIG is a DNA-binding module in the Rel family, where it is the RHD's
# C-terminal half, and a protein-interaction module elsewhere. Both arms of the
# external validation convict it independently in the two families where it is
# testable: EBF1-4 and RBPJ/RBPJL rest their whole span on an IPT/TIG match,
# Pfam puts their DNA-binding domain somewhere else entirely (COE1 DNA-binding
# domain; Beta-trefoil/LAG1), and the two solved complexes that exist -- 3MLN
# "DNA binding domain of Early B-cell Factor 1 bound to DNA" and 3IAG "CSL
# bound to HES-1" -- put ZERO of their DNA contacts inside the span.
#
# Removing it changes 18 records, in three different ways:
#   * EBF1-4, RBPJ, RBPJL (12, both species) lose their only keep-feature. The
#     Pfam fallback below then gives them the right span, so they end up
#     EXTRACTED rather than abstained -- EBF1 12-253, RBPJ 42-334.
#   * CAMTA1/2 (4) keep their CG-1 DNA_BIND and shed a span that ran to the
#     IPT/TIG 700 aa downstream: CAMTA1 goes from 57-959 (903 aa, flagged) to
#     57-194, and from `bipartite` to `dna_bind` -- the CG-1 is the only
#     DNA-binding feature left. Pfam names no DNA-binding family in CAMTA at
#     all, so this is the annotated DBD standing on its own.
#   * NFATC4 (2) shrinks from 395-689 to 395-585, dropping the IPT/TIG at
#     586-683 (585, not 588, because the flank is pulled back to the excluded
#     domain's edge). This is the cost, and it is real: arm 2 found NFATC2 contacting
#     DNA at 663-666, in the region its paralog annotates as IPT/TIG. It also
#     makes NFATC4 consistent with NFATC1/2/3, whose UniProt entries annotate
#     only RHD and whose spans never reached that region. The NFAT family needs
#     a keep rule of its own, not an accident that covered 1 of 4 paralogs;
#     recorded as a limitation in docs/DBD_EXTERNAL_VALIDATION.md.

# Short domain names that must match the WHOLE description, not a substring of
# it: 'tea' is inside 'Protease' and 'rel' inside 'Related', so substring
# matching here would keep domains that have nothing to do with DNA. The note is
# compared with its trailing repeat index stripped ('Myb 2' -> 'myb').
KEEP_DOMAIN_EXACT = frozenset({
    "tea", "ets", "rhd", "mh1", "mbd", "cut", "pax", "sand", "csd",
    "cbf", "runt", "stat", "tbp", "gcm", "arid", "psq", "hth psq-type",
})

# F5 recommendation 2, part 1. UniProt files some DNA-binding domains as
# REGION rather than DNA_BIND or DOMAIN, and this script read three of the four
# feature types -- so the DBD was in the table and the extraction abstained.
# HSF1's is `REGION 15..120; /note="DNA-binding domain"` with experimental
# evidence (ECO:0000269|PubMed:7935471), and TP73's is `REGION 131..310;
# /note="DNA-binding"` -- which is why TP53 extracts cleanly from a DNA_BIND
# line and its paralog next door did not.
#
# Matched against the WHOLE cleaned note, never as a substring. UniProt also
# writes "No DNA binding activity or transactivation activity..." and p53's
# "Basic (repression of DNA-binding)", both of which contain the substring and
# assert the opposite.
#
# MIN_PLAUSIBLE_SPAN_AA is applied as well, because the same note text does two
# different jobs: a single substantial region naming the DBD (HSF1 106 aa,
# TP73 180 aa, PURB 227 aa), and a per-residue contact footnote. RBPJ carries
# three of the latter -- `REGION 57..67`, `165..170`, `192..197` -- and keeping
# those as boundaries produced a 401-aa span covering 80% of the protein.
# "Interaction with DNA" (76 occurrences) is excluded for the same reason it is
# excluded from the hazard list: those are contacts inside a DBD already found.
#
# Only the two hyphenated forms are kept. Measured over all 2,818 cached
# feature tables: "dna-binding" appears as a whole note 45 times (17 at >= 40
# aa) and "dna-binding domain" twice (both >= 40 aa), while the unhyphenated
# "dna binding"/"dna binding domain" never appear as a whole note at all --
# unhyphenated spelling does occur (10 times) but only inside a longer
# sentence ("Required for DNA binding"), which an exact match refuses. Every
# "dna-binding; minor groove" (4) and "dna-binding; major groove" (2) is a
# per-residue footnote under the 40-aa floor. Carrying them would have implied
# a coverage this rule does not have.
KEEP_REGION_EXACT = frozenset({
    "dna-binding", "dna-binding domain",
})

# Not every zinc finger binds DNA. PHD, RING, ZZ, SP-RING, RanBP2 and MYND
# fingers are protein- or histone-interaction modules, and treating them as part
# of the array is not a near miss -- NSD2's real DBD is the 69-aa HMG box at
# 453-521, but its four PHD fingers (667-1286) dragged the span to 447-1292,
# 846 aa carrying the whole SET methyltransferase. That is precisely the
# self-antagonism this script exists to remove, reintroduced by the extractor.
#
# So ZN_FING features are WHITELISTED by their /note, not taken wholesale.
# Everything not named here is recorded as an excluded feature instead.
#
# F5 recommendation 2, part 2: "matrin" is on the list. A matrin-type finger is
# a C2H2-class finger -- Pfam calls ZNF385A's and Zmat3's "Zinc-finger of C2H2
# type" and CIZ1's "CIZ1-like, C2H2-type zinc finger" -- and omitting the name
# abstained on 8 records whose only annotated features were matrin fingers
# (CIZ1, ZNF318, ZNF385A/B, Zfp318, Zfp385a, Zmat1, Zmat3). Adding it costs no
# existing span: no record loses or changes one.
ZNF_KEEP_TYPES = ("c2h2", "nr c4", "gata", "thap", "bed", "matrin")

# Domains that must never ride along on a fusion. KRAB/SCAN/BTB recruit their own
# corepressors; carrying one into an opener is the failure this script exists for.
#
# NOTE ON WHAT EXCLUSION CAN AND CANNOT DO: the span is min..max over the KEPT
# features, so naming a domain here never shortens it. Exclusion drives the
# `dbd_excluded_domains` record and the review flag -- it is bookkeeping, and
# the real protection is that an excluded domain is not a keep-feature and so
# cannot anchor a boundary. KRAB is absent from the spans because it sits
# N-terminal to the ZnF array, not because this tuple removed it.
EXCLUDE_DOMAIN_HINTS = ("krab", "scan", "btb", "poz", "bromo", "chromo", "phd",
                        "jmjc", "swirm", "elm2", "pwwp", "myb-binding")

# Short excluded names, matched against the WHOLE description like the keep side:
# 'set' is inside 'Basic motif'-free text but 'ank' is inside 'flanking', and an
# exclusion that fires by accident silently discards a real DNA-binding domain.
EXCLUDE_DOMAIN_EXACT = frozenset({
    "set", "pre-set", "post-set", "sam", "ank", "wd", "ring", "zz", "ddt",
    "bah", "lim", "sh2", "pnt", "orange", "nr lbd", "pas", "pac", "mh2",
})

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
    (RAW_DIR / "uniprot").mkdir(parents=True, exist_ok=True)
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
            # Do NOT continue. main() cannot tell "fetch failed" from "entry has
            # no features", so skipping a batch writes dbd_extraction='none' for
            # up to 100 accessions -- erasing a good prior extraction on one
            # transient HTTP error. Fail loudly and let the run be repeated;
            # everything already fetched is cached, so a retry is cheap.
            raise RuntimeError(
                f"UniProt feature fetch failed for batch {i // BATCH + 1} "
                f"({len(chunk)} accessions): {exc}. Nothing was written; re-run "
                f"to resume from the cache."
            ) from exc

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


def _family(row) -> str:
    """`row["subtype"] or "?"` never falls back: a pandas NULL is float('nan'),
    which is truthy, so the family tally logged a literal 'nan' key."""
    v = row.get("subtype")
    return str(v) if v is not None and str(v).lower() != "nan" and str(v) else "?"


def _merge(spans: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping/adjacent (start, end) intervals, so coverage counts each
    residue once -- ZnF arrays and DNA_BIND features routinely overlap."""
    out: list[tuple[int, int]] = []
    for s, e in sorted(spans):
        if out and s <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def _clean_label(note: str) -> str:
    """One region label, safe to join into a ';'-delimited 'start-end:label' field."""
    lab = re.sub(r"\s+", " ", (note or "").replace(";", ",").replace(":", ",")).strip()
    return lab or "unnamed"


def _repeat_motif(segment: str) -> str | None:
    """'W-x5-n10' when `segment` carries a periodic single-residue ladder, else None."""
    n = len(segment)
    if n < MIN_REPEAT_RUNGS:
        return None
    cap = n * MAX_REPEAT_RESIDUE_FRACTION
    for aa in sorted(set(segment)):
        pos = [i for i, c in enumerate(segment) if c == aa]
        if len(pos) < MIN_REPEAT_RUNGS or len(pos) > cap:
            continue
        for k in REPEAT_SPACINGS:
            best = run = 1
            for p, q in zip(pos, pos[1:]):
                if q - p in (k, 2 * k):
                    run += 1
                    best = max(best, run)
                else:
                    run = 1
            if best >= MIN_REPEAT_RUNGS:
                return f"{aa}-x{k}-n{best}"
    return None


def _labelled_gaps(seq: str, lo: int, hi: int,
                   feats: dict) -> list[tuple[int, int, str]]:
    """Every residue OUTSIDE the kept span [lo, hi], split into labelled runs.

    Annotated features supply the label wherever they cover a dropped residue --
    REGION included, which `dbd_excluded_domains` never records, and which is
    how Nanog's annotated tryptophan repeat came to be dropped silently. Runs no
    feature covers are labelled 'unannotated'. Where UniProt gives a residue two
    annotations (Nanog's 198-247 is both the Trp repeat and a transactivation
    region) both are kept, joined with ' + ', so neither label can hide behind
    the other."""
    named: list[tuple[int, int, str]] = []
    for key in ("ZN_FING", "DNA_BIND", "DOMAIN", "REGION"):
        for s, e, note in feats.get(key, []):
            named.append((s, e, _clean_label(note)))

    out: list[tuple[int, int, str]] = []
    for a, b in ((1, lo - 1), (hi + 1, len(seq))):
        if b < a:
            continue
        marks: list[frozenset[str]] = [frozenset()] * (b - a + 1)
        for s, e, lab in named:
            for i in range(max(s, a), min(e, b) + 1):
                marks[i - a] = marks[i - a] | {lab}
        start = a
        for i in range(1, len(marks) + 1):
            if i == len(marks) or marks[i] != marks[i - 1]:
                out.append((start, a + i - 1,
                            " + ".join(sorted(marks[i - 1])) or "unannotated"))
                start = a + i
    return out


def _dropped(seq: str, lo: int, hi: int, feats: dict) -> tuple[list[str], list[str]]:
    """(every dropped region as 'start-end:label', those that need review).

    A named region needs review when its annotation reads as self-association
    (SELF_ASSOCIATION_HINTS), or as a tandem repeat at least
    MIN_REPEAT_REGION_AA long; an unnamed one when the repeat screen finds a
    ladder in it."""
    regions, flagged = [], []
    for s, e, lab in _labelled_gaps(seq, lo, hi, feats):
        if lab == "unannotated":
            motif = _repeat_motif(seq[s - 1:e])
            if motif:
                lab = f"unannotated+repeat({motif})"
        entry = f"{s}-{e}:{lab}"
        low = lab.lower()
        if ("+repeat(" in lab
                or any(h in low for h in SELF_ASSOCIATION_HINTS)
                or (REPEAT_REGION_HINT in low
                    and e - s + 1 >= MIN_REPEAT_REGION_AA)):
            flagged.append(entry)
        regions.append(entry)
    return regions, flagged


def classify_domain(note: str) -> str:
    """'keep' | 'exclude' | 'other' for a DOMAIN description.

    Exclusion wins over keeping: a description that reads as both should not be
    carried into a fusion on the strength of the weaker signal."""
    low = note.lower().strip()
    bare = _TRAILING_INDEX.sub("", low)
    if any(h in low for h in EXCLUDE_DOMAIN_HINTS) or bare in EXCLUDE_DOMAIN_EXACT:
        return "exclude"
    if any(h in low for h in KEEP_DOMAIN_HINTS) or bare in KEEP_DOMAIN_EXACT:
        return "keep"
    return "other"


def extract_one(seq: str, feats: dict, acc: str | None = None) -> dict:
    """
    Compute the DNA-binding span for one record. Returns a dict of the DBD
    columns; `dbd_extraction='none'` when no keep-feature was found.

    `acc` is only used for the Pfam fallback: with no accession, or with no
    cached Pfam entry for it, a record with no UniProt keep-feature abstains
    exactly as before.
    """
    zn, dropped = [], []
    for s, e, note in feats.get("ZN_FING", []):
        if any(t in note.lower() for t in ZNF_KEEP_TYPES):
            zn.append((s, e))
        else:
            dropped.append((f"ZnF:{note or '?'}", s, e))
    dna = [(s, e) for s, e, _ in feats.get("DNA_BIND", [])]
    keep_dom = []
    for s, e, note in feats.get("DOMAIN", []):
        if classify_domain(note) == "keep":
            keep_dom.append((s, e))
        else:
            dropped.append((note or "?", s, e))

    spans = zn + dna + keep_dom

    # F5 recommendation 2, part 1: a REGION whose note NAMES it as the
    # DNA-binding module counts as a keep-feature (see KEEP_REGION_EXACT).
    # Treated as DNA_BIND, because it is the same assertion UniProt makes on a
    # DNA_BIND line, filed under a feature type this script did not read.
    #
    # ONLY on abstention, for the same reason the Pfam fallback below is. Run
    # unconditionally it edits spans that were never abstentions -- 17 records,
    # against the 9 it recovers -- and because a REGION and a DOMAIN routinely
    # annotate the SAME module, `dna and keep_dom` then relabels one module as
    # `bipartite`: ARNT's REGION 88..128 "DNA-binding" sits inside its bHLH
    # 89..142, TFCP2's inside its CP2 domain, PRDM9's inside its own C2H2
    # array. `bipartite` drives the pitch check, the review-queue rank and the
    # per-kind validation tables, so a false one is not cosmetic.
    region_promoted = False
    if not spans:
        promoted = [(s, e) for s, e, note in feats.get("REGION", [])
                    if " ".join((note or "").split()).lower().rstrip(".")
                    in KEEP_REGION_EXACT
                    and (e - s + 1) >= MIN_PLAUSIBLE_SPAN_AA]
        if promoted:
            dna = promoted
            spans = promoted
            region_promoted = True

    pfam_families: list[str] = []
    if not spans:
        # F5 recommendation 2, part 3: the UniProt feature table names nothing
        # usable. Ask Pfam before abstaining. This recovers 30 records -- the
        # whole STAT family (14), EBF/COE (8), GCFC2/PAXBP1 (4) and RBPJ/CSL (4)
        # -- and for STAT1 there is nothing else to go on: its entire feature
        # table is one SH2 domain at 573-670.
        #
        # ONLY on abstention. Pfam never overrides a UniProt keep-feature, so
        # every span that exists today is unaffected, and `dbd_extraction` is
        # set to 'pfam' so no consumer can mistake the provenance. The Pfam
        # family names are written to data/processed/dbd_extraction_report.tsv;
        # they are not a library column.
        hit = dbd_domains(acc) if acc else None
        if hit is None:
            # No kept span, so the extraction dropped nothing: fuse.py falls
            # back to the full-length record, which is F1's problem, not F7's.
            return {"dbd_extraction": "none", "dbd_start": None, "dbd_end": None,
                    "dbd_sequence_aa": None, "dbd_length_aa": None,
                    "dbd_excluded_domains": None, "dbd_review_flag": 0,
                    "dbd_dropped_regions": None, "dbd_dropped_review_flag": 0,
                    "dbd_zinc_finger_count": None, "n_znf": len(zn)}
        spans, pfam_families = hit

    first, last = min(s for s, _ in spans), max(e for _, e in spans)

    # Rule 3: the feature coordinates have to fit the stored sequence. They
    # come from a different place than the sequence does -- the feature table
    # and, since F5, InterPro, which locates domains on UniProt's canonical
    # isoform -- while the span is cut from the stored `sequence_aa`. An
    # isoform mismatch would otherwise be published as a silently truncated
    # span. Measured over all 2,818 extractable records across all three span
    # sources (feature table, promoted REGION, Pfam): 0 trip this, so it is
    # prophylactic rather than a fix for live data.
    #
    # It replaces an `if not sub:` test that could never fire. `lo` clamped to
    # min(first - FLANK, len(seq)) and `hi` to len(seq), so a span beginning
    # past the end of the sequence produced a ONE-residue `sub`, not an empty
    # one, and was emitted as a 1-aa extraction instead of an abstention.
    if last > len(seq):
        return {"dbd_extraction": "none", "dbd_start": None, "dbd_end": None,
                "dbd_sequence_aa": None, "dbd_length_aa": None,
                "dbd_excluded_domains": None, "dbd_review_flag": 1,
                "dbd_dropped_regions": None, "dbd_dropped_review_flag": 0,
                "dbd_zinc_finger_count": None, "n_znf": len(zn),
                "_pfam_families": pfam_families,
                "_region_promoted": region_promoted,
                "_inside": [f"feature coordinates {first}-{last} run past the "
                            f"{len(seq)}-aa stored sequence"]}

    lo = max(1, min(first - FLANK, len(seq)))
    hi = min(len(seq), last + FLANK)

    # The flank must not reach into an excluded domain. ZKSCAN8's first C2H2
    # finger starts at 322 and its KRAB ends at 316, so a 6-residue flank put
    # one KRAB residue in the construct -- 26 records did this. Pull the
    # boundary back to the excluded domain's edge; the span stays contiguous
    # and only the flank is shortened. A domain that is genuinely INSIDE the
    # kept region is untouched here and still raises the review flag below.
    for _n, s, e in dropped:
        if e < first and e >= lo:
            lo = e + 1
        if s > last and s <= hi:
            hi = s - 1
    # lo <= first <= last <= hi holds by construction: the pull-back above only
    # moves `lo` to e + 1 for an excluded domain ending BEFORE `first`, and
    # `hi` to s - 1 for one starting AFTER `last`. So `sub` is never empty and
    # needs no guard of its own -- the guard that matters is the coordinate
    # check above.
    sub = seq[lo - 1:hi]

    # Measuring the gates
    # -------------------
    # `ext` is the annotated feature extent: first kept feature to last, with
    # no flank. Every gate below measures it rather than len(sub), so FLANK
    # cannot move a verdict. The pull-back above only ever shortens the flank
    # (it fires when an excluded domain lies OUTSIDE [first, last]), so the
    # stored span always contains `ext` and `ext <= len(sub)` holds.
    #
    # The two gates then divide the work cleanly, which is what the first
    # attempt at FLANK=15 got wrong by moving one and not the other:
    #   * coverage = covered / ext asks "are these features one module, or two
    #     distant ones with a protein between them" -- the MGA case. For a
    #     single merged interval it is 1.000, which is correct and not a
    #     silence, because there is no gap to measure.
    #   * ext < MIN_PLAUSIBLE_SPAN_AA asks "is there enough annotated feature
    #     here to be a domain at all" -- the KMT2C case, a 13-aa A.T hook.
    #     That question is about the feature, never about the padding.
    ext = last - first + 1

    # Rule 5 (F2): the span is the protein. Report the ratio and write no span,
    # rather than storing a near-full-length sequence under a DBD column name.
    # dbd_review_flag is set as well, so a `dispersed` record cannot be mistaken
    # for a clean abstention by a consumer that only looks at the flag.
    frac = ext / len(seq)
    if (frac >= MAX_SPAN_FRACTION and ext >= MIN_DISPERSED_SPAN_AA) \
            or ext >= MAX_SPAN_AA:
        return {"dbd_extraction": "dispersed", "dbd_start": None, "dbd_end": None,
                "dbd_sequence_aa": None, "dbd_length_aa": None,
                "dbd_excluded_domains": ";".join(n for n, _, _ in dropped) or None,
                "dbd_review_flag": 1,
                "dbd_dropped_regions": None, "dbd_dropped_review_flag": 0,
                "dbd_zinc_finger_count": None, "n_znf": len(zn),
                "_inside": [f"feature extent {ext} aa = {frac:.0%} of the "
                            f"{len(seq)}-aa protein ({first}-{last}); dispersed"]}

    # Rule 2: a dropped feature the span carries is a conflict. The test is
    # OVERLAP, not containment: a domain straddling a boundary is still partly
    # in the construct, and containment let TADA2A carry SWIRM unflagged.
    def _overlaps(s, e):
        return s <= hi and e >= lo

    inside = [n for n, s, e in dropped if _overlaps(s, e)]
    for s, e, note in feats.get("REGION", []):
        low = (note or "").lower()
        if any(h in low for h in HAZARD_REGION_HINTS) and _overlaps(s, e):
            inside.append(note)

    # A span that is mostly NOT DNA-binding feature is a boundary failure even
    # when nothing named lands inside it: MGA extracts 78-2480 (2,403 aa of a
    # 3,065-aa protein) because two distant keep-features bracket the whole
    # middle. Gate 1 only catches spans that are too SHORT.
    #
    # Denominator is `ext`, so this is a fragmentation measure: how much of
    # the region the features bracket is actually feature. It reads 1.000 for
    # a single merged interval by construction (1,361 of the span-bearing
    # records), which is why the floor below moved to `ext` in the same
    # change -- on len(sub) these two gates overlapped and both drifted with
    # FLANK; on `ext` they ask different questions and neither does.
    covered = sum(e - s + 1 for s, e in _merge(spans))
    coverage = covered / ext
    if coverage < MIN_SPAN_COVERAGE:
        inside.append(f"bracketed region {coverage:.0%} DNA-binding feature")

    # F2b: too short to hold one complete DNA-binding module, whatever the kind.
    # The acceptance gate checks ZnF pitch only, and only for znf_array and
    # bipartite, so dna_bind spans of 25 aa passed with the flag at 0.
    if ext < MIN_PLAUSIBLE_SPAN_AA:
        inside.append(f"feature extent {ext} aa, under the "
                      f"{MIN_PLAUSIBLE_SPAN_AA}-aa minimum for a DNA-binding module")

    if pfam_families:
        kind = "pfam"
    elif zn and (dna or keep_dom):
        kind = "bipartite"
    elif zn:
        kind = "znf_array"
    elif dna and keep_dom:
        kind = "bipartite"
    elif dna:
        kind = "dna_bind"
    else:
        kind = "domain"

    # F7: what fell outside the span, named or not. `dbd_excluded_domains` above
    # records only annotated ZN_FING and DOMAIN features, so an annotated REGION
    # -- Nanog's tryptophan repeat -- was dropped without appearing anywhere.
    regions, dropped_flagged = _dropped(seq, lo, hi, feats)

    return {
        "dbd_extraction": kind,
        "dbd_start": lo,
        "dbd_end": hi,
        "dbd_sequence_aa": sub,
        "dbd_length_aa": len(sub),
        "dbd_excluded_domains": ";".join(n for n, _, _ in dropped) or None,
        "dbd_review_flag": 1 if inside else 0,
        "dbd_dropped_regions": ";".join(regions) or None,
        "dbd_dropped_review_flag": 1 if dropped_flagged else 0,
        # F3: fingers the span actually KEPT. `zinc_finger_count` counts every
        # ZN_FING feature over the whole protein, PHD/RING/FOG-type included, so
        # af3.py added Zn ions for fingers the construct does not contain --
        # ZFPM2 claims 8 for a 102-aa span of a 1,151-aa protein.
        #
        # None, not 0, when the span did not come from the UniProt feature
        # table: a Pfam-derived span is positioned by a source that says
        # nothing about ZN_FING, and PFAM_DBD_HINTS does name C2H2 families, so
        # a 0 here would assert "no fingers in this span" on no evidence and
        # suppress af3.py's warn-and-fall-back to the protein-wide count.
        "dbd_zinc_finger_count": None if pfam_families else len(zn),
        "n_znf": len(zn),
        "_inside": inside,
        "_dropped_flagged": dropped_flagged,
        "_pfam_families": pfam_families,
        "_region_promoted": region_promoted,
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
        fam = _family(r)
        families_total[fam] += 1
        feats = parse_features(features.get(acc, ""))
        for _, _, note in feats.get("DOMAIN", []):
            v = classify_domain(note)
            by_verdict[v] += 1
            if v == "other":
                unknown[note] += 1
        # Must name every keep path extract_one has, or --survey reports
        # abstentions the run does not make. The Pfam fallback is not checked
        # here: it needs the interpro cache, and survey() is a feature-table
        # report, so STAT/EBF/RBPJ/GCFC2 still appear below as "no keep
        # feature" -- true of the UniProt table, which is what this surveys.
        if not (feats.get("ZN_FING") or feats.get("DNA_BIND")
                or any(classify_domain(n) == "keep" for _, _, n in feats.get("DOMAIN", []))
                or any(" ".join((n or "").split()).lower().rstrip(".") in KEEP_REGION_EXACT
                       and (e - st + 1) >= MIN_PLAUSIBLE_SPAN_AA
                       for st, e, n in feats.get("REGION", []))):
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
# to be re-justified, not a silent drift. The coordinates are FLANK-dependent,
# which makes this a drift detector rather than independent evidence -- arm 1
# and arm 2 of docs/DBD_EXTERNAL_VALIDATION.md supply that. The point here is
# that these nine spans cannot change without someone noticing. The bracketed
# values are what they were at FLANK=6; the boundaries move with the flank,
# the VERDICTS (flag, dispersed) no longer do.
ANCHORS = {
    "HNF1A":  (72, 294),      # was (81, 285) at FLANK=6
    "EGR3":   (260, 370),     # was (269, 361)
    "ZNF793": (212, 404),     # was (221, 395)
    # F2 reclassified this one deliberately, so the anchor is the VERDICT, not a
    # span: ZNF850's extraction ran 191-1090, i.e. 900 aa = 83% of a 1,090-aa
    # protein, ending at the C-terminus. That is what MAX_SPAN_FRACTION is for.
    # It is also the fib->beta opener whose native KRAB motivated this script,
    # and a 900-aa chassis was never usable anyway -- ~980 aa of fusion, ~2.9 kb
    # against a ~4.7 kb AAV cargo (DBD_EXTRACTION_PLAN.md "Still open").
    "ZNF850": "dispersed",
    "ZNF540": (172, 660),     # was (181, 660); the end is the C-terminus
    "ZNF444": (164, 315),     # was (173, 306)
    "HINFP":  (1, 383),       # was (9, 374); the start clamps to residue 1
    # The two rules F5 added, each with the record that motivated it, so neither
    # can go silently inert -- which is the real risk for the Pfam fallback,
    # since an empty data/raw/interpro/ turns it off without any error.
    # The NFAT family is RHD-only by decision, not by accident: NFATC4 is the
    # one paralog that annotates the IPT/TIG its siblings leave unannotated, so
    # a rule reaching that region would fire on one of four. See
    # docs/DBD_EXTERNAL_VALIDATION.md Limitations. 585, not 597, because the
    # pull-back keeps the flank out of the IPT/TIG starting at 586.
    "NFATC4": (386, 585),
    "HSF1":   (1, 135),       # REGION 15..120, flanked and clamped to 1
    "EBF1":   (3, 261),       # Pfam fallback: COE1 DNA-binding domain
}

KRAB_HINT = "krab"


def validate(lib, log) -> int:
    """Three acceptance gates. Returns the number of failures."""
    df = lib.to_dataframe("DBD")
    done = df[df["dbd_extraction"].notna()
               & ~df["dbd_extraction"].isin(SPANLESS_VERDICTS)]
    fails = 0

    # Gate 1 -- length sanity against the independently stored finger count.
    # A C2H2 array runs about 23 aa per finger plus ~5 aa TGEKP linkers, so
    # 28 aa/finger is the expected pitch; allow generous slack for degenerate
    # fingers and split arrays, and only flag spans that are too SHORT to hold
    # the fingers claimed (too long is legitimate -- split arrays span a gap).
    #
    # The reference is the number of fingers the span actually KEPT, re-derived
    # from the cached feature table rather than read back from the extraction --
    # `zinc_finger_count` counts every ZN_FING feature, PHD and RING included,
    # so it is an upper bound and comparing against it fails records whose
    # non-DNA-binding fingers were (correctly) left out.
    znf = done[done["dbd_extraction"].isin(["znf_array", "bipartite"])
               & done["uniprot_id"].notna()]
    short, overcount = [], []
    for _, r in znf.iterrows():
        p = _cache_path(str(r["uniprot_id"]).strip())
        if not p.exists():
            continue
        notes = parse_features(p.read_text(encoding="utf-8")).get("ZN_FING", [])
        kept = sum(1 for _, _, n in notes
                   if any(t in n.lower() for t in ZNF_KEEP_TYPES))
        if kept and r["dbd_length_aa"] < 20 * kept:
            short.append((r["gene_symbol"], r["dbd_length_aa"], kept))
        if r["zinc_finger_count"] and r["zinc_finger_count"] > len(notes):
            overcount.append(r["gene_symbol"])
    log.info(f"[gate 1/length] {len(znf):,} ZnF records checked, "
             f"{len(short):,} shorter than 20 aa per KEPT finger")

    # F2b: the pitch test above covers znf_array and bipartite only, so a 25-aa
    # dna_bind span was never length-checked at all. This one applies to every
    # kind and needs no feature cache.
    tiny = done[done["dbd_length_aa"].notna()
                & (done["dbd_length_aa"] < MIN_PLAUSIBLE_SPAN_AA)]
    log.info(f"[gate 1/length] {len(tiny):,} spans under "
             f"{MIN_PLAUSIBLE_SPAN_AA} aa across all kinds "
             f"({dict(tiny['dbd_extraction'].value_counts())})")
    for _, r in tiny.head(10).iterrows():
        log.info(f"    {r['gene_symbol']}: {int(r['dbd_length_aa'])} aa "
                 f"({r['dbd_extraction']}), flag {r['dbd_review_flag']}")
    # A short span is not itself a failure -- KMT2C's 25-aa homeobox is what
    # UniProt annotates, and the honest response is to flag it, not to invent a
    # boundary. The failure is a short span that is NOT flagged, because that
    # one reaches a design run as if it were a clean chassis.
    unflagged_tiny = tiny[tiny["dbd_review_flag"].fillna(0).astype(float) != 1]
    log.info(f"[gate 1/length] {len(unflagged_tiny):,} of those carry no "
             f"review flag")
    if len(unflagged_tiny):
        fails += 1
        for _, r in unflagged_tiny.head(10).iterrows():
            log.warning(f"    UNFLAGGED {r['gene_symbol']}: "
                        f"{int(r['dbd_length_aa'])} aa ({r['dbd_extraction']})")
    if short:
        fails += 1
        for gene, length, kept in short[:10]:
            log.warning(f"    {gene}: {length} aa for {kept} kept fingers")
    if overcount:
        log.warning(f"[gate 1/length] {len(overcount)} records where "
                    f"zinc_finger_count exceeds the annotated ZN_FING features")

    # Gate 2 -- known anchors.
    bad = []
    for gene, expected in ANCHORS.items():
        hit = df[(df["gene_symbol"] == gene) & (df["organism"].astype(str).str.contains("sapiens", case=False, na=False))]
        if hit.empty:
            hit = df[df["gene_symbol"] == gene]
        if hit.empty:
            bad.append(f"{gene}: not in library"); continue
        r = hit.iloc[0]
        if isinstance(expected, str):
            if r["dbd_extraction"] != expected:
                bad.append(f"{gene}: got {r['dbd_extraction']}, "
                           f"expected {expected}")
            continue
        s, e = expected
        if (r["dbd_start"], r["dbd_end"]) != (s, e):
            bad.append(f"{gene}: got {r['dbd_start']}-{r['dbd_end']}, expected {s}-{e}")
    log.info(f"[gate 2/anchors] {len(ANCHORS) - len(bad)}/{len(ANCHORS)} match")
    for b in bad:
        log.warning(f"    {b}")
    if bad:
        fails += 1

    # Gate 3 -- KRAB audit. This is the specific failure that motivated the work,
    # so it gets its own gate rather than being assumed away by gate 2.
    #
    # The KRAB span is re-read from the feature table and intersected with the
    # kept span, rather than trusting dbd_review_flag -- the flag is set by any
    # of several conditions, so testing it here would pass or fail for reasons
    # that have nothing to do with KRAB.
    #
    # Two outcomes are distinguished. A KRAB the span reaches into on its own is
    # a boundary BUG and fails the gate. A KRAB whose UniProt annotation
    # OVERLAPS the first kept finger is an annotation conflict with no correct
    # trim: pulling the boundary back would cut into the finger. Those are
    # reported, not failed, because the honest response is to exclude the
    # chassis at design time, which dbd_review_flag already asks for.
    krab = done[done["dbd_excluded_domains"].fillna("").str.lower().str.contains(KRAB_HINT)
                & done["uniprot_id"].notna()]
    bugs, conflicts = [], []
    for _, r in krab.iterrows():
        p = _cache_path(str(r["uniprot_id"]).strip())
        if not p.exists():
            continue
        feats = parse_features(p.read_text(encoding="utf-8"))
        keeps = [(s, e) for s, e, n in feats.get("ZN_FING", [])
                 if any(t in n.lower() for t in ZNF_KEEP_TYPES)]
        keeps += [(s, e) for s, e, _ in feats.get("DNA_BIND", [])]
        first_keep = min((s for s, _ in keeps), default=None)
        for s, e, n in feats.get("DOMAIN", []):
            if KRAB_HINT not in n.lower():
                continue
            if s <= r["dbd_end"] and e >= r["dbd_start"]:
                entry = (r["gene_symbol"], f"{int(r['dbd_start'])}-{int(r['dbd_end'])}",
                         f"KRAB {s}-{e}")
                (conflicts if first_keep is not None and e >= first_keep
                 else bugs).append(entry)
    log.info(f"[gate 3/KRAB] {len(krab):,} records dropped a KRAB domain; "
             f"{len(bugs)} span reaches into one, {len(conflicts)} annotation "
             f"overlaps a kept finger")
    if bugs:
        fails += 1
        for gene, span, krab_span in bugs[:10]:
            log.warning(f"    BUG {gene}: span {span} vs {krab_span}")
    for gene, span, krab_span in conflicts[:10]:
        log.info(f"    conflict {gene}: span {span} vs {krab_span} "
                 f"(flagged for review; no correct trim exists)")

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
            # F9: the fetch returned nothing for this accession, which is NOT
            # the same as "the entry has features, none of them usable". The
            # fetch logic above argues exactly this and the old code then wrote
            # both as 'none', so the 262 `none` records were two causes mixed.
            # `no_features` accessions are retryable; `none` ones are not.
            res = {"dbd_extraction": "no_features", "dbd_start": None,
                   "dbd_end": None, "dbd_sequence_aa": None,
                   "dbd_length_aa": None,
                   "dbd_excluded_domains": None, "dbd_review_flag": 0,
                   "dbd_dropped_regions": None, "dbd_dropped_review_flag": 0,
                   "dbd_zinc_finger_count": None,
                   "n_znf": 0, "_inside": []}
        else:
            res = extract_one(seq, parse_features(line), acc)

        kinds[res["dbd_extraction"]] += 1
        by_family[_family(r)][res["dbd_extraction"]] += 1

        updates.append((
            seq, res["dbd_start"], res["dbd_end"], res["dbd_sequence_aa"],
            res["dbd_length_aa"], res["dbd_extraction"],
            res["dbd_excluded_domains"], res["dbd_review_flag"],
            res["dbd_dropped_regions"], res["dbd_dropped_review_flag"],
            res["dbd_zinc_finger_count"],
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
            "dbd_zinc_finger_count": res["dbd_zinc_finger_count"],
            "excluded_domains": res["dbd_excluded_domains"],
            "review_flag": res["dbd_review_flag"],
            "inside_kept_span": ";".join(res.get("_inside", [])) or None,
            "dropped_regions": res["dbd_dropped_regions"],
            "dropped_review_flag": res["dbd_dropped_review_flag"],
            "dropped_needs_review": ";".join(res.get("_dropped_flagged", [])) or None,
            # Provenance for the Pfam fallback: `dbd_extraction='pfam'` says a
            # span came from Pfam, this says which families it came from. Kept
            # here rather than in a library column.
            "pfam_families": "; ".join(res.get("_pfam_families", [])) or None,
            # Same purpose for the other new source: `dbd_extraction` cannot
            # carry it, because a promoted REGION is recorded as `dna_bind` --
            # which is exactly what UniProt asserted, on the wrong feature line.
            "region_promoted": 1 if res.get("_region_promoted") else 0,
        })

    log.info(f"[extract] by kind: {dict(kinds)}")
    log.info("[extract] coverage by family (extracted / total):")
    for fam, c in sorted(by_family.items(), key=lambda kv: -sum(kv[1].values()))[:20]:
        tot = sum(c.values())
        log.info(f"    {tot - c['none']:5d} / {tot:5d}  {fam}")
    flagged = sum(1 for x in report if x["review_flag"])
    log.info(f"[extract] {flagged:,} records need review "
             f"(an excluded domain or effector region lies inside the kept span)")
    log.info(f"[extract] spanless verdicts: "
             f"none {kinds['none']:,} (rule-3 abstention), "
             f"no_features {kinds['no_features']:,} (fetch returned nothing, "
             f"retryable), dispersed {kinds['dispersed']:,} (span would have "
             f"been >= {MAX_SPAN_FRACTION:.0%} of the protein and "
             f">= {MIN_DISPERSED_SPAN_AA} aa, or >= {MAX_SPAN_AA} aa)")
    pfam_rows = [x for x in report if x["dbd_extraction"] == "pfam"]
    if pfam_rows:
        fams = Counter(x["pfam_families"] for x in pfam_rows)
        log.info(f"[extract] {len(pfam_rows):,} spans came from the Pfam "
                 f"fallback (no UniProt keep-feature); by family: "
                 f"{dict(fams.most_common(8))}")
    else:
        log.info("[extract] no Pfam fallback spans. If that is unexpected, the "
                 "cache is missing: run scripts/07_validate_dbd_external.py "
                 "--fetch to populate data/raw/interpro/.")
    short = sum(1 for x in report
                if x["inside_kept_span"] and "-aa minimum" in x["inside_kept_span"])
    log.info(f"[extract] {short:,} spans under {MIN_PLAUSIBLE_SPAN_AA} aa flagged")
    dropped_flagged = sum(1 for x in report if x["dropped_review_flag"])
    log.info(f"[extract] {dropped_flagged:,} records dropped a self-association "
             f"or repeat module OUTSIDE the kept span (dbd_dropped_review_flag; "
             f"under-inclusion, a different failure class from the above)")

    import pandas as pd
    out = PROCESSED_DIR / "dbd_extraction_report.tsv"
    pd.DataFrame(report).to_csv(out, sep="\t", index=False)
    log.info(f"[extract] wrote {out}")

    # F5 recommendation 4. Both arms of the external validation rank the kinds
    # the same way, and `domain` is last on both: recall >= 0.95 for 93% of
    # records against 98-99% for the others, and 40 of 51 with full contact
    # coverage against 91/97 for dna_bind. That is the kind whose boundaries
    # rest on KEEP_DOMAIN_HINTS alone, with no DNA_BIND or ZN_FING line under
    # them. `pfam` ranks ahead of it only because arm 1 cannot check it at all
    # -- those spans ARE Pfam, so the comparison is circular, and they have
    # never been reviewed. Flagged records only; the order is the whole point,
    # so this is a separate file rather than a column on the report.
    priority = {"pfam": 1, "domain": 2, "bipartite": 3, "znf_array": 4,
                "dna_bind": 5, "dispersed": 6, "none": 7, "no_features": 8}
    # Flagged records, PLUS every span whose provenance is a source this
    # script did not use before, flagged or not. Ranking `pfam` first because
    # arm 1 cannot check it is empty if 21 of its 30 spans never enter the
    # queue; the same holds for a promoted REGION, whose largest effect
    # (MLLT10, a 364-aa REGION) clears every gate.
    queue = [dict(x, review_priority=priority.get(x["dbd_extraction"], 9))
             for x in report
             if x["review_flag"] or x["dropped_review_flag"]
             or x["dbd_extraction"] == "pfam" or x["region_promoted"]]
    queue.sort(key=lambda x: (x["review_priority"], x["gene_symbol"] or ""))
    qout = PROCESSED_DIR / "dbd_review_queue.tsv"
    pd.DataFrame(queue).to_csv(qout, sep="\t", index=False)
    log.info(f"[extract] wrote {qout}: {len(queue):,} records to review, "
             f"worst kind first "
             f"({dict(Counter(x['dbd_extraction'] for x in queue))})")

    if args.dry_run:
        log.info("[dry-run] library not modified")
        return

    with ModuleLibrary(args.db) as lib:
        with lib.transaction() as cur:
            cur.executemany(
                "UPDATE modules SET sequence_full_aa=?, dbd_start=?, dbd_end=?, "
                "dbd_sequence_aa=?, dbd_length_aa=?, dbd_extraction=?, "
                "dbd_excluded_domains=?, dbd_review_flag=?, "
                "dbd_dropped_regions=?, dbd_dropped_review_flag=?, "
                "dbd_zinc_finger_count=?, "
                "date_modified=? WHERE module_id=?", updates)
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
                      f"{flagged} flagged for review; {dropped_flagged} flagged "
                      f"for a dropped repeat or self-association module"),
        })
    log.info(f"[extract] updated {len(updates):,} records in {args.db}")
    log.info("Next: python scripts/06_extract_dbd.py --validate")


if __name__ == "__main__":
    main()
