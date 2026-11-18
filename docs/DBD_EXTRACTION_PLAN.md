# DBD extraction plan for eCR_mod_lib

**Status:** **implemented 2026-09-04** — `scripts/06_extract_dbd.py`, schema columns
live, all three validation gates pass. Results at the end of this document.
Originally proposed after the fib→β eCR run (2026-09-03), which could not
build correct fusions from the library as it stands. Seven chassis were extracted
by hand for that run; this is the plan to do it for all 3,248.

## The problem

`scripts/01_fetch_dbd.py` writes the **entire UniProt sequence** into
`sequence_aa` for every `type='DBD'` record:

```python
"sequence_aa": seq,
"length_aa": len(seq) if seq else None,
```

So the `DBD` type holds whole transcription factors. `DBD_HO_HNF1A` is all 631 aa
of HNF1A; its homeodomain is 81. Nothing in the schema records where the
DNA-binding part starts or ends.

For motif scanning this is harmless — `target.py` and `offtarget.py` use the
JASPAR PWM, not the protein. It becomes serious at `fuse.py`, which concatenates
`sequence_aa` to an effector domain. That builds a full-length-TF fusion carrying
every native domain the TF has.

### Why this is not cosmetic

Measured on the seven chassis the fib→β run selected:

| chassis | full | DNA-binding | fraction | also carried |
|---|---|---|---|---|
| HNF1A | 631 | 205 | 32% | HNF-p1 dimerization, C-terminal TAD |
| EGR3 | 387 | 93 | 24% | disordered regions |
| **ZNF793** | 406 | 175 | 43% | **KRAB** |
| **ZNF850** | 1090 | 900 | 83% | **KRAB** |
| ZNF540 | 660 | 480 | 73% | KRAB |
| ZNF444 | 327 | 134 | 41% | SCAN box |
| HINFP | 517 | 366 | 71% | NPAT interaction, activation region |

**ZNF793 and ZNF850 were selected as openers.** KRAB recruits KAP1/TRIM28 →
SETDB1 → H3K9me3. Fusing a BAF recruiter to a KRAB-bearing protein builds a
construct that opens and silences the same locus at once — self-antagonism inside
a single chain, which no design-level rule (such as the opener/closer DBD
disjointness rule) can detect or prevent.

**~39% of the library is `zf-C2H2` (1,262 of 3,248).** KRAB-ZnFs are a large part
of that family, so this is a systematic hazard, not seven unlucky picks.

## Scope

| | count |
|---|---|
| `type='DBD'` rows | 3,248 |
| with `uniprot_id` and `sequence_aa` | **2,818** |
| with `jaspar_id` | 1,177 |
| with `zinc_finger_count` | 1,010 |
| human / mouse | 1,637 / 1,611 |

430 rows have no UniProt accession and no sequence; they cannot be extracted and
should be marked, not guessed at.

Families: zf-C2H2 1262, Homeobox 423, bHLH 214, Others 165, bZIP 108, HMG 104,
ZBTB 96, Fork_head 96, MYB 62, ETS 57.

## Method

Source of truth is the **UniProt feature table**: `ft_zn_fing`, `ft_dna_bind`,
`ft_domain`, `ft_region`.

> **Correction (2026-09-04).** The draft said these were already cached, so the
> pass would be fully offline. They were not. `data/raw/uniprot/` holds `.fasta`
> sequences (3,253) and `_znf.tsv` zinc-finger tables (1,008) — `01_fetch_dbd.py`
> fetches `ft_zn_fing` only. `ft_dna_bind` / `ft_domain` / `ft_region` had to be
> fetched. This turned out not to matter: the search endpoint takes 100
> accessions per request, so the whole library was 29 requests and a few seconds,
> and results are now cached per accession as `<acc>_feat.tsv`.

**Keep** the DNA-contacting apparatus:

- `Zinc finger` spans → first start to last end **as one contiguous block**.
  The inter-finger TGEKP linkers are part of the array; splicing fingers together
  individually destroys the register and the DNA contacts.
- `DNA binding` spans (Homeobox, HMG box, Fork head, T-box, ...).
- Partner domains that form part of the interface, matched by description:
  `POU-specific`, `POU-like`, `Homeo`. HNF1A binds as a bipartite POU-specific +
  homeodomain unit; taking the homeobox alone removes half the interface.

**Drop** everything else, explicitly including `KRAB`, `SCAN`, `BTB`/`POZ`, and
interaction/activation `Region`s.

Then take `min(start) - FLANK` to `max(end) + FLANK`, with **FLANK = 6**.

### Rules that matter

1. **Contiguous span, never concatenated fragments.** Splicing across an excluded
   internal domain would create a junction that exists in no real protein — a
   novel neoepitope surface and an unpredictable fold.
2. **A dropped domain lying INSIDE the kept span is a conflict, not a silent
   pass.** Flag it for review; the current implementation reports these under
   `inside_kept_span`. For the seven chassis there were none — every KRAB/SCAN
   sat N-terminal to the ZnF array, which is the canonical arrangement.
3. **Never fabricate boundaries.** No UniProt features → leave `sequence_aa`
   as-is and set `dbd_extraction = 'none'`. A guessed domain is worse than an
   honest full-length record.
4. **Keep the full sequence.** Extraction is additive; nothing is destroyed.

## Schema change

Add to `modules`, all nullable so existing consumers keep working:

| column | meaning |
|---|---|
| `dbd_start`, `dbd_end` | 1-indexed span in `sequence_full_aa` |
| `dbd_sequence_aa` | the extracted region |
| `dbd_length_aa` | its length |
| `sequence_full_aa` | the original full-length sequence, preserved |
| `dbd_extraction` | `znf_array` / `dna_bind` / `bipartite` / `none` |
| `dbd_excluded_domains` | `;`-joined, e.g. `KRAB` |
| `dbd_review_flag` | true when an excluded domain fell inside the kept span |

**Do not overwrite `sequence_aa`.** Several call sites read it; a silent change of
meaning would be the same class of bug as the one being fixed. Migrate consumers
deliberately: `fuse.py` should prefer `dbd_sequence_aa` and fall back with a
warning.

## Steps

1. New `scripts/06_extract_dbd.py` — offline over the cached UniProt JSON, since
   `01_fetch_dbd.py` already retrieves it. No re-fetch, no rate limit.
2. Extend `schema.py` with the columns above; `04_build_library.py` carries them
   through to the TSV export.
3. Run over all 2,818 extractable records; write a coverage report by family.
4. **Validate before trusting it**, in three ways:
   - *Length sanity* — a C2H2 array should be about 23 aa per finger plus linkers.
     `zinc_finger_count` is already stored for 1,010 records; cross-check.
   - *Known anchors* — a hand-checked set (HNF1A 81-285, EGR3 269-361,
     ZNF793 221-395, ZNF850 191-1090, ZNF540 181-660, ZNF444 173-306,
     HINFP 9-374) as a regression fixture.
   - *KRAB audit* — after extraction, assert no `dbd_sequence_aa` contains a
     KRAB span. That is the specific failure that motivated this work, so it
     deserves its own test rather than being assumed.
5. Point `fuse.py` at `dbd_sequence_aa`; warn and fall back where absent.
6. Re-run the fib→β fusions and confirm the by-hand and pipeline extractions agree.

## Known edge cases

- **Very large arrays.** ZNF850 is 32 fingers, 900 aa extracted — barely smaller
  than the 1,090-aa protein. Extraction does not make every chassis compact, and
  a ~980-aa fusion has real delivery consequences (AAV cargo is about 4.7 kb
  total). Flag `dbd_length_aa` in construct design rather than assuming trimming
  solved it.
- **Boundary regions that do double duty.** HINFP 374-407 is annotated *"Required
  for activation ... and contributes to DNA-binding"*. Excluded here, because
  carrying an activation region into a NuRD-recruiting closer is exactly the
  self-antagonism being removed — but the DNA-binding contribution is genuinely
  lost. Such cases need a per-record note, not a silent rule.
- **Degenerate fingers.** UniProt marks some C2H2 fingers `degenerate`; they are
  kept, since they sit inside the array and removing them would break the register.
- **Non-C2H2, non-classical binders** (ZBTB, bHLH dimerization, HMG) need their
  own keep-lists; the description-matching keep list here covers POU/homeo only.

## Cost

**Measured 2026-09-04:** 29 batched UniProt requests (~4 s) plus a single-threaded
local pass; the whole run, extraction and library update included, was **under a
minute**. No GPU. The computational resources that deferred this work were never
the constraint.

---

## As implemented (2026-09-04)

`scripts/06_extract_dbd.py`, with `--survey` (domain vocabulary, no writes),
`--dry-run` (extract + report, no writes), `--validate` (the three gates).

### Coverage

2,818 of 3,248 DBD rows are extractable; **2,556 (91%) got a span**.

| `dbd_extraction` | n | meaning |
|---|---|---|
| `znf_array` | 1,055 | contiguous first-finger → last-finger block |
| `dna_bind` | 793 | a `DNA_BIND` feature |
| `domain` | 463 | a DNA-binding `DOMAIN` (bHLH, bZIP, MBD, …) with no `DNA_BIND` |
| `bipartite` | 245 | binding feature + partner domain (e.g. HNF1A POU + homeobox) |
| `none` | 262 | no usable feature — full-length record left honest, per rule 3 |
| *(NULL)* | 430 | no accession or no sequence; not touched |

`domain` was added to the four kinds the draft listed: bHLH, bZIP, MBD, RHD, MH1
and others annotate their binding module as `DOMAIN` with no `DNA_BIND` feature
at all, so without it 214 bHLH and 108 bZIP records would have fallen to `none`.

By family: zf-C2H2 911/927, Homeobox 378/379, bHLH 208/208, Fork_head 94/94,
ZBTB 83/83, MYB 52/52, T-box 34/34. The families that stay at `none` are
`Others`, MH1, STAT, Tub, AP-2, GTF2I, AF-4 — each needs a curated keep rule,
and guessing one is worse than reporting the gap.

### Keep list is data-driven, and short names match exactly

`--survey` prints every unclassified `DOMAIN` description, which is how the keep
list was built rather than assumed. It also caught a real hazard: substring
matching on short domain names is wrong — `tea` is inside *Protease* and `rel`
inside *Related*. Short names live in `KEEP_DOMAIN_EXACT` / `EXCLUDE_DOMAIN_EXACT`
and must match the whole description (with any trailing repeat index stripped).

The survey confirms the drops are right too: the most common dropped domain is
**NR LBD** (96) — the nuclear-receptor *ligand*-binding domain, correctly not
carried into a fusion — followed by SET, PAS, Orange, ELM2, LIM, SH2.

## Post-review corrections (2026-09-04)

A code review of the first implementation found a defect that reintroduced the
exact failure the script exists to prevent, plus five smaller ones. All are
fixed; the numbers above are the corrected run.

**Not every zinc finger binds DNA.** `ZN_FING` features were kept
unconditionally, so PHD, RING, ZZ, SP-RING and RanBP2 fingers — protein- and
histone-interaction modules — anchored the span. **NSD2**'s real DBD is the
69-aa HMG box at 453-521, but its four PHD fingers at 667-1286 stretched the
span to 447-1292: **846 aa carrying the entire SET methyltransferase**, on a
record the gates passed. `ZNF_KEEP_TYPES` now whitelists C2H2, NR C4, GATA,
THAP and BED fingers by their `/note`; NSD2 extracts to 447-527. TADA2A
(anchored at 6 by a ZZ finger, 436 of 443 aa) and ZMIZ2 (labelled `znf_array`
on an SP-RING SUMO-ligase finger) are likewise fixed.

**The review flag tested containment, not overlap** — a dropped domain
straddling a boundary was carried but never flagged, and Gate 3 passed
vacuously on those. Now an overlap.

**The span had no upper bound.** MGA extracted 78-2480 — 2,403 aa of a 3,065-aa
protein — because two distant keep-features bracketed the middle, and nothing
named landed inside, so `dbd_review_flag` stayed 0. Spans whose DNA-binding
feature coverage falls below `MIN_SPAN_COVERAGE` (0.5) are now flagged: MGA at
10%, TRPS1 19%, ADNP 37%. Review flags went 32 → 287, which is the honest count.

**The flank could reach into an excluded domain.** ZKSCAN8's first finger starts
at 322 and its KRAB ends at 316, so `FLANK = 6` put a KRAB residue in the
construct — 26 records. The boundary is now pulled back to the excluded domain's
edge; only the flank shortens, the span stays contiguous.

**A failed fetch batch was skipped**, and `main()` could not tell "fetch failed"
from "no features" — one transient HTTP error would have written
`dbd_extraction='none'` over a good extraction for up to 100 accessions. It now
raises; everything already fetched is cached, so a re-run resumes.

**`04_build_library.py` silently erased the whole extraction.**
`insert_or_replace` filled every schema column from the record and
INSERT-OR-REPLACEd; no source TSV carries `dbd_*`, so an ordinary pipeline run
NULLed all 3,248 DBD rows — and `fuse.py` would then quietly fall back to
full-length sequences. An update now touches only the columns the record
mentions. Verified: after a full `04` run, 2,556 spans survive.

Two smaller ones: `lo` was not clamped to the sequence length (an isoform
mismatch would emit a zero-length "extraction"), and `row["subtype"] or "?"`
never fell back because a pandas NULL is `float('nan')`, which is truthy, so the
family tally logged a literal `nan` key.

### Validation

All three gates pass (`--validate`, exit 0):

- **Length sanity** — 1,300 ZnF records cross-checked; **0** spans too short to
  hold the fingers they kept (threshold 20 aa/finger; a C2H2 finger is ~23 aa
  plus a ~5 aa TGEKP linker). The reference is the count of fingers the span
  actually kept, re-derived from the cached feature table — `zinc_finger_count`
  counts every `ZN_FING` including PHD and RING, so it is an upper bound.
  *(That over-count is a live bug elsewhere: `af3.py` adds one Zn ion per
  `zinc_finger_count`, so a PHD-bearing chassis gets too many ions.)*
- **Known anchors** — **7/7** exact, and byte-identical to the by-hand
  extraction used for the fib→β constructs (plan step 6).
- **KRAB audit** — **422** records dropped a KRAB; **0** where the span reaches
  into one. Three (ZNF132, ZNF418, ZNF256) have a KRAB annotation that overlaps
  their first kept finger — an annotation conflict with no correct trim, since
  pulling the boundary back would cut into the finger. Those are reported and
  flagged for review rather than failing the gate; the honest response is to
  exclude the chassis at design time.

### Review flags

287 records set `dbd_review_flag` — an excluded domain or effector-annotated
region overlaps the kept span, or the span is mostly not DNA-binding feature.
Real ones: **ZNF541/TRERF1** (ELM2, a NuRD recruiter — an unintended closer),
**ATF7** (transactivation domain), **ZHX2/ZHX3**, **MGA** and **TRPS1** (span
coverage). Conservative ones: **REST**, **YY1**, **FOS**, **HINFP**, where the
flagged region is itself required for DNA binding.

**[corrected 2026-11-18]** This list named **NSD2/KMT2C** first. Both carry
`dbd_review_flag = 0` — NSD2 extracts `dna_bind` 447-527, KMT2C 28-52, and
their SET/PHD/PWWP domains fall *outside* the kept span, so flag 0 is correct
behaviour. The error propagated from here into the first draft of
`DBD_EXTRACTION_FIXES.md` (F4.2). KMT2C has a different problem: its span is
25 aa of a 4,911-aa protein, which F2b now flags.
The flag is a review prompt, not a verdict.

### Consumers

`sequence_aa` is untouched; the original is also copied to `sequence_full_aa`.
`fuse.py` now prefers `dbd_sequence_aa` and prints a loud warning when it has to
fall back (per-gene if the column exists but a record has no span; a
whole-input warning if the column is missing entirely, i.e. a library exported
before this change). `output.py --include-sequence` carries `dbd_sequence_aa`,
`dbd_start`, `dbd_end`, `dbd_extraction` and `dbd_review_flag` through Step 1.

### Still open

- The **262** `none` records and the `Others` family need per-family keep rules.
  (**[corrected 2026-11-18]** this line said 196, a figure from an earlier run;
  `:188` in this same file already tabulated 262. F9 splits these into
  `none` abstentions and retryable `no_features` fetch failures.)
- The **287** flagged records need curation; ELM2-bearing chassis should
  probably be excluded from opener designs outright. (**[corrected
  2026-11-18]** this line said 32, the pre-rework count that `:236` records as
  "32 → 287". As of F1 in eCR_predictor they are excluded by default rather
  than awaiting curation.)
- `dbd_length_aa` is not yet a design-time constraint. ZNF850 extracts to 900 aa,
  so the fusion is ~980 aa (~2.9 kb) — real against a ~4.7 kb AAV cargo.
