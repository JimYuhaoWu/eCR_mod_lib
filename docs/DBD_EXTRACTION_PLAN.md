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

2,818 of 3,248 DBD rows are extractable; **2,622 (93%) got a span**.

| `dbd_extraction` | n | meaning |
|---|---|---|
| `znf_array` | 1,121 | contiguous first-finger → last-finger block |
| `dna_bind` | 787 | a `DNA_BIND` feature |
| `domain` | 440 | a DNA-binding `DOMAIN` (bHLH, bZIP, MBD, …) with no `DNA_BIND` |
| `bipartite` | 274 | binding feature + partner domain (e.g. HNF1A POU + homeobox) |
| `none` | 196 | no usable feature — full-length record left honest, per rule 3 |
| *(NULL)* | 430 | no accession or no sequence; not touched |

`domain` was added to the four kinds the draft listed: bHLH, bZIP, MBD, RHD, MH1
and others annotate their binding module as `DOMAIN` with no `DNA_BIND` feature
at all, so without it 214 bHLH and 108 bZIP records would have fallen to `none`.

By family: zf-C2H2 921/927, Homeobox 378/379, bHLH 208/208, bZIP 100/105,
HMG 98/101, Fork_head 94/94, ZBTB 83/83, MYB 52/52, T-box 34/34. The families
that stay at `none` are `Others` (74/160), MH1, STAT, Tub, AP-2, GTF2I, AF-4 —
each needs a curated keep rule, and guessing one is worse than reporting the gap.

### Keep list is data-driven, and short names match exactly

`--survey` prints every unclassified `DOMAIN` description, which is how the keep
list was built rather than assumed. It also caught a real hazard: substring
matching on short domain names is wrong — `tea` is inside *Protease* and `rel`
inside *Related*. Short names now live in `KEEP_DOMAIN_EXACT` and must match the
whole description (with any trailing repeat index stripped).

The survey confirms the drops are right too: the most common dropped domain is
**NR LBD** (96) — the nuclear-receptor *ligand*-binding domain, correctly not
carried into a fusion — followed by SET, PAS, Orange, ELM2, LIM, SH2.

### Validation

All three gates pass (`--validate`, exit 0):

- **Length sanity** — 1,004 ZnF records cross-checked against the independently
  stored `zinc_finger_count`; **0** spans too short to hold the fingers claimed
  (threshold 20 aa/finger; a C2H2 finger is ~23 aa plus a ~5 aa TGEKP linker).
  Only *too short* fails: a span longer than the pitch is legitimate, because a
  split array legitimately spans the gap between two clusters.
- **Known anchors** — **7/7** exact, and the extracted sequences are byte-identical
  to the by-hand extraction used for the fib→β constructs (plan step 6).
- **KRAB audit** — **422** records dropped a KRAB domain; **0** retain one inside
  the kept span. This is the failure that motivated the work, and it is now
  measured rather than assumed: 422 chassis, not the two the run happened to pick.

### Review flags

32 records set `dbd_review_flag`, i.e. an excluded domain or an effector-annotated
region lies *inside* the kept span and is therefore carried along. Real ones:
**NSD2/KMT2C** (a SET methyltransferase inside the span), **ZNF541/TRERF1/Zfp541**
(ELM2, a NuRD-recruiting module — an unintended closer), **ATF7** (transactivation
domain), **ZHX2/ZHX3** ("required for repressor activity"), **BAZ2A/BAZ2B** (DDT).
Conservative ones: **REST**, **YY1**, **FOS**, where the flagged region is itself
required for DNA binding. The flag is a review prompt, not a verdict.

Note the shape of this result: the KRAB gate is clean because KRAB is always
N-terminal to the ZnF array, but ELM2 and SET sit *inside* their spans, so they
cannot be removed by trimming at all. Those chassis need excluding at design
time, which is what the flag is for.

### Consumers

`sequence_aa` is untouched; the original is also copied to `sequence_full_aa`.
`fuse.py` now prefers `dbd_sequence_aa` and prints a loud warning when it has to
fall back (per-gene if the column exists but a record has no span; a
whole-input warning if the column is missing entirely, i.e. a library exported
before this change). `output.py --include-sequence` carries `dbd_sequence_aa`,
`dbd_start`, `dbd_end`, `dbd_extraction` and `dbd_review_flag` through Step 1.

### Still open

- The 196 `none` records and the `Others` family need per-family keep rules.
- The 32 flagged records need curation; ELM2- and SET-bearing chassis should
  probably be excluded from opener designs outright.
- `dbd_length_aa` is not yet a design-time constraint. ZNF850 extracts to 900 aa,
  so the fusion is ~980 aa (~2.9 kb) — real against a ~4.7 kb AAV cargo.
