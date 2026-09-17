# Next modification: expand the DBD library with Codebook motifs

**Status:** proposed 2026-09-04, **revised 2026-09-17** after the paper was located
and read and the consumer code was traced. Follows the completed DBD extraction
([DBD_EXTRACTION_PLAN.md](DBD_EXTRACTION_PLAN.md)).

**Source:** Codebook Consortium, *An expanded codebook of human transcription
factor DNA-binding specificity*, **Nature, August 2026**,
[10.1038/s41586-026-10798-9](https://www.nature.com/articles/s41586-026-10798-9)
(open access: [PMC13538051](https://pmc.ncbi.nlm.nih.gov/articles/PMC13538051/)).
Companion benchmarking paper: *Cross-platform motif discovery and benchmarking…*,
Commun Biol 2025, [s42003-025-08909-9](https://www.nature.com/articles/s42003-025-08909-9).
Motif explorer: [mex.autosome.org](https://mex.autosome.org). Bulk data: Zenodo
concept DOI [10.5281/zenodo.8327372](https://doi.org/10.5281/zenodo.8327372),
current version [15667805](https://zenodo.org/records/15667805) (2025-06-15).

> **Revision note.** The first draft of this plan carried a "verify before
> implementing" caveat on the paper's numbers. They have now been verified, and
> three of the assumptions built on them were wrong. Corrections are marked
> **[corrected]** below. The headline change: the Codebook cannot close the gap
> this plan was written to close — it can close at most a quarter of it.

---

## 1. What the paper actually contains  **[corrected]**

| | |
|---|---|
| proteins assayed | 393 (332 putative/poorly-characterised TFs + 61 controls) |
| experiments | 4,237, across five platforms |
| platforms | HT-SELEX, GHT-SELEX, SMiLE-seq, PBM, ChIP–seq (inducible tagged transgenes, HEK293) |
| **Codebook TFs with a motif** | **177 of 332 (53%)** |
| distinct motifs among them | 129 |
| control TFs with a motif | 58 of 61 |
| motifs per TF | **one** representative PWM, selected by an expert committee |
| confidence tiers | **none** |

Two consequences the earlier draft got wrong:

**[corrected] The ceiling is 177, not 673.** The previous draft said the
Codebook's purpose "is exactly this set" and implied it closes the 673-record
motif gap. It cannot. 177 TFs received motifs, so even if every one of them were
a library record we currently carry motif-less, **the gap closes by at most 26%**
— and the real figure is the intersection, which is smaller. Sizing the work as
"close the 673" was wrong; sizing it as "recover up to 177, measure which" is
right. This does not make the work unattractive — 177 new designable chassis is
a large fraction of what we could design today (964 human motif-bearing DBDs) —
but it must be planned as an increment, not a fix.

**[corrected] There are no confidence tiers, so `validation_level` cannot be
mapped to them.** The earlier draft instructed "decide this explicitly against
the paper's own confidence tiers"; the paper has none. Curation was an expert
committee judging AUROC, Jaccard, peak overlap and family plausibility, and
publishing a single representative PWM. The mapping rule must therefore be built
from something that does exist in the data — see §5.

**[resolved risk] One motif per TF.** A serious schema worry going in was
one-to-many: if a TF carried several motifs, a single `jaspar_id`-shaped column
could not hold them and the join key would have to become a separate table. The
paper's "single representative PWM" settles it — **one motif per TF, so the
existing one-column shape survives.** Recorded because it was the most expensive
thing that could have been true.

**The 155 TFs that yielded no motif are themselves a result worth ingesting.**
They were assayed by five platforms across 4,237 experiments and still produced
nothing reproducible. Recording that as an explicit "assayed, no motif"
annotation stops us, and anyone reading the library later, from treating those
records as merely un-annotated and waiting for a motif that the field has already
looked hard for. This is cheap and it is the sort of negative result a library
should carry.

---

## 2. The gap, re-measured (2026-09-17, live library)

| | count |
|---|---|
| `type='DBD'` rows, all organisms | 3,248 |
| human `type='DBD'` | 1,637 |
| — with a `jaspar_id` | 964 |
| — **without any motif** | **673** |
| ⤷ of those, with an extracted DBD span | 553 |
| &nbsp;&nbsp;&nbsp;&nbsp;⤷ span clean (`dbd_review_flag=0`) | **468** ← immediately actionable |
| &nbsp;&nbsp;&nbsp;&nbsp;⤷ span review-flagged | 85 |
| ⤷ with no `uniprot_id` at all | 7 |
| mouse `type='DBD'` | 1,611 |
| — with a `jaspar_id` | 213 |

A DBD with no motif is not merely under-annotated — it is **undesignable**.
`target.py:307` drops motif-less DBDs outright (unlike `cli.py`, which keeps them
with `motif_score = NA`), because a DBD that cannot be placed on a locus by
sequence cannot be assigned a binding site. So those 673 records cost nothing and
contribute nothing today.

**468**, not 553, is the immediately actionable number: motif + a *clean*
extracted boundary is what `target.py` → `cli.py` → `fuse.py` needs end to end,
and an 85-record flagged tail should not be counted as ready.

**The 7 records with no `uniprot_id`** are the concrete unmappable set. The
earlier draft promised to report unmappable TFs rather than guess; this is the
floor of that count, known in advance.

### The organism problem the earlier draft never mentioned

**The Codebook is human-only, and our flagship runs are mouse.** Mouse DBDs
outnumber human in the library (1,611 vs 1,637, near-parity) and are far worse
covered — only **213 of 1,611 (13%)** carry a motif, against 59% for human. The
MEF→mES work runs on mouse.

So: does a human Codebook motif transfer onto the mouse orthologue's row? For a
conserved DBD it is scientifically defensible — the recognition helix is what
determines specificity, and for most families it is identical between mouse and
human. But it is a **decision, not a default**, and it must be recorded per row
(`motif_source='codebook-human-ortholog'` or similar) so that a design resting
on a transferred motif is visibly resting on one. Doing it silently would put
un-marked cross-species inference into `design_score`.

Recommendation: **ingest human-to-human first, ship it, and treat orthologue
transfer as a separate follow-up with its own gate** (percent identity over the
extracted DBD span, which we now have for both organisms — that is exactly the
measurement that licenses or refuses the transfer, and we could not have made it
before the extraction work).

---

## 3. Which artefact is authoritative — decide before downloading  **[corrected]**

This was not in the earlier draft at all, and it is the easiest way to get the
whole exercise subtly wrong.

The Zenodo record is the **GRECO-BIT/MEX benchmarking** dataset, and its motif
files rank motifs *by benchmark metric*:

| file | size | content |
|---|---|---|
| `MEX_top1.zip` | 0.21 MB | best-scoring motif per dataset/TF **by metric** |
| `MEX_top20.zip` | 3.58 MB | top 20 per TF |
| `MEX_motifs_complete.tgz` | 51.87 MB | all PFMs and PCMs |
| `metadata_complete_motif.zip` | 2.58 MB | tool + dataset provenance per motif |
| `metadata_complete_dataset.tsv` | 6.97 MB | data-slice ID → Codebook experiment ID |

**`MEX_top1` is not the paper's representative PWM.** "Top by benchmark metric"
and "selected by the expert curation committee" are different selection rules
over the same candidate pool, and they disagree for some TFs by construction —
that disagreement is a large part of what the companion benchmarking paper is
about. The authoritative set is **Supplementary Data 1 of the Nature paper**
(representative PWMs for all 1,421 TFs). Ingest that; use MEX for provenance, for
the per-motif metrics, and as a cross-check.

Two further traps in the same place:

- **Pin the Zenodo *version*, not the concept DOI.** A 2025 update replaced the
  files for ten TFs — CHCHD3, CPXCR1, MYRFL, PRDM5, TIGD4, TIGD5, TTF1, ZBED5,
  ZBTB5, ZNF578 — whose datasets were partially mislabelled in the original
  upload. An unpinned reference silently resolves to a different dataset than the
  one we validated against. Record the version DOI and the file checksums in
  `download_manifest.json` like every other source.
- **The README lists `MEX.GHTS.tar` and `MEX.HTS.tar`, but the record's file list
  does not contain them.** Either the record is incomplete or the archive names
  changed. Irrelevant to us (we want PWMs, not raw reads) but it means the README
  is not a reliable index of what is actually there — check the file listing
  through the API, not the README.

- **`Supplementary Data 1` covers 1,421 TFs, not 177.** It is a full collection,
  most of which duplicates motifs we already resolve through JASPAR. The ingest
  must filter to the Codebook set (and the controls, see §6) rather than
  wholesale-importing 1,421 motifs and creating a second, competing motif source
  for TFs that already have one.

---

## 4. The architectural decision: `jaspar_id` is the join key

Unchanged in substance, but the **blast radius is larger than the earlier draft
stated**, now that the consumers have been traced.

A motif today *is* a `jaspar_id`, resolved by `scan.fetch_pwm` in a fixed order:
local `jaspar_cache/<id>.jaspar` → BioPython `JASPAR2020` → JASPAR REST. A
Codebook motif has no JASPAR accession, and **the paper does not state any JASPAR
or CIS-BP deposition** — distribution is Zenodo + MEX. That removes the evidence
for option (a).

**(a) Wait for JASPAR deposition.** Zero code, zero risk, motifs arrive through
the existing path on a future release. But there is no announced deposition, so
the delay is not merely unknown, it is unevidenced; and partial deposition gives
partial coverage with no way to tell which is which. **Now the weakest option**,
not the safe one.

**(b) Add `motif_id` + `motif_source`, keep `jaspar_id` as-is.** `jaspar_id`
keeps its current meaning exactly, so nothing reading it changes behaviour.
**Recommended.**

**(c) Write Codebook motifs into `jaspar_id` with a prefix.** Cheapest to write
and worst to live with: it makes `jaspar_id` mean two things — the same class of
mistake as `sequence_aa` holding whole proteins, which the last cycle was spent
undoing. Rejected, and recorded so it is not rediscovered as a shortcut.

### What (b) actually touches  **[expanded]**

| file | why |
|---|---|
| `eCR_mod_lib/scripts/schema.py` | two columns via `MIGRATIONS` (the mechanism is already there) |
| `eCR_mod_lib/scripts/01_fetch_dbd.py`, `05_validate.py` | read/validate `jaspar_id` |
| `ecr_predictor/scan.py` | `fetch_pwm` gains a `motif_source` branch — **the single resolution point, and the only place where the change is clean** |
| `ecr_predictor/prefetch.py` | **re-implements the JASPAR REST fetch itself** rather than delegating to `scan.fetch_pwm`. Either the branch is duplicated here or prefetch is refactored to delegate. **Refactor** — a second copy of resolution logic is how the two paths drift. |
| `ecr_predictor/fimo.py` | a **second, independent motif path**: `_jaspar_to_meme_block` emits `MOTIF <jaspar_id>` and `fimo_pvalue` is mapped back by `jaspar_id`. Untouched, FIMO validation silently skips every Codebook motif. |
| `ecr_predictor/target.py` | `SITE_COLS`/`GAP_COLS` carry `jaspar_id`/`best_jaspar_id`; `load_burdens` joins the off-target TSV on `jaspar_id` |
| `ecr_predictor/offtarget.py` | writes `jaspar_id` as a **TSV column header** — a file-format contract, so old and new bundles must both parse |
| `ecr_predictor/output.py` | `--include-sequence` column list (already extended once for the `dbd_*` columns) |

Six modules across two repos, two of them motif paths the earlier draft did not
know about. Still the right option — but it is a coordinated release, and
`fimo.py` and `prefetch.py` are where an incomplete job will hide.

---

## 5. `validation_level`, without tiers to lean on  **[corrected]**

`validation_level` feeds `annotation_confidence`, one of the two independent
scores the predictor reports, so mapping every Codebook row to one value for
convenience would flatten a signal we deliberately keep separate.

With no tiers in the paper, derive the level from what the data does contain —
**which platforms supported the motif**:

- motif from *in vitro* binding assays only (HT-SELEX / GHT-SELEX / SMiLE-seq /
  PBM) → **`motif-only`**. This is a direct binding measurement and nothing more.
- motif additionally supported by the **ChIP–seq** arm (bound in cells) →
  **`ChIP-validated`**, which is what that level means everywhere else in the
  library.

`metadata_complete_motif.zip` carries the per-motif tool and dataset provenance
needed to make this assignment mechanically rather than by hand. Do not reach for
`screen-validated`: no screen was run.

---

## 6. Step 0 — measure the overlap before writing ingest code

Unchanged in intent, sharpened by §1: the number that decides the size of the
work is *not* 673, it is |Codebook 177 ∩ library|. Map the Codebook TF list to
UniProt accessions and intersect. Report five numbers:

1. Codebook TFs already in the library **with** a motif — do we agree, or does
   the Codebook motif differ from the JASPAR one? A disagreement is a finding,
   not a merge conflict to paper over. **Needs a precedence rule decided in
   advance** (see below).
2. Codebook TFs in the library **without** a motif — **the target set, and the
   real size of this project.**
3. Codebook TFs **not in the library** — new DBD records, which also need a DBD
   span extracting, so they land back in `06_extract_dbd.py`.
4. Codebook TFs we cannot map to an accession — reported, not guessed.
5. Of the target set, how many have a **clean extracted span** — the
   end-to-end-buildable count, the honest headline.

**Precedence rule, to decide now rather than at merge time:** where a TF has both
a JASPAR and a Codebook motif, **JASPAR wins for `jaspar_id` and the Codebook
motif is still recorded** in `motif_id`/`motif_source`. Nothing that reads
`jaspar_id` changes behaviour, existing designs do not move, and the disagreement
stays measurable instead of being resolved by import order. An ingest that
overwrites a curated `jaspar_id` is the `04_build_library.py` erasure bug in a
new costume.

This step is about half a day and it determines everything after it. Do not skip
to the ingest.

---

## 6b. Step 0 — RESULTS (measured 2026-09-17)

Step 0 is done. It changes the recommendation, so read this before §7.

**Source artefacts used.** Zenodo 15667805 → `MEX_top1.zip` (0.21 MB, 236 files,
one representative motif per TF) and `metadata_complete_motif.zip` (2.58 MB →
`complete_metadata_motif_15-06-2025.tsv`, 202,868 candidate motifs). The metadata
carries **`Approved`** (dataset passed curation, yes/no), `Overall rank`, `Tool`
and `Train data` per motif. 394 TFs assayed; 159,063 approved motif rows;
**236 TFs with an approved motif** (= 177 Codebook + 58 controls + 1).

> **[corrected] The paper has no confidence *tiers*, but the data does carry a
> binary `Approved` curation flag plus a per-motif `Train data` platform.** §5's
> rule is therefore implementable mechanically after all, straight off this table
> — better than the draft assumed.

> **[corrected] The motifs ship as `.ppm` — position *probability* matrices.**
> This confirms the §7 normalisation hazard as fact rather than caution:
> `scan.py` does `motif.counts.normalize(pseudocounts=0.5)`, which assumes
> counts. Applied to probabilities summing to 1.0 per column, a 0.5 pseudocount
> **swamps the matrix entirely** and every score is garbage — silently.
> `MEX_motifs_complete.tgz` ships PCMs as well; use those, or restore counts from
> the PPM against a stated effective count.

### The intersection

| | count |
|---|---|
| Codebook TFs with an approved motif | 236 |
| **in our library, already carrying a motif** | **153** |
| ⤷ of which the motif is **already the Codebook motif**, arrived via JASPAR | **96** |
| ⤷ independent motif → genuine comparison case | 57 |
| **in our library, motif-less → TARGET SET** | **52** |
| ⤷ clean extracted span (buildable end to end) | **35** |
| ⤷ span review-flagged | 12 |
| ⤷ no span | 5 |
| ⤷ **already in JASPAR today** (free, no ingest) | 5 |
| ⤷ genuinely absent from JASPAR | **47** |
| not in the library at all → new DBD records | 31 |

### The finding: the Codebook is already reaching us through JASPAR

**96 of the 153 overlap motifs are Codebook-derived**, sourced to PMIDs
39605368 / 39605320 / 39605530. And **101 of the 143 human `MA≥2000`
(JASPAR-2024-era) motifs in the library — 71% of that entire bin — belong to
Codebook TFs.** JASPAR has been absorbing the Codebook, and our ordinary JASPAR
path has been quietly collecting it for us.

**[corrected] This reverses §4's verdict on option (a).** The revision called
"wait for JASPAR deposition" the weakest option because the paper announces no
deposition. The paper doesn't, but JASPAR did it anyway — and has already
delivered 96 of the 236. Option (a) is not a gamble on a future release; it is a
process **already two-thirds complete**.

So the true value of building the `motif_id` / `motif_source` architecture is not
673 records, and not 177 TFs. It is the **47-TF residual** JASPAR has not yet
absorbed — of which ~26 have a clean span — and that residual **shrinks with
every JASPAR release**. Six modules across two repos, a threshold recalibration
and a specificity re-scaling, to capture 47 TFs that are arriving anyway, is poor
value. **Recommendation: defer the custom ingest.**

### What to do instead, in order of value per unit of work

1. **Pay down JASPAR annotation debt first — no architecture at all.** Five of
   the 52 target TFs have had JASPAR motifs for years: RXRA (`MA0512`), SOX5
   (`MA0087`), HMGA1 (`MA2124`), PRDM5 (`MA1999`), ZFAT (`UN0676`). `MA0512` and
   `MA0087` are not obscure. Their absence from our rows is **our annotation
   debt, not the field's gap** — part of the 673 is simply un-recorded. Sizing
   this across all 673 is running now (§6c); if it is large, it is by far the
   cheapest designability win available and needs no schema change, no new motif
   source, and no recalibration.
2. **Mouse orthologue transfer — the biggest payoff, and still no new motif
   source.** 30 of the 52 target TFs have mouse rows and **25 of those are
   motif-less**. Mouse coverage is 213/1,611 (13%) against 59% for human, and the
   flagship runs are mouse. Transfer reuses existing `jaspar_id` values on mouse
   rows, gated on percent identity over the extracted DBD span — a measurement
   the extraction work has only just made possible. This is worth more to the
   actual pipeline than the entire human Codebook ingest.
3. **Re-run `prefetch.py` against current JASPAR** and re-measure the residual.
   Cheap, and it tells us the absorption rate — which is the number that decides
   whether the residual is ever worth custom code.
4. **Only then reconsider the ingest**, if the residual is both large and stable.

### Gate 2's fixture is 57 TFs, not 153 — and the difference is circularity

The 96 TFs whose JASPAR motif *is* the Codebook motif cannot validate an ingest
of Codebook motifs: comparing a motif to itself agrees by construction, and a
gate that cannot fail is not a gate. The usable fixture is the **57 TFs with an
independently derived motif** (ENCODE, ReMap, older JASPAR). Anyone running
Gate 2 on all 153 would get a spectacular agreement statistic that means nothing
— and this is precisely the trap the "gates that can actually fail" principle
exists to catch.

### Caveats on these numbers

- Matching is by **gene symbol**, uppercased, against `gene_symbol` falling back
  to `name`. Symbol aliases will cause some undercount in the "31 not in library"
  row; several of those (KDM2A, PHF21A, TET3, TCF20, MBNL2, FBXL19, GLYR1,
  DNTTIP1) are chromatin proteins that AnimalTFDB plausibly does not class as
  TFs, so they are genuinely new records rather than alias misses. A UniProt
  accession join would firm this up and is the one refinement worth making.
- `AC092835`, `ZNF705EP`, `ZNF788P` are a clone-based locus and two pseudogenes.
  Do not create DBD records for them without checking they encode a real protein.
- `Approved` is a **dataset-level** curation flag, not a per-motif quality score;
  do not read it as the expert committee's representative-PWM selection.
- MEX top-1 selection is **by benchmark metric**, so this intersection is a TF
  *list* measurement. It is valid for sizing — the TF roster is the same — but
  the PWMs to ingest must still come from the paper's Supplementary Data 1 (§3).

## 6c. The real finding: `01_fetch_dbd.py` assigns motifs by the wrong taxon

Step 0's cheapest side-check turned out to matter more than the Codebook. All 673
motif-less human DBD rows were queried against current JASPAR by gene name:

| | count |
|---|---|
| motif-less human DBD genes | 673 |
| **JASPAR has a motif we never recorded** | **94** |
| ⤷ CORE collection (not UNVALIDATED) | 89 |
| ⤷ **with a clean extracted span → designable *and* buildable today** | **74** |
| ⤷ span flagged / no span | 13 / 7 |
| genuinely absent from JASPAR | 579 |

**94 is twice the Codebook's 47-TF residual, and it costs no new architecture at
all** — no `motif_id`, no `motif_source`, no threshold recalibration, no
coordinated two-repo release. 74 of them are immediately buildable.

### Why they are missing, proven

`fetch_jaspar_id` (`scripts/01_fetch_dbd.py:230`) queries
`/matrix/?name={gene}&tax_id={taxon_id}`. **JASPAR tags a profile with the
species of the protein the experiment used, not the species the profile applies
to.** CORE profiles are vertebrate-wide. So:

```
/matrix/?name=RFX6&tax_id=9606              → count: 0
/matrix/?name=RFX6&tax_id=10090             → MA1724.1, MA1724.2  (CORE)
/matrix/?name=RFX6&tax_group=vertebrates    → MA1724.1, MA1724.2  (CORE)
```

`MA1724` was derived from mouse Rfx6, so the human query returns nothing and
human RFX6 is recorded as motif-less — while the mouse row carries `MA1724.1`.
The library confirms this pattern exactly, same gene, same profile:

| gene | human `jaspar_id` | mouse `jaspar_id` |
|---|---|---|
| RFX6 | — | MA1724.1 |
| ISL1 | — | MA1608.1 |
| PTF1A | — | MA1619.1 |
| NKX2-1 | — | MA1994.1 |
| SOX17 | — | MA0078.1 |
| RXRA | — | MA0512.1 |
| FOXO1 | — | MA0480.1 |

**93 of the 94 have a mouse row, and 91 of those already carry the motif.** The
gap is not missing data — the same profile is already in the library on the other
species' row. Motif coverage has been split by which species the JASPAR
experiment happened to use (human 59%, mouse 13%), which is an artefact of the
query, not biology.

### This has been silently constraining the flagship deliverable

The 94 include **RFX6, ISL1, PTF1A, NKX2-1, NR5A2, FOXO1, SOX17, MECOM, WT1** —
pancreatic and endocrine identity TFs. `target.py:307` drops motif-less DBDs, so
**every fib→β design run so far has been unable to consider RFX6, ISL1 or PTF1A
as a DBD chassis**, for a reason that has nothing to do with their biology. Two
of them (RFX6, ISL1) are on the carried TODO list of anchors for the model
stratification re-run, where the same absence would have skewed the comparison.

### Two further defects in the same function

1. **`results[0]` takes an arbitrary version, and the library is mostly stale.**
   RFX6 returns `MA1724.1` *and* `MA1724.2`; the script keeps whichever the API
   lists first. Measured on a 120-motif sample of the human rows: **82 of 120
   (68%) have a newer version in JASPAR than the one we hold** — e.g.
   `MA0028.1 → .3`, `MA0798.1 → .3`, `MA0761.1 → .3`, `MA0500.2 → .3`. Library-wide,
   1,074 of 1,177 `jaspar_id` values end in `.1`, which is the signature of
   taking the first result rather than the current one.

   This is **not cosmetic**: a JASPAR version bump is a re-derived matrix, so
   every `motif_score`, every off-target burden and every `design_score` in the
   project has been computed from superseded PWMs for roughly two thirds of the
   library. Fixing it will move existing scores — which is exactly the rank-churn
   problem §8 Gate 4 describes, arriving from a different direction and before
   any Codebook work. Pin the specificity scale first (§8, option (i)), then
   re-score, so the shift is measurable rather than silent.
2. **The empty result is cached and never revisited.** The cache key is
   `data/raw/jaspar/{gene}_{taxon_id}.json`, and a cached `results: []` returns
   `None` forever. The bug froze itself in place, which is why it survived every
   subsequent re-run and why `--force` alone will not clear it.

### Fix

Query `tax_group=vertebrates` (or query both taxa and prefer same-species),
match the gene name **exactly** rather than trusting the API's filter, prefer the
highest version, and invalidate the negative cache entries. Then re-run 01 → 04
and re-measure. This is a small, self-contained change to one function with a
ready-made regression fixture: the 94 genes must gain a motif, and the 7 rows
with no `uniprot_id` must stay empty.

**This supersedes the Codebook as the next piece of work.** It is smaller,
strictly cheaper, twice the yield, and it repairs a defect rather than adding a
source on top of one.

### Revised order of work

1. **Fix `fetch_jaspar_id`** — taxon filter, exact name match, highest version,
   negative-cache invalidation. Yield: 94 human rows gain a motif (74 buildable),
   plus a corrected matrix for ~68% of the motifs we already hold.
2. **Pin the specificity scale** (§8 Gate 4, option (i)) *before* re-scoring, so
   the ranking shift from step 1 is measured instead of absorbed.
3. **Re-run 01 → 04, re-score, report rank churn** against the existing fib→β
   bundle. RFX6/ISL1/PTF1A becoming available may change which chassis the
   design run picks — that is the point, and it needs to be visible.
4. **Mouse orthologue transfer**, gated on span identity (§6b item 2).
5. **Reconsider the Codebook ingest** against whatever residual remains. On
   today's numbers that is 47 TFs and falling.

## 7. Ingest

Follows the established `02_seed_ed.py` pattern — a declarative source entry,
supplementary files in `data/manual/`, a `literature/papers.yaml` record, a
provenance row.

1. `data/manual/codebook_2026_*` — Supplementary Data 1 plus the MEX metadata,
   checksummed into `download_manifest.json`, **with the pinned Zenodo version
   DOI and the paper DOI**.
2. `scripts/07_fetch_codebook.py` — parse to PFM, write
   `data/raw/codebook/<motif_id>.pfm`, emit `data/processed/codebook_raw.tsv`.
   A separate script, not a branch in `01_fetch_dbd.py`: different shape,
   different failure modes, and 01 is already slow.
3. `04_build_library.py` gains `codebook_raw.tsv` as a source.

**[corrected] `motif_id` has to be minted — the Codebook has no motif ID
scheme.** The earlier draft wrote `<motif_id>.pfm` as though ids came with the
data; they do not. Motifs are identified by **TF protein name** (TPRX1, DACH1).
So mint one deterministically and namespaced — `CB.<TF>.<assay-or-rev>` — with
three properties: it must not collide with a JASPAR accession pattern, it must be
stable across re-ingests (so re-running does not renumber every motif and dirty
every downstream bundle), and it must be derivable from the source row alone
rather than from a counter. A gene symbol alone is not enough: symbols change,
and the same symbol will need a second motif the moment an orthologue transfer or
a revised PWM arrives.

**PFM normalisation is a real hazard, not bookkeeping.** `scan.py` does
`motif.counts.normalize(pseudocounts=0.5)`, which assumes `counts` are *counts*.
A Codebook matrix distributed as probabilities makes a 0.5 pseudocount
overwhelming; one distributed as raw SELEX read counts makes it negligible. Both
silently change every score. Normalise to a stated effective count at ingest and
record which convention the source used. `MEX_motifs_complete.tgz` ships both PFM
and PCM, which is what makes this checkable.

---

## 8. Validation gates

Same principle as the extraction work: gates that can actually fail. Two of the
four have changed materially.

**Gate 1 — Round-trip.** Every ingested PFM re-reads to the same matrix, sums to
a consistent count per column, and has the length the source claims. Extended:
the matrix must survive `scan.fetch_pwm` → `motif.counts.normalize` → `log_odds`
without NaN or inf, and `best_site_score` on the motif's own consensus must
return ≈1.0. A motif that cannot score its own consensus is broken, and this is
the cheapest possible detector.

**Gate 2 — Known-motif agreement, anchored on the paper's own controls.**
**[expanded]** The Codebook assayed **61 control TFs and obtained motifs for 58**
— TFs whose specificity was already known, run through the same five pipelines.
That is a purpose-built positive control, and it is the regression fixture this
plan was missing: the extraction work's single most useful feature was the seven
hand-checked `ANCHORS`, and there was no equivalent here. Use the controls that
also have a JASPAR motif as a **committed fixture**, comparing by Pearson on
aligned PFM columns or Tomtom. High agreement validates the whole ingest path
end to end; disagreement on a *control* is a bug in our ingest, not a biological
finding — which is exactly what makes it a usable gate. Disagreement on a
Codebook TF, by contrast, is a real result about assay differences and should be
reported per TF rather than averaged away.

**Gate 3 — Threshold comparability.** **[new — this is the gap that matters
most.]** `offtarget.DEFAULT_THRESHOLD = 0.8` is a **fixed** cut on
`best_site_score`, which normalises log-odds by the motif's own theoretical
maximum. That normalisation makes scores dimensionless but **not comparable
across motifs**: a short, low-information motif reaches 0.8 far more often by
chance than a long, high-information one. Every motif in the library today came
from JASPAR, so this bias has been roughly constant and has cancelled out of our
comparisons. Codebook motifs come from different assays with a different length
and information-content distribution, so **importing them at the same fixed
threshold mixes two scales in one ranking** — and it does so invisibly, because
nothing in the current output reports motif length or IC.

The gate: before any design run, report the length and information-content
distribution of the Codebook motifs against the incumbent JASPAR set, and the
per-motif hit rate on the accessible background at threshold 0.8. If the
distributions differ, a fixed normalised cut must be replaced by a **per-motif
threshold at matched expected false-positive rate** — `fimo.py` already carries
the p-value machinery (`FIMO_PVALUE_THRESHOLD = 1e-4`) to calibrate against.
Skipping this gate is how the designability gate below returns a beautiful number
for the wrong reason: low-IC motifs that "hit" everywhere.

**Gate 4 — Designability delta, and the ranking churn it causes.**
**[corrected — the earlier version of this gate was wrong.]** Re-run `target.py`
on an existing run bundle and report how many previously undesignable regions
gain a candidate DBD. That is the number that says whether the work paid off, and
it is the honest headline: if the gap report barely shrinks, the motifs are not
reaching the loci we nominate.

The earlier draft then said to "confirm the burden distribution of existing DBDs
is unchanged (it must be — burden is per motif)". **Half right, and the wrong
half is load-bearing.** Burden *is* per motif and will not move. But
`target.compute_specificity` divides burden by **the median positive burden
across the DBDs in that run**, so adding several hundred motifs **moves the
median, changes every existing DBD's `specificity`, changes `design_score`
(= motif_score × specificity), and reshuffles the Pareto front** — including for
designs we have already built and sent to AF3. The gate as originally written
would have passed while the ranking silently moved underneath it.

So the gate is: report **rank churn** on a fixed bundle — how many existing DBDs
change Pareto status, and the largest `design_score` change — and decide
deliberately whether to (i) pin the specificity scale to a recorded constant so
rankings stay comparable across library versions, or (ii) accept the re-scaling
and re-rank everything at once, re-running the affected designs. Either is
defensible; drifting into (ii) without noticing is not. **Option (i) is the
recommendation**: a library-relative denominator that moves whenever the library
grows means no two runs are ever comparable, which is a defect independent of the
Codebook.

---

## 9. Cost

Ingest and PFM handling: local, minutes. `Supplementary Data 1` and
`MEX_top1.zip` are sub-megabyte; only `MEX_motifs_complete.tgz` (52 MB) is worth
a deliberate download, and the multi-gigabyte raw-read archives are not needed
at all.

The real cost is `offtarget.py`, which scans the accessible background per motif
at roughly **45 minutes per 32-way shard per motif set**. Bounded here by the
Codebook's own ceiling — **at most 177 new motifs, and only those in the
intersection** — so this is a single overnight batch, not an open-ended one.
Note it buys nothing for the mouse runs until the orthologue question in §2 is
decided.

---

## 10. Sequencing against the open extraction work

- **262 `none` records** (the earlier draft said 196 — **[corrected]**, measured
  2026-09-17) and the `Others` family need per-family keep rules. If Codebook TFs
  land there, the two tasks share a fix. Likely: the Codebook's TFs are
  *poorly characterised by selection*, which is precisely the population whose
  UniProt feature tables are thinnest — so expect the intersection with `none`
  and with the `Others` family to be **enriched, not average**. A newly ingested
  Codebook TF with a motif but no extractable span is designable by `target.py`
  and **not buildable** by `fuse.py`; report that count separately, because it is
  a different kind of gap from the one being closed.
- **287 review-flagged records** need curation; ELM2- and SET-bearing chassis
  should probably be excluded from opener designs outright.
- `dbd_length_aa` is still not a design-time constraint (ZNF850 → a ~980-aa
  fusion, ~2.9 kb against ~4.7 kb of AAV cargo). Independent of this work.

## 11. Also worth fixing while in this code

`zinc_finger_count` counts every `ZN_FING` feature, PHD and RING fingers
included, and `af3.py` adds one Zn ion per count — so a PHD-bearing chassis gets
too many ions in its AF3 input. The extraction work now has the whitelist
(`ZNF_KEEP_TYPES`) that fixes it. Small, unrelated to the Codebook, surfaced here
so it is not lost.

`aggregation_risk` is still collapsed to binary at `ranking.py:63` (every opener
1.0, every closer 0.0) while `max_junction_score` varies 0.15–0.55 — and it is
now load-bearing on the prune. Also unrelated, also worth not losing.
