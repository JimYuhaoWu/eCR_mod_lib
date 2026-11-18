# DBD extraction + library metadata — fix plan

**Status:** proposed 2026-11-12, **revised 2026-11-17** after an independent
re-verification against the code and the shipped `library/module_library.tsv`.
Several premises in the first draft were wrong and are corrected in place; the
corrections are marked **[rev]**. Nothing here is implemented yet.

Audit finding in one line: the extraction work does what it was built to do
(remove the KRAB self-antagonism hazard), but **its safety flag is inert** — 287
records are flagged, the flag is read in two places and acted on in neither —
**the span-length hazard it does flag is larger than the flag reports** (112 of
the 256 spans over 500 aa, and 80 of the 151 spans covering ≥80% of their
protein, carry `dbd_review_flag = 0`), and **692 DBD records have no span at all
and fuse full-length today**, which no flag-based refusal can catch.

**[rev3] The largest live defect is not in the extraction at all:** `fuse.py`
keys effector domains by name with `setdefault`, so **1,439 of 6,969 ED rows
never reach a design run** (F8a). It is also the cheapest fix here and needs no
rebuild.

Two repos, so two PRs. Order matters — see **Sequencing**, which changed in this
revision because F2 as originally sequenced would have opened the exact hazard
the extractor exists to close.

---

## F1 — make `dbd_review_flag` a refusal, not a note

**Repo:** eCR_predictor · **Files:** `fuse.py` (repo root, **not**
`ecr_predictor/fuse.py`) (+ test)

`fuse.py` prefers `dbd_sequence_aa` and warns on fallback — twice: a
whole-input warning when the column is absent (`fuse.py:111-115`) and a
per-gene warning naming up to 10 genes (`fuse.py:95-100`, `116-120`).

**[rev]** The first draft said `dbd_review_flag` "appears only in
`ecr_predictor/output.py`". That is wrong. It is read at `fuse.py:108` into the
`meta` dict (alongside `dbd_extraction` at `:107`) and then **never consulted**;
`ecr_predictor/output.py:47` carries it as a passthrough column; and
`ecr_predictor/fusion/structure.py:56-62` builds the AF3/FoldX frame from
`dbd_meta` but copies only `tf_family` and `zinc_finger_count`, dropping the
flag. So the defect is "read and ignored", not "never read" — a weaker claim,
and the one to state in the paper.

**[rev]** The first draft named "SET-bearing NSD2/KMT2C" among the 287 flagged
chassis. Both are `dbd_review_flag = 0` (NSD2 `dna_bind` 447-527; KMT2C
`dna_bind` 28-52). Their SET/PHD/PWWP domains fall *outside* the kept span, so
flag 0 is correct extractor behaviour — the error is in the documentation, and
it is inherited verbatim from `DBD_EXTRACTION_PLAN.md:284` (see F4).

Verified flagged (`dbd_review_flag = 1`): ZNF541 and TRERF1 (ELM2, i.e. NuRD
recruiters and so unintended closers), ATF7, REST, YY1, FOS, HINFP, ZHX2, ZHX3,
TRPS1, MGA.

**[rev2] The flag is the wrong thing to gate on by itself.** Found on
re-verification 2026-11-17: **692 of the 3,248 DBD records have no
`dbd_sequence_aa`** — 262 with `dbd_extraction = 'none'` (written with
`dbd_review_flag = 0` at `06:298` and `06:586`) and 430 never attempted
(`dbd_review_flag` NULL). Every one of them already falls through
`fuse.py:95-100` to a **full-length fusion today**, which is the KRAB
self-antagonism the extractor exists to prevent (`06:16-21`). A refusal keyed on
`dbd_review_flag == 1` catches **none** of them, because none is flagged. So F1
as first drafted would not have closed the hazard it was written for.

The refusal therefore needs two independent conditions:

- **Absence of a span** (`dbd_sequence_aa` null) → refuse. This is the one that
  matters: 692 records, live today, flag 0 or NULL.
- `dbd_review_flag == 1` → refuse. 287 records, the over-inclusion hazard.
- `--allow-flagged-dbd` to override the second, since some flags are
  conservative (REST, YY1, FOS, HINFP, where the flagged region is itself
  required for binding). A separate `--allow-fulllength-dbd` for the first, since
  falling back to a full-length TF is a different decision from accepting a span
  with a hazard in it.
- Count both skip classes in the run summary so a design run cannot silently
  lose chassis or silently fuse a whole TF.

**[rev] Non-code consumer to update.** `docs/CODEBOOK_EXPANSION_PLAN.md:83`
uses `dbd_review_flag = 0` to define a 468-record "immediately actionable"
cohort. Making the flag a refusal changes what that cohort means, and F7 will
change its membership again.

**Test:** a flagged record is skipped by default and fused under the override.

---

## F2 — gate on span-to-protein ratio, not on feature coverage

**[rev] Retitled and re-scoped.** The first draft blamed the `bipartite` rule
and proposed reusing `MIN_SPAN_COVERAGE` as the detector. Both were wrong.

**Repo:** eCR_mod_lib · **Files:** `scripts/06_extract_dbd.py`, then re-run
`06` + `04_build_library.py`

`min(start) → max(end)` (`06:301-303`) is applied to **every** extraction kind,
not to `bipartite`. Of the 256 spans over 500 aa: **`znf_array` 203**,
`bipartite` 41, `dna_bind` 9, `domain` 3. Of the 618 over 300 aa: `znf_array`
538, `bipartite` 51, **`dna_bind` 18, `domain` 11** (**[rev2]** the last two were
omitted in the previous revision, which made 538 + 51 read as exhaustive when it
sums to 589, not 618). Median span is 106 aa.

**[rev]** The doc's own example table mis-attributed the worst case —
`ZNF462` (2,442 of 2,506 aa, 97%) is `znf_array`, not `bipartite`:

| gene | rule | span | full | % | flag |
|---|---|---|---|---|---|
| ZNF462 | `znf_array` | 2,442 | 2,506 | 97% | 1 |
| Zfhx3 | `bipartite` | 3,510 | 3,726 | 94% | 1 |
| ZFHX4 | `bipartite` | 2,822 | 3,567 | 79% | 1 |
| MGA | `bipartite` | 2,403 | 3,065 | 78% | 1 |
| **ZNF521** | `znf_array` | **1,271** | **1,311** | **97%** | **0** |
| **ZNF423** | `znf_array` | **1,224** | **1,284** | **95%** | **0** |

**[rev] `MIN_SPAN_COVERAGE = 0.5` does not detect these.** `06:342-343`
computes coverage as *merged keep-feature residues ÷ span length* — feature
**density inside the span**, not span ÷ protein. A 97%-of-protein span densely
packed with C2H2 fingers passes that gate comfortably. Measured consequence:

- **112 of the 256 spans >500 aa carry `dbd_review_flag = 0`.**
- **80 of the 151 spans covering ≥80% of their protein carry flag 0** — ZNF521
  (96.9%), ZNF423 (95.3%), ZNF664 (97.7%), **Zfp239** (100.0%, mouse — **[rev3]**
  earlier drafts wrote "ZFP239"; there is no human record), ZNF146 (96.9%) —
  all `znf_array`, all with `dbd_excluded_domains` null.

So the original fix ("below `MIN_SPAN_COVERAGE` → `dispersed`, write no span")
would have left untouched exactly the records the argument condemns.

- Add a **new** `MAX_SPAN_FRACTION` gate on `span ÷ len(sequence_full_aa)`,
  applied to every extraction kind. Threshold to be chosen from the measured
  distribution, not guessed.
- Above it, set `dbd_extraction = 'dispersed'` and write **no**
  `dbd_start/end/sequence_aa/length_aa`; keep the ratio in a note.
- **[rev2] `dispersed` records must also set `dbd_review_flag = 1`** (or F1 must
  gate on span absence, per F1's `[rev2]`). Without one of the two, a
  `dispersed` record is indistinguishable from the 262 existing `none` records:
  flag 0, no span, straight through to a full-length fusion. The previous
  revision's sequencing argument ("F1 must land first") silently assumed this
  and did not state it.
- **[rev]** No `schema.py` change is needed: `dbd_extraction` is plain `TEXT`
  with no CHECK constraint (`schema.py:58`, migration at `schema.py:114`), so
  there is no vocabulary to extend.
- Effect on the shipped span count is **not predicted here** — the per-accession
  feature cache is not on this machine (`data/raw/uniprot/` holds 571 files, all
  `*.fasta`, zero `*_feat.tsv`), so report the real post-run numbers.
- **[rev2]** If the ratio is reported, note that `coverage` can exceed 1.0, so it
  is not the fraction its name implies. Harmless for a `< 0.5` test, but do not
  reuse it as a percentage. **[rev3]** The mechanism given earlier — the flank
  pull-back — was wrong: `lo`/`hi` always contain `[first, last]`, so that cannot
  do it. The real cause is `hi = min(len(seq), last + FLANK)` at `06:303`,
  clipping the span to the stored sequence when feature coordinates run past it
  (isoform mismatch), while the numerator still sums the un-clipped features.

**[rev] Claim dropped.** The first draft said a near-full-length span "violates
the plan's own rule 3 (*never fabricate boundaries*)". Rule 3 (`06:48-50`) is
about emitting a span where no usable feature exists; `min..max` over two
genuine annotated features is not a fabricated boundary. The objection is that
the output is mislabelled as an extraction, which is the argument to make.

**Open question for the paper:** a dispersed multi-homeodomain or long-array TF
may simply not be a usable eCR chassis. If so that is a *design-time exclusion
rule*, which is a more interesting result than a trimming heuristic.

---

## F2b — **[new]** implausibly short spans are unchecked

**Repo:** eCR_mod_lib · **File:** `scripts/06_extract_dbd.py`

Found during the 2026-11-17 re-verification; the first draft audited only the
over-long end. **27 spans are under 40 aa and 23 of those carry
`dbd_review_flag = 0`**; 8 are under 30 aa. `KMT2C` extracts **28-52 = 25 aa
from a 4,911-aa protein** (`dna_bind`, flag 0), `AHCTF1` 25 aa of 2,266,
`TPRX1` 30 aa of 411. The length-sanity gate at `06:443-444` filters to
`znf_array`/`bipartite` only, so `dna_bind` and `domain` spans get **no** length
check at all — and `fuse.py` will happily build a construct around a 25-aa
fragment.

- Extend the length-sanity gate to all four extraction kinds.
- Minimum plausible span per family (a homeodomain is ~60 aa, a C2H2 finger
  ~28 aa); below it, flag rather than silently emit.
- **[rev2] The threshold shape matters as much as the coverage.** Gate 1's test
  is `dbd_length_aa < 20 * kept` (`06:453`), so a 37-aa `znf_array` span with one
  kept finger passes — which is why 16 of the 27 under-40-aa spans are
  `znf_array` and unflagged despite being inside the gate's current scope.
  Extending the gate to all four kinds is necessary but not sufficient.
- **[rev3] Two precisions.** Of the 27 short spans, 16 are `znf_array` and 11 are
  `dna_bind`; **none is `domain`**, so "all four kinds" is a precaution rather
  than a measured gap for that kind. And gate 1 never writes
  `dbd_review_flag` — it logs and counts failures — so "unflagged" above means
  the stored flag is 0, which is a separate thing from the gate's report.

---

## F3 — zinc ion count over-counts fingers dropped by the extraction

**Repo:** eCR_mod_lib (**[rev]** moved — see below) · **Consumer:**
`eCR_predictor/ecr_predictor/af3.py:162`

`zinc_count = int(pd.to_numeric(row.get("zinc_finger_count"), ...) or 0)` adds
one `ZN` ion per count (`af3.py:109-113`), gated to
`ZINC_FINGER_FAMILIES = {"zf-C2H2", "ZBTB"}` (`af3.py:50`, `:161`). The count
itself is `text.count("ZN_FING")` (`01_fetch_dbd.py:194`) — every finger type,
no `/note` filter — computed over the **full protein**, while the fusion path
passes the **span-only** sequence (`fuse.py:95-97` → `fusion/structure.py:59`)
with the full-protein count alongside it (`structure.py:61`).

**[rev] The original example was impossible.** The first draft said "a
PHD-bearing chassis gets ions that are not in its DNA-binding apparatus" and
named NSD2. `zinc_finger_count` is populated *only* for `zf-C2H2`/`ZBTB`
(`01_fetch_dbd.py:171`) and `af3.py:161` gates on the same two families —
**zero** records outside them carry a count, and NSD2/KMT2C are `HMG` with a
null count. The real over-count is non-C2H2 fingers *inside* zf-C2H2 members.

**[rev3] Retracting a `[rev2]` retraction.** `[rev2]` claimed only
`ZnF:CCHC FOG-type` and `ZnF:CCHC HIVEP-type` appear, and that UBZ/C2HC "cannot
be checked without the feature cache". **Both wrong.** The dropped fingers are
recorded in the shipped TSV's own `dbd_excluded_domains` — no cache needed —
and parsing the 16 affected records gives **five** types: `ZnF:CCHC FOG-type`
(20 features), `ZnF:CCHC HIVEP-type` (4), `ZnF:CCHC-type` (4),
**`ZnF:UBZ-type` (2 — ZBTB1/Zbtb1)**, **`ZnF:C2HC-type` (2 — BCL11A/Bcl11a)**.
The original "(CCHC/FOG/UBZ/C2HC)" was correct. Asserting unverifiability about
a column this document mines elsewhere was the error, not the original list.

**[rev] Measured incidence** (the first draft gave none): 1,010 records carry a
count; 994 have both a count and a span; **16** demonstrably drop a finger the
count still includes (32 dropped finger features in total). Worst case
`ZFPM2`/`Zfpm2`: **count 8, span 102 aa of 1,151**, five FOG-type CCHC fingers
excluded — 8 Zn ions for a construct that can hold about 3. Also `ZFPM1`
(9, 600 aa) and `HIVEP1`/`HIVEP3` (5 each). So this is a real bug on a narrow
class, not a broad one.

- Count only `ZNF_KEEP_TYPES` (C2H2, NR C4, GATA, THAP, BED), and only fingers
  inside the kept span.
- Route: have `06_extract_dbd.py` write a `dbd_zinc_finger_count` column — the
  number it already re-derives for the length-sanity gate — and have `af3.py`
  prefer it, falling back to `zinc_finger_count` with a warning.
- **[rev3] That new column needs the same schema work F7 spells out**, which
  this fix omitted: a `schema.py` DDL line **and** a `MIGRATIONS` entry, plus the
  **hardcoded UPDATE column list at `06:636-639`**, which neither F3 nor F7
  mentions. Any new `dbd_*` column needs all three or it will not reach a
  deployed library.
- **[rev] This is therefore a mod_lib change requiring a rebuild**, so it moves
  into sequencing step 3, not the predictor PR. The first draft filed it under
  eCR_predictor and scheduled it *after* the only rebuild, which would not work.
- **[rev] Scope of affected structures narrows.** Only fusion-path structures
  are wrong. `refine.py:150` calls the same `run_af3_prediction` but passes the
  **full-length** TF (`output.py:45`), where the full-protein ion count was
  correct.

---

## F4 — stale and incorrect numbers in the plan doc

**Repo:** eCR_mod_lib · **File:** `docs/DBD_EXTRACTION_PLAN.md`

1. `:301` says *"The 196 `none` records"* and `:302` *"the 32 flagged records
   need curation"*. Shipped values are **262** and **287**. Both stale figures
   are contradicted earlier in the same file — `:236` records "Review flags went
   32 → 287" and `:188` already tabulates `none` 262.
2. **[rev]** `:284` names NSD2/KMT2C as flagged chassis. Both are flag 0 (F1).
   This is the source of the same error in the first draft of this document.

Fix before anything is quoted into a methods table.

---

## F5 — independent validation (the one that is publishable)

**Repo:** eCR_mod_lib · **New:** `scripts/07_validate_dbd_external.py`,
report in `docs/`

**[rev] Narrowed.** The first draft said "all three current gates are
self-consistency checks" and that the validation "is circular". Only gate 2 is
circular in the strict sense: its 7 `ANCHORS` (`06:413-421`) are commented as
"Hand-checked spans from the fib->beta run (2026-09-03)", and
`DBD_EXTRACTION_PLAN.md:271-272` calls the agreement "byte-identical to the
by-hand extraction used for the fib→β constructs" — so it shows the script
reproduces the hand method, not that either is right. Gates 1 (`06:443-465`) and
3 (`06:498-527`) do test something real — span against annotated features — but
re-read the same cached feature table the extractor used, so none of the three
reaches outside the pipeline's own inputs. That is the reviewer objection.

- Pull **InterPro** (or Pfam) domain boundaries via the InterPro REST API,
  cached per accession like `*_feat.tsv`. **[rev2]** The fetch is **2,803
  distinct `uniprot_id`**, not 2,818 — 2,818 is the record count, and some
  human/mouse orthologue pairs share an accession.
- Report per-record boundary agreement (start/end offset, Jaccard over residues)
  stratified by `dbd_extraction` kind and by family.
- Second, smaller arm: for chassis with a solved protein–DNA complex in the PDB,
  check that the span covers the residues making DNA contacts.
- Deliverable is a **distribution, not a pass/fail gate** — disagreement with
  InterPro is informative, since the two disagree legitimately (InterPro splits
  C2H2 fingers individually; this method keeps the array contiguous by design).

This is also the figure panel for the extraction part of the methods section.

---

## F6 — library metadata errors (independent of extraction)

**Repo:** eCR_mod_lib

1. **ED count wrong in docs.** `README.md:47` and `CLAUDE.md:21` both state ED
   ≈ 7,700; the snapshot has **6,973** — a shortfall of 727, i.e. **−10.4%** of
   the actual count (−9.4% of the quoted one). Total is 10,758. **[rev]** Also:
   README's own table sums to 11,490 against the "~10,800 records" at
   `README.md:42`. **[rev2]** The DBD and CR rows are *not* errors: that table's
   header is "Approx. count" and its entries are `~3,250` and `~540` against
   3,248 and 537, which is what an approximation should look like. Only the ED
   row is wrong, and it is wrong far outside rounding.
2. **582 viral EDs labelled `organism = Homo sapiens`** — the whole Ludwig 2023
   set: HHV 261, vTR 186, CoV 135. Only `ED_VP16_TAD` and `ED_RTA_TAD` carry
   `organism = viral`. The labels live in `source`, not `name`, so a
   viral-effector row cannot be produced from `organism` today. Root cause:
   `scripts/02_seed_ed.py:731` hardcodes `"organism": "Homo sapiens"` for every
   screen-derived ED and the `SCREEN_FILES` entries carry no organism key — so a
   `subtype_override`-style per-source organism override is the right fix.
   **[rev]** Also: `README.md:215` records Ludwig as vTR 195 / CoV 135 / HHV 268
   = 598, which disagrees with the data's 186 / 135 / 261 = 582.
3. **`structurally-resolved` is documented as a validation level and used by zero
   records** (`validation_level`: screen-validated 6,167, predicted 2,603,
   motif-only 1,177, ChIP-validated 811). **[rev]** It appears not only in
   `README.md:56` but in the schema CHECK constraint (`schema.py:71`, inside the
   CHECK at `:66-72`) and the allowed set in `scripts/05_validate.py:51`, so
   dropping it from the README alone would leave the schema accepting a level no
   document describes. **[rev3]** That enumeration was still incomplete — there
   are two more live references, both in the *predictor*:
   **`eCR_predictor/ecr_predictor/score.py:11`**, which maps
   `"structurally-resolved"` to confidence `"high"`, and
   `eCR_predictor/CLAUDE.md:186`. Removing the level without the `score.py` entry
   leaves a dangling mapping. **[rev2] This is therefore not doc-only, and the
   schema half cannot ride the migration machinery.** `MIGRATIONS` only performs
   `ALTER TABLE … ADD COLUMN` (`schema.py:172-174`), SQLite cannot drop or alter
   a CHECK constraint in place, and `CREATE TABLE IF NOT EXISTS` is a no-op on a
   deployed library — so removing the level needs a full table rebuild (new
   table, copy, drop, rename) or it will not take effect on any existing
   database. Simplest honest options: drop it from the README and
   `05_validate.py` now and leave the CHECK alone with a comment, or schedule the
   table rebuild deliberately.

**[rev] Consumer impact of (2).** `organism` is read for DBD species matching
(`query.py:39,46`) and ED organism is unused, so predictor risk is low. The
visible effect is on the public API: `load(type="ED", organism="Homo sapiens")`
(`scripts/__init__.py:21-25`) **today returns 582 rows too many**, and would
return 582 fewer once the labels are corrected. **[rev3]** An earlier draft
stated this as "silently returns 582 fewer rows", which reads as current
behaviour and inverts it. That draft also cited `target.py:76` for
`_SPECIES_BY_LABEL` — it is at **`ecr_predictor/target.py:77`** (the root
`target.py` is a 6-line shim), it is used at `:403`, and it drives **DBD**
species selection only. EDs are never species-matched, so the claim that the new
organism value "must be chosen deliberately" because of that map was overstated;
choose it deliberately for the API's sake, not for `target.py`'s.

(1) is doc-only. (2) and (3) change data or schema and ride with the F2 rebuild.

---

## F8 — **[rev2 new]** the same ED sequence appears as both activator and repressor

**Repo:** eCR_mod_lib · **Files:** `scripts/02_seed_ed.py`, `README.md`

Found on re-verification 2026-11-17; the ED metadata is contradictory, and the
collision it exposes (F8a below) is the most consequential live defect found in
this audit. Earlier revisions missed both entirely.

**Mukund 2023 was ingested twice, with opposite labels.** `source` counts give
`Mukund_2023_activators` 52 and `Mukund_2023_repressors` 52 — **104 rows for 52
unique names and 52 unique sequences**. Every domain is in both pools:

| module_id | name | subtype | chromatin_state_effect | score |
|---|---|---|---|---|
| `…_ACTIVATORS_LARP7_HUMAN_RRM_1` | LARP7_HUMAN_RRM_1 | activator | open / H3K27ac deposition | +0.593 |
| `…_REPRESSORS_LARP7_HUMAN_RRM_1` | LARP7_HUMAN_RRM_1 | repressor | close / H3K9me2-3 or H3K27me3 | −0.268 |

The cause is structural, not a typo: the dual-readout screen was imported
through two `SCREEN_FILES` entries with `subtype_override` set both ways
(`02_seed_ed.py:496` and `:519`; the generic application site is `:709-710`), so each domain's activation score and its repression
score became two records with contradictory `subtype`.

**It is not confined to Mukund.** Library-wide: **577 ED rows share a sequence
with another ED row** (272 distinct sequences), and **68 ED sequences carry both
an `activator` and a `repressor` label**.

**[rev3] The consequence an earlier draft gave for this was invented.** It
claimed "`fuse.py` picks an opener or a closer by `subtype`, so the opener and
closer pools are not disjoint". **`fuse.py` never reads `subtype`** — a grep
returns no match — and `chromatin_state_effect` is read nowhere in the
predictor. `_load_eds` (`fuse.py:150-162`) loads every ED, optionally filters by
an explicit `--ed` **name** list, and `assemble_fusions` (`:193-195`) builds the
full DBD × linker × ED cross product. **There are no opener and closer ED pools
in the predictor at all.** ("opener"/"closer" in `preview.py` refers to OCR
regions, not effectors.)

The data defect is real and the measurements above stand. The real consequence
is the one below, and it is larger.

### F8a — the name collision silently drops 1,439 ED rows

`_load_eds` builds `eds` as `{name: sequence}` with **`setdefault`**
(`fuse.py:159-161`), so the first row reached for a name wins and the rest are
discarded. The library has **6,969 sequence-bearing ED rows under 5,530 unique
names**, and **1,148 names map to more than one distinct sequence** — so
**1,439 rows never reach a design run**, and which sequence survives depends on
pandas row order. The Mukund pairs are one instance of this; the collision class
is **25× larger** (1,439 rows vs 104).

This is the defect worth fixing first, and it is a predictor-side fix
independent of any rebuild.

- Make the key unique (`module_id`, or `name` + `source`), or fail loudly on a
  collision instead of `setdefault`. Report how many rows collide.
- Then decide the Mukund representation: one record per sequence carrying *both*
  scores (an `activation_score` and a `repression_score`), or two records with an
  explicit link and a `readout` field. The second is closer to the current shape;
  the first is what a dual-readout assay actually measured.
- **[rev3] The representation fix does not cover everything.** Only 52 of the 68
  both-labelled sequences are Mukund. The other 16 are **cross-source**
  activator/repressor contradictions — `DelRosso_2023` (8 repressor / 7
  activator), `HiTEff` (4/4), `Ludwig_2023_vTR` (4), `Ludwig_2023_HHV` (4),
  `Tycko_2025_ADs` (2), `TENet_2024_RDs` (1) — not a dual-readout artefact, so
  they need a separate adjudication rule.
- Also wrong in the docs: `README.md:101` and `:218` claim "55 unique entries" /
  "55 non-HiTEff short nuclear domains" against 52 unique / 104 rows, and
  **[rev3]** `CLAUDE.md`'s "Pending work" section carries a third copy
  ("55 activator + 55 repressor entries (110 total)").

### F8b — a Staller 2022 record is a parse failure

`ED_STALLER_2022_WT_DOMAINS_NONE` has `name` literally `'none'`, `subtype`
activator, 39 aa — and it is the **only** `Staller_2022_wt_domains` row, where
`README.md:213` claims "+ 5 WT domains". So 1 of 5 was ingested and the one that
survived has a garbage name. The 50 `Staller_2022_predictions` rows are fine
(though their names retain raw FASTA headers, `>sp|P01100|FOS_HUMAN_208-246`).

---

## F9 — **[rev2 new]** `dbd_extraction = 'none'` conflates two different causes

**Repo:** eCR_mod_lib · **File:** `scripts/06_extract_dbd.py`

At `06:583-587`, when `features.get(acc)` is `None` — i.e. UniProt returned
nothing for that accession — the record is written as `'none'` with
`dbd_review_flag = 0` and no excluded domains, which is **indistinguishable from
a genuine rule-3 abstention** ("entry has features, none of them usable").

The script's own fetch logic argues at length (`06:197-206`) that a failed fetch
must never masquerade as "no features", and then this path does exactly that.
So the 262 `none` records are a mixture of two causes with no column separating
them, and F4's and F5's discussion of the `none` cohort treats it as one thing.

- Separate the verdicts: `'none'` for a real abstention, `'no_features'` (or a
  `dbd_extraction_reason` column) when the fetch returned nothing.
- Report the split. If a material share of the 262 is actually fetch failure,
  those accessions are retryable and the `none` count is overstated.

---

## F7 — record what was dropped, even when UniProt does not name it

**Repo:** eCR_mod_lib · **File:** `scripts/06_extract_dbd.py`
**Status: code landed 2026-11-18** (`scripts/06_extract_dbd.py`,
`scripts/schema.py`, and the passthrough in `eCR_predictor/ecr_predictor/output.py`).
The columns reach the library only on the next `06` → `04` pass, step 3 of the
sequencing below. Until then `dbd_dropped_regions` does not exist in the
deployed DB, and the `MIGRATIONS` entry is what adds it there. `06:NNN` line
references in this entry are to the file as it stood **before** the fix; the
additions shifted them.

Found 2026-11-16 while checking Nanog for the manuscript's step 3. `DBD_MU_NANOG`
(Q80Z64) extracts to `dna_bind` **90-161 (72 aa)** of 305, `dbd_review_flag = 0`,
`dbd_excluded_domains = NaN`. The mouse sequence has tryptophans at 143, 154,
159, then **198, 203, 208, 213, 218, 223, 228, 233, 238, 243** — ten at exact
5-residue spacing, i.e. the **tryptophan-repeat (WR) domain**, through which
Nanog dimerizes. It lies entirely outside the kept span, and the pipeline
recorded nothing in `dbd_excluded_domains` and did not flag the record.

**[corrected 2026-11-18] UniProt does name it.** The first draft of this entry
said no feature was annotated there. Fetching Q80Z64's feature table shows
otherwise:

```
REGION 198..247; /note="10 X repeats starting with a Trp in each unit"
REGION 198..247; /note="Sufficient for transactivation activity"
REGION 248..305; /note="Sufficient for strong transactivation activity"
```

So this was not missing knowledge, it was a **discarded annotation** — which
makes the defect worse, and the fix cheaper. Mechanism in code:
`dbd_excluded_domains` is written only from `dropped` (`06:364`), populated
solely from annotated `ZN_FING` and `DOMAIN` features (`06:281-292`). `REGION`
features enter the extraction at exactly one place — a hazard-hint test for
regions lying *inside* the kept span (`06:333-335`) — and are otherwise read
and thrown away. Everything outside the span, named or not, goes unrecorded.
Corroborating: **199 of the 287 flagged records have a null
`dbd_excluded_domains`**, so even a flag's reason is often unrecorded.

**[rev] Two qualifications.** (i) The total unrecorded drop is **233 of 305
residues (76%)**, not 144 — the extraction also drops 1-89, which the first
draft omitted; 162-305 = 144 aa (47%) is only the C-terminal half.
(ii) `DBD_HO_NANOG` is not the identical case: human tryptophans are at 142,
153, 158, 196, 201, 206, 216, 221, 226, 231, 236 — eleven, with gaps of
5, 5, **10**, 5, 5, 5, 5, not a clean every-5 ladder. Verify per record rather
than generalising across the pair. (The human record drops the same 233 aa /
76% in total, since its span is also 72 aa of 305.)

This is a different failure class from the KRAB case the extraction was built
to fix. The existing gates audit **over-inclusion** of named effector domains;
nothing audits **under-inclusion** of function-essential modules that carry no
UniProt feature. F1 would not catch Nanog, because its flag is 0.

- Record every residue range outside the kept span in a new
  `dbd_dropped_regions` column — named features and unannotated gaps alike,
  as `start-end:label` with `label = unannotated` where UniProt is silent.
  Where UniProt gives a residue two annotations they are joined with ` + `, so
  Nanog's 198-247 reads as both the Trp repeat and a transactivation region
  rather than one hiding the other.
- Flag in a **separate `dbd_dropped_review_flag` column**, never in
  `dbd_review_flag`. **[decided 2026-11-17]** See the F7/F1 resolution under
  Sequencing.
- **[corrected 2026-11-18] Drop size cannot be the detector.** This entry said
  to flag "when an unannotated dropped region is large (threshold to pick from
  the distribution)". The distribution rules that out. Over the 2,556 DBD
  records with a span and a sequence, the dropped fraction runs p25 **45.3%**,
  p50 **70.5%**, p75 **80.1%** — Nanog's 76.4% sits between the median and p75,
  so it is a *typical* record, not an outlier. A threshold low enough to catch
  it flags 1,306 records at ≥70% and 1,827 at ≥50% (898 and 1,266 of them with
  no `dbd_excluded_domains`). This is not a tuning problem: a TF's DBD is
  *supposed* to be a small fraction of the protein, so drop size measures "is
  this a normal transcription factor", not "was something essential lost".
  Size is still recorded, as the span and the regions; it just cannot gate.
- **What the flag is instead.** Two tests, both on dropped regions only:
  (i) the region's annotation reads as **self-association** —
  `SELF_ASSOCIATION_HINTS` = dimeriz / oligomeriz / trimeriz / tetrameriz /
  self-associ / coiled coil — or as a **tandem repeat ≥20 aa**
  (`MIN_REPEAT_REGION_AA`, which excludes compositional annotations such as
  MESP1's 4-residue "2 X 2 AA tandem repeats of G-Q"); or
  (ii) for a stretch UniProt does **not** name, a **repeat screen** finds a
  periodic single-residue ladder — ≥6 rungs at a constant spacing k ∈ 3..8,
  tolerating one skipped rung (gap 2k, which is what human Nanog's W ladder
  has: 5,5,**10**,5,5,5,5), with the residue capped at 8% of the segment so
  poly-Q and poly-G composition does not swamp it.
  Activation and interaction annotations are deliberately **not** flagged:
  removing those is the purpose of the extraction.
  **Measured against live UniProt feature tables for all 2,818 extractable DBD
  records: 97 flagged** — 63 by a self-association annotation (55 of those
  naming dimerization), 31 by a named tandem repeat ≥20 aa (the Nanog family
  among them), and 3 by the unnamed-ladder screen alone. Reviewable by hand,
  which a 1,300-record size threshold would not be.
- **A second case the same flag catches: p53.** `DBD_HO_TP53` extracts 96-298
  of 393 and drops **325-356, annotated `Oligomerization`** — the
  tetramerization domain p53 binds DNA through. Nanog is not an isolated
  example of a dropped assembly module.
- No keep-rule change: whether WR belongs in the construct is a biology
  question, and step 3 of the manuscript is designed to answer it. The fix is
  that the pipeline must **say** it dropped it.
- **[rev] Implementation surface is wider than one file.** A new column needs a
  `schema.py` DDL line **and** a `MIGRATIONS` entry (`schema.py:108-117`, whose
  comment at `:104-107` explains why omitting it leaves deployed libraries
  silently without the column), plus `ecr_predictor/output.py:45-47` if the
  predictor is to see it. **[decided 2026-11-17]** With the flag in its own
  column this is **two** new columns — `dbd_dropped_regions` and
  `dbd_dropped_review_flag` — so two DDL lines and two `MIGRATIONS` entries.

**Why this matters beyond the data.** Nanog is a better motivating example for
the methods section than the KRAB case: the KRAB failure was at least labelled,
this one was silent. It also gives the manuscript a concrete prediction — a
72-aa Nanog chassis lacks the dimerization module, so DBD+ED may underperform
for a recorded reason rather than an unexplained one.

---

## Sequencing

**[rev2] Revised twice.** The first order had three defects; the second had one
more, noted in step 2.

1. **F4, F6.1 — doc-only, no rebuild.** (F6.3 moved out of this step: its schema
   and `05_validate.py` halves are code, and the CHECK cannot be migrated — see
   F6.3's `[rev2]`.)
2. **F1 — predictor PR, *before* the rebuild, and gating on span absence as well
   as the flag.** The original order put F2 first: F2 nulls `dbd_sequence_aa` for
   `dispersed` records, and until F1 lands those records fall straight through
   `fuse.py:95-100` to a **full-length fusion**, reintroducing the KRAB
   self-antagonism the extractor exists to prevent (`06:16-21`).
   **[rev2]** But the flag-only version of F1 would not have helped, because
   nothing says `dispersed` records get flagged — so they would have behaved
   exactly like the 262 existing `none` records and fused full-length anyway.
   F1 must gate on **span absence**; that also closes the 692-record hazard that
   is live today, independently of F2.
3. **F2 + F2b + F3 + F6.2 + F6.3 (code half) + F7 + F9 — one pass of `06` → `04`.**
   F3 moves here from the predictor PR because its preferred route writes a new
   mod_lib column; the `af3.py` side follows in step 4. **[rev2]** Do **not**
   pass `04 --rebuild`: that unlinks the DB (`04:63-65`) and would discard the
   extraction. Plain `04` is safe, because `insert_or_replace` writes only the
   columns a record mentions (`schema.py:209-220`). The TSV is exported by `04`
   itself (`04:108-109`), so no separate re-export step is needed.
4. **F3 consumer side (`af3.py` prefers `dbd_zinc_finger_count`) — predictor PR,
   against the rebuilt library.**
5. **[rev3] F8a — the `_load_eds` name collision — predictor PR, independent of
   everything above and of any rebuild.** 1,439 ED rows never reach a design run
   today; this is the largest live defect in the document and the cheapest to
   fix. It can go first if convenient.
6. **F8 representation + F8b — decide before implementing.** A data-model
   question (one record with two scores, or two linked records with a `readout`
   field), not a patch: it changes the 104 Mukund rows, the 16 cross-source
   both-labelled sequences, and the one Staller parse failure. Runs in parallel
   with steps 2–4 since it touches `02_seed_ed.py`, not the DBD path.
7. **F5 — after the data is final**, so the concordance numbers are the ones the
   paper reports.

**[decided 2026-11-17] Interaction between F7 and F1 — resolved.** F7 flags
dropped self-association and repeat modules. *If* it reused `dbd_review_flag`,
`DBD_MU_NANOG` would become flag 1 — its dropped 198-247 is an annotated
tryptophan repeat — and F1 would skip it by default. The manuscript's step-3
Nanog chassis would disappear from design runs unless `--allow-flagged-dbd` were
passed, directly against F7's "no keep-rule change". (**[corrected 2026-11-18]**
This clause previously justified the flag by the drop being "233 aa / 76% of the
protein, past any plausible threshold". That reasoning was wrong — 76% is near
the median drop — and is replaced by the detector described under F7. The
decision it supports is unaffected: the flag is in its own column either way.)

**Decision: F7 writes a separate `dbd_dropped_review_flag`, and F1 does not gate
on it.** Rationale: over-inclusion (a kept span containing an effector domain)
and under-inclusion (a dropped region containing a function-essential module)
are different hazards with different remedies, and a single flag cannot express
which one fired. `dbd_review_flag` keeps its current meaning, so F1's refusal
logic and the `--allow-flagged-dbd` override are unaffected, and the 468-record
`dbd_review_flag = 0` cohort at `docs/CODEBOOK_EXPANSION_PLAN.md:83` does not
change membership from F7 — only from F1's re-reading of the flag.
`DBD_MU_NANOG` therefore stays in design runs with its 233-aa drop recorded,
which is what step 3 needs: the chassis is testable *and* the pipeline says what
it removed.

**[rev2]** F2 would *not* also catch Nanog — its span is 23.6% of the protein,
far below any `MAX_SPAN_FRACTION` — so F7 is the only fix that changes Nanog's
status, and the interaction above is the only one to resolve for it.

## Not in scope

No functional validation that a trimmed DBD still binds. Flank = 6 aa stays an
assumption. The manuscript should scope the claim as *removes a known
self-antagonism hazard*, not *produces functional minimal DBDs*.

## Environment caveat

The local `library/module_library.db` that `eCR_predictor/query.py:12` reads by
default holds 38 DBD rows and **no `dbd_*` columns**, so nothing in this plan
can be cross-validated against what the predictor reads on a dev machine. The
authoritative library is the server rebuild (`CLAUDE.md`). All counts in this
document come from the committed `library/module_library.tsv` (10,758 records;
note the file has 10,798 *lines* — 39 fields contain embedded newlines, in
`target_or_mechanism`, `linker_notes` and `engineering_compatibility`).

**[rev3] Correcting `[rev2]`, which asserted a local outcome without running it.**
`[rev2]` said gates 1 and 3 "pass locally by vacuity" and that a local "all gates
pass" is meaningless. Actually running `python scripts/06_extract_dbd.py
--validate` here gives: gate 1 `0 ZnF records checked`, gate 2 **`0/7 match`**
(all seven anchors "not in library"), gate 3 `0 records`, and
**`[validate] 1 gate(s) failed`** — a local run *fails*, so the state `[rev2]`
warned about cannot arise. The mechanism was also wrong: gates 1 and 3 report
zero because `done` is empty (the local DB has no `dbd_extraction` column at
all), not because of the cache check at `06:448-449` / `06:503-504`.

What remains true at the code level: given a database that *does* carry spans
but no `*_feat.tsv`, those two gates would skip every record. Run validation on
the server, where both the library and the cache are authoritative.

The local DB is also missing `zinc_finger_count`, so **F3's consumer behaviour
is untestable here too** — not only the `dbd_*` columns.
