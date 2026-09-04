# Next modification: expand the DBD library with Codebook motifs

**Status:** proposed, 2026-09-04. Follows the completed DBD extraction
([DBD_EXTRACTION_PLAN.md](DBD_EXTRACTION_PLAN.md)).

**Source:** *An expanded codebook of human transcription factor DNA-binding
specificity* (Codebook consortium). **Verify before implementing** — the exact
title, DOI, assay set, TF count and supplementary file layout must be read off
the paper itself, not from this document. Everything below that depends on those
numbers is marked as an assumption.

## The gap this fills

Measured on the current library:

| | count |
|---|---|
| human `type='DBD'` rows | 1,637 |
| with a `jaspar_id` | 964 |
| **without any motif** | **673** |
| of those, with an extracted DBD span | 553 |

A DBD with no motif is not merely under-annotated — it is **undesignable**.
`target.py` drops motif-less DBDs outright, because a DBD that cannot be placed
on a locus by sequence cannot be assigned a binding site; `offtarget.py` scores
burden per motif and joins on `jaspar_id`. So those 673 records are inert in the
design path: they cost nothing and contribute nothing.

That is the direct measurement of what the Codebook adds. Its stated purpose is
motifs for human TFs that previously had none, which is exactly this set. The
553 with an extracted span are the immediately actionable subset: motif plus a
correct DBD boundary is everything `target.py` → `cli.py` → `fuse.py` needs.

**Assumption to check first:** that the Codebook's TFs actually intersect our
673. They may already be in AnimalTFDB (then this is pure enrichment), or they
may be TFs we do not carry at all (then it is enrichment *and* new records). The
first task below decides which, and the answer changes the size of the work.

## Step 0 — measure the overlap before writing any ingest code

Pull the supplementary motif table, map its identifiers to UniProt accessions,
and intersect with the library. Report four numbers:

1. Codebook TFs already in the library **with** a motif — do we agree, or does
   the Codebook motif differ from the JASPAR one? A disagreement is a finding,
   not a merge conflict to paper over.
2. Codebook TFs in the library **without** a motif — the target set.
3. Codebook TFs **not in the library** — new DBD records.
4. Codebook TFs we cannot map to an accession — reported, not guessed.

This is a half-day of work and it determines everything after it. Do not skip to
the ingest.

## The architectural decision: `jaspar_id` is the join key

This is the part that needs a decision before code, because it touches both repos.

Today, a motif is a `jaspar_id`, and `scan.py` resolves it in a fixed order:
local `jaspar_cache/<id>.jaspar` → BioPython `JASPAR2020` → JASPAR REST.
`offtarget.py` and `target.py` both join on `jaspar_id`. A Codebook motif has no
JASPAR accession unless and until it is deposited there.

Three options:

**(a) Wait for JASPAR deposition.** Zero code, zero risk, and the motifs arrive
through the existing path on the next JASPAR release. Costs an unknown delay,
and partial deposition means partial coverage with no way to tell which is which.

**(b) Add `motif_id` + `motif_source`, keep `jaspar_id` as-is.** New columns;
`scan.fetch_pwm` gains a resolution branch for `motif_source='codebook'` reading
a local PWM cache. `jaspar_id` keeps its current meaning exactly, so nothing that
reads it changes behaviour. Both repos need a coordinated release, and the
run-bundle/site contracts need the new column threaded through.

**(c) Write Codebook motifs into `jaspar_id` with a prefix.** Cheapest to
implement and worst to live with: it makes `jaspar_id` mean two things, which is
the same class of mistake as `sequence_aa` holding whole proteins — the bug the
last cycle was spent undoing. Rejected, and recorded here so it is not
rediscovered as a shortcut.

**Recommendation: (b)**, with (a) as the fallback if the Codebook motifs turn out
to be JASPAR-deposited already — which Step 0 will reveal. (b) is more work but
it is the only option that keeps provenance honest: a design built on a Codebook
motif should be traceable to a Codebook motif.

## Ingest

Follows the established `02_seed_ed.py` pattern — a declarative source entry
plus supplementary files in `data/manual/`, a `literature/papers.yaml` record,
and a provenance row. Concretely:

1. `data/manual/codebook_<year>_<file>` — the supplementary motif and metadata
   tables, checksummed into `download_manifest.json` like every other source.
2. `scripts/07_fetch_codebook.py` — parse motifs to PFM, write
   `data/raw/codebook/<motif_id>.pfm`, emit `data/processed/codebook_raw.tsv`.
   Separate script, not a branch inside `01_fetch_dbd.py`: the sources have
   different shapes and different failure modes, and 01 is already slow.
3. `04_build_library.py` gains `codebook_raw.tsv` as a source for the new/updated
   DBD rows.
4. `validation_level` — the Codebook's assays are direct binding measurements, so
   these are at least `motif-only`, and TFs with orthogonal assay support are
   arguably `screen-validated`. **Decide this explicitly against the paper's own
   confidence tiers**; do not map every row to one level for convenience. This
   column feeds `annotation_confidence`, which is one of the two independent
   scores the predictor reports.

## Validation gates

Same principle as the extraction work: gates that can actually fail.

1. **Round-trip** — every ingested PFM re-reads to the same matrix, sums to a
   consistent count per column, and has the length the source claims.
2. **Known-motif agreement** — for Codebook TFs that *also* have a JASPAR motif,
   compare (Pearson on aligned PFM columns, or Tomtom). High agreement validates
   the ingest; disagreement is a real result about assay differences and should
   be reported per TF rather than averaged away.
3. **Designability delta** — re-run `target.py` on an existing run bundle and
   report how many previously undesignable regions gain a candidate DBD. This is
   the number that says whether the work paid off, and it is the honest headline:
   if the gap report barely shrinks, the motifs are not reaching the loci we
   nominate.
4. **No silent widening of off-target burden** — new motifs mean new genome-wide
   scans in `offtarget.py`. Confirm the burden distribution of existing DBDs is
   unchanged (it must be — burden is per motif), and budget the scan time for
   the new ones.

## Cost

Ingest and PFM handling: local, minutes. The real cost is `offtarget.py`, which
scans genome-wide per motif — at the current sharded rate that is roughly
**45 minutes per 32-way shard per motif set**, so a few hundred new motifs is an
overnight job, not an interactive one. Budget it as a batch run.

## Sequencing against the open extraction work

The extraction cycle left three things open, and two of them interact with this:

- The **196 `none` records** and the `Others` family need per-family keep rules.
  If Codebook TFs land there, the two tasks share a fix.
- The **287 review-flagged records** need curation; ELM2- and SET-bearing
  chassis should probably be excluded from opener designs outright.
- `dbd_length_aa` is not yet a design-time constraint (ZNF850 → a ~980-aa
  fusion, ~2.9 kb against a ~4.7 kb AAV cargo). Independent of this work.

A newly ingested Codebook TF with a motif but no extractable DBD span is
designable by `target.py` and **not** buildable by `fuse.py` — worth reporting as
its own count, since it is a different kind of gap from the one being closed.

## Also worth fixing while in this code

`zinc_finger_count` counts every `ZN_FING` feature, PHD and RING fingers
included, and `af3.py` adds one Zn ion per count — so a PHD-bearing chassis gets
too many ions in its AF3 input. The extraction work now has the whitelist
(`ZNF_KEEP_TYPES`) that would fix it. Small, unrelated to the Codebook, and
surfaced here so it is not lost.
