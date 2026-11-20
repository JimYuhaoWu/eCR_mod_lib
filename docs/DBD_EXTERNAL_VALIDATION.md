# External validation of the extracted DBD spans

**F5 in [`DBD_EXTRACTION_FIXES.md`](DBD_EXTRACTION_FIXES.md) · run 2026-11-18 ·
script [`scripts/07_validate_dbd_external.py`](../scripts/07_validate_dbd_external.py)**

The three acceptance gates in `06_extract_dbd.py --validate` all re-read the
same cached UniProt feature table the extraction used, and gate 2's anchors are
the by-hand spans the script was written to reproduce. So they show the method
is self-consistent, not that it is right. This report answers that with two
arms that reach outside the pipeline's inputs.

**It is a distribution, not a pass/fail gate.** Pfam and this extraction
disagree by design — Pfam matches each C2H2 finger separately, while the
extraction keeps a finger array contiguous because the inter-finger linkers are
part of the array — so *precision* below 1 is expected and *recall* is the
number that has to be high.

---

## Headline

| | result |
|---|---|
| Records compared against Pfam | **2,816** of 2,818 extractable (2 accessions failed to fetch) |
| Recall of the best-matching Pfam DNA-binding domain | **median 1.000**, p25 1.000, min 0.488 |
| Spans overlapping **no** Pfam DNA-binding domain | **20** — and they are 5 gene families, not noise |
| Records checked against DNA contacts in a solved complex | **220** across **207** structures |
| Of those, spans containing **every** DNA-contacting residue | **185 (84%)** |
| Abstentions (`dbd_extraction='none'`) where Pfam names a DBD | **28** of 262 — recoverable chassis |

**The two arms independently convict the same two gene families.** Pfam says
EBF1–4 and RBPJ/RBPJL have their DNA-binding domain somewhere the span does not
cover; the crystal structures of Ebf1 (3MLN) and Rbpj (3IAG) put **0 of 30** and
**0 of 16** DNA-contacting residues inside the span. Different evidence, same
verdict, on records that carry `dbd_review_flag = 0`.

---

## Arm 1 — Pfam domain boundaries via InterPro

### Method

Pfam locations per accession from the InterPro REST API, cached in
`data/raw/interpro/{acc}.json` (2,801 of 2,803 accessions; 2 failed).

Pfam is HMM-based. The independence claim is measurable, so it was measured:
of the **1,156 `DNA_BIND` features** in the cached UniProt feature tables,
**1,040 (90%) carry PROSITE-ProRule evidence** (`ECO:0000255|PROSITE-ProRule`),
29 are by-similarity, 4 cite a PubMed experiment and 50 carry no evidence
string. So the annotation this extraction follows is overwhelmingly one
profile-based method, and Pfam is a different one on the same sequence — not an
independent experiment, but not a copy of the same curation either.

All **513** distinct Pfam family names the fetch returned were classified by
hand into **98 DNA-binding families (2,575 hits)** and **419 others (2,038
hits)**. The DNA-binding side is led by Zinc finger C2H2 (788), Homeodomain
(335), bHLH (162), the nuclear-receptor C4 double treble clef (86), Forkhead
(82) and HMG box (79); the other side by KRAB box (356), nuclear-receptor
ligand-binding (90), BTB/POZ (67), SCAN (60), PAS (40) and OAR (27) — which are
exactly the modules the extraction is built to drop.

Recall is scored against the **single** Pfam DNA-binding domain the span
overlaps most, not their union: ZSCAN20 carries a Myb/SANT-like domain hundreds
of residues from its C2H2 array, and scoring against the union would mark a
correct span as having missed a module it was never meant to contain. Omitted
DNA-binding domains are counted separately instead.

### Verdicts

| verdict | n |
|---|---|
| span overlaps a Pfam DNA-binding domain | 2,302 |
| no span (`none` / `no_features` / `dispersed`) | 379 |
| span exists, Pfam names no DNA-binding family anywhere in the protein | 115 |
| **span overlaps no Pfam DNA-binding domain** | **20** |

### Agreement where the span overlaps a domain (n = 2,302)

| metric | median | p25 | p75 | min |
|---|---|---|---|---|
| recall (domain residues inside the span) | **1.000** | 1.000 | 1.000 | 0.488 |
| precision (span residues in that domain) | 0.770 | 0.565 | 0.802 | 0.021 |
| Jaccard | 0.768 | 0.565 | 0.800 | 0.021 |
| start offset (span − domain) | −7 aa | −11 | −6 | −873 |
| end offset (span − domain) | +7.5 aa | +6 | +11 | −51 |

The ±7 aa halo is `FLANK = 6` plus the boundary pull-back, visible as expected.

| kind | n | recall median | recall ≥0.95 | precision median |
|---|---|---|---|---|
| `dna_bind` | 756 | 1.000 | **99%** | 0.792 |
| `znf_array` | 899 | 1.000 | **99%** | 0.613 |
| `bipartite` | 211 | 1.000 | **98%** | 0.697 |
| `domain` | 436 | 1.000 | **93%** | 0.800 |

`znf_array`'s low precision is the predicted legitimate disagreement — the
array is contiguous by design and Pfam matches only the fingers.
`domain` is the weakest kind on recall, which is the kind whose spans rest
entirely on `KEEP_DOMAIN_HINTS` rather than on an annotated `DNA_BIND` feature.

By family, median recall is 1.000 for every one of the 15 largest (zf-C2H2 803,
Homeobox 363, bHLH 204, TF_bZIP 100, Fork_head 94, HMG 92, ZBTB 74, ETS 52, MYB
50, THR-like 38, T-box 32, Pou 30, zf-GATA 25, RXR-like 24, E2F 21), except
zf-GATA at 0.956.

### The 20 spans that overlap no Pfam DNA-binding domain

| records | span is | Pfam says the DBD is |
|---|---|---|
| EBF1, EBF2, EBF3, EBF4 (+ mouse) — 8 | the IPT/TIG domain | Transcription factor COE1 DNA-binding domain |
| RBPJ, RBPJL (+ mouse) — 4 | the IPT/TIG domain | Beta-trefoil DNA-binding domain / LAG1 DNA binding |
| ZNF280C, ZNF280D (+ mouse) — 4 | a C2H2 stretch | Z280C/D-like C2H2 zinc finger |
| ZNF618 (+ mouse) — 2 | a C2H2 stretch | Hermes transposase DNA-binding domain |
| RAG1 (+ mouse) — 2 | 386-465 | Recombination-activation protein 1 DNA-binding domain |

**18 of the 20 carry `dbd_review_flag = 0`.** The EBF and RBPJ cases share a
cause: `IPT/TIG` is in `KEEP_DOMAIN_HINTS`, and in these proteins it is the
protein-interaction module, not the DNA-binding one. 12 library records rest on
an `IPT/TIG` match (of 491 whose span rests on the hint list alone), so that
hint is the single highest-value thing in this report to revisit.

### The 28 abstentions Pfam disagrees with

Of the 262 `dbd_extraction='none'` records, **28 have a Pfam DNA-binding
domain** — chassis the method is losing rather than correctly declining:

| Pfam family | n | records |
|---|---|---|
| STAT protein, DNA binding domain | 14 | STAT1–6, STAT5A/5B, human + mouse |
| Zinc-finger of C2H2 type | 5 | ZNF385A/B, Zfp385a, Zmat1, Zmat3 |
| GC-rich sequence DNA-binding factor-like | 4 | GCFC2, PAXBP1, + mouse |
| HSF-type DNA-binding | 2 | HSF1, Hsf1 |
| P53 DNA-binding domain | 2 | TP73, Trp73 |
| CIZ1-like, C2H2-type zinc finger | 1 | CIZ1 |

The STAT family is the striking one: 14 records, a major TF family, abstained on
because UniProt annotates no feature the keep rules recognise. TP73 is the same
story next door to TP53, which extracts fine.

The remaining 111 of the 139 spanless-but-Pfam-positive records are `dispersed`
— F2 removed their spans on purpose — and the Pfam hits confirm what F2 said:
RREB1, ZNF236, TRPS1, ZHX3, MGA and HIVEP3 do carry real DNA-binding domains,
spread across proteins too large to use as chassis.

---

## Arm 2 — DNA contacts in solved complexes

### Method

For each accession, UniProt's own PDB cross-references are read in batches, one
RCSB GraphQL query per 200 entries returns each entry's DNA chain count and
resolution, and the best-resolution protein–DNA entry per accession is selected
— deterministic, so a capped run is reproducible. For the top 220 by resolution
the mmCIF is downloaded, every protein residue with a heavy atom within
**4.0 Å** of a DNA heavy atom is collected, and those residues are mapped to
UniProt numbering through SIFTS. Both author and label numbering are handled:
SIFTS leaves the author fields null for some entries (1DSZ/RARA), and requiring
them silently dropped three quarters of the structures in the first run.

### Results

| | n |
|---|---|
| records measured | **220** (207 structures) |
| span contains **every** DNA-contacting residue | **185 (84%)** |
| partial | 30 |
| zero coverage | 5 |

| kind | n | mean coverage | full |
|---|---|---|---|
| `dna_bind` | 97 | 0.986 | 91/97 |
| `bipartite` | 46 | 0.959 | 33/46 |
| `znf_array` | 26 | 0.945 | 21/26 |
| `domain` | 51 | 0.910 | 40/51 |

**The 30 partial records fail at the boundary, not in the middle.** 93 of the
97 missed contacts lie within 30 aa of a span edge; the other 4 are NFATC2's
second contact region at 663-666, 83 residues past a span ending at 580:

| extra flank | missed contacts recovered | partial records completed |
|---|---|---|
| +6 aa | 49/97 | 10/30 |
| +10 aa | 66/97 | 15/30 |
| **+15 aa** | **80/97** | **24/30** |
| +20 aa | 92/97 | 28/30 |

**[corrected 2026-11-18] That table is a numeric sweep and overstates the
gain.** It shifts the recorded boundaries arithmetically, so it applies neither
the excluded-domain pull-back, nor the dispersed gate, nor the clamp to the
sequence ends. Re-running `extract_one` at each value and re-scoring the same
220 records gives:

| FLANK | contacts recovered | full coverage | partial | zero | span lost entirely |
|---|---|---|---|---|---|
| 6 (then) | — | 185 | 30 | 5 | 0 |
| 10 | 33/97 | 191 | 23 | 5 | 1 |
| **15 (now)** | **63/97** | **199** | 15 | 5 | 1 |
| 20 | 78/97 | 208 | 6 | 5 | 1 |

The pull-back accounts for the whole difference: where a longer flank would
reach toward an excluded domain it is shortened again, which is the behaviour
that makes a larger FLANK safe in the first place. (These rows are the flank
alone. The figures for all four changes together are in "What the four changes
did, end to end" below.)

Examples: GATA3 (3DFX) misses 348–365 with a span ending at 347; THRB (2NLL)
misses 188–198 ending at 187; NR5A2 (5L0M) misses 161–173 ending at 160.

**The 5 zero-coverage records split two ways**, and the distinction matters:

| record | structure | reading |
|---|---|---|
| Rbpj | 3IAG — *CSL (RBP-Jk) bound to HES-1 nonconsensus site* | **span is wrong** (IPT/TIG, not the DNA-binding NTD/BTD) |
| Ebf1 | 3MLN — *DNA binding domain of Early B-cell Factor 1 bound to DNA* | **span is wrong** (IPT/TIG) |
| MBD4 | 4OFA — *thymine glycosylase activity on T:O6-methylG mismatch* | span is the MBD; the structure is of the **glycosylase** domain |
| MBD1 | 5W9Q — *zinc finger region of MBD1 in complex with CpG DNA* | span is the MBD; the structure is of the **CXXC fingers** |
| NSD2 | 5VC8 — *WHSC1 PWWP1 domain* | span is the annotated DNA-binding region; the structure is of **PWWP** |

So three of the five are not extraction errors at all: the protein has a second
DNA-contacting module and the available structure happens to be of that one.
Only Rbpj and Ebf1 are genuine — the same two families arm 1 convicts.

---

## What this suggests, as tradeoffs rather than changes

**All four were implemented on 2026-11-18.** Each entry below keeps what the
recommendation said, then what implementing it actually showed — and in two
cases the recommendation's own reasoning turned out to be wrong.

1. **Revisit `IPT/TIG` as a keep hint.** 12 records rest on it, and both
   families where it is testable (EBF, RBPJ) are wrong by two independent
   measures. The cheapest fix is to stop treating `IPT/TIG` as DNA-binding on
   its own; the cost is that any protein where it genuinely participates loses
   its span and falls to `dbd_extraction='none'`.

   **Done: removed from `KEEP_DOMAIN_HINTS`.** It changes **16** records, not
   12, in three different ways, and the stated cost did not materialise.
   - EBF1—4, RBPJ and RBPJL (12) lose their only keep-feature and are then
     picked up by the Pfam fallback of item 2, so they end **extracted on the
     right domain** rather than abstained. Recomputed from the deposited
     coordinates: 3MLN puts **30 of 30** of Ebf1's DNA-contacting residues
     inside the new 3—261 span, against **0 of 30** inside 256—351, and
     3IAG puts **16 of 16** of Rbpj's inside 59—369 against **0 of 16**
     inside 375—477. That is the same coordinate evidence that convicted
     these spans, re-run; the fix does not rest on the Pfam family names that
     motivated it. (Recomputed rather than read off the per-record table, whose
     `missed` column is capped at 12 residues — these two and MBD4 are the
     only records where that cap bites.)
   - CAMTA1/2 (4) keep their annotated CG-1 `DNA_BIND` and shed a span that ran
     all the way to an IPT/TIG 700 aa downstream: CAMTA1 goes from **57—959
     (903 aa, flagged) to 48—203 (156 aa)**. Pfam names no DNA-binding family
     anywhere in CAMTA, so neither arm had an opinion here; this is the
     annotated DBD standing on its own once the hint stops bracketing past it.
   - NFATC4 (2) shrinks from 395—689 to 386—585, dropping the IPT/TIG at
     586—683. **This is the real cost.** Arm 2 found NFATC2 contacting DNA at
     663—666, in the region its paralog annotates as IPT/TIG, so that
     material is not inert. It is recorded as a limitation below rather than
     preserved: keeping the hint meant keeping something wrong in 12 records in
     order to be accidentally right in 2, and NFATC1/2/3 never had it —
     their UniProt entries annotate only RHD.
2. **28 chassis are recoverable** by adding the six Pfam families above to the
   keep rules — STAT most of all. The cost is that keep rules keyed on Pfam
   names need Pfam locations at extraction time, which is a new dependency for
   `06_extract_dbd.py` (the cache now exists, so it is a small one).

   **Done, and the premise was wrong for 10 of the 28.** Checking each of the 28
   against its own UniProt feature table, rather than against Pfam, splits them
   three ways.
   - **4 were already in the UniProt table, filed under a feature type this
     script did not read.** HSF1 is `REGION 15..120; /note="DNA-binding
     domain"` with experimental evidence; TP73 is `REGION 131..310;
     /note="DNA-binding"`. That is why TP73 abstained while TP53 next door
     extracts cleanly — TP53's is a `DNA_BIND` line. Fixed with
     `KEEP_REGION_EXACT`, matched against the whole note and only at ≥ 40 aa:
     UniProt writes the same words for per-residue contact footnotes (RBPJ has
     three, 6—11 aa each) and for negations ("No DNA binding activity...",
     p53's "Basic (repression of DNA-binding)"), so a substring rule would keep
     all of them. Recovers 9 records in total — HSF1, TP73, PURB, LRRFIP1,
     MLLT10 and orthologs.
   - **6 were a naming gap in the zinc-finger whitelist.** A matrin-type finger
     is C2H2-class — Pfam calls ZNF385A's "Zinc-finger of C2H2 type" — and
     "matrin" was simply absent from `ZNF_KEEP_TYPES`. Recovers 8 records (the 6
     plus ZNF318/Zfp318, which Pfam gives no DBD family at all) and costs
     nothing: no existing span changes.
   - **18 genuinely needed Pfam**, and for those the recommendation was right.
     STAT1's entire feature table is one SH2 domain at 573—670. Implemented
     as a fallback that runs only on abstention, never overriding a UniProt
     keep-feature, writing `dbd_extraction='pfam'` so the provenance sits on the
     record. It recovers 30 records: STAT 14, EBF/COE 8, GCFC2/PAXBP1 4,
     RBPJ/CSL 4.

   **Both new paths run only on abstention.** The REGION rule did not, in the
   first attempt, and that was wrong in a way measurement caught and reasoning
   had not: run unconditionally it edits 17 records that were never
   abstentions, against the 9 it recovers. A REGION and a DOMAIN routinely
   annotate the *same* module, so `dna and keep_dom` then fires and relabels
   one module as `bipartite` — ARNT's `REGION 88..128 "DNA-binding"` sits
   inside its own bHLH at 89..142, TFCP2's inside its CP2 domain, PRDM9's
   inside its own C2H2 array. `bipartite` drives the ZnF pitch check, the
   review-queue rank and the per-kind tables in this report, so a false one is
   not cosmetic. ARID4A was the worst single case: `REGION 4..121` plus
   `DOMAIN 309..401 "ARID"`, 188 disordered residues apart, bracketed into a
   416-aa span at 53% coverage — three points above the gate, so
   `dbd_review_flag = 0`. Gated to abstention, all 17 keep their old spans.

   The rule's whitelist is also smaller than first written. Of six entries only
   two can ever fire: measured across all 2,818 cached feature tables,
   `"dna-binding"` appears as a whole note 45 times (17 at ≥ 40 aa) and
   `"dna-binding domain"` twice, while the unhyphenated `"dna binding"` and
   `"dna binding domain"` never appear as a whole note at all. The
   unhyphenated spelling does occur, 10 times, but only inside a longer
   sentence ("Required for DNA binding", "No DNA binding activity...") that an
   exact match correctly refuses — and every
   `"dna-binding; minor groove"` (4) and `"dna-binding; major groove"` (2) is a
   per-residue footnote under the 40-aa floor. The four dead entries were
   removed rather than left implying a coverage the rule does not have.

   One cost stands, recorded rather than fixed: the rule has a lower bound and
   no upper one, so MLLT10's `REGION 311..674` — 364 aa, UniProt's own
   annotation — becomes a 376-aa `dna_bind` span that clears every gate with
   `dbd_review_flag = 0`. It is the rule's largest effect and the justifying
   comment cited only HSF1 (106 aa), TP73 (180 aa) and PURB (227 aa). Rather
   than invent a threshold, every promoted-REGION span now enters the review
   queue whether flagged or not.

   All 28 are now extracted. The dependency is a cached read rather than a
   fetch: with `data/raw/interpro/` empty the fallback is inert and those
   records abstain exactly as before — which is why `ANCHORS` now pins EBF1,
   since an empty cache would otherwise switch the rule off in silence.
3. **`FLANK = 6` trims real DNA contacts.** Raising it to 15 completes 24 of the
   30 partial records; 20 completes 28. The cost is construct length and a
   greater chance of reaching toward an excluded domain — the existing
   pull-back already guards the latter, and the former is real for AAV cargo.

   **Done: `FLANK = 15`, after moving all three gates off the flanked span.**
   The first attempt raised FLANK alone and had to be reverted; the
   measurement that forced the revert is the reason this version works. FLANK
   sat in the *numerator* of all three gates, so raising it loosened every one
   at once, and correcting one of the three produced a library that was
   simultaneously too permissive and too strict.

   - **Too permissive.** `MIN_PLAUSIBLE_SPAN_AA = 40` measures the flanked
     span, so 30 aa of padding carries a span over it with no change in
     evidence. Spans under the 40-aa floor fell **27 —> 5**, so 22 records
     left the set acceptance gate 1 inspects, and 31 records flagged on `main`
     come back unflagged. KMT2C is the clearest: a 4,911-aa protein whose sole
     DNA-binding annotation is `DNA_BIND 34..46 "A.T hook"`, 13 residues,
     emitted at FLANK=15 as a 43-aa span with `dbd_review_flag = 0`. AHCTF1 is
     the same shape. This also hollows out acceptance gate 1, which fails the
     run if any span under 40 aa is unflagged: those 22 records are now 43—56
     aa, so the gate passed while reporting "5 spans under 40 aa".
   - **Too strict.** `MAX_SPAN_FRACTION` / `MAX_SPAN_AA` measure the flanked
     span too, so 14 records lost their span to padding alone — BCL11A,
     TSHZ2/3, ZNF420, ZNF433, ZNF469, ZNF48, ZNF595, ZNF780B, ZNF99, Zfp553,
     Zfp865 and orthologs. ZNF469 shows what that is worth: a 983-aa span at
     **24.9%** of its 3,953-aa protein, dropped because 983 + 18 crosses the
     1,000-aa cap, by a gate whose stated rationale is "the span is the
     protein".
   - **And the one gate that was corrected went too far.** The denominator was
     changed from the flanked span to the feature-bracketed region, which
     removes the flank dependency but is degenerate: `first`/`last` are
     min/max over the same intervals the numerator sums, so a single merged
     interval makes coverage **1.000 by construction** — at least 1,361 of
     2,460 records by direct count. It stops measuring how much of the span is feature and measures
     only whether the features are fragmented. MECP2, HMGA2 and INSM1 (48—50%)
     lose the flag with the artifact, and so do AKNA (42.9%) and NFYB (36.8%),
     which are not artifacts.

   **The fix is one rule: the flank is construct padding, so no gate may see
   it.** Each now measures `ext`, the annotated feature extent — first kept
   feature to last, no flank. The pull-back only ever shortens the flank, so
   the stored span always contains `ext` and the substitution is well defined.
   The two gates then divide the work instead of overlapping:

   | gate | measures | asks |
   |---|---|---|
   | `coverage = covered / ext` | fragmentation of the bracketed region | one module, or two distant ones with a protein between them (MGA) |
   | `ext < MIN_PLAUSIBLE_SPAN_AA` | the feature itself | is there enough annotated feature to be a domain at all (KMT2C) |
   | `ext / len(seq)`, `ext >= MAX_SPAN_AA` | the feature extent | is the span the protein (ZNF462) |

   Coverage still reads 1.000 for a single merged interval, and that is now
   correct rather than a silence — there is no gap to measure, and "is this
   enough feature" is the floor's question, asked of `ext` rather than of the
   padding.

   Measured against the three shipped changes (recommendations 1, 2 and 4 at
   FLANK=6):

   | | #13 | + FLANK=15 on `ext` |
   |---|---|---|
   | records with a span | 2,474 | **2,487** |
   | `dbd_extraction='dispersed'` | 117 | **104** |
   | `dbd_review_flag = 1` | 361 | 354 |
   | stored spans under 40 aa | 27 | 5 |
   | feature extents under 40 aa | 31 | **31, all flagged** |
   | median span | 104 aa | 122 aa |

   **No span is lost.** 13 records that were `dispersed` under the flanked
   measure now keep a span, because their feature extent was never near the
   cap: Repin1, ZNF30, ZNF319, ZNF512B, ZNF573, ZNF585A/B, ZNF749, ZNF799,
   ZNF841, Zfp319, Zfp35, Zfp62. ZNF469 keeps its 983-aa extent at 24.9% of
   the protein; BCL11A keeps its 79.8% span and gains the 6 residues to its
   C-terminus. The records the first attempt silenced stay flagged, now with a
   message that names the real quantity: KMT2C "feature extent 13 aa", AKNA
   "9 aa", NFYB "7 aa".

   **The cost, measured: 5 records lose a coverage flag and nothing replaces
   it** — MECP2/Mecp2, HMGA2/Hmga2 and INSM1, all previously flagged at
   48—50% span coverage. Their features merge into one interval, so the
   fragmentation measure reads 1.000, and their extents (MECP2 188 aa, HMGA2
   59 aa, INSM1 226 aa) are all above the floor. These are the records where "the span is half padding"
   was a true statement about a *small* feature, and no gate asks that question
   any more. They are not errors — MECP2's MBD and INSM1's finger array are
   correctly located — but the flag was doing something, and this is what the
   rule costs. 10 records gain a flag in exchange (CIZ1, FOXP3, LEUTX, PLAG1,
   SOX18, ZNF503 and orthologs), all on a feature extent under the floor.
4. **`domain` is the weakest kind** on both arms (93% recall ≥0.95; 40/51 full
   contact coverage), which is the kind whose spans rest on `KEEP_DOMAIN_HINTS`
   alone. If review effort is finite, spend it there.

   **Done:** `06_extract_dbd.py` writes `data/processed/dbd_review_queue.tsv`,
   ordered `pfam` —> `domain` —> `bipartite` —> `znf_array` —>
   `dna_bind` —> `dispersed`. `pfam` goes first because arm 1 cannot check it:
   those spans *are* Pfam, so the comparison is circular, and nothing has ever
   reviewed them.

   The queue holds flagged records **plus every span whose provenance is a
   source this script did not use before**, flagged or not. Filtering to
   flagged records only — the first attempt — left 21 of the 30 `pfam` spans
   out of the queue entirely, which empties the reason for ranking `pfam`
   first. 472 records: 30 `pfam`, 36 `domain`, 63 `bipartite`, 137
   `znf_array`, 89 `dna_bind`, 117 `dispersed`.

## What all four changes did, end to end

Measured by running the new `extract_one` over the same 2,818 extractable DBD
records.

| | before | after |
|---|---|---|
| records with a span | 2,439 | **2,487** |
| `dbd_extraction='none'` | 262 | **227** |
| `dbd_extraction='dispersed'` | 117 | **104** |
| `dbd_review_flag = 1` | 352 | 354 |
| spans under 40 aa, unflagged | 0 | 0 |
| feature extents under 40 aa, unflagged | 0 | 0 |
| median span | 103 aa | 122 aa |
| total span residues | 448,411 | 500,214 |

**No record loses a span.** 35 abstentions gain one (18 `pfam`, 9 `dna_bind`
via the REGION rule, 8 `znf_array` via `"matrin"`); 13 more keep a span the
old `dispersed` gate took on padding; 12 move onto a different domain
(EBF1—4, RBPJ/RBPJL and orthologs, via the Pfam fallback); 6 shrink
(CAMTA1/2, NFATC4 and orthologs, from dropping `ipt/tig`); every other span
moves by the flank alone.

All three acceptance gates pass against the new spans: 9/9 anchors (`ANCHORS`
gained HSF1 for the REGION rule and EBF1 for the Pfam fallback, so neither can
go inert in silence), 27 spans under 40 aa of which **0 unflagged**, 0 spans
reaching into a KRAB, 0 spans too short for their finger count, and 0 `pfam`
spans carrying a zinc-finger count they have no evidence for.

**Arm 2 has not been re-run end to end for this configuration**, and the
aggregate `185 —> 201` should be treated as provisional until it is. Two
parts of it are directly verified, because the gates do not move boundaries:
the spans scored by hand from the cached coordinates are the spans this
configuration produces. 3MLN puts **30 of 30** Ebf1 DNA contacts inside
`3—261` against 0 of 30 in the old `256—351`, and 3IAG **16 of 16** for
Rbpj inside `33—343` against 0 of 16 in `349—451`. The flank sweep's 63 of
97 recovered contacts also still holds, since those boundaries are unchanged;
what changed is that 14 fewer spans are nulled by the dispersed gate, so the
aggregate can only improve on the figure measured during the first attempt.

## Limitations

- Pfam is not ground truth. It splits C2H2 arrays, so `znf_array` precision of
  0.613 is a property of the comparison, not a defect.
- 220 structures is a biased sample: TFs that crystallise with DNA. The 84%
  full-coverage figure describes those, not the library.
- A structure can capture a *different* DNA-binding module of the same protein,
  as 3 of the 5 zero-coverage cases show. Zero coverage is a prompt to look, not
  a verdict.
- Two accessions failed the InterPro fetch and are absent from arm 1.
- **NFAT needs a keep rule of its own.** Arm 2 puts NFATC2 contacts at
  663—666, past a span ending at 580, and UniProt annotates that region as
  IPT/TIG in NFATC4 and as nothing at all in NFATC1/2/3. All four now stop at
  the RHD, so the family is at least consistent and all four are missing the
  same material. No rule here reaches it: the region is unannotated in three of
  the four entries, and Pfam calls it "Rel homology dimerisation domain", which
  the fallback would keep — but the fallback runs only on abstention, and
  these records have a span from RHD.
- **`MAX_SPAN_FRACTION = 0.80` sits at the boundary for at least one real
  chassis.** BCL11A's span is 666 aa, **79.8%** of the protein, and arm 2 found
  all 17 of its mapped DNA contacts (6U9Q, 1.83 Å) inside it — a verified
  span two tenths of a percent under the cap. It keeps that span here, and an
  earlier draft of this report was wrong about why it might not: it listed
  BCL11A as a cost of `MAX_SPAN_FRACTION` and said the remedy was to revert
  FLANK or lift the threshold. In fact *both* of the first attempt's changes
  destroyed it independently — the flank pushed it to 81.9%, and the
  then-ungated REGION rule promoted `REGION 737..835 "DNA-binding"`, dragging
  the span to the C-terminus and over 80% even at FLANK=6. Neither is a
  property of the threshold. With the gates measuring the feature extent it
  keeps its span (`155—835`) and its 17 verified contacts, and the threshold
  is back to meaning what its name says. It remains two tenths of a percent
  from the cut, so it is still the record to watch.
- Contacts are computed on the deposited coordinates of one chain mapping, with
  no symmetry expansion, so a contact made only by a crystallographic neighbour
  is not counted.

## Reproduce

```bash
python scripts/07_validate_dbd_external.py --fetch     # ~2,800 cached GETs
python scripts/07_validate_dbd_external.py --survey    # Pfam family census
python scripts/07_validate_dbd_external.py --report    # arm 1
python scripts/07_validate_dbd_external.py --pdb --limit 220   # arm 2
```

Per-record tables land in `data/processed/dbd_external_validation.tsv` and
`data/processed/dbd_pdb_contact_validation.tsv` (both gitignored; this report is
the committed artefact).

Note the ordering, now that the extraction has a Pfam fallback: `--fetch`
populates `data/raw/interpro/`, and `06_extract_dbd.py` reads that cache. Run
the fetch **before** the extraction on a fresh checkout, or 18 records (the
STAT family among them) will abstain and gate 2's EBF1 anchor will fail. The
numbers in this report describe the spans as they were *before* the four
changes above; re-running `--report` and `--pdb` after an extraction re-scores
them.
