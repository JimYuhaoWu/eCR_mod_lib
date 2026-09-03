#!/usr/bin/env python3
"""Extract the DNA-binding region from full-length TFs, for the 8 fib->beta eCRs.

eCR_mod_lib stores whole UniProt sequences under type='DBD' (01_fetch_dbd.py
takes the entire sequence), so a fusion built from those records is a
full-length-TF fusion, not an eCR. That matters here beyond tidiness: ZNF793 and
ZNF850 are OPENERS carrying native KRAB domains, which recruit KAP1/SETDB1 and
would silence the very loci SS18-N70 is meant to open -- self-antagonism inside
one construct.

Rule: keep the DNA-contacting apparatus, drop everything else.
  - C2H2 arrays -> first ZnF start .. last ZnF end (inter-finger TGEKP linkers are
    part of the array and must stay), plus FLANK residues either side.
  - Non-ZnF -> the union of DNA_BIND and DNA-binding-partner domains, e.g. HNF1A
    binds as a bipartite POU-specific + homeodomain unit, so taking the homeobox
    alone would remove half the interface.
  - Explicitly excluded: KRAB, SCAN, BTB/POZ, and interaction/activation regions.
"""
import json, sys, urllib.request

FLANK = 6
EXCLUDE_DOMAIN = ("krab", "scan", "btb", "poz")
KEEP_DOMAIN_HINTS = ("pou", "homeo")     # DNA-contacting partner domains

TARGETS = [   # (gene, uniprot, role)
    ("HNF1A",  "P20823", "opener"),
    ("EGR3",   "Q06889", "opener"),
    ("ZNF793", "Q6ZN11", "opener"),
    ("ZNF850", "A8MQ14", "opener"),
    ("ZNF540", "Q8NDQ6", "closer"),
    ("ZNF444", "Q8N0Y2", "closer"),
    ("HINFP",  "Q9BQA5", "closer"),
]

def fetch(uid):
    url = (f"https://rest.uniprot.org/uniprotkb/{uid}.json"
           "?fields=ft_zn_fing,ft_dna_bind,ft_domain,sequence")
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.load(r)

out = []
for gene, uid, role in TARGETS:
    d = fetch(uid)
    seq = d["sequence"]["value"]
    zn, dna, keep_dom, dropped = [], [], [], []
    for f in d.get("features", []):
        s, e = f["location"]["start"]["value"], f["location"]["end"]["value"]
        desc = (f.get("description") or "").lower()
        if f["type"] == "Zinc finger":
            (zn if "degenerate" not in desc else zn).append((s, e))
        elif f["type"] == "DNA binding":
            dna.append((s, e))
        elif f["type"] == "Domain":
            if any(k in desc for k in EXCLUDE_DOMAIN):
                dropped.append((f.get("description"), s, e))
            elif any(k in desc for k in KEEP_DOMAIN_HINTS):
                keep_dom.append((s, e))
            else:
                dropped.append((f.get("description"), s, e))
    spans = zn + dna + keep_dom
    if not spans:
        print(f"{gene}: NO DNA-binding feature found -- skipping"); continue
    lo = max(1, min(s for s, _ in spans) - FLANK)
    hi = min(len(seq), max(e for _, e in spans) + FLANK)
    # a dropped domain must not sit inside the kept span, or we would carry it anyway
    inside = [d0 for d0, s, e in dropped if s >= lo and e <= hi]
    sub = seq[lo-1:hi]
    out.append(dict(gene_name=gene, uniprot_id=uid, role=role,
                    full_len=len(seq), dbd_start=lo, dbd_end=hi, dbd_len=len(sub),
                    n_znf=len(zn), dropped=";".join(d0 or "?" for d0, _, _ in dropped) or "-",
                    inside_kept_span=";".join(inside) or "-", sequence_aa=sub))
    print(f"{gene:7s} {role:6s} full {len(seq):5d} -> DBD {lo}-{hi} ({len(sub)} aa, "
          f"{100*len(sub)/len(seq):.0f}%)  ZnF={len(zn)}  dropped: "
          f"{', '.join(d0 or '?' for d0,_,_ in dropped) or 'none'}"
          + (f"  !! STILL INSIDE: {inside}" if inside else ""))

import csv
with open("dbd_extracted.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(out)
print(f"\nwrote dbd_extracted.tsv ({len(out)} records)")
