"""
scripts/02_seed_ed.py
---------------------
Populate Effector Domain (ED) entries in the module library from two sources:

  1. data/manual/ed_curated.yaml
     – Manually curated canonical EDs (VP64, VPR, KRAB, p65-TAD, etc.)
     – Sequences fetched from UniProt for entries with a uniprot_id;
       sequence_manual used directly for synthetic constructs.

  2. Published activation-domain screens (Alerasool et al. 2022, DelRosso et al. 2023)
     – These require manual download of supplementary tables from the journal
       website; the script checks for them and processes them if present.
     – Instructions for manual download are printed if files are missing.

Outputs
-------
  data/raw/uniprot/{id}.fasta         cached sequences
  data/processed/ed_raw.tsv           merged ED records
  library/build_manifest.json         updated build record

Run
---
  python scripts/02_seed_ed.py [--config config/sources.yaml] [--dry-run]
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests
import yaml
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))
from utils import get_logger, write_build_manifest

ROOT = Path(__file__).resolve().parents[1]
LOGS_DIR = ROOT / "logs"
RAW_DIR = ROOT / "data" / "raw"
PROCESSED_DIR = ROOT / "data" / "processed"
MANUAL_DIR = ROOT / "data" / "manual"

# Mapping: curated subtype string → canonical group for chromatin_state_effect
_EFFECT_MAP = {
    "activator": "open / H3K27ac deposition",
    "repressor": "close / H3K9me2-3 or H3K27me3",
    "dual": "context-dependent",
}


# ── UniProt sequence fetch (shared with 01_fetch_dbd.py) ─────────────────────

def _get(url: str, session: requests.Session, timeout: int = 20) -> requests.Response | None:
    import time
    for attempt in range(3):
        try:
            r = session.get(url, timeout=timeout)
            if r.status_code == 200:
                return r
            if r.status_code == 404:
                return None
        except requests.RequestException:
            pass
        time.sleep(2 ** attempt)
    return None


def fetch_domain_sequence(uniprot_id: str, domain_residues: str | None,
                           session: requests.Session, log) -> str | None:
    """
    Fetch the full protein sequence from UniProt, then slice to domain_residues
    if provided (e.g. '411-490').  Returns the domain AA sequence or None.
    """
    if not uniprot_id:
        return None

    cache = RAW_DIR / "uniprot" / f"{uniprot_id}.fasta"
    if cache.exists():
        fasta = cache.read_text()
    else:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
        resp = _get(url, session, timeout=20)
        if resp is None:
            log.warning(f"[UniProt] could not fetch {uniprot_id}")
            return None
        cache.write_text(resp.text, encoding="utf-8")
        fasta = resp.text

    # Parse FASTA
    lines = fasta.strip().splitlines()
    full_seq = "".join(l.strip() for l in lines if not l.startswith(">"))
    if not full_seq:
        return None

    # Slice domain
    if domain_residues:
        try:
            start, end = (int(x) for x in domain_residues.split("-"))
            return full_seq[start - 1 : end]   # 1-indexed, inclusive
        except Exception:
            log.warning(f"[UniProt] could not parse domain_residues='{domain_residues}'")
    return full_seq


# ── Process curated YAML ──────────────────────────────────────────────────────

def process_curated_eds(cfg_path: Path, session: requests.Session,
                         log, dry_run: bool) -> list[dict]:
    with open(cfg_path, encoding="utf-8") as f:
        raw = yaml.safe_load(f)

    records = []
    now = datetime.now(timezone.utc).isoformat()

    for category in ("activators", "repressors"):
        entries = raw.get(category, [])
        for entry in tqdm(entries, desc=f"Curated {category}", ncols=80):
            name = entry["name"]
            subtype = entry.get("subtype", "activator")
            uniprot_id = entry.get("uniprot_id") or None
            domain_residues = entry.get("domain_residues") or None
            seq_manual = entry.get("sequence_manual") or None

            # Sequence resolution order: manual > UniProt slice
            seq = seq_manual
            if seq is None and uniprot_id and not dry_run:
                seq = fetch_domain_sequence(uniprot_id, domain_residues, session, log)
            length = len(seq) if seq else None

            module_id = f"ED_{name.upper().replace('-','_').replace(' ','_')}"

            record = {
                "module_id": module_id,
                "type": "ED",
                "subtype": subtype,
                "name": name,
                "organism": entry.get("organism", ""),
                "source_species": entry.get("organism", ""),
                "gene_symbol": None,
                "uniprot_id": uniprot_id,
                "sequence_aa": seq,
                "length_aa": length,
                "target_or_mechanism": entry.get("mechanism", ""),
                "quantitative_metric": entry.get("activity_score"),
                "quantitative_metric_label": entry.get("activity_metric"),
                "quantitative_metric_source": entry.get("activity_source"),
                "validation_level": entry.get("validation_level", "screen-validated"),
                "source": "Curated: " + entry.get("doi", ""),
                "source_doi": entry.get("doi"),
                "source_version": "manual_v1",
                "source_date": now[:10],
                "known_interactors": entry.get("known_interactors", ""),
                "linker_notes": entry.get("linker_notes", ""),
                "engineering_compatibility": entry.get("engineering_compatibility", ""),
                "chromatin_state_effect": _EFFECT_MAP.get(subtype, ""),
                "notes": str(entry.get("notes", "")).strip(),
                "date_added": now,
                "date_modified": now,
            }
            records.append(record)
            log.debug(f"[curated ED] {name}  seq={'yes' if seq else 'no'}  len={length}")

    return records


# ── Process Alerasool / DelRosso screen data ──────────────────────────────────

SCREEN_FILES = {
    "Alerasool_2022": {
        "path": MANUAL_DIR / "Alerasool_2022_SupTable.xlsx",
        "doi": "10.1016/j.molcel.2021.12.005",
        "sheet": "tAD-seq",
        "col_name": None,           # built from Gene + Fragment below
        "col_gene": "Gene",
        "col_fragment": "Fragment",
        "col_sequence": "Sequence",
        "col_score": "logFC high GFP",
        "col_hit": "Hit high GFP",
        "hit_value": "hit",
        "notes": (
            "Download mmc2.xlsx (Supplementary Table 2) from "
            "https://doi.org/10.1016/j.molcel.2021.12.005 "
            "and save as data/manual/Alerasool_2022_SupTable.xlsx"
        ),
    },
    "DelRosso_2023_activators": {
        # F8: this paper scored the same domains in both arms, ingested as two
        # SCREEN_FILES entries. `readout` says which arm a row is, so the
        # opposite `subtype` on a shared sequence reads as two measurements
        # rather than a contradiction.
        "readout": "activation",
        "path": MANUAL_DIR / "DelRosso_2023_SupTable.xlsx",
        "doi": "10.1038/s41586-023-05906-y",
        "sheet": "Activation Domains",
        "col_name": None,
        "col_gene": "HGNC symbol",
        "col_fragment": "Domain",
        "col_sequence": "Sequence",
        "col_score": "Max avg act",
        "col_hit": None,
        "hit_value": None,
        "subtype_override": "activator",
        "notes": (
            "Download Supplementary Table 6 from "
            "https://doi.org/10.1038/s41586-023-05906-y "
            "and save as data/manual/DelRosso_2023_SupTable.xlsx"
        ),
    },
    "DelRosso_2023_repressors": {
        # F8: this paper scored the same domains in both arms, ingested as two
        # SCREEN_FILES entries. `readout` says which arm a row is, so the
        # opposite `subtype` on a shared sequence reads as two measurements
        # rather than a contradiction.
        "readout": "repression",
        "path": MANUAL_DIR / "DelRosso_2023_SupTable.xlsx",
        "doi": "10.1038/s41586-023-05906-y",
        "sheet": "Repression Domains",
        "col_name": None,
        "col_gene": "HGNC symbol",
        "col_fragment": "Domain",
        "col_sequence": "Sequence",
        "col_score": "Max avg pEF",
        "col_score_fallback": "Max avg PGK",
        "col_hit": None,
        "hit_value": None,
        "subtype_override": "repressor",
        "notes": (
            "Download Supplementary Table 6 from "
            "https://doi.org/10.1038/s41586-023-05906-y "
            "and save as data/manual/DelRosso_2023_SupTable.xlsx"
        ),
    },
    # Staller 2022 — directed mutational scan + novel predicted ADs tested in human cells.
    # mmc5 (activities) is merged on-the-fly with mmc4 (UniProt IDs) inside process_screen_data.
    "Staller_2022_predictions": {
        "path": MANUAL_DIR / "Staller_2022_mmc5.csv",
        "doi": "10.1016/j.cels.2022.01.002",
        "sheet": None,
        "col_name": None,
        "col_gene": "GeneName",
        "col_fragment_start": "Start",
        "col_fragment_end": "End",
        "col_sequence": "ProteinRegionSeq",
        "col_score": "Activity_Zscore_mean",
        "min_score_threshold": 0.5,
        "row_prefilter": {"col": "RegionType", "value": "Prediction"},
        "col_hit": None,
        "hit_value": None,
        "subtype_override": "activator",
        "validation_level_override": "screen-validated",
        # Merge mmc4 to add UniProt IDs; joined on ProteinRegionSeq
        "merge_file": MANUAL_DIR / "Staller_2022_mmc4.csv",
        "merge_on": "ProteinRegionSeq",
        "merge_extract": {
            "col_uniprot": "uniprotID",
            "col_gene_symbol": ("GeneName", r"\|([^|]+)_HUMAN"),
        },
        "notes": (
            "Save mmc5.csv and mmc4.csv from https://doi.org/10.1016/j.cels.2022.01.002 "
            "as data/manual/Staller_2022_mmc5.csv and Staller_2022_mmc4.csv"
        ),
    },
    "Staller_2022_wt_domains": {
        "path": MANUAL_DIR / "Staller_2022_mmc2.csv",
        "doi": "10.1016/j.cels.2022.01.002",
        "sheet": None,
        "col_name": "Variant_Name",
        "col_gene": None,
        "col_fragment": None,
        "col_sequence": "ADseq",
        "col_score": "Activity_Mean_MSS18",
        "col_hit": "WT",
        "col_hit_values": [True],
        "subtype_override": "activator",
        "validation_level_override": "screen-validated",
        "notes": (
            "Download mmc2.csv from https://doi.org/10.1016/j.cels.2022.01.002 "
            "and save as data/manual/Staller_2022_mmc2.csv"
        ),
    },
    # Compendium of human TF effector domains — experimentally characterised ADs/RDs.
    # Filter: Activity H or M (drop L, NaN, S). Validation: ChIP-validated.
    "Compendium_2021_activators": {
        "path": MANUAL_DIR / "Compendium_2021_SupTable.xlsx",
        "doi": "10.1016/j.molcel.2021.11.007",
        "sheet": "Table S2",
        "col_name": None,
        "col_gene": "TF name",
        "col_fragment": "Effector domain ID",
        "col_sequence": "Sequence",
        "col_score": "Activity (H, M or L)",
        "score_map": {"H": 1.0, "M": 0.5},
        "col_hit": "Activity (H, M or L)",
        "col_hit_values": ["H", "M"],
        "row_prefilter": {"col": "Domain type", "value": "AD"},
        "subtype_override": "activator",
        "validation_level_override": "ChIP-validated",
        "notes": (
            "Download Supplementary Table 2 (mmc8.xlsx) from "
            "https://doi.org/10.1016/j.molcel.2021.11.007 "
            "and save as data/manual/Compendium_2021_SupTable.xlsx"
        ),
    },
    "Compendium_2021_repressors": {
        "path": MANUAL_DIR / "Compendium_2021_SupTable.xlsx",
        "doi": "10.1016/j.molcel.2021.11.007",
        "sheet": "Table S2",
        "col_name": None,
        "col_gene": "TF name",
        "col_fragment": "Effector domain ID",
        "col_sequence": "Sequence",
        "col_score": "Activity (H, M or L)",
        "score_map": {"H": 1.0, "M": 0.5},
        "col_hit": "Activity (H, M or L)",
        "col_hit_values": ["H", "M"],
        "row_prefilter": {"col": "Domain type", "value": "RD"},
        "subtype_override": "repressor",
        "validation_level_override": "ChIP-validated",
        "notes": (
            "Download Supplementary Table 2 (mmc8.xlsx) from "
            "https://doi.org/10.1016/j.molcel.2021.11.007 "
            "and save as data/manual/Compendium_2021_SupTable.xlsx"
        ),
    },
    # High-throughput nuclear domain screens (Pfam-domain library).
    # DOI: TODO — confirm from paper. Validation: screen-validated.
    "HiTEff_activators": {
        # F8: this paper scored the same domains in both arms, ingested as two
        # SCREEN_FILES entries. `readout` says which arm a row is, so the
        # opposite `subtype` on a shared sequence reads as two measurements
        # rather than a contradiction.
        "readout": "activation",
        "path": MANUAL_DIR / "HiTEff_SupTable.xlsx",
        "doi": "10.1016/j.cell.2020.11.024",
        "sheet": "NucAct_data",
        "col_name": None,
        "col_gene": "Gene entry name",
        "col_fragment": "Domain ID",
        "col_sequence": "Extended Domain sequence",
        "col_score": "Avg Act",
        "col_hit": "Hit",
        "col_hit_values": [True],
        "subtype_override": "activator",
        "notes": (
            "Download mmc4.xlsx from https://doi.org/10.1016/j.cell.2020.11.024 "
            "and save as data/manual/HiTEff_SupTable.xlsx"
        ),
    },
    "HiTEff_repressors": {
        # F8: this paper scored the same domains in both arms, ingested as two
        # SCREEN_FILES entries. `readout` says which arm a row is, so the
        # opposite `subtype` on a shared sequence reads as two measurements
        # rather than a contradiction.
        "readout": "repression",
        "path": MANUAL_DIR / "HiTEff_SupTable.xlsx",
        "doi": "10.1016/j.cell.2020.11.024",
        "sheet": "NucRepr_data",
        "col_name": None,
        "col_gene": "Gene entry name",
        "col_fragment": "Domain ID",
        "col_sequence": "Extended Domain sequence",
        "col_score": "Avg ReprD13",
        "col_hit": "Hit",
        "col_hit_values": [True],
        "subtype_override": "repressor",
        "notes": (
            "Download mmc4.xlsx from https://doi.org/10.1016/j.cell.2020.11.024 "
            "and save as data/manual/HiTEff_SupTable.xlsx"
        ),
    },

    # ── Tycko 2025 NatBiotech — compact effectors, multi-context HT-recruit ──────
    "Tycko_2025_ADs": {
        "path": MANUAL_DIR / "Tycko_2025_MOESM3.xlsx",
        "doi": "10.1038/s41587-024-02442-6",
        "sheet": "ST3 CD2 ADs",
        "col_name": None,
        "col_gene": "HGNC symbol",
        "col_fragment": "Domain",
        "col_uniprot": "UniProt ID",
        "col_sequence": "Sequence",
        "col_score": "Max Avg Act at CD2 log2(OFF:ON)",
        "col_hit": None,
        "subtype_override": "activator",
        "validation_level_override": "screen-validated",
        "notes": (
            "Download 41587_2024_2442_MOESM3_ESM.xlsx from "
            "https://doi.org/10.1038/s41587-024-02442-6 "
            "and save as data/manual/Tycko_2025_MOESM3.xlsx"
        ),
    },
    "Tycko_2025_RDs": {
        "path": MANUAL_DIR / "Tycko_2025_MOESM3.xlsx",
        "doi": "10.1038/s41587-024-02442-6",
        "sheet": "ST3 CD43 RDs",
        "col_name": None,
        "col_gene": "HGNC symbol",
        "col_fragment": "Domain",
        "col_uniprot": "UniProt ID",
        "col_sequence": "Domain sequence",
        "col_score": "Max Avg repression at CD43 log2(OFF:ON)",
        "col_hit": None,
        "subtype_override": "repressor",
        "validation_level_override": "screen-validated",
        "notes": (
            "Download 41587_2024_2442_MOESM3_ESM.xlsx from "
            "https://doi.org/10.1038/s41587-024-02442-6 "
            "and save as data/manual/Tycko_2025_MOESM3.xlsx"
        ),
    },

    # ── Ludwig 2023 Cell Systems — viral transcriptional effectors ────────────────
    # Three sheets (vTR=adenovirus/other, CoV=coronavirus, HHV=herpesvirus).
    # Effect column determines subtype. Max Sequence = 80aa core domain.
    # organism_override: these are viral proteins, not human (F6.2).
    "Ludwig_2023_vTR": {
        "path": MANUAL_DIR / "Ludwig_2023_mmc4.xlsx",
        "doi": "10.1016/j.cels.2023.05.008",
        "sheet": "vTR Domains",
        "col_name": None,
        "col_gene": "Gene",
        "col_fragment": "Tile ID",
        "col_sequence": "Max Sequence",
        "col_score": "Max Score",
        "col_hit": None,
        "col_subtype": "Effect",
        "subtype_map": {"Activation": "activator", "Repression": "repressor"},
        "validation_level_override": "screen-validated",
        "organism_override": "viral",
        "notes": (
            "Download mmc4.xlsx from https://doi.org/10.1016/j.cels.2023.05.008 "
            "and save as data/manual/Ludwig_2023_mmc4.xlsx"
        ),
    },
    "Ludwig_2023_CoV": {
        "path": MANUAL_DIR / "Ludwig_2023_mmc4.xlsx",
        "doi": "10.1016/j.cels.2023.05.008",
        "sheet": "CoV Domains",
        "col_name": None,
        "col_gene": "Gene",
        "col_fragment": "Tile ID",
        "col_sequence": "Max Sequence",
        "col_score": "Max Score",
        "col_hit": None,
        "col_subtype": "Effect",
        "subtype_map": {"Activation": "activator", "Repression": "repressor"},
        "validation_level_override": "screen-validated",
        "organism_override": "viral",
        "notes": (
            "Download mmc4.xlsx from https://doi.org/10.1016/j.cels.2023.05.008 "
            "and save as data/manual/Ludwig_2023_mmc4.xlsx"
        ),
    },
    "Ludwig_2023_HHV": {
        "path": MANUAL_DIR / "Ludwig_2023_mmc4.xlsx",
        "doi": "10.1016/j.cels.2023.05.008",
        "sheet": "HHV Domains",
        "col_name": None,
        "col_gene": "Gene",
        "col_fragment": "Tile ID",
        "col_sequence": "Max Sequence",
        "col_score": "Max Score",
        "col_hit": None,
        "col_subtype": "Effect",
        "subtype_map": {"Activation": "activator", "Repression": "repressor"},
        "validation_level_override": "screen-validated",
        "organism_override": "viral",
        "notes": (
            "Download mmc4.xlsx from https://doi.org/10.1016/j.cels.2023.05.008 "
            "and save as data/manual/Ludwig_2023_mmc4.xlsx"
        ),
    },

    # ── TENet 2024 preprint — 54 WT repressor domains from DMS study ─────────────
    # media-6 has all clinical variants; deduplicate on repressor_domain to get
    # one WT entry per domain. log2_offon_wt is the WT activity from HT-recruit.
    "TENet_2024_RDs": {
        "path": MANUAL_DIR / "TENet_2024_media-6.xlsx",
        "doi": "10.1101/2024.09.21.614253",
        "sheet": "RD-DMS",
        "col_name": "repressor_domain",
        "col_gene": "gene_hgnc_id",
        "col_fragment": None,
        "col_uniprot": "uniprot_id",
        "col_sequence": "wt_aa",
        "col_score": "log2_offon_wt",
        "col_hit": None,
        "deduplicate_on": "repressor_domain",
        "subtype_override": "repressor",
        "validation_level_override": "screen-validated",
        "notes": (
            "Download media-6.xlsx from https://doi.org/10.1101/2024.09.21.614253 "
            "and save as data/manual/TENet_2024_media-6.xlsx. "
            "Preprint — update DOI when published."
        ),
    },

    # ── Mukund 2023 Cell Systems — Pfam nuclear domain screen ────────────────────
    # mmc5 contains 4,168 short nuclear domains; 4,113 are identical labels to
    # HiTEff (same assay platform, same oligo library). Only the 55 unique entries
    # are loaded by excluding any label present in the HiTEff reference sheets.
    # Two entries: one for activator scores (Avg Act), one for repressor (Avg ReprD13).
    "Mukund_2023_activators": {
        # F8: this paper scored the same domains in both arms, ingested as two
        # SCREEN_FILES entries. `readout` says which arm a row is, so the
        # opposite `subtype` on a shared sequence reads as two measurements
        # rather than a contradiction.
        "readout": "activation",
        "path": MANUAL_DIR / "Mukund_2023_mmc5.xlsx",
        "doi": "10.1016/j.cels.2023.07.001",
        "sheet": "Sheet1",
        "col_name": None,
        "col_gene": "Gene entry name",
        "col_fragment": "Domain ID",
        "col_label": "label",
        "col_sequence": "Domain sequence",
        "col_score": "Avg Act",
        "col_hit": None,
        "row_prefilter": {"col": "Category", "value": "Short nuclear domain"},
        "subtype_override": "activator",
        "validation_level_override": "screen-validated",
        "exclude_labels_from": [
            {"path": MANUAL_DIR / "HiTEff_SupTable.xlsx", "sheet": "NucAct_data",  "col": "label"},
            {"path": MANUAL_DIR / "HiTEff_SupTable.xlsx", "sheet": "NucRepr_data", "col": "label"},
        ],
        "notes": (
            "Download mmc5.xlsx from https://doi.org/10.1016/j.cels.2023.07.001 "
            "and save as data/manual/Mukund_2023_mmc5.xlsx"
        ),
    },
    "Mukund_2023_repressors": {
        # F8: this paper scored the same domains in both arms, ingested as two
        # SCREEN_FILES entries. `readout` says which arm a row is, so the
        # opposite `subtype` on a shared sequence reads as two measurements
        # rather than a contradiction.
        "readout": "repression",
        "path": MANUAL_DIR / "Mukund_2023_mmc5.xlsx",
        "doi": "10.1016/j.cels.2023.07.001",
        "sheet": "Sheet1",
        "col_name": None,
        "col_gene": "Gene entry name",
        "col_fragment": "Domain ID",
        "col_label": "label",
        "col_sequence": "Domain sequence",
        "col_score": "Avg ReprD13",
        "col_hit": None,
        "row_prefilter": {"col": "Category", "value": "Short nuclear domain"},
        "subtype_override": "repressor",
        "validation_level_override": "screen-validated",
        "exclude_labels_from": [
            {"path": MANUAL_DIR / "HiTEff_SupTable.xlsx", "sheet": "NucAct_data",  "col": "label"},
            {"path": MANUAL_DIR / "HiTEff_SupTable.xlsx", "sheet": "NucRepr_data", "col": "label"},
        ],
        "notes": (
            "Download mmc5.xlsx from https://doi.org/10.1016/j.cels.2023.07.001 "
            "and save as data/manual/Mukund_2023_mmc5.xlsx"
        ),
    },

    # ── Kristof 2025 Genome Biology — engineered CRISPRi repressors ──────────────
    # ST1 lists repressor domain constructs with AA sequences used in the paper.
    # Treated as ChIP-validated (experimentally characterised, not a pooled screen).
    "Kristof_2025_RDs": {
        "path": MANUAL_DIR / "Kristof_2025_MOESM1.xlsx",
        "doi": "10.1186/s13059-025-03640-4",
        "sheet": "Supplementary Table S1",
        "header_row": 2,          # row 0=blank, row 1=section title, row 2=column headers
        "col_name": "Repressor Domain",
        "col_gene": None,
        "col_fragment": None,
        "col_sequence": "AA Sequence",
        "col_score": None,
        "col_hit": None,
        "subtype_override": "repressor",
        "validation_level_override": "ChIP-validated",
        "notes": (
            "Download 13059_2025_3640_MOESM1_ESM.xlsx from "
            "https://doi.org/10.1186/s13059-025-03640-4 "
            "and save as data/manual/Kristof_2025_MOESM1.xlsx"
        ),
    },
}


def process_screen_data(log) -> list[dict]:
    """
    Parse available screen supplementary tables. Prints instructions for
    missing files; skips silently if absent.
    """
    records = []
    now = datetime.now(timezone.utc).isoformat()

    for screen_name, info in SCREEN_FILES.items():
        if not info["path"].exists():
            log.warning(
                f"[{screen_name}] file not found: {info['path']}\n"
                f"  → {info['notes']}"
            )
            continue

        log.info(f"[{screen_name}] parsing {info['path']}")
        try:
            sheet = info.get("sheet")
            header_row = info.get("header_row", 0)
            if info["path"].suffix == ".xlsx":
                df = pd.read_excel(info["path"], sheet_name=sheet, header=header_row)
            elif info["path"].suffix == ".csv":
                df = pd.read_csv(info["path"])
            else:
                df = pd.read_csv(info["path"], sep="\t")
        except Exception as e:
            log.error(f"[{screen_name}] parse error: {e}")
            continue

        # Exclude rows whose label already appears in reference file(s)
        exclude_cfg = info.get("exclude_labels_from")
        col_label = info.get("col_label")
        if exclude_cfg and col_label and col_label in df.columns:
            excluded_labels: set = set()
            for exc in exclude_cfg:
                try:
                    edf = pd.read_excel(exc["path"], sheet_name=exc.get("sheet"))
                    excluded_labels |= set(edf[exc["col"]].dropna())
                except Exception as e:
                    log.warning(f"[{screen_name}] exclude_labels_from load failed: {e}")
            before = len(df)
            df = df[~df[col_label].isin(excluded_labels)]
            log.info(f"[{screen_name}] excluded {before - len(df)} duplicate labels, {len(df)} unique remain")

        # Optional merge with a second file (e.g. to add UniProt IDs)
        merge_file = info.get("merge_file")
        if merge_file and Path(merge_file).exists():
            try:
                mdf = pd.read_csv(merge_file) if str(merge_file).endswith(".csv") else pd.read_excel(merge_file)
                merge_on = info["merge_on"]
                extract = info.get("merge_extract", {})
                # Extract/rename columns from the merge file before joining
                for target_key, src in extract.items():
                    if isinstance(src, tuple):
                        col, pattern = src
                        mdf[target_key] = mdf[col].str.extract(pattern)
                    else:
                        mdf = mdf.rename(columns={src: target_key})
                keep = [merge_on] + list(extract.keys())
                mdf = mdf[[c for c in keep if c in mdf.columns]].drop_duplicates(merge_on)
                df = df.merge(mdf, on=merge_on, how="left")
                # Update col_uniprot to the merged column name
                if "col_uniprot" in extract:
                    info = {**info, "col_uniprot": "col_uniprot"}
            except Exception as e:
                log.warning(f"[{screen_name}] merge failed: {e}")

        # Apply row pre-filter (e.g. restrict to a specific domain type)
        pf = info.get("row_prefilter")
        if pf and pf["col"] in df.columns:
            df = df[df[pf["col"]].astype(str).str.strip() == str(pf["value"])]

        # Filter to hits — supports single value (str) or list of values
        col_hit = info.get("col_hit")
        hit_val = info.get("hit_value")
        hit_vals = info.get("col_hit_values")
        if col_hit and col_hit in df.columns:
            if hit_vals is not None:
                # Coerce to str for comparison so bool True matches string "True"
                str_vals = [str(v) for v in hit_vals]
                df = df[df[col_hit].astype(str).isin(str_vals)]
            elif hit_val:
                df = df[df[col_hit].astype(str).str.strip().str.lower() == hit_val.lower()]

        # Deduplicate on a column (e.g. keep one WT row per domain)
        dedup_col = info.get("deduplicate_on")
        if dedup_col and dedup_col in df.columns:
            df = df.drop_duplicates(subset=[dedup_col], keep="first")

        # Numeric score threshold filter
        min_score = info.get("min_score_threshold")
        col_q_thresh = info.get("col_score")
        if min_score is not None and col_q_thresh and col_q_thresh in df.columns:
            df = df[pd.to_numeric(df[col_q_thresh], errors="coerce") >= min_score]

        col_n = info.get("col_name")
        col_s = info.get("col_sequence")
        col_q = info.get("col_score")
        col_gene = info.get("col_gene")
        col_frag = info.get("col_fragment")
        col_frag_start = info.get("col_fragment_start")
        col_frag_end = info.get("col_fragment_end")
        # col_uniprot: direct column name, or from merge_extract
        merge_extract = info.get("merge_extract", {})
        col_uniprot_col = "col_uniprot" if "col_uniprot" in merge_extract else info.get("col_uniprot")

        for _, row in df.iterrows():
            # Build name: gene+fragment, gene+start-end, col_name, or just gene
            if col_n and col_n in df.columns:
                name = str(row.get(col_n, "")).strip()
            elif col_gene and col_frag and col_frag in df.columns:
                name = f"{str(row.get(col_gene, '')).strip()}_{str(row.get(col_frag, '')).strip()}"
            elif col_gene and col_frag_start and col_frag_end:
                g = str(row.get(col_gene, "")).strip()
                s = str(row.get(col_frag_start, "")).strip().split(".")[0]
                e = str(row.get(col_frag_end, "")).strip().split(".")[0]
                name = f"{g}_{s}-{e}"
            elif col_gene:
                name = str(row.get(col_gene, "")).strip()
            else:
                name = ""
            if not name or name.lower() in ("nan", "_", "nan_nan", "nan_nan-nan"):
                continue

            seq = str(row.get(col_s, "")).strip() if col_s and col_s in df.columns else None
            if seq in (None, "nan", ""):
                seq = None
            def _parse_score(val):
                try:
                    f = float(val)
                    return None if pd.isna(f) else f
                except (ValueError, TypeError):
                    return None

            score = None
            col_score_used = col_q
            score_map = info.get("score_map")
            if col_q and col_q in df.columns:
                raw_val = row[col_q]
                if score_map and str(raw_val).strip() in score_map:
                    score = score_map[str(raw_val).strip()]
                else:
                    score = _parse_score(raw_val)
            # Fall back to secondary score column if primary is null
            if score is None:
                col_fb = info.get("col_score_fallback")
                if col_fb and col_fb in df.columns:
                    fb_score = _parse_score(row[col_fb])
                    if fb_score is not None:
                        score = fb_score
                        col_score_used = col_fb

            if info.get("subtype_override"):
                subtype = info["subtype_override"]
            elif info.get("col_subtype") and info.get("subtype_map"):
                raw_st = str(row.get(info["col_subtype"], "")).strip()
                subtype = info["subtype_map"].get(raw_st, "activator")
            else:
                is_repressor = "repression" in screen_name.lower() or "repression" in (col_q or "")
                subtype = "repressor" if is_repressor else "activator"

            # UniProt from dedicated column if specified
            uniprot_from_col = None
            if col_uniprot_col and col_uniprot_col in df.columns:
                v = str(row.get(col_uniprot_col, "")).strip()
                if v and v.lower() not in ("nan", ""):
                    uniprot_from_col = v

            module_id = f"ED_{screen_name.upper()}_{name.upper().replace(' ','_')[:40]}"
            record = {
                "module_id": module_id,
                "type": "ED",
                "subtype": subtype,
                "name": name,
                # F6.2: screen-derived EDs were all written as Homo sapiens,
                # which labelled the 582 Ludwig 2023 viral effectors (HHV 261,
                # vTR 186, CoV 135) as human. `organism` is what the public
                # `load(type="ED", organism=...)` filter reads, so the viral
                # sets carry an override; "viral" matches ED_VP16_TAD and
                # ED_RTA_TAD, the two records that were already right.
                "organism": info.get("organism_override", "Homo sapiens"),
                "uniprot_id": uniprot_from_col,
                "sequence_aa": seq,
                "length_aa": len(seq) if seq else None,
                "quantitative_metric": score,
                "quantitative_metric_label": col_score_used,
                "quantitative_metric_source": info["doi"],
                "validation_level": info.get("validation_level_override", "screen-validated"),
                "source": screen_name,
                "source_doi": info["doi"],
                "source_version": "supplementary_table",
                "source_date": now[:10],
                "chromatin_state_effect": _EFFECT_MAP.get(subtype, ""),
                "readout": info.get("readout"),
                "notes": f"From {screen_name}; {info['doi']}",
                "date_added": now,
                "date_modified": now,
            }
            records.append(record)

        log.info(f"[{screen_name}] loaded {len([r for r in records if screen_name in r['source']]):,} records")

    return records


# ── Main ──────────────────────────────────────────────────────────────────────

def _row_digest(rec: dict) -> str:
    """Six hex characters over the content that distinguishes two rows sharing
    a name: the sequence, and the metric label (DelRosso 2023 reports the same
    domain at two promoters, `Max avg pEF` and `Max avg PGK`). Deterministic,
    so a re-ingest produces the same id -- which an insertion-order suffix
    would not."""
    key = f"{rec.get('sequence_aa') or ''}|{rec.get('quantitative_metric_label') or ''}"
    return hashlib.sha256(key.encode()).hexdigest()[:6]


def resolve_ids(records: list[dict], log) -> list[dict]:
    """Make `module_id` unique, and link dual-readout siblings (F10, F8).

    `module_id` was f"ED_{screen}_{name[:40]}", and a library `name` is not
    unique within a screen: 725 ids carried more than one row, so 04's upsert
    overwrote 730 of them -- 713 with a sequence the survivor did not have.
    Truncation was never the cause (no colliding id had two distinct names);
    the rows are genuinely different tiles or conditions under one label.

    Exact duplicates collapse with a count (Ludwig 2023 ships 16 such rows).
    Anything still sharing an id gets `~<digest>` appended to EVERY member, so
    the id is a function of content rather than of which row came first."""
    by_id: dict[str, list[dict]] = {}
    for rec in records:
        by_id.setdefault(rec["module_id"], []).append(rec)

    kept, dropped_dup, rekeyed = [], 0, 0
    for mid, group in by_id.items():
        if len(group) == 1:
            kept.extend(group)
            continue
        seen, unique = set(), []
        for rec in group:
            fingerprint = (rec.get("sequence_aa"), rec.get("quantitative_metric_label"),
                           rec.get("readout"), rec.get("quantitative_metric"))
            if fingerprint in seen:
                dropped_dup += 1
                continue
            seen.add(fingerprint)
            unique.append(rec)
        if len(unique) == 1:
            kept.extend(unique)
            continue
        for rec in unique:
            rec["module_id"] = f"{mid}~{_row_digest(rec)}"
            rekeyed += 1
        kept.extend(unique)

    log.info(f"[ids] {len(by_id):,} distinct ids from {len(records):,} rows; "
             f"{dropped_dup:,} exact duplicates collapsed, {rekeyed:,} rows "
             f"re-keyed with a content digest")

    # A digest over sequence and metric label still leaves Ludwig 2023's
    # within-sheet contradictions colliding: BALF2_P03227_010 appears twice in
    # the HHV sheet, same 80 aa, once as Activation (3.49) and once as
    # Repression (1.93). Those are two directional calls on one tile, so the
    # subtype is the discriminator -- content, not row order. Dropping either
    # would be the silent loss this fix exists to stop.
    final: dict[str, dict] = {}
    remaining: dict[str, list[dict]] = {}
    for rec in kept:
        remaining.setdefault(rec["module_id"], []).append(rec)
    for mid, group in remaining.items():
        if len(group) == 1:
            final[mid] = group[0]
            continue
        # EVERY member takes the suffix, not just the ones that arrived later,
        # so the id stays a function of content rather than of row order.
        for rec in group:
            rec["module_id"] = f"{mid}_{(rec.get('subtype') or 'na')[:4].upper()}"
        regrouped: dict[str, list[dict]] = {}
        for rec in group:
            regrouped.setdefault(rec["module_id"], []).append(rec)
        for mid2, grp in regrouped.items():
            if len(grp) == 1:
                final[mid2] = grp[0]
                continue
            log.error(f"[ids] {len(grp)} rows share name, sequence, metric "
                      f"label and subtype under {mid2}; numbering them by "
                      f"score so the ids are at least stable")
            ordered = sorted(grp, key=lambda r: (str(r.get("quantitative_metric")),
                                                 str(r.get("sequence_aa"))))
            for i, rec in enumerate(ordered, start=1):
                rec["module_id"] = f"{mid2}_{i}"
                final[rec["module_id"]] = rec

    # F8: a sequence one sheet labels both ways is the same situation as a
    # dual-readout paper, so give those rows a `readout` too and let the
    # pairing below link them. The contradiction is then recorded rather than
    # resolved by picking a winner from a sheet that disagrees with itself.
    contradictory: dict[tuple, set] = {}
    for rec in final.values():
        if rec.get("sequence_aa"):
            contradictory.setdefault((rec.get("source"), rec.get("name"),
                                      rec["sequence_aa"]), set()).add(rec.get("subtype"))
    marked = 0
    for rec in final.values():
        if rec.get("readout") or not rec.get("sequence_aa"):
            continue
        subs = contradictory.get((rec.get("source"), rec.get("name"),
                                  rec["sequence_aa"]), set())
        if {"activator", "repressor"} <= subs:
            rec["readout"] = ("activation" if rec.get("subtype") == "activator"
                              else "repression")
            marked += 1
    if marked:
        log.info(f"[ids] {marked:,} rows carry a within-source directional "
                 f"contradiction; `readout` set from their own subtype so the "
                 f"two calls stay visible and linked")

    # F8: link the arms. Two rows are siblings when the same paper measured the
    # same sequence under the same name in a different arm.
    groups: dict[tuple, list[dict]] = {}
    for rec in final.values():
        if not rec.get("readout") or not rec.get("sequence_aa"):
            continue
        groups.setdefault((rec.get("source_doi"), rec.get("name"),
                           rec["sequence_aa"]), []).append(rec)
    paired = 0
    for group in groups.values():
        if len({r["readout"] for r in group}) < 2:
            continue
        for rec in group:
            siblings = [r["module_id"] for r in group
                        if r["module_id"] != rec["module_id"]]
            rec["paired_module_id"] = ";".join(sorted(siblings))
            paired += 1
    log.info(f"[ids] {paired:,} rows linked to a dual-readout sibling")
    return list(final.values())


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="config/sources.yaml")
    parser.add_argument("--dry-run", action="store_true",
                        help="Skip network calls; parse local files only")
    args = parser.parse_args()

    log = get_logger("02_seed_ed", LOGS_DIR)
    log.info("=== 02_seed_ed.py started ===")

    session = requests.Session()
    session.headers.update({"User-Agent": "module_library/1.0 (research)"})

    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Canonical curated EDs
    curated = process_curated_eds(MANUAL_DIR / "ed_curated.yaml", session, log, args.dry_run)
    log.info(f"Curated EDs: {len(curated)} entries")

    # 2. Screen data (skipped if files missing)
    screen = process_screen_data(log)
    log.info(f"Screen EDs: {len(screen)} entries")

    all_records = resolve_ids(curated + screen, log)
    df = pd.DataFrame(all_records)
    out = PROCESSED_DIR / "ed_raw.tsv"
    df.to_csv(out, sep="\t", index=False)
    log.info(f"Wrote {len(df):,} ED records to {out}")

    write_build_manifest(
        ROOT / "library",
        script_name="02_seed_ed.py",
        sources_used=[
            {"name": "ed_curated.yaml", "version": "manual_v1"},
            {"name": "Alerasool_2022", "doi": SCREEN_FILES["Alerasool_2022"]["doi"]},
            {"name": "DelRosso_2023",  "doi": SCREEN_FILES["DelRosso_2023_activators"]["doi"]},
        ],
        records_added=len(all_records),
    )

    log.info("=== 02_seed_ed.py complete ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
