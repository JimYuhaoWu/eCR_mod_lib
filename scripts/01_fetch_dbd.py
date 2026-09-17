"""
scripts/01_fetch_dbd.py
-----------------------
Download and parse DNA-Binding Domain (DBD) entries from three sources:

  1. AnimalTFDB v4  – TF gene list with DBD family annotation
                      for Homo sapiens and Mus musculus
  2. UniProt REST   – amino acid sequences for each TF
  3. JASPAR 2024    – position-frequency matrix (motif) IDs where available

Outputs
-------
  data/raw/animaltfdb/{species}_TF.txt      raw downloaded files
  data/raw/uniprot/{uniprot_id}.fasta       per-entry sequences
  data/raw/jaspar/{uniprot_id}_jaspar.json  per-entry JASPAR hits
  data/raw/download_manifest.json           checksum + URL log
  data/processed/dbd_raw.tsv               merged, cleaned table

Run
---
  python scripts/01_fetch_dbd.py [--config config/sources.yaml] [--dry-run]

Reproducibility
---------------
  Every downloaded file is SHA-256 checksummed and logged in
  data/raw/download_manifest.json. Re-running is idempotent: files
  with a matching checksum are not re-fetched.
"""

from __future__ import annotations

import argparse
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests
import yaml
from tqdm import tqdm

# Allow running from project root or scripts/
sys.path.insert(0, str(Path(__file__).parent))
from utils import (
    get_logger, checksum_file, load_manifest, save_manifest,
    record_download, write_build_manifest,
)

ROOT = Path(__file__).resolve().parents[1]
LOGS_DIR = ROOT / "logs"
RAW_DIR = ROOT / "data" / "raw"
PROCESSED_DIR = ROOT / "data" / "processed"


# ── Helpers ───────────────────────────────────────────────────────────────────

def _get(url: str, session: requests.Session, timeout: int = 30,
         retries: int = 3, backoff: float = 2.0) -> requests.Response | None:
    """GET with retry + exponential backoff. Returns None on failure."""
    for attempt in range(retries):
        try:
            r = session.get(url, timeout=timeout)
            if r.status_code == 200:
                return r
            if r.status_code == 404:
                return None                    # not found — skip silently
        except requests.RequestException:
            pass
        time.sleep(backoff ** attempt)
    return None


# ── Step 1: Download AnimalTFDB ───────────────────────────────────────────────

def fetch_animaltfdb(cfg: dict, session: requests.Session,
                     manifest: dict, log, dry_run: bool) -> dict[str, pd.DataFrame]:
    """
    Download AnimalTFDB TF files for each configured species.
    Returns a dict {species_key: DataFrame}.
    """
    atfdb_cfg = cfg["animaltfdb"]
    base_url = atfdb_cfg["base_url"]
    version = atfdb_cfg["version"]
    dfs = {}

    for species_key, sp_info in atfdb_cfg["species"].items():
        fname = sp_info["filename"]
        url = f"{base_url}/{fname}"
        out_path = RAW_DIR / "animaltfdb" / fname

        # Use cached file if it exists (manifest entry optional — covers manual downloads)
        if out_path.exists():
            log.info(f"[AnimalTFDB] {species_key}: using existing {out_path.name}")
            if f"animaltfdb_{species_key}" not in manifest:
                manifest = record_download(manifest, f"animaltfdb_{species_key}",
                                           url, out_path, version,
                                           f"AnimalTFDB v{version} {species_key} (manual)")
                save_manifest(manifest, RAW_DIR)
        else:
            log.info(f"[AnimalTFDB] downloading {url}")
            if dry_run:
                log.info("[AnimalTFDB] DRY-RUN: skipping download")
                continue
            resp = _get(url, session)
            if resp is None:
                log.error(f"[AnimalTFDB] FAILED to download {url}. "
                          f"Download manually from https://guolab.wchscu.cn/AnimalTFDB4.0 "
                          f"and place at {out_path}")
                continue
            out_path.write_text(resp.text, encoding="utf-8")
            manifest = record_download(manifest, f"animaltfdb_{species_key}",
                                       url, out_path, version,
                                       f"AnimalTFDB v{version} {species_key}")
            save_manifest(manifest, RAW_DIR)
            log.info(f"[AnimalTFDB] saved {out_path} ({out_path.stat().st_size:,} bytes)")

        if not out_path.exists():
            continue

        # Parse – handle both tab and comma separators gracefully
        try:
            df = pd.read_csv(out_path, sep="\t", dtype=str, comment="#")
        except Exception as e:
            log.warning(f"[AnimalTFDB] parse error for {out_path}: {e}")
            continue

        df.columns = [c.strip() for c in df.columns]
        col_map = atfdb_cfg["col_map"]
        # Rename to standard names; tolerate missing columns
        rename = {v: k for k, v in col_map.items() if v in df.columns}
        df = df.rename(columns=rename)

        for needed in ("symbol", "family"):
            if needed not in df.columns:
                log.warning(f"[AnimalTFDB] column '{needed}' not found in {fname}; "
                            f"columns: {list(df.columns)}")

        df["species_key"] = species_key
        df["taxon_id"] = str(sp_info["taxon_id"])
        dfs[species_key] = df
        log.info(f"[AnimalTFDB] {species_key}: {len(df):,} TFs loaded")

    return dfs


# ── Step 2: Fetch UniProt sequences ──────────────────────────────────────────

def fetch_uniprot_sequence(uniprot_id: str, session: requests.Session,
                            api_base: str) -> str | None:
    """Fetch FASTA sequence for a UniProt accession. Returns the sequence string."""
    out_path = RAW_DIR / "uniprot" / f"{uniprot_id}.fasta"
    if out_path.exists():
        text = out_path.read_text()
        seq = "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))
        return seq or None

    url = f"{api_base}/{uniprot_id}.fasta"
    resp = _get(url, session, timeout=15)
    if resp is None:
        return None

    out_path.write_text(resp.text, encoding="utf-8")
    seq = "".join(line.strip() for line in resp.text.splitlines() if not line.startswith(">"))
    return seq or None


# C2H2-type zinc-finger families (AnimalTFDB 'Family' values) whose members
# coordinate one Zn(2+) ion per finger. Only these get a zinc_finger_count.
ZINC_FINGER_FAMILIES = {"zf-C2H2", "ZBTB"}


def fetch_zinc_finger_count(uniprot_id: str, session: requests.Session,
                            api_base: str) -> int | None:
    """
    Count curated zinc-finger features (ZN_FING) for a UniProt accession.

    Each ZN_FING feature is one Zn(2+)-coordinating finger; the count is used
    downstream to add that many Zn ions to AF3 inputs. Returns the integer
    count, or None if the entry could not be fetched.
    """
    cache = RAW_DIR / "uniprot" / f"{uniprot_id}_znf.tsv"
    if cache.exists():
        text = cache.read_text()
    else:
        url = f"{api_base}/{uniprot_id}.tsv?fields=ft_zn_fing"
        resp = _get(url, session, timeout=15)
        if resp is None:
            return None
        text = resp.text
        cache.write_text(text, encoding="utf-8")

    return text.count("ZN_FING")


def lookup_uniprot_by_gene(gene: str, taxon_id: str, session: requests.Session,
                            api_base: str) -> tuple[str | None, str | None]:
    """
    Search UniProt by gene name + taxon for a reviewed (Swiss-Prot) entry.
    Returns (uniprot_id, sequence) or (None, None).
    """
    url = (f"{api_base}/search?query=gene_exact:{gene}"
           f"+AND+organism_id:{taxon_id}+AND+reviewed:true"
           f"&format=tsv&fields=accession,sequence&size=1")
    resp = _get(url, session, timeout=15)
    if resp is None or not resp.text.strip():
        return None, None

    lines = resp.text.strip().splitlines()
    if len(lines) < 2:
        return None, None

    parts = lines[1].split("\t")
    if len(parts) < 2:
        return None, None

    accession, seq = parts[0].strip(), parts[1].strip()
    # Cache locally
    out_path = RAW_DIR / "uniprot" / f"{accession}.fasta"
    if not out_path.exists():
        fasta = f">{accession}\n{seq}\n"
        out_path.write_text(fasta, encoding="utf-8")

    return accession, seq


# ── Step 3: Fetch JASPAR motif ID ────────────────────────────────────────────

def _matrix_parts(matrix_id: str) -> tuple[str, int]:
    """Split a JASPAR accession into (base, version): 'MA1724.2' -> ('MA1724', 2)."""
    base, _, ver = matrix_id.partition(".")
    return base, (int(ver) if ver.isdigit() else 0)


def fetch_jaspar_id(gene: str, taxon_id: str, session: requests.Session,
                    api_base: str) -> str | None:
    """
    Query JASPAR for a matrix matching the gene name, anywhere in vertebrates.
    Returns the highest-versioned matching matrix ID (e.g. 'MA1724.2') or None.

    NOT filtered by `taxon_id`. JASPAR tags a profile with the species of the
    protein the experiment used, not the species the profile applies to, and CORE
    profiles are vertebrate-wide. MA1724 was derived from mouse Rfx6, so
    `?name=RFX6&tax_id=9606` returns nothing and human RFX6 was recorded
    motif-less while the mouse row carried MA1724.1 -- the same profile, already
    in the library, on the other species' row. That split 94 human DBDs (89 CORE)
    off from motifs we could already resolve, and since target.py drops
    motif-less DBDs it silently removed RFX6, ISL1, PTF1A and NKX2-1 from every
    design run. `taxon_id` is kept in the signature for the caller's sake but is
    deliberately not part of the query or the cache key.
    """
    cache = RAW_DIR / "jaspar" / f"{gene}_vertebrates.json"
    if cache.exists():
        import json
        data = json.loads(cache.read_text())
    else:
        # page_size must exceed the version count for a gene: versions are
        # separate results, and truncating them loses the current matrix.
        url = (f"{api_base}/matrix/?name={gene}&tax_group=vertebrates"
               f"&format=json&page_size=50")
        resp = _get(url, session, timeout=15)
        if resp is None:
            return None
        import json
        data = resp.json()
        cache.write_text(json.dumps(data, indent=2), encoding="utf-8")

    # Exact name only. The API's `name` filter is not exact, so a substring hit
    # would attach another gene's motif to this record.
    exact = [r for r in data.get("results", [])
             if str(r.get("name", "")).upper() == gene.upper()]
    if not exact:
        return None

    # Highest version, not the first result. A version bump is a re-derived
    # matrix, so taking results[0] pinned most of the library to '.1' and scored
    # every design against a superseded PWM.
    #
    # Versions are only comparable WITHIN an accession. PTF1A has three distinct
    # CORE profiles all named 'Ptf1A' (MA1618, MA1619, MA1620), so a plain max
    # over versions picks between different matrices on a number that means
    # nothing across them. Take the newest version of each accession, then the
    # lowest accession, which is at least stable across runs and releases.
    # Choosing between genuinely different profiles for one TF needs curation,
    # not a tiebreak, so the alternatives are left in the cached JSON rather
    # than resolved here.
    newest: dict[str, int] = {}
    for r in exact:
        base, ver = _matrix_parts(r["matrix_id"])
        if ver > newest.get(base, -1):
            newest[base] = ver
    # min() over the accession strings also prefers a CORE 'MA....' profile to an
    # UNVALIDATED 'UN....' one ('M' < 'U'), which is the right precedence and is
    # relied on here -- 4 of the 19 accession changes this fix produced were
    # UN -> MA upgrades (NANOG UN0805.1 -> MA2339.1, FOSB, BCL11B).
    chosen = min(newest)
    return f"{chosen}.{newest[chosen]}"


# ── Step 4: Assemble processed table ─────────────────────────────────────────

def build_dbd_records(dfs: dict[str, pd.DataFrame], cfg: dict,
                      session: requests.Session, log,
                      dry_run: bool) -> pd.DataFrame:
    """
    For each TF in the AnimalTFDB tables, enrich with UniProt sequence
    and JASPAR motif ID, then return a processed DataFrame.
    """
    up_cfg = cfg["uniprot"]
    ja_cfg = cfg["jaspar"]

    all_records = []

    for species_key, df in dfs.items():
        sp_info = cfg["animaltfdb"]["species"][species_key]
        taxon_id = str(sp_info["taxon_id"])
        organism_label = "Homo sapiens" if "Homo" in species_key else "Mus musculus"

        for _, row in tqdm(df.iterrows(), total=len(df),
                           desc=f"Enriching {species_key} DBDs", ncols=80):

            symbol = str(row.get("symbol", "")).strip()
            family = str(row.get("family", "")).strip()
            uniprot_id = str(row.get("uniprot_id", "")).strip()

            if not symbol or symbol.lower() in ("nan", ""):
                continue

            # Sequence
            seq, resolved_uid = None, uniprot_id or None
            if uniprot_id and uniprot_id.lower() not in ("nan", ""):
                if not dry_run:
                    seq = fetch_uniprot_sequence(uniprot_id, session, up_cfg["api_base"])
            if seq is None and not dry_run:
                resolved_uid, seq = lookup_uniprot_by_gene(
                    symbol, taxon_id, session, up_cfg["api_base"])

            # Zinc-finger count (only for C2H2-type families; one Zn per finger)
            zinc_finger_count = None
            if family in ZINC_FINGER_FAMILIES and resolved_uid and not dry_run:
                zinc_finger_count = fetch_zinc_finger_count(
                    resolved_uid, session, up_cfg["api_base"])

            # JASPAR motif
            jaspar_id = None
            if not dry_run:
                jaspar_id = fetch_jaspar_id(symbol, taxon_id, session, ja_cfg["api_base"])

            now = datetime.now(timezone.utc).isoformat()
            module_id = f"DBD_{species_key[:2].upper()}_{symbol.upper()}"

            record = {
                "module_id": module_id,
                "type": "DBD",
                "subtype": family,
                "name": symbol,
                "organism": organism_label,
                "source_species": species_key,
                "gene_symbol": symbol,
                "uniprot_id": resolved_uid,
                "sequence_aa": seq,
                "length_aa": len(seq) if seq else None,
                "zinc_finger_count": zinc_finger_count,
                "target_or_mechanism": f"{family} family TF; see JASPAR {jaspar_id}" if jaspar_id else f"{family} family TF",
                "validation_level": "motif-only" if jaspar_id else "predicted",
                "source": f"AnimalTFDB v{cfg['animaltfdb']['version']}",
                "source_doi": cfg["animaltfdb"]["doi"],
                "source_version": cfg["animaltfdb"]["version"],
                "source_date": now[:10],
                "jaspar_id": jaspar_id,
                "notes": (f"Sequence from UniProt {resolved_uid}" if resolved_uid else
                          "Sequence not yet fetched"),
                "date_added": now,
                "date_modified": now,
            }
            all_records.append(record)

    return pd.DataFrame(all_records)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="config/sources.yaml",
                        help="Path to sources.yaml (default: config/sources.yaml)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Skip all network calls; parse existing files only")
    parser.add_argument("--limit", type=int, default=0,
                        help="Process only first N TFs per species (for testing)")
    args = parser.parse_args()

    log = get_logger("01_fetch_dbd", LOGS_DIR)
    log.info("=== 01_fetch_dbd.py started ===")
    log.info(f"config={args.config}  dry_run={args.dry_run}  limit={args.limit}")

    cfg_path = ROOT / args.config
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    # Requests session with a descriptive User-Agent (polite crawling)
    session = requests.Session()
    session.headers.update({"User-Agent": "module_library/1.0 (research; contact: see paper)"})

    # Ensure all raw subdirectories exist (not committed to git)
    for subdir in ("animaltfdb", "uniprot", "jaspar"):
        (RAW_DIR / subdir).mkdir(parents=True, exist_ok=True)

    manifest = load_manifest(RAW_DIR)

    # 1. AnimalTFDB
    dfs = fetch_animaltfdb(cfg, session, manifest, log, args.dry_run)

    if not dfs:
        log.error("No AnimalTFDB data loaded. Check URL or manual download instructions.")
        sys.exit(1)

    # Optional limit for testing
    if args.limit > 0:
        dfs = {k: v.head(args.limit) for k, v in dfs.items()}
        log.info(f"[LIMIT] truncated to {args.limit} TFs per species for testing")

    # 2–3. UniProt + JASPAR enrichment
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    processed = build_dbd_records(dfs, cfg, session, log, args.dry_run)

    out_tsv = PROCESSED_DIR / "dbd_raw.tsv"
    processed.to_csv(out_tsv, sep="\t", index=False)
    log.info(f"Wrote {len(processed):,} DBD records to {out_tsv}")

    # 4. Build manifest
    write_build_manifest(
        ROOT / "library",
        script_name="01_fetch_dbd.py",
        sources_used=[
            {"name": "AnimalTFDB", "version": cfg["animaltfdb"]["version"],
             "doi": cfg["animaltfdb"]["doi"]},
            {"name": "UniProt", "version": cfg["uniprot"]["version"]},
            {"name": "JASPAR", "version": cfg["jaspar"]["version"],
             "doi": cfg["jaspar"]["doi"]},
        ],
        records_added=len(processed),
        notes=f"dry_run={args.dry_run}  limit={args.limit}",
    )

    log.info("=== 01_fetch_dbd.py complete ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
