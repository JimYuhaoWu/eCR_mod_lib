"""
scripts/utils.py
----------------
Shared helpers: logging, checksums, download manifests, build manifests.
"""

import hashlib
import json
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path


# ── Logging ───────────────────────────────────────────────────────────────────

def get_logger(name: str, logs_dir: Path) -> logging.Logger:
    logs_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = logs_dir / f"{name}_{ts}.log"

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")

    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


# ── Checksums ─────────────────────────────────────────────────────────────────

def checksum_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


# ── Download manifest ─────────────────────────────────────────────────────────

def load_manifest(raw_dir: Path) -> dict:
    path = raw_dir / "download_manifest.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return {}


def save_manifest(manifest: dict, raw_dir: Path) -> None:
    raw_dir.mkdir(parents=True, exist_ok=True)
    path = raw_dir / "download_manifest.json"
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def record_download(manifest: dict, key: str, url: str,
                    file_path: Path, version: str, description: str) -> dict:
    manifest[key] = {
        "url": url,
        "version": version,
        "description": description,
        "sha256": checksum_file(file_path),
        "downloaded_at": datetime.now(timezone.utc).isoformat(),
        "local_path": str(file_path),
    }
    return manifest


# ── Build manifest ────────────────────────────────────────────────────────────

def write_build_manifest(library_dir: Path, script_name: str,
                         sources_used: list[dict], records_added: int,
                         notes: str = "") -> None:
    library_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = library_dir / "build_manifest.json"

    if manifest_path.exists():
        data = json.loads(manifest_path.read_text(encoding="utf-8"))
    else:
        data = {"runs": []}

    data["runs"].append({
        "script": script_name,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "records_added": records_added,
        "sources": sources_used,
        "notes": notes,
    })

    manifest_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
