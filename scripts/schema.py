"""
scripts/schema.py
-----------------
Defines the SQLite schema for the module library and exposes a ModuleLibrary
context-manager class for thread-safe read/write access.

Schema overview
---------------
  modules      – one row per DBD / ED / CR entry
  provenance   – one row per script run that added/modified records

Usage
-----
    from schema import ModuleLibrary

    with ModuleLibrary("library/module_library.db") as lib:
        lib.insert_or_replace(records)   # list[dict]
        df = lib.to_dataframe()
"""

from __future__ import annotations

import sqlite3
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

import pandas as pd


# ── DDL ──────────────────────────────────────────────────────────────────────

CREATE_MODULES = """
CREATE TABLE IF NOT EXISTS modules (
    module_id               TEXT    PRIMARY KEY,
    type                    TEXT    NOT NULL
                                    CHECK(type IN ('DBD','ED','CR')),
    subtype                 TEXT,
    name                    TEXT    NOT NULL,
    organism                TEXT,
    source_species          TEXT,
    gene_symbol             TEXT,
    uniprot_id              TEXT,
    sequence_aa             TEXT,
    length_aa               INTEGER,
    zinc_finger_count       INTEGER,
    -- DBD extraction (06_extract_dbd.py). `sequence_aa` holds the WHOLE UniProt
    -- protein for type='DBD' rows, so a fusion built from it is a full-length-TF
    -- fusion carrying every native domain (KRAB, SCAN, TADs). These columns
    -- record the DNA-binding span without destroying the original: consumers
    -- migrate to `dbd_sequence_aa` deliberately rather than by a silent change
    -- of meaning under `sequence_aa`.
    sequence_full_aa        TEXT,
    dbd_start               INTEGER,
    dbd_end                 INTEGER,
    dbd_sequence_aa         TEXT,
    dbd_length_aa           INTEGER,
    dbd_extraction          TEXT,
    dbd_excluded_domains    TEXT,
    dbd_review_flag         INTEGER,
    -- F7: every residue outside the kept span, as 'start-end:label'; the flag
    -- is set when a dropped region reads as a self-association or repeat
    -- module. Separate from dbd_review_flag, which reports the opposite
    -- failure (an excluded domain INSIDE the span).
    dbd_dropped_regions     TEXT,
    dbd_dropped_review_flag INTEGER,
    target_or_mechanism     TEXT,
    quantitative_metric     REAL,
    quantitative_metric_label TEXT,
    quantitative_metric_source TEXT,
    validation_level        TEXT    NOT NULL
                                    CHECK(validation_level IN (
                                        'predicted',
                                        'motif-only',
                                        'ChIP-validated',
                                        'screen-validated',
                                        'structurally-resolved'
                                    )),
    source                  TEXT    NOT NULL,
    source_doi              TEXT,
    source_version          TEXT,
    source_date             TEXT,
    jaspar_id               TEXT,
    known_interactors       TEXT,
    linker_notes            TEXT,
    engineering_compatibility TEXT,
    chromatin_state_effect  TEXT,
    notes                   TEXT,
    date_added              TEXT    NOT NULL,
    date_modified           TEXT    NOT NULL
);
"""

CREATE_PROVENANCE = """
CREATE TABLE IF NOT EXISTS provenance (
    run_id          TEXT    PRIMARY KEY,
    script          TEXT    NOT NULL,
    run_date        TEXT    NOT NULL,
    source_name     TEXT    NOT NULL,
    source_url      TEXT,
    source_version  TEXT,
    file_path       TEXT,
    file_checksum   TEXT,
    records_added   INTEGER DEFAULT 0,
    records_updated INTEGER DEFAULT 0,
    notes           TEXT
);
"""

# Columns added after the table shipped. CREATE TABLE IF NOT EXISTS is a no-op on
# an existing database, so a new column in the DDL above would never reach one --
# every deployed library would silently lack it. Each entry is applied only when
# absent, which makes _create_schema() idempotent and safe to run repeatedly.
MIGRATIONS = [
    ("sequence_full_aa",     "TEXT"),
    ("dbd_start",            "INTEGER"),
    ("dbd_end",              "INTEGER"),
    ("dbd_sequence_aa",      "TEXT"),
    ("dbd_length_aa",        "INTEGER"),
    ("dbd_extraction",       "TEXT"),
    ("dbd_excluded_domains", "TEXT"),
    ("dbd_review_flag",      "INTEGER"),
    ("dbd_dropped_regions",  "TEXT"),
    ("dbd_dropped_review_flag", "INTEGER"),
]

CREATE_INDEXES = [
    "CREATE INDEX IF NOT EXISTS idx_modules_type ON modules(type);",
    "CREATE INDEX IF NOT EXISTS idx_modules_organism ON modules(organism);",
    "CREATE INDEX IF NOT EXISTS idx_modules_gene ON modules(gene_symbol);",
    "CREATE INDEX IF NOT EXISTS idx_modules_validation ON modules(validation_level);",
]


# ── ModuleLibrary class ───────────────────────────────────────────────────────

class ModuleLibrary:
    """
    Context manager wrapping a SQLite connection to the module library.

    Example
    -------
    with ModuleLibrary(db_path) as lib:
        lib.insert_or_replace(records)
        df = lib.to_dataframe("DBD")
    """

    REQUIRED_FIELDS = {"module_id", "type", "name", "validation_level", "source",
                       "date_added", "date_modified"}

    def __init__(self, db_path: str | Path):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self.conn: sqlite3.Connection | None = None

    def __enter__(self):
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        self.conn.execute("PRAGMA journal_mode=WAL;")
        self.conn.execute("PRAGMA foreign_keys=ON;")
        self._create_schema()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.conn:
            if exc_type is None:
                self.conn.commit()
            else:
                self.conn.rollback()
            self.conn.close()
        return False

    # ── schema ────────────────────────────────────────────────────────────────

    def _create_schema(self):
        cur = self.conn.cursor()
        cur.executescript(CREATE_MODULES + CREATE_PROVENANCE)
        cur.execute("PRAGMA table_info(modules);")
        have = {row["name"] for row in cur.fetchall()}
        for col, decl in MIGRATIONS:
            if col not in have:
                cur.execute(f"ALTER TABLE modules ADD COLUMN {col} {decl};")
        for idx in CREATE_INDEXES:
            cur.execute(idx)
        self.conn.commit()

    # ── write ─────────────────────────────────────────────────────────────────

    def insert_or_replace(self, records: list[dict]) -> tuple[int, int]:
        """
        Upsert a list of record dicts.  Returns (inserted, updated) counts.
        Missing optional fields are silently set to NULL.
        Raises ValueError if any required field is absent.
        """
        inserted, updated = 0, 0
        cur = self.conn.cursor()

        # Discover all column names from schema
        cur.execute("PRAGMA table_info(modules);")
        columns = [row["name"] for row in cur.fetchall()]

        for rec in records:
            missing = self.REQUIRED_FIELDS - rec.keys()
            if missing:
                raise ValueError(f"Record '{rec.get('module_id','?')}' missing: {missing}")

            # Check if already exists
            cur.execute("SELECT 1 FROM modules WHERE module_id=?", (rec["module_id"],))
            exists = cur.fetchone() is not None

            # Only the columns the record actually mentions. Filling every
            # schema column and INSERT-OR-REPLACEing would NULL any column the
            # caller does not know about: the dbd_* span columns exist only in
            # the database (no source TSV carries them), so one ordinary
            # 04_build_library.py run would silently erase the whole extraction
            # for all 3,248 DBD rows. An update touches what it was given.
            row = {col: rec[col] for col in columns if col in rec}
            cols = list(row.keys())
            if exists:
                assignments = ", ".join(f"{c}=?" for c in cols)
                cur.execute(f"UPDATE modules SET {assignments} WHERE module_id=?",
                            list(row.values()) + [rec["module_id"]])
            else:
                placeholders = ", ".join("?" * len(cols))
                cur.execute(
                    f"INSERT INTO modules ({', '.join(cols)}) VALUES ({placeholders})",
                    list(row.values()))

            if exists:
                updated += 1
            else:
                inserted += 1

        self.conn.commit()
        return inserted, updated

    def insert_provenance(self, record: dict) -> None:
        """Insert a provenance record."""
        cur = self.conn.cursor()
        cur.execute("PRAGMA table_info(provenance);")
        columns = [row["name"] for row in cur.fetchall()]
        row = {col: record.get(col) for col in columns}
        cols = list(row.keys())
        sql = (f"INSERT OR REPLACE INTO provenance ({', '.join(cols)}) "
               f"VALUES ({', '.join('?'*len(cols))})")
        cur.execute(sql, list(row.values()))
        self.conn.commit()

    # ── read ──────────────────────────────────────────────────────────────────

    def to_dataframe(self, module_type: str | None = None) -> pd.DataFrame:
        """Return all (or type-filtered) modules as a DataFrame."""
        if module_type:
            df = pd.read_sql_query(
                "SELECT * FROM modules WHERE type=?", self.conn,
                params=(module_type,))
        else:
            df = pd.read_sql_query("SELECT * FROM modules", self.conn)
        return df

    def counts(self) -> dict:
        """Return per-type and per-validation_level record counts."""
        cur = self.conn.cursor()
        cur.execute("SELECT type, COUNT(*) FROM modules GROUP BY type")
        by_type = dict(cur.fetchall())
        cur.execute("SELECT validation_level, COUNT(*) FROM modules GROUP BY validation_level")
        by_level = dict(cur.fetchall())
        return {"by_type": by_type, "by_validation_level": by_level}

    # ── export ────────────────────────────────────────────────────────────────

    def to_tsv(self, out_path: str | Path) -> None:
        """Export the full library to a TSV file (human-readable / git-trackable)."""
        df = self.to_dataframe()
        df.to_csv(out_path, sep="\t", index=False)

    # ── utils ─────────────────────────────────────────────────────────────────

    @contextmanager
    def transaction(self) -> Iterator[sqlite3.Cursor]:
        """Explicit transaction block for batched writes."""
        cur = self.conn.cursor()
        try:
            yield cur
            self.conn.commit()
        except Exception:
            self.conn.rollback()
            raise
