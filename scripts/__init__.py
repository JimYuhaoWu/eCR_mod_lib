"""
ecr_mod_lib.scripts
-------------------
Exposes a zero-setup load() function for users who want to query the
pre-built library without running the pipeline from scratch.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd


# Repo root is three levels up from this file (scripts/__init__.py → scripts/ → repo root)
_REPO_ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_TSV = _REPO_ROOT / "library" / "module_library.tsv"


def load(
    type: Optional[str] = None,
    organism: Optional[str] = None,
    tsv_path: Optional[str | Path] = None,
) -> pd.DataFrame:
    """
    Load the pre-built module library as a pandas DataFrame.

    Parameters
    ----------
    type : str, optional
        Filter by module type: 'DBD', 'ED', or 'CR'.
    organism : str, optional
        Filter by organism, e.g. 'Homo sapiens' or 'Mus musculus'.
    tsv_path : str or Path, optional
        Path to a custom TSV file. Defaults to library/module_library.tsv
        in the repository root.

    Returns
    -------
    pd.DataFrame

    Examples
    --------
    >>> from scripts import load
    >>> df = load()                              # all ~10,800 records
    >>> eds = load(type="ED")                   # effector domains only
    >>> human_dbds = load(type="DBD", organism="Homo sapiens")
    """
    path = Path(tsv_path) if tsv_path else _DEFAULT_TSV

    if not path.exists():
        raise FileNotFoundError(
            f"Library TSV not found at: {path}\n"
            "To build the library from scratch, run the pipeline:\n"
            "  python scripts/01_fetch_dbd.py\n"
            "  python scripts/02_seed_ed.py\n"
            "  python scripts/03_fetch_cr.py\n"
            "  python scripts/04_build_library.py\n"
            "  python scripts/05_validate.py"
        )

    df = pd.read_csv(path, sep="\t", low_memory=False)

    if type is not None:
        type_upper = type.upper()
        if type_upper not in {"DBD", "ED", "CR"}:
            raise ValueError(f"type must be 'DBD', 'ED', or 'CR'; got '{type}'")
        df = df[df["type"] == type_upper]

    if organism is not None:
        df = df[df["organism"] == organism]

    return df.reset_index(drop=True)
