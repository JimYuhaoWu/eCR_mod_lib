"""
Merge Staller 2022 mmc4 (UniProt IDs) + mmc5 (activities) into one clean file.
Run once from project root: python scripts/prepare_staller2022.py
"""
import pandas as pd
from pathlib import Path

MANUAL = Path("data/manual")

mmc4 = pd.read_csv(MANUAL / "Staller_2022_mmc4.csv")
mmc5 = pd.read_csv(MANUAL / "Staller_2022_mmc5.csv")

# Extract UniProt ID and gene symbol from mmc4 GeneName (">sp|Q8N587|ZN561_HUMAN")
mmc4["uniprot_id"] = mmc4["uniprotID"]
mmc4["gene_symbol"] = mmc4["GeneName"].str.extract(r"\|([^|]+)_HUMAN")
mmc4_slim = mmc4[["ProteinRegionSeq", "uniprot_id", "gene_symbol", "Start", "End"]]

# Keep only Prediction rows from mmc5
preds = mmc5[mmc5["RegionType"] == "Prediction"].copy()
merged = preds.merge(mmc4_slim, on="ProteinRegionSeq", how="left")

out = MANUAL / "Staller_2022_predictions.csv"
merged.to_csv(out, index=False)

n_uniprot = merged["uniprot_id"].notna().sum()
n_pos = (merged["Activity_Zscore_mean"] > 0.5).sum()
print(f"Prediction rows: {len(merged)}")
print(f"With UniProt ID: {n_uniprot}")
print(f"Z-score > 0.5:   {n_pos}")
print(f"Saved: {out}")
