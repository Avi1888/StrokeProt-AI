"""Input loading and validation for Model 1."""

from pathlib import Path

import pandas as pd

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
REQUIRED_COLUMNS = {"uniprot_id", "sequence", "label"}


def load_dataset(path: str | Path) -> pd.DataFrame:
    """Load a clean protein table and fail early on data-quality problems."""
    data = pd.read_csv(path, dtype={"uniprot_id": "string", "sequence": "string"})
    missing_columns = REQUIRED_COLUMNS - set(data.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")
    if data.empty:
        raise ValueError("Input dataset has no rows.")
    for column in REQUIRED_COLUMNS:
        if data[column].isna().any():
            raise ValueError(f"Column {column!r} contains missing values.")

    data["uniprot_id"] = data["uniprot_id"].astype(str).str.strip()
    data["sequence"] = data["sequence"].astype(str).str.upper().str.replace(r"\s+", "", regex=True)
    if (data["uniprot_id"] == "").any():
        raise ValueError("Protein accessions cannot be empty.")
    if data["uniprot_id"].duplicated().any():
        ids = data.loc[data["uniprot_id"].duplicated(keep=False), "uniprot_id"].tolist()[:10]
        raise ValueError(f"Duplicate UniProt accessions found (examples: {ids}).")

    invalid = data["sequence"].map(lambda seq: sorted(set(seq) - STANDARD_AA))
    invalid_rows = invalid[invalid.map(bool)]
    if invalid_rows.any():
        examples = [(str(data.loc[i, "uniprot_id"]), invalid.loc[i]) for i in invalid_rows.index[:5]]
        raise ValueError(f"Sequences contain symbols outside the 20 standard amino acids: {examples}")
    if data["sequence"].duplicated().any():
        ids = data.loc[data["sequence"].duplicated(keep=False), "uniprot_id"].astype(str).tolist()[:10]
        raise ValueError(f"Duplicate amino-acid sequences must be resolved before splitting (examples: {ids}).")

    numeric_labels = pd.to_numeric(data["label"], errors="coerce")
    if numeric_labels.isna().any() or not set(numeric_labels.unique()).issubset({0, 1}):
        raise ValueError("Labels must be binary integers 0 or 1.")
    data["label"] = numeric_labels.astype("int8")
    if data["label"].nunique() != 2:
        raise ValueError("Supervised binary training requires both labels 0 and 1; no labels will be fabricated.")
    if "organism" in data and not data["organism"].astype(str).str.contains("Homo sapiens", case=False).all():
        raise ValueError("All rows must be human (Homo sapiens).")
    if "cluster_id" in data:
        if data["cluster_id"].isna().any() or (data["cluster_id"].astype(str).str.strip() == "").any():
            raise ValueError("cluster_id is present but contains missing or empty values.")
    return data.reset_index(drop=True)

