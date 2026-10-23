"""Sequence-derived features for the Model 1 baseline."""

from collections import Counter

import numpy as np
import pandas as pd

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
FEATURE_NAMES = [f"aac_{aa}" for aa in AMINO_ACIDS] + ["sequence_length"]


def amino_acid_composition(sequence: str) -> np.ndarray:
    """Return frequencies of the 20 standard residues in a validated sequence."""
    seq = sequence.upper()
    if not seq:
        raise ValueError("Protein sequence is empty after validation.")
    counts = Counter(seq)
    return np.asarray([counts[aa] / len(seq) for aa in AMINO_ACIDS], dtype=np.float32)


def extract_features(data: pd.DataFrame) -> pd.DataFrame:
    """Build AAC fractions plus sequence length, in a stable column order."""
    rows = []
    for sequence in data["sequence"]:
        composition = amino_acid_composition(sequence)
        rows.append(np.concatenate((composition, np.asarray([len(sequence)], dtype=np.float32))))
    return pd.DataFrame(rows, columns=FEATURE_NAMES, index=data.index)

