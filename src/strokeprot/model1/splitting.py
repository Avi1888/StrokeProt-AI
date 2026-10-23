"""Reproducible train/validation/test splits with optional cluster grouping."""

import warnings

import numpy as np
from sklearn.model_selection import GroupShuffleSplit, train_test_split


def _has_both_classes(labels: np.ndarray, indices: np.ndarray) -> bool:
    return len(np.unique(labels[indices])) == 2


def split_indices(data, random_pilot: bool = False, seed: int = 42):
    """Return 70/15/15 indices; use similarity clusters unless pilot is requested."""
    labels = data["label"].to_numpy(dtype=int)
    all_indices = np.arange(len(data))
    if random_pilot:
        warnings.warn(
            "Using a row-random pilot split. Homologous proteins may cross partitions; do not treat this as final evaluation.",
            RuntimeWarning,
        )
        train_val, test = train_test_split(
            all_indices, test_size=0.15, random_state=seed, stratify=labels
        )
        train, validation = train_test_split(
            train_val, test_size=0.15 / 0.85, random_state=seed, stratify=labels[train_val]
        )
        split_type = "stratified_random_pilot"
    else:
        if "cluster_id" not in data:
            raise ValueError("Leakage-aware splitting requires cluster_id. Add sequence similarity clusters or explicitly pass --random-pilot for a preliminary run.")
        groups = data["cluster_id"].astype(str).to_numpy()
        found = None
        for offset in range(300):
            outer = GroupShuffleSplit(n_splits=1, test_size=0.15, random_state=seed + offset)
            train_val, test = next(outer.split(all_indices, labels, groups))
            if not _has_both_classes(labels, train_val) or not _has_both_classes(labels, test):
                continue
            inner = GroupShuffleSplit(n_splits=1, test_size=0.15 / 0.85, random_state=seed + 1000 + offset)
            train_rel, val_rel = next(inner.split(train_val, labels[train_val], groups[train_val]))
            train, validation = train_val[train_rel], train_val[val_rel]
            if _has_both_classes(labels, train) and _has_both_classes(labels, validation):
                found = (train, validation, test)
                break
        if found is None:
            raise ValueError("Could not create three cluster-disjoint partitions containing both labels. Review cluster counts/class distribution or use --random-pilot for a clearly preliminary split.")
        train, validation, test = found
        split_type = "cluster_grouped"

    partitions = {"train": train, "validation": validation, "test": test}
    for name, indices in partitions.items():
        if not _has_both_classes(labels, indices):
            raise ValueError(f"The {name} split does not contain both labels.")
    if not random_pilot:
        group_sets = {name: set(data.iloc[idx]["cluster_id"].astype(str)) for name, idx in partitions.items()}
        names = list(group_sets)
        for i, left in enumerate(names):
            for right in names[i + 1:]:
                if group_sets[left] & group_sets[right]:
                    raise RuntimeError(f"Cluster leakage detected between {left} and {right}.")
    return partitions, split_type

