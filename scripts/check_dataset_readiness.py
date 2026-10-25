"""Summarize evidence curation and stop premature binary training."""
import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CURATION = ROOT / "data/raw/evidence_curation.csv"
DATASET = ROOT / "data/processed/proteins_esm.csv"


def main():
    with CURATION.open(newline="", encoding="utf-8-sig") as handle:
        curated = list(csv.DictReader(handle))
    with DATASET.open(newline="", encoding="utf-8-sig") as handle:
        records = list(csv.DictReader(handle))

    counts = Counter(row["status"] for row in curated)
    print(f"Reviewed evidence candidates: {len(curated)}")
    print(f"  positive: {counts['positive']}")
    print(f"  candidate_negative: {counts['candidate_negative']}")
    print(f"  unresolved: {counts['unresolved']}")
    print(f"Sequence records in training table: {len(records)}")
    sequence_counts = Counter(row.get("label_status", "") for row in records)
    print(f"Sequence-backed examples by status: {dict(sequence_counts)}")
    groups = Counter()
    for row in records:
        status = row.get("label_status", "")
        cluster = row.get("cluster_id", "")
        if status in {"positive", "candidate_negative"} and cluster:
            groups[(status, cluster)] += 1
    independent = Counter(status for status, _ in groups)
    print(f"Independent similarity groups by training status: {dict(independent)}")
    # At least six independent groups per class are needed even for a coarse
    # train/validation/test partition, and this is only a floor, not a guarantee.
    ready = all(independent[status] >= 6 for status in ("positive", "candidate_negative"))
    if ready:
        print("A minimal group-count check passed; inspect split feasibility and class balance before training.")
    else:
        print("dataset not ready for binary training")
        print("Do not create artificial negatives or report training accuracy as performance.")


if __name__ == "__main__":
    main()
