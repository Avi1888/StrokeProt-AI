"""Command-line entry point for Model 1."""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from strokeprot.model1.data import load_dataset
from strokeprot.model1.train import run_model1


def main() -> None:
    parser = argparse.ArgumentParser(description="Train Model 1: AAC + sequence length Random Forest")
    parser.add_argument("--data", type=Path, required=True, help="CSV containing uniprot_id, sequence, label, and preferably cluster_id")
    parser.add_argument("--output", type=Path, default=Path("results/model1"), help="Directory for model, metrics, predictions, and figures")
    parser.add_argument("--random-pilot", action="store_true", help="Allow a preliminary stratified random split when no usable cluster_id is available")
    parser.add_argument("--seed", type=int, default=42, help="Reproducible random seed")
    args = parser.parse_args()
    data = load_dataset(args.data)
    metrics = run_model1(data, args.output, random_pilot=args.random_pilot, seed=args.seed)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()

