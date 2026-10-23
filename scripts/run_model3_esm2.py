"""Run the frozen ESM-2 Model 3 pipeline."""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from strokeprot.model3.esm import embed_dataset, fit_evaluate, read_proteins


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="data/processed/proteins_esm.csv")
    parser.add_argument("--output", default="results/model3_esm2")
    parser.add_argument("--cache", default="data/processed/esm2_t6_embeddings.npz")
    parser.add_argument("--device", default="auto", help="auto, cpu, or a PyTorch device such as cuda")
    parser.add_argument("--random-pilot", action="store_true", help="allow a preliminary row-stratified split without sequence clusters")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    data = read_proteins(args.data)
    embeddings = embed_dataset(data, args.cache, args.device)
    metrics = fit_evaluate(data, embeddings, args.output, args.random_pilot, args.seed)
    print(f"Saved ESM-2 evaluation to {args.output}; held-out F1={metrics['f1']:.3f}")


if __name__ == "__main__":
    main()
