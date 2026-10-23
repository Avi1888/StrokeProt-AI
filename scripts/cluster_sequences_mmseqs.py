"""Create leakage-control sequence clusters with MMseqs2 and add cluster_id to CSV."""
import argparse
import csv
import shutil
import subprocess
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="data/processed/proteins_esm.csv")
    parser.add_argument("--min-seq-id", type=float, default=0.30)
    parser.add_argument("--coverage", type=float, default=0.80)
    args = parser.parse_args()
    executable = shutil.which("mmseqs")
    if not executable:
        raise SystemExit("MMseqs2 is not installed or not on PATH. Install it before final sequence-aware evaluation.")
    path = Path(args.data)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
        columns = list(rows[0]) if rows else []
    if not rows or not {"protein_accession", "sequence"}.issubset(columns):
        raise SystemExit("Input must contain protein_accession and sequence columns.")
    with tempfile.TemporaryDirectory(prefix="strokeprot_mmseqs_") as temporary:
        temp = Path(temporary)
        fasta, prefix, tmp_dir = temp / "seqs.fasta", temp / "cluster", temp / "tmp"
        with fasta.open("w", encoding="utf-8") as handle:
            for row in rows:
                handle.write(f">{row['protein_accession']}\n{row['sequence']}\n")
        command = [executable, "easy-cluster", str(fasta), str(prefix), str(tmp_dir),
                   "--min-seq-id", str(args.min_seq_id), "-c", str(args.coverage), "--cov-mode", "0"]
        subprocess.run(command, check=True)
        mapping = {}
        with open(f"{prefix}_cluster.tsv", encoding="utf-8") as handle:
            for line in handle:
                representative, member = line.rstrip("\n").split("\t")[:2]
                mapping[member] = representative
    missing = [row["protein_accession"] for row in rows if row["protein_accession"] not in mapping]
    if missing:
        raise SystemExit(f"MMseqs2 did not assign {len(missing)} proteins; refusing to create incomplete clusters.")
    for row in rows:
        row["cluster_id"] = mapping[row["protein_accession"]]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns if "cluster_id" in columns else columns + ["cluster_id"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Added MMseqs2 clusters to {path} (identity >= {args.min_seq_id}, coverage >= {args.coverage}).")


if __name__ == "__main__":
    main()
