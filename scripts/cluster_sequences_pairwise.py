"""Cluster proteins by pairwise local-alignment identity (portable fallback)."""
import argparse
import csv
from pathlib import Path

from Bio.Align import PairwiseAligner


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="data/processed/proteins_esm.csv")
    parser.add_argument("--min-seq-id", type=float, default=0.30)
    parser.add_argument("--coverage", type=float, default=0.80)
    args = parser.parse_args()
    path = Path(args.data)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
        fields = list(rows[0]) if rows else []
    if not rows:
        raise SystemExit("No sequence rows found.")
    if len({r["protein_accession"] for r in rows}) != len(rows):
        raise SystemExit("Duplicate protein accessions; refusing to cluster.")
    sequences = [r["sequence"].upper() for r in rows]
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -8
    aligner.extend_gap_score = -0.5
    parent = list(range(len(rows)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        a, b = find(i), find(j)
        if a != b:
            parent[b] = a

    links = []
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            aln = aligner.align(sequences[i], sequences[j])[0]
            coords = aln.coordinates
            matches = aligned1 = aligned2 = columns = 0
            for k in range(coords.shape[1] - 1):
                a0, a1 = int(coords[0, k]), int(coords[0, k + 1])
                b0, b1 = int(coords[1, k]), int(coords[1, k + 1])
                da, db = a1 - a0, b1 - b0
                if da and db:
                    n = min(da, db)
                    matches += sum(sequences[i][a0 + t] == sequences[j][b0 + t] for t in range(n))
                    aligned1 += da
                    aligned2 += db
                    columns += max(da, db)
                else:
                    aligned1 += da
                    aligned2 += db
                    columns += da + db
            identity = matches / columns if columns else 0.0
            cov1 = aligned1 / len(sequences[i])
            cov2 = aligned2 / len(sequences[j])
            if identity >= args.min_seq_id and min(cov1, cov2) >= args.coverage:
                union(i, j)
                links.append((rows[i]["protein_accession"], rows[j]["protein_accession"], identity, cov1, cov2))

    roots = {}
    for i in range(len(rows)):
        root = find(i)
        roots.setdefault(root, rows[root]["protein_accession"])
    for i, row in enumerate(rows):
        row["cluster_id"] = roots[find(i)]
    out_fields = fields if "cluster_id" in fields else fields + ["cluster_id"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote pairwise-alignment clusters to {path}: {len(set(r['cluster_id'] for r in rows))} groups; {len(links)} qualifying pair links.")
    print(f"Thresholds: identity >= {args.min_seq_id:.0%}; aligned coverage >= {args.coverage:.0%} on both proteins; connected components keep transitive groups together.")


if __name__ == "__main__":
    main()
