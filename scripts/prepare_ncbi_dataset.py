"""Fetch curated human protein accessions from NCBI; labels must be supplied with evidence."""
import argparse
import csv
import time
from io import StringIO
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from Bio import SeqIO

REQUIRED = {"protein_accession", "label", "gene_symbol", "label_source", "evidence_id", "evidence_url"}
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fetch(accession: str):
    query = urlencode({"db": "protein", "id": accession, "rettype": "fasta", "retmode": "text"})
    request = Request(f"{BASE}?{query}", headers={"User-Agent": "StrokeProt-AI/0.1 (academic reproducibility)"})
    with urlopen(request, timeout=45) as response:
        records = list(SeqIO.parse(StringIO(response.read().decode("utf-8")), "fasta"))
    if len(records) != 1:
        raise ValueError(f"NCBI returned {len(records)} FASTA records for {accession}; provide an exact protein accession.")
    return records[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--curation", default="data/raw/curated_labels.csv")
    parser.add_argument("--out", default="data/processed/proteins_esm.csv")
    parser.add_argument("--fasta", default="data/processed/proteins_esm.fasta")
    parser.add_argument("--xlsx", default="data/processed/proteins_esm.xlsx")
    args = parser.parse_args()
    with open(args.curation, newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or REQUIRED - set(rows[0]):
        raise SystemExit(f"Curation CSV needs rows and columns: {', '.join(sorted(REQUIRED))}")
    output_rows, fasta_records, seen = [], [], set()
    for row in rows:
        accession = row["protein_accession"].strip()
        if not accession or accession in seen:
            raise SystemExit(f"Empty or duplicate protein accession: {accession!r}")
        seen.add(accession)
        try:
            label = int(row["label"])
        except ValueError:
            raise SystemExit(f"Label must be 0 or 1 for {accession}")
        if label not in (0, 1) or not all(row[key].strip() for key in ("label_source", "evidence_id", "evidence_url")):
            raise SystemExit(f"{accession}: use binary label and cite its label source, evidence ID and URL")
        if label == 0 and "candidate" not in row["label_source"].lower():
            raise SystemExit(f"{accession}: label 0 must be explicitly described as a candidate negative")
        record = fetch(accession)
        sequence = str(record.seq).upper()
        if "[Homo sapiens]" not in record.description:
            raise SystemExit(f"{accession}: NCBI header does not confirm Homo sapiens: {record.description}")
        if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
            raise SystemExit(f"{accession}: returned sequence is empty or includes non-standard residues; review manually")
        header_accession = record.id.split("|")[-2] if "|" in record.id else record.id
        if not header_accession.startswith(accession.split(".")[0]):
            raise SystemExit(f"Requested {accession}, NCBI returned {record.id}; verify the accession mapping")
        output_rows.append({"protein_accession": record.id, "gene_symbol": row["gene_symbol"].strip(),
                            "sequence": sequence, "length": len(sequence), "organism": "Homo sapiens",
                            "label": label, "label_type": "literature-supported association" if label == 1 else "candidate negative; source-limited",
                            "label_source": row["label_source"].strip(), "evidence_id": row["evidence_id"].strip(),
                            "evidence_url": row["evidence_url"].strip(), "ncbi_fasta_id": record.id,
                            "sequence_source": "NCBI Protein EFetch", "cluster_id": ""})
        fasta_records.append((record.id, sequence))
        time.sleep(0.34)
    if {row["label"] for row in output_rows} != {0, 1}:
        raise SystemExit("Dataset must contain evidence-curated examples of both labels; no negatives are auto-generated.")
    out_path, fasta_path = Path(args.out), Path(args.fasta)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    Path(args.xlsx).parent.mkdir(parents=True, exist_ok=True)
    columns = list(output_rows[0])
    with open(out_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(output_rows)
    with open(fasta_path, "w", encoding="utf-8") as handle:
        for accession, sequence in fasta_records:
            handle.write(f">{accession}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + "\n")
    import pandas as pd
    pd.DataFrame(output_rows).to_excel(args.xlsx, index=False)
    print(f"Wrote {len(output_rows)} evidence-linked NCBI proteins to CSV, XLSX and FASTA.")
    print("cluster_id is intentionally blank: perform sequence clustering before a final evaluation.")


if __name__ == "__main__":
    main()
