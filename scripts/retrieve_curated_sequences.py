"""Fetch human RefSeq sequences for positive evidence rows from NCBI."""
import argparse
import csv
import time
import xml.etree.ElementTree as ET
from datetime import date
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from Bio import SeqIO
from io import StringIO

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
AA = set("ACDEFGHIKLMNPQRSTVWY")


def ncbi(endpoint, params):
    params = {**params, "tool": "StrokeProtAI"}
    request = Request(BASE + endpoint + "?" + urlencode(params),
                      headers={"User-Agent": "StrokeProt-AI/1.0 (reproducible research)"})
    with urlopen(request, timeout=12) as response:
        return response.read()


def protein_candidates(gene_id):
    result = ET.fromstring(ncbi("esearch.fcgi", {
        "db": "protein", "term": f"{gene_id}[Gene ID] AND srcdb_refseq[PROP]", "retmax": "100"}))
    ids = [node.text for node in result.findall(".//Id")]
    if not ids:
        raise ValueError(f"No RefSeq proteins returned for NCBI Gene {gene_id}")
    summary = ET.fromstring(ncbi("esummary.fcgi", {"db": "protein", "id": ",".join(ids)}))
    rows = []
    for doc in summary.findall("DocSum"):
        items = {item.attrib.get("Name"): item.text or "" for item in doc.findall("Item")}
        accession, title = items.get("AccessionVersion", ""), items.get("Title", "")
        if accession.startswith("NP_") and "[Homo sapiens]" in title:
            rows.append((accession, title, int(items.get("Length", "0"))))
    if not rows:
        raise ValueError(f"No curated NP_ RefSeq protein returned for Gene {gene_id}")
    return rows


def retrieve(accession):
    payload = ncbi("efetch.fcgi", {"db": "protein", "id": accession,
                                   "rettype": "fasta", "retmode": "text"}).decode("utf-8")
    records = list(SeqIO.parse(StringIO(payload), "fasta"))
    if len(records) != 1:
        raise ValueError(f"NCBI returned {len(records)} records for {accession}")
    record = records[0]
    sequence = str(record.seq).upper()
    returned = record.id.split("|")[-2] if "|" in record.id else record.id
    if returned != accession or "[Homo sapiens]" not in record.description:
        raise ValueError(f"Record mismatch or non-human sequence for {accession}: {record.description}")
    if not sequence or set(sequence) - AA:
        raise ValueError(f"{accession} is empty or contains non-standard amino-acid symbols")
    return record, sequence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--curation", default="data/raw/evidence_curation.csv")
    parser.add_argument("--out", default="data/processed/proteins_esm.csv")
    parser.add_argument("--fasta", default="data/processed/proteins_esm.fasta")
    args = parser.parse_args()
    with open(args.curation, newline="", encoding="utf-8-sig") as handle:
        curated = list(csv.DictReader(handle))
    existing_by_gene = {}
    if Path(args.out).exists():
        with open(args.out, newline="", encoding="utf-8-sig") as handle:
            existing_by_gene = {row["gene_symbol"]: row for row in csv.DictReader(handle)}
    rows, fasta, seen_accessions, seen_sequences = [], [], set(), set()
    verified_metadata = {}
    for evidence in curated:
        if evidence["status"] not in {"positive", "candidate_negative"}:
            continue
        label = 1 if evidence["status"] == "positive" else 0
        gene_id = evidence["ncbi_gene_id"].strip()
        previous = existing_by_gene.get(evidence["gene_symbol"])
        if previous and previous.get("organism") == "Homo sapiens" and previous.get("verification_status", "").startswith("NCBI RefSeq"):
            accession = previous["protein_accession"]
            sequence = previous["sequence"].upper()
            if accession.startswith("NP_") and sequence and not (set(sequence) - AA):
                record = None
                candidates = [(
                    accession, f"{previous['protein_name']} [Homo sapiens]", len(sequence))]
                row = dict(previous)
                row["label"] = label
                row["label_status"] = evidence["status"]
                row["evidence_category"] = evidence["evidence_category"]
            else:
                previous = None
        if not previous:
            try:
                candidates = protein_candidates(gene_id)
                requested = evidence["refseq_accession"].strip()
                options = [item for item in candidates if item[0] == requested]
                # If no accession was supplied, use RefSeq isoform 1; otherwise the
                # longest curated NP_ record. Record the policy and remaining ambiguity.
                if not options:
                    options = [item for item in candidates if "isoform 1 " in item[1].lower()]
                chosen = max(options or candidates, key=lambda item: item[2])
                accession, title, _ = chosen
                record, sequence = retrieve(accession)
                row = {
                    "protein_accession": accession, "gene_symbol": evidence["gene_symbol"],
                    "protein_name": evidence["protein_name"], "ncbi_gene_id": gene_id,
                    "sequence": sequence, "sequence_length": len(sequence),
                    "organism": "Homo sapiens", "verification_status": "NCBI RefSeq; Homo sapiens verified",
                    "retrieval_date": date.today().isoformat(),
                    "sequence_source_url": f"https://www.ncbi.nlm.nih.gov/protein/{accession}",
                    "accession_mapping_uncertainty": "", "label": label,
                    "label_status": evidence["status"], "evidence_category": evidence["evidence_category"],
                    "cluster_id": ""}
            except Exception as error:
                evidence["sequence_status"] = "retrieval_failed"
                evidence["curator_notes"] = f"Sequence retrieval failed; retry when NCBI is reachable. Error: {error}"
                print(f"{evidence['gene_symbol']}: retrieval failed ({error})")
                continue
        if accession in seen_accessions:
            raise ValueError(f"Duplicate RefSeq accession encountered: {accession}")
        if sequence in seen_sequences:
            raise ValueError(f"Identical sequence encountered for {accession}; review the gene/isoform mapping")
        seen_accessions.add(accession)
        seen_sequences.add(sequence)
        if len(candidates) > 1 and not row.get("accession_mapping_uncertainty"):
            row["accession_mapping_uncertainty"] = "One RefSeq protein isoform selected (isoform 1 when available, else longest); alternatives not modeled."
        rows.append(row)
        verified_metadata[evidence["gene_symbol"]] = {
            "refseq_accession": accession,
            "sequence_status": "retrieved_verified_human_refseq",
            "sequence_length": len(sequence), "organism": "Homo sapiens",
            "retrieval_date": date.today().isoformat(),
            "verification_status": "NCBI RefSeq; Homo sapiens verified"}
        fasta.append((accession, sequence))
        if record is None:
            print(f"{evidence['gene_symbol']}: reused validated local sequence {accession}; {len(sequence)} aa")
        else:
            print(f"{evidence['gene_symbol']}: {accession}; {len(sequence)} aa; {record.description}")
            time.sleep(0.36)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else [
            "protein_accession", "gene_symbol", "protein_name", "ncbi_gene_id", "sequence",
            "sequence_length", "organism", "verification_status", "retrieval_date",
            "sequence_source_url", "accession_mapping_uncertainty", "label", "label_status",
            "evidence_category", "cluster_id"])
        writer.writeheader()
        writer.writerows(rows)
    Path(args.fasta).parent.mkdir(parents=True, exist_ok=True)
    with open(args.fasta, "w", encoding="utf-8") as handle:
        for accession, sequence in fasta:
            handle.write(f">{accession} [Homo sapiens]\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + "\n")
    for row in curated:
        if row["gene_symbol"] in verified_metadata:
            row.update(verified_metadata[row["gene_symbol"]])
    curation_fields = list(curated[0])
    for field in ("sequence_length", "organism", "retrieval_date", "verification_status"):
        if field not in curation_fields:
            curation_fields.append(field)
    with open(args.curation, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=curation_fields)
        writer.writeheader()
        writer.writerows(curated)
    from collections import Counter
    print(f"Retrieved {len(rows)} evidence-curated sequences: {dict(Counter(row['label_status'] for row in rows))}.")


if __name__ == "__main__":
    main()
