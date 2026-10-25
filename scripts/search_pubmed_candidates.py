"""Log reproducible PubMed searches for candidate-negative review; never assigns labels."""
import argparse
import csv
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
STROKE = '("ischemic stroke"[Title/Abstract] OR "cerebral ischemia"[Title/Abstract] OR "cerebral infarction"[Title/Abstract] OR stroke[Title/Abstract] OR "Brain Ischemia"[MeSH Terms] OR "Stroke"[MeSH Terms])'


def search(query):
    url = BASE + "esearch.fcgi?" + urlencode({
        "db": "pubmed", "term": query, "retmax": "1000", "retmode": "xml", "tool": "StrokeProtAI"})
    request = Request(url, headers={"User-Agent": "StrokeProt-AI/1.0 (reproducible research)"})
    with urlopen(request, timeout=20) as response:
        result = ET.fromstring(response.read())
    return result.findtext("Count", "0"), [node.text or "" for node in result.findall(".//Id")]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="data/raw/negative_search_pool.csv")
    parser.add_argument("--output", default="data/raw/pubmed_negative_search_results.csv")
    parser.add_argument("--workers", type=int, default=2,
                        help="Concurrent requests (keep low for NCBI's unauthenticated rate limit).")
    args = parser.parse_args()
    with open(args.input, newline="", encoding="utf-8-sig") as handle:
        pool = list(csv.DictReader(handle))
    existing = {}
    if Path(args.output).exists():
        with open(args.output, newline="", encoding="utf-8-sig") as handle:
            existing = {row["gene_symbol"]: row for row in csv.DictReader(handle)}
    pending = []
    for row in pool:
        if row["gene_symbol"] in existing and existing[row["gene_symbol"]]["search_status"] == "completed":
            continue
        terms = [term.strip() for term in row["protein_aliases"].split("|") if term.strip()]
        protein = " OR ".join(f'"{term}"[Title/Abstract]' for term in terms)
        query = f"({protein}) AND {STROKE}"
        pending.append((row, query))

    def run(item):
        row, query = item
        # Retry transient HTTP/network failures once, without turning a failure
        # into an apparent zero-hit result.
        error_text = ""
        for attempt in range(2):
            try:
                count, pmids = search(query)
                return {
                    "gene_symbol": row["gene_symbol"], "aliases": row["protein_aliases"],
                    "study_univariate_result": row["study_univariate_result"],
                    "pubmed_query": query, "search_date": date.today().isoformat(),
                    "pubmed_result_count": count, "pmids_returned_first_1000": ";".join(pmids),
                    "search_status": "completed", "error": "",
                    "curation_status": "unreviewed_search_results",
                    "review_note": "Search hits require human review; zero hits alone do not assign a label."}
            except Exception as error:
                error_text = str(error)
                time.sleep(1.0 + attempt)
        return {
            "gene_symbol": row["gene_symbol"], "aliases": row["protein_aliases"],
            "study_univariate_result": row["study_univariate_result"],
            "pubmed_query": query, "search_date": date.today().isoformat(),
            "pubmed_result_count": "", "pmids_returned_first_1000": "",
            "search_status": "search_failed", "error": error_text,
            "curation_status": "unreviewed_search_failed",
            "review_note": "Search failed; retry required. Failure is not a zero-hit result."}

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    fields = ["gene_symbol", "aliases", "study_univariate_result", "pubmed_query", "search_date",
              "pubmed_result_count", "pmids_returned_first_1000", "search_status", "error",
              "curation_status", "review_note"]
    results = list(existing.values())
    completed = 0
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool_exec:
        futures = [pool_exec.submit(run, item) for item in pending]
        for future in as_completed(futures):
            record = future.result()
            results = [r for r in results if r["gene_symbol"] != record["gene_symbol"]]
            results.append(record)
            with open(args.output, "w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(sorted(results, key=lambda r: r["gene_symbol"].casefold()))
            completed += 1
            print(f"{record['gene_symbol']}: {record['pubmed_result_count'] or 'not searched'} ({record['search_status']})", flush=True)
            # Unauthenticated NCBI E-utilities limit is 3 requests/second.
            time.sleep(0.36)
    print(f"Logged {len(results)} queries to {args.output}; no labels were assigned.")


if __name__ == "__main__":
    main()
