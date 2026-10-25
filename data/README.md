# Dataset location

Place the cleaned, documented protein table at `processed/proteins.csv`. See the repository README for the required columns and label rules. Do not put raw downloads or unreviewed labels in this directory. The model refuses all-positive input such as a small biomarker-candidate list; valid supervised training needs sourced examples from both classes.

For ESM-2, `raw/evidence_curation.csv` is the source-reviewed evidence table. `scripts/retrieve_curated_sequences.py` fetches verified human RefSeq sequences only for reviewed positive and candidate-negative rows and writes `processed/proteins_esm.csv` and `processed/proteins_esm.fasta`; it does not assign labels. The current table has 56 positives, 10 provisional candidate negatives and 2 unresolved cases; unresolved sequences are excluded. Search logs and source supplements in `raw/` support the candidate-negative review. See `../docs/evidence_and_labels.md`, `../docs/candidate_negative_search.md`, and `../docs/dataset_manifest.md`.

