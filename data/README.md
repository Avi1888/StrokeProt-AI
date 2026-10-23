# Dataset location

Place the cleaned, documented protein table at `processed/proteins.csv`. See the repository README for the required columns and label rules. Do not put raw downloads or unreviewed labels in this directory. The model refuses all-positive input such as a small biomarker-candidate list; valid supervised training needs sourced examples from both classes.

For ESM-2, `raw/curated_labels_template.csv` shows the evidence fields required by `scripts/prepare_ncbi_dataset.py`. The script fetches exact NCBI accessions and writes CSV, XLSX, and FASTA; it does not assign labels or generate negatives. See `../docs/model3_esm2_dataset.md`.

