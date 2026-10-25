# Model 3: frozen ESM-2 and dataset protocol

## Model definition

Model 3 uses the existing pretrained ESM-2 checkpoint `esm2_t6_8M_UR50D`. The encoder is frozen. The mean of final-layer (layer 6) residue representations becomes one 320-dimensional vector per protein. A class-balanced Logistic Regression head is trained on those vectors. The downloaded checkpoint is an existing general protein language model; the project contribution is this downstream application and its evaluation, not a new pretrained model.

The embedding cache records a SHA-256 digest of the ordered sequences, checkpoint, layer, and pooling policy. Proteins longer than 1,022 residues are processed as overlapping windows (128-residue overlap), with window means combined by residue count. The script never silently truncates long proteins.

## Label and sequence data

The active evidence source is `data/raw/evidence_curation.csv`; the sequence-backed model table is `data/processed/proteins_esm.csv`. The older `curated_labels_template.csv` and `prepare_ncbi_dataset.py` describe the original pipeline and are not the current curation workflow. Use `label=1` only for primary-human evidence meeting the project's documented stroke-association rule. `label=0` means only that qualifying evidence was not found under the documented candidate search; it is not a biologically confirmed negative. Exclude unresolved records from training.

NCBI sequence records were retrieved for 66 reviewed candidates and checked for organism, accession, amino-acid alphabet, and duplicates. There are 56 positives, 10 provisional candidate negatives, and 2 unresolved evidence rows. The dataset supports only an exploratory model demonstration; the candidate-negative class and evaluation are too small for dependable scientific performance claims.

```powershell
python -m pip install -r requirements.txt
# Fill data/raw/curated_labels_template.csv and save it as curated_labels.csv
python scripts/prepare_ncbi_dataset.py --curation data/raw/curated_labels.csv
# Install MMseqs2 separately, then create clusters for leakage-aware splits
python scripts/cluster_sequences_mmseqs.py --data data/processed/proteins_esm.csv
```

Review each accession and label source independently before scientific reuse. Pin the disease definition, source release/date, query/filter rules, and manual decisions in the manifest. Keep pretrained weights, caches and model artifacts local and ignored by Git.

## Run Model 3

```powershell
python scripts/run_model3_esm2.py --data data/processed/proteins_esm.csv --output results/model3_esm2
```

By default, the existing Model 1 split helper makes a 70/15/15 cluster-grouped train/validation/test split (seed 42), rejects class-deficient partitions, and checks that clusters do not cross partitions. Validation F1 selects the classification threshold; the test partition is used once for reporting. `--random-pilot` enables a warned-about stratified random pilot if clusters cannot be generated. That pilot is not a leakage-controlled result.

Saved outputs: cached NumPy embeddings, fitted classifier, test-set predictions, and metrics JSON (ignored by Git). The recorded grouped split is train 46 (40/6 rows, 36/6 groups), validation 9 (7/2 rows, 7/2 groups), test 11 (9/2 rows, 8/2 groups). The model missed both test candidate-negative examples, with ROC-AUC 0.17. See `docs/dataset_manifest.md` for all metrics. This run demonstrates code execution only and is not evidence of dependable biological classification or clinical performance.

## Dataset sources

- NCBI Protein E-utilities EFetch: https://www.ncbi.nlm.nih.gov/books/NBK25500/
- Open Targets explains that disease-to-target lists can include indirect ontology-descendant evidence, so do not copy a broad list into positives without filtering and review: https://platform-docs.opentargets.org/associations
- Open Targets recommends bulk data downloads for systematic queries: https://platform-docs.opentargets.org/data-access

NCBI sequences alone do not establish stroke association. The evidence sheet remains the source of the labels.
