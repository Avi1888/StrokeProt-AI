# Model 3: frozen ESM-2 and dataset protocol

## Model definition

Model 3 uses the existing pretrained ESM-2 checkpoint `esm2_t6_8M_UR50D`. The encoder is frozen. The mean of final-layer (layer 6) residue representations becomes one 320-dimensional vector per protein. A class-balanced Logistic Regression head is trained on those vectors. The downloaded checkpoint is an existing general protein language model; the project contribution is this downstream application and its evaluation, not a new pretrained model.

The embedding cache records a SHA-256 digest of the ordered sequences, checkpoint, layer, and pooling policy. Proteins longer than 1,022 residues are processed as overlapping windows (128-residue overlap), with window means combined by residue count. The script never silently truncates long proteins.

## Label and sequence data

`data/raw/curated_labels_template.csv` is a header-only template. Populate it with exact NCBI Protein accessions and labels reviewed by the team. Each row must include gene symbol, label source, evidence ID, and a durable evidence URL. Use `label=1` only for evidence that meets the team's predefined ischemic-stroke association rule. A `label=0` is a **candidate negative**, meaning only that the agreed search did not find qualifying evidence; it is not a biologically confirmed negative. Include `candidate` in the label-source text for label 0. Exclude uncertain records rather than guessing.

NCBI retrieves sequence records by the exact accession provided. The output includes accession, FASTA, label provenance, and source metadata. Review NCBI's human organism annotation and accession mapping before accepting the records. This script does not choose positives, infer negatives, or treat missing evidence as a negative. No training-ready CSV is included until both labels have defensible, source-linked records.

```powershell
python -m pip install -r requirements.txt
# Fill data/raw/curated_labels_template.csv and save it as curated_labels.csv
python scripts/prepare_ncbi_dataset.py --curation data/raw/curated_labels.csv
# Install MMseqs2 separately, then create clusters for leakage-aware splits
python scripts/cluster_sequences_mmseqs.py --data data/processed/proteins_esm.csv
```

Review each accession and label source as a team before training. Pin the disease definition, source release/date, query/filter rules, and manual curation decisions in the dataset manifest. Keep sensitive or very large raw downloads outside Git.

## Run Model 3

```powershell
python scripts/run_model3_esm2.py --data data/processed/proteins_esm.csv --output results/model3_esm2
```

By default, the existing Model 1 split helper makes a 70/15/15 cluster-grouped train/validation/test split (seed 42), rejects class-deficient partitions, and checks that clusters do not cross partitions. Validation F1 selects the classification threshold; the test partition is used once for reporting. `--random-pilot` enables a warned-about stratified random pilot if clusters cannot be generated. That pilot is not a leakage-controlled result.

Saved outputs: cached NumPy embeddings, fitted classifier, test-set predictions, and metrics JSON. Report test-set counts alongside precision, recall, F1, balanced accuracy, confusion matrix, ROC-AUC, and average precision. Metrics are not evidence of causal biology or clinical performance.

## Dataset sources

- NCBI Protein E-utilities EFetch: https://www.ncbi.nlm.nih.gov/books/NBK25500/
- Open Targets explains that disease-to-target lists can include indirect ontology-descendant evidence, so do not copy a broad list into positives without filtering and review: https://platform-docs.opentargets.org/associations
- Open Targets recommends bulk data downloads for systematic queries: https://platform-docs.opentargets.org/data-access

NCBI sequences alone do not establish stroke association. The evidence sheet remains the source of the labels.
