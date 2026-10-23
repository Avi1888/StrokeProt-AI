# StrokeProt-AI

Starter repository for sequence-based classification of human proteins by curated stroke-association evidence. Model 1 is a traditional baseline using amino-acid composition (AAC), sequence length, and a Random Forest.

This is a research prototype, not a clinical tool. A positive label means the chosen disease-association source contains a qualifying record under documented rules. A candidate negative means that the same source has no qualifying record under those rules; it does not mean the protein has no biological role in stroke.

## Model 1

Model 1 extracts the fraction of each of the 20 standard amino acids and the sequence length, then fits a `RandomForestClassifier` with class balancing. It selects hyperparameters by validation average precision and selects the decision threshold by validation F1. The test set is evaluated once after those choices are fixed.

### Input data

Put a cleaned CSV at `data/processed/proteins.csv` with these required columns:

| Column | Meaning |
|---|---|
| `uniprot_id` | Stable UniProt protein accession; one canonical sequence per accession |
| `sequence` | Amino-acid sequence using the 20 standard amino-acid symbols |
| `label` | `1` for qualifying stroke association, `0` for candidate negative under the same source/rules |
| `cluster_id` | Optional but strongly recommended sequence-similarity cluster for leakage-aware splitting |

An optional `gene_symbol` column is retained only for outputs. The loader rejects missing values, duplicate accessions or sequences, non-standard sequence symbols, invalid labels, and datasets with only one class. No sample dataset, negative labels, or performance values are fabricated or included. The seven all-positive biomarker candidates from the separate starter package are not suitable for supervised binary training.

Record the disease scope, association database and release, query/filters, evidence rules, download date, accession mapping, and sequence source using the [dataset manifest template](docs/dataset_manifest_template.md) before training. Keep raw downloads out of Git unless licensing and repository size allow them.

### Install and run

Use Python 3.10 or 3.11:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python scripts/run_model1.py --data data/processed/proteins.csv --output results/model1
```

The default is a 70/15/15 train/validation/test split grouped by `cluster_id`. If clusters are unavailable, a stratified random split is available only as an explicitly marked preliminary pilot:

```powershell
python scripts/run_model1.py --data data/processed/proteins.csv --output results/model1 --random-pilot
```

The split code checks that every partition contains both classes and, for cluster-based splitting, that clusters do not cross partitions. It stops with an explanation if the data are too small or too few clusters are available. It does not silently fall back to row-random splitting.

Outputs include the fitted model, held-out test predictions, JSON metrics, feature-importance CSV, class-count and sequence-length plots, confusion matrix, and precision-recall curve. Model outputs are ignored by Git by default.

## Model 3: ESM-2 transfer learning

Model 3 uses the frozen pretrained `esm2_t6_8M_UR50D` encoder, mean-pooled final-layer residue embeddings, and a class-balanced Logistic Regression head. It uses the same cluster-aware split helper and seed as Model 1 when both models receive the same rows and `cluster_id` values. See [the Model 3 dataset and method protocol](docs/model3_esm2_dataset.md).

Prepare a source-linked curation CSV from `data/raw/curated_labels_template.csv`, then retrieve its exact human protein accessions from NCBI:

```powershell
python scripts/prepare_ncbi_dataset.py --curation data/raw/curated_labels.csv
python scripts/cluster_sequences_mmseqs.py --data data/processed/proteins_esm.csv
python scripts/run_model3_esm2.py --data data/processed/proteins_esm.csv --output results/model3_esm2
```

The NCBI script does not assign labels or auto-generate negatives. A candidate negative means no qualifying record was found under the specified evidence rule; it is not a confirmed biological negative. This repository intentionally does not include a training-ready labeled dataset or fabricated performance results. The first ESM-2 run downloads pretrained weights through Fair-ESM and may take time on CPU.

## Team workflow

Keep `main` as the shared integration branch. Each of the two contributors should create a focused branch from the latest `main`, commit only their work, and open a pull request back to `main`. The other contributor reviews the pull request; merge only after conflicts and data/provenance questions are resolved. Do not both commit directly to `main` or use one shared working branch.

Suggested ownership: one contributor owns dataset construction, evidence, and the dataset manifest; the other owns model implementation and evaluation. Agree on the label rules and final split before comparing models. See [CONTRIBUTING.md](CONTRIBUTING.md).

## Limitations

This baseline ignores residue order, and its learned feature importance is predictive association rather than biological causation. Sequence similarity and weak labels can inflate apparent performance. A high score would not establish causal biology or clinical usefulness. Report the test-set size and metrics honestly; do not fill results with expected values.
