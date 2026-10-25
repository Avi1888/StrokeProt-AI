# StrokeProt-AI

A reproducible, evidence-curated research project asking whether one **human protein sequence** meets the project's documented definition of ischemic-stroke association. It is a protein-level research classification task. It does not diagnose a person, predict an individual's stroke, or prove that a protein causes stroke. The contribution is the evidence-curated question and dataset, not a claim that the model architecture is new.

## Current dataset status

The targeted review contains 56 positive candidates, 10 provisional `candidate_negative` candidates, and 2 unresolved candidates (F8 and CST5). All 66 trainable candidates have verified human RefSeq sequences. This permits an exploratory model run, but the evaluation set is too small for a reliable performance claim. The [seven-slide project overview](presentations/StrokeProt-AI_project_overview.pptx) predates this dataset expansion; use this README and manifest for current results. See the [presentation-ready project brief](docs/project_brief.md), [decisions and label rules](docs/evidence_and_labels.md), [candidate-negative search status](docs/candidate_negative_search.md), [evidence table](data/raw/evidence_curation.csv), and [dataset manifest](docs/dataset_manifest.md).

## Label definition

The project defines `positive` from primary human research directly relating a protein to ischemic stroke through one or more separately recorded categories: `risk_before_stroke`, `post_event_measurement`, or `severity_or_outcome`. Reviews, animal-only studies, database entries, candidate lists, and expression changes without qualifying human stroke evidence are not sufficient. `candidate_negative` means no qualifying evidence was found under a documented search; it is not a confirmed biological negative. `unresolved` cases stay in the curation table and are excluded from training. Full rules and the candidate-specific decisions are in [docs/evidence_and_labels.md](docs/evidence_and_labels.md).

## Local demo UI

Run `.venv\Scripts\python.exe app.py` and open `http://127.0.0.1:8000` to explore a pasted sequence or FASTA record. The demo returns a score from the trained exploratory model and prominently explains its limits. See [docs/demo_ui.md](docs/demo_ui.md).

## Primary model and current result

Use the existing Model 3 pipeline: frozen `esm2_t6_8M_UR50D`, mean-pooled per-residue embeddings, then regularized logistic regression. The encoder is not fine-tuned. This small pretrained model is used as a sequence representation; the project contribution remains evidence curation and evaluation, not architecture novelty. Model 1 is retained as an older amino-acid-composition/Random Forest baseline and is not the proposed primary model.

The method does not use patient clinical features, AlphaFold structure features, SMOTE, CNNs, ensembles, or a neural-network classification head. Plain-language method details and limits are in [docs/model3_esm2_dataset.md](docs/model3_esm2_dataset.md).

The seed-42 grouped split had train 46 (40 positive/6 candidate-negative), validation 9 (7/2), and test 11 (9/2), with 8 positive and 2 candidate-negative groups in test. Test precision was 0.80, recall 0.89, F1 0.84, balanced accuracy 0.44, PR-AUC 0.72, ROC-AUC 0.17; confusion matrix `[[0,2],[1,8]]`. The model missed both test candidate-negatives. Despite the F1, this is not evidence of useful discrimination: metrics are extremely uncertain and outputs are for software demonstration only. Full details are in the [manifest](docs/dataset_manifest.md).

## Evidence-first workflow

1. Review and extend `data/raw/evidence_curation.csv`. Record a primary paper title, DOI or PMID, evidence category, human study population, and short rationale for every decision. Do not assign candidate negatives without a logged, comparable literature search.
2. Run `python scripts/retrieve_curated_sequences.py` to fetch the positive candidates' human RefSeq sequences by NCBI Gene ID, validate them, and update the curation metadata. It does not create labels or negatives. Network access to NCBI is required.
3. The sequence-backed table is `data/processed/proteins_esm.csv`; verified sequences are also written to `data/processed/proteins_esm.fasta`. It contains 56 positive and 10 provisional candidate-negative rows; unresolved proteins are excluded from training.
4. Run `python scripts/check_dataset_readiness.py`. Similarity groups were computed locally with the documented pairwise-alignment fallback after the Windows MMseqs2 binary failed; 51 positive and 10 candidate-negative groups were found. The group floor is only a coarse split feasibility check.
5. Run the frozen model workflow:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python scripts/run_model3_esm2.py --data data/processed/proteins_esm.csv --output results/model3_esm2
```

The Model 3 run is saved locally under ignored `results/model3_esm2/`. Current seed-42 grouped split: train 46 (40 positive/6 candidate-negative), validation 9 (7/2), test 11 (9/2); test groups are 8 positive and 2 candidate-negative. Metrics: precision 0.80, recall 0.89, F1 0.84, balanced accuracy 0.44, PR-AUC 0.72, confusion matrix `[[0,2],[1,8]]`. The model missed both test candidate-negatives and ROC-AUC was 0.17. Despite the F1, this is not evidence of useful discrimination: with only two provisional negatives, results are highly uncertain and likely unrepresentative. Do not use training accuracy or the demo score as evidence of performance.

## Existing Model 1 baseline

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

## Sequence retrieval dependencies

The NCBI retrieval helper writes sequence CSV, XLSX, and FASTA for evidence-reviewed positive and candidate-negative rows when network access is available; it does not create labels. The first ESM-2 run downloads pretrained weights via Fair-ESM. Weights, embedding caches, and trained model outputs are not committed.

## Team workflow

Keep `main` as the shared integration branch. Each of the two contributors should create a focused branch from the latest `main`, commit only their work, and open a pull request back to `main`. The other contributor reviews the pull request; merge only after conflicts and data/provenance questions are resolved. Do not both commit directly to `main` or use one shared working branch.

Suggested ownership: one contributor owns dataset construction, evidence, and the dataset manifest; the other owns model implementation and evaluation. Agree on the label rules and final split before comparing models. See [CONTRIBUTING.md](CONTRIBUTING.md).

## Limitations

The evidence review is targeted, not systematic. Candidate negatives are provisional absence-of-evidence labels, not true negatives. The test set contains only two candidate-negative groups and both were classified as positive; ROC-AUC was 0.17. Any metric is extremely uncertain and does not establish generalization, causality, clinical usefulness, or person-level risk. The model is suitable only as an explicitly caveated software demonstration.
