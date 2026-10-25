# StrokeProt-AI: project brief

## Research question

Does a human protein sequence belong to a protein with a documented primary-human-research association to ischemic stroke? The three separately tracked evidence types are pre-event risk, post-event measurement, and severity/outcome. This is a protein-level literature classification question. It does not diagnose people, predict an individual's stroke, or establish causation.

## What has been assembled

- 68 reviewed rows: 56 positive, 10 provisional candidate negatives, and 2 unresolved (F8 and CST5).
- The evidence includes a targeted review of the five carried-forward candidates and the 20 FDR-significant proteins reported in Supplementary Table S2 of a 534-case ischemic-stroke outcome cohort.
- All 66 labeled proteins have retrieved human RefSeq sequences that passed accession, organism, sequence alphabet, and duplicate checks.
- The evidence table records the primary paper title, DOI/PMID, human study population, evidence type, rationale, and sequence provenance/status.

## Model status

An exploratory binary classifier was trained using a grouped split. This run is suitable only to demonstrate the software: the test set had 9 positives and 2 provisional candidate negatives, and both negatives were misclassified. The metrics are highly uncertain and do not demonstrate useful classification.

If future evidence curation supports both classes with enough independent sequence groups, the planned primary method is frozen ESM-2 (`esm2_t6_8M_UR50D`), mean pooling across residues, and regularized logistic regression. It will not fine-tune the encoder or use clinical, structural, resampling, CNN, or ensemble features.

## Main limitations

The review is targeted rather than systematic. Candidate negatives mean no evidence was found within the documented search and are not biological negatives. The positive class mixes pre-event risk evidence with acute post-event and outcome associations. Only two candidate-negative groups occur in test, and both were missed. These limitations prevent a reliable performance claim.

## Recommended next work

1. Independently re-review the 10 provisional candidate-negative decisions and expand primary human literature searches; keep CST5 unresolved unless evidence supports adjudication.
2. Add a larger and more representative set of reviewed proteins and independent sequence groups, then rerun grouped evaluation. The present result is not a dependable predictor.
3. Add uncertainty/calibration analysis only when the test set is large enough to support it.

## Reproduction commands

```powershell
python scripts/check_dataset_readiness.py
python scripts/check_dataset_readiness.py
python scripts/run_model3_esm2.py --data data/processed/proteins_esm.csv --output results/model3_esm2
python app.py
```

The first model run downloads public ESM-2 weights; local trained outputs and caches are ignored by Git. See [the manifest](dataset_manifest.md), [evidence decisions](evidence_and_labels.md), and [README](../README.md).

## Demo interface

A local exploratory scoring UI is available with `python app.py`; it accepts one sequence and returns the model score and threshold-based class, with a warning about the model's test-set limitations. See [demo instructions](demo_ui.md).
