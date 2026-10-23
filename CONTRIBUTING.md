# Contributing

## Two-contributor workflow

1. Start from the latest `main` and create a branch for one task, such as `data/curated-labels` or `model/model-1-random-forest`.
2. Keep changes focused and commit them with a short descriptive message.
3. Open a pull request to `main` and explain the scope, data provenance, and how the change was checked.
4. The other contributor reviews the change, especially label rules, split leakage, and whether reported results were measured.
5. Resolve comments and conflicts on the feature branch, then merge the pull request. Do not force-push over another contributor's work.

Both contributors should agree on the stroke subtype definition, association source and release, inclusion/evidence threshold, candidate-negative interpretation, and final test split before model comparisons. Never use the test set for feature choices, hyperparameter selection, threshold selection, or repeated model iteration.

## Data handling

- Keep raw source downloads read-only and document query details, versions, and retrieval dates.
- Preserve stable identifiers and source evidence alongside labels.
- Candidate negatives mean “no qualifying record found in the selected source under these rules,” not confirmed biological negatives.
- Do not commit large raw downloads, trained checkpoints, secrets, or patient-level data without checking permission and licensing.
- Do not report unmeasured or placeholder metrics as results.
