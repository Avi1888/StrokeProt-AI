"""Train and evaluate the AAC + length Random Forest baseline."""

import json
import time
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    average_precision_score,
    confusion_matrix,
    precision_recall_curve,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
)

from .features import extract_features
from .splitting import split_indices

PARAMETER_GRID = [
    {"n_estimators": trees, "min_samples_leaf": leaf, "max_features": max_features}
    for trees in (300, 400)
    for leaf in (1, 2, 4)
    for max_features in ("sqrt", 0.5)
]


def _best_f1_threshold(labels: np.ndarray, scores: np.ndarray) -> float:
    precision, recall, thresholds = precision_recall_curve(labels, scores)
    if len(thresholds) == 0:
        return 0.5
    f1 = 2 * precision[:-1] * recall[:-1] / np.maximum(precision[:-1] + recall[:-1], 1e-12)
    best = np.flatnonzero(f1 == np.nanmax(f1))
    # On ties, choose the more conservative (higher) score threshold.
    return float(thresholds[best[-1]])


def _specificity(labels: np.ndarray, predictions: np.ndarray) -> float:
    tn, fp, _, _ = confusion_matrix(labels, predictions, labels=[0, 1]).ravel()
    return float(tn / (tn + fp)) if tn + fp else float("nan")


def _save_plots(data: pd.DataFrame, features: pd.DataFrame, y_test, scores, predictions, model, out: Path) -> None:
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    counts = data["label"].value_counts().reindex([0, 1], fill_value=0)
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(["Candidate negative", "Positive association"], counts.values, color=["#8097ad", "#d66b4d"])
    ax.set_ylabel("Protein count")
    ax.set_title("Class distribution")
    fig.tight_layout()
    fig.savefig(figures / "class_distribution.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    for label, name, color in [(0, "Candidate negative", "#8097ad"), (1, "Positive association", "#d66b4d")]:
        ax.hist(features.loc[data["label"] == label, "sequence_length"], bins=25, alpha=0.65, label=name, color=color)
    ax.set_xlabel("Sequence length (amino acids)")
    ax.set_ylabel("Protein count")
    ax.set_title("Sequence length by label")
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures / "sequence_length_by_label.png", dpi=160)
    plt.close(fig)

    cm = confusion_matrix(y_test, predictions, labels=[0, 1])
    fig, ax = plt.subplots(figsize=(4.5, 4))
    image = ax.imshow(cm, cmap="Blues")
    ax.set_xticks([0, 1], labels=["Candidate negative", "Positive"])
    ax.set_yticks([0, 1], labels=["Candidate negative", "Positive"])
    ax.set_xlabel("Predicted label")
    ax.set_ylabel("True label")
    ax.set_title("Held-out test confusion matrix")
    for (row, col), value in np.ndenumerate(cm):
        ax.text(col, row, str(value), ha="center", va="center")
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(figures / "confusion_matrix.png", dpi=160)
    plt.close(fig)

    precision, recall, _ = precision_recall_curve(y_test, scores)
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(recall, precision, color="#315f87")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Held-out test precision-recall curve")
    fig.tight_layout()
    fig.savefig(figures / "precision_recall_curve.png", dpi=160)
    plt.close(fig)

    importance = pd.DataFrame({"feature": features.columns, "importance": model.feature_importances_})
    importance = importance.sort_values("importance", ascending=False)
    importance.to_csv(out / "feature_importance.csv", index=False)
    top = importance.head(15).sort_values("importance")
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.barh(top["feature"], top["importance"], color="#4d83aa")
    ax.set_xlabel("Random Forest feature importance")
    ax.set_title("Most important predictive features")
    fig.tight_layout()
    fig.savefig(figures / "feature_importance.png", dpi=160)
    plt.close(fig)


def run_model1(data: pd.DataFrame, output_dir: str | Path, random_pilot: bool = False, seed: int = 42) -> dict:
    """Tune on validation data, then report once on the held-out test partition."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    features = extract_features(data)
    y = data["label"].to_numpy(dtype=int)
    partitions, split_type = split_indices(data, random_pilot=random_pilot, seed=seed)
    train_idx, val_idx, test_idx = (partitions[k] for k in ("train", "validation", "test"))

    best_model, best_ap, best_params = None, -np.inf, None
    tuning_rows = []
    start = time.perf_counter()
    for params in PARAMETER_GRID:
        model = RandomForestClassifier(class_weight="balanced", random_state=seed, n_jobs=-1, **params)
        model.fit(features.iloc[train_idx], y[train_idx])
        val_scores = model.predict_proba(features.iloc[val_idx])[:, 1]
        val_ap = float(average_precision_score(y[val_idx], val_scores))
        tuning_rows.append({**params, "validation_pr_auc": val_ap})
        if val_ap > best_ap:
            best_model, best_ap, best_params = model, val_ap, params
    train_seconds = time.perf_counter() - start
    pd.DataFrame(tuning_rows).to_csv(output / "validation_tuning.csv", index=False)

    val_scores = best_model.predict_proba(features.iloc[val_idx])[:, 1]
    threshold = _best_f1_threshold(y[val_idx], val_scores)
    test_scores = best_model.predict_proba(features.iloc[test_idx])[:, 1]
    test_predictions = (test_scores >= threshold).astype(int)

    test_labels = y[test_idx]
    metrics = {
        "model": "Model 1: AAC + sequence length + Random Forest",
        "split_type": split_type,
        "split_counts": {name: int(len(idx)) for name, idx in partitions.items()},
        "split_class_counts": {
            name: {str(label): int((y[idx] == label).sum()) for label in (0, 1)}
            for name, idx in partitions.items()
        },
        "seed": seed,
        "selected_parameters": best_params,
        "validation_pr_auc_for_selected_parameters": best_ap,
        "decision_threshold_selected_on_validation_f1": threshold,
        "test_metrics": {
            "pr_auc_average_precision": float(average_precision_score(test_labels, test_scores)),
            "roc_auc": float(roc_auc_score(test_labels, test_scores)),
            "precision": float(precision_score(test_labels, test_predictions, zero_division=0)),
            "recall": float(recall_score(test_labels, test_predictions, zero_division=0)),
            "f1": float(f1_score(test_labels, test_predictions, zero_division=0)),
            "specificity": _specificity(test_labels, test_predictions),
            "confusion_matrix_tn_fp_fn_tp": confusion_matrix(test_labels, test_predictions, labels=[0, 1]).ravel().tolist(),
            "test_protein_count": int(len(test_idx)),
        },
        "training_and_tuning_seconds": train_seconds,
        "interpretation_note": "Feature importance and predictions describe patterns associated with the dataset labels; they do not establish causality or clinical validity.",
    }
    (output / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    pred = data.iloc[test_idx][[c for c in ("uniprot_id", "gene_symbol", "label", "cluster_id") if c in data]].copy()
    pred["score"] = test_scores
    pred["predicted_label"] = test_predictions
    pred["threshold"] = threshold
    pred.to_csv(output / "test_predictions.csv", index=False)
    joblib.dump(best_model, output / "random_forest.joblib")
    _save_plots(data, features, test_labels, test_scores, test_predictions, best_model, output)
    return metrics

