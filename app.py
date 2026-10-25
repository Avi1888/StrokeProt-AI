"""Local research demo for the exploratory StrokeProt-AI classifier."""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import Counter, defaultdict
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse

ROOT = Path(__file__).resolve().parent
os.environ.setdefault("TORCH_HOME", str(ROOT / "data" / "processed" / ".torch_cache"))
CURATION = ROOT / "data" / "raw" / "evidence_curation.csv"
SEQUENCES = ROOT / "data" / "processed" / "proteins_esm.csv"
HTML = ROOT / "web" / "index.html"
MODEL_DIR = ROOT / "results" / "model3_esm2"
sys.path.insert(0, str(ROOT / "src"))
AA = set("ACDEFGHIKLMNPQRSTVWY")
MAX_BODY = 250_000
_PREDICTOR = None


def dataset_status() -> dict:
    with CURATION.open(newline="", encoding="utf-8-sig") as handle:
        evidence = list(csv.DictReader(handle))
    with SEQUENCES.open(newline="", encoding="utf-8-sig") as handle:
        sequences = list(csv.DictReader(handle))
    status_counts = Counter(row.get("status", "") for row in evidence)
    sequence_counts = Counter(row.get("label_status", "") for row in sequences)
    sequence_genes = {row.get("gene_symbol", "") for row in sequences}
    positive_genes = {row.get("gene_symbol", "") for row in evidence
                      if row.get("status") == "positive"}
    missing_positive = sorted(positive_genes - sequence_genes)
    groups: dict[str, set[str]] = defaultdict(set)
    for row in sequences:
        label = row.get("label_status", "")
        group = (row.get("cluster_id") or "").strip()
        if label in {"positive", "candidate_negative"} and group:
            groups[label].add(group)
    positive_groups = len(groups["positive"])
    negative_groups = len(groups["candidate_negative"])
    ready = status_counts["candidate_negative"] > 0 and positive_groups >= 6 and negative_groups >= 6
    metrics_path = MODEL_DIR / "metrics.json"
    metrics = json.loads(metrics_path.read_text(encoding="utf-8")) if metrics_path.exists() else None
    return {
        "reviewed_candidates": len(evidence),
        "positive": status_counts["positive"],
        "candidate_negative": status_counts["candidate_negative"],
        "unresolved": status_counts["unresolved"],
        "verified_sequence_rows": len(sequences),
        "positive_sequence_rows": sequence_counts["positive"],
        "missing_positive_sequences": missing_positive,
        "positive_similarity_groups": positive_groups,
        "negative_similarity_groups": negative_groups,
        "minimum_group_count_floor_passed": ready,
        "model_trained": metrics is not None and (MODEL_DIR / "logistic_regression.joblib").exists(),
        "prediction_available": metrics is not None and (MODEL_DIR / "logistic_regression.joblib").exists(),
        "model_metrics": metrics,
        "demo_warning": "Exploratory research demo only. The test set has just two candidate-negative proteins, and both were misclassified; scores are not reliable for biological conclusions.",
    }


def clean_single_sequence(raw: str) -> tuple[str | None, str | None]:
    """Accept a raw sequence or one FASTA record; reject multi-record input."""
    lines = [line.strip() for line in raw.splitlines() if line.strip()]
    if not lines:
        return None, "Paste one amino-acid sequence or one FASTA record."
    headers = [i for i, line in enumerate(lines) if line.startswith(">")]
    if len(headers) > 1:
        return None, "Paste one protein at a time; this input contains multiple FASTA records."
    if headers and headers[0] != 0:
        return None, "A FASTA header must appear before the sequence."
    if headers:
        lines = lines[1:]
    sequence = "".join(lines).upper()
    if not sequence:
        return None, "No amino-acid sequence was found."
    invalid = sorted(set(sequence) - AA)
    if invalid:
        return None, "Unsupported character(s): " + ", ".join(invalid) + ". Use the 20 standard amino-acid letters."
    return sequence, None


class Handler(BaseHTTPRequestHandler):
    server_version = "StrokeProtDemo/1.0"

    def _send(self, status: int, body: bytes, content_type: str) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:
        route = urlparse(self.path).path
        if route == "/":
            try:
                body = HTML.read_bytes()
            except OSError:
                self._send(500, b"Demo page is missing.", "text/plain; charset=utf-8")
                return
            self._send(200, body, "text/html; charset=utf-8")
        elif route == "/api/status":
            try:
                body = json.dumps(dataset_status()).encode("utf-8")
            except (OSError, csv.Error, KeyError) as error:
                body = json.dumps({"error": f"Could not read project data: {error}"}).encode("utf-8")
                self._send(500, body, "application/json; charset=utf-8")
                return
            self._send(200, body, "application/json; charset=utf-8")
        else:
            self._send(404, b"Not found.", "text/plain; charset=utf-8")

    def do_POST(self) -> None:
        if urlparse(self.path).path not in {"/api/validate", "/api/predict"}:
            self._send(404, b"Not found.", "text/plain; charset=utf-8")
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
        except ValueError:
            length = 0
        if length <= 0 or length > MAX_BODY:
            self._send(413, b'{"error":"Input is empty or too large."}', "application/json; charset=utf-8")
            return
        try:
            payload = json.loads(self.rfile.read(length))
            raw = payload.get("sequence", "")
            if not isinstance(raw, str):
                raise ValueError("Sequence must be text.")
        except (json.JSONDecodeError, ValueError, AttributeError) as error:
            body = json.dumps({"valid": False, "message": str(error), "prediction_available": False}).encode("utf-8")
            self._send(400, body, "application/json; charset=utf-8")
            return
        sequence, error = clean_single_sequence(raw)
        response = {
            "valid": sequence is not None,
            "length": len(sequence) if sequence else 0,
            "message": error or "Sequence format is valid.",
            "prediction_available": False,
        }
        if sequence is not None and urlparse(self.path).path == "/api/predict":
            model_file = MODEL_DIR / "logistic_regression.joblib"
            metrics_file = MODEL_DIR / "metrics.json"
            if not model_file.exists() or not metrics_file.exists():
                response["prediction_message"] = "No trained model files are available."
            else:
                try:
                    import joblib
                    import esm
                    import torch
                    import numpy as np
                    from strokeprot.model3.esm import _sequence_embedding, CHECKPOINT
                    metrics = json.loads(metrics_file.read_text(encoding="utf-8"))
                    global _PREDICTOR
                    if _PREDICTOR is None:
                        model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
                        model.eval()
                        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
                        model.to(device)
                        _PREDICTOR = (model, alphabet, device, joblib.load(model_file))
                    model, alphabet, device, classifier = _PREDICTOR
                    vector = _sequence_embedding(sequence, model, alphabet,
                                                 alphabet.get_batch_converter(), device)
                    probability = float(classifier.predict_proba(np.asarray(vector).reshape(1, -1))[0, 1])
                    threshold = float(metrics["validation_selected_threshold"])
                    response.update({
                        "prediction_available": True,
                        "predicted_label": "stroke-associated under the project definition" if probability >= threshold else "not classified as stroke-associated under the project definition",
                        "model_score": probability,
                        "threshold": threshold,
                        "prediction_message": "Exploratory score only; this model is not reliable for biological conclusions.",
                    })
                except Exception as error:
                    response["prediction_message"] = f"Could not load or run the local model: {error}"
        self._send(200, json.dumps(response).encode("utf-8"), "application/json; charset=utf-8")

    def log_message(self, fmt: str, *args) -> None:
        print("[demo] " + fmt % args)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args()
    if not 1 <= args.port <= 65535:
        parser.error("--port must be between 1 and 65535")
    server = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    print(f"StrokeProt-AI demo: http://127.0.0.1:{args.port}")
    print("Local-only exploratory research demo; not a diagnostic or clinical tool.")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nDemo server stopped.")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
