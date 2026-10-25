# Local demo UI

The browser page reads the curation and sequence tables, displays their counts and held-out metrics, and accepts one pasted sequence or FASTA record. It runs the locally trained frozen ESM-2 plus logistic regression model and displays its score and threshold-based class. This is an exploratory software demo: the test set had only two provisional candidate-negative proteins, and the model missed both (ROC-AUC 0.17). Do not interpret the output as reliable biological classification, diagnosis, personal stroke risk, or causation.

## Start the demo

From the repository folder in a PowerShell terminal:

```powershell
.venv\Scripts\Activate.ps1
python app.py
```

Open `http://127.0.0.1:8000` in your browser. Keep the terminal open while presenting. Press `Ctrl+C` in that terminal to stop the local server. If port 8000 is already in use, run `python app.py --port 8001` and open `http://127.0.0.1:8001`.

The first prediction loads ESM-2 weights and may take several seconds; later sequences reuse the in-memory encoder. Keep the training dependencies and local model artifacts under `results/model3_esm2/` available. The first training run downloads weights and caches embeddings under ignored `data/processed/`. A raw sequence must use the 20 standard amino-acid letters. Its displayed score is not a calibrated clinical probability.
