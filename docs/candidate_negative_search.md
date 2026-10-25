# Candidate-negative search status

## Search rule

A candidate is not negative merely because it was in a panel without a statistically significant result. A `candidate_negative` label requires a documented search of primary human ischemic-stroke evidence over the agreed scope, review of any hits, and a conclusion that no qualifying association was found. The status means “no supporting evidence found under this search”; it is not a biological negative.

## Candidate pool and query

`data/raw/negative_search_pool.csv` now lists 1 unresolved protein from the primary human outcome panel in Angerfors et al. whose univariable p-values were above 0.05 in that cohort. The starting pool had 23 candidates. CCL19, CXCL1, CXCL5, IL18, MMP10, and CCL28 were reclassified positive after targeted primary-human literature review; nonsignificant panel results were not used as negative labels. This one cohort is context only and does not establish negative labels.

The helper `scripts/search_pubmed_candidates.py` constructs an alias-aware PubMed query for each candidate. It searches title/abstract aliases together with ischemic stroke, cerebral ischemia, cerebral infarction, stroke, and relevant MeSH headings. Run it with:

```powershell
python scripts/search_pubmed_candidates.py
```

It writes query strings, result counts, and returned PMIDs to `data/raw/pubmed_negative_search_results.csv`; each hit still needs manual review for primary research, human population, ischemic-stroke relevance, and a protein-specific association. Zero hits alone never assign a label.

## Current run and decisions

The carried-forward query log and subsequent searches now support 10 provisional candidate-negative decisions. They are not true negatives and remain subject to broader independent review. CST5 returned zero records for its logged aliases; this alone is insufficient, so CST5 remains unresolved. The repository has **10 candidate-negative labels**.

Targeted literature discovery did identify positive human evidence for proteins whose result was nonsignificant in the 2023 panel study: CXCL9 (PMID 36744204; DOI 10.2147/NDT.S391578) and CCL2/MCP-1 (PMID 39119971; DOI 10.1161/JAHA.124.035820). A separate prospective cohort supports FGF19 as an outcome-associated protein (PMID 40117962; DOI 10.1016/j.clnu.2025.03.013). Their full study details and cautions are in `data/raw/evidence_curation.csv`. These examples show why a nonsignificant result in one panel study cannot be used as a negative label.

## Additional relevant candidate source identified

Lagging et al. (PMID 36418382; DOI 10.1038/s41598-022-23288-5) is a primary human cohort relevant to the project. The three residual candidates KYNU, CSF2RA, and FCRL2 were searched and their hits manually reviewed. No qualifying protein-specific human ischemic-stroke evidence was confirmed in this limited search, so they are provisional candidate negatives. Query logs are in `data/raw/2022_91panel_pubmed_results.csv`; the supplement is retained at `data/raw/lagging_2022_91protein_supplement.xlsx`. This is not a comprehensive review of all 91 proteins. Sources: [primary article](https://www.nature.com/articles/s41598-022-23288-5), [PubMed record](https://pubmed.ncbi.nlm.nih.gov/36418382/).

## 2026 human HIV-associated stroke source

Chen et al. (PMID 42265127; DOI 10.1038/s41467-026-74258-8) report a primary human targeted-plasma proteomics study in people living with HIV. Its ischemic-stroke supplement has 27 Up/Down and 155 Ns entries. Alias-aware PubMed searches were completed and expanded for initial zero-hit entries. Seven proteins remained zero-hit under this search scope and were provisionally marked candidate_negative: CPA1, CTRC, SCGB3A2, HAO1, IL17D, PRSS27, and PIGR. Search pools/logs and hit abstracts are retained in `data/raw/2026_ischemic_candidate_search_pool.csv`, `data/raw/2026_ischemic_pubmed_search_results.csv`, `data/raw/2026_ischemic_zero_hit_expanded_pool.csv`, and `data/raw/2026_ischemic_pubmed_expanded_zero_results.csv`. This population is HIV-specific and not broadly representative. See the [primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC13402305/) and retained supplements.

## Training consequence

Dataset counts are 56 positive, 10 candidate_negative, and 2 unresolved (F8/CST5). Sequences exist for all 66 trainable rows. This is enough to execute an exploratory grouped split (51 positive and 10 candidate-negative sequence groups), but not enough to make a reliable classifier: test has only two candidate-negative proteins, both misclassified. Binary model files and honest metrics are available for demo purposes, with severe uncertainty. Broaden and independently review the negative evidence before treating the dataset as training-ready for research claims.
