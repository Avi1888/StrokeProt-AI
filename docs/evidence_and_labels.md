# Evidence and label rules

## Task definition

For one human protein sequence, this project asks whether the protein meets the documented definition of **stroke-associated**. This is a protein-level research classification label. It does not diagnose a person, predict an individual's stroke, or show that a protein causes stroke.

The scope is ischemic stroke. Assign `positive` only when a primary human research study directly reports the protein in at least one of these separately recorded categories:

1. `risk_before_stroke`: protein measured before an incident ischemic stroke and associated with later stroke occurrence/risk.
2. `post_event_measurement`: protein measured in people after ischemic stroke in relation to the event (including case/control or temporal measurements).
3. `severity_or_outcome`: protein measured after stroke in relation to severity, complications, recovery, or functional/cognitive outcome.

These categories describe different questions and are not interchangeable. A post-event signal is not evidence of pre-event risk. A reported association is not proof of causation.

## Statuses

- `positive`: a qualifying primary human ischemic-stroke study is identified, and the rationale cites what it measured and reported.
- `candidate_negative`: only after a documented search of the agreed sources and query scope finds no qualifying human primary evidence. This means no supporting evidence found under that search; it is not a biological negative. Current provisional candidate negatives are described in `docs/candidate_negative_search.md`.
- `unresolved`: evidence, protein identity, or protein-specific interpretation is ambiguous. Exclude it from training but retain its rationale in `data/raw/evidence_curation.csv`.

Review articles, animal-only experiments, database annotations, and candidate lists do not establish a positive label. They may help locate primary studies but cannot substitute for one. Expression changes alone do not qualify unless the primary human study relates the protein to one of the categories above.

## Current review scope and limitations

The review began with carried-forward human protein candidates and added targeted primary human evidence. It includes published papers and one clearly identified preprint (AGRP, PMID 40832370). The current evidence table has 68 rows: 56 positive, 10 provisional candidate_negative, and 2 unresolved (F8 and CST5). All 66 labeled candidates have verified human RefSeq sequences. The 10 provisional negatives arose from completed alias-aware PubMed searches of seven 2026-supplement proteins and three residual candidates in a 2022 human-proteomics panel; details and saved query logs are in `docs/candidate_negative_search.md`. This literature review remains targeted, not systematic. See each curation row for the paper title, DOI/PMID, evidence type, human population, rationale, and uncertainty.

CST5 remains unresolved after a zero-hit alias-aware PubMed query, which is not treated as proof of absence; their nonsignificant result in one cohort is not treated as a negative. The table records title, DOI, PMID, population, evidence type, rationale, sequence status, and source URL. `unresolved` cases remain visible. Forty-eight positive accessions are currently available as verified sequences; F8 and CST5 remain unresolved and unsequenced.

## Current decisions

| Protein | Status | Category | Decision summary |
|---|---|---|---|
| CRP | positive | risk before stroke | WHI nested case-control study: baseline CRP associated with incident ischemic stroke and remained an independent single biomarker. |
| IL6 | positive | risk before stroke | Same WHI study: IL-6 was among individually predictive biomarkers, but it was not the only independent single predictor after mutual biomarker adjustment. |
| MMP9 | positive | post-event measurement; severity/outcome | CATIS cohort: acute serum MMP-9 associated with death or major disability at 3 months. |
| VWF | positive | post-event measurement; severity/outcome | Stroke-registry study reports adverse outcome for joint FVIII+/VWF+ category; retain caveat that independent VWF effect is not isolated. |
| F8 | unresolved | post-event measurement; severity/outcome | Same registry's joint category does not isolate an F8 effect, and the paper attributes the observed relationship to VWF. |
| TNFSF10, S100A12, OSM, FLT3LG, TNFSF14, CSF1, CCL25, HGF | positive | severity/outcome | Each was selected by LASSO in the 534-case primary cohort for separating favorable and unfavorable 3-month outcome; directions are recorded in the evidence table. |
| CCL7, CD6, KITLG, TNFSF11, IL10, CCL23, EIF4EBP1, SIRT2, DNER, CASP8, STAMBP | positive | severity/outcome | Supplementary Table S2 reports FDR-significant univariable outcome associations in the same cohort; this is an association label, not evidence of causality. |
| CXCL9, CCL2 (MCP-1), FGF19 | positive | post-event measurement; severity/outcome | Separate primary human studies report protein measurements related to ischemic stroke or its outcomes; the curation table records populations and cautions. |
| CCL19 | positive | post-event measurement; severity/outcome | Two prospective human ischemic-stroke cohorts (PMID 37399099) associated plasma CCL19 with 3-month outcome. |
| CXCL1, CXCL5, CCL28 | positive | post-event measurement; severity/outcome | A controlled human large-vessel stroke comparison (PMID 36179808) reported post-event plasma differences; cohort is selected for thrombectomy. |
| IL18 | positive | post-event measurement | Acute ischemic-stroke patients had higher serum IL-18 than controls (PMID 21424610; 45 patients, 40 controls). |
| MMP10 | positive | post-event measurement; severity/outcome | Human acute ischemic-stroke cohorts linked serum proMMP-10 with brain damage and 3-month outcome (PMID 23742289). |
| CD244, CCL13, CXCL6, CCL8 | positive | post-event measurement | Primary human acute ischemic-stroke study measured these proteins in pericyte-derived plasma microvesicles; CD244, MCP-4/CCL13, and CXCL6 increased at 12–24 h versus controls, while MCP-2/CCL8 showed a nominal increase (PMID 38323422). Interpret as microvesicle cargo signals, not free plasma concentrations. |
| CXCL11 | positive | post-event measurement | Human thrombectomy study reported local CXCL11 elevation within the occluded arterial compartment in hyperacute ischemic stroke (PMID 31899551); selected large-vessel cohort. |
| FGF23 | positive | post-event measurement | In 262 first-ever ischemic-stroke patients, plasma FGF23 was associated with intracranial cerebral atherosclerosis (PMID 31942755). |
| FGF21 | positive | severity/outcome | CATIS cohort of 3,412 acute ischemic-stroke patients: higher plasma FGF21 associated with death or major disability at one year (PMID 33320402). |
| IFNG | positive | post-event measurement; severity/outcome | Small case-control study (22 ischemic stroke, 22 controls) found higher salivary IFN-γ; level correlated with cognitive score (PMID 40520895). Replication needed. |
| CX3CL1 | positive | post-event measurement; severity/outcome | Serial plasma study of 55 acute ischemic-stroke patients associated CX3CL1 temporal dynamics with severity and 90-day outcome (PMID 24722201). |
| ADA | positive | post-event measurement | Retrospective human case-control study of 3,968 acute cerebral infarction cases and 3,945 controls found lower serum ADA activity in cases (PMID 36260871). |
| CDCP1 | positive | severity/outcome | In the 534-case primary human cohort, higher CDCP1 was associated with unfavorable prognosis in the cryptogenic-stroke subgroup and selected by LASSO (PMID 37794467); exploratory subgroup finding. |
| SLAMF1, IL10RB | positive | severity/outcome | Both contributed to LASSO outcome separation in women in the same 534-case human cohort (PMID 37794467); sex-stratified exploratory findings. |
| NTF3 (NT-3) | positive | severity/outcome | Higher NT-3 was associated with favorable outcome in the cryptogenic subgroup in the same cohort (PMID 37794467); exploratory subgroup finding. |
| LIFR | positive | severity/outcome | Human thrombectomy cohort (81 subjects) found higher LIFR associated with lower discharge MoCA in non-Appalachian subjects (PMID 39154085; DOI 10.1186/s12974-024-03201-9); selected cohort and subgroup-specific. |
| CST5 | unresolved | — | Alias-aware PubMed query returned zero hits; the nonsignificant result in one panel cohort does not make it a biological negative. CST5/cystatin D is distinct from CST3/cystatin C. |
| ALCAM | positive | severity/outcome | Human acute ischemic-stroke cohort of 244 patients; admission ALCAM related to long-term mortality (PMID 21757661). |
| CD84 | positive | post-event measurement | Human arterial blood from ischemic cerebral circulation showed local CD84 shedding (PMID 32762491). |
| LGALS4 | positive | risk before stroke | Population cohort (1,688 total; 59 prevalent ischemic strokes) reported adjusted association with plasma galectin-4 (PMID 37373212). |
| PRTN3 | positive | post-event measurement | Human serum proteomics discovery and validation cohorts reported protein as a component of an ischemic-stroke multianalyte model (PMID 39167481). |
| BLMH | positive | post-event measurement | Differential abundance in 24 acute ischemic-stroke thrombectomy clots by HCAS status (PMID 36878687). |
| CCL15, IL16 | positive | post-event measurement | Human acute ischemic-stroke cytokine cohort (60 cases, 20 controls) reported differences by stroke subtype (PMID 39831358). |
| AGRP | positive | risk before stroke | SIREN study (90 ischemic-stroke cases, 90 controls) measured Olink plasma; source is a preprint and needs peer-reviewed replication (PMID 40832370). |
| CPA1, CTRC, SCGB3A2, HAO1, IL17D, PRSS27, PIGR | candidate_negative | — | No qualifying human protein-specific ischemic-stroke evidence found in logged alias-aware PubMed searches of the 2026 HIV-associated stroke study's nonsignificant panel subset. Provisional absence-of-evidence statuses, not biological negatives. |
| KYNU, CSF2RA, FCRL2 | candidate_negative | — | Limited alias-aware searches and manual hit review seeded from the 2022 human 91-protein panel found no qualifying protein-specific human ischemic-stroke association. Not exhaustive; provisional statuses only. |
| IL12B | positive | post-event measurement | A primary human large-vessel acute ischemic-stroke cohort measured the IL12B-encoded IL-12/IL-23 p40 analyte in thrombi and plasma compartments (PMID 35673992); this is not a specific measure of intact IL-12p70. |


The 11 additional FDR-significant proteins are CCL7, CD6, KITLG, TNFSF11, IL10, CCL23, EIF4EBP1, SIRT2, DNER, CASP8, and STAMBP; individual paper details are in the evidence table. Similarity groups were clustered at 30% identity and 80% aligned coverage for both sequences. Model evaluation is exploratory only: the 11-protein test split included two candidate negatives and missed both. See `docs/dataset_manifest.md` for metrics and exact split counts. Re-review all decisions before scientific reuse.
