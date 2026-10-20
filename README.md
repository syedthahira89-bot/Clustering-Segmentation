# Genomic CRM Segmentation Demo

A small, reproducible Python project built from `Untitled0.ipynb`. It turns the
notebook's synthetic genomic/clinical cohort, PCA visualization, patient
profiling, CRM mapping, CSV export, and ROI scenario into a testable pipeline
and local dashboard.

> **Synthetic demo only.** The generated cohort is artificial. This project is
> not a medical device, does not diagnose or recommend treatment, and has not
> been clinically validated. Do not upload real patient records or use its
> segment labels for care decisions. CRM output is a demonstration CSV; there
> is no live EHR or CRM connection.

## Presentation chart

![Genomic CRM segmentation summary chart](docs/images/segmentation-summary.png)

The one-slide summary shows the PCA cohort map, segment mix, and the
illustrative cost-versus-net-savings scenario. Regenerate it (optionally from
your own synthetic CSV) with:

```powershell
genomic-segmentation chart --output docs/images/segmentation-summary.png
```

The silhouette score is computed in the full scaled feature space, so it can be
low even when the 2-D PCA view shows clearly separated clusters.

## Run locally

From this folder in PowerShell:

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
streamlit run dashboard.py
```

The dashboard starts with deterministic synthetic data. You can also upload a
synthetic CSV with the schema below, choose a cluster count, inspect the PCA
map and cluster profiles, download a pseudonymous CRM roster, and adjust the
illustrative financial scenario.

## Command line

Generate a synthetic cohort:

```powershell
genomic-segmentation demo --output outputs/demo_patients.csv --patients 150 --genes 500
```

Cluster a CSV and export CRM campaign assignments:

```powershell
genomic-segmentation segment --input outputs/demo_patients.csv --output outputs/crm_roster.csv --clusters 3
```

The CRM export contains only `patient_id`, `cluster_id`, `segment_name`,
`next_best_action`, and `campaign_id`. It does not export genomic features,
clinical features, names, or the synthetic cohort label.

## Input CSV contract

Supply one row per pseudonymous patient. Required columns:

- `patient_id`: unique, non-empty pseudonymous identifier (not a name, email,
  medical record number, or other direct identifier).
- `age`, `bmi`, `cardiovascular_prs`, `diabetes_prs`,
  `portal_logins_30d`, `wellness_participation`: numeric clinical/engagement
  fields. Precomputed PRS fields are expected to be on a comparable 0–1 scale
  for the demo's profile labels.
- `pgx_status`: categorical value such as `normal`, `intermediate`, or `poor`.
- At least two numeric `gene_*` expression fields.

Numeric missing values are median-imputed; missing PGx categories are marked
`unknown`. If present, `demo_cohort` is ignored for clustering and excluded
from the CRM export. It exists only in generated data for demo validation.
Other fields are not used by the model.

## Method

1. Validate the input schema and pseudonymous IDs.
2. Impute and scale numeric genomic, clinical, and engagement features; encode
   PGx status.
3. Fit K-means to the complete scaled feature space and PCA to display the
   first two dimensions.
4. Summarize each cluster using clinical/engagement averages and the top five
   differentiating synthetic gene features.
5. Translate cohort profiles into the three notebook CRM campaigns.

The segment names are descriptive rules over unsupervised cluster profiles,
not diagnoses or individualized medical assessments. PCA is a visualization,
and the silhouette score is a mathematical separation metric, not a clinical
quality measure.

## Illustrative ROI assumptions

The calculator defaults reproduce the notebook's scenario: 20% baseline
readmission rate, 25% relative reduction, $15,200 per readmission, $75,000
fixed setup cost, and $150 variable cost per patient. At 1,000 patients this
produces $760,000 gross scenario savings, $225,000 project cost, and $535,000
net scenario savings. These are assumptions—not validated estimates, study
results, or guarantees of financial/clinical impact.

## Tests

```powershell
python -m pytest
```
