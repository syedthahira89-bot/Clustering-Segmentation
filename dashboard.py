"""Streamlit dashboard for synthetic genomic patient segmentation."""

from __future__ import annotations

import pandas as pd
import streamlit as st

from genomic_segmentation.pipeline import cluster_patients, generate_demo_data, roi_projection

st.set_page_config(
    page_title="Precision Care | Cohort Explorer",
    page_icon="P",
    layout="wide",
)
st.title("Precision Care · Cohort Explorer")
st.caption(
    "Explore synthetic genomic and clinical cohorts, inspect cluster profiles, "
    "and create a CRM-ready CSV."
)
st.warning(
    "Demonstration only. Use synthetic data here; this tool does not provide "
    "diagnoses, treatment recommendations, or validated clinical risk estimates."
)

with st.sidebar:
    st.header("Data and clustering")
    data_source = st.radio("Patient data", ["Generate synthetic demo", "Upload CSV"])
    if data_source == "Generate synthetic demo":
        patient_count = st.slider("Synthetic patients", 60, 600, 150, step=30)
        gene_count = st.slider("Synthetic gene features", 20, 500, 200, step=20)
        frame = generate_demo_data(patient_count, gene_count)
    else:
        upload = st.file_uploader("Synthetic patient CSV", type=["csv"])
        frame = pd.read_csv(upload) if upload is not None else None
    cluster_count = st.slider("Number of clusters", 2, 8, 3)

    st.divider()
    st.header("Illustrative pilot scenario")
    roi_patients = st.number_input("Pilot cohort size", min_value=0, value=1000, step=100)
    baseline_rate = st.slider("Assumed baseline readmission rate", 0.0, 1.0, 0.20)
    reduction = st.slider("Assumed relative reduction", 0.0, 1.0, 0.25)
    readmission_cost = st.number_input(
        "Assumed cost per readmission ($)", min_value=0.0, value=15_200.0, step=500.0
    )
    fixed_cost = st.number_input(
        "Assumed fixed setup cost ($)", min_value=0.0, value=75_000.0, step=5_000.0
    )
    variable_cost = st.number_input(
        "Assumed variable cost per patient ($)",
        min_value=0.0,
        value=150.0,
        step=25.0,
    )

if frame is None:
    st.info("Generate the synthetic demo or upload a CSV to begin.")
    st.stop()

st.subheader("Input preview")
st.dataframe(
    frame.drop(columns=["demo_cohort"], errors="ignore").head(8),
    width="stretch",
)

try:
    result = cluster_patients(frame, n_clusters=cluster_count)
except (ValueError, TypeError) as error:
    st.error(f"Could not segment this file: {error}")
    st.stop()

metric_columns = st.columns(3)
metric_columns[0].metric("Patients segmented", f"{len(result.assignments):,}")
metric_columns[1].metric("PCA variance shown", f"{result.explained_variance:.1%}")
metric_columns[2].metric("Silhouette score", f"{result.silhouette:.3f}")
st.caption(
    "The silhouette score summarizes separation in this dataset; it is not a measure "
    "of clinical validity."
)

left, right = st.columns([1.2, 1])
with left:
    st.subheader("PCA cohort map")
    st.scatter_chart(
        result.assignments,
        x="pca_1",
        y="pca_2",
        color="segment_name",
        width="stretch",
    )
with right:
    st.subheader("Cluster profiles")
    st.dataframe(
        result.profiles[
            [
                "cluster_id",
                "segment_name",
                "patient_count",
                "cardiovascular_prs_mean",
                "diabetes_prs_mean",
                "pgx_concern_rate",
                "top_differentiating_genes",
            ]
        ].round(3),
        width="stretch",
        hide_index=True,
    )

st.subheader("CRM assignment preview")
crm_export = result.crm_export()
st.dataframe(crm_export, width="stretch", hide_index=True)
st.download_button(
    "Download CRM roster CSV",
    data=crm_export.to_csv(index=False).encode("utf-8"),
    file_name="genomic_crm_roster.csv",
    mime="text/csv",
)

scenario = roi_projection(
    int(roi_patients),
    baseline_readmission_rate=baseline_rate,
    target_relative_reduction=reduction,
    average_readmission_cost=readmission_cost,
    fixed_setup_cost=fixed_cost,
    variable_cost_per_patient=variable_cost,
)
st.subheader("Illustrative financial scenario")
st.caption(
    "These notebook assumptions are scenario inputs, not evidence of clinical impact "
    "or guaranteed savings. Replace them with locally validated, approved estimates."
)
roi_columns = st.columns(4)
roi_columns[0].metric(
    "Prevented readmissions (expected)", f"{scenario['prevented_readmissions']:,.1f}"
)
roi_columns[1].metric("Gross savings", f"${scenario['gross_savings']:,.0f}")
roi_columns[2].metric("Estimated project cost", f"${scenario['project_cost']:,.0f}")
roi_columns[3].metric("Net scenario", f"${scenario['net_savings']:,.0f}")
