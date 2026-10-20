import numpy as np
import pandas as pd
import pytest

from genomic_segmentation.pipeline import (
    cluster_patients,
    generate_demo_data,
    roi_projection,
)


def test_demo_generation_is_deterministic_and_excludes_direct_identifiers() -> None:
    first = generate_demo_data(n_patients=60, n_genes=20, seed=10)
    second = generate_demo_data(n_patients=60, n_genes=20, seed=10)

    pd.testing.assert_frame_equal(first, second)
    assert len([column for column in first if column.startswith("gene_")]) == 20
    assert "name" not in first.columns
    assert first["patient_id"].is_unique


def test_clustering_creates_crm_exports_and_profiles_without_label_leakage() -> None:
    patients = generate_demo_data(n_patients=90, n_genes=30, seed=12)
    result = cluster_patients(patients, n_clusters=3)
    changed_labels = patients.copy()
    changed_labels["demo_cohort"] = "changed label"
    result_with_changed_labels = cluster_patients(changed_labels, n_clusters=3)

    assert len(result.assignments) == len(patients)
    pd.testing.assert_series_equal(
        result.assignments["cluster_id"],
        result_with_changed_labels.assignments["cluster_id"],
    )
    assert result.assignments["cluster_id"].nunique() == 3
    assert set(result.crm_export().columns) == {
        "patient_id",
        "cluster_id",
        "segment_name",
        "next_best_action",
        "campaign_id",
    }
    assert set(result.assignments["campaign_id"]) <= {
        "CAMP-LONGEVITY",
        "CAMP-CARDIO-2026",
        "CAMP-PGX-ALERT",
    }
    assert "demo_cohort" not in result.crm_export().columns
    assert np.isfinite(result.explained_variance)
    assert -1 <= result.silhouette <= 1
    assert len(result.profiles) == 3


def test_clustering_coerces_numeric_strings_and_imputes_missing_values() -> None:
    patients = generate_demo_data(n_patients=30, n_genes=10)
    patients["bmi"] = patients["bmi"].astype(str)
    patients.loc[0, "bmi"] = np.nan

    result = cluster_patients(patients)

    assert result.profiles["bmi_mean"].notna().all()


def test_clustering_rejects_duplicate_patient_ids_and_missing_features() -> None:
    patients = generate_demo_data(n_patients=30, n_genes=10)
    patients.loc[1, "patient_id"] = patients.loc[0, "patient_id"]
    with pytest.raises(ValueError, match="unique"):
        cluster_patients(patients)

    patients = generate_demo_data(n_patients=30, n_genes=10).drop(columns="pgx_status")
    with pytest.raises(ValueError, match="missing required features"):
        cluster_patients(patients)


def test_presentation_chart_is_written(tmp_path) -> None:
    from genomic_segmentation.charts import build_presentation_chart

    result = cluster_patients(generate_demo_data(n_patients=60, n_genes=20), n_clusters=3)
    output = build_presentation_chart(result, tmp_path / "charts" / "summary.png")

    assert output.is_file()
    assert output.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_roi_projection_matches_notebook_assumptions() -> None:
    result = roi_projection(1000)

    assert result["baseline_readmissions"] == 200
    assert result["prevented_readmissions"] == 50
    assert result["gross_savings"] == 760_000
    assert result["project_cost"] == 225_000
    assert result["net_savings"] == 535_000


@pytest.mark.parametrize(
    "kwargs",
    [
        {"cohort_size": -1},
        {"cohort_size": 1.5},
        {"cohort_size": 10, "baseline_readmission_rate": 1.1},
        {"cohort_size": 10, "target_relative_reduction": -0.1},
        {"cohort_size": 10, "average_readmission_cost": np.inf},
    ],
)
def test_roi_projection_rejects_invalid_assumptions(kwargs: dict[str, object]) -> None:
    with pytest.raises(ValueError):
        roi_projection(**kwargs)
