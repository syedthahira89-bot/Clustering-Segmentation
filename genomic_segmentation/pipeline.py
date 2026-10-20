"""Synthetic cohort generation, clustering, CRM routing, and ROI scenarios."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.metrics import silhouette_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

PATIENT_ID = "patient_id"
DEMO_LABEL = "demo_cohort"
CLINICAL_NUMERIC = (
    "age",
    "bmi",
    "cardiovascular_prs",
    "diabetes_prs",
    "portal_logins_30d",
    "wellness_participation",
)
CLINICAL_CATEGORICAL = ("pgx_status",)
PGX_CONCERN_STATUSES = {"intermediate", "poor", "suboptimal", "slow"}

CRM_ACTIONS = {
    "Proactive Wellness": (
        "Offer optional wellness resources; obtain consent before outreach.",
        "CAMP-LONGEVITY",
    ),
    "High-Risk Cardiological": (
        "Route for qualified preventive-care review; do not infer a diagnosis.",
        "CAMP-CARDIO-2026",
    ),
    "PGx Suboptimal Metabolizer": (
        "Request qualified pharmacist review; do not change medication automatically.",
        "CAMP-PGX-ALERT",
    ),
}


@dataclass(frozen=True)
class SegmentationResult:
    """Patient assignments, clinical cluster profiles, and validation metrics."""

    assignments: pd.DataFrame
    profiles: pd.DataFrame
    explained_variance: float
    silhouette: float

    def crm_export(self) -> pd.DataFrame:
        """Return the pseudonymous patient-to-campaign payload, without biomarkers."""
        return self.assignments[
            [PATIENT_ID, "cluster_id", "segment_name", "next_best_action", "campaign_id"]
        ].copy()


def generate_demo_data(
    n_patients: int = 150, n_genes: int = 500, seed: int = 42
) -> pd.DataFrame:
    """Create a reproducible, entirely synthetic three-cohort dataset."""
    if n_patients < 6:
        raise ValueError("n_patients must be at least 6.")
    if n_genes < 2:
        raise ValueError("n_genes must be at least 2.")

    rng = np.random.default_rng(seed)
    cohort_sizes = np.full(3, n_patients // 3, dtype=int)
    cohort_sizes[: n_patients % 3] += 1
    expression = rng.normal(2.0, 0.5, size=(n_patients, n_genes))
    labels = np.repeat(
        ["Proactive Wellness", "Cardiometabolic Risk", "PGx Risk"], cohort_sizes
    )

    first_block = max(1, n_genes // 10)
    second_start = first_block
    second_end = min(n_genes, 2 * first_block)
    first_cohort_end = cohort_sizes[0]
    second_cohort_end = first_cohort_end + cohort_sizes[1]
    expression[first_cohort_end:second_cohort_end, :first_block] += 1.5
    expression[second_cohort_end:, second_start:second_end] += 2.0

    cardiovascular_prs = rng.beta(2, 5, size=n_patients)
    diabetes_prs = rng.beta(2, 5, size=n_patients)
    bmi = rng.normal(26, 4, size=n_patients).clip(16, 48)
    pgx_status = rng.choice(
        ["normal", "intermediate", "poor"],
        size=n_patients,
        p=[0.85, 0.12, 0.03],
    )

    cardiovascular_slice = slice(first_cohort_end, second_cohort_end)
    pgx_slice = slice(second_cohort_end, n_patients)
    cardiovascular_prs[cardiovascular_slice] = rng.uniform(
        0.72, 0.99, size=cohort_sizes[1]
    )
    diabetes_prs[cardiovascular_slice] = rng.uniform(0.62, 0.95, size=cohort_sizes[1])
    bmi[cardiovascular_slice] += 5
    pgx_status[pgx_slice] = rng.choice(
        ["intermediate", "poor"], size=cohort_sizes[2], p=[0.45, 0.55]
    )

    data: dict[str, object] = {
        PATIENT_ID: [f"DEMO-{index:05d}" for index in range(1, n_patients + 1)],
        "age": rng.integers(18, 86, size=n_patients),
        "bmi": bmi.round(1),
        "cardiovascular_prs": cardiovascular_prs.round(4),
        "diabetes_prs": diabetes_prs.round(4),
        "pgx_status": pgx_status,
        "portal_logins_30d": rng.poisson(4, size=n_patients),
        "wellness_participation": rng.binomial(1, 0.3, size=n_patients),
        DEMO_LABEL: labels,
    }
    data.update(
        {
            f"gene_{index + 1:03d}": expression[:, index].round(4)
            for index in range(n_genes)
        }
    )
    return pd.DataFrame(data)


def _validate_and_select_features(frame: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    if frame.empty:
        raise ValueError("Input data contains no patient rows.")
    if frame.columns.duplicated().any():
        raise ValueError("Input data contains duplicate column names.")
    if PATIENT_ID not in frame:
        raise ValueError(f"Input data must include a pseudonymous '{PATIENT_ID}' column.")
    patient_ids = frame[PATIENT_ID].astype("string").str.strip()
    if patient_ids.isna().any() or patient_ids.eq("").any():
        raise ValueError("Patient IDs must be present and non-empty.")
    if patient_ids.duplicated().any():
        raise ValueError("Patient IDs must be unique within the input file.")

    missing = (
        set(CLINICAL_NUMERIC + CLINICAL_CATEGORICAL) - set(frame.columns)
    )
    if missing:
        raise ValueError(f"Input data is missing required features: {sorted(missing)}")
    gene_columns = sorted(column for column in frame if column.startswith("gene_"))
    if len(gene_columns) < 2:
        raise ValueError("Input data must contain at least two numeric gene_* features.")

    features = frame[list(CLINICAL_NUMERIC + CLINICAL_CATEGORICAL) + gene_columns].copy()
    for column in CLINICAL_NUMERIC + tuple(gene_columns):
        features[column] = pd.to_numeric(features[column], errors="coerce")
        if features[column].notna().sum() == 0:
            raise ValueError(f"Feature '{column}' contains no numeric values.")
    if features[list(CLINICAL_NUMERIC) + gene_columns].isna().all(axis=1).any():
        raise ValueError("Each patient must have at least one observed numeric feature.")
    features["pgx_status"] = features["pgx_status"].astype("string").fillna("unknown")
    if features.empty:
        raise ValueError("No usable model features remain.")
    return features, gene_columns


def _describe_segment(profile: pd.Series) -> str:
    if profile["pgx_concern_rate"] >= 0.35:
        return "PGx Suboptimal Metabolizer"
    if (
        profile["cardiovascular_prs_mean"] >= 0.70
        or profile["diabetes_prs_mean"] >= 0.70
    ):
        return "High-Risk Cardiological"
    return "Proactive Wellness"


def cluster_patients(
    frame: pd.DataFrame, n_clusters: int = 3, random_state: int = 42
) -> SegmentationResult:
    """Cluster patients using scaled genomic, clinical, and engagement features."""
    if n_clusters < 2:
        raise ValueError("n_clusters must be at least 2.")
    if len(frame) <= n_clusters:
        raise ValueError("Input data must contain more patients than clusters.")

    features, gene_columns = _validate_and_select_features(frame)
    numeric_columns = list(CLINICAL_NUMERIC) + gene_columns
    preprocessor = ColumnTransformer(
        transformers=[
            (
                "numeric",
                Pipeline(
                    steps=[
                        ("imputer", SimpleImputer(strategy="median")),
                        ("scaler", StandardScaler()),
                    ]
                ),
                numeric_columns,
            ),
            (
                "categorical",
                Pipeline(
                    steps=[
                        ("imputer", SimpleImputer(strategy="most_frequent")),
                        ("encoder", OneHotEncoder(handle_unknown="ignore")),
                    ]
                ),
                list(CLINICAL_CATEGORICAL),
            ),
        ],
        remainder="drop",
    )
    matrix = preprocessor.fit_transform(features)
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()

    pca = PCA(n_components=2, random_state=random_state)
    embedding = pca.fit_transform(matrix)
    cluster_ids = KMeans(
        n_clusters=n_clusters, n_init=10, random_state=random_state
    ).fit_predict(matrix)
    if np.unique(cluster_ids).size < 2:
        raise ValueError("Clustering produced fewer than two distinct clusters.")

    profile_frame = frame[[PATIENT_ID, *CLINICAL_NUMERIC, *CLINICAL_CATEGORICAL]].copy()
    for column in CLINICAL_NUMERIC:
        values = features[column]
        profile_frame[column] = values.fillna(values.median()).to_numpy()
    profile_frame["pgx_status"] = features["pgx_status"].to_numpy()
    profile_frame["cluster_id"] = cluster_ids
    profiles = profile_frame.groupby("cluster_id", sort=True).agg(
        patient_count=(PATIENT_ID, "size"),
        cardiovascular_prs_mean=("cardiovascular_prs", "mean"),
        diabetes_prs_mean=("diabetes_prs", "mean"),
        bmi_mean=("bmi", "mean"),
        portal_logins_mean=("portal_logins_30d", "mean"),
        wellness_participation_rate=("wellness_participation", "mean"),
        pgx_concern_rate=(
            "pgx_status",
            lambda values: values.astype("string")
            .str.strip()
            .str.lower()
            .isin(PGX_CONCERN_STATUSES)
            .mean(),
        ),
    )
    profiles["segment_name"] = profiles.apply(_describe_segment, axis=1)

    top_gene_names: dict[int, str] = {}
    expression_frame = features[gene_columns]
    overall_means = expression_frame.mean()
    for cluster_id in sorted(np.unique(cluster_ids)):
        group = expression_frame.loc[cluster_ids == cluster_id]
        differences = (group.mean() - overall_means).abs()
        top_gene_names[cluster_id] = ", ".join(differences.nlargest(5).index)
    profiles["top_differentiating_genes"] = pd.Series(top_gene_names)
    profiles = profiles.reset_index()

    segment_names = profiles.set_index("cluster_id")["segment_name"]
    assignments = pd.DataFrame(
        {
            PATIENT_ID: frame[PATIENT_ID].astype("string").to_numpy(),
            "cluster_id": cluster_ids,
            "segment_name": pd.Series(cluster_ids).map(segment_names).to_numpy(),
            "pca_1": embedding[:, 0],
            "pca_2": embedding[:, 1],
        }
    )
    assignments["next_best_action"] = assignments["segment_name"].map(
        {segment: action for segment, (action, _) in CRM_ACTIONS.items()}
    )
    assignments["campaign_id"] = assignments["segment_name"].map(
        {segment: campaign for segment, (_, campaign) in CRM_ACTIONS.items()}
    )
    if assignments[["next_best_action", "campaign_id"]].isna().any().any():
        raise ValueError("A segment did not map to a configured CRM action.")

    return SegmentationResult(
        assignments=assignments,
        profiles=profiles,
        explained_variance=float(pca.explained_variance_ratio_.sum()),
        silhouette=float(silhouette_score(matrix, cluster_ids)),
    )


def roi_projection(
    cohort_size: int,
    baseline_readmission_rate: float = 0.20,
    target_relative_reduction: float = 0.25,
    average_readmission_cost: float = 15_200.0,
    fixed_setup_cost: float = 75_000.0,
    variable_cost_per_patient: float = 150.0,
) -> dict[str, float]:
    """Calculate an illustrative scenario; this is not a clinical or ROI forecast."""
    if isinstance(cohort_size, bool) or cohort_size < 0:
        raise ValueError("cohort_size must be a non-negative integer.")
    if int(cohort_size) != cohort_size:
        raise ValueError("cohort_size must be a non-negative integer.")
    if not np.isfinite(baseline_readmission_rate) or not 0 <= baseline_readmission_rate <= 1:
        raise ValueError("baseline_readmission_rate must be between 0 and 1.")
    if not np.isfinite(target_relative_reduction) or not 0 <= target_relative_reduction <= 1:
        raise ValueError("target_relative_reduction must be between 0 and 1.")
    costs = (average_readmission_cost, fixed_setup_cost, variable_cost_per_patient)
    if any(not np.isfinite(cost) or cost < 0 for cost in costs):
        raise ValueError("Cost assumptions must be finite and non-negative.")

    baseline = cohort_size * baseline_readmission_rate
    prevented = baseline * target_relative_reduction
    gross_savings = prevented * average_readmission_cost
    project_cost = fixed_setup_cost + cohort_size * variable_cost_per_patient
    return {
        "cohort_size": float(cohort_size),
        "baseline_readmissions": float(baseline),
        "target_readmissions": float(baseline - prevented),
        "prevented_readmissions": float(prevented),
        "gross_savings": float(gross_savings),
        "project_cost": float(project_cost),
        "net_savings": float(gross_savings - project_cost),
        "roi_ratio": float(gross_savings / project_cost) if project_cost else 0.0,
    }
