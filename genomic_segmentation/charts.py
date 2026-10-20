"""Presentation chart summarizing the synthetic segmentation demo."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from .pipeline import SegmentationResult, roi_projection

SEGMENT_COLORS = {
    "Proactive Wellness": "#1f77b4",
    "High-Risk Cardiological": "#2ca02c",
    "PGx Suboptimal Metabolizer": "#9467bd",
}
COHORT_SIZES = (500, 1000, 5000, 10000, 25000)


def build_presentation_chart(result: SegmentationResult, output: str | Path) -> Path:
    """Save a three-panel summary figure and return its path."""
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)

    assignments = result.assignments
    counts = assignments["segment_name"].value_counts()
    colors = [SEGMENT_COLORS.get(name, "#7f7f7f") for name in counts.index]

    fig = plt.figure(figsize=(16, 9), facecolor="white")
    grid = fig.add_gridspec(2, 2, width_ratios=[1.3, 1], hspace=0.38, wspace=0.5)
    scatter_ax = fig.add_subplot(grid[:, 0])
    mix_ax = fig.add_subplot(grid[0, 1])
    roi_ax = fig.add_subplot(grid[1, 1])

    for name, group in assignments.groupby("segment_name"):
        scatter_ax.scatter(
            group["pca_1"],
            group["pca_2"],
            s=70,
            alpha=0.8,
            label=f"{name} (n={len(group)})",
            color=SEGMENT_COLORS.get(name, "#7f7f7f"),
            edgecolor="white",
            linewidth=0.6,
        )
    scatter_ax.set_title("Patient cohorts in genomic + clinical space (PCA)", fontsize=14, loc="left")
    scatter_ax.set_xlabel("Principal component 1")
    scatter_ax.set_ylabel("Principal component 2")
    scatter_ax.grid(True, linestyle="--", alpha=0.4)
    scatter_ax.legend(loc="best", frameon=True)
    scatter_ax.text(
        0.01,
        0.01,
        f"2 components explain {result.explained_variance:.1%} of variance · "
        f"silhouette {result.silhouette:.2f}",
        transform=scatter_ax.transAxes,
        fontsize=9,
        color="#555555",
    )

    bars = mix_ax.barh(counts.index[::-1], counts.values[::-1], color=colors[::-1])
    mix_ax.bar_label(
        bars,
        labels=[f"{v} ({v / counts.sum():.0%})" for v in counts.values[::-1]],
        padding=4,
    )
    mix_ax.set_title("Segment mix", fontsize=14, loc="left")
    mix_ax.set_xlabel("Patients")
    mix_ax.set_xlim(0, counts.max() * 1.25)
    mix_ax.grid(True, axis="x", linestyle="--", alpha=0.4)

    scenarios = [roi_projection(size) for size in COHORT_SIZES]
    positions = range(len(scenarios))
    width = 0.38
    roi_ax.bar(
        [p - width / 2 for p in positions],
        [s["project_cost"] for s in scenarios],
        width,
        label="Project cost",
        color="#d95f5f",
    )
    roi_ax.bar(
        [p + width / 2 for p in positions],
        [s["net_savings"] for s in scenarios],
        width,
        label="Net savings",
        color="#3aa76d",
    )
    roi_ax.set_xticks(list(positions))
    roi_ax.set_xticklabels([f"{size:,}" for size in COHORT_SIZES])
    roi_ax.yaxis.set_major_formatter(
        FuncFormatter(lambda value, _: f"${value / 1_000_000:.1f}M")
    )
    roi_ax.set_title("Illustrative scenario: cost vs. net savings", fontsize=14, loc="left")
    roi_ax.set_xlabel("Pilot cohort size")
    roi_ax.grid(True, axis="y", linestyle="--", alpha=0.4)
    roi_ax.legend(loc="upper left")

    fig.suptitle(
        "Genomic CRM Segmentation — Synthetic Demo", fontsize=20, fontweight="bold", x=0.05, ha="left"
    )
    fig.text(
        0.05,
        0.015,
        "Synthetic data only. Not clinically validated. ROI uses notebook assumptions "
        "(20% baseline readmission, 25% reduction, $15,200/readmission), not evidence of impact.",
        fontsize=9,
        color="#555555",
    )
    fig.subplots_adjust(top=0.88, bottom=0.08, left=0.06, right=0.97)
    fig.savefig(destination, dpi=150)
    plt.close(fig)
    return destination
