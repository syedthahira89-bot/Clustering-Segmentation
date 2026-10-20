"""Synthetic genomic and clinical segmentation demo."""

from .pipeline import (
    SegmentationResult,
    cluster_patients,
    generate_demo_data,
    roi_projection,
)

__all__ = [
    "SegmentationResult",
    "cluster_patients",
    "generate_demo_data",
    "roi_projection",
]
