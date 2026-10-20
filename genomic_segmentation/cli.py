"""Command-line tools for generating demo data and exporting CRM segments."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .pipeline import cluster_patients, generate_demo_data


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Genomic CRM segmentation demo")
    commands = parser.add_subparsers(dest="command", required=True)

    demo = commands.add_parser("demo", help="Generate synthetic patient data")
    demo.add_argument("--output", default="outputs/demo_patients.csv")
    demo.add_argument("--patients", type=int, default=150)
    demo.add_argument("--genes", type=int, default=500)
    demo.add_argument("--seed", type=int, default=42)

    segment = commands.add_parser("segment", help="Cluster patients and export CRM rows")
    segment.add_argument("--input", required=True, help="Patient features CSV")
    segment.add_argument("--output", default="outputs/crm_roster.csv")
    segment.add_argument("--clusters", type=int, default=3)
    chart = commands.add_parser("chart", help="Create the presentation chart (PNG)")
    chart.add_argument("--input", help="Patient features CSV; defaults to synthetic demo data")
    chart.add_argument("--output", default="docs/images/segmentation-summary.png")
    chart.add_argument("--clusters", type=int, default=3)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    destination = Path(args.output)
    destination.parent.mkdir(parents=True, exist_ok=True)

    if args.command == "chart":
        from .charts import build_presentation_chart

        patients = pd.read_csv(args.input) if args.input else generate_demo_data()
        result = cluster_patients(patients, n_clusters=args.clusters)
        print(f"Wrote presentation chart to {build_presentation_chart(result, destination)}")
        return 0

    if args.command == "demo":
        generate_demo_data(args.patients, args.genes, args.seed).to_csv(
            destination, index=False
        )
        print(f"Wrote synthetic demo patients to {destination}")
        return 0

    result = cluster_patients(
        pd.read_csv(args.input), n_clusters=args.clusters
    )
    result.crm_export().to_csv(destination, index=False)
    print(
        f"Clustered {len(result.assignments)} patients "
        f"(silhouette={result.silhouette:.3f}); CRM roster: {destination}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
