#!/usr/bin/env python3
"""Plot empirical ARG double-recombination coverage across an inversion."""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot the fraction of logged ARG double-recombination intervals "
            "covering each position and compare it with the expected triangle."
        )
    )
    parser.add_argument(
        "input",
        type=Path,
        help="Run directory or arg_recombination_events.csv path.",
    )
    parser.add_argument("--left", type=float, required=True, help="Inversion left boundary in bp.")
    parser.add_argument("--right", type=float, required=True, help="Inversion right boundary in bp.")
    parser.add_argument("--points", type=int, default=201, help="Number of plotted positions.")
    parser.add_argument(
        "--output",
        type=Path,
        help="Output PNG path; defaults to arg_dr_triangle.png in the run directory.",
    )
    return parser.parse_args()


def read_dr_intervals(path):
    starts = []
    ends = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"event_type", "bp1_bp", "bp2_bp"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing CSV columns: {', '.join(sorted(missing))}")

        for row in reader:
            if row["event_type"].strip().lower() != "dr":
                continue
            starts.append(float(row["bp1_bp"]))
            ends.append(float(row["bp2_bp"]))

    return np.asarray(starts), np.asarray(ends)


def empirical_coverage(positions, starts, ends, midpoint):
    sorted_starts = np.sort(starts)
    sorted_ends = np.sort(ends)
    n_events = starts.size
    coverage = np.empty_like(positions)

    left_mask = positions <= midpoint
    coverage[left_mask] = (
        np.searchsorted(sorted_starts, positions[left_mask], side="right") / n_events
    )
    coverage[~left_mask] = (
        n_events
        - np.searchsorted(sorted_ends, positions[~left_mask], side="left")
    ) / n_events
    return coverage


def main():
    args = parse_args()
    if args.right <= args.left:
        raise SystemExit("--right must be greater than --left")
    if args.points < 3:
        raise SystemExit("--points must be at least 3")

    csv_path = args.input
    if csv_path.is_dir():
        csv_path = csv_path / "arg_recombination_events.csv"
    if not csv_path.is_file():
        raise SystemExit(f"Input CSV not found: {csv_path}")

    starts, ends = read_dr_intervals(csv_path)
    if starts.size == 0:
        raise SystemExit(f"No DR events found in {csv_path}")

    midpoint = 0.5 * (args.left + args.right)
    positions = np.linspace(args.left, args.right, args.points)
    empirical = empirical_coverage(positions, starts, ends, midpoint)

    invalid = int(
        np.count_nonzero(
            (starts < args.left)
            | (starts > midpoint)
            | (ends < midpoint)
            | (ends > args.right)
            | (starts >= ends)
        )
    )

    output_path = args.output or csv_path.with_name("arg_dr_triangle.png")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(positions, empirical, color="#3366aa", linewidth=2, label="Observed DR coverage")
    ax.set(
        xlabel="Position in inversion (bp)",
        ylabel="Fraction of DR events covering position",
        ylim=(-0.02, 1.02),
        title=f"ARG double-recombination coverage (n = {starts.size:,})",
    )
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    print(f"DR events: {starts.size:,}")
    print(f"Invalid breakpoint pairs: {invalid:,}")
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
