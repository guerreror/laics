#!/usr/bin/env python3
"""Plot observed ARG gene-conversion start positions and tract lengths."""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input",
        type=Path,
        help="Run directory or arg_recombination_events.csv path.",
    )
    parser.add_argument("--left", type=float, required=True)
    parser.add_argument("--right", type=float, required=True)
    parser.add_argument("--tract-length", type=float, default=200.0)
    parser.add_argument("--start-bins", type=int, default=50)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def read_gc_events(path):
    starts = []
    ends = []
    lengths = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"event_type", "bp1_bp", "bp2_bp", "tract_length_bp"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing CSV columns: {', '.join(sorted(missing))}")

        for row in reader:
            if row["event_type"].strip().lower() != "gc":
                continue
            starts.append(float(row["bp1_bp"]))
            ends.append(float(row["bp2_bp"]))
            lengths.append(float(row["tract_length_bp"]))

    return np.asarray(starts), np.asarray(ends), np.asarray(lengths)


def main():
    args = parse_args()
    if args.right <= args.left:
        raise SystemExit("--right must be greater than --left")
    if args.tract_length <= 0:
        raise SystemExit("--tract-length must be positive")

    csv_path = args.input
    if csv_path.is_dir():
        csv_path = csv_path / "arg_recombination_events.csv"
    if not csv_path.is_file():
        raise SystemExit(f"Input CSV not found: {csv_path}")

    starts, ends, lengths = read_gc_events(csv_path)
    if starts.size == 0:
        raise SystemExit(f"No GC events found in {csv_path}")

    tolerance = 0.01
    truncated = lengths < args.tract_length - tolerance
    invalid = (
        (starts < args.left)
        | (starts >= args.right)
        | (ends > args.right + tolerance)
        | (ends <= starts)
        | (lengths > args.tract_length + tolerance)
    )

    output_path = args.output or csv_path.with_name("arg_gc_diagnostics.png")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

    axes[0].hist(
        starts,
        bins=args.start_bins,
        range=(args.left, args.right),
        color="#3366aa",
        edgecolor="white",
        linewidth=0.4,
    )
    axes[0].set(
        xlabel="GC tract start (bp)",
        ylabel="Number of GC events",
        title="Observed GC start positions",
        xlim=(args.left, args.right),
    )

    length_bins = np.linspace(0.0, args.tract_length, 41)
    axes[1].hist(
        lengths,
        bins=length_bins,
        color="#228833",
        edgecolor="white",
        linewidth=0.4,
    )
    axes[1].set_yscale("log")
    axes[1].set(
        xlabel="Observed GC tract length (bp)",
        ylabel="Number of GC events (log scale)",
        title="Observed GC tract lengths",
        xlim=(0.0, args.tract_length),
    )
    axes[1].text(
        0.03,
        0.95,
        f"Full length: {starts.size - np.count_nonzero(truncated):,}\n"
        f"Boundary-truncated: {np.count_nonzero(truncated):,}",
        transform=axes[1].transAxes,
        va="top",
    )

    fig.suptitle(f"ARG gene-conversion diagnostics (n = {starts.size:,})")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    print(f"GC events: {starts.size:,}")
    print(f"Invalid GC intervals: {np.count_nonzero(invalid):,}")
    print(f"Full-length tracts: {starts.size - np.count_nonzero(truncated):,}")
    print(f"Boundary-truncated tracts: {np.count_nonzero(truncated):,}")
    print(f"Mean observed tract length: {np.mean(lengths):.6g} bp")
    print(f"Minimum observed tract length: {np.min(lengths):.6g} bp")
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
