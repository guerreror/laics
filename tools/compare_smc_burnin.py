#!/usr/bin/env python3
"""Compare early-position SMC TMRCA transients across tree-shape outputs.

The input files are smc_tree_shape.csv diagnostics. Each row is one SMC hop
state, so this script reports the transient three ways:

  1. by SMC hop number,
  2. by physical bp position, weighted by interval length,
  3. by fraction of the simulated region, also weighted by interval length.

The plateau is estimated from a late part of the chromosome, by default the
last half of the region. The reported crossing points are the first hop/bin
where mean root_time reaches a chosen fraction of that plateau.
"""

import argparse
import csv
import math
import textwrap
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_rows(path):
    rows = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            parsed = {}
            for key, value in row.items():
                if key in {"run", "hop"}:
                    parsed[key] = int(value)
                else:
                    parsed[key] = float(value)
            rows.append(parsed)
    return rows


def group_by_run(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["run"]].append(row)
    for run_rows in grouped.values():
        run_rows.sort(key=lambda row: row["hop"])
    return dict(sorted(grouped.items()))


def mean(values):
    values = [value for value in values if math.isfinite(value)]
    return sum(values) / len(values) if values else float("nan")


def median(values):
    values = sorted(value for value in values if math.isfinite(value))
    if not values:
        return float("nan")
    middle = len(values) // 2
    if len(values) % 2:
        return values[middle]
    return (values[middle - 1] + values[middle]) / 2.0


def quantile(values, probability):
    values = sorted(value for value in values if math.isfinite(value))
    if not values:
        return float("nan")
    if len(values) == 1:
        return values[0]
    index = probability * (len(values) - 1)
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return values[lower]
    fraction = index - lower
    return values[lower] * (1.0 - fraction) + values[upper] * fraction


def mean_by_hop(rows):
    by_hop = defaultdict(list)
    for row in rows:
        by_hop[row["hop"]].append(row["root_time"])
    return [
        {"hop": hop, "mean_root_time": mean(values), "n": len(values)}
        for hop, values in sorted(by_hop.items())
    ]


def root_time_by_hop_values(rows, max_hop):
    by_hop = defaultdict(list)
    for row in rows:
        if row["hop"] <= max_hop:
            by_hop[row["hop"]].append(row["root_time"])
    return dict(sorted(by_hop.items()))


def summarize_hop_distributions(label, rows, max_hop):
    summaries = []
    for hop, values in root_time_by_hop_values(rows, max_hop).items():
        summaries.append({
            "label": label,
            "hop": hop,
            "n": len(values),
            "mean_root_time": mean(values),
            "median_root_time": median(values),
            "q10_root_time": quantile(values, 0.10),
            "q25_root_time": quantile(values, 0.25),
            "q75_root_time": quantile(values, 0.75),
            "q90_root_time": quantile(values, 0.90),
            "min_root_time": min(values) if values else float("nan"),
            "max_root_time": max(values) if values else float("nan"),
        })
    return summaries


def plot_hop_ecdf(label, rows, max_hop, out_path, run_label=None):
    by_hop = root_time_by_hop_values(rows, max_hop)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    colors = ["#3b6f8f", "#b279a2", "#59a14f", "#e15759", "#9c755f"]
    for index, (hop, values) in enumerate(by_hop.items()):
        values = sorted(value for value in values if math.isfinite(value))
        if not values:
            continue
        y = [(i + 1) / len(values) for i in range(len(values))]
        ax.step(values, y, where="post", linewidth=1.7,
                color=colors[index % len(colors)], label=f"hop {hop}")
    ax.set_xlabel("Root time / TMRCA (generations)")
    ax.set_ylabel("Empirical cumulative probability")
    ax.set_title(label)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="lower right")
    if run_label:
        fig.text(
            0.5,
            0.015,
            textwrap.fill(run_label, width=110),
            ha="center",
            va="bottom",
            fontsize=9,
            color="#444444",
        )
        fig.tight_layout(rect=(0.0, 0.14, 1.0, 1.0))
    else:
        fig.tight_layout()
    fig.savefig(out_path, dpi=300, facecolor="white")
    plt.close(fig)


def weighted_bin_series(grouped, region_length, n_bins):
    bin_width = region_length / n_bins
    per_run = []

    for rows in grouped.values():
        weighted = [0.0] * n_bins
        covered = [0.0] * n_bins
        for index, row in enumerate(rows):
            start = max(0.0, row["current_x"])
            end = rows[index + 1]["current_x"] if index + 1 < len(rows) else region_length
            end = min(region_length, end)
            if end <= start:
                continue

            first_bin = max(0, min(n_bins - 1, int(start / bin_width)))
            last_bin = max(0, min(n_bins - 1, int(math.nextafter(end, start) / bin_width)))
            for bin_index in range(first_bin, last_bin + 1):
                bin_start = bin_index * bin_width
                bin_end = bin_start + bin_width
                overlap = max(0.0, min(end, bin_end) - max(start, bin_start))
                if overlap > 0.0:
                    weighted[bin_index] += row["root_time"] * overlap
                    covered[bin_index] += overlap

        per_run.append([
            weighted[i] / covered[i] if covered[i] > 0.0 else float("nan")
            for i in range(n_bins)
        ])

    series = []
    for bin_index in range(n_bins):
        values = [run[bin_index] for run in per_run if math.isfinite(run[bin_index])]
        bin_start = bin_index * bin_width
        bin_end = bin_start + bin_width
        series.append({
            "bin_index": bin_index,
            "bp_start": bin_start,
            "bp_end": bin_end,
            "bp_center": (bin_start + bin_end) / 2.0,
            "fraction_start": bin_start / region_length,
            "fraction_end": bin_end / region_length,
            "fraction_center": (bin_start + bin_end) / (2.0 * region_length),
            "mean_root_time": mean(values),
            "median_root_time": median(values),
            "n_runs_with_coverage": len(values),
        })
    return series


def first_crossing(series, x_key, y_key, threshold):
    for row in series:
        value = row[y_key]
        if math.isfinite(value) and value >= threshold:
            return row[x_key], value
    return float("nan"), float("nan")


def summarize_run(label, path, region_length, n_fraction_bins, n_bp_bins, threshold_fraction, plateau_start_fraction, max_distribution_hop):
    rows = read_rows(path)
    grouped = group_by_run(rows)
    last_rows = [run_rows[-1] for run_rows in grouped.values()]

    fraction_series = weighted_bin_series(grouped, region_length, n_fraction_bins)
    bp_series = weighted_bin_series(grouped, region_length, n_bp_bins)
    hop_series = mean_by_hop(rows)

    plateau_bins = [
        row["mean_root_time"]
        for row in fraction_series
        if row["fraction_center"] >= plateau_start_fraction
    ]
    plateau_mean = mean(plateau_bins)
    threshold = threshold_fraction * plateau_mean

    hop_cross, hop_cross_value = first_crossing(hop_series, "hop", "mean_root_time", threshold)
    bp_cross, bp_cross_value = first_crossing(bp_series, "bp_center", "mean_root_time", threshold)
    fraction_cross, fraction_cross_value = first_crossing(
        fraction_series, "fraction_center", "mean_root_time", threshold
    )

    summary = {
        "label": label,
        "path": str(path),
        "n_rows": len(rows),
        "n_runs": len(grouped),
        "mean_final_hop": mean(row["hop"] for row in last_rows),
        "median_final_hop": median(row["hop"] for row in last_rows),
        "mean_final_x": mean(row["current_x"] for row in last_rows),
        "median_final_x": median(row["current_x"] for row in last_rows),
        "hop0_mean_root_time": hop_series[0]["mean_root_time"] if hop_series else float("nan"),
        "hop1_mean_root_time": hop_series[1]["mean_root_time"] if len(hop_series) > 1 else float("nan"),
        "late_plateau_mean_root_time": plateau_mean,
        "threshold_fraction": threshold_fraction,
        "threshold_root_time": threshold,
        "first_hop_at_threshold": hop_cross,
        "first_hop_mean_root_time_at_threshold": hop_cross_value,
        "first_bp_center_at_threshold": bp_cross,
        "first_bp_mean_root_time_at_threshold": bp_cross_value,
        "first_fraction_center_at_threshold": fraction_cross,
        "first_fraction_mean_root_time_at_threshold": fraction_cross_value,
        "fraction_bin_count": n_fraction_bins,
        "fraction_bin_width_bp": region_length / n_fraction_bins,
        "bp_bin_count": n_bp_bins,
        "bp_bin_width_bp": region_length / n_bp_bins,
    }
    for hop in range(10):
        row = hop_series[hop] if hop < len(hop_series) and hop_series[hop]["hop"] == hop else None
        summary[f"hop{hop}_mean_root_time"] = row["mean_root_time"] if row else float("nan")
        summary[f"hop{hop}_n_runs"] = row["n"] if row else 0

    hop_distribution_summaries = summarize_hop_distributions(label, rows, max_distribution_hop)
    return summary, hop_series, bp_series, fraction_series, hop_distribution_summaries, rows


def parse_labeled_input(value):
    if "=" not in value:
        path = Path(value)
        return path.parent.name, path
    label, path = value.split("=", 1)
    return label, Path(path)


def write_dicts(path, rows):
    if not rows:
        return
    fieldnames = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", help="Inputs as LABEL=/path/to/smc_tree_shape.csv or just /path/to/csv")
    parser.add_argument("--region-length", type=float, required=True, help="Simulated region length in bp")
    parser.add_argument("--out-dir", default=None, help="Optional directory for comparison CSV outputs")
    parser.add_argument("--threshold-fraction", type=float, default=0.90, help="Fraction of late plateau used to call transient resolved")
    parser.add_argument("--plateau-start-fraction", type=float, default=0.50, help="Late-region fraction used to estimate plateau")
    parser.add_argument("--fraction-bins", type=int, default=100, help="Number of equal-fraction bins across the chromosome")
    parser.add_argument("--bp-bin-size", type=float, default=1000.0, help="Physical bin width in bp for bp-position crossing")
    parser.add_argument("--distribution-hops", type=int, default=2, help="Highest hop number to include in hop root-time distribution summaries/ECDFs")
    parser.add_argument("--run-label", default=None, help="Short parameter caption printed below each ECDF plot")
    args = parser.parse_args()

    n_bp_bins = max(1, int(math.ceil(args.region_length / args.bp_bin_size)))
    summaries = []
    per_bin_rows = []
    distribution_rows = []
    ecdf_inputs = []

    for input_value in args.inputs:
        label, path = parse_labeled_input(input_value)
        summary, hop_series, bp_series, fraction_series, hop_distribution_summaries, raw_rows = summarize_run(
            label,
            path,
            args.region_length,
            args.fraction_bins,
            n_bp_bins,
            args.threshold_fraction,
            args.plateau_start_fraction,
            args.distribution_hops,
        )
        summaries.append(summary)
        distribution_rows.extend(hop_distribution_summaries)
        ecdf_inputs.append((label, raw_rows))

        for row in hop_series:
            per_bin_rows.append({"label": label, "axis": "hop", **row})
        for row in bp_series:
            per_bin_rows.append({"label": label, "axis": "bp", **row})
        for row in fraction_series:
            per_bin_rows.append({"label": label, "axis": "fraction", **row})

    hop_headers = ",".join(f"hop{hop}_mean" for hop in range(10))
    print(f"label,n_runs,{hop_headers},plateau_mean,threshold,first_hop,first_bp,first_fraction,median_final_hops")
    for row in summaries:
        hop_values = ",".join(f"{row[f'hop{hop}_mean_root_time']:.6g}" for hop in range(10))
        print(
            f"{row['label']},{row['n_runs']},"
            f"{hop_values},"
            f"{row['late_plateau_mean_root_time']:.6g},{row['threshold_root_time']:.6g},"
            f"{row['first_hop_at_threshold']},{row['first_bp_center_at_threshold']:.6g},"
            f"{row['first_fraction_center_at_threshold']:.6g},{row['median_final_hop']:.6g}"
        )

    if args.out_dir:
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        write_dicts(out_dir / "smc_burnin_summary.csv", summaries)
        write_dicts(out_dir / "smc_burnin_series.csv", per_bin_rows)
        write_dicts(out_dir / "smc_hop_root_time_distributions.csv", distribution_rows)
        for label, raw_rows in ecdf_inputs:
            safe_label = "".join(char if char.isalnum() or char in "-_" else "_" for char in label)
            plot_hop_ecdf(
                label,
                raw_rows,
                args.distribution_hops,
                out_dir / f"{safe_label}_hop_root_time_ecdf.png",
                run_label=args.run_label,
            )
        print(f"Wrote {out_dir / 'smc_burnin_summary.csv'}")
        print(f"Wrote {out_dir / 'smc_burnin_series.csv'}")
        print(f"Wrote {out_dir / 'smc_hop_root_time_distributions.csv'}")


if __name__ == "__main__":
    main()
