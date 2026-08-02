#!/usr/bin/env python3
"""Plot SMC standard-neutral-model diagnostics from smc_tree_shape.csv.

This is meant for no-inversion control runs. It checks Rafael's expectation:
local TMRCA/root time should be around 2N and should not systematically grow
as the SMC moves along the chromosome.
"""

import argparse
import csv
import math
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
                    try:
                        parsed[key] = float(value)
                    except ValueError:
                        parsed[key] = value
            rows.append(parsed)
    return rows


def group_by_run(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["run"]].append(row)
    for run in grouped:
        grouped[run].sort(key=lambda row: row["hop"])
    return dict(sorted(grouped.items()))


def quantile(values, probability):
    values = sorted(values)
    if not values:
        return float("nan")
    if len(values) == 1:
        return values[0]
    index = probability * (len(values) - 1)
    lower = int(math.floor(index))
    upper = int(math.ceil(index))
    if lower == upper:
        return values[lower]
    fraction = index - lower
    return values[lower] * (1.0 - fraction) + values[upper] * fraction


def mean_by_hop(rows, fields):
    by_hop = defaultdict(list)
    for row in rows:
        by_hop[row["hop"]].append(row)

    means = []
    for hop in sorted(by_hop):
        out = {"hop": hop, "n": len(by_hop[hop])}
        for field in fields:
            vals = [row[field] for row in by_hop[hop] if field in row and math.isfinite(row[field])]
            out[field] = sum(vals) / len(vals) if vals else float("nan")
            out[f"{field}_q10"] = quantile(vals, 0.10)
            out[f"{field}_q90"] = quantile(vals, 0.90)
        means.append(out)
    return means


def add_run_label(fig, run_label):
    if run_label:
        fig.text(
            0.5,
            0.015,
            run_label,
            ha="center",
            va="bottom",
            fontsize=9,
            color="#444444",
        )


def finish_figure(fig, out_path, run_label=None):
    if run_label:
        fig.subplots_adjust(bottom=0.18)
    else:
        fig.tight_layout()
    fig.savefig(out_path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def nice_ceiling(value):
    if value <= 0.0:
        return value
    exponent = math.floor(math.log10(value))
    base = 10 ** exponent
    scaled = value / base
    if scaled <= 1:
        nice = 1
    elif scaled <= 2:
        nice = 2
    elif scaled <= 5:
        nice = 5
    else:
        nice = 10
    return nice * base


def choose_bin_count(region_length, bin_size_bp=None):
    if bin_size_bp and bin_size_bp > 0.0:
        return max(1, int(math.ceil(region_length / bin_size_bp)))

    target_bins = 200
    min_bin_width = 1000.0
    raw_width = max(min_bin_width, region_length / target_bins)
    bin_width = nice_ceiling(raw_width)
    return max(1, int(math.ceil(region_length / bin_width)))


def position_weighted_bins(grouped, field, region_length, bin_size_bp=None):
    if not region_length or region_length <= 0.0:
        return [], [], [], []

    n_bins = choose_bin_count(region_length, bin_size_bp)
    bin_width = region_length / n_bins
    per_run = []
    for rows in grouped.values():
        weighted_sum = [0.0] * n_bins
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
                    weighted_sum[bin_index] += row[field] * overlap
                    covered[bin_index] += overlap
        per_run.append([
            weighted_sum[i] / covered[i] if covered[i] > 0.0 else float("nan")
            for i in range(n_bins)
        ])

    centers = [(i + 0.5) * bin_width for i in range(n_bins)]
    means, lows, highs = [], [], []
    for bin_index in range(n_bins):
        values = [run[bin_index] for run in per_run if math.isfinite(run[bin_index])]
        means.append(sum(values) / len(values) if values else float("nan"))
        lows.append(quantile(values, 0.10))
        highs.append(quantile(values, 0.90))
    return centers, means, lows, highs


def write_summary(rows, out_path, reference_time):
    fields = [
        "run",
        "first_hop",
        "last_hop",
        "first_x",
        "last_x",
        "first_root_time",
        "last_root_time",
        "first_standard_branch_length",
        "last_standard_branch_length",
        "first_rho",
        "last_rho",
        "last_unary_parent_branch_length",
        "last_root_time_over_reference",
    ]
    grouped = group_by_run(rows)
    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for run, run_rows in grouped.items():
            first = run_rows[0]
            last = run_rows[-1]
            writer.writerow(
                {
                    "run": run,
                    "first_hop": first["hop"],
                    "last_hop": last["hop"],
                    "first_x": first["current_x"],
                    "last_x": last["current_x"],
                    "first_root_time": first["root_time"],
                    "last_root_time": last["root_time"],
                    "first_standard_branch_length": first.get("standard_branch_length", float("nan")),
                    "last_standard_branch_length": last.get("standard_branch_length", float("nan")),
                    "first_rho": first["rho"],
                    "last_rho": last["rho"],
                    "last_unary_parent_branch_length": last.get("unary_parent_branch_length", float("nan")),
                    "last_root_time_over_reference": last["root_time"] / reference_time if reference_time else float("nan"),
                }
            )


def reference_label(ref_y, ref_label):
    if ref_label:
        return f"{ref_label} ({ref_y:g})"
    return f"Reference ({ref_y:g})"


def plot_by_hop(grouped, means, field, ylabel, out_path, ref_y=None, ref_label=None, log_y=False, run_label=None):
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    add_run_label(fig, run_label)
    min_replicates = max(2, math.ceil(len(grouped) / 2))
    supported = [row for row in means if row["n"] >= min_replicates]
    hops = [row["hop"] for row in supported]
    ax.fill_between(
        hops,
        [row[f"{field}_q10"] for row in supported],
        [row[f"{field}_q90"] for row in supported],
        color="#1f4e79",
        alpha=0.14,
        linewidth=0,
        label="Variation among replicates (10–90%)",
    )
    ax.plot(
        hops,
        [row[field] for row in supported],
        color="#1f4e79",
        linewidth=2.0,
        label=f"Mean across replicates (≥{min_replicates} runs)",
    )
    if ref_y is not None:
        ax.axhline(ref_y, color="#b13b2e", linestyle="--", linewidth=1.4, label=reference_label(ref_y, ref_label))

    ax.set_xlabel("SMC hop")
    ax.set_ylabel(ylabel)
    if log_y:
        ax.set_yscale("log")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="upper left")
    finish_figure(fig, out_path, run_label)


def plot_by_x(grouped, field, ylabel, out_path, ref_y=None, ref_label=None, log_y=False, region_length=None, run_label=None, bin_size_bp=None):
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    add_run_label(fig, run_label)
    centers, bin_means, bin_lows, bin_highs = position_weighted_bins(
        grouped, field, region_length, bin_size_bp
    )
    if centers:
        bin_width = region_length / len(centers)
        bin_label = f"{bin_width / 1000:g}-kb bins" if bin_width >= 1000 else f"{bin_width:g}-bp bins"
        ax.fill_between(
            centers, bin_lows, bin_highs, color="#1f4e79", alpha=0.14,
            linewidth=0, label="Variation among replicates (10–90%)"
        )
        ax.plot(
            centers, bin_means, color="#1f4e79", linewidth=2.0,
            label=f"Mean TMRCA across replicates ({bin_label})"
        )

    if ref_y is not None:
        ax.axhline(ref_y, color="#b13b2e", linestyle="--", linewidth=1.4, label=reference_label(ref_y, ref_label))
    ax.legend(frameon=False, loc="upper left")

    ax.set_xlabel("Current position (bp, scientific notation)")
    ax.set_ylabel(ylabel)
    if log_y:
        ax.set_yscale("log")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    finish_figure(fig, out_path, run_label)


def plot_by_fraction_traversed(grouped, field, ylabel, out_path, ref_y=None, ref_label=None, log_y=False, region_length=None, run_label=None):
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    add_run_label(fig, run_label)
    # Fraction plots should compare the same relative location across runs, so
    # keep bins fixed at 1% of the region instead of adapting to bp length.
    fraction_bin_size_bp = region_length / 100.0 if region_length else None
    centers, bin_means, bin_lows, bin_highs = position_weighted_bins(
        grouped, field, region_length, fraction_bin_size_bp
    )
    if centers:
        frac_centers = [center / region_length for center in centers]
        bin_width = region_length / len(centers)
        bin_label = f"{bin_width / 1000:g}-kb bins" if bin_width >= 1000 else f"{bin_width:g}-bp bins"
        ax.fill_between(
            frac_centers, bin_lows, bin_highs, color="#1f4e79", alpha=0.14,
            linewidth=0, label="Variation among replicates (10–90%)"
        )
        ax.plot(
            frac_centers, bin_means, color="#1f4e79", linewidth=2.0,
            label=f"Mean TMRCA across replicates ({bin_label})"
        )

    if ref_y is not None:
        ax.axhline(ref_y, color="#b13b2e", linestyle="--", linewidth=1.4, label=reference_label(ref_y, ref_label))
    ax.legend(frameon=False, loc="best")

    ax.set_xlabel("Fraction of region traversed")
    ax.set_ylabel(ylabel)
    if log_y:
        ax.set_yscale("log")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    finish_figure(fig, out_path, run_label)


def plot_example_by_x(grouped, field, ylabel, out_path, ref_y=None, ref_label=None, region_length=None, run_label=None):
    run, rows = next(iter(grouped.items()))
    xs = [row["current_x"] for row in rows]
    ys = [row[field] for row in rows]
    if region_length and region_length > xs[-1]:
        xs.append(region_length)
        ys.append(ys[-1])

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    add_run_label(fig, run_label)
    ax.step(xs, ys, where="post", color="#3b6f8f", linewidth=1.15,
            label=f"Replicate {run} local TMRCA")
    if ref_y is not None:
        ax.axhline(ref_y, color="#b13b2e", linestyle="--", linewidth=1.4,
                   label=reference_label(ref_y, ref_label))
    ax.set_xlabel("Current position (bp, scientific notation)")
    ax.set_ylabel(ylabel)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="best")
    finish_figure(fig, out_path, run_label)


def plot_unary_fraction_by_hop(rows, out_path, run_label=None):
    by_hop = defaultdict(list)
    for row in rows:
        total = row.get("total_branch_length", float("nan"))
        unary = row.get("internal_unary_chain_branch_length", row.get("unary_parent_branch_length", float("nan")))
        if total and total > 0.0 and math.isfinite(total) and math.isfinite(unary):
            by_hop[row["hop"]].append(unary / total)

    hops = sorted(by_hop)
    means = [sum(by_hop[hop]) / len(by_hop[hop]) for hop in hops]

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    add_run_label(fig, run_label)
    ax.plot(hops, means, color="#7b2f5c", linewidth=2.2)
    ax.set_xlabel("SMC hop")
    ax.set_ylabel("Internal unary-chain branch length / total branch length")
    ax.set_ylim(0, 1.02)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    finish_figure(fig, out_path, run_label)


def plot_branch_length_composition(rows, out_path, run_label=None):
    exclusive_fields = [
        ("terminal_from_unary_parent_branch_length", "Terminal from unary parent", "#9c755f"),
        ("terminal_from_branching_parent_branch_length", "Terminal from branching parent", "#4c78a8"),
        ("internal_unary_chain_branch_length", "Internal unary chain", "#b279a2"),
        ("internal_branching_exclusive_branch_length", "Internal branching", "#59a14f"),
    ]
    fallback_fields = [
        ("terminal_leaf_branch_length", "Terminal sample branches", "#4c78a8"),
        ("unary_parent_branch_length", "Unary/event-history branches", "#b279a2"),
        ("internal_branching_branch_length", "Internal branching branches", "#59a14f"),
    ]
    fields = exclusive_fields if exclusive_fields[0][0] in rows[0] else fallback_fields
    by_hop = defaultdict(list)
    for row in rows:
        by_hop[row["hop"]].append(row)

    hops = sorted(by_hop)
    series = []
    for field, _, _ in fields:
        vals = []
        for hop in hops:
            hop_vals = [row.get(field, float("nan")) for row in by_hop[hop]]
            hop_vals = [v for v in hop_vals if math.isfinite(v)]
            vals.append(sum(hop_vals) / len(hop_vals) if hop_vals else float("nan"))
        series.append(vals)

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    add_run_label(fig, run_label)
    ax.stackplot(
        hops,
        series,
        labels=[label for _, label, _ in fields],
        colors=[color for _, _, color in fields],
        alpha=0.86,
    )
    ax.set_xlabel("SMC hop")
    ax.set_ylabel("Mean active-tree branch length")
    ax.set_yscale("log")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="upper left")
    finish_figure(fig, out_path, run_label)


def plot_hop_locations(grouped, out_path, region_length=None, run_label=None):
    hop_xs = []
    hops_per_run = []
    for run, rows in grouped.items():
        run_hops = [row["current_x"] for row in rows if row["hop"] > 0]
        hop_xs.extend(run_hops)
        hops_per_run.append(len(run_hops))

    median_hops = quantile(hops_per_run, 0.50) if hops_per_run else float("nan")
    mean_hops = sum(hops_per_run) / len(hops_per_run) if hops_per_run else float("nan")
    example_run, example_rows = next(iter(grouped.items()))
    example_hops = [row["current_x"] for row in example_rows if row["hop"] > 0]

    fig, (ax_hist, ax_rug) = plt.subplots(
        2, 1, figsize=(7.2, 5.4), sharex=True,
        gridspec_kw={"height_ratios": [2.0, 1.15], "hspace": 0.08},
    )
    add_run_label(fig, run_label)

    upper_x = region_length if region_length and region_length > 0.0 else max(hop_xs, default=1.0)
    bins = min(80, max(10, int(math.sqrt(max(len(hop_xs), 1)))))
    ax_hist.hist(
        hop_xs,
        bins=bins,
        range=(0.0, upper_x),
        color="#3b6f8f",
        edgecolor="white",
        linewidth=0.45,
    )
    ax_hist.text(
        0.99,
        0.95,
        f"Mean/run: {mean_hops:.1f}\nMedian/run: {median_hops:.1f}",
        transform=ax_hist.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "none", "alpha": 0.88},
    )
    ax_hist.set_ylabel("Hop count")
    ax_hist.spines["top"].set_visible(False)
    ax_hist.spines["right"].set_visible(False)

    ax_rug.eventplot(
        example_hops,
        orientation="horizontal",
        lineoffsets=0,
        linelengths=0.75,
        linewidths=0.9,
        colors="#7b2f5c",
    )
    ax_rug.text(
        0.99,
        0.86,
        f"Example replicate {example_run}: {len(example_hops)} hops",
        transform=ax_rug.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "none", "alpha": 0.88},
    )
    ax_rug.set_yticks([])
    ax_rug.set_ylabel(f"Run {example_run}")
    ax_rug.set_xlabel("Hop location along region (bp)")
    ax_rug.spines["top"].set_visible(False)
    ax_rug.spines["right"].set_visible(False)
    ax_rug.set_xlim(0.0, upper_x)

    finish_figure(fig, out_path, run_label)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", help="SMC output directory containing smc_tree_shape.csv")
    parser.add_argument("--out-dir", default=None, help="Directory for plots; default: run_dir/snm_diagnostics")
    parser.add_argument("--pop-size", type=float, default=10000.0, help="Deprecated: retained for old commands; reference lines now require --reference-time")
    parser.add_argument("--reference-time", type=float, default=None, help="Expected TMRCA/root-time reference line in generations")
    parser.add_argument("--reference-label", default=None, help="Label for the reference line, e.g. 'N haploid expectation'")
    parser.add_argument("--prefix", default="smc_snm", help="Output filename prefix")
    parser.add_argument("--region-length", type=float, default=None, help="Region length in bp for fraction-traversed plots")
    parser.add_argument("--bin-size-bp", type=float, default=None, help="Optional fixed bp bin width for position-weighted plots; default chooses an adaptive width from region length")
    parser.add_argument("--run-label", default=None, help="Short context label/caption printed at the top of each plot")
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    input_path = run_dir / "smc_tree_shape.csv"
    if not input_path.exists():
        raise SystemExit(f"Missing {input_path}")

    out_dir = Path(args.out_dir) if args.out_dir else run_dir / "snm_diagnostics"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = read_rows(input_path)
    grouped = group_by_run(rows)
    fields = ["root_time", "standard_branch_length", "unary_parent_branch_length", "rho"]
    means = mean_by_hop(rows, fields)
    reference_time = args.reference_time
    reference_label_text = args.reference_label or "Reference"

    write_summary(rows, out_dir / f"{args.prefix}_summary.csv", reference_time)
    plot_by_hop(grouped, means, "root_time", "Local TMRCA / root time (generations)", out_dir / f"{args.prefix}_root_time_by_hop.png", ref_y=reference_time, ref_label=reference_label_text, run_label=args.run_label)
    plot_by_x(grouped, "root_time", "Local TMRCA / root time (generations)", out_dir / f"{args.prefix}_root_time_by_x.png", ref_y=reference_time, ref_label=reference_label_text, region_length=args.region_length, run_label=args.run_label, bin_size_bp=args.bin_size_bp)
    plot_by_fraction_traversed(grouped, "root_time", "Local TMRCA / root time (generations)", out_dir / f"{args.prefix}_root_time_by_fraction_traversed.png", ref_y=reference_time, ref_label=reference_label_text, region_length=args.region_length, run_label=args.run_label)
    plot_example_by_x(grouped, "root_time", "Local TMRCA / root time (generations)", out_dir / f"{args.prefix}_root_time_example_run0_by_x.png", ref_y=reference_time, ref_label=reference_label_text, region_length=args.region_length, run_label=args.run_label)
    plot_by_hop(grouped, means, "standard_branch_length", "Standard active-tree branch length", out_dir / f"{args.prefix}_standard_branch_length_by_hop.png", log_y=True, run_label=args.run_label)
    plot_by_hop(grouped, means, "unary_parent_branch_length", "Unary-parent branch length", out_dir / f"{args.prefix}_unary_branch_length_by_hop.png", log_y=True, run_label=args.run_label)
    plot_by_hop(grouped, means, "rho", "Horizontal rate rho", out_dir / f"{args.prefix}_rho_by_hop.png", log_y=True, run_label=args.run_label)
    plot_unary_fraction_by_hop(rows, out_dir / f"{args.prefix}_unary_fraction_by_hop.png", run_label=args.run_label)
    plot_branch_length_composition(rows, out_dir / f"{args.prefix}_branch_length_composition_by_hop.png", run_label=args.run_label)
    plot_hop_locations(grouped, out_dir / f"{args.prefix}_hop_locations.png", region_length=args.region_length, run_label=args.run_label)

    print(f"Wrote SNM diagnostics to {out_dir}")


if __name__ == "__main__":
    main()
