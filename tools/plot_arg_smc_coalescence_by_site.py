#!/usr/bin/env python3
"""Compare ARG and SMC pairwise coalescence-time distributions by position."""

import argparse
import bisect
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, StrMethodFormatter


ARG_COLOR = "#2A6FBB"
SMC_COLOR = "#E66B4E"


def quantile(values, probability):
    ordered = sorted(values)
    index = probability * (len(ordered) - 1)
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return ordered[lower]
    fraction = index - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def summarize(values):
    n = len(values)
    mean = statistics.mean(values)
    sd = statistics.stdev(values)
    se = sd / math.sqrt(n)
    return {
        "n": n,
        "mean": mean,
        "sd": sd,
        "se": se,
        "ci_low": mean - 1.96 * se,
        "ci_high": mean + 1.96 * se,
        "min": min(values),
        "q25": quantile(values, 0.25),
        "median": quantile(values, 0.50),
        "q75": quantile(values, 0.75),
        "max": max(values),
    }


def read_arg_by_position(path, bases_per_morgan, total_population_size):
    values = defaultdict(list)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            position_bp = int(round(float(row["site_position_morgan"]) * bases_per_morgan))
            values[position_bp].append(
                float(row["tmrca_scaled"]) * total_population_size
            )
    return dict(sorted(values.items()))


def read_smc_by_position(path, positions):
    by_run = defaultdict(list)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"run", "current_x", "root_time"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "SMC input must be smc_tree_shape.csv; "
                f"missing: {', '.join(sorted(missing))}"
            )
        for row in reader:
            by_run[int(row["run"])].append(
                (
                    float(row["current_x"]),
                    float(row["root_time"]),
                )
            )

    values = {position: [] for position in positions}
    for run in sorted(by_run):
        rows = sorted(by_run[run], key=lambda row: row[0])
        xs = [row[0] for row in rows]
        for position in positions:
            index = bisect.bisect_right(xs, position) - 1
            if index < 0:
                raise ValueError(
                    f"SMC run {run} has no tree at or before position {position}"
                )
            _, root_time = rows[index]
            values[position].append(root_time)
    return values


def write_summary(path, positions, arg_values, smc_values):
    fields = [
        "position_bp", "simulator", "n", "mean_generations", "sd_generations",
        "se_generations", "ci_low_generations", "ci_high_generations",
        "min_generations", "q25_generations", "median_generations",
        "q75_generations", "max_generations",
    ]
    summaries = {}
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for position in positions:
            for simulator, source in (("ARG", arg_values), ("SMC", smc_values)):
                summary = summarize(source[position])
                summaries[(position, simulator)] = summary
                writer.writerow({
                    "position_bp": position,
                    "simulator": simulator,
                    "n": summary["n"],
                    **{
                        f"{key}_generations": summary[key]
                        for key in (
                            "mean", "sd", "se", "ci_low", "ci_high", "min",
                            "q25", "median", "q75", "max",
                        )
                    },
                })
    return summaries


def distribution_panel(ax, position, arg, smc, panel_label, x_limit):
    bins = 42
    ax.hist(arg, bins=bins, range=(0, x_limit), density=True, histtype="stepfilled",
            alpha=0.28, color=ARG_COLOR, edgecolor=ARG_COLOR, linewidth=1.4,
            label=f"ARG (n={len(arg):,})")
    ax.hist(smc, bins=bins, range=(0, x_limit), density=True, histtype="step",
            color=SMC_COLOR, linewidth=2.0, label=f"SMC (n={len(smc):,})")
    ax.set_title(f"{panel_label}  Position = {position:,} bp",
                 loc="left", fontweight="bold")
    ax.set_xlabel("Coalescence time (generations)")
    ax.set_ylabel("Density")
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.xaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
    ax.grid(axis="y", alpha=0.2)


def make_figure(path_base, positions, arg_values, smc_values, summaries,
                region_length, analytical_expectation, title,
                show_analytical_expectation=True,
                mean_panel_title="Expected coalescence time",
                density_positions=None):
    if density_positions is None:
        density_positions = [positions[0], positions[1], positions[-1]]
    density_values = [
        value
        for position in density_positions
        for source in (arg_values, smc_values)
        for value in source[position]
    ]
    density_x_limit = quantile(density_values, 0.995)

    fig = plt.figure(figsize=(16, 10), constrained_layout=True)
    grid = fig.add_gridspec(2, 6, height_ratios=[1, 1.18])
    density_axes = [
        fig.add_subplot(grid[0, 0:2]),
        fig.add_subplot(grid[0, 2:4]),
        fig.add_subplot(grid[0, 4:6]),
    ]
    for ax, position, label in zip(density_axes, density_positions, ("(a)", "(b)", "(c)")):
        distribution_panel(
            ax, position, arg_values[position], smc_values[position],
            label, density_x_limit,
        )
    density_axes[0].legend(frameon=False)

    box_ax = fig.add_subplot(grid[1, 0:4])
    centers = list(range(len(positions)))
    offsets = (-0.18, 0.18)
    widths = 0.30
    for offset, simulator, source, color in (
        (offsets[0], "ARG", arg_values, ARG_COLOR),
        (offsets[1], "SMC", smc_values, SMC_COLOR),
    ):
        artists = box_ax.boxplot(
            [source[position] for position in positions],
            positions=[center + offset for center in centers],
            widths=widths,
            whis=(0, 100),
            showmeans=True,
            showfliers=False,
            patch_artist=True,
            meanprops={"marker": "D", "markerfacecolor": "white", "markeredgecolor": "black", "markersize": 4},
            medianprops={"color": "black", "linewidth": 1.5},
            boxprops={"facecolor": color, "edgecolor": color, "alpha": 0.58},
            whiskerprops={"color": color},
            capprops={"color": color},
        )
        artists["boxes"][0].set_label(
            f"{simulator} (n={len(source[positions[0]]):,})"
        )
    box_ax.set_title("(d)  Distribution summaries", loc="left", fontweight="bold")
    box_ax.set_xticks(centers, [f"{position:,}" for position in positions])
    box_ax.set_xlabel("Genomic position (bp)")
    box_ax.set_ylabel("Coalescence time (generations)")
    box_ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
    box_ax.grid(axis="y", alpha=0.2)
    box_ax.legend(frameon=False, loc="upper right")
    box_ax.text(
        0.01, 0.98, "Boxes: quartiles and median; whiskers: min–max; diamonds: means",
        transform=box_ax.transAxes, va="top", fontsize=9, color="#444444",
    )

    mean_ax = fig.add_subplot(grid[1, 4:6])
    for simulator, color, marker in (("ARG", ARG_COLOR, "o"), ("SMC", SMC_COLOR, "x")):
        means = [summaries[(position, simulator)]["mean"] for position in positions]
        errors = [1.96 * summaries[(position, simulator)]["se"] for position in positions]
        mean_ax.errorbar(
            positions, means, yerr=errors, color=color, marker=marker,
            linewidth=1.6, capsize=3,
            label=f"{simulator} (n={summaries[(positions[0], simulator)]['n']:,})",
        )
    if show_analytical_expectation:
        mean_ax.axhline(analytical_expectation, color="#666666", linewidth=1.6,
                        label=f"Analytical expectation ({analytical_expectation:,.0f})")
    mean_ax.set_title(f"(e)  {mean_panel_title}", loc="left", fontweight="bold")
    mean_ax.set_xlabel("Genomic position (bp)")
    mean_ax.xaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
    mean_ax.set_ylabel("Mean coalescence time (generations)")
    mean_ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
    mean_ax.grid(alpha=0.2)
    mean_ax.legend(frameon=False, fontsize=9)

    fig.suptitle(
        title,
        fontsize=18,
        fontweight="bold",
    )
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(path_base.with_suffix(f".{suffix}"), dpi=300, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("arg_tmrca_by_site", type=Path)
    parser.add_argument("smc_tree_shape", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument(
        "--smc-x0-tree-shape",
        type=Path,
        help=(
            "Optional smc_tree_shape.csv used only for the exact x=0 tree; "
            "positive positions continue to use the primary SMC input."
        ),
    )
    parser.add_argument("--bases-per-morgan", type=float, default=1.0e8)
    parser.add_argument("--total-population-size", type=float, default=10833.0)
    parser.add_argument("--region-length", type=float, default=100000.0)
    parser.add_argument("--analytical-expectation", type=float, default=150000.0)
    parser.add_argument("--hide-analytical-expectation", action="store_true")
    parser.add_argument("--mean-panel-title", default="Expected coalescence time")
    parser.add_argument(
        "--density-positions",
        type=lambda value: [int(item) for item in value.split(",")],
        help="Comma-separated genomic positions for the three density panels.",
    )
    parser.add_argument(
        "--positions",
        type=lambda value: [int(item) for item in value.split(",")],
        help="Comma-separated genomic positions to include in the comparison.",
    )
    parser.add_argument("--prefix", default="arg_vs_smc_nomigration_coalescence")
    parser.add_argument(
        "--title", default="ARG versus SMC coalescence times — one pop SNM"
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    arg_values = read_arg_by_position(
        args.arg_tmrca_by_site,
        args.bases_per_morgan,
        args.total_population_size,
    )
    if args.positions is not None:
        missing = [position for position in args.positions if position not in arg_values]
        if missing:
            raise SystemExit(
                "Requested positions absent from ARG data: "
                + ", ".join(str(position) for position in missing)
            )
        arg_values = {position: arg_values[position] for position in args.positions}
    positions = sorted(arg_values)
    if args.density_positions is not None:
        if len(args.density_positions) != 3:
            raise SystemExit("--density-positions requires exactly three positions")
        missing = [
            position for position in args.density_positions
            if position not in arg_values
        ]
        if missing:
            raise SystemExit(
                "Density-panel positions absent from ARG data: "
                + ", ".join(str(position) for position in missing)
            )
    smc_values = read_smc_by_position(args.smc_tree_shape, positions)
    if args.smc_x0_tree_shape is not None and 0 in positions:
        smc_values[0] = read_smc_by_position(args.smc_x0_tree_shape, [0])[0]
    summary_path = args.output_dir / f"{args.prefix}_summary.csv"
    summaries = write_summary(
        summary_path, positions, arg_values, smc_values
    )
    make_figure(
        args.output_dir / args.prefix,
        positions,
        arg_values,
        smc_values,
        summaries,
        args.region_length,
        args.analytical_expectation,
        args.title,
        not args.hide_analytical_expectation,
        args.mean_panel_title,
        args.density_positions,
    )
    print(f"Wrote {summary_path}")
    for suffix in ("png", "pdf", "svg"):
        print(f"Wrote {args.output_dir / (args.prefix + '.' + suffix)}")


if __name__ == "__main__":
    main()
