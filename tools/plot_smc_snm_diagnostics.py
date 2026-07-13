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


def mean_by_hop(rows, fields):
    by_hop = defaultdict(list)
    for row in rows:
        by_hop[row["hop"]].append(row)

    means = []
    for hop in sorted(by_hop):
        out = {"hop": hop}
        for field in fields:
            vals = [row[field] for row in by_hop[hop] if field in row and math.isfinite(row[field])]
            out[field] = sum(vals) / len(vals) if vals else float("nan")
        means.append(out)
    return means


def write_summary(rows, out_path, pop_size):
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
        "last_root_time_over_2N",
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
                    "last_root_time_over_2N": last["root_time"] / (2.0 * pop_size),
                }
            )


def plot_by_hop(grouped, means, field, ylabel, out_path, ref_y=None, log_y=False):
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for _, rows in grouped.items():
        ax.plot([r["hop"] for r in rows], [r[field] for r in rows], color="#8aa0b8", alpha=0.18, linewidth=0.8)

    ax.plot([r["hop"] for r in means], [r[field] for r in means], color="#1f4e79", linewidth=2.2, label="Mean across replicates")
    if ref_y is not None:
        ax.axhline(ref_y, color="#b13b2e", linestyle="--", linewidth=1.4, label=f"2N reference ({ref_y:g})")

    ax.set_xlabel("SMC hop")
    ax.set_ylabel(ylabel)
    if log_y:
        ax.set_yscale("log")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_by_x(grouped, field, ylabel, out_path, ref_y=None, log_y=False):
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for _, rows in grouped.items():
        ax.plot([r["current_x"] for r in rows], [r[field] for r in rows], color="#3b6f8f", alpha=0.22, linewidth=0.8)

    if ref_y is not None:
        ax.axhline(ref_y, color="#b13b2e", linestyle="--", linewidth=1.4, label=f"2N reference ({ref_y:g})")
        ax.legend(frameon=False, loc="best")

    ax.set_xlabel("Current position (bp, scientific notation)")
    ax.set_ylabel(ylabel)
    if log_y:
        ax.set_yscale("log")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_by_fraction_traversed(grouped, field, ylabel, out_path, ref_y=None, log_y=False, region_length=None):
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for _, rows in grouped.items():
        xs = [r["current_x"] for r in rows]
        run_start = min(xs) if xs else 0.0
        denom = region_length if region_length and region_length > 0.0 else (max(xs) - run_start)
        if denom <= 0.0:
            frac = [0.0 for _ in rows]
        else:
            frac = [(x - run_start) / denom for x in xs]
        ax.plot(frac, [r[field] for r in rows], color="#3b6f8f", alpha=0.22, linewidth=0.8)

    if ref_y is not None:
        ax.axhline(ref_y, color="#b13b2e", linestyle="--", linewidth=1.4, label=f"2N reference ({ref_y:g})")
        ax.legend(frameon=False, loc="best")

    ax.set_xlabel("Fraction of region traversed")
    ax.set_ylabel(ylabel)
    if log_y:
        ax.set_yscale("log")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_unary_fraction_by_hop(rows, out_path):
    by_hop = defaultdict(list)
    for row in rows:
        total = row.get("total_branch_length", float("nan"))
        unary = row.get("internal_unary_chain_branch_length", row.get("unary_parent_branch_length", float("nan")))
        if total and total > 0.0 and math.isfinite(total) and math.isfinite(unary):
            by_hop[row["hop"]].append(unary / total)

    hops = sorted(by_hop)
    means = [sum(by_hop[hop]) / len(by_hop[hop]) for hop in hops]

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(hops, means, color="#7b2f5c", linewidth=2.2)
    ax.set_xlabel("SMC hop")
    ax.set_ylabel("Internal unary-chain branch length / total branch length")
    ax.set_ylim(0, 1.02)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_branch_length_composition(rows, out_path):
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
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", help="SMC output directory containing smc_tree_shape.csv")
    parser.add_argument("--out-dir", default=None, help="Directory for plots; default: run_dir/snm_diagnostics")
    parser.add_argument("--pop-size", type=float, default=10000.0, help="Population size N for the 2N reference line")
    parser.add_argument("--prefix", default="smc_snm", help="Output filename prefix")
    parser.add_argument("--region-length", type=float, default=None, help="Region length in bp for fraction-traversed plots")
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
    two_n = 2.0 * args.pop_size

    write_summary(rows, out_dir / f"{args.prefix}_summary.csv", args.pop_size)
    plot_by_hop(grouped, means, "root_time", "Local TMRCA / root time (generations)", out_dir / f"{args.prefix}_root_time_by_hop.png", ref_y=two_n)
    plot_by_x(grouped, "root_time", "Local TMRCA / root time (generations)", out_dir / f"{args.prefix}_root_time_by_x.png", ref_y=two_n)
    plot_by_fraction_traversed(grouped, "root_time", "Local TMRCA / root time (generations)", out_dir / f"{args.prefix}_root_time_by_fraction_traversed.png", ref_y=two_n, region_length=args.region_length)
    plot_by_hop(grouped, means, "standard_branch_length", "Standard active-tree branch length", out_dir / f"{args.prefix}_standard_branch_length_by_hop.png", log_y=True)
    plot_by_hop(grouped, means, "unary_parent_branch_length", "Unary-parent branch length", out_dir / f"{args.prefix}_unary_branch_length_by_hop.png", log_y=True)
    plot_by_hop(grouped, means, "rho", "Horizontal rate rho", out_dir / f"{args.prefix}_rho_by_hop.png", log_y=True)
    plot_unary_fraction_by_hop(rows, out_dir / f"{args.prefix}_unary_fraction_by_hop.png")
    plot_branch_length_composition(rows, out_dir / f"{args.prefix}_branch_length_composition_by_hop.png")

    print(f"Wrote SNM diagnostics to {out_dir}")


if __name__ == "__main__":
    main()
