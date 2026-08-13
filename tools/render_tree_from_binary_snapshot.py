#!/usr/bin/env python3
import argparse
import os
import re
import struct
import sys
from collections import defaultdict
from typing import Dict, List, Tuple

MAGIC = b"SMCTREE1"
RECORD = struct.Struct("<iiddQqdIH")


def clean_x(value: float) -> str:
    text = f"{value:.6f}".rstrip("0").rstrip(".")
    if not text:
        text = "0"
    return text.replace(".", "p").replace("-", "m")


def read_snapshots(path: str):
    rows_by_snapshot = defaultdict(list)
    with open(path, "rb") as f:
        magic = f.read(len(MAGIC))
        if magic != MAGIC:
            raise ValueError("Invalid snapshot binary header.")
        while True:
            chunk = f.read(RECORD.size)
            if not chunk:
                break
            if len(chunk) != RECORD.size:
                raise ValueError("Truncated snapshot binary record.")
            run, hop, x_start, x_end, node_id, parent_id, time, pop, inversion = RECORD.unpack(chunk)
            key = (run, hop, x_start, x_end)
            rows_by_snapshot[key].append({
                "node_id": node_id,
                "parent_id": parent_id,
                "time": time,
                "pop": pop,
                "inversion": inversion,
            })
    if not rows_by_snapshot:
        raise ValueError("No tree snapshots found.")
    return rows_by_snapshot


def select_run(snapshots, requested_run):
    runs = sorted({key[0] for key in snapshots})
    if requested_run is not None:
        if requested_run not in runs:
            raise ValueError(f"run {requested_run} not found. Available runs: {runs}")
        return requested_run
    if len(runs) == 1:
        return runs[0]
    print(f"Available runs: {runs[0]}..{runs[-1]} ({len(runs)} total)")
    text = input(f"Run to render [{runs[0]}]: ").strip()
    return runs[0] if text == "" else int(text)


def select_snapshot(snapshots, run, target_x):
    candidates = sorted([key for key in snapshots if key[0] == run], key=lambda k: (k[2], k[3], k[1]))
    if not candidates:
        raise ValueError(f"No snapshots found for run {run}.")

    eps = 1e-12
    exact = [key for key in candidates if abs(key[2] - target_x) <= eps and abs(key[3] - target_x) <= eps]
    if exact:
        return exact[0]

    containing = [key for key in candidates if key[2] + eps < target_x < key[3] - eps]
    if containing:
        key = containing[0]
        prev_key = None
        for candidate in candidates:
            if candidate[0] == run and candidate[1] == key[1] - 1:
                prev_key = candidate
                break
        options = []
        if prev_key is not None:
            options.append((f"tree at x_start={key[2]} before hop {key[1]}", prev_key))
        options.append((f"tree at x_end={key[3]} after hop {key[1]}", key))
        print(f"x={target_x} falls inside hop {key[1]} interval: {key[2]} -> {key[3]}")
        for i, (label, option_key) in enumerate(options, start=1):
            print(f"{i}. {label} (stored hop={option_key[1]})")
        choice = input(f"Choose 1-{len(options)}: ").strip()
        idx = int(choice) - 1
        if idx < 0 or idx >= len(options):
            raise ValueError("Invalid choice.")
        return options[idx][1]

    boundary = [key for key in candidates if abs(key[2] - target_x) <= eps or abs(key[3] - target_x) <= eps]
    if boundary:
        return boundary[-1]

    lower = None
    higher = None
    for key in candidates:
        if key[3] <= target_x:
            lower = key
        if key[2] >= target_x and higher is None:
            higher = key

    options = []
    if lower is not None:
        options.append(lower)
    if higher is not None and higher != lower:
        options.append(higher)
    if not options:
        return min(candidates, key=lambda k: min(abs(k[2] - target_x), abs(k[3] - target_x)))

    print(f"No snapshot interval contains x={target_x}.")
    for i, key in enumerate(options, start=1):
        print(f"{i}. hop={key[1]}, x_start={key[2]}, x_end={key[3]}")
    choice = input(f"Choose 1-{len(options)}: ").strip()
    idx = int(choice) - 1
    if idx < 0 or idx >= len(options):
        raise ValueError("Invalid choice.")
    return options[idx]


def build_nodes_edges(rows):
    nodes: Dict[int, dict] = {}
    edges: List[Tuple[int, int]] = []
    for row in rows:
        node_id = int(row["node_id"])
        parent_id = int(row["parent_id"])
        nodes[node_id] = {
            "id": node_id,
            "time": float(row["time"]),
            "pop": int(row["pop"]),
            "inv": int(row["inversion"]),
            "label": f"{node_id}\\nt={float(row['time']):.6g}\\npop={row['pop']} inv={row['inversion']}",
        }
        if parent_id >= 0:
            edges.append((parent_id, node_id))
    if not nodes:
        raise ValueError("Selected snapshot has no nodes.")
    return nodes, edges


def build_tskit(nodes: Dict[int, dict], edges: List[Tuple[int, int]], sequence_length: float):
    try:
        import tskit  # type: ignore
    except Exception as exc:
        raise RuntimeError("tskit is required. Install with: pip install tskit") from exc

    children = {child for _, child in edges}
    parents = {parent for parent, _ in edges}
    leaves = children - parents if edges else set(nodes)

    max_pop = max(node["pop"] for node in nodes.values())
    tables = tskit.TableCollection(sequence_length=sequence_length)
    tables.nodes.metadata_schema = tskit.MetadataSchema.permissive_json()
    for _ in range(max_pop + 1):
        tables.populations.add_row()

    eps = 1e-9
    changed = True
    while changed:
        changed = False
        for parent, child in edges:
            if nodes[parent]["time"] <= nodes[child]["time"]:
                nodes[parent]["time"] = nodes[child]["time"] + eps
                changed = True

    id_to_row = {}
    for node_id, node in nodes.items():
        flags = tskit.NODE_IS_SAMPLE if node_id in leaves else 0
        id_to_row[node_id] = tables.nodes.add_row(
            flags=flags,
            time=node["time"],
            population=node["pop"],
            metadata={"orig_id": node_id, "inv": node["inv"], "label": node["label"]},
        )

    for parent, child in edges:
        tables.edges.add_row(left=0.0, right=sequence_length, parent=id_to_row[parent], child=id_to_row[child])

    tables.sort()
    return tables.tree_sequence()


def render_svg(ts, out_svg: str, width: int, height: int, show_labels: bool):
    labels = None
    if show_labels:
        labels = {u: str(ts.node(u).metadata.get("orig_id", u)) for u in range(ts.num_nodes)}
    svg = ts.first().draw_svg(size=(width, height), node_labels=labels, y_axis=True, x_axis=True)
    svg = re.sub(r"<g class=\"background\">.*?</g>", "", svg, flags=re.DOTALL)
    svg = svg.replace("<svg ", "<svg style=\"background:#ffffff\" ", 1)
    svg = svg.replace(".axes line, .edge {stroke: black; fill: none}", ".axes line, .edge {stroke: #444; fill: none}")
    svg = svg.replace(".axes, .tree {font-size: 14px; text-anchor: middle}", ".axes, .tree {font-size: 16px; text-anchor: middle}")
    with open(out_svg, "w", encoding="utf-8") as f:
        f.write(svg)


def svg_to_png(svg_path: str, png_path: str):
    try:
        import cairosvg  # type: ignore
    except Exception as exc:
        raise RuntimeError("PNG export requires cairosvg. Install with: pip install cairosvg") from exc
    cairosvg.svg2png(url=svg_path, write_to=png_path)


def main() -> int:
    parser = argparse.ArgumentParser(description="Render one SMC tree from smc_tree_snapshots.bin.")
    parser.add_argument("--snapshot-bin", default="smc_tree_snapshots.bin")
    parser.add_argument("--target-x", type=float)
    parser.add_argument("--run", type=int)
    parser.add_argument("--sequence-length", type=float, default=1.0)
    parser.add_argument("--width", type=int, default=2400)
    parser.add_argument("--height", type=int, default=1600)
    parser.add_argument("--show-labels", action="store_true")
    args = parser.parse_args()

    if not os.path.isfile(args.snapshot_bin):
        print(f"Error: file not found: {args.snapshot_bin}", file=sys.stderr)
        return 1

    try:
        snapshots = read_snapshots(args.snapshot_bin)
        run = select_run(snapshots, args.run)
        target_x = args.target_x
        if target_x is None:
            target_x = float(input("Chromosome x position to render: ").strip())
        key = select_snapshot(snapshots, run, target_x)
        nodes, edges = build_nodes_edges(snapshots[key])
        ts = build_tskit(nodes, edges, args.sequence_length)
        base = f"selected_tree_run{key[0]}_hop{key[1]}_x{clean_x(target_x)}"
        out_svg = base + ".tskit.svg"
        out_png = base + ".tskit.png"
        render_svg(ts, out_svg, args.width, args.height, args.show_labels)
        svg_to_png(out_svg, out_png)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Selected run={key[0]}, hop={key[1]}, x_start={key[2]}, x_end={key[3]}")
    print(f"Wrote {out_svg}")
    print(f"Wrote {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
