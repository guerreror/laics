#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import re
import sys
from typing import Dict, List

import cairosvg
import tskit


def read_nodes(path: str) -> Dict[int, dict]:
    nodes = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            nid = int(row["id"])
            nodes[nid] = {
                "time": float(row["time"]),
                "pop": int(row["pop"]),
                "inv": int(row["inv"]),
                "is_sample": int(row["is_sample"]),
            }
    return nodes


def read_edges(path: str) -> List[dict]:
    edges = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            edges.append({
                "parent": int(row["parent"]),
                "child": int(row["child"]),
                "left": float(row["left"]),
                "right": float(row["right"]),
            })
    return edges


def build_ts(nodes: Dict[int, dict], edges: List[dict], sequence_length: float | None):
    positive_rights = [e["right"] for e in edges if e["right"] > e["left"]]
    points = sorted({e["left"] for e in edges if e["right"] <= e["left"]})
    if sequence_length is None:
        sequence_length = max(positive_rights) if positive_rights else (max(points) if points else 1.0)
    if sequence_length <= 0:
        sequence_length = 1.0

    point_intervals = {}
    if points:
        for i, p in enumerate(points):
            left = 0.0 if i == 0 else 0.5 * (points[i - 1] + p)
            right = sequence_length if i == len(points) - 1 else 0.5 * (p + points[i + 1])
            if right <= left:
                right = left + 1e-9
            point_intervals[p] = (left, right)

    for e in edges:
        if e["right"] <= e["left"]:
            e["left"], e["right"] = point_intervals[e["left"]]

    edge_map = {}
    for e in edges:
        if e["right"] <= e["left"]:
            continue
        key = (e["parent"], e["child"])
        edge_map.setdefault(key, []).append((e["left"], e["right"]))
    edges = []
    for (parent, child), intervals in edge_map.items():
        intervals = sorted(set(intervals))
        merged = []
        for left, right in intervals:
            if merged and left <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], right))
            else:
                merged.append((left, right))
        for left, right in merged:
            edges.append({"parent": parent, "child": child, "left": left, "right": right})

    # tskit requires parent times to be strictly older than child times.
    eps = 1e-9
    changed = True
    while changed:
        changed = False
        for e in edges:
            p = e["parent"]
            c = e["child"]
            if p in nodes and c in nodes and nodes[p]["time"] <= nodes[c]["time"]:
                nodes[p]["time"] = nodes[c]["time"] + eps
                changed = True

    tables = tskit.TableCollection(sequence_length=sequence_length)
    tables.nodes.metadata_schema = tskit.MetadataSchema.permissive_json()

    max_pop = max((n["pop"] for n in nodes.values()), default=-1)
    for _ in range(max_pop + 1):
        tables.populations.add_row()

    id_to_row = {}
    for nid in sorted(nodes):
        n = nodes[nid]
        flags = tskit.NODE_IS_SAMPLE if n["is_sample"] else 0
        id_to_row[nid] = tables.nodes.add_row(
            flags=flags,
            time=n["time"],
            population=n["pop"],
            metadata={"arg_id": nid, "inv": n["inv"]},
        )

    for e in edges:
        if e["parent"] not in id_to_row or e["child"] not in id_to_row:
            continue
        if e["right"] <= e["left"]:
            continue
        tables.edges.add_row(
            left=e["left"],
            right=e["right"],
            parent=id_to_row[e["parent"]],
            child=id_to_row[e["child"]],
        )

    tables.sort()
    return tables.tree_sequence()


def write_outputs(ts, out_base: str, width: int, height: int) -> None:
    trees_path = out_base + ".trees"
    svg_path = out_base + ".svg"
    png_path = out_base + ".png"

    ts.dump(trees_path)
    svg = ts.draw_svg(size=(width, height), y_axis=True, x_axis=True)
    svg = svg.replace("<svg ", "<svg style=\"background:#ffffff\" ", 1)
    with open(svg_path, "w", encoding="utf-8") as f:
        f.write(svg)
    cairosvg.svg2png(url=svg_path, write_to=png_path)

    print(f"Wrote {trees_path}")
    print(f"Wrote {svg_path}")
    print(f"Wrote {png_path}")


def run_id_from_path(path: str) -> str:
    m = re.search(r"arg_nodes_run(.+)\.tsv$", os.path.basename(path))
    return m.group(1) if m else "0"


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert laBP ARG TSV outputs to tskit .trees and images.")
    parser.add_argument("--input-dir", default="ARG_outputs")
    parser.add_argument("--output-dir", default="ARG_outputs")
    parser.add_argument("--sequence-length", type=float, default=None)
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--height", type=int, default=900)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    node_files = sorted(glob.glob(os.path.join(args.input_dir, "arg_nodes_run*.tsv")))
    if not node_files:
        print(f"No ARG node TSV files found in {args.input_dir}", file=sys.stderr)
        return 1

    try:
        for node_path in node_files:
            run_id = run_id_from_path(node_path)
            edge_path = os.path.join(args.input_dir, f"arg_edges_run{run_id}.tsv")
            if not os.path.exists(edge_path):
                raise FileNotFoundError(edge_path)
            nodes = read_nodes(node_path)
            edges = read_edges(edge_path)
            ts = build_ts(nodes, edges, args.sequence_length)
            write_outputs(ts, os.path.join(args.output_dir, f"arg_run{run_id}"), args.width, args.height)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
