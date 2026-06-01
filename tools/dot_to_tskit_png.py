#!/usr/bin/env python3
import argparse
import os
import re
import sys
from typing import Dict, List, Tuple


NODE_RE = re.compile(r'^\s*node(\d+)\s+\[label="([^"]+)"\];\s*$')
EDGE_RE = re.compile(r'^\s*node(\d+)\s*->\s*node(\d+)')
INV_RE = re.compile(r"inv=(\d+)")
POP_RE = re.compile(r"pop=(\d+)")
TIME_RE = re.compile(r"t=([0-9eE+\-\.]+)")


def parse_dot(dot_path: str):
    nodes: Dict[int, dict] = {}
    edges: List[Tuple[int, int]] = []

    with open(dot_path, "r", encoding="utf-8") as f:
        for line in f:
            n = NODE_RE.match(line)
            if n:
                node_id = int(n.group(1))
                label = n.group(2)
                inv_m = INV_RE.search(label)
                pop_m = POP_RE.search(label)
                time_m = TIME_RE.search(label)
                if inv_m is None or pop_m is None or time_m is None:
                    raise ValueError(f"Could not parse inv/pop/time for node {node_id}")
                nodes[node_id] = {
                    "id": node_id,
                    "inv": int(inv_m.group(1)),
                    "pop": int(pop_m.group(1)),
                    "time": float(time_m.group(1)),
                    "label": label,
                }
                continue

            e = EDGE_RE.match(line)
            if e:
                parent = int(e.group(1))
                child = int(e.group(2))
                edges.append((parent, child))

    if not nodes:
        raise ValueError("No nodes found in DOT file.")
    if not edges:
        raise ValueError("No edges found in DOT file.")

    return nodes, edges


def build_tskit(nodes: Dict[int, dict], edges: List[Tuple[int, int]], sequence_length: float):
    try:
        import tskit  # type: ignore
    except Exception as exc:  # pragma: no cover
        raise RuntimeError(
            "tskit is required. Install with: pip install tskit"
        ) from exc

    children = set(c for _, c in edges)
    parents = set(p for p, _ in edges)
    leaves = children - parents

    max_pop = max(n["pop"] for n in nodes.values())
    tables = tskit.TableCollection(sequence_length=sequence_length)
    tables.nodes.metadata_schema = tskit.MetadataSchema.permissive_json()

    for _ in range(max_pop + 1):
        tables.populations.add_row()

    # Ensure strict parent>child times required by tskit.
    # Some DOT exports can contain equal times on connected nodes.
    eps = 1e-9
    changed = True
    while changed:
        changed = False
        for parent, child in edges:
            pt = nodes[parent]["time"]
            ct = nodes[child]["time"]
            if pt <= ct:
                nodes[parent]["time"] = ct + eps
                changed = True

    # tskit row index per original node id
    id_to_row: Dict[int, int] = {}
    for nid, nd in nodes.items():
        flags = tskit.NODE_IS_SAMPLE if nid in leaves else 0
        row = tables.nodes.add_row(
            flags=flags,
            time=nd["time"],
            population=nd["pop"],
            metadata={"orig_id": nd["id"], "inv": nd["inv"], "label": nd["label"]},
        )
        id_to_row[nid] = row

    for parent, child in edges:
        if parent not in id_to_row or child not in id_to_row:
            continue
        tables.edges.add_row(
            left=0.0,
            right=sequence_length,
            parent=id_to_row[parent],
            child=id_to_row[child],
        )

    tables.sort()
    return tables.tree_sequence(), id_to_row


def render_svg(ts, out_svg: str, width: int, height: int, show_labels: bool):
    node_labels = None
    if show_labels:
        node_labels = {
            u: str(ts.node(u).metadata.get("orig_id", u))
            for u in range(ts.num_nodes)
        }
    tree = ts.first()
    svg = tree.draw_svg(
        size=(width, height),
        node_labels=node_labels,
        y_axis=True,
        x_axis=True,
    )
    # Make output readable on any viewer/theme:
    # 1) remove checkered tskit background group
    # 2) force white canvas
    # 3) use dark-gray edges and larger fonts
    svg = re.sub(r"<g class=\"background\">.*?</g>", "", svg, flags=re.DOTALL)
    svg = svg.replace(
        "<svg ",
        "<svg style=\"background:#ffffff\" ",
        1,
    )
    svg = svg.replace(".axes line, .edge {stroke: black; fill: none}",
                      ".axes line, .edge {stroke: #444; fill: none}")
    svg = svg.replace(".axes, .tree {font-size: 14px; text-anchor: middle}",
                      ".axes, .tree {font-size: 16px; text-anchor: middle}")

    with open(out_svg, "w", encoding="utf-8") as f:
        f.write(svg)


def svg_to_png(svg_path: str, png_path: str):
    try:
        import cairosvg  # type: ignore
    except Exception as exc:  # pragma: no cover
        raise RuntimeError(
            "PNG export requires cairosvg. Install with: pip install cairosvg"
        ) from exc
    cairosvg.svg2png(url=svg_path, write_to=png_path)


def main() -> int:
    p = argparse.ArgumentParser(description="Convert DOT tree to tskit-rendered PNG/SVG.")
    p.add_argument("dot_file", help="Input DOT file")
    p.add_argument("--sequence-length", type=float, default=1.0, help="Tree sequence length")
    p.add_argument("--svg", help="Output SVG path")
    p.add_argument("--png", help="Output PNG path")
    p.add_argument("--width", type=int, default=2400, help="SVG/PNG width in pixels")
    p.add_argument("--height", type=int, default=1600, help="SVG/PNG height in pixels")
    p.add_argument(
        "--show-labels",
        action="store_true",
        help="Show original node-id labels (off by default for readability)",
    )
    args = p.parse_args()

    if not os.path.isfile(args.dot_file):
        print(f"Error: file not found: {args.dot_file}", file=sys.stderr)
        return 1

    base, _ = os.path.splitext(args.dot_file)
    out_svg = args.svg or (base + ".tskit.svg")
    out_png = args.png or (base + ".tskit.png")

    try:
        nodes, edges = parse_dot(args.dot_file)
        ts, _ = build_tskit(nodes, edges, args.sequence_length)
        render_svg(ts, out_svg, args.width, args.height, args.show_labels)
        svg_to_png(out_svg, out_png)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote {out_svg}")
    print(f"Wrote {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
