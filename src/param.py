#!/usr/bin/env python3
import sys
import yaml
import demes
import subprocess
import argparse
import json
import re
import random
from collections import defaultdict, deque

# ---------- File paths ----------
OTHER_YAML      = "src/other.yaml"
DEMES_YAML      = "src/demes.yaml"
DEMOGRAPHY_JSON = "src/demo_dict.json"   
EXECUTABLE      = "./executables/labp_v19"

# ---------- Defaults ----------
parameters = {
    "seed":         "1",
    "nruns":        "1",
    "kingman_coal": "1",
    "drift_sim":    "0",
    "msOutput":     "1",
    "popSizeVec":   "10000 10000 10000 10000",
    "inv_freq":     "0.2 0.3 0.4",
    "speciation":   "1 0 2 1000 0.3",
    "demography":   "0 0 0",             # placeholder; may be overridden by JSON
    "inv_age":      "0",
    "migRate":      "0.02",
    "BasesPerMorgan":"1e8",
    "randPhi":      "0",
    "phi":          "0.2",
    "invRange":     "0 1e3",
    "fixedSNPs":    "1 10 1e-6",
    "n_SNPs":       "0",
    "snpPositions": "500 550 600 650 700 750 800 850 900 950",
    "randomSample": "0",
    "tempRead":     "10 0",
    "nCarriers":    "10 0"
}

# ---------- Helpers for demes ----------
def epoch_start_size(epoch):
    return int(epoch.start_size)

def leaf_demes(graph):
    return [d for d in graph.demes if d.epochs[-1].end_time == 0]

def pop_index_order(leaves):
    def keyfn(name):
        m = re.fullmatch(r"pop(\d+)", name)
        if m:
            return (0, int(m.group(1)))
        return (1, name)
    names = sorted([d.name for d in leaves], key=keyfn)
    return {name: i for i, name in enumerate(names)}

def build_children_map(graph):
    children = defaultdict(list)
    for d in graph.demes:
        if getattr(d, "ancestors", None):
            if len(d.ancestors) > 1:
                raise NotImplementedError(
                    f"Admixture not supported for speciation parsing (deme={d.name}, ancestors={d.ancestors})"
                )
            parent = d.ancestors[0]
            children[parent].append(d.name)
    # deterministic ordering
    for p in children:
        children[p].sort(key=lambda x: (x.startswith("pop"), x))
    return children

def leaves_under(node, children, leafset):
    out = []
    q = deque([node])
    while q:
        u = q.popleft()
        if u in leafset:
            out.append(u)
        for v in children.get(u, []):
            q.append(v)
    return out

def find_roots(graph):
    return [d.name for d in graph.demes if len(getattr(d, "ancestors", [])) == 0]

def pick_root(graph):
    roots = find_roots(graph)
    if not roots:
        raise ValueError("No root deme (a deme with no ancestors) found in demes.yaml")
    if len(roots) == 1:
        return roots[0]
    by_name = {d.name: d for d in graph.demes}
    roots.sort(key=lambda n: by_name[n].start_time, reverse=True)
    return roots[0]

def fmt_time(x):
    try:
        if x == float("inf"):
            return "inf"
        return f"{x:g}"
    except Exception:
        return str(x)

def ancestor_size_at(deme, T):
    """
    Return the parent's epoch start_size at the split time T.
    Prefer an epoch whose end_time == T; otherwise choose the first
    epoch with end_time > T.
    """
    eps = 1e-9
    for ep in deme.epochs:
        if abs(ep.end_time - T) < eps:
            return int(ep.start_size)
    for ep in deme.epochs:
        if ep.end_time > T:
            return int(ep.start_size)
    # Fallback to the last epoch's start_size (shouldn't happen in well-formed graphs)
    return int(deme.epochs[-1].start_size)

def render_tree_ascii(graph, only_these_leaves=None):
    """
    Pretty-print tree (past -> present) with connectors.
    """
    children = build_children_map(graph)
    root = pick_root(graph)
    by_name = {d.name: d for d in graph.demes}

    if only_these_leaves:
        needed = set()
        parents = {c: (by_name[c].ancestors[0] if getattr(by_name[c], "ancestors", None) else None)
                   for c in by_name}
        for leaf in only_these_leaves:
            u = leaf
            while u is not None and u not in needed:
                needed.add(u)
                u = parents[u]
        pruned = {p: [c for c in cs if c in needed] for p, cs in children.items() if p in needed}
        children = pruned

    lines = []
    def dfs(u, prefix="", depth=0, is_last=True):
        connector = "" if depth == 0 else ("└─ " if is_last else "├─ ")
        label = f"{u}  (start={fmt_time(by_name[u].start_time)})"
        lines.append(prefix + connector + label)
        kids = children.get(u, [])
        for i, v in enumerate(kids):
            last = (i == len(kids) - 1)
            next_prefix = prefix if depth == 0 else prefix + ("   " if is_last else "│  ")
            dfs(v, next_prefix, depth + 1, last)
    dfs(root)
    return "\n".join(lines)

def build_speciation_from_demes(graph, ancestor_freqs=None):
    """
    Build:
      - pop_sizes_str: "N0 N1 N2 ..."
      - speciation_str: "1 A B T F R A B T F R ..."
      - debug_info: dict
    If ancestor_freqs is provided (dict: parent_name -> F), use it for F.
    Otherwise fall back to random choices in {0.1,0.2,0.3,0.4,0.5}.
    """
    leaves = leaf_demes(graph)
    if not leaves:
        raise ValueError("No present-day demes (end_time==0) found.")
    idx_map = pop_index_order(leaves)
    leaf_names_in_order = [name for name, _ in sorted(idx_map.items(), key=lambda kv: kv[1])]

    # Present-day sizes in our index order
    sizes = []
    for name in leaf_names_in_order:
        d = next(dd for dd in leaves if dd.name == name)
        sizes.append(str(epoch_start_size(d.epochs[-1])))

    children = build_children_map(graph)
    by_name = {d.name: d for d in graph.demes}

    # active_groups: list of lists of present-day leaf names
    active_groups = [[name] for name in leaf_names_in_order]

    def find_group_idx(leaf):
        for i, g in enumerate(active_groups):
            if leaf in g:
                return i
        raise RuntimeError(f"Leaf {leaf} not found in current groups")

    internal_nodes = list(children.keys())
    # Earliest (present-ward) splits first
    internal_nodes.sort(key=lambda n: min((by_name[c].start_time for c in children[n])), reverse=False)

    default_freqs = [0.1, 0.2, 0.3, 0.4, 0.5]
    events = []       # (A, B, T, F, R)
    plan_log = []
    groups_log = []

    # Validation helper
    def get_parent_F(parent):
        if ancestor_freqs is None:
            return random.choice(default_freqs)
        if parent not in ancestor_freqs:
            raise ValueError(f"ancestor_frequencies is missing a value for internal node '{parent}'")
        F = float(ancestor_freqs[parent])
        if not (0.0 <= F <= 1.0):
            raise ValueError(f"ancestor_frequencies['{parent}'] must be in [0,1], got {F}")
        return F

    for parent in internal_nodes:
        child_sets = []
        for child in children[parent]:
            under = leaves_under(child, children, set(leaf_names_in_order))
            if under:
                child_sets.append((child, sorted(under)))
        if len(child_sets) < 2:
            continue

        # Event time and resulting size (R comes from the parent’s epoch at that time)
        T = min(by_name[ch].start_time for ch, _ in child_sets)
        R = ancestor_size_at(by_name[parent], T)
        F_parent = get_parent_F(parent)

        # Map each child to current group index
        child_group_idxs = []
        for ch, under in child_sets:
            gi = find_group_idx(under[0])
            child_group_idxs.append((gi, ch, under))

        # Use the smallest current index as sink to match our C++ index behavior
        child_group_idxs.sort(key=lambda x: x[0])
        sink_idx = child_group_idxs[0][0]
        sink_child = child_group_idxs[0][1]

        # Pairwise fold-in of the remaining children into the sink
        for gi, ch, under in child_group_idxs[1:]:
            A = sink_idx
            B = gi
            F = F_parent
            events.append((A, B, T, F, R))
            plan_log.append(
                f"t={fmt_time(T)}: merge group {B} ({ch}:{under}) → group {A} ({sink_child}); F={F:g}, R={R}"
            )

            # Apply merge to groups and log the new mapping
            active_groups[A].extend(active_groups[B])
            del active_groups[B]
            groups_log.append([list(g) for g in active_groups])

    pop_sizes_str = " ".join(sizes)
    if events:
        speciation_flat = ["1"]
        for (A, B, T, F, R) in events:
            speciation_flat += [str(A), str(B), f"{T:g}", f"{F:g}", str(R)]
        speciation_str = " ".join(speciation_flat)
    else:
        speciation_str = "0 0 0"

    debug_info = {
        "leaf_order": leaf_names_in_order,
        "initial_groups": [[x] for x in leaf_names_in_order],
        "events": events,
        "plan_log": plan_log,
        "groups_log": groups_log,
    }
    return pop_sizes_str, speciation_str, debug_info

# ---------- Load other.yaml ----------
try:
    with open(OTHER_YAML, 'r') as f:
        other_params = yaml.safe_load(f) or {}
except Exception as e:
    print(f"Error loading {OTHER_YAML}: {e}", file=sys.stderr)
    sys.exit(1)

# Copy simple keys; leave popSizeVec/speciation to be derived here
for key in other_params:
    if key in parameters:
        parameters[key] = str(other_params[key])

# Pull ancestor_frequencies (dict) if present
ancestor_freqs = other_params.get("ancestor_frequencies", None)
if ancestor_freqs is not None and not isinstance(ancestor_freqs, dict):
    print("Error: 'ancestor_frequencies' in other.yaml must be a mapping (dict).", file=sys.stderr)
    sys.exit(1)

# ---------- Load demes.yaml & derive sizes/speciation/migration ----------
try:
    graph = demes.load(DEMES_YAML)
except Exception as e:
    print(f"Error loading {DEMES_YAML}: {e}", file=sys.stderr)
    sys.exit(1)

try:
    pop_sizes_str, speciation_str, debug_info = build_speciation_from_demes(graph, ancestor_freqs=ancestor_freqs)
    parameters["popSizeVec"] = pop_sizes_str
    parameters["speciation"] = speciation_str
except NotImplementedError as e:
    print(f"[speciation parser] {e}", file=sys.stderr)
    sys.exit(1)
except ValueError as e:
    print(f"[speciation parser] {e}", file=sys.stderr)
    sys.exit(1)

# Migration (first rate if present)
if graph.migrations:
    try:
        parameters["migRate"] = str(graph.migrations[0].rate)
    except Exception:
        pass

# ---------- Optional: sparse demography JSON ----------
try:
    with open(DEMOGRAPHY_JSON, 'r') as f:
        demog_map = json.load(f)
except FileNotFoundError:
    demog_map = {}

if demog_map:
    demog_entries = ["1"]
    for t_str, sizes in sorted(demog_map.items(), key=lambda kv: float(kv[0])):
        demog_entries.append(t_str)
        demog_entries += [str(int(x)) for x in sizes]
    parameters["demography"] = " ".join(demog_entries)

# ---------- ASCII tree + plan print ----------
print("\n=== Parsed speciation tree (past → present) ===")
ascii_tree = render_tree_ascii(graph, only_these_leaves=debug_info["leaf_order"])
print(ascii_tree)

print("\n=== Present-day leaves (order → indices) ===")
for i, name in enumerate(debug_info["leaf_order"]):
    print(f"  {i}: {name}")

print("\n=== Speciation plan (after reindexing between events) ===")
if debug_info["events"]:
    for line in debug_info["plan_log"]:
        print("  " + line)
    if debug_info["groups_log"]:
        print("\n--- Group reindex snapshots ---")
        for step, groups in enumerate(debug_info["groups_log"], 1):
            print(f"  after merge {step}: {groups}")
else:
    print("  (none)")

# ---------- Allow user tweaks EXCEPT popSizeVec/speciation (derived) ----------
def print_parameters(params):
    for k, v in params.items():
        print(f"{k}: {v}")

def modify_params(parameters):
    print("Current Parameters:")
    print_parameters(parameters)
    ans = input("\nDo you want to modify any parameters? (Y/N): ").strip().lower()
    if ans == 'y':
        print("(Note: 'popSizeVec' and 'speciation' are derived from demes.yaml and will be ignored if changed.)")
        user_input = input("\nEnter parameters to modify (e.g., --seed 2 --nruns 3): ")
        parser = argparse.ArgumentParser()
        for key in parameters:
            if key in ["popSizeVec", "speciation"]:
                continue
            if key in ["inv_freq","demography","invRange","fixedSNPs","tempRead","nCarriers","snpPositions"]:
                parser.add_argument(f"--{key}", nargs='+')
            else:
                parser.add_argument(f"--{key}")
        if user_input.strip():
            args = parser.parse_args(user_input.split())
            for key, val in vars(args).items():
                if val is not None:
                    parameters[key] = " ".join(val) if isinstance(val, list) else str(val)
        print("\nUpdated Parameters:")
        print_parameters(parameters)
    else:
        print("\nNo changes made.")

modify_params(parameters)

# Re-assert derived values (in case user tried to change them)
parameters["popSizeVec"] = pop_sizes_str
parameters["speciation"] = speciation_str

# ---------- Final assembly ----------
keys_order = [
    "seed","nruns","kingman_coal","drift_sim","msOutput",
    "popSizeVec","inv_freq","speciation","demography","inv_age",
    "migRate","BasesPerMorgan","randPhi","phi","invRange","fixedSNPs",
    "n_SNPs","snpPositions","randomSample","tempRead","nCarriers"
]
args_list = [parameters[k] for k in keys_order]

print("\nFinal parameters being passed:")
for k in keys_order:
    print(f"{k}: {parameters[k]}")

# ---------- Run ----------
try:
    result = subprocess.run([EXECUTABLE] + args_list, capture_output=True, text=True)
except FileNotFoundError:
    print(f"\nERROR: Executable not found at {EXECUTABLE}", file=sys.stderr)
    sys.exit(1)

print("\nC++ Program Output:\n", result.stdout)
if result.stderr:
    print("\nErrors:\n", result.stderr)
    with open("Output_log.txt", "w") as log_file:
        log_file.write(result.stderr)