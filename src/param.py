#!/usr/bin/env python3
import sys
import yaml
import demes
import subprocess
import argparse
import re
import random
from collections import defaultdict, deque

# ---------- File paths ----------
OTHER_YAML      = "src/other.yaml"
DEMES_YAML      = "src/demes.yaml"
EXECUTABLE      = "./executables/labp_v19"

# ---------- Defaults ----------
parameters = {
    "seed":         "1",
    "nruns":        "1",
    "kingman_coal": "1",
    "drift_sim":    "0",
    "msOutput":     "1",
    "popSizeVec":   "1000 1000",
    "inv_freq":     "0.2 0.3",              
    "speciation":   "0 0 0",
    "demography":   "0 0 0",
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
    "nCarriers":    "10 0",
}

# ---------- Helpers for demes ----------
def epoch_start_size(epoch):
    return int(epoch.start_size)

def leaf_demes(graph):
    # present-day demes: last epoch must end at 0
    return [d for d in graph.demes if d.epochs and d.epochs[-1].end_time == 0]

def pop_index_order(leaves):
    # sort pop-like names numerically: pop0, pop1, …; others after
    def keyfn(name):
        m = re.fullmatch(r"pop(\d+)", name)
        if m: return (0, int(m.group(1)))
        return (1, name)
    names = sorted([d.name for d in leaves], key=keyfn)
    return {name: i for i, name in enumerate(names)}

def build_children_map(graph):
    children = defaultdict(list)
    for d in graph.demes:
        if getattr(d, "ancestors", None):
            if len(d.ancestors) > 1:
                raise NotImplementedError(
                    f"Admixture is not supported for this parser (deme={d.name}, ancestors={d.ancestors})."
                )
            parent = d.ancestors[0]
            children[parent].append(d.name)
    for p in children:
        # deterministic: leaf-like labels last component numeric ordering
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
    Return the parent's epoch start_size at split time T.
    Prefer an epoch whose end_time == T; otherwise the first with end_time > T.
    """
    eps = 1e-9
    for ep in deme.epochs:
        if abs(ep.end_time - T) < eps:
            return int(ep.start_size)
    for ep in deme.epochs:
        if ep.end_time > T:
            return int(ep.start_size)
    # fallback
    return int(deme.epochs[-1].start_size)

def render_tree_ascii(graph, only_these_leaves=None):
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

# ---------- Speciation (sizes + events) ----------
def build_speciation_from_demes(graph, ancestor_freqs=None, default_F=0.2):
    """
    Build:
      - pop_sizes_str: "N0 N1 ... (current leaf sizes at end_time=0)"
      - speciation_str: "1 A B T F R A B T F R ..."
      - debug_info: dict with leaf order, events, and plan logs

    A,B are CURRENT group indices at time T (pre-merge state); R is the resulting
    size from the parent's epoch at T; F is the inversion frequency for the merged group.
    """
    leaves = leaf_demes(graph)
    if not leaves:
        raise ValueError("No present-day demes (end_time==0) found.")
    idx_map = pop_index_order(leaves)
    leaf_names_in_order = [name for name, _ in sorted(idx_map.items(), key=lambda kv: kv[1])]

    # Present-day sizes in that order (from the last epoch)
    sizes = []
    for name in leaf_names_in_order:
        d = next(dd for dd in leaves if dd.name == name)
        sizes.append(str(epoch_start_size(d.epochs[-1])))

    children = build_children_map(graph)
    by_name = {d.name: d for d in graph.demes}
    leafset = set(leaf_names_in_order)

    # Active groups evolve as we "rewind" backward splits to forward merges
    active_groups = [[name] for name in leaf_names_in_order]

    def find_group_idx(leaf):
        for i, g in enumerate(active_groups):
            if leaf in g:
                return i
        raise RuntimeError(f"Leaf {leaf} not found in current groups")

    # Internal nodes ordered by earliest child start_time (present-ward first)
    internal_nodes = list(children.keys())
    internal_nodes.sort(key=lambda n: min((by_name[c].start_time for c in children[n])), reverse=False)

    events = []       # list of (A,B,T,F,R)
    plan_log = []
    groups_log = []

    def get_parent_F(parent):
        if ancestor_freqs is None:
            return float(default_F)
        if parent not in ancestor_freqs:
            raise ValueError(f"ancestor_frequencies is missing a value for internal node '{parent}'")
        F = float(ancestor_freqs[parent])
        if not (0.0 <= F <= 1.0):
            raise ValueError(f"ancestor_frequencies['{parent}'] must be in [0,1], got {F}")
        return F

    for parent in internal_nodes:
        child_sets = []
        for child in children[parent]:
            under = leaves_under(child, children, leafset)
            if under:
                child_sets.append((child, sorted(under)))
        if len(child_sets) < 2:
            continue

        # Event time and resulting pop size (R from parent at this split)
        T = min(by_name[ch].start_time for ch, _ in child_sets)
        R = ancestor_size_at(by_name[parent], T)
        F_parent = get_parent_F(parent)

        # Map each child to its current group
        child_group_idxs = []
        for ch, under in child_sets:
            gi = find_group_idx(under[0])
            child_group_idxs.append((gi, ch, under))

        # Smallest index becomes sink (to mirror C++ index behavior)
        child_group_idxs.sort(key=lambda x: x[0])
        sink_idx = child_group_idxs[0][0]
        sink_child = child_group_idxs[0][1]

        # Fold the rest into the sink
        for gi, ch, under in child_group_idxs[1:]:
            A = sink_idx
            B = gi
            events.append((A, B, T, F_parent, R))
            plan_log.append(
                f"t={fmt_time(T)}: merge group {B} ({ch}:{under}) → group {A} ({sink_child}); F={F_parent:g}, R={R}"
            )
            # Update active groups & log snapshot
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

# ---------- Demography (built from all demes) ----------
def build_demography_from_demes_full(graph, leaf_order, spec_events):
    """
    Build demography from ALL demes (leaves + ancestors).

    For each deme D and each epoch with end_time > 0, at t=end_time emit a size
    change to the CURRENT group index that represents D at that time (if it exists),
    setting it to that epoch's start_size.

    Important: demography applies BEFORE speciation at the same time t.
    So we compute changes using the pre-merge mapping, then apply merges at t.
    """
    by_name = {d.name: d for d in graph.demes}
    npops   = len(leaf_order)

    children = build_children_map(graph)
    leafset  = set(leaf_order)

    def leaves_under_name(name):
        return sorted(leaves_under(name, children, leafset))

    # Collect all (t, deme_name, size) for end_time>0
    demog_records = defaultdict(list)  # t -> list[(deme_name, size)]
    for d in graph.demes:
        for ep in d.epochs:
            if ep.end_time and ep.end_time > 0:
                demog_records[ep.end_time].append((d.name, int(ep.start_size)))

    if not demog_records:
        return "0 0 0", []

    # Build merges at time lookup from spec events
    merges_at_time = defaultdict(list)  # t -> list[(A,B)]
    for (A, B, T, _, _) in spec_events:
        merges_at_time[T].append((A, B))

    # Track evolving groups
    active_groups = [[name] for name in leaf_order]  # list[list[str]]

    def group_index_for_deme(name):
        target = set(leaves_under_name(name))
        for i, g in enumerate(active_groups):
            if set(g) == target:
                return i
        return None  # doesn't exist yet in the current mapping

    demog_entries = ["1"]
    dbg_lines = []

    for t in sorted(demog_records.keys()):
        sizes = [0] * npops
        applied = []

        # (1) Apply demography at t using PRE-MERGE mapping
        for name, sz in demog_records[t]:
            gi = group_index_for_deme(name)
            if gi is not None and gi < len(active_groups):
                sizes[gi] = sz
                applied.append((gi, name, sz))

        if any(v != 0 for v in sizes[:len(active_groups)]):
            demog_entries.append(f"{t:g}")
            demog_entries += [str(v) for v in sizes]
            # Debug line: show only nonzero assignments
            applied_str = ", ".join([f"g{gi}({name})={sz}" for gi, name, sz in applied])
            dbg_lines.append(f"  t={fmt_time(t)}: {applied_str}")

        # (2) Apply all merges scheduled at time t (to match C++ ordering: demography before speciation on ties)
        for (A, B) in sorted(merges_at_time.get(t, []), key=lambda x: (x[0], x[1])):
            if A < len(active_groups) and B < len(active_groups):
                active_groups[A].extend(active_groups[B])
                del active_groups[B]

    return " ".join(demog_entries), dbg_lines

# ---------- Load other.yaml ----------
try:
    with open(OTHER_YAML, 'r') as f:
        other_params = yaml.safe_load(f) or {}
except Exception as e:
    print(f"Error loading {OTHER_YAML}: {e}", file=sys.stderr)
    sys.exit(1)

# Copy simple keys (the derived ones will be overwritten later)
for key in other_params:
    if key in parameters:
        parameters[key] = str(other_params[key])

# ancestor_frequencies (dict: internal node name -> F)
ancestor_freqs = other_params.get("ancestor_frequencies", None)
if ancestor_freqs is not None and not isinstance(ancestor_freqs, dict):
    print("Error: 'ancestor_frequencies' in other.yaml must be a mapping (dict).", file=sys.stderr)
    sys.exit(1)

# Optional default F when ancestor_frequencies not provided
DEFAULT_ANCESTOR_F = 0.2

# ---------- Load demes.yaml & derive sizes/speciation/migration/demography ----------
try:
    graph = demes.load(DEMES_YAML)
except Exception as e:
    print(f"Error loading {DEMES_YAML}: {e}", file=sys.stderr)
    sys.exit(1)

# Speciation & present-day pop sizes
try:
    pop_sizes_str, speciation_str, debug_info = build_speciation_from_demes(
        graph,
        ancestor_freqs=ancestor_freqs,
        default_F=DEFAULT_ANCESTOR_F,
    )
    parameters["popSizeVec"] = pop_sizes_str
    parameters["speciation"] = speciation_str
except (NotImplementedError, ValueError) as e:
    print(f"[speciation parser] {e}", file=sys.stderr)
    sys.exit(1)

# Demography from demes.yaml (all demes)
demog_str, demog_dbg = build_demography_from_demes_full(
    graph,
    debug_info["leaf_order"],
    debug_info["events"]
)
parameters["demography"] = demog_str

def _extract_default_mig_rate_from_yaml(path):
    try:
        with open(path, "r") as _f:
            raw = yaml.safe_load(_f) or {}
        dflt = (raw.get("defaults") or {}).get("migration") or {}
        rate = dflt.get("rate", None)
        return float(rate) if rate is not None else None
    except Exception:
        return None

mig_rate = None

try:
    if graph.migrations and hasattr(graph.migrations[0], "rate"):
        mig_rate = float(graph.migrations[0].rate)
except Exception:
    mig_rate = None

if mig_rate is None:
    mig_rate = _extract_default_mig_rate_from_yaml(DEMES_YAML)

if mig_rate is not None:
    parameters["migRate"] = f"{mig_rate:g}"


# ---------- Debug prints ----------
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

print("\n=== Demography plan (from demes.yaml) ===")
if parameters["demography"].startswith("1"):
    if demog_dbg:
        for line in demog_dbg:
            print(line)
    else:
        print("  (events exist, but all mapped to times where affected groups didn't yet exist pre-merge)")
else:
    print("  (none)")

# ---------- Allow user tweaks EXCEPT derived fields ----------
def print_parameters(params):
    for k, v in params.items():
        print(f"{k}: {v}")

def modify_params(parameters):
    print("Current Parameters:")
    print_parameters(parameters)
    ans = input("\nDo you want to modify any parameters? (Y/N): ").strip().lower()
    if ans == 'y':
        print("(Note: 'popSizeVec', 'speciation', and 'demography' are derived from demes.yaml and will be ignored if changed.)")
        user_input = input("\nEnter parameters to modify (e.g., --seed 2 --nruns 3): ")
        parser = argparse.ArgumentParser()
        for key in parameters:
            if key in ["popSizeVec", "speciation", "demography"]:
                continue
            if key in ["inv_freq","invRange","fixedSNPs","tempRead","nCarriers","snpPositions"]:
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
parameters["demography"] = demog_str

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
