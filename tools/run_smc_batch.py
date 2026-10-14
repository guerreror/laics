#!/usr/bin/env python3
import argparse
import copy
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_DIR = ROOT / "SMC_input"
DEFAULT_OUTPUT_DIR = ROOT / "batch_outputs"
PARAM_RUNNER = ROOT / "src" / "param.py"
GENERATED_PATHS = [
    "smc_trees",
    "smc_hop_trace.csv",
    "smc_hop_events.csv",
    "Output_log.txt",
    "src/migration_matrices.json",
    "src/original_migration_matrices.json",
]
GRID_KEYS = ("GridSearch", "grid_search", "gridSearch")


def remove_existing_outputs():
    for rel in GENERATED_PATHS:
        path = ROOT / rel
        if path.is_dir():
            shutil.rmtree(path)
        elif path.exists():
            path.unlink()


def collect_outputs(result_dir: Path):
    result_dir.mkdir(parents=True, exist_ok=True)
    for rel in GENERATED_PATHS:
        path = ROOT / rel
        if not path.exists():
            continue
        dest = result_dir / Path(rel).name
        if dest.exists():
            if dest.is_dir():
                shutil.rmtree(dest)
            else:
                dest.unlink()
        if rel.startswith("src/"):
            shutil.copy2(path, dest)
        else:
            shutil.move(str(path), str(dest))


def load_yaml(path: Path):
    with open(path, "r") as f:
        return yaml.safe_load(f) or {}


def strip_grid(doc):
    doc = copy.deepcopy(doc)
    grid = {}
    for key in GRID_KEYS:
        if key in doc:
            grid = doc.pop(key) or {}
            break
    if not isinstance(grid, dict):
        raise ValueError("GridSearch must be a mapping of yaml_path: [values].")
    return doc, grid


def path_part(container, part):
    if isinstance(container, list):
        if part.isdigit():
            return int(part)
        for i, item in enumerate(container):
            if isinstance(item, dict) and item.get("name") == part:
                return i
    return part


def set_yaml_path(doc, dotted_path, value):
    parts = dotted_path.split(".")
    cur = doc
    for part in parts[:-1]:
        cur = cur[path_part(cur, part)]
    cur[path_part(cur, parts[-1])] = value


def safe_name(value):
    text = str(value)
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_")
    return text or "value"


def grid_variants(base_doc, grid, source_name):
    variants = []
    for dotted_path, values in grid.items():
        if not isinstance(values, list):
            values = [values]
        for value in values:
            doc = copy.deepcopy(base_doc)
            set_yaml_path(doc, str(dotted_path), value)
            name = f"{source_name}__{str(dotted_path).replace('.', '_')}__{safe_name(value)}"
            variants.append((name, doc, dotted_path, value))
    return variants


def write_yaml(path: Path, data):
    with open(path, "w") as f:
        yaml.safe_dump(data, f, sort_keys=False)


def doc_key(demes_doc, params_doc):
    return (
        yaml.safe_dump(demes_doc, sort_keys=True),
        yaml.safe_dump(params_doc, sort_keys=True),
    )


def suffix(path: Path, prefix: str):
    m = re.fullmatch(rf"{re.escape(prefix)}(?:_(.+))?\.ya?ml", path.name)
    return m.group(1) if m else None


def find_pairs(input_dir: Path):
    demes = {}
    params = {}
    for path in input_dir.glob("*.y*ml"):
        d = suffix(path, "demes")
        p = suffix(path, "parameters")
        if d is not None:
            demes[d or "default"] = path
        elif p is not None:
            params[p or "default"] = path

    pairs = []
    for key in sorted(demes, key=lambda x: (x != "default", x)):
        if key not in params:
            raise FileNotFoundError(f"Missing parameters file for {demes[key].name}; expected parameters_{key}.yaml")
        pairs.append((key, demes[key], params[key]))

    extra = sorted(set(params) - set(demes))
    if extra:
        names = ", ".join(params[k].name for k in extra)
        raise FileNotFoundError(f"Parameter file(s) without matching demes file: {names}")
    return pairs


def main():
    parser = argparse.ArgumentParser(description="Run SMC for every demes/parameters YAML pair in an input folder.")
    parser.add_argument("--input-dir", default=str(DEFAULT_INPUT_DIR), help="Folder containing demes_N.yaml and parameters_N.yaml pairs.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR), help="Folder where per-pair outputs are saved.")
    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()

    try:
        pairs = find_pairs(input_dir)
    except Exception as e:
        print(f"Batch input error: {e}", file=sys.stderr)
        return 1
    if not pairs:
        print(f"No demes/parameters YAML pairs found in {input_dir}", file=sys.stderr)
        return 1

    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Found {len(pairs)} demes/parameters pair(s) in {input_dir}", flush=True)

    with tempfile.TemporaryDirectory(prefix="laics_grid_") as tmp:
        tmp_dir = Path(tmp)
        for idx, (key, demes_file, params_file) in enumerate(pairs, 1):
            demes_base, demes_grid = strip_grid(load_yaml(demes_file))
            params_base, params_grid = strip_grid(load_yaml(params_file))
            runs = []
            if demes_grid:
                for name, doc, _, _ in grid_variants(demes_base, demes_grid, "demes"):
                    runs.append((name, doc, params_base))
            if params_grid:
                for name, doc, _, _ in grid_variants(params_base, params_grid, "parameters"):
                    runs.append((name, demes_base, doc))
            if not runs:
                runs.append(("base", demes_base, params_base))
            deduped_runs = []
            seen_runs = set()
            base_run_key = doc_key(demes_base, params_base)
            for run in runs:
                key_for_run = doc_key(run[1], run[2])
                if key_for_run in seen_runs:
                    continue
                seen_runs.add(key_for_run)
                if key_for_run == base_run_key:
                    run = ("base", run[1], run[2])
                deduped_runs.append(run)
            runs = deduped_runs

            for run_idx, (variant_name, demes_doc, params_doc) in enumerate(runs, 1):
                run_label = demes_file.stem if variant_name == "base" else f"{demes_file.stem}__{variant_name}"
                result_dir = output_dir / run_label
                tmp_demes = tmp_dir / f"{run_label}_demes.yaml"
                tmp_params = tmp_dir / f"{run_label}_parameters.yaml"
                write_yaml(tmp_demes, demes_doc)
                write_yaml(tmp_params, params_doc)

                print(
                    f"\n[{idx}/{len(pairs)}:{run_idx}/{len(runs)}] Running {demes_file.name} with {params_file.name}"
                    + ("" if variant_name == "base" else f" ({variant_name})"),
                    flush=True,
                )
                remove_existing_outputs()
                cmd = [sys.executable, str(PARAM_RUNNER), "--mode", "smc", "--config", str(tmp_params), "--demes", str(tmp_demes)]
                proc = subprocess.run(cmd, cwd=ROOT)
                if proc.returncode != 0:
                    print(
                        f"Run failed for {demes_file.name}/{params_file.name} ({variant_name}) with exit code {proc.returncode}",
                        file=sys.stderr,
                    )
                    return proc.returncode
                collect_outputs(result_dir)
                print(f"Saved outputs to {result_dir}", flush=True)

    print(f"\nBatch complete. Results are in {output_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
