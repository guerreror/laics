#!/usr/bin/env python3
"""Submit one SMC parameter/demes pair as a Slurm array."""

import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path

import yaml

if __package__:
    from .run_smc_batch import doc_key, grid_variants, strip_grid
else:
    from run_smc_batch import doc_key, grid_variants, strip_grid


ROOT = Path(__file__).resolve().parents[1]
DEFAULTS = {
    "Cores": 10,
    "MemoryPerTask": "8G",
    "TimeLimit": "00:20:00",
    "Partition": "",
    "Account": "",
    "OutputRoot": str(ROOT / "batch_outputs"),
}


# check if slurm mode is on
def _boolean(value, name):
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"1", "true"}:
        return True
    if text in {"0", "false"}:
        return False
    raise ValueError(f"{name} must be 0, 1, false, or true.")

# type/error check for yaml value, stupid
def _positive_integer(value, name):
    if isinstance(value, bool) or (isinstance(value, float) and not value.is_integer()):
        raise ValueError(f"{name} must be a positive integer.")
    try:
        value = int(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a positive integer.")
    if value <= 0:
        raise ValueError(f"{name} must be a positive integer.")
    return value

# yaml adv file path 
def _yaml_mapping(path, description):
    path = Path(path).resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{description} does not exist: {path}")
    with path.open() as handle:
        document = yaml.safe_load(handle) or {}
    if not isinstance(document, dict):
        raise ValueError(f"{description} must contain a YAML mapping: {path}")
    return document

# check default or adv settings for slurm submissions
def load_slurm_settings(parameters_path):
    """Load the Slurm toggle and optional advanced settings for one run."""
    parameters_path = Path(parameters_path).resolve()
    parameters = _yaml_mapping(parameters_path, "Parameter file")
    if not _boolean(parameters.get("SlurmMode", False), "SlurmMode"):
        return None

    settings = dict(DEFAULTS)
    advanced_name = str(parameters.get("SlurmConfig", "")).strip()
    if advanced_name:
        advanced_path = Path(os.path.expandvars(advanced_name)).expanduser()
        if not advanced_path.is_absolute():
            advanced_path = parameters_path.parent / advanced_path
        advanced = _yaml_mapping(advanced_path, "SlurmConfig")
        unknown = sorted(set(advanced) - set(DEFAULTS))
        if unknown:
            raise ValueError("Unknown Slurm setting(s): " + ", ".join(unknown))
        settings.update(advanced)

    settings["Cores"] = _positive_integer(settings["Cores"], "Cores")

    for name in ("MemoryPerTask", "TimeLimit"):
        settings[name] = str(settings[name]).strip()
        if not settings[name]:
            raise ValueError(f"{name} must not be empty.")
    for name in ("Partition", "Account"):
        settings[name] = str(settings[name] or "").strip()

    output_root = Path(os.path.expandvars(str(settings["OutputRoot"]))).expanduser()
    if not output_root.is_absolute():
        output_root = ROOT / output_root
    settings["OutputRoot"] = str(output_root.resolve())
    return settings


# auto resource/replicate allocation for slurm array tasks
def plan_work(number_of_replicates, settings):
    """Divide replicates among one-CPU Slurm array tasks."""
    number_of_replicates = _positive_integer(
        number_of_replicates, "NumberOfReplicates"
    )
    task_count = min(settings["Cores"], number_of_replicates)
    replicates_per_task, extra_replicates = divmod(
        number_of_replicates, task_count
    )
    tasks = [
        replicates_per_task + 1
        if task_index < extra_replicates
        else replicates_per_task
        for task_index in range(task_count)
    ]
    return {
        "NumberOfReplicates": number_of_replicates,
        "TaskCount": task_count,
        "Tasks": tasks,
    }

# create a directory for each slurm task, and save the param/demes files for those tasks
def materialize_task_inputs(parameters_path, demes_path, run_directory, plan):
    """Write one parameter file per task and one shared demes snapshot."""
    parameters_path = Path(parameters_path).resolve()
    demes_path = Path(demes_path).resolve()
    parameters = _yaml_mapping(parameters_path, "Parameter file")
    if not demes_path.is_file():
        raise FileNotFoundError(f"Demes file does not exist: {demes_path}")

    try:
        base_seed = int(parameters.get("Seed", 1))
    except (TypeError, ValueError):
        raise ValueError("Seed must be an integer.")

    run_directory = Path(run_directory).resolve()
    run_directory.mkdir(parents=True, exist_ok=True)
    shared_demes = run_directory / "demes.yaml"
    shutil.copy2(demes_path, shared_demes)

    task_inputs = []
    for task_index, replicate_count in enumerate(plan["Tasks"]):
        task_directory = run_directory / f"task_{task_index:04d}"
        task_directory.mkdir(parents=True, exist_ok=True)

        task_parameters = dict(parameters)
        task_parameters.pop("SlurmMode", None)
        task_parameters.pop("SlurmConfig", None)
        task_parameters["NumberOfReplicates"] = str(replicate_count)
        task_parameters["Seed"] = str(base_seed + task_index)

        task_parameter_path = task_directory / "parameters.yaml"
        with task_parameter_path.open("w") as handle:
            yaml.safe_dump(task_parameters, handle, sort_keys=False)

        task_inputs.append(
            {
                "TaskIndex": task_index,
                "Replicates": replicate_count,
                "Seed": base_seed + task_index,
                "Directory": str(task_directory),
                "Parameters": str(task_parameter_path),
            }
        )

    return {"Demes": str(shared_demes), "Tasks": task_inputs}

# allow the functionallity of grid search mode with slurm mode
def materialize_scenarios(parameters_path, demes_path, run_directory, settings):
    """Expand existing GridSearch scenarios and prepare their array tasks."""
    parameters_path = Path(parameters_path).resolve()
    demes_path = Path(demes_path).resolve()
    parameters = _yaml_mapping(parameters_path, "Parameter file")
    demes = _yaml_mapping(demes_path, "Demes file")
    runs = expand_grid_runs(demes, parameters)

    run_directory = Path(run_directory).resolve()
    scenarios = []
    with tempfile.TemporaryDirectory(prefix="laics_slurm_grid_") as temporary:
        temporary = Path(temporary)
        for variant_name, demes_document, parameters_document in runs:
            run_label = (
                demes_path.stem
                if variant_name == "base"
                else f"{demes_path.stem}__{variant_name}"
            )
            scenario_directory = run_directory / run_label
            temporary_demes = temporary / f"{run_label}_demes.yaml"
            temporary_parameters = temporary / f"{run_label}_parameters.yaml"
            with temporary_demes.open("w") as handle:
                yaml.safe_dump(demes_document, handle, sort_keys=False)
            with temporary_parameters.open("w") as handle:
                yaml.safe_dump(parameters_document, handle, sort_keys=False)

            replicate_count = parameters_document.get("NumberOfReplicates", 1)
            plan = plan_work(replicate_count, settings)
            materialized = materialize_task_inputs(
                temporary_parameters,
                temporary_demes,
                scenario_directory,
                plan,
            )
            scenarios.append(
                {
                    "Name": run_label,
                    "Variant": variant_name,
                    "Directory": str(scenario_directory),
                    "Plan": plan,
                    "Inputs": materialized,
                }
            )
    return scenarios

    
def expand_grid_runs(demes_doc, params_doc):
    """Adapt the existing grid helpers without changing the local batch runner."""
    demes_base, demes_grid = strip_grid(demes_doc)
    params_base, params_grid = strip_grid(params_doc)
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
    return deduped_runs



# before submitting slurm array, make sure compiled execs are availible
def validate_runtime(root=ROOT, sbatch_command="sbatch"):
    """Verify the installed LAiCS runtime and Slurm client are available."""
    root = Path(root).resolve()
    executable = root / "executables" / "labp_smc"
    parameter_runner = root / "src" / "param.py"

    if not executable.is_file():
        raise FileNotFoundError(f"SMC executable does not exist: {executable}")
    if not os.access(executable, os.X_OK):
        raise PermissionError(f"SMC executable is not executable: {executable}")
    if not parameter_runner.is_file():
        raise FileNotFoundError(f"Parameter runner does not exist: {parameter_runner}")

    sbatch = shutil.which(sbatch_command)
    if sbatch is None:
        raise FileNotFoundError(f"Slurm command is not available: {sbatch_command}")

    return {
        "Root": str(root),
        "Executable": str(executable),
        "ParameterRunner": str(parameter_runner),
        "Python": sys.executable,
        "Sbatch": sbatch,
    }

# Write the Slurm array script for all tasks in all scenarios
def write_array_script(scenarios, settings, run_directory, runtime):
    """Write the one-CPU Slurm array script for all prepared scenarios."""
    task_directories = [
        task["Directory"]
        for scenario in scenarios
        for task in scenario["Inputs"]["Tasks"]
    ]
    if not task_directories:
        raise ValueError("Cannot write a Slurm array script without tasks.")

    run_directory = Path(run_directory).resolve()
    run_directory.mkdir(parents=True, exist_ok=True)
    if any(character.isspace() for character in str(run_directory)):
        raise ValueError("The Slurm run directory must not contain whitespace.")

    concurrency = min(settings["Cores"], len(task_directories))
    directives = [
        "#SBATCH --job-name=laics-smc",
        "#SBATCH --cpus-per-task=1",
        f"#SBATCH --mem={settings['MemoryPerTask']}",
        f"#SBATCH --time={settings['TimeLimit']}",
        f"#SBATCH --array=0-{len(task_directories) - 1}%{concurrency}",
        f"#SBATCH --output={run_directory}/slurm-%A_%a.out",
    ]
    if settings["Partition"]:
        directives.append(f"#SBATCH --partition={settings['Partition']}")
    if settings["Account"]:
        directives.append(f"#SBATCH --account={settings['Account']}")

    quoted_tasks = "\n".join(
        f"    {shlex.quote(str(Path(task).resolve()))}"
        for task in task_directories
    )
    script = f"""#!/bin/bash
{chr(10).join(directives)}

set -euo pipefail

ROOT={shlex.quote(runtime['Root'])}
PYTHON={shlex.quote(runtime['Python'])}
TASK_DIRS=(
{quoted_tasks}
)

task_id="${{SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set}}"
if (( task_id < 0 || task_id >= ${{#TASK_DIRS[@]}} )); then
    echo "Invalid Slurm array task ID: $task_id" >&2
    exit 2
fi

TASK_DIR="${{TASK_DIRS[$task_id]}}"
SCENARIO_DIR="$(dirname "$TASK_DIR")"

cd "$ROOT"
"$PYTHON" "$ROOT/src/param.py" \\
    --mode smc \\
    --config "$TASK_DIR/parameters.yaml" \\
    --demes "$SCENARIO_DIR/demes.yaml" \\
    --output-dir "$TASK_DIR" \\
    > "$TASK_DIR/simulation.log" 2>&1
"""

    script_path = run_directory / "run_array.slurm"
    script_path.write_text(script)
    script_path.chmod(0o755)
    return script_path


# output directory needed for each slurm sub
def create_run_directory(parameters_path, output_root):
    """Create a unique directory for one Slurm-mode invocation."""
    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", Path(parameters_path).stem).strip("_")
    label = label or "smc"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    run_directory = Path(output_root).resolve() / f"{label}_{timestamp}"
    run_directory.mkdir(parents=True, exist_ok=False)
    return run_directory

# submit indiv slurm arrays
def submit_array(script_path, runtime):
    """Submit one prepared Slurm array and return its job ID."""
    completed = subprocess.run(
        [runtime["Sbatch"], "--parsable", str(Path(script_path).resolve())],
        cwd=runtime["Root"],
        text=True,
        capture_output=True,
    )
    if completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError(f"Slurm array submission failed: {message}")

    submission = completed.stdout.strip()
    if not submission:
        raise RuntimeError("Slurm array submission returned no job ID.")
    return submission.split(";", 1)[0]


def prepare_and_submit(parameters_path, demes_path, dry_run=False):
    """Prepare a complete Slurm run and optionally submit its array."""
    settings = load_slurm_settings(parameters_path)
    if settings is None:
        raise ValueError("SlurmMode is not enabled in the parameter file.")

    runtime = validate_runtime()
    run_directory = create_run_directory(parameters_path, settings["OutputRoot"])
    scenarios = materialize_scenarios(
        parameters_path,
        demes_path,
        run_directory,
        settings,
    )
    script_path = write_array_script(
        scenarios,
        settings,
        run_directory,
        runtime,
    )
    task_count = sum(
        len(scenario["Inputs"]["Tasks"])
        for scenario in scenarios
    )

    job_id = None if dry_run else submit_array(script_path, runtime)
    return {
        "RunDirectory": str(run_directory),
        "Script": str(script_path),
        "ScenarioCount": len(scenarios),
        "TaskCount": task_count,
        "JobID": job_id,
    }


# run everythin 
# there's a dry run option if you want
def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Prepare and submit one SMC run as a Slurm array."
    )
    parser.add_argument("--config", required=True, help="SMC parameter YAML.")
    parser.add_argument("--demes", required=True, help="Demes YAML.")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Prepare all files but do not call sbatch.",
    )
    args = parser.parse_args(argv)

    try:
        result = prepare_and_submit(args.config, args.demes, args.dry_run)
    except Exception as error:
        print(f"Slurm preparation failed: {error}", file=sys.stderr)
        return 1

    print(f"Prepared {result['ScenarioCount']} scenario(s) and {result['TaskCount']} task(s).")
    print(f"Run directory: {result['RunDirectory']}")
    print(f"Array script: {result['Script']}")
    if result["JobID"] is None:
        print("Dry run complete; no Slurm job was submitted.")
    else:
        print(f"Submitted Slurm array job {result['JobID']}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
