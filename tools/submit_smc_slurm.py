#!/usr/bin/env python3
"""Submit one SMC parameter/demes pair as a Slurm array."""

import os
import shutil
import tempfile
from pathlib import Path

import yaml

if __package__:
    from .run_smc_batch import expand_grid_runs
else:
    from run_smc_batch import expand_grid_runs


ROOT = Path(__file__).resolve().parents[1]
DEFAULTS = {
    "Cores": 10,
    "MemoryPerTask": "8G",
    "TimeLimit": "00:20:00",
    "Partition": "",
    "Account": "",
    "OutputRoot": str(ROOT / "batch_outputs"),
    "Compile": True,
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
    settings["Compile"] = _boolean(settings["Compile"], "Compile")

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
