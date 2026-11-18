#!/usr/bin/env python3
"""Submit one SMC parameter/demes pair as a Slurm array."""

import os
from pathlib import Path

import yaml


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


def _boolean(value, name):
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"1", "true"}:
        return True
    if text in {"0", "false"}:
        return False
    raise ValueError(f"{name} must be 0, 1, false, or true.")


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


def _yaml_mapping(path, description):
    path = Path(path).resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{description} does not exist: {path}")
    with path.open() as handle:
        document = yaml.safe_load(handle) or {}
    if not isinstance(document, dict):
        raise ValueError(f"{description} must contain a YAML mapping: {path}")
    return document


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
