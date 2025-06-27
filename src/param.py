#!/usr/bin/env python3
import sys
import yaml
import demes
import subprocess
import argparse

# Default file names for the YAML files.
OTHER_YAML = "src/other.yaml"
DEMES_YAML = "src/demes.yaml"

# The default parameters dictionary (its structure remains unchanged)
parameters = {
    "seed": "1",
    "nruns": "1",
    "kingman_coal": "1",
    "drift_sim": "0",
    "msOutput": "1",
    "popSizeVec": "10000 10000",
    "inv_freq": "0.2 0.4",
    "speciation": "1 10000 0.2",
    "demography": "0 0 0",
    "inv_age": "0",
    "migRate": "0.02",
    "BasesPerMorgan": "1e8",
    "randPhi": "0",
    "phi": "0.2",
    "invRange": "0 1e3",
    "fixedSNPs": "1 10 1e-6",
    "n_SNPs": "0",
    "snpPositions": "500 550 600 650 700 750 800 850 900 950",
    "randomSample": "0",
    "tempRead": "10 0",
    "nCarriers": "10 0"
}

# Load parameters from the YAML files.
try:
    with open(OTHER_YAML, 'r') as f:
        other_params = yaml.safe_load(f)
except Exception as e:
    print(f"Error loading {OTHER_YAML}: {e}")
    sys.exit(1)

# Update the parameters dictionary with values from other.yaml (only for keys not updated from demes.yaml)
for key in other_params:
    if key in parameters and key not in ["popSizeVec", "speciation", "migRate"]:
        parameters[key] = str(other_params[key])

try:
    graph = demes.load(DEMES_YAML)
except Exception as e:
    print(f"Error loading {DEMES_YAML}: {e}")
    sys.exit(1)

# Update popSizeVec: For each deme except the one named "ancestor", extract the start_size from the first epoch.
pop_sizes = []
for deme in graph.demes:
    if deme.name != "ancestor":
        pop_sizes.append(str(deme.epochs[0].start_size))
parameters["popSizeVec"] = " ".join(pop_sizes)

# Update speciation: The format is "1 10000 0". Replace the middle value with the "ancestor" deme’s end_time.
ancestor_end_time = None
for deme in graph.demes:
    if deme.name == "ancestor":
        ancestor_end_time = deme.epochs[0].end_time
        break
if ancestor_end_time is not None:
    speciation_parts = parameters["speciation"].split()
    if len(speciation_parts) >= 3:
        speciation_parts[1] = str(ancestor_end_time)
        parameters["speciation"] = " ".join(speciation_parts)

# Update migRate: Set it to the migration rate from the first migration (if present) in the YAML.
if graph.migrations:
    parameters["migRate"] = str(graph.migrations[0].rate)

# --- Reintroduce the modify/update parameters function ---
def print_parameters(params):
    for key, value in params.items():
        print(f"{key}: {value}")

def modify_params(parameters):
    print("Current Parameters:")
    print_parameters(parameters)
    
    modify = input("\nDo you want to modify any parameters? (Y/N): ")
    if modify.lower() == 'y':
        user_input = input("\nEnter parameters to modify (e.g., --seed 2 --nruns 3): ")
        parser = argparse.ArgumentParser()
        for key in parameters.keys():
            if key in ["popSizeVec", "inv_freq", "speciation", "demography", "invRange", "fixedSNPs", "tempRead", "nCarriers", "snpPositions"]:
                parser.add_argument(f"--{key}", nargs='+')
            else:
                parser.add_argument(f"--{key}")
        args = parser.parse_args(user_input.split())
        
        for key, value in vars(args).items():
            if value is not None:
                if key in ["popSizeVec", "inv_freq", "speciation", "demography", "invRange", "fixedSNPs", "tempRead", "nCarriers", "snpPositions"]:
                    parameters[key] = " ".join(value)
                else:
                    parameters[key] = value
        print("\nUpdated Parameters:")
        print_parameters(parameters)
    else:
        print("\nNo changes made.")

# Call modify_params to allow interactive parameter updates.
modify_params(parameters)
# --- End modify/update function ---

# Prepare the ordered list of parameters to pass to the C++ executable.
keys_order = ["seed", "nruns", "kingman_coal", "drift_sim", "msOutput",
              "popSizeVec", "inv_freq", "speciation", "demography", "inv_age",
              "migRate", "BasesPerMorgan", "randPhi", "phi", "invRange", "fixedSNPs",
              "n_SNPs", "snpPositions", "randomSample", "tempRead", "nCarriers"]

args_list = [parameters[key] for key in keys_order]

print("\nFinal parameters being passed:")
for key in keys_order:
    print(f"{key}: {parameters[key]}")

# Run the C++ executable (assumed to be located at "./executables/labp_v19").
command = ["./executables/labp_v19"] + args_list
result = subprocess.run(command, capture_output=True, text=True)

print("\nC++ Program Output:\n", result.stdout)
if result.stderr:
    print("\nErrors:\n", result.stderr)
    with open("Output_log.txt", "w") as log_file:
        log_file.write(result.stderr)
