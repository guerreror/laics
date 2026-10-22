# InvLaBP

Coalescent simulation of locally adapted inversions

InvLaBP simulates coalescent histories for sequences within a polymorphic chromosomal inversion under migration–selection balance (Guerrero et al., 2012). The inversion primarily alters recombination: heterokaryotypes exchange genetic material (“gene flux”) at a different rate than homokaryotypes.

---

## Quickstart

    # 1) Install deps (macOS/Homebrew)
    brew install yaml-cpp python
    python -m pip install demes numpy pyyaml

    # 2) Add matching demes and parameter YAML files
    #    - SMC_input/ for SMC
    #    - ARG_input/ for ARG
    #    param.py automatically writes migration_matrices.json into each result directory

    # 4) Build (from src/)
    g++ -std=c++14 -O3 \
      -I /opt/homebrew/opt/yaml-cpp/include \
      -L /opt/homebrew/opt/yaml-cpp/lib \
      argnode.cpp chromosome.cpp poisevents.cpp simulate.cpp sitenode.cpp migprob.cpp world.cpp parameters.cpp main.cpp chromrecomb.cpp \
      -lyaml-cpp \
      -o labp_v19

    # 5) Place the executable inside the executables folder

    # 6) Run (from repo root)
    make run-smc  # SMC
    make run-arg  # ARG

---

## Newly Added Features

- Separate YAML input directories (`SMC_input/`, `ARG_input/`)
- Demes-aware, time-varying migration via `demes.yaml` → `migration_matrices.json`
- Multiple speciation merges and demography changes on explicit schedules
- Python launcher (`src/param.py`) that parses configs and prints summaries

---

## Requirements

### OS / Tooling

- macOS or Linux (tested primarily on macOS)
- C++14 compiler (`clang++` / `g++`)
- macOS: Xcode Command Line Tools (for profiling purposes)

### Libraries

- `yaml-cpp`
- Python 3.9+
  - `demes`
  - `numpy`
  - `pyyaml`

> Adjust include/lib paths if your yaml-cpp installation lives elsewhere.

---

## Profiling Build

    cd src
    g++ -std=c++14 -O3 -g -fno-omit-frame-pointer …….. -o labp_v19

---

## Run (Standard Workflow)

Run (standard workflow):

1. **Add YAML inputs**
   Put matching `demes_*.yaml` and `parameters_*.yaml` pairs in `SMC_input/`. ARG uses `ARG_input/demes_arg.yaml` and `ARG_input/parameters_arg.yaml`.

2. **Let `param.py` generate migration matrices**

   `param.py` reads the selected demes file and writes `migration_matrices.json` into the run's output directory.

3. **Run from the root folder**

   Set `TreeOutputDirectory` in an SMC parameter YAML to an optional scratch or
   output location. Trees are written to `<TreeOutputDirectory>/<run-name>/smc_trees`;
   an empty value keeps them in the normal batch result directory.

       make run-smc

   Or use:

       make run-arg

---

## Configuration

SMC configurations live in `SMC_input/`; ARG configuration lives in `ARG_input/`.

- `Seed`: RNG seed (`0` = random; the effective seed is printed)
- `NumberOfReplicates`: number of simulation replicates

- `KingmanCoal`:
  - `1` = event-driven Kingman approximation
  - `0` = generation-by-generation

- `DriftSimulation`:
  - (Drift simulation mode; see code / `param.py` for details)

- `InversionFrequency`:  
  Initial inversion frequency per population

- `InversionAge`:  
  If `> 0`, step the inversion to loss outside origin context at that time

- `BasesPerMorgan`:
  - Bases per Morgan (scaling from physical bp to recombination units)

- `RandPhi`:
  - ARG only. `1` to sample `log10(Phi)` uniformly in a range.
  - `0` to use fixed `Phi`

- `Phi`:  
  ARG fixed gene-flux when `RandPhi = 0`

- `InversionRange`:  
  Inversion span in bp (scaled internally)

  ```yaml
  InversionRange: "L_bp R_bp"
  ```

- `FixedSNPs`:
  ARG only. `FixedSNPs = 1 ...` means exactly fixed marker count.

- `RandomSNPs`:
  ARG only. `1` means random marker placement.

- `SiteNodePositions`:  
  Site positions in bp.

- `RandomSample`: 
  `1` → one per-pop total sample; carriers drawn binomially by `InversionFrequency`  
  `0` → one pair per pop: `<standard> <inverted>` (exact)

- `Samples`: 
  Per-pop sample definitions.

  ```yaml
  Samples:
    pop0: "2 0"
    pop1: "2 0"
  ```

- `GeneConversionRate`:
  SMC local GC rectangle height.

- `DoubleRecombinationRate`:
  SMC DR triangle peak height.
