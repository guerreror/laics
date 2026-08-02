# InvLaBP

Coalescent simulation of locally adapted inversions

InvLaBP simulates coalescent histories for sequences within a polymorphic chromosomal inversion under migration–selection balance (Guerrero et al., 2012). The inversion primarily alters recombination: heterokaryotypes exchange genetic material (“gene flux”) at a different rate than homokaryotypes.

---

## Quickstart

    # 1) Install deps (macOS/Homebrew)
    brew install boost yaml-cpp python
    python -m pip install demes numpy pyyaml

    # 2) Edit demes.yaml
    #    param.py automatically writes src/migration_matrices.json from demes.yaml

    # 3) Edit simulation parameters
    #    - src/parameters_smc.yaml for SMC
    #    - src/parameters_arg.yaml for ARG

    # 4) Build (from src/)
    g++ -std=c++14 -O3 \
      -I /opt/homebrew/opt/boost/include \
      -I /opt/homebrew/opt/yaml-cpp/include \
      -L /opt/homebrew/opt/boost/lib \
      -L /opt/homebrew/opt/yaml-cpp/lib \
      argnode.cpp chromosome.cpp poisevents.cpp simulate.cpp sitenode.cpp migprob.cpp world.cpp parameters.cpp main.cpp chromrecomb.cpp \
      -lboost_random -lboost_system -lboost_math_c99 -lyaml-cpp \
      -o labp_v19

    # 5) Place the executable inside the executables folder

    # 6) Run (from repo root)
    make run      # SMC
    make run-arg  # ARG

---

## Newly Added Features

- YAML-driven configuration (`src/parameters_smc.yaml`, `src/parameters_arg.yaml`)
- Demes-aware, time-varying migration via `demes.yaml` → `migration_matrices.json`
- Multiple speciation merges and demography changes on explicit schedules
- Python launcher (`src/param.py`) that parses configs, prints summaries, and can emit ms-style data

---

## Requirements

### OS / Tooling

- macOS or Linux (tested primarily on macOS)
- C++14 compiler (`clang++` / `g++`)
- macOS: Xcode Command Line Tools (for profiling purposes)

### Libraries

- Boost (headers + libs): `random`, `system`, `math_c99`
- `yaml-cpp`
- Python 3.9+
  - `demes`
  - `numpy`
  - `pyyaml`

> Adjust include/lib paths if your Boost or yaml-cpp lives elsewhere.

---

## Profiling Build

    cd src
    g++ -std=c++14 -O3 -g -fno-omit-frame-pointer …….. -o labp_v19

---

## Run (Standard Workflow)

Run (standard workflow):

1. **Edit `demes.yaml`**  
   Define populations, sizes, epochs, and migrations.

2. **Let `param.py` generate migration matrices**

   `param.py` automatically reads `demes.yaml` and writes `src/migration_matrices.json`.

3. **Edit the parameter YAML accordingly**  
   Use `src/parameters_smc.yaml` for SMC or `src/parameters_arg.yaml` for ARG.

4. **Run the Python launcher from the root folder**

       python src/param.py --config src/parameters_smc.yaml

   Or use:

       make run
       make run-arg

---

## Configuration

SMC knobs live in `src/parameters_smc.yaml`; ARG knobs live in `src/parameters_arg.yaml`.

- `Seed`: RNG seed (`0` = random; the effective seed is printed)
- `NumberOfReplicates`: number of simulation replicates

- `KingmanCoal`:
  - `1` = event-driven Kingman approximation
  - `0` = generation-by-generation

- `DriftSimulation`:
  - (Drift simulation mode; see code / `param.py` for details)

- `MSOutput`:  
  `1` to emit `outLABP*.sites` and `outLABP*.stats`

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

- `TargetSNPs`:
  SMC only. Emit diagnostics for selected horizontal positions.

- `GeneConversionRate`:
  SMC local GC rectangle height.

- `DoubleRecombinationRate`:
  SMC DR triangle peak height.
