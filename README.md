# InvLaBP

Coalescent simulation of locally adapted inversions

InvLaBP simulates coalescent histories for sequences within a polymorphic chromosomal inversion under migration–selection balance (Guerrero et al., 2012). The inversion primarily alters recombination: heterokaryotypes exchange genetic material (“gene flux”) at a different rate than homokaryotypes.

---

## Quickstart

    # 1) Install deps (macOS/Homebrew)
    brew install boost yaml-cpp python
    python -m pip install demes numpy pyyaml

    # 2) Generate migration matrices from your demes graph
    python src/mig_matrix_demes.py  # Outputs the migration_matrices.json file

    # 3) Edit simulation parameters
    #    - parameters.yaml (points to migration_matrices.json and sets all sim knobs)

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
    python src/param.py

---

## Newly Added Features

- YAML-driven configuration (`parameters.yaml`)
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

2. **Generate migration matrices using the provided `demes.yaml` file**

       python src/mig_matrix_demes.py

   This produces a time-stamped `migration_matrices.json` which is used by the simulation software.

3. **Edit `parameters.yaml` accordingly**  
   Make sure it points to `migration_matrices.json` and sets all simulation parameters.

4. **Run the Python launcher from the root folder**

       python src/param.py

---

## Configuration (Key Fields in `parameters.yaml` File)

All knobs live in `parameters.yaml`. Names map directly to internal parameters.

- `seed`: RNG seed (`0` = random; the effective seed is printed)
- `nruns`: number of replicates of the simulation

- `kingman_coal`:
  - `1` = event-driven Kingman approximation
  - `0` = generation-by-generation

- `drift_sim`:
  - (Drift simulation mode; see code / `param.py` for details)

- `msOutput`:  
  `1` to emit `outLABP*.sites` and `outLABP*.stats`

- `inv_freq`:  
  Initial inversion frequency per population

- `inv_age`:  
  If `> 0`, step the inversion to loss outside origin context at that time

- `BasesPerMorgan`:
  - Bases per Morgan (scaling from physical bp to recombination units)

- `randPhi`:
  - `1` to sample `log10(phi)` uniformly in `phi_range = [min, max]`
  - `0` to use fixed `phi`

- `phi`:  
  Fixed gene-flux when `randPhi = 0`

- `invRange`:  
  Inversion span in bp (scaled internally)

  ```yaml
  invRange: [L_bp, R_bp]
- `fixedS`:
  `fixedS = 1` → exactly `n_SNPs` markers

- `n_SNPs`: 
  Number of site node positions.

- `randSNP`:
  `1` → uniform placement in `[snpPositions[0], snpPositions[-1]]`

- `snpPositions`:  
  Site positions in bp.

- `randomSample`: 
  `1` → one per-pop total sample; carriers drawn binomially by `inv_freq`  
  `0` → one pair per pop: `<standard> <inverted>` (exact)

- `nCarriers / tempRead`: 
  As printed by the summary from `src/param.py`.