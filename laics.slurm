#!/bin/bash
#SBATCH --job-name=test
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

source $(conda info --base)/etc/profile.d/conda.sh
conda activate /usr/local/usrapps/guerrerolab/baturne5

# cd src/
# g++ -std=c++14 -O3 -g -fno-omit-frame-pointer \
# -I $CONDA_PREFIX/include \
# -L $CONDA_PREFIX/lib \
# -Wl,-rpath,$CONDA_PREFIX/lib \
# argnode.cpp chromosome.cpp poisevents.cpp simulate.cpp sitenode.cpp migprob.cpp world.cpp parameters.cpp smc.cpp chromrecomb.cpp smc_helpers.cpp treemod.cpp simulate_smc.cpp \
# -lboost_random \
# -lboost_system \
# -lboost_math_c99 \
# -lyaml-cpp \
# -o labp_smc 

# mv labp_smc ../executables/


# g++ -std=c++14 -O3 -g -fno-omit-frame-pointer \
#   -I $CONDA_PREFIX/include \
#   -L $CONDA_PREFIX/lib \
#   -Wl,-rpath,$CONDA_PREFIX/lib \
#   argnode.cpp chromosome.cpp poisevents.cpp simulate.cpp sitenode.cpp migprob.cpp world.cpp parameters.cpp main.cpp chromrecomb.cpp \
#   -lboost_random \
#   -lboost_system \
#   -lboost_math_c99 \
#   -lyaml-cpp \
#   -o labp_v21

# mv labp_v21 ../executables/
# cd ../

#python src/mig_matrix_demes.py
python src/param.py --config src/parameters_smc.yaml
