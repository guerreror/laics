ROOT := $(CURDIR)
SRC := $(ROOT)/src
EXEC := $(ROOT)/executables/labp_smc
PYTHON := $(if $(wildcard .venv/bin/python),.venv/bin/python,python3)
.PHONY: run run-smc run-arg compile delete-op clean-generated run-tree

run-smc:
	@cd "$(ROOT)" && "$(PYTHON)" src/param.py --config src/parameters_smc.yaml

run-arg:
	@cd "$(ROOT)" && "$(PYTHON)" src/param.py --config src/parameters_arg.yaml

run-tree:
	@cd "$(ROOT)" && "$(PYTHON)" tools/render_tree_from_csv_snapshot.py

compile:
	@cd "$(SRC)" && g++ -std=c++14 -O3 -g -fno-omit-frame-pointer \
		-I /opt/homebrew/opt/boost/include \
		-I /opt/homebrew/opt/yaml-cpp/include \
		-L /opt/homebrew/opt/boost/lib \
		-L /opt/homebrew/opt/yaml-cpp/lib \
		argnode.cpp chromosome.cpp poisevents.cpp simulate.cpp sitenode.cpp migprob.cpp world.cpp parameters.cpp smc.cpp chromrecomb.cpp smc_helpers.cpp treemod.cpp simulate_smc.cpp \
		-lboost_random -lboost_system -lboost_math_c99 -lyaml-cpp \
		-o labp_smc
	@mv -f "$(SRC)/labp_smc" "$(EXEC)"

delete-op:
	@cd "$(ROOT)" && find . -maxdepth 1 -type f -name 'genetree_*' -delete
