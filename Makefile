ROOT := $(CURDIR)
SRC := $(ROOT)/src
EXEC_SMC := $(ROOT)/executables/labp_smc
EXEC_ARG := $(ROOT)/executables/labp_v21
PYTHON := $(if $(wildcard .venv/bin/python),.venv/bin/python,python3)
TSKIT_C := $(ROOT)/third_party/tskit/c
TSKIT_BUILD := $(ROOT)/build/tskit
TSKIT_OBJS := $(TSKIT_BUILD)/core.o $(TSKIT_BUILD)/tables.o $(TSKIT_BUILD)/kastore.o
DEPS_FLAGS :=
ifneq ($(strip $(CONDA_PREFIX)),)
DEPS_FLAGS := -I "$(CONDA_PREFIX)/include" -L "$(CONDA_PREFIX)/lib" -Wl,-rpath,"$(CONDA_PREFIX)/lib"
endif
.PHONY: run run-smc run-arg compile compile-smc compile-arg delete-op clean-generated

run-smc:
	@cd "$(ROOT)" && "$(PYTHON)" src/param.py --config src/parameters_smc.yaml

run-arg:
	@cd "$(ROOT)" && "$(PYTHON)" src/param.py --config src/parameters_arg.yaml

compile: compile-arg compile-smc

compile-smc:
	@mkdir -p "$(TSKIT_BUILD)"
	@cc -std=c99 -O3 -I "$(TSKIT_C)" -I "$(TSKIT_C)/subprojects/kastore" -c "$(TSKIT_C)/tskit/core.c" -o "$(TSKIT_BUILD)/core.o"
	@cc -std=c99 -O3 -I "$(TSKIT_C)" -I "$(TSKIT_C)/subprojects/kastore" -c "$(TSKIT_C)/tskit/tables.c" -o "$(TSKIT_BUILD)/tables.o"
	@cc -std=c99 -O3 -I "$(TSKIT_C)" -I "$(TSKIT_C)/subprojects/kastore" -c "$(TSKIT_C)/subprojects/kastore/kastore.c" -o "$(TSKIT_BUILD)/kastore.o"
	@cd "$(SRC)" && g++ -std=c++14 -O3 -g -fno-omit-frame-pointer \
		$(DEPS_FLAGS) \
		-I "$(TSKIT_C)" \
		-I "$(TSKIT_C)/subprojects/kastore" \
		argnode.cpp chromosome.cpp poisevents.cpp simulate.cpp sitenode.cpp migprob.cpp world.cpp parameters.cpp smc.cpp chromrecomb.cpp smc_helpers.cpp treemod.cpp simulate_smc.cpp \
		"$(TSKIT_BUILD)/core.o" "$(TSKIT_BUILD)/tables.o" "$(TSKIT_BUILD)/kastore.o" \
		-lboost_random -lboost_system -lboost_math_c99 -lyaml-cpp \
		-o labp_smc
	@mv -f "$(SRC)/labp_smc" "$(EXEC_SMC)"

compile-arg:
	@cd "$(SRC)" && g++ -std=c++14 -O3 -g -fno-omit-frame-pointer \
		$(DEPS_FLAGS) \
		argnode.cpp chromosome.cpp poisevents.cpp simulate.cpp sitenode.cpp migprob.cpp world.cpp parameters.cpp main.cpp chromrecomb.cpp \
		-lboost_random -lboost_system -lboost_math_c99 -lyaml-cpp \
		-o labp_v21
	@mv -f "$(SRC)/labp_v21" "$(EXEC_ARG)"

delete-op:
	@cd "$(ROOT)" && find . -maxdepth 1 -type f -name 'genetree_*' -delete
