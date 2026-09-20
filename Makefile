PYTHON ?= python

.PHONY: install test classical quantum figures smoke smoke-quantum all

install:
	$(PYTHON) -m pip install -e ".[dev]"

test:
	PYTHONPATH=src $(PYTHON) -m pytest -v

classical:
	PYTHONPATH=src MPLCONFIGDIR=.mplconfig $(PYTHON) -m experiments.run_classical --config configs/classical.toml

quantum:
	PYTHONPATH=src MPLCONFIGDIR=.mplconfig $(PYTHON) -m experiments.run_quantum_single_pauli --config configs/quantum_single_pauli.toml

figures:
	PYTHONPATH=src MPLCONFIGDIR=.mplconfig $(PYTHON) -m experiments.make_figures --input results/data/classical_fer.csv --output results/figures/classical_fer_vs_ebn0.png
	PYTHONPATH=src MPLCONFIGDIR=.mplconfig $(PYTHON) -m experiments.make_figures --kind quantum-single-pauli --input results/data/quantum_single_pauli_failure.csv --output results/figures/quantum_single_pauli_failure_vs_sigma.png

smoke:
	PYTHONPATH=src MPLCONFIGDIR=.mplconfig $(PYTHON) -m experiments.run_classical --config configs/classical.toml --trials 1000

smoke-quantum:
	PYTHONPATH=src MPLCONFIGDIR=.mplconfig $(PYTHON) -m experiments.run_quantum_single_pauli --config configs/quantum_single_pauli.toml --trials 250 --no-figure

all: test classical quantum
