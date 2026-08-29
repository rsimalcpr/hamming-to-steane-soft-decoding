PYTHON ?= python

.PHONY: install test classical figures smoke all

install:
	$(PYTHON) -m pip install -e ".[dev]"

test:
	PYTHONPATH=src $(PYTHON) -m unittest discover -s tests -v

classical:
	PYTHONPATH=src MPLCONFIGDIR=.mplconfig $(PYTHON) -m experiments.run_classical --config configs/classical.toml

figures:
	PYTHONPATH=src MPLCONFIGDIR=.mplconfig $(PYTHON) -m experiments.make_figures --input results/data/classical_fer.csv --output results/figures/classical_fer_vs_ebn0.png

smoke:
	PYTHONPATH=src MPLCONFIGDIR=.mplconfig $(PYTHON) -m experiments.run_classical --config configs/classical.toml --trials 1000

all: test classical
