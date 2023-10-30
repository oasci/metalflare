SHELL := /usr/bin/env bash
PYTHON_VERSION := 3.11
PYTHON_VERSION_CONDENSED := 311
PACKAGE_NAME := metalflare
REPO_PATH := $(shell git rev-parse --show-toplevel)
PACKAGE_PATH := $(REPO_PATH)/02-methods/$(PACKAGE_NAME)
TESTS_PATH := $(REPO_PATH)/02-methods/99-tests/
CONDA_NAME := $(PACKAGE_NAME)-dev
CONDA_BASE_PATH = $(shell conda info --base)
CONDA_PATH := $(CONDA_BASE_PATH)/envs/$(CONDA_NAME)
CONDA := conda run -n $(CONDA_NAME)
DOCS_URL := https://metalflare.oasci.org

###   ENVIRONMENT   ###

.PHONY: conda-setup
conda-setup:
	- conda deactivate
	conda remove -y --name $(CONDA_NAME) --all
	conda create -y -n $(CONDA_NAME) python=$(PYTHON_VERSION)
	conda install -y conda-lock -n $(CONDA_NAME)
	conda install -y -c conda-forge poetry pre-commit tomli tomli-w -n $(CONDA_NAME)
	$(CONDA) pip install conda_poetry_liaison

.PHONY: write-conda-lock
write-conda-lock:
	- rm $(REPO_PATH)/conda-lock.yml
	$(CONDA) conda env export --from-history | grep -v "^prefix" > environment.yml
	$(CONDA) conda-lock -f environment.yml -p linux-64 -p osx-64 -p win-64
	$(CONDA) cpl-deps $(REPO_PATH)/pyproject.toml --env_path $(CONDA_PATH)
	$(CONDA) cpl-clean $(CONDA_PATH)

.PHONY: from-conda-lock
from-conda-lock:
	$(CONDA) conda-lock install -n $(CONDA_NAME) $(REPO_PATH)/conda-lock.yml
	$(CONDA) pip install conda_poetry_liaison
	$(CONDA) cpl-clean $(CONDA_PATH)

.PHONY: pre-commit-install
pre-commit-install:
	$(CONDA) pre-commit install

# Reads `pyproject.toml`, solves environment, then writes lock file.
.PHONY: poetry-lock
poetry-lock:
	$(CONDA) poetry lock --no-interaction
	$(CONDA) poetry export --without-hashes > requirements.txt

.PHONY: install
install:
	$(CONDA) poetry install --no-interaction

.PHONY: refresh
refresh: conda-setup from-conda-lock pre-commit-install install formatting validate


###   FORMATTING   ###

.PHONY: validate
validate:
	- $(CONDA) pre-commit run --all-files

.PHONY: formatting
formatting:
	- $(CONDA) isort $(PACKAGE_PATH)
	- $(CONDA) black --config pyproject.toml $(PACKAGE_PATH)



###   LINTING   ###

.PHONY: test
test:
	$(CONDA) pytest -c pyproject.toml --cov=$(PACKAGE_PATH) --cov-report=xml $(TESTS_PATH)

.PHONY: check-codestyle
check-codestyle:
	$(CONDA) isort --diff --check-only $(PACKAGE_PATH)
	$(CONDA) black --diff --check --config pyproject.toml $(PACKAGE_PATH)
	-$(CONDA) pylint $(PACKAGE_PATH)

.PHONY: mypy
mypy:
	-$(CONDA) mypy --config-file pyproject.toml ./

.PHONY: check-safety
check-safety:
	$(CONDA) poetry check
	$(CONDA) safety check --full-report
	$(CONDA) bandit -ll --recursive $(PACKAGE_PATH) $(TESTS_PATH)

.PHONY: lint
lint: check-codestyle mypy check-safety



###   CLEANING   ###

.PHONY: pycache-remove
pycache-remove:
	find . | grep -E "(__pycache__|\.pyc|\.pyo$$)" | xargs rm -rf

.PHONY: dsstore-remove
dsstore-remove:
	find . | grep -E ".DS_Store" | xargs rm -rf

.PHONY: mypycache-remove
mypycache-remove:
	find . | grep -E ".mypy_cache" | xargs rm -rf

.PHONY: ipynbcheckpoints-remove
ipynbcheckpoints-remove:
	find . | grep -E ".ipynb_checkpoints" | xargs rm -rf

.PHONY: pytestcache-remove
pytestcache-remove:
	find . | grep -E ".pytest_cache" | xargs rm -rf

.PHONY: coverage-remove
coverage-remove:
	find . | grep -E ".coverage" | xargs rm -rf

.PHONY: build-remove
build-remove:
	rm -rf build/

.PHONY: cleanup
cleanup: pycache-remove dsstore-remove mypycache-remove ipynbcheckpoints-remove pytestcache-remove psi-remove coverage-remove



###   WEBSITE   ###

.PHONY: website
website:
	rm -rf ./website/html/
	$(CONDA) sphinx-build -nT ./ ./website/html/
	touch ./website/html/.nojekyll

.PHONY: website-versioned
website-versioned:
	rm -rf ./website/html/
	$(CONDA) sphinx-multiversion -nT ./ ./website/html/
	touch ./website/html/.nojekyll

	# Create html redirect to main
	echo "<head>" > ./website/html/index.html
	echo "  <meta http-equiv='refresh' content='0; URL=$(DOCS_URL)/main/index.html'>" >> ./website/html/index.html
	echo "</head>" >> ./website/html/index.html

.PHONY: open-website
open-website:
	xdg-open ./website/html/index.html 2>/dev/null
