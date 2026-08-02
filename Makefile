# ExonDomainCompare - convenience targets

PYTHON ?= ./.venv/bin/python

.PHONY: setup serve doctor test lint build release-check check architecture architecture-diagrams architecture-guide architecture-check

setup:
	./scripts/setup_local.sh

serve:
	./scripts/start_local.sh

doctor:
	$(PYTHON) scripts/edc.py doctor --redact-paths

test:
	$(PYTHON) -m pytest -p no:cacheprovider -q

lint:
	npm --prefix webapp/frontend run lint

build:
	npm --prefix webapp/frontend run build

release-check:
	$(PYTHON) scripts/release/check_public_release.py

check: release-check test lint build

## Build the full Architecture Atlas (all diagram SVGs + PDFs and the guide).
architecture:
	$(PYTHON) scripts/docs/build_architecture_atlas.py

## Build only the diagram SVGs + PDFs.
architecture-diagrams:
	$(PYTHON) scripts/docs/build_architecture_atlas.py --diagrams

## Build only the Architecture Guide PDF.
architecture-guide:
	$(PYTHON) scripts/docs/build_architecture_atlas.py --guide

## Verify the diagram toolchain (d2, rsvg-convert, typst) is installed.
architecture-check:
	$(PYTHON) scripts/docs/build_architecture_atlas.py --check
