# ExonDomainCompare - convenience targets

PYTHON ?= ./.venv/bin/python
EDC ?= ./.venv/bin/edc

.PHONY: setup serve doctor test lint build release-check check

setup:
	./scripts/setup_local.sh

serve:
	./scripts/start_local.sh

doctor:
	$(EDC) doctor --redact-paths

test:
	$(PYTHON) -m pytest -p no:cacheprovider -q

lint:
	npm --prefix webapp/frontend run lint

build:
	npm --prefix webapp/frontend run build

release-check:
	$(PYTHON) scripts/release/check_public_release.py

check: release-check test lint build
