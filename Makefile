.PHONY: clean test black tox venv

VENV_DIR := .venv
PYTHON := ${VENV_DIR}/bin/python

venv:
	test -d $(VENV_DIR) || python3 -m venv $(VENV_DIR)
	$(PYTHON) -m pip install --upgrade pip setuptools wheel
	$(PYTHON) -m pip install --upgrade tox black isort

tox: venv
	. $(VENV_DIR)/bin/activate & $(PYTHON) -m tox

black:
	. $(VENV_DIR)/bin/activate & $(PYTHON) -m black naval tests

isort:
	. $(VENV_DIR)/bin/activate & $(PYTHON) -m isort naval tests

test: tox

clean:
	rm -rf .pytest_cache
	rm -rf $(VENV_DIR)
	rm -rf .tox .mypy_cache .coverage
