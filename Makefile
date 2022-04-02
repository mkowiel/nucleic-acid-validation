VENV_DIR := .venv
PYTHON := ${VENV_DIR}/bin/python

venv:
	test -d $(VENV_DIR) || python3 -m venv $(VENV_DIR)
	. $(VENV_DIR)/bin/activate & $(PYTHON) -m pip install --upgrade pip setuptools wheel tox

tox: venv
	. $(VENV_DIR)/bin/activate & $(PYTHON) -m tox

test: tox

black:
	. $(VENV_DIR)/bin/activate & $(PYTHON) -m black naval tests

clean:
	rm -rf .pytest_cache
	rm -rf $(VENV_DIR)
	rm -rf .tox .mypy_cache .coverage
