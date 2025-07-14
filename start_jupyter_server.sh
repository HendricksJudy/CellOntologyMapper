#!/bin/bash
# This script sets up a Python virtual environment and launches JupyterLab so the notebooks
# in this repository can be used via a web browser.

set -euo pipefail

ENV_DIR="${ENV_DIR:-venv}"

# Create the virtual environment if it doesn't exist
if [ ! -d "$ENV_DIR" ]; then
  python3 -m venv "$ENV_DIR"
fi

# Activate the environment
source "$ENV_DIR/bin/activate"

# Install required packages
pip install --upgrade pip
# Install all dependencies from requirements.txt
if [ -f "requirements.txt" ]; then
  pip install -r requirements.txt
else
  pip install git+https://github.com/Starlitnightly/CellOntologyMapper.git
  pip install jupyterlab
fi

# Launch JupyterLab accessible on port 8888
exec jupyter lab --ip 0.0.0.0 --port 8888 --no-browser
