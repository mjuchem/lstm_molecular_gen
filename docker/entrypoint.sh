#!/bin/bash -xe

source /srv/venv/bin/activate

exec jupyter-lab \
  --allow-root --no-browser \
  --ip="0.0.0.0" \
  --port="${LAB_PORT:-8080}" \
  --notebook-dir="/src/notebooks" \
  --NotebookApp.token='' \
  --NotebookApp.password='' \
  "$@"
