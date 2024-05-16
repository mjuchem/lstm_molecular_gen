#!/bin/bash -e

SRC_DIR="${SRC_DIR:-"$(git rev-parse --show-toplevel)"}"
NOTEBOOKS_DIR="${NOTEBOOKS_DIR:-"${SRC_DIR}"}"
LAB_ADDR="${LAB_ADDR:-127.0.0.1}"
LAB_PORT="${LAB_PORT:-8080}"

docker run -it --rm \
  --env LAB_PORT="${LAB_PORT}" \
  --publish "${LAB_ADDR}:${LAB_PORT}:${LAB_PORT}" \
  --volume "${NOTEBOOKS_DIR}:/src/notebooks" \
  "$@" \
  lstm_molecular_gen
