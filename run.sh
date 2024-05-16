#!/bin/bash -e

SRC_DIR="${SRC_DIR:-"$(git rev-parse --show-toplevel)"}"
NOTEBOOKS_DIR="${NOTEBOOKS_DIR:-"${SRC_DIR}"}"
LAB_ADDR="${LAB_ADDR:-127.0.0.1}"
LAB_PORT="${LAB_PORT:-52019}"

DOCKER_ARGS=( \
  --env LAB_PORT="${LAB_PORT}" \
  --publish "${LAB_ADDR}:${LAB_PORT}:${LAB_PORT}" \
  --volume "${NOTEBOOKS_DIR}:/src/notebooks" \
)

if [ -e /dev/dri ]; then
  DOCKER_ARGS+=(--device /dev/dri:/dev/dri)
fi

docker run -it --rm \
  "${DOCKER_ARGS[@]}" \
  "$@" \
  lstm_molecular_gen
