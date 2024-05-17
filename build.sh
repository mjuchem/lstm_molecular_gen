#!/bin/bash -e

SRC_DIR="${SRC_DIR:-"$(git rev-parse --show-toplevel)"}"
REPO_NAME="$(basename "${SRC_DIR}")"

image_name="${REPO_NAME}"

(set -x; \
  DOCKER_BUILDKIT=1 \
    docker build \
      -t "${image_name}" \
      docker \
)
