#!/bin/bash -e

DOCKER_BUILDKIT=1 \
  docker build \
    -t lstm_molecular_gen \
    docker
