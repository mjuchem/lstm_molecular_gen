#!/bin/bash -e

SRC_DIR="${SRC_DIR:-"$(git rev-parse --show-toplevel)"}"
REPO_NAME="$(basename "${SRC_DIR}")"

image_name="${REPO_NAME}"
container_name="${REPO_NAME}"

NOTEBOOKS_DIR="${NOTEBOOKS_DIR:-"${SRC_DIR}"}"
LAB_ADDR="${LAB_ADDR:-127.0.0.1}"
LAB_PORT="${LAB_PORT:-52019}"

docker_args=()

mode=run

while [ "$#" -gt 0 ]; do
  arg="$1"; shift
  case "${arg}" in
    --run|--start|--restart|--stop|--kill)
      mode="${arg:2}"
      ;;

    *)
      docker_args+=("${arg}")
      ;;
  esac
done

case "${mode}" in
  run)
    docker_args+=(run -it --rm)
    ;;

  start)
    docker_args+=( \
      run
      --name "${container_name}"
      --detach
      --restart unless-stopped
    )
    ;;

  restart)
    docker_args+=(container restart "${container_name}")
    ;;

  stop)
    docker_args+=(container stop "${container_name}")
    ;;

  kill)
    docker_args+=(container rm -f "${container_name}")
    ;;

  *)
    echo "invalid mode: ${mode}"
    exit 1
    ;;
esac


if [ "${mode}" = "start" ] || [ "${mode}" = "run" ]; then
  docker_args+=( \
    --env LAB_PORT="${LAB_PORT}"
    --publish "${LAB_ADDR}:${LAB_PORT}:${LAB_PORT}"
    --volume "${NOTEBOOKS_DIR}:/src/notebooks"
  )

  if [ -e /dev/dri ]; then
    docker_args+=(--device /dev/dri:/dev/dri)
  fi

  docker_args+=("${image_name}")
fi

(set -x; docker "${docker_args[@]}")
