LAB_ADDR ?= 127.0.0.1
LAB_PORT ?= 52019

.PHONY: all build run interactive

all: build

build:
	./build.sh

run: build
	LAB_ADDR="${LAB_ADDR}" \
	LAB_PORT="${LAB_PORT}" \
		./run.sh

interactive: build
	LAB_ADDR="${LAB_ADDR}" \
	LAB_PORT="${LAB_PORT}" \
		./run.sh \
			--entrypoint bash
