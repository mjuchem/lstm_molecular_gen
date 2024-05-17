LAB_ADDR ?= 127.0.0.1
LAB_PORT ?= 52019

.PHONY: all build run interactive start restart stop kill

all: build

build:
	./build.sh

interactive: build
	LAB_ADDR="${LAB_ADDR}" \
	LAB_PORT="${LAB_PORT}" \
		./run.sh \
			--entrypoint bash

run: build
	LAB_ADDR="${LAB_ADDR}" \
	LAB_PORT="${LAB_PORT}" \
		./run.sh

start: build
	LAB_ADDR="${LAB_ADDR}" \
	LAB_PORT="${LAB_PORT}" \
		./run.sh --start

restart:
	./run.sh --restart

stop:
	./run.sh --stop

kill:
	./run.sh --kill
