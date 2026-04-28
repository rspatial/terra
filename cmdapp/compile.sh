#!/bin/bash
# Wrapper: builds the standalone library + ./terra via Makefile.
# See Makefile for targets (make lib, make clean, make -j4).

set -e
cd "$(dirname "$0")"
exec make "$@"
