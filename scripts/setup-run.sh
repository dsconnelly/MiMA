#!/bin/bash -e

RUN_DIR=$(realpath "$1")
RUN_NAME=$(basename "$RUN_DIR")

mkdir -p "$RUN_DIR/RESTART"
cp -r "$PWD/input/"* "$RUN_DIR"
cp "$PWD/build/mima.x" "$RUN_DIR"

sed -i "s|RUN_DIR|$RUN_DIR|" "$RUN_DIR/submit.slurm"
sed -i "s|RUN_NAME|$RUN_NAME|" "$RUN_DIR/submit.slurm"
sed -i "s|MIMA_DIR|$PWD|" "$RUN_DIR/submit.slurm"
sed -i "s|USER|$USER|" "$RUN_DIR/submit.slurm"

