#!/bin/bash

RUN_DIR="$1"
mkdir -p "${RUN_DIR}/RESTART"
MIMA_DIR="/home/dsc7746/MiMA"

cp -r "${MIMA_DIR}/input/"* "${RUN_DIR}"
cp "${MIMA_DIR}/exp/exec.greene/mima.x" "${RUN_DIR}"
cp "${MIMA_DIR}/submit.slurm" "${RUN_DIR}"

run_name=$(basename "${RUN_DIR}")
sed -i "s|RUN_DIR|${RUN_DIR}|" "${RUN_DIR}/submit.slurm"
sed -i "s|RUN_NAME|${run_name}|" "${RUN_DIR}/submit.slurm"
