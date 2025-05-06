#!/bin/bash

RUN_DIR="$1"
mkdir -p "${RUN_DIR}/RESTART"
MIMA_DIR="/home/dsc7746/MiMA"

cp -r "${MIMA_DIR}/input/"* "${RUN_DIR}"
cp "${MIMA_DIR}/exp/exec.greene/mima.x" "${RUN_DIR}"
sed -i "s|RUN_DIR|${RUN_DIR}|" "${RUN_DIR}/recompile.sh"

RUN_NAME=$(basename "${RUN_DIR}")
sed -i "s|RUN_DIR|${RUN_DIR}|" "${RUN_DIR}/submit.slurm"
sed -i "s|RUN_NAME|${RUN_NAME}|" "${RUN_DIR}/submit.slurm"

shift
for arg in "$@"; do
    if [[ "$arg" =~ ^--([^=]+)=(.*)$ ]]; then
        name="${BASH_REMATCH[1]}"
        value="${BASH_REMATCH[2]}"
        num=$(grep -n "\<${name}\>" "${RUN_DIR}/input.nml" | cut -d: -f1)

        if [ -n "$num" ]; then
            orig=$(sed -n "${num}p" "${RUN_DIR}/input.nml")
            ending=$(echo "$orig" | \
                sed -E 's/.*[^[:space:]]([,/])[[:space:]]*$/\1/')

            line="${name} = ${value}${ending}"
            sed -i "${num}s/.*/${line}/" "${RUN_DIR}/input.nml"
        fi
    fi
done