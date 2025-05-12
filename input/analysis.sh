#!/bin/bash

python ()
{
    singularity exec --overlay ${overlay}:ro ${image} \
        /bin/bash -c "source /ext3/env.sh; python ${*}"
}

MIMA_DIR="/home/dsc7746/MiMA"
python "${MIMA_DIR}/analysis.py" "$@"