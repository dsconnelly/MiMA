#!/bin/bash -e

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00
#SBATCH --mem=48G
#SBATCH --job-name=RUN_NAME
#SBATCH --out=RUN_DIR/slurm.out
#SBATCH --mail-user=USER@nyu.edu
#SBATCH --mail-type=END,FAIL

cd RUN_DIR
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

source MIMA_DIR/scripts/modules.sh
ccomb=MIMA_DIR/build/mppnccombine

# This variable sets the number of years to integrate for.
# Be sure to keep it as two digits - that is, you need to pad
# it with a zero if you integrate for fewer than ten years.
N_YEARS="05"

for yy in {01..$N_YEARS}; do
    mpirun -np $SLURM_NTASKS ./mima.x

    for fname in *.nc.0000; do
        head=${fname%.0000}
        mpirun -np 1 $ccomb -r $head $head.*
    done

    mkdir -p ${yy}
    mv *.nc ${yy}

    tar -czf ${yy}/${yy}.restart.tar.gz -C RESTART .
    mv RESTART/* INPUT
done

