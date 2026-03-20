#!/bin/bash -e

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=3
#SBATCH --time=24:00:00
#SBATCH --mem=64G
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
N_YEARS="25"

for yy in $(seq -f "%02g" 1 $N_YEARS); do
    mpirun --bind-to none \
        -np $SLURM_NTASKS \
	env OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK \
	./mima.x

    for fname in *.nc.0000; do
        head=${fname%.0000}
        $ccomb -r $head $head.*
    done

    mkdir -p ${yy}
    mv *.nc ${yy}

    tar -czf ${yy}/${yy}.restart.tar.gz -C RESTART .
    mv RESTART/* INPUT
done

module load cdo/2.5.0
cdo -zonmean -timselmean,28, -selname,u -mergetime ??/atmos_out.nc qbo-weekly.nc
