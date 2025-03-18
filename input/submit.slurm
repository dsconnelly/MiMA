#!/bin/bash -e

#SBATCH --nodes=8
#SBATCH --ntasks-per-node=4
#SBATCH --time=2:00:00
#SBATCH --mem=64G
#SBATCH --job-name=RUN_NAME
#SBATCH --out=RUN_DIR/slurm.out
#SBATCH --mail-user=dsconnelly@nyu.edu
#SBATCH --mail-type=END,FAIL

cd RUN_DIR
./modules.sh

for yy in {01..01}; do
    mpirun -n ${SLURM_NTASKS} ./mima.x

    ccomb=/home/dsc7746/MiMA/bin/mppnccombine.greene
    for fname in *.nc.0000; do
        head=${fname%.0000}
        $ccomb -r $head $head.*
    done

    mkdir -p ${yy}
    mv *.nc ${yy}

    tar -czf ${yy}/${yy}.restart.tar.gz -C RESTART .
    mv RESTART/* INPUT
    
done
