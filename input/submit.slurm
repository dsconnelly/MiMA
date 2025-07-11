#!/bin/bash -e

#SBATCH --nodes=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00
#SBATCH --mem=64G
#SBATCH --job-name=RUN_NAME
#SBATCH --out=RUN_DIR/slurm.out
#SBATCH --mail-user=dsconnelly@nyu.edu
#SBATCH --mail-type=END,FAIL

cd RUN_DIR
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load netcdf-c/intel/4.7.4
module load netcdf-fortran/intel/4.5.3
module load openmpi/intel/4.1.1
module load cdo/intel/1.9.10

for yy in {01..05}; do
    srun ./mima.x

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

cdo -timselmean,28 -mergetime ??/atmos_4xdaily.nc weekly.nc
srun -n 1 ./analysis.sh weekly.nc plot-climatology:DJF plot-climatology:JJA plot-qbo