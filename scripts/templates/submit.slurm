#!/bin/bash -e

#SBATCH --nodes=8
#SBATCH --ntasks-per-node=4
#SBATCH --time=96:00:00
#SBATCH --mem=64G
#SBATCH --job-name=CASE_NAME-online
#SBATCH --out=CASE_DIR/slurm.out
#SBATCH --mail-user=dsconnelly@nyu.edu
#SBATCH --mail-type=END,FAIL

module purge
module load netcdf-c/intel/4.7.4
module load netcdf-fortran/intel/4.5.3
module load openmpi/intel/4.0.5
module load cdo/intel/1.9.10

cd CASE_DIR
if HAS_BASE; then
    cp BASE_DIR/INPUT/* CASE_DIR/INPUT
fi

cat << "EOF" > .model.run
#!/bin/bash
module purge
module load netcdf-c/intel/4.7.4
module load netcdf-fortran/intel/4.5.3
module load openmpi/intel/4.0.5
module load cdo/intel/1.9.10
source /ext3/activate.sh
export LD_LIBRARY_PATH=/lib:/lib64:${LD_LIBRARY_PATH}
export PYTHONWARNINGS=ignore::UserWarning
./mima.x
EOF

chmod +x .model.run
for yy in {01..N_YEARS}; do
    mpirun -n $SLURM_NTASKS singularity exec \
        --bind /usr,/lib,/lib64 \
        --overlay ${overlay}:ro \
        $image \
        ./.model.run

    ccomb=/home/dsc7746/MiMA/bin/mppnccombine.greene
    for fname in *.nc.0000; do
        head=${fname%.0000}
        $ccomb -r $head $head.*
    done

    mkdir -p ${yy}
    mv *.nc ${yy}

    tar -czf ${yy}/${yy}.restart.tar.gz -C RESTART .
    mv RESTART/* INPUT

    if SAVE_DATA; then

        cdo --reduce_dim -daymean -zonmean \
            -selname,ps,u_gwf,v_gwf,t_gwf,bf_cgwd,gwfu_cgwd,gwfv_cgwd \
            ${yy}/atmos_4xdaily.nc ${yy}/zonal_mean.nc

        cdo -sellonlatbox,0.0,360.0,1.0,2.0 \
            -selname,u_gwf,gwfu_cgwd \
            ${yy}/atmos_4xdaily.nc ${yy}/tropics.nc

        cdo -sellonlatbox,0.0,360.0,59.0,61.0 \
            -selname,u_gwf,gwfu_cgwd \
            ${yy}/atmos_4xdaily.nc ${yy}/vortex.nc

    fi

    rm ${yy}/atmos_4xdaily.nc

done
rm .model.run

if SAVE_DATA; then

    cdo mergetime ??/zonal_mean.nc zonal_mean.nc
    cdo mergetime ??/tropics.nc tropics.nc
    cdo mergetime ??/vortex.nc vortex.nc
    cdo collgrid tropics.nc vortex.nc covariance.nc

    rm ??/zonal_mean.nc
    rm ??/tropics.nc tropics.nc
    rm ??/vortex.nc vortex.nc

fi

for yy in {01..N_YEARS}; do
    rm -rf ${yy}
done
