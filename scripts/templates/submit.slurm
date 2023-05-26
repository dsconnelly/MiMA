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

        cdo --reduce_dim -mermean -daymean -zonmean \
            -sellonlatbox,0.0,360.0,-5.0,5.0 -selname,u_gwf \
            ${yy}/atmos_4xdaily.nc ${yy}/qbo.nc

        cdo -mermean -daymean -zonmean \
            -sellonlatbox,0.0,360.0,59.0,61.0 -selname,u_gwf \
            ${yy}/atmos_4xdaily.nc ${yy}/ssw.nc

        cdo --reduce_dim -timmean -zonmean \
            -selname,u_gwf,v_gwf,t_gwf,gwfu_cgwd,gwfv_cgwd \
            ${yy}/atmos_4xdaily.nc ${yy}/clim.nc

        cdo -sellonlatbox,0.0,360.0,1.0,2.0 \
            -selname,u_gwf,t_gwf \
            ${yy}/atmos_4xdaily.nc ${yy}/tropical.nc

        cdo -sellonlatbox,0.0,360.0,39.0,41.0 \
            -selname,u_gwf,t_gwf \
            ${yy}/atmos_4xdaily.nc ${yy}/midlatitude.nc

        cdo -sellonlatbox,0.0,360.0,59.0,61.0 \
            -selname,u_gwf,t_gwf \
            ${yy}/atmos_4xdaily.nc ${yy}/polar.nc

    fi

    rm ${yy}/atmos_4xdaily.nc

done
rm .model.run

if SAVE_DATA; then

    cdo mergetime ??/qbo.nc qbo.nc
    cdo mergetime ??/ssw.nc ssw.nc
    cdo timmean -mergetime ??/clim.nc clim.nc

    cdo mergetime ??/tropical.nc tropical.nc
    cdo mergetime ??/midlatitude.nc midlatitude.nc
    cdo mergetime ??/polar.nc polar.nc
    cdo collgrid tropical.nc midlatitude.nc polar.nc covariance.nc

    rm ??/qbo.nc ??/ssw.nc ??/clim.nc
    rm ??/tropical.nc tropical.nc
    rm ??/midlatitude.nc midlatitude.nc
    rm ??/polar.nc polar.nc

fi

for yy in {01..N_YEARS}; do
    rm -rf ${yy}
done
