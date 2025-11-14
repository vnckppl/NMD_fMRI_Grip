#!/usr/bin/env bash
#SBATCH --account=koppelmans-np
#SBATCH --mail-user=u6012627@utah.edu
#SBATCH --partition=koppelmans-shared-np
#SBATCH --job-name=build_fmriprep
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=24G
#SBATCH --time=3:00:00
#SBATCH -o build_fmriprep_sif-%j.out-%N
#SBATCH -e build_fmriprep_sif-%j.err-%N

# * Environment
odir="/uufs/chpc.utah.edu/common/home/koppelmans-group1/20230809_Software\
/SingularityImages"

# * Load singularity/apptainer
module load apptainer

# * Version
container_version="25.2.3"

# * Environment
apptainer \
    build \
    "${odir}/fmriprep-${container_version}.simg" \
    "docker://nipreps/fmriprep:${container_version}"

exit
