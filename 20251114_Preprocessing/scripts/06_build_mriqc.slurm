#!/usr/bin/env bash
#SBATCH --account=owner-guest
#SBATCH --mail-user=u6012627@utah.edu
#SBATCH --partition=notchpeak-shared-guest
#SBATCH --job-name=build_mriqc
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=24G
#SBATCH --time=3:00:00
#SBATCH -o build_mriqc_sif-%j.out-%N
#SBATCH -e build_mriqc_sif-%j.err-%N

# * Environment
odir="/uufs/chpc.utah.edu/common/home/koppelmans-group1/20230809_Software\
/SingularityImages"

# * Load singularity/apptainer
module load apptainer

# * Version
container_version="25.0.0rc0"

# * Environment
apptainer \
    build \
    "${odir}/mriqc-${container_version}.simg" \
    "docker://nipreps/mriqc:${container_version}"

exit
