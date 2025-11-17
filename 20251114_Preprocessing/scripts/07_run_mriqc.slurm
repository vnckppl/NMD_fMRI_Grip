#!/usr/bin/env bash
#SBATCH --account=owner-guest
#SBATCH --mail-user=u6012627@utah.edu
#SBATCH --partition=notchpeak-shared-guest
#SBATCH --job-name=run_mriqc
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH -o run_mriqc_sif-%j.out-%N
#SBATCH -e run_mriqc_sif-%j.err-%N

# * Environment
base="/uufs/chpc.utah.edu/common/home/koppelmans-group1"
# ** MRIQC version and container location
mriqc_version="25.0.0rc0"
mriqc="${base}/20230809_Software/SingularityImages/mriqc-${mriqc_version}.simg"
# ** IO paths
idir="${base}/20230809_Kladblok/20251114_NMD_fMRI/20251114_BIDS_prepped"
odir="${base}/20230809_Kladblok/20251114_NMD_fMRI/20251116_MRIQC/mriqc_${mriqc_version}"
mkdir -p "${odir}"


# * Load singularity/apptainer
module load apptainer


# * Run MRIQC
apptainer \
   run \
   --cleanenv \
   --bind "${idir}":/data:ro \
   --bind "${odir}":/out \
   ${mriqc} \
   /data \
   /out \
   participant \
   --nprocs 32 \
   --omp-nthreads 4 \
   --float32 \
   --write-graph \
   --no-sub \
   --notrack


exit
