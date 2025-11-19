#!/usr/bin/env bash
# * Environment
rid="u6012627"
rip="notchpeak.chpc.utah.edu"
rdir="/uufs/chpc.utah.edu/common/home/koppelmans-group1/20230809_Kladblok\
/20251114_NMD_fMRI"
ldir="/datadisk/Utah/Kladblok/20201120_Margolis_NMD/20220504_Data_Processing/20251113_BayesfMRI/20251114_Preprocessing/20251117_QC"
mkdir -p "${ldir}"


# * Copy over the data
# ** fMRIprep
rsync \
    --include="*/" \
    --include="*.html" \
    --include="*.svg" \
    --exclude="*" \
    -avmhe ssh \
    "${rid}@${rip}:${rdir}/20251115_fMRIprep/derivatives/fMRIPrep_25.2.3" \
    "${ldir}"/

# ** MRIQC
rsync \
    --include="*/" \
    --include="*.html" \
    --include="*.svg" \
    --include="*.json" \
    --include="*.tsv" \
    --exclude="*" \
    -avmhe ssh \
    "${rid}@${rip}:${rdir}/20251116_MRIQC/mriqc_25.0.0rc0" \
    "${ldir}"/
