#!/usr/bin/env bash
# * Environment
rid="u6012627"
rip="notchpeak.chpc.utah.edu"
rdir="/uufs/chpc.utah.edu/common/home/koppelmans-group1/20230809_Kladblok\
/20251114_NMD_fMRI/20251115_fMRIprep/derivatives/fMRIPrep_25.2.3"
ldir="/datadisk/tmp/20251119_NMD_ClenchTask/data"
mkdir -p "${ldir}"


# * Copy over the data
rsync \
    --include="*/" \
    --include="*_ses-01_task-handdom_desc-confounds_timeseries.tsv" \
    --include="*_ses-01_task-handdom_space-fsLR_den-91k_bold.dtseries.nii" \
    --exclude="*" \
    -avmhe ssh \
    "${rid}@${rip}:${rdir}/**/ses-01/func/" \
    "${ldir}"/
