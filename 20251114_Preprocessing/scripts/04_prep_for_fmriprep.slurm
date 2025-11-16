#!/usr/bin/env bash
#SBATCH --account=koppelmans-np
#SBATCH --mail-user=u6012627@utah.edu
#SBATCH --partition=koppelmans-shared-np
#SBATCH --job-name=prep_for_fmriprep
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=12G
#SBATCH --time=2:00:00
#SBATCH -o prep_for_fmriprep-%j.out-%N
#SBATCH -e prep_for_fmriprep-%j.err-%N

# This script:
# a) reconstructs the T1 image from the multi-echo MPRAGE images
# b) copies relevant fMRI (task + rest) data and field maps from the BIDS folder

# * Environment
base="/uufs/chpc.utah.edu/common/home/koppelmans-group1/20230809_Kladblok\
/20251114_NMD_fMRI"
idir="${base}/20240715_BIDS"
odir="${base}/20251114_BIDS_prepped"
mkdir -p "${odir}"


# * Load FSL module
module load fsl


# * Function for creating the T1 image
mkt1 () {
    # ** Environment
    sub="${1}"
    iidir="${2}/${sub}/ses-01/anat"
    oodir="${3}/${sub}/ses-01/anat"
    mkdir -p "${oodir}"
    tfile=$(mktemp --suffix=".nii.gz")

    # ** Find the four MEMPRAGE images
    images=$(
        find "${iidir}" \
             -mindepth 1 \
             -maxdepth 1 \
             -type f \
             -iname "*multiecho*echo-*nii.gz" \
             | sort
    )

    # ** Merge the four MEMPRAGE images into a single 4D file
    fslmerge -t "${tfile}" $(echo ${images})

    # ** Convert the 4D image to a single T1w image
    fslmaths \
             "${tfile}" \
             -sqr -Tmean -sqrt \
             "${oodir}/${sub}_ses-01_T1w.nii.gz"
}


# * List subjects
subjects=($(
    find "${idir}" \
         -mindepth 1 \
         -maxdepth 1 \
         -type d \
         -iname "sub-*" \
         -exec basename {} \; \
         | sort
))

# ** [for testing]
# subjects=("sub-NMD001")


# * Loop over subjects to select data and create the T1w images
for sub in ${subjects[@]}; do

    # ** Reconstruct T1
    mkt1 "${sub}" "${idir}" "${odir}"

    # ** Create output folders
    mkdir -p "${odir}/${sub}/ses-01"/{func,fmap}

    # ** Some shortcuts
    andir_i="${idir}/${sub}/ses-01/anat"
    andir_o="${odir}/${sub}/ses-01/anat"
    fudir_i="${idir}/${sub}/ses-01/func"
    fudir_o="${odir}/${sub}/ses-01/func"
    fmdir_i="${idir}/${sub}/ses-01/fmap"
    fmdir_o="${odir}/${sub}/ses-01/fmap"

    # ** Copy over resting state fMRI data
    cp \
        "${fudir_i}/${sub}_ses-01_task-rest_part-mag_bold.nii.gz" \
        "${fudir_o}/${sub}_ses-01_task-rest_bold.nii.gz"
    cp \
        "${fudir_i}/${sub}_ses-01_task-rest_part-mag_bold.json" \
        "${fudir_o}/${sub}_ses-01_task-rest_bold.json"

    # ** Copy over task-based fMRI data
    # This finds all task data, sorts it in reverse order (so, last run on top),
    # and then picks the top item.
    nfile=$(
        find "${fudir_i}" -iname "*handdom*nii.gz*" \
             | sort -r \
             | head -1
         )
    cp "${nfile}" "${fudir_o}/${sub}_ses-01_task-handdom_bold.nii.gz"
    cp "${nfile/nii.gz/json}" "${fudir_o}/${sub}_ses-01_task-handdom_bold.json"

    # ** Copy over the field maps
    # Loop over reversed phase encoding directions and tasks
    for dir in AP PA; do
        for task in rest handdom; do

            # *** Define input file
            nfile=$(
                find "${fmdir_i}" -iname "*acq-${task}*${dir}_epi.nii.gz*" \
                   | sort -r \
                   | head -1
                 )

            # *** Copy over the nifti file
            cp "${nfile}" \
                "${fmdir_o}/${sub}_ses-01_acq-${task}_dir-${dir}_epi.nii.gz"

            # *** Copy over the json file
            cp "${nfile/nii.gz/json}" \
                "${fmdir_o}/${sub}_ses-01_acq-${task}_dir-${dir}_epi.json"

        done
    done

    # ** Copy over the T2w data
    # Only 39 out of 40 subjects have T2w data.
    cp \
        "${andir_i}/${sub}_ses-01_rec-bias_T2w.nii.gz" \
        "${andir_o}/${sub}_ses-01_T2w.nii.gz"
    cp \
        "${andir_i}/${sub}_ses-01_rec-bias_T2w.json" \
        "${andir_o}/${sub}_ses-01_T2w.json"

done


# * Copy some BIDS data meta files
for i in participants.json participants.tsv; do
    cp "${idir}/${i}" "${odir}"
done


# * Fill out the datset_description.json file
cat <<EOF>${odir}/dataset_description.json
{
  "Name": "Neuro Motor Depression study: hand clench task",
  "BIDSVersion": "1.7.0",
  "DatasetType": "raw",
  "License": "",
  "Authors": [
    "Vincent Koppelmans, (PI)"
  ],
  "Funding": [
    "Ben B & Iris Margolis Foundation"
  ],
  "EthicsApprovals": ["University of Utah, IRB_00139658"],
  "ReferencesAndLinks": [
    "https://doi.org/10.3389/fpsyt.2025.1624776",
    "https://doi.org/10.1123/jmld.2024-0105"
  ],
  "DatasetDOI": "doi:"
}
EOF


# * Add .bidsignore for the scripts folder
echo "scripts/" > "${odir}/.bidsignore"


exit
