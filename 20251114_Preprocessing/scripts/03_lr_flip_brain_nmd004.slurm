#!/usr/bin/env bash
#SBATCH --account=koppelmans-np
#SBATCH --mail-user=u6012627@utah.edu
#SBATCH --partition=koppelmans-shared-np
#SBATCH --job-name=lr_flip_NMD004
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=12G
#SBATCH --time=1:00:00
#SBATCH -o lr_flip_NMD004-%j.out-%N
#SBATCH -e lr_flip_NMD004-%j.err-%N

# This script makes a backup of the data of subject NMD004, the only left
# handed subject in the data set. It then flips the MRI data left-to-right to
# ensure that we will be able to do group level analyses for the dominant hand.

# * Environment
dir="/uufs/chpc.utah.edu/common/home/koppelmans-group1/20230809_Kladblok\
/20251114_NMD_fMRI/20240715_BIDS"
cd "${dir}" || exit


# * Load FSL module
module load fsl


# * Backup the original data
# Only do this if the backup file does not exists. This will prevent overwriting
# the true original data.
bkpfile="${dir}/sub-NMD004_non_flipped_bkp.gz"
if [ ! -f "${bkpfile}" ]; then
    echo "Backing up data for sub-NMD004"
    tar -zcf "${bkpfile}" "sub-NMD004"
    echo $(basename "${bkpfile}") >> .bidsignore
    echo "scripts/" >> .bidsignore
else
    echo "Data for sub-NMD004 already backed up. Skipping..."
fi


# * List NMD004's nifti files
ddir="${dir}/sub-NMD004"
files=($(find "${ddir}" -iname "*.nii.gz" | sort))


# * Flip each nifti file
for file in ${files[@]}; do

    fbase=$(echo "${file}" | sed "s#${ddir}#..#g")
    echo "$(tput setaf 3)Left-Right flipping:$(tput sgr0) ${fbase}"
    fslswapdim "${file}" -x y z "${file}"

done


exit
