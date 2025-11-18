#!/usr/bin/env bash
# * Environment
# ** IO paths
base="/uufs/chpc.utah.edu/common/home/koppelmans-group1"
idir="${base}/20230809_Kladblok/20251114_NMD_fMRI/20251114_BIDS_prepped"
odir="${base}/20230809_Kladblok/20251114_NMD_fMRI/20251115_fMRIprep"
sdir="${odir}/scripts"
mkdir -p "${sdir}"/{logs,jobs}
# ** fMRIprep version, license, and container location
fmriprep_version="25.2.3"
simg="${base}/20230809_Software/SingularityImages/fmriprep-${fmriprep_version}.simg"
fslicense="${base}/20230809_Software/license/freesurfer_license.txt"
# ** Go to log folder
cd "${sdir}/logs" || exit


# * List subjects to be processed
subjects=(
    $(find \
          "${idir}" \
          -mindepth 1 \
          -maxdepth 1 \
          -type d \
          -iname "sub-*" \
          -exec basename {} \; \
          | sort)
)


# * Create the fMRIprep script for each subject
# Loop over all subjects and create a script for fMRIprep
# for sub in sub-NMD001; do
for sub in ${subjects[@]}; do

    # ** Define job
    subject_job="${sdir}/jobs/fMRIprep_${sub}.slurm"

    # ** Create job via heredoc
    cat <<-EOF2> "${subject_job}"
	#!/bin/bash
	#SBATCH --account=koppelmans-np
	#SBATCH --mail-user=u6012627@utah.edu
	#SBATCH --partition=koppelmans-shared-np
	#SBATCH --job-name=${sub//sub-/}_fmriprep
	#SBATCH --nodes=1
	#SBATCH --ntasks=12
	#SBATCH --mem=64G
	#SBATCH --time=24:00:00
	#SBATCH -o log_${sub}_fmriprep-%j.out-%N
	#SBATCH -e log_${sub}_fmriprep-%j.err-%N


	# * Load module(s)
	module load apptainer


	# * Environment
	sub="${sub}"
	base="${base}"
	idir="\${base}/20230809_Kladblok/20251114_NMD_fMRI/20251114_BIDS_prepped"
	odir="\${base}/20230809_Kladblok/20251114_NMD_fMRI/20251115_fMRIprep"
	sdir="\${odir}/scripts/jobs"
	fmriprep_version="${fmriprep_version}"
	fmriprep="\${base}/20230809_Software/SingularityImages/fmriprep-\${fmriprep_version}.simg"
	fslicense="\${base}/20230809_Software/license/freesurfer_license.txt"


	# * Run fMRIprep
	# ** Create output folders
	mkdir -p \${odir}/derivatives/fMRIPrep_\${fmriprep_version}
	workdir=\$(mktemp -d)
	echo "Working directory: \${workdir}"

	# ** Run fMRIprep
	apptainer \\
	   run \\
	   --cleanenv \\
	   --bind \${idir}:/bids:ro \\
	   --bind \${odir}/derivatives/fMRIPrep_\${fmriprep_version}:/out \\
	   --bind \${fslicense}:/license.txt:ro \\
	   \${fmriprep} \\
	   /bids \\
	   /out \\
	   participant \\
	   --participant-label \${sub} \\
	   --output-spaces MNI152NLin6Asym:res-2 \\
	   --cifti-output 91k \\
	   --notrack \\
	   --fs-license-file /license.txt \\
	   --skip_bids_validation \\
	   --nprocs 12 \\
	   --fs-subjects-dir /out/\${sub}/freesurfer \\
	   --work-dir="\${workdir}" \\
	   --clean-workdir \\
	   --stop-on-first-crash \\
	   --write-graph \\
	   --random-seed 357 \\
	   -vv

exit
EOF2

    # ** Submit job to the CHPC slurm via sbatch
    sbatch "${subject_job}"

done
