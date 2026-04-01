# ** Create screenshots
# * Environment
base="/Users/u1478302/Library/CloudStorage/Box-Box/Seungmi Lee UROP/20251123_1stlevel"
tdir=$(mktemp -d)


# * Subjects
subjects=$(
    find "${base}" \
         -mindepth 1 \
         -maxdepth 1 \
         -type d \
         -iname "sub-NMD0*" \
         -exec basename {} \; \
        | sort
        )
subjects=(sub-NMD001) # Test


# * Function for display
# ** 02_BayesfMRI output
myss_02_bayes() {

    # * Environment
    sub="${1}"
    idir="${base}/${sub}"
    iidir="${idir}/02_BayesfMRI"
    ifile_nii="${iidir}/${sub}_bglm_1st_level.sep.nii"
    ofile_nii="${iidir}/${sub}_02_BayesfMRI_subcortical.png"

    # ** Only run create screenshots if the imaging 1st level data exists
    if [ -f "${ifile_nii}" ]; then

        # *** Announce
        echo "- [X] BAYES fMRI data for ${sub} found"

        # *** Create screenshot of subcortical results
        fsleyes \
            render \
            -of "${ofile_nii}" \
            --size 1680 800 \
            --scene lightbox \
            --zrange 0.05 0.45 \
            --sliceSpacing 0.02 \
            --hideCursor \
            -std1mm \
            "${ifile_nii}" \
            -cm red-yellow \
            -nc blue-lightblue \
            -dr 1 2 \
            --gamma +0.5

        # *** Create screenshots of cortical results
        HS=("lh" "rh")            # Hemisphere
        CR=("-90" "90")           # Camera Rotation: degrees
        CP=("lat" "med")          # Camera Postion: label
        LP=("0 0 -100" "0 0 100") # Light position

        # *** Loop over hemispheres
        for H in ${HS[@]}; do

            # **** Hemisphere letter
            HL=$(echo ${H:0:1} | tr [:lower:] [:upper:])
            #HL=$(echo "${H^^}" | cut -c 1)

            # **** Loop over camera positions
            for D in 0 1; do

                # ***** Output file
                ofile_gii="${tdir}/${sub}_02_BayesfMRI_surface_${H}_${CP[${D}]}.png"

                # ***** Render screenshot
                fsleyes \
                    render \
                    -of "${ofile_gii}" \
                    --size 900 800 \
                    -ds world \
                    --scene 3d \
                    --bgColour 1 1 1 \
                    --performance 3 \
                    --colourBarLocation bottom \
                    --cameraRotation "${CR[${D}]}" 0 0 \
                    --lightPos ${LP[${D}]} \
                    --zoom 100 \
                    --hideLegend \
                    --hideCursor \
                    "${iidir}/S1200.${HL}.inflated_MSMAll.10k_fs_LR.surf.gii" \
                    --colour 0.5 0.5 0.5 \
                    -vd "${iidir}/${sub}_bglm_1st_level.sep.${HL}.func.gii" \
                    -dr 1 2 \
                    -cr 1 100000 \
                    -cm red-yellow \
                    -nc blue-lightblue

            done
        done

        # *** Combine surface figures
        magick \
            "${tdir}/${sub}_02_BayesfMRI_surface_lh_med.png" \
            "${tdir}/${sub}_02_BayesfMRI_surface_lh_lat.png" \
            "${tdir}/${sub}_02_BayesfMRI_surface_rh_med.png" \
            "${tdir}/${sub}_02_BayesfMRI_surface_rh_lat.png" \
	    +append \
	    +repage \
            "${iidir}/${sub}_02_BayesfMRI_surface.png"

    else \

        echo "- [ ] BAYES fMRI data for ${sub} *NOT* found"

    fi
}


# ** 03_Activation output
myss_03_activation() {

    # * Environment
    sub="${1}"
    idir="${base}/${sub}"
    iidir="${idir}/03_Activation"
    surfdir="${idir}/02_BayesfMRI"
    ifile_nii="${iidir}/${sub}_bglm_1st_level.sep.nii"
    ofile_nii="${iidir}/${sub}_03_Activation_subcortical.png"

    # ** Only run create screenshots if the imaging 1st level data exists
    if [ -f "${ifile_nii}" ]; then

        # *** Announce
        echo "- [X] ACTIVATION data for ${sub} found"

        # *** Prepare images for display
        # **** Create nifti mask
        tdir=$(mktemp -d)
        fslmaths "${ifile_nii}" -bin "${tdir}/mask.nii.gz"
        fslmaths "${ifile_nii}" -thr 2 -bin "${tdir}/activation.nii.gz"

        # *** Create screenshot of subcortical results
        fsleyes \
            render \
            -of "${ofile_nii}" \
            --size 1680 800 \
            --scene lightbox \
            --zrange 0.05 0.45 \
            --sliceSpacing 0.02 \
            --hideCursor \
            "${FSLDIR}/data/standard/MNI152_T1_1mm.nii.gz" \
            "${tdir}/mask.nii.gz" -a 90 \
            "${tdir}/activation.nii.gz" -ot mask --maskColour 1 0 0

        # *** Create screenshots of cortical results
        HS=("lh" "rh")            # Hemisphere
        CR=("-90" "90")           # Camera Rotation: degrees
        CP=("lat" "med")          # Camera Postion: label
        LP=("0 0 -100" "0 0 100") # Light position

        # *** Loop over hemispheres
        for H in ${HS[@]}; do

            # **** Hemisphere letter
            HL=$(echo ${H:0:1} | tr [:lower:] [:upper:])
            #HL=$(echo "${H^^}" | cut -c 1)

            # **** Loop over camera positions
            for D in 0 1; do

                # ***** Output file
                ofile_gii="${tdir}/${sub}_03_Activation_surface_${H}_${CP[${D}]}.png"

                # ***** Render screenshot
                fsleyes \
                    render \
                    -of "${ofile_gii}" \
                    --size 900 800 \
                    -ds world \
                    --scene 3d \
                    --bgColour 1 1 1 \
                    --performance 3 \
                    --colourBarLocation bottom \
                    --cameraRotation "${CR[${D}]}" 0 0 \
                    --lightPos ${LP[${D}]} \
                    --zoom 100 \
                    --hideLegend \
                    --hideCursor \
                    "${surfdir}/S1200.${HL}.inflated_MSMAll.10k_fs_LR.surf.gii" \
                    --colour 0.5 0.5 0.5 \
                    -vd "${iidir}/${sub}_bglm_1st_level.sep.${HL}.label.gii" \
                    -cr 1 100000 \
                    -cm red
            done
        done

        # *** Combine surface figures
        magick \
            "${tdir}/${sub}_03_Activation_surface_lh_med.png" \
            "${tdir}/${sub}_03_Activation_surface_lh_lat.png" \
            "${tdir}/${sub}_03_Activation_surface_rh_med.png" \
            "${tdir}/${sub}_03_Activation_surface_rh_lat.png" \
	    +append \
	    +repage \
            "${iidir}/${sub}_03_Activation_surface.png"

    else \

        echo "- [ ] ACTIVATION data for ${sub} *NOT* found"

    fi
}


# * Loop over subjects and create screenshots
for sub in ${subjects[@]}; do

    myss_02_bayes "${sub}"
    myss_03_activation "${sub}"

done


** Create Montage
# * Environment
base="/Users/u1478302/Library/CloudStorage/Box-Box/Seungmi Lee UROP/20251123_1stlevel"
odir="${base}/montage"
mkdir -p "${odir}"


# * Montage BayesfMRI
# ** Glob files
files_bayes=$(
    find "${base}" \
         -mindepth 3 \
         -maxdepth 3 \
         -type f \
         -iname "*BayesfMRI_surface.png" \
        | sort
           )

# ** Create montage
montage \
    ${files_bayes} \
    -geometry '900x200+50+0' \
    -tile 2x^ \
    "${odir}/BayesfMRI_surface.png"


# * Montage Activation
# ** Glob files
files_activation=$(
    find "${base}" \
         -mindepth 3 \
         -maxdepth 3 \
         -type f \
         -iname "*Activation_surface.png" \
        | sort
     )

# ** Create montage
montage \
    ${files_activation} \
    -geometry '900x200+50+0' \
    -tile 2x^ \
    "${odir}/Activation_surface.png"