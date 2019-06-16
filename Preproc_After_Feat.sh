#!/bin/bash
#
# Preprocessing of fMRI after Feat
#

wkpath=~/Desktop/Data/RestfMRI/Dopar
# full-path-to-subject list 
subjlist=~/Desktop/Data/RestfMRI/Dopar/subject_list.txt

# spatial smoothing paramters
sigma=3   

# band-pass filtering frequencies, in Hz
lp=0.1
hp=0.008

while read subj
do
    featdir=${wkpath}/${subj}/rest/rest.feat
    fsmaksdir=${wkpath}/${subj}/anat/NuiMask

    start_time=$( date )
    echo "  Starting batch at ${start_time}"

    # filename after despike and detrend: filtered_func_data_rlt
    printf "    Batch: despike and detrend ...\n"
    ./DespikeAndDetrend.sh ${featdir}/filtered_func_data

    # extract nuisance time series
    #vol=${featdir}/filtered_func_data_rlt
    #roidir=${fsmaksdir}
    #refimg=${featdir}/example_func
    #reg=${featdir}/reg/highres2example_func.mat
    #roisavedir=${featdir}/ROIs
    tsdir=${featdir}/TimeSeries/Nuisance

    volname=filtered_func_data_rlt
    printf "    Batch: Extracting nuisance regressors ...\n"
    ./Extract_Brain_Nuisance.sh ${featdir} ${volname} ${fsmaksdir}

    # combine motion and CSF, WM, whole-brain
    paste -d ' ' ${featdir}/mc/prefiltered_func_data_mcf.par ${tsdir}/wholebrain_wm_csf.txt > ${tsdir}/nuisance_regressor.txt

    # regres out nuisance
    printf "    Batch: Regressing out nuisance regressor ...\n"
    vol=${featdir}/filtered_func_data_rlt
    design=${tsdir}/nuisance_regressor.txt

    vol2save=${featdir}/filtered_func_data_rlt_nui
    fsl_glm -i $vol -d ${design} -o ${vol2save}_beta --out_res=${vol2save}_res --demean
    fslmaths ${vol} -Tmean -add ${vol2save}_res ${vol2save}_nuicorr -odt float
    [ -f ${vol2save}_res.nii.gz ] && rm ${vol2save}_res.nii.gz
	[ -f ${vol}.nii.gz ] && rm ${vol}.nii.gz

    # global intensity normalization
    printf "    Batch: Global intensity normalization ...\n"
    fslmaths ${vol2save}_nuicorr -ing 10000 ${vol2save}_nuicorr_volsorm -odt float

    # register to MNI space
    printf "    Batch: Register to MNI space ...\n"
    reg=${featdir}/reg/example_func2standard.mat
    
    # band-pass filtering and smooth for functional connectivity analysis
    printf "    Batch: Band-pass filtering and smoothing ...\n"
    [ -f ${vol2save}_nuicorr_volsorm_filt.nii.gz ] && rm ${vol2save}_nuicorr_volsorm_filt.nii.gz
    3dFourier -lowpass ${lp} -highpass ${hp} -prefix ${vol2save}_nuicorr_volsorm_filt.nii.gz ${vol2save}_nuicorr_volsorm.nii.gz

    # registration to standard brain
    flirt -in ${vol2save}_nuicorr_volsorm_filt -ref ${FSLDIR}/data/standard/MNI152_T1_2mm_brain -applyxfm -init ${reg} -out ${vol2save}_nuicorr_volsorm_filt2std
    fslmaths ${vol2save}_nuicorr_volsorm_filt2std -kernel gauss ${sigma} -fmean ${vol2save}_nuicorr_volsorm_filt2std_s${sigma} -odt float

    clean_data=0
    if [[ clean_data -eq 0 ]]; then
        [ -f ${vol2save}_nuicorr.nii.gz ] && rm ${vol2save}_nuicorr.nii.gz
        [ -f ${vol2save}_nuicorr_volsorm.nii.gz ] && rm ${vol2save}_nuicorr_volsorm.nii.gz
        [ -f ${vol2save}_nuicorr_volsorm_filt.nii.gz ] && rm ${vol2save}_nuicorr_volsorm_filt.nii.gz
		[ -f ${vol2save}_nuicorr_volsorm_filt2std.nii.gz ] && rm ${vol2save}_nuicorr_volsorm_filt2std.nii.gz
    fi

    end_time=$( date )
    echo "Batch ended at ${end_time}"
	
done < ${subjlist}

