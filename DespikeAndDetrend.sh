#!/usr/bin/env bash
# 
# despike and detrend after slice timing and motion correction
#
#
# full-path-to-NIFTI-image, without extension
#vol=/Volumes/Untitled/RestOlfAnaly/SubjectData/${subj}/${sess}/${sess}.feat/filtered_func_data
# despiked data ${vol}_dspk; detrended and despiked data ${vol}_rlt
vol=$1

# AFNI environment variable setup
export PATH=$PATH:/Users/gz/Downloads/afni 
export DYLD_FALLBACK_LIBRARY_PATH=/Users/gz/Downloads/afni 
export AFNI_3dDespike_NEW=YES

# 1 - do not delete existing files and re-run; 0 - overwrite existing files
rerun=1
vol=$( ${FSLDIR}/bin/remove_ext ${vol} )
volname=$( basename $vol )
voldir=$( dirname $vol )

## 1. Despike
if [[ $rerun -eq 0 ]]; then
    printf "Despiking ...\n"
    [ -f ${vol}_dspk.nii.gz ] && rm ${vol}_dspk.nii.gz
    3dDespike -prefix ${vol}_dspk.nii.gz ${vol}.nii.gz
    
else
    if [[ -f ${vol}_dspk.nii.gz ]]; then
        printf "Despike has alreay been run, skipped.\n"
    else
        printf "Despiking ...\n"
        3dDespike -prefix ${vol}_dspk.nii.gz ${vol}.nii.gz
    fi    
fi

## 2. Detrending for linear trends
if [[ $rerun -eq 0 ]]; then
    [ -f ${vol}_rlt.nii.gz ] && rm ${vol}_rlt.nii.gz
    printf "Removing linear trend in timeseries\n"
    3dTcat -rlt+ -prefix ${vol}_rlt.nii.gz ${vol}_dspk.nii.gz
    
else
    if [[ -f ${vol}_rlt.nii.gz ]]; then
        printf "Linear trend has already been removed, skipped.\n"
    else
        printf "Removing linear trend in timeseries\n"
        3dTcat -rlt+ -prefix ${vol}_rlt.nii.gz ${vol}_dspk.nii.gz
    fi    
fi
