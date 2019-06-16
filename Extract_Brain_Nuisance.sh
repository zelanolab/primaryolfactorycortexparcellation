#!/bin/bash
#
# Extract white-matter, csf and whole-brain EPI time series using anatomical masks
#
# featdir: feat directory
# 

#rest=rfmriNasal
#subjdir=/Users/gzc4893/Desktop/Tmp/TestGLM/subj002
#roidir=/Volumes/Untitled/RestOlfAnaly/SubjectData/ro_subj002/anat/FSMasks

#featdir=${subjdir}/${rest}/${rest}+++.feat
#vol=${featdir}/filtered_func_data_rlt
#roisavedir=${featdir}/${rest}/ROIs
#refimg=${featdir}/example_func
#reg=${featdir}/reg/highres2example_func.mat
#tsdir=${featdir}/TimeSeries/Nuisance

featdir=$1
vol4ts=$2
fsmaskdir=$3

vol=${featdir}/${vol4ts}
roidir=${fsmaskdir}
roisavedir=${featdir}/ROIs
refimg=${featdir}/example_func
reg=${featdir}/reg/highres2example_func.mat
tsdir=${featdir}/TimeSeries/Nuisance

[ -d $tsdir ] || mkdir -p ${tsdir}
[ -d $roisavedir ] || mkdir -p ${roisavedir}

echo ${vol} > ${tsdir}/vol4nuisance.txt

flirt -in ${roidir}/brainmask -ref ${refimg} -applyxfm -init ${reg} -interp nearestneighbour -out ${roisavedir}/brainmask
flirt -in ${roidir}/csf_erode -ref ${refimg} -applyxfm -init ${reg} -interp nearestneighbour -out ${roisavedir}/csf_erode
flirt -in ${roidir}/wm_erode -ref ${refimg} -applyxfm -init ${reg} -interp nearestneighbour -out ${roisavedir}/wm_erode 


#fslmeants -i $vol -m ${roisavedir}/brainmask -o ${tsdir}/wholebrain.txt
fslmeants -i $vol -m ${featdir}/mask -o ${tsdir}/wholebrain.txt
fslmeants -i $vol -m ${roisavedir}/csf_erode -o ${tsdir}/csf_erode.txt
fslmeants -i $vol -m ${roisavedir}/wm_erode -o ${tsdir}/wm_erode.txt
paste ${tsdir}/wholebrain.txt ${tsdir}/wm_erode.txt ${tsdir}/csf_erode.txt > ${tsdir}/wholebrain_wm_csf.txt
[ -f ${tsdir}/wm_erode.txt ] && rm ${tsdir}/wm_erode.txt
[ -f ${tsdir}/csf_erode.txt ] && rm ${tsdir}/csf_erode.txt
[ -f ${tsdir}/wholebrain.txt ] && rm ${tsdir}/wholebrain.txt


