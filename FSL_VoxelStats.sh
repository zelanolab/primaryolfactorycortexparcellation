#!/usr/bin/env bash
#
# count the number of significant voxels within each region using HarvardOxford cortical and subcortical atlases
#
# Input
#   a p value image
#   p value threshold
#   atlas name
#   roi radius

################################################################
## Start of parameter setup
# p value file, e.g.: *_tfce_corrp_*.nii.gz
#pfile=/Volumes/Untitled/RestOlfAnaly/Results/sup_par_s6mmR2Z/MouthVsNasal/PairedT_tfce_corrp_tstat1.nii.gz

pfile=$1

# p value threshold
#threshold=0.95
threshold=$2

# figure folder
stats_fold=$3

# subcortical regions index 2-10 13-20
# cortical regions index 0-47
# To check the labels and its corresponding index
# cp ${FSLDIR}/data/atlases/HarvardOxford-Cortical.xml ~/Desktop/Tmp/
# the corresponding subcortical atlas will be set automatically
cort_atlas=${FSLDIR}/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr50-2mm.nii.gz

# create sphere ROIs centered at each voxel, radius is in mm
roi_radius=6

## End of parameter setup
################################################################


# No error should happen if fsl has been installed
command -v fsl 1>/dev/null
if [[ $? -eq 1 ]]; then
    printf '\e[1;31m %s\n \e[m' "ERROR: FSL must be installed and configured properly."
    exit 1
fi

tmp=$( remove_ext ${cort_atlas})
[ -f ${tmp}.nii.gz ] || { printf "${cort_atlas} was not found.\n"; exit 1; }

tmp=$( remove_ext ${pfile})
[ -f ${tmp}.nii.gz ] || { printf "${pfile} was not found.\n"; exit 1; }

statsdir=$( dirname ${pfile} )
statsdir=${statsdir}/${stats_fold}
[ -d ${statsdir} ] || mkdir ${statsdir}
statsfile=${statsdir}/SignificantVoxels.txt

subcort_atlas=$( echo $cort_atlas | sed -e "s/-cort-/-sub-/g")
# standard MNI brain, same resolution version to the atlas
if [[ `echo ${cort_atlas} | grep -c "1mm" ` -gt 0 ]]; then
    stdbrain=${FSLDIR}/data/standard/MNI152_T1_1mm_brain
    
elif [[ `echo ${cort_atlas} | grep -c "2mm" ` -gt 0 ]]; then
    stdbrain=${FSLDIR}/data/standard/MNI152_T1_2mm_brain
    
else
    printf '\e[1;31m %s\n \e[m' "ERROR: Unknown atlas."
    exit 1
fi

## remove overlapping subcortical labels from cortical label
printf "Remove subcortical overlaps from cortical label\n"
fslmaths ${cort_atlas} -mul 0 ${statsdir}/tmp_subcort
cnt=0
while read line
do
    if [[ `echo ${line} | grep -c "index" ` -gt 0 ]]; then
        label=$( echo $line | cut -d \> -f2 )
        label=${label%</label}
        index=$( echo $line | cut -d " "  -f2 )
        index=$( echo "$index" | sed -e "s/^index=//" )
        # echo $index $label
        if [[ ( ${cnt} -ge 2  && ${cnt} -le 10 ) || (${cnt} -ge 13  &&  ${cnt} -le 20) ]]
        then
            val=`expr ${cnt} + 1`
            fslmaths ${subcort_atlas} -thr ${val} -uthr ${val} -bin -add ${statsdir}/tmp_subcort ${statsdir}/tmp_subcort
        fi
        
        cnt=`expr ${cnt} + 1`
    fi
done < ${FSLDIR}/data/atlases/HarvardOxford-Subcortical.xml
fslmaths ${cort_atlas} -mul 0 -add 1 -sub ${statsdir}/tmp_subcort -mul ${cort_atlas} ${statsdir}/tmp_cort


info=( $( fslinfo ${cort_atlas} ))
nb_left_slice=`echo ${info[5]} / 2 | bc` 
printf "Stats Image: %s\n Threshold: %s\nCortical Atlas: %s\nSubcortical Atlas: %s\n\n" ${pfile} ${threshold} ${cort_atlas} ${subcort_atlas} > ${statsfile}
printf "Label; Max_Img(X); Max_Img(Y); Max_Img(Z); Max_STD(X); Max_STD(Y); Max_STD(Z); Number of voxels; Volume\n" >> ${statsfile}
printf "Label; Max_Img(X); Max_Img(Y); Max_Img(Z); Max_STD(X); Max_STD(Y); Max_STD(Z); Number of voxels; Volume\n"

printf "Thresholding ${pfile} at ${threshold}\n"
fslmaths ${pfile} -thr ${threshold} ${statsdir}/tmp_pfile
fslmaths ${statsdir}/tmp_cort -bin -mul ${statsdir}/tmp_pfile  ${statsdir}/tmp_mul
nbvoxs=( $( fslstats  ${statsdir}/tmp_mul -V))
if [[ ${nbvoxs} -gt 0 ]]; then
    #statements
    printf "Significant cortical voxels \n"
    
    raw_clustdir=${statsdir}/RawClust
    [ -d ${raw_clustdir} ] || mkdir ${raw_clustdir}
    peak_maskdir=${statsdir}/PeakVox
    [ -d ${peak_maskdir} ] || mkdir ${peak_maskdir}
    peak_spheredir=${statsdir}/PeakSphere
    [ -d ${peak_spheredir} ] || mkdir ${peak_spheredir}
    
    # cortical labels
    cnt=0
    while read line
    do
        if [[ `echo ${line} | grep -c "index" ` -gt 0 ]]; then
            label=$( echo $line | cut -d \> -f2 )
            label=${label%</label}
            index=$( echo $line | cut -d " "  -f2 )
            index=$( echo "$index" | sed -e "s/^index=//" )
            # echo $index $label
            if [ ${cnt} -ge 0 ] && [ ${cnt} -le 47 ]
            then
                # within-cluster voxels
                val=`expr ${cnt} + 1`
                fslmaths ${statsdir}/tmp_cort -thr ${val} -uthr ${val} -mul ${statsdir}/tmp_pfile -bin ${statsdir}/tmp_cluster
            
                # left and right
                fslmaths ${statsdir}/tmp_cluster -roi 0 ${nb_left_slice} -1 -1 -1 -1 0 1 ${statsdir}/tmp_cluster_right
                fslmaths ${statsdir}/tmp_cluster -roi ${nb_left_slice} ${nb_left_slice} -1 -1 -1 -1 0 1 ${statsdir}/tmp_cluster_left
            
                nbvoxs_left=( $( fslstats ${statsdir}/tmp_cluster_left -V))
                if [[ ${nbvoxs_left[0]} -gt 0 ]]; then
                    # echo "Left_$label" ${nbvoxs_left[0]}
                    stats=($( fslstats ${statsdir}/tmp_cluster_left -x -V ))
                    mni_coord=($( printf "${stats[0]} ${stats[1]} ${stats[2]}" | img2stdcoord -img ${stdbrain} -std ${stdbrain} - ))                
                    printf "    Left ${label}; ${stats[0]}; ${stats[1]}; ${stats[2]}; ${mni_coord[0]}; ${mni_coord[1]}; ${mni_coord[2]}; ${stats[3]}; ${stats[4]}\n"
                    printf "Left ${label}; ${stats[0]}; ${stats[1]}; ${stats[2]}; ${mni_coord[0]}; ${mni_coord[1]}; ${mni_coord[2]}; ${stats[3]}; ${stats[4]}\n" >> ${statsfile}

                    # save cluster mask
                    label_name=$(echo $label | sed -e "s/[ ,]/_/g")
                    fslmaths ${statsdir}/tmp_cluster_left ${raw_clustdir}/Left_${label_name}
                    fslmaths ${cort_atlas} -mul 0 -add 1 -roi ${stats[0]} 1 ${stats[1]} 1 ${stats[2]} 1 0 1 ${peak_maskdir}/Left_${label_name}_peak
                    fslmaths ${peak_maskdir}/Left_${label_name}_peak -kernel sphere ${roi_radius} -fmean -bin ${peak_spheredir}/Left_${label_name}_peak${roi_radius}mm
                fi
            
                nbvoxs_right=( $( fslstats ${statsdir}/tmp_cluster_right -V))
                if [[ ${nbvoxs_right[0]} -gt 0 ]]; then
                    # echo "Right_$label" ${nbvoxs_right[0]}
                    stats=($( fslstats ${statsdir}/tmp_cluster_right -x -V ))
                    mni_coord=($( printf "${stats[0]} ${stats[1]} ${stats[2]}" | img2stdcoord -img ${stdbrain} -std ${stdbrain} - ))  
                
                    printf "    Right ${label}; ${stats[0]}; ${stats[1]}; ${stats[2]}; ${mni_coord[0]}; ${mni_coord[1]}; ${mni_coord[2]}; ${stats[3]}; ${stats[4]}\n"                
                    printf "Right ${label}; ${stats[0]}; ${stats[1]}; ${stats[2]}; ${mni_coord[0]}; ${mni_coord[1]}; ${mni_coord[2]}; ${stats[3]}; ${stats[4]}\n"  >> ${statsfile}
                
                    label_name=$(echo $label | sed -e "s/[ ,]/_/g")
                    fslmaths ${statsdir}/tmp_cluster_right ${raw_clustdir}/Right_${label_name}
                    fslmaths ${cort_atlas} -mul 0 -add 1 -roi ${stats[0]} 1 ${stats[1]} 1 ${stats[2]} 1 0 1 ${peak_maskdir}/Right_${label_name}_peak
                    fslmaths ${peak_maskdir}/Right_${label_name}_peak -kernel sphere ${roi_radius} -fmean -bin ${peak_spheredir}/Right_${label_name}_peak${roi_radius}mm
                fi
            
                # stats
                #fslstats ${statsdir}/tmp_pfile -k ${statsdir}/tmp_cluster -r
            fi
            cnt=`expr ${cnt} + 1`
        
        fi
    done < ${FSLDIR}/data/atlases/HarvardOxford-Cortical.xml


else
    printf "No significant cortical voxel found\n."
fi

# subcortical labels
printf "\n"
fslmaths ${statsdir}/tmp_subcort -bin -mul ${statsdir}/tmp_pfile ${statsdir}/tmp_mul
nbvoxs=( $( fslstats ${statsdir}/tmp_mul -V))
if [[ ${nbvoxs} -gt 0 ]]; then
    printf "Significant subcortical voxels \n"
    cnt=0
    
    raw_clustdir=${statsdir}/RawClust
    [ -d ${raw_clustdir} ] || mkdir ${raw_clustdir}
    peak_maskdir=${statsdir}/PeakVox
    [ -d ${peak_maskdir} ] || mkdir ${peak_maskdir}
    peak_spheredir=${statsdir}/PeakSphere
    [ -d ${peak_spheredir} ] || mkdir ${peak_spheredir}
    
    while read line
    do
        if [[ `echo ${line} | grep -c "index" ` -gt 0 ]]; then
            label=$( echo $line | cut -d \> -f2 )
            label=${label%</label}
            index=$( echo $line | cut -d " "  -f2 )
            index=$( echo "$index" | sed -e "s/^index=//" )
            # echo $index $label
            if [[ ( ${cnt} -ge 2  && ${cnt} -le 10 ) || (${cnt} -ge 13  &&  ${cnt} -le 20) ]]
            then
                val=`expr ${cnt} + 1`
                fslmaths ${subcort_atlas} -thr ${val} -uthr ${val} -mul ${statsdir}/tmp_pfile -bin ${statsdir}/tmp_cluster
                nbvoxs=( $( fslstats ${statsdir}/tmp_cluster -V))
                if [[ ${nbvoxs[0]} -gt 0 ]]; then
                    stats=($( fslstats ${statsdir}/tmp_cluster -x -V))
                    mni_coord=($( printf "${stats[0]} ${stats[1]} ${stats[2]}" | img2stdcoord -img ${stdbrain} -std ${stdbrain} - ))   
                    printf "    ${label}; ${stats[0]}; ${stats[1]}; ${stats[2]}; ${mni_coord[0]}, ${mni_coord[1]}; ${mni_coord[2]}; ${stats[3]}; ${stats[4]}\n"  
                    printf "${label}; ${stats[0]}; ${stats[1]}; ${stats[2]}; ${mni_coord[0]}; ${mni_coord[1]}; ${mni_coord[2]}; ${stats[3]}; ${stats[4]}\n"   >> ${statsfile}
                
                    label_name=$(echo $label | sed -e "s/[ ,]/_/g")
                    fslmaths ${statsdir}/tmp_cluster ${raw_clustdir}/${label_name}
                    fslmaths ${cort_atlas} -mul 0 -add 1 -roi ${stats[0]} 1 ${stats[1]} 1 ${stats[2]} 1 0 1 ${peak_maskdir}/${label_name}_peak
                    fslmaths ${peak_maskdir}/${label_name}_peak -kernel sphere ${roi_radius} -fmean  -bin ${peak_spheredir}/${label_name}_peak${roi_radius}mm
                fi
            fi
        
            cnt=`expr ${cnt} + 1`
        fi
    done < ${FSLDIR}/data/atlases/HarvardOxford-Subcortical.xml

else
    printf "No significant subcortical voxel found.\n"
fi

# clear temporal data
ls ${statsdir}/tmp* 1>/dev/null 2>/dev/null
if [[ $? -eq 0 ]]; then
    printf "Clearing temporary files\n"
    rm ${statsdir}/tmp*
fi
