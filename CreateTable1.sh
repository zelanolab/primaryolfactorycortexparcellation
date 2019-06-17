#!/usr/bin/env bash

## Make Table 1

# Subregion specific connectivity
wkpath=~/Documents/RestfMRI/Revision_1/SourceData/Figure_4_sourcedata
FSL_VoxelStats.sh ${wkpath}/unique_pos_LR_Aon 0.999 unique_pos_LR_Aon
FSL_VoxelStats.sh ${wkpath}/unique_pos_LR_PirF 0.999 unique_pos_LR_PirF
FSL_VoxelStats.sh ${wkpath}/unique_pos_LR_PirT 0.999 unique_pos_LR_PirT
FSL_VoxelStats.sh ${wkpath}/unique_pos_LR_Tub 0.999 unique_pos_LR_Tub

# Connectivity common to all subregions
wkpath=~/Documents/RestfMRI/Revision_1/SourceData/Figure_5_sourcedata
FSL_VoxelStats.sh ${wkpath}/Figure_5_pos_Common2AllROIs 0.999 pos_Common2AllROIs


# The number of voxels within each brain region are given in SignificantVoxels.txt
