%% Figure 1B, plot the primary olfactory cortex ROI

% full-path to roi file
roi = '~/Documents/RestfMRI/SourceData/Figure_1_sourcedata/Figure_1B_MNI152_ROI_1mm.nii.gz';

% directory to save figures
figdir = '~/Documents/RestfMRI/Figures';


if ~exist( figdir, 'dir')
    mkdir( figdir);
end

% coronal slice
figure;
SetPrintProp( gcf, 0.3, 0.3);
subplot( 2, 4, 1);
MNIOverlay( roi, [93, 0, 0],...
    'brain_threshold', [3, 8]*1000,...
    'cm', [1 0 0],...
    'hold', 'on',...
    'refline', [125:3:141],...
    'refline_view', 'y',...
    'refline_style', '-',...
    'refline_width', 0.5,...
    'refline_color', 'k');

cnt = 3;
for s = [125:3:141]
    subplot( 2, 4, cnt);
    cnt = cnt + 1;
    MNIOverlay( roi, [0, s, 0],...
        'brain_threshold', [3, 8]*1000,...
        'cm', [1 0 0],...
        'hold', 'on');
end

print( gcf, '-dpdf', '-fillpage', fullfile( figdir, 'ROI_Sagittal.pdf'));


% axial slice
axial_slice = [52 : 2 : 62];
figure;
SetPrintProp( gcf, 0.3, 0.3);
subplot( 2, 4, 1);
MNIOverlay( roi, [93, 0, 0],...
    'brain_threshold', [3, 8]*1000,...
    'cm', [1 0 0],...
    'hold', 'on',...
    'refline', axial_slice,...
    'refline_view', 'z',...
    'refline_style', '-',...
    'refline_width', 0.5,...
    'refline_color', 'k');

cnt = 2;
for s = axial_slice
    subplot( 2, 4, cnt);
    cnt = cnt + 1;
    MNIOverlay( roi, [0, 0, s],...
        'brain_threshold', [3, 8]*1000,...
        'cm', [1 0 0],...
        'hold', 'on');
end

print( gcf, '-dpdf', '-fillpage', fullfile( figdir, 'ROI_Axial.pdf'));



%% Figure 2, plot parcellation

% plot correlation
roi = '~/Documents/RestfMRI/SourceData/Figure_2_sourcedata1/Figure_2A_1mm.nii.gz';
figdir = '~/Documents/RestfMRI/Figures';

stdbrain = '/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz';
figure;
SetPrintProp( gcf, 0.3, 0.3);
subplot( 2, 4, 1);
MNIOverlay( roi, [93, 0, 0],...
    'brain_threshold', [3, 8]*1000,...
    'threshold', [eps, 0.4],...
    'hold', 'on',...
    'refline', [125:3:141],...
    'refline_view', 'y',...
    'refline_style', '-',...
    'refline_width', 0.5,...
    'refline_color', 'k',...
    'brainfile', stdbrain);

cnt = 3;
for s = [125:3:141]
    subplot( 2, 4, cnt);
    cnt = cnt + 1;
    MNIOverlay( roi, [0, s, 0],...
        'brain_threshold', [3, 8]*1000,...
        'threshold', [eps, 0.4],...
        'hold', 'on',...
        'brainfile', stdbrain,...
        'show_colorbar', 'yes');
end

print( gcf, '-dpdf', '-fillpage', fullfile( figdir, 'matplot_stabilityMap.pdf'));
close
print( gcf, '-dpdf', '-fillpage', fullfile( figdir, 'matplot_stabilityMap_cb.pdf'));


% plot correlation histogram
clear;clc
cd('~/Documents/RestfMRI/SourceData/Figure_2_sourcedata1');
figdir = '~/Documents/RestfMRI/Figures';

h1 = load( 'Figure_2A_correlation_LeftHemisphere.mat'); 
h2 = load( 'Figure_2A_correlation_RightHemisphere.mat');
v = cat( 2, h1.leave_one_out_rval, h2.leave_one_out_rval); 
figure; 
SetPrintProp( gcf, 0.05, 0.1);
hist( v, 20);
set( gca, 'Box', 'off',...
    'TickLength', [1, 1]*0.03, ...
    'TickDir', 'out',...
    'XTick', 0:0.1:0.4,...
    'XTickLabel', {'0', '', '0.2', '', '0.4'},...
    'ytick', 0:20:60,...
    'ylim', [0, 60]);

print( gcf, '-dpdf', '-fillpage', fullfile( figdir, 'matplot_stabilityMap_histogram.pdf'));

% mean +/- std
[mean(v), std(v)/sqrt( length(v))]

% 95% confidence interval
SEM = std(v) / sqrt( length( v));   
ts = tinv([0.025  0.975],length(v)-1); 
CI = mean(v) + ts*SEM


% plot parcellation result
img2plot = '~/Documents/RestfMRI/SourceData/Figure_2_sourcedata2/Figure_2_Parcellation_1mm.nii.gz';
figdir = '~/Documents/RestfMRI/Figures';

% coronal slice
figure;
SetPrintProp( gcf, 0.3, 0.3);
subplot( 2, 4, 1);
MNIOverlay( img2plot, [93, 0, 0],...
    'brain_threshold', [3, 8]*1000,...
    'cm', [1 0 0; 0 1 0; 0 0 1; [204, 51, 255]/255 ],...
    'hold', 'on',...
    'refline', [125:3:141],...
    'refline_view', 'y',...
    'refline_style', '-',...
    'refline_width', 0.5,...
    'refline_color', 'k');

cnt = 3;
for s = [125:3:141]
    subplot( 2, 4, cnt);
    cnt = cnt + 1;
    MNIOverlay( img2plot, [0, s, 0],...
        'brain_threshold', [3, 8]*1000,...
        'cm', [1 0 0; 0 1 0; 0 0 1; [204, 51, 255]/255 ],...
        'hold', 'on');
end

print( gcf, '-dpdf', '-fillpage', fullfile( figdir, 'Parcellation_Sagittal.pdf'));

% axial slice
figure;
SetPrintProp( gcf, 0.3, 0.3);
cnt = 3;
for s = [56 60]
    subplot( 2, 4, cnt);
    cnt = cnt + 1;
    MNIOverlay( img2plot, [0, 0, s],...
        'brain_threshold', [3, 8]*1000,...
        'cm', [1 0 0; 0 1 0; 0 0 1; [204, 51, 255]/255 ],...
        'hold', 'on');
end

print( gcf, '-dpdf', '-fillpage', fullfile( figdir, 'Parcellation_Axial.pdf'));


% overlap bettween parcellation and anatomical sub-region
% see Overlap_Manual_Parc.m


%% Figure 3. parcellation of left and right combined ROI
clear;
clc;
close all

figdir = '~/Documents/RestfMRI/Figures';

figure;
SetPrintProp( gcf, 0.2, 0.6);

% colormap
cm = [0 0 1;...
    0 1 0;...
    1 1 0];

subplot( 421);
roi = '~/Documents/RestfMRI/SourceData/Figure_3_sourcedata/Figure_3_K3_MNI152_ROI_1mm.nii.gz';
MNIOverlay( roi, [0, 137, 0],...
    'brain_threshold', [3, 8]*1000,...
    'cm', cm,...
    'hold', 'on');

subplot( 422);
MNIOverlay( roi, [0, 129, 0],...
    'brain_threshold', [3, 8]*1000,...
    'cm', cm,...
    'hold', 'on');


roi = '~/Documents/RestfMRI/SourceData/Figure_3_sourcedata/Figure_3_K4_MNI152_ROI_1mm.nii.gz';
cm = [0 0 1;...
    1 0 0;...
    0 1 0;...
    [204, 51, 255]/255 ];

subplot( 423);
MNIOverlay( roi, [0, 137, 0],...
    'brain_threshold', [3, 8]*1000,...
    'cm', cm,...
    'hold', 'on');
subplot( 424);
MNIOverlay( roi, [0, 129, 0],...
    'brain_threshold', [3, 8]*1000,...
    'cm', cm,...
    'hold', 'on');


roi = '~/Documents/RestfMRI/SourceData/Figure_3_sourcedata/Figure_3_K5_MNI152_ROI_1mm.nii.gz';
cm = [1 1 0;...    
    0 0 1;...
    [204, 51, 255]/255;...
    0 1 0;...  
    1 0 0];

subplot( 425);
MNIOverlay( roi, [0, 137, 0],...
    'brain_threshold', [3, 8]*1000,...
    'cm', cm,...
    'hold', 'on');

subplot( 426);
MNIOverlay( roi, [0, 129, 0],...
    'brain_threshold', [3, 8]*1000,...
    'cm', cm,...
    'hold', 'on');


roi = '~/Documents/RestfMRI/SourceData/Figure_3_sourcedata/Figure_3_K6_MNI152_ROI_1mm.nii.gz';
cm = [1 0 0;...    
    0 0 1;...
    0 1 0;...
    0.6 0.3 0;...
    [204, 51, 255]/255;...     
    1 0 1];

subplot( 427);
MNIOverlay( roi, [0, 137, 0],...
    'brain_threshold', [3, 8]*1000,...
    'cm', cm,...
    'hold', 'on');

subplot( 428);
MNIOverlay( roi, [0, 129, 0],...
    'brain_threshold', [3, 8]*1000,...
    'cm', cm,...
    'hold', 'on');


print( gcf, '-dpdf', '-fillpage', fullfile( figdir, 'Manual_Vs_Parcellation_LRMerged_C3-6.pdf'));


% overlap bettween parcellation and anatomical sub-region
% see Overlap_Manual_Parc.m


%% Figure 4 & 5. sub-region unique functional connectivity and common connection
% Plot_ConnectivityOnSurface_Tstat.m

%% Figure 6.
% Plot_Figure6.m


%% Suppl. Figure 1. supplementary parcellation analysis

roi = '~/Documents/RestfMRI/SourceData/Supplementary_Figure_1_sourcedata/Suppl_Figure1_Parcellation_v2_1mm.nii.gz';

label_color = [1 0 0;...
    [204, 51, 255]/255;...
    0 0 1;...
    0 1 0 ];
figure;
SetPrintProp( gcf, 0.3, 0.3);
subplot( 2, 4, 1);
MNIOverlay( roi, [93, 0, 0], 'brain_threshold', [3, 8]*1000, 'cm', label_color, 'hold', 'on',...
    'refline', [125:3:141], 'refline_view', 'y', 'refline_style', '-', 'refline_width', 1,'refline_color', 'k');

cnt = 3;
for s = [125:3:141]
    subplot( 2, 4, cnt);
    cnt = cnt + 1;
    MNIOverlay( roi, [0, s, 0], 'brain_threshold', [3, 8]*1000, 'cm', label_color, 'hold', 'on');
end
print( gcf, '-dpdf', '-fillpage', '~/Documents/RestfMRI/Figures/Parcellation_V2_Sagittal.pdf');

% axial slice
figure;
SetPrintProp( gcf, 0.3, 0.3);
cnt = 3;
for s = [56 60]
    subplot( 2, 4, cnt);
    cnt = cnt + 1;
    MNIOverlay( roi, [0, 0, s], 'brain_threshold', [3, 8]*1000, 'cm', label_color, 'hold', 'on');
end
print( gcf, '-dpdf', '-fillpage', '~/Documents/RestfMRI/Figures/Parcellation_V2_Axial.pdf');




roi = '~/Documents/RestfMRI/SourceData/Supplementary_Figure_1_sourcedata/SupplFigure1_TN_Parcellation_1mm.nii.gz';

label_color = [1 0 0;...
    [204, 51, 255]/255;...
    0 0 1;...
    0 1 0 ];
figure;
SetPrintProp( gcf, 0.3, 0.3);
subplot( 2, 4, 1);
MNIOverlay( roi, [93, 0, 0], 'brain_threshold', [3, 8]*1000, 'cm', label_color, 'hold', 'on',...
    'refline', [125:3:141], 'refline_view', 'y', 'refline_style', '-', 'refline_width', 1,'refline_color', 'k');

cnt = 3;
for s = [125:3:141]
    subplot( 2, 4, cnt);
    cnt = cnt + 1;
    MNIOverlay( roi, [0, s, 0], 'brain_threshold', [3, 8]*1000, 'cm', label_color, 'hold', 'on');
end
print( gcf, '-dpdf', '-fillpage', '~/Documents/RestfMRI/Figures/TN_Parcellation_Sagittal.pdf');

% axial slice
figure;
SetPrintProp( gcf, 0.3, 0.3);
cnt = 3;
for s = [56 60]
    subplot( 2, 4, cnt);
    cnt = cnt + 1;
    MNIOverlay( roi, [0, 0, s], 'brain_threshold', [3, 8]*1000, 'cm', label_color, 'hold', 'on');
end
print( gcf, '-dpdf', '-fillpage', '~/Documents/RestfMRI/Figures/TN_Parcellation_Axial.pdf');


% overlap bettween parcellation and anatomical sub-region
% see Overlap_Manual_Parc.m


%% Suppl. Figure 2. functional connecivity laterality
wkpath = '~/Documents/RestfMRI/SourceData/Supplementary_Figure_2_sourcedata';

rois = {'Aon', 'Tub', 'PirF', 'PirT'};
stdbrain = '/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz';

threshold = [0.9989, 1];
r = autumn;
b = winter;

figure;
SetPrintProp( gcf, 0.5, 0.3);
for roi_idx = 1 : length( rois)
    roi = fullfile( wkpath, rois{ roi_idx}, 'pos_s3_OneSampT_tfce_corrp_tstat1_1mm.nii.gz');
        
    subplot( 2, 5, roi_idx);
    if roi_idx == 2
        s = 61;
    else
        s = 57;
    end
    MNIOverlay( roi, s,...
        'coord_style', 'z',...
        'brain_threshold', [3, 8]*1000,...
        'cm', r,...
        'threshold', threshold,...
        'hold', 'on',...
        'brainfile', stdbrain,...
        'fontsize', 6,...
        'show_colorbar', 'yes');

    roi = fullfile( wkpath, rois{ roi_idx}, 'neg_s3_OneSampT_tfce_corrp_tstat1_1mm.nii.gz');
    subplot( 2, 5, roi_idx+5);
    if roi_idx == 2
        s = 61;
    else
        s = 57;
    end
    MNIOverlay( roi, s,...
        'coord_style', 'z',...
        'brain_threshold', [3, 8]*1000,...
        'cm', b,...
        'threshold', threshold,...
        'hold', 'on',...
        'brainfile', stdbrain, ...
        'fontsize', 6, ...
        'show_colorbar', 'yes');
end


print( gcf, '-dpdf', '-fillpage', '~/Documents/RestfMRI/Figures/lateralization_brain2.pdf');
close
print( gcf, '-dpdf', '-fillpage', '~/Documents/RestfMRI/Figures/lateralization_cb1.pdf');
close
print( gcf, '-dpdf', '-fillpage', '~/Documents/RestfMRI/Figures/lateralization_cb2.pdf');


