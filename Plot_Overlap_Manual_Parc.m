% overlap between manual ROIs and automatic parcellations

%% parameter setup for the permutation test
% Figure 2B
parfile = '~/Documents/RestfMRI/SourceData/Figure_2_sourcedata2/Figure_2_Parcellation_2mm.nii.gz';

roi = '~/Documents/RestfMRI/SourceData/Figure_2_sourcedata2/Figure_2_MNI152_ROI_ManualSegmentation_2mm.nii.gz';
roi_label = {'Aon', 'Tub', 'PirF', 'PirT'};
roi_label = fliplr( roi_label);

% 1-LR 2-L 3-R
hemi = 2:3;
% 'first' | 'last'
ref_label = 'last';
% 'normal' | 'reverse'
ydir = 'normal';


% Figure 3, K=3-6
roi = '~/Documents/RestfMRI/SourceData/Figure_3_source_data/Figure_3_MNI152_ROI_ManualSegmentation_2mm_LRMerged.nii.gz';
parfile = '~/Documents/RestfMRI/SourceData/Figure_3_source_data/Figure_3_K6_MNI152_ROI_2mm.nii.gz';
roi_label = {'LAon', 'RAon', 'LTub', 'TRub', 'LPirF', 'RPirF', 'LPirT', 'RPirT'};

hemi = 1;
% 'first' | 'last'
ref_label = 'first';
% 'normal' | 'reverse'
ydir = 'reverse';


% Supplementary Figure 1B 
roi = '~/Documents/RestfMRI/SourceData/Supplementary_Figure_1_sourcedata/Suppl_Figure1_Manual_Segmentation_v2_2mm.nii.gz';
roi_label = {'Aon', 'Tub', 'PirF', 'PirT'};
parfile = '~/Documents/RestfMRI/SourceData/Supplementary_Figure_1_sourcedata/Suppl_Figure1_Parcellation_v2_2mm.nii.gz';
hemi = 2:3;
% 'first' | 'last'
ref_label = 'first';
% 'normal' | 'reverse'
ydir = 'reverse';


% Supplementary Figure 1D 
roi = '~/Documents/RestfMRI/SourceData/Figure_2_sourcedata2/Figure_2_MNI152_ROI_ManualSegmentation_2mm.nii.gz';
roi_label = {'Aon', 'Tub', 'PirF', 'PirT'};
roi_label = fliplr( roi_label);
parfile = '~/Documents/RestfMRI/SourceData/Supplementary_Figure_1_sourcedata/SupplFigure1_TN_Parcellation_2mm.nii.gz';
hemi = 2:3;
% 'first' | 'last'
ref_label = 'last';
% 'normal' | 'reverse'
ydir = 'normal';




%% code to do the permutation test
figdir = '~/Documents/RestfMRI/Figures';

% permutation test
nbpermuts = 10000;

% % 'first' | 'last'
% ref_label = 'first'; 
% % 'normal' | 'reverse'
% ydir = 'reverse';

roi_vol = MRIread( roi);
roi_vol.vol = permute( roi_vol.vol, [2, 1, 3]);
par_vol = MRIread( parfile);
par_vol.vol = permute( par_vol.vol, [2, 1, 3]);

% 1:N, left; N:2*N, right.
N = floor( size( roi_vol.vol, 1) / 2);

roi_vals = unique( roi_vol.vol(:));
roi_vals( roi_vals == 0) = [];
nb_roi_vals = length( roi_vals);

par_vals = unique( par_vol.vol(:));
par_vals( par_vals == 0 ) = [];
nb_par_vals = length( par_vals);


% switch parcellation orders, hopefully this works
% antomical x parcellation, roi x par x hemi (all, left, right)
nb_overlaps = zeros( nb_roi_vals, nb_par_vals, 3);
reorder_nb_overlap = zeros( nb_roi_vals, nb_par_vals, 3);

hemi_name = {'LR', 'Left', 'Right'};

% this automic sorting method only works when the number of anatomical ROIs
% is smaller than the number of parcellations
figure;
SetPrintProp( gcf, 0.4, 0.6);
SetPrintProp( gcf, 0.2, 0.4)
sort_par_vals = cell( 3, 1);
for hemi_idx = hemi
    tmp_roi_vol = roi_vol.vol;
    tmp_par_vol = par_vol.vol;    
    if hemi_idx == 2
        tmp_roi_vol( 1:N, :, :) = 0;
        tmp_par_vol( 1:N, :, :) = 0;
    elseif hemi_idx == 3
        tmp_roi_vol( N+1:end, :, :) = 0;
        tmp_par_vol( N+1:end, :, :) = 0;
    end
    
    if sum( tmp_roi_vol(:) ~= 0) ~= sum( tmp_par_vol(:) ~=0)
        error( 'The anatomical template and parcellation results are not matched.');
    end
    
    for par_idx = 1 : nb_par_vals 
        cur_par_val = par_vals( par_idx);
        nb_par_voxs = sum( tmp_par_vol(:) == cur_par_val);
        for roi_idx = 1 : nb_roi_vals
            cur_roi_val = roi_vals( roi_idx);
            tmp = tmp_par_vol == cur_par_val & tmp_roi_vol == cur_roi_val;
            nb_overlaps( roi_idx, par_idx, hemi_idx) = sum( tmp(:)) / nb_par_voxs;
        end
    end
  
    switch lower( ref_label)
        case 'last'
            ref_row = nb_roi_vals;
        case 'first'
            ref_row = 1;
        otherwise
            error('Invalid reference label.');
    end
    
    [mval, midx] = max( nb_overlaps( :, :, hemi_idx), [], 1);
    [~, start_col] = min( abs( ref_row - midx));

    % old -> new column index
    last_row = [];
    last_row( 1, 1:2) = [ref_row, 1];
    last_col = [];
    last_col( 1, 1:2) = [start_col, 1];

    for k = 2 : nb_roi_vals
        tmp = nb_overlaps( :, :, hemi_idx);
        tmp( :, last_col( :, 1)) = -1;
        [tmp_mval, tmp_midx] = max( tmp, [], 1);        
        % row next to the last one
        tmp_midx_diff = abs( last_row( k-1, 1) - tmp_midx);
        tmp_midx_diff( last_col( :, 1)) = nb_par_vals * 2;
        [min_val, min_ind] = min( tmp_midx_diff);
        
        last_row( k, 1) = tmp_midx( min_ind);
        last_row( k, 2) = k;

        if k <= nb_par_vals
            last_col( k, 1) = min_ind;
            last_col( k, 2) = k;
        end
    end

    % if nb_roi_vals > nb_par_vals
    if nb_roi_vals < nb_par_vals        
        v = (1 : nb_par_vals)';
        v( ismember( v, last_col( :, 1))) = [];        
        loc = [];
        loc( :, 1) = v;
        loc( :, 2) = size( last_col, 1) + 1 : length( tmp);
        last_col = cat( 1, last_col, loc);
    end
    
    sort_par_vals{ hemi_idx} = last_col( :, 1);
    
    reorder_nb_overlap( :, :, hemi_idx) = nb_overlaps( :, last_col( :, 1), hemi_idx);

    xticklabel = num2cell( 1:nb_par_vals);
    xticklabel = cellfun( @(x) ['C', num2str(x)], xticklabel, 'UniformOutput', false);
    
    subplot( 3, 2, 2*(hemi_idx-1) + 1);
    imagesc( nb_overlaps( :, :, hemi_idx));
    title( ['Unsorted ', hemi_name{ hemi_idx}]);
    set( gca, 'ytick', 1:nb_roi_vals,...
        'yticklabels', roi_label,...
        'xtick', 1:nb_par_vals,...
        'xticklabel', xticklabel,...
        'ydir', ydir, 'ticklength', [0, 0]);
    colorbar;
    caxis( [0, 1]);
    
    subplot( 3, 2, 2*(hemi_idx-1) + 2);
    imagesc( 100* reorder_nb_overlap( :, :, hemi_idx));
    title( ['Sorted ', hemi_name{ hemi_idx}]);
    set( gca, 'ytick', 1:nb_roi_vals,...
        'yticklabels', roi_label,...
        'xtick', 1:nb_par_vals,...
        'xticklabel', xticklabel, 'ydir', ydir, 'ticklength', [0, 0]);
    ylabel( 'Manual ROI'); xlabel( 'Parcellation');
    
    c = colorbar;
    caxis( [0, 100]);
    set( c, 'ticks', 0:25:100);
    xlabel( c, 'Percentage (%)');
end

[rb, r,b] = CustomColormap;
colormap( r);


permut_nb_overlaps = zeros( nb_roi_vals, nb_par_vals, 3, nbpermuts+1);
zscore = nan( nb_roi_vals, nb_par_vals, 3);
figure;
SetPrintProp( gcf, 0.3, 0.1);
SetPrintProp( gcf, 0.2, 0.4)
for hemi_idx = hemi
    tmp_roi_vol = roi_vol.vol;
    tmp_par_vol = par_vol.vol;    
    if hemi_idx == 2
        tmp_roi_vol( 1:N, :, :) = 0;
        tmp_par_vol( 1:N, :, :) = 0;
    elseif hemi_idx == 3
        tmp_roi_vol( N+1:end, :, :) = 0;
        tmp_par_vol( N+1:end, :, :) = 0;
    end
    
    for pidx = 1 : nbpermuts + 1
        fprintf( 'Run %d/%d, Permutation: %d/%d           \r', hemi_idx, 3, pidx, nbpermuts+1);
        if pidx > 1
            tmp_loc = find( tmp_roi_vol(:));
            tmp_val = tmp_roi_vol( tmp_loc);
            tmp_roi_vol( tmp_loc) = tmp_val( randperm( length( tmp_loc)));
        end        
        
        for par_idx = 1 : nb_par_vals
            cur_par_val = sort_par_vals{ hemi_idx}( par_idx);
            nb_par_voxs = sum( tmp_par_vol(:) == cur_par_val);
            for roi_idx = 1 : nb_roi_vals
                cur_roi_val = roi_vals( roi_idx);
                tmp = tmp_par_vol == cur_par_val & tmp_roi_vol == cur_roi_val;
                permut_nb_overlaps( roi_idx, par_idx, hemi_idx, pidx) = sum( tmp(:)) / nb_par_voxs;
            end
        end
    end    
    fprintf( '\n');
    
    [avg, sd] = normfit( permute( squeeze( permut_nb_overlaps( :, :, hemi_idx, 2:end)), [3 1 2]));
    avg = squeeze( avg);
    sd = squeeze( sd);
    zscore( :, :, hemi_idx) = ( permut_nb_overlaps( :, :, hemi_idx, 1) - avg) ./ sd;
    
    subplot( 1, 3, hemi_idx);
    v = squeeze( permut_nb_overlaps( randperm( nb_roi_vals, 1), randperm( nb_par_vals, 1), hemi_idx, 2:end));
    histfit( v, 20);
    xlabel( 'Overlap percentage');
    ylabel( 'Count');
    title( [hemi_name{ hemi_idx}, ' example hist']);
    set( gca, 'fontsize', 10);
end


% FDR on left and right 
tmp_z = zscore( :, :, hemi);
tmp_pval = G_z2p( tmp_z);
fdrp = fdr( tmp_pval(:), 0.001);
    
figure;
SetPrintProp( gcf, 0.4, 0.6);
SetPrintProp( gcf, 0.2, 0.4)
for hemi_idx = hemi
    subplot( 3, 2, 2*(hemi_idx-1) + 1);
    imagesc( 100* permut_nb_overlaps( :, :, hemi_idx, 1));
    title( ['Sorted ', hemi_name{ hemi_idx}]);
    set( gca, 'ytick', 1:nb_roi_vals, 'yticklabels', roi_label, 'xtick', 1:nb_par_vals, 'xticklabel', xticklabel, 'ydir', ydir,...
        'ticklength', [0, 0], 'fontsize', 6);
    ylabel( 'Manual ROI'); xlabel( 'Parcellation');
    xtickangle( 45);
    
    c = colorbar;
    caxis( [0, 100]);
    set( c, 'Limits', [0, 100], 'ticks', 0:25:100);
    xlabel( c, 'Percentage (%)');    
    
    % show positive only
    tmp_pval = G_z2p( zscore( :, :, hemi_idx)) + (1 - (zscore( :, :, hemi_idx) > 0));

    subplot( 3, 2, 2*(hemi_idx-1) + 2);
    G_Imagesc( zscore( :, :, hemi_idx), 'aterisks', tmp_pval, 'aterisks_level', fdrp);
    c = colorbar;
    caxis( [-10, 10]);
    set( c, 'Limits', [-10, 10], 'ticks', -10:5:10);
    set( gca, 'ytick', 1:nb_roi_vals, 'yticklabels', roi_label, 'xtick', 1:nb_par_vals, 'xticklabel', xticklabel, 'ydir', ydir,...
        'ticklength', [0, 0], 'fontsize', 6);
    ylabel( 'Manual ROI'); xlabel( 'Parcellation');
    xlabel( c, 'z score');      
    xtickangle( 45);
end


% jet colormap for z score
colormap( jet);
print( gcf, '-dpdf', '-fillpage', fullfile( figdir, [name_prefix, '_1.pdf']));

% red colormap for percentage figure
[~, r] = CustomColormap;
colormap( r);
print( gcf, '-dpdf', '-fillpage', fullfile( figdir, [name_prefix, '_2.pdf']));



