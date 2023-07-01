function [pre_pad, post_pad] = MNIOverlay( statsfile, coord2plot, varargin)
% Overlay statistical images onto MNI standard brain
% Use this function with FSL's fslview
% Note. Image index in fslview starts from 0, here in Matlab, it starts from 1.
% 
% Requires FSL and Freesurfer's matlab toolbox
% 
% Usage
%   % plot all slices with a p > 0.95, use skull-stripped MNI brain as background
%   % this would generate a bunch of figures
%   MNIOverlay( '~/tfce_corrp_tstat2.nii.gz', []); 
% 
%   % plot sagittal slice x=48 (voxel, starting from 1), the voxel index starts from 0 in fslview!
%   MNIOverlay( '~/tfce_corrp_tstat2.nii.gz', [48 0 0]); % or
%   MNIOverlay( '~/tfce_corrp_tstat2.nii.gz', [48], 'coord_style', 'x');
% 
%   % plot multiple sagittal slices
%   MNIOverlay( '~/tfce_corrp_tstat2.nii.gz', [48 60 62], 'coord_style', 'x');
% 
%   % specfiy threshold, brain image etc. as key-value input
%   MNIOverlay( '~/tfce_corrp_tstat2.nii.gz', [48 0 0], 'key', value, ...);
% 
%   % plot brain slices without stats image
%   MNIOverlay( '', [48 0 60], 'brainfile', 'MNI152_T1_1mm_brain'); 
% 
%   % plot t stats image with both positive and negative values  
%   MNIOverlay( '~/Desktop/tstat2.nii.gz', [0 55 0], 'cm', cat( 3,
%   autumn(100), pink(100)), 'threshold', [2, 3; -2  -5]);
% 
%   % plot different slices on the same figure, there can only be one non-zero number in the coordinates
%   figure;
%   subplot( 231);
%   MNIOverlay( '~/Desktop/tstat2.nii.gz',...
%       [0 55 0], 'cm', cat( 3, autumn(100), pink(100)), 'threshold', [2, 3; -2  -5], 'hold', 'on');
% 
%   subplot( 232);
%   MNIOverlay( '~/Desktop/tstat2.nii.gz',...
%       [60 0 0], 'cm', cat( 3, autumn(100), pink(100)), 'threshold', [2, 3; -2  -5], 'hold', 'on');
% 
%   subplot( 233);
%   MNIOverlay( '~/Desktop/tstat2.nii.gz',...
%       [0 0 50], 'cm', cat( 3, autumn(100), pink(100)), 'threshold', [2, 3; -2  -5], 'hold', 'on');
% 
%   % I want to add a marker at coordinate (x, y, z)
%    MNIOverlay( '~/tfce_corrp_tstat2.nii.gz', [48 60 62], 'coord_style',
%    'xyz', 'show_marker', 'yes', 'marker_type', 'dot'); % 'dot' | 'cross'
%    | 'dot+cross'
%  
%   % I am plotting only sagittal, axial or cronal slices, but I still want to add a marker
%    MNIOverlay( '~/tfce_corrp_tstat2.nii.gz', [48 60 62], 'coord_style',
%    'x', 'cross', [48 60 62], 'cross_type', 'dot+cross'); % 'dot' | 'cross' | 'dot+cross'
% 
%   see MNIOverlay_ExampleUsage.m for more examples
% 
%  Outline options
%       To outline clusters, use 
%           'outline_overlay', 'yes' | 'no'
%           'outline_only', 'yes' | 'no'
%           'outline_color', 'r' | 'b' etc. or rgb value [1 0 0]
% 
% Input
%   statsfile,  full path to stats image (including extension), [] to plot brain images only 
%   coord2plot, image coordinates for plotting (index starts from 1)
% 
%   Key-Value pairs
%       'brainfile',       '' standard brain image, e.g. 'MNI152_T1_1mm_brain.nii.gz', use full path if FSL is not installed
%       'figdir',          '', directory to save figures, default to stats image directory
%       'threshold',       [0.95, 1], threshold of stats image
%       'cm',              autumn( 200), colormap Nx3 matrix
%       'nbticks',         5, number of tick for colorbar
%       'tick_precision',  2,
%       'brain_threshold', e.g. [3000, 8000], threhsold for standard brain
%       'brain_background', 'k'|'w', background color, works only for skull-stripped brain images
%       'outline_overlay', 'yes' | 'no' (default), outline clusters
%       'outline_only', 'yes' | 'no' (default), only show outlines
%       'outline_color', rgb color (default: []), outline color
%       'rect', [],... % to make a zoomed-in version, 3x2 matrix
%       'rect_trim', 'no',...
%       'rect_color', 'k',...
%       'rect_linewidth', 1.5,...  
%       'atlas', '',... % full-path-to-nifti-image, outline atlas
%       'atlas_value', [],...
%       'atlas_linewidth', 1,...
%       'atlas_hemi', 'lr',...
%       'atlas_color', [],... % Nx3 matrix of RGB color
%       'text_color', 'b', ... % label text color
%       'fontsize', 6,...
%       'refline', [],... % plot vertical or horizontal lines to indicates
%                       the location of slices being plotted
%       'refline_view', 'x',...%'x' | 'y' | 'z'
%       'refline_style', '-',...
%       'refline_color', 'k',...
%       'refline_width', 0.5' ); 
% 
%   If the stats image has both positive and negative values and you want
%   to diplay them with different colormaps, see below
%       'threshold', [2, 3; -2, -3]
%       'cm', cat( 3, autumn( 100), winter( 100))
%   
% Bash scripts to generate white background standard brain
%   stdbrain=${FSLDIR}/data/standard/MNI152_T1_2mm_brain
%   savefile=~/Desktop/Tmp/bk_MNI512_2mm_brain
%   fslmaths ${stdbrain} -add 0.001 -uthr 0.01 -bin -mul 10000 -add ${stdbrain} ${savefile}
%  See also MNIWhiteBackground.m
%  The idea is to make background has a much higher value than the brain
% 
% naturalzhou@gmail.com
% Zelano Lab @Northwestern University
% https://sites.northwestern.edu/zelano/
% 

warning( 'Use with caution');

if nargin < 2
    error( 'Not enought input arguments.');
end

if ~exist( 'MRIread', 'file')
    error( 'Freesurfer toolbox is not available for MRIread.');
end


options = struct( 'brainfile', '',...
    'figdir', '',...
    'close_fig', 'no',... % close figure if it's saved
    'threshold', [0.95, 1],...
    'cm', autumn( 200),...
    'nbticks', 5,...
    'tick_precision', 2,...
    'brain_threshold', [],... % 1 x 2 matrix
    'brain_background', 'w',... % only works with 1mm and 2mm brain
    'img_sz', 0.15,... % figure size 0-1, relative to your screen size
    'marker_size', 10,...
    'show_marker', 'no',...% only when x,y&z are all set
    'marker_type', 'dot',... % 'dot' | 'cross' | 'dot+cross
    'hold', 'off',...
    'outline_overlay', 'no',...
    'outline_only', 'no',...
    'outline_color', [],...
    'outline_linewidth', 1,...
    'overlay_transparency', 1,...
    'coord_style', 'xyz',...%'xyz', 'x', 'y', 'z'
    'cross', [],...% N x 3, this can be used when coordinates are not given as [x, y, z] which is a prerequest for show_marker
    'cross_type', 'cross',...% 'dot' | 'cross' | 'dot+cross
    'rect', [],... % to make a zoomed-in version, 3x2 matrix
    'rect_trim', 'no',...
    'rect_color', 'k',...
    'rect_linewidth', 1.5,...    
    'atlas', '',... % each region is indicated by an integer number
    'atlas_value', [],...
    'atlas_linewidth', 1,...
    'atlas_hemi', 'lr',...
    'atlas_color', [],... % Nx3 matrix of RGB color
    'text_color', 'b', ... % label text color
    'fontsize', 6,...
    'refline', [],...
    'refline_view', 'x',...%'x' | 'y' | 'z'
    'refline_style', '-',...
    'refline_color', 'k',...
    'refline_width', 0.5',...
    'show_orientlabel', 'yes',... 'yes'|'no'
    'show_colorbar', 'no'); % 'yes' | 'no' 
     
options.outline_color = {'k', 'k'};
options = G_SparseArgs( options, varargin);

% make sure all input arguments are valid
if ~isscalar( options.img_sz) || options.img_sz < 0 || options.img_sz > 1
    % use default image size if user-specified value is out of range
    options.img_sz = 0.15;
end

% directory to save the figures
if ~isempty( options.figdir) && ~exist( options.figdir, 'dir')
    [status, msg] = mkdir( options.figdir);
    if status == 0
        warning( ['Failed to create figure saving directory:', msg]);
        options.figdir = '';
    end
end

switch options.close_fig
    case {'yes', 'no'}
        % do nothing
    otherwise
        options.close_fig = 'no';
end

% transparency of overlays
if ~isscalar( options.overlay_transparency) || options.overlay_transparency < 0 || options.overlay_transparency > 1
    options.overlay_transparency = 1;
end

% colorbar ticks
if ~isscalar( options.nbticks) || abs( options.nbticks - round( options.nbticks)) > eps
    options.nbticks = 5;
end
if ~isscalar( options.tick_precision) || abs( options.tick_precision - round( options.tick_precision)) > eps || ...
        options.tick_precision < 0 || options.tick_precision > 10
    options.tick_precision = 2;
end

cm = options.cm;
cm_len = size( cm, 1);
nb_cm = size( cm, 3);
nbticks = options.nbticks;
nbticks = min( [cm_len, nbticks]);
tick_precision = options.tick_precision;
if ~(isreal( tick_precision) && rem(tick_precision, 1) == 0)
    tick_precision = 2;
end

% threshold must be 1 x 2 or 2 x 2 matrix
% [negative_lo, negative_hi], [positive_lo, positive_hi], [negative_lo, negative_hi; positive_lo, positive_hi]
threshold = options.threshold;
if isvector( threshold)
    threshold = reshape( threshold(:), [1, length( threshold)]);
end

sz_thresh = size( threshold);
nb_threshold = sz_thresh( 1);
if ~isempty( threshold)
    if sz_thresh( 2) ~= 2 || sz_thresh( 1) > 2
        error( 'threshold must be a 1x2 or 2x2 matrix.');
    end

    % threshold is given as [lo, hi]
    if nb_cm ~= nb_threshold
        error( 'Number of colormaps must be equal to the number of thresholds.');
    end
end

for k = 1 : nb_threshold
    if threshold( k, 1) * threshold( k, 2) < 0
        error( 'Threshold must be set for positive and negative separately.');
    end
    
    if threshold( k, 1) < 0
        threshold( k, :) = sort( threshold( k, :), 'descend');
    else
        threshold( k, :) = sort( threshold( k, :));
    end
end

% read stats image
stats = [];
if ~isempty( statsfile)
    stats = MRIread( statsfile);
    stats.vol = permute( stats.vol, [2, 1, 3]);
else
    % you may just want to plot some brain slices
    if ~isempty( coord2plot) && isempty( options.brainfile)
        error( 'To plot plain brain slices, ''brainfile'' can not be empty.');
    end
end

% coordinates could be given as [x1, y1, z1; x2, y2, z2;...] and [x1 x2 x3]
% , [y1, y2, y3, ...] or [z1 z2 z3 ...],'coord_style must be specified in
% the latter case
if ~isempty( coord2plot) && ~strcmpi( options.coord_style, 'xyz')
    if ~isvector( coord2plot)
        error( 'Axial, sagittal, or cronal slices index must be given as a vector.');
    else
        nb_coords = length( coord2plot);
        switch lower( options.coord_style)
            case 'x'
                coord2plot = [coord2plot(:), zeros( nb_coords, 2)];
            case 'y'
                coord2plot = [zeros( nb_coords, 1), coord2plot(:), zeros( nb_coords, 1)];
            case 'z'
                coord2plot = [zeros( nb_coords, 2), coord2plot(:)];
            otherwise
                error( 'Unsupported coordinates input style.');
        end
    end
end

% determine which standard brain to use
fsldir = getenv( 'FSLDIR');
if isempty( fsldir) && exist( '/usr/local/fsl', 'dir')
    fsldir = '/usr/local/fsl';
end

if ~isempty( options.brainfile)
    if ~strcmpi( options.brainfile( max([1, end-6]) : end), '.nii.gz')
        % .nii.gz 
        options.brainfile = [options.brainfile, '.nii.gz'];
    end
    [p, brain_name] = fileparts( options.brainfile);
    [~, brain_name] = fileparts( brain_name);
    if isempty( p)
        if ~isempty( fsldir)
            options.brainfile = fullfile( fsldir, 'data', 'standard', [brain_name, '.nii.gz']);
        else
            options.brainfile = fullfile( pwd, [brain_name, '.nii.gz']);
        end        
        if ~exist( options.brainfile, 'file')
            error( [brain_name, ' was not found in standard brain atlas or current working path.']);
        end
    end 
end

if isempty( options.brainfile) && ~isempty( stats)
    % determine which MNI brain to use automatically
    if ~isempty( fsldir)
        if abs( stats.volres(1) - 0.5) < eps
            options.brainfile = fullfile( fsldir, 'data', 'standard', 'MNI152_T1_0.5mm.nii.gz');            
        elseif stats.volres(1) == 1
            options.brainfile = fullfile( fsldir, 'data', 'standard', 'MNI152_T1_1mm_brain.nii.gz');            
        elseif stats.volres(1) == 2
            options.brainfile = fullfile( fsldir, 'data', 'standard', 'MNI152_T1_2mm_brain.nii.gz');
        else
            error( ['No ', num2str(stats.volres(1)), ' mm MNI brain was found']);
        end
        
    else
        error( '''FSLDIR'' was not set. No MNI standard brain was available');
    end 
    
elseif ~isempty( options.brainfile) && ~isempty( stats)
    % background image was set, but with  a different resolution
    % resample the stats image using FSL
    [~, brain_name] = fileparts( options.brainfile);
    tmp_loc = strfind( brain_name, 'mm');
    if ~isempty( tmp_loc)
        brain_resol = brain_name( tmp_loc-1 : tmp_loc+1);
        if abs( stats.volres(1) - str2double( brain_resol(1))) > eps
            [p, n, ext] = fileparts( statsfile);
            statsfile_new_resol = fullfile( p, [brain_resol, '_', n, ext]);
            cmd = ['source ${FSLDIR}/etc/fslconf/fsl.sh && export PATH=$PATH:${FSLDIR}/bin && flirt -in ', statsfile,...
            ' -ref ', options.brainfile, ' -interp nearestneighbour ', ' -applyxfm -usesqform -out ', statsfile_new_resol, ' -omat ', [statsfile_new_resol, '.mat 1>/dev/null']];
            status = system( cmd);

            % re-load stats image
            if status == 0
                if ~isempty( coord2plot)
                    mni_coord2plot = stats.vox2ras1 * [coord2plot'; ones( 1, size( coord2plot, 1))];
                end
                
                % cross point
                if ~isempty( options.cross)
                    options.cross = stats.vox2ras1 * [options.cross'; ones( 1, size( options.cross, 1))];
                end
                
                % rect
                if ~isempty( options.rect)
                    options.rect = stats.vox2ras1 * [options.rect; ones( 1, size( options.rect, 2))];
                end
                
                statsfile = statsfile_new_resol;
                stats = MRIread( statsfile);
                stats.vol = permute( stats.vol, [2, 1, 3]);

                % coordinates in new space
                if ~isempty( coord2plot)
                    coord2plot = stats.vox2ras1 \ mni_coord2plot;
                    coord2plot = coord2plot( 1:3, :)';
                end
                
                % cross point
                if ~isempty( options.cross)
                    options.cross = stats.vox2ras1 \ options.cross;
                    options.cross = round( options.cross(1:3))';
                end
                
                % rect
                if ~isempty( options.rect)
                    options.rect = stats.vox2ras1 \ options.rect;
                    options.rect = round( options.rect( 1:3, :));
                end  

            else
                error( ['Failed to resample the stats image to a new resolution of ', brain_resol]);
            end
        end
    end
    
else
    % do nothing
end

% % directory to save figures
% % all figures will be saved to file and closed, otherwise they will be left open
% if isempty( options.figdir) && ~isempty( statsfile)    
%     p = fileparts( statsfile);
%     p = fullfile( p, 'Figures');
%     if ~exist( p, 'dir')
%         status = mkdir( p);
%         if status == 1
%             options.figdir = p;
%         end
%     else
%         options.figdir = p;
%     end
% end

% read template brain image
brain = MRIread( options.brainfile);
brain.vol = permute( brain.vol, [2, 1, 3, 4]);
    
% make sure brain threshold is valid
brain_threshold = options.brain_threshold;
if isempty( brain_threshold) || ~isvector( brain_threshold) || length( brain_threshold) ~= 2
    if ~isempty( strfind( options.brainfile, 'MNI152_'))
        brain_threshold = [3, 8] * 1000; 
    else
        brain_threshold = prctile( brain.vol(:), [5, 95]);
    end
    
else
    lo_intensity = min( brain.vol(:));
    hi_intensity = max( brain.vol(:));
    brain_threshold = sort( brain_threshold);
    if brain_threshold(1) > hi_intensity || brain_threshold(2) < lo_intensity
        warning( 'Invalid brain threshold, automatic threshold will be used.\n');
        brain_threshold = prctile( brain.vol(:), [5, 95]);
    end
end

% read atlas
if ~isempty( options.atlas)
    atlas = MRIread( options.atlas);
    atlas_transf = atlas.vox2ras1;
    atlas.vol = permute( atlas.vol, [2, 1, 3, 4]);
    % make sure atlas and brain of the same dimension  
    if ndims( brain.vol) ~= ndims( atlas.vol) || ( ~all( size( brain.vol) == size( atlas.vol)) )
        error( 'Brain and atlas may not be in the same space.');
    end
    
    % set non-atlas voxels to 0
    if ~isempty( options.atlas_value)
        atlas.vol( ~ismember( atlas.vol, options.atlas_value)) = 0;
    end
    
    mid_slice = atlas_transf \ [0, 0, 0, 1]';
    mid_slice = mid_slice(1);
    left_hemi_slice = atlas_transf \ [-1, 0, 0, 1]';
    left_hemi_slice = left_hemi_slice( 1);
    
    switch options.atlas_hemi        
        case 'l'
            if left_hemi_slice > mid_slice
                atlas.vol( 1: round( mid_slice), :, :) = 0;
            else
                atlas.vol( round( mid_slice) : end, :, :) = 0;
            end
        case 'r'
            if left_hemi_slice > mid_slice
                atlas.vol( round( mid_slice) : end, :, :) = 0;
            else
                atlas.vol( 1: round( mid_slice), :, :) = 0;
            end
        case 'lr'
            % do nothing
        otherwise
            % do nothing
    end
        
    unique_atlas = unique( atlas.vol(:));
    unique_atlas( unique_atlas == 0) = [];
    
    if isempty( options.atlas_color)
        if length( unique_atlas) < 4
            options.atlas_color = ...
                [1 0 0;...
                0 0 1;...
                0 1 0;...
                1 1 1];
        else
            options.atlas_color = ...
                cat( 1,...
                [1 0 0; 0 0 1; 0 1 0; 1 1 1],...
                rand( length( unique_atlas) - 4, 3));
        end
    end
end

% pad matrix
sz = size( brain.vol);
max_sz = max( sz);
pre_pad = zeros( 3, 1);
post_pad = zeros( 3, 1);
for k = 1 : 3
    total_pad = max_sz - sz( k);
    pre_pad( k) = floor( total_pad/2);
    post_pad( k) = total_pad - pre_pad( k);
end

pad_brain = padarray( brain.vol, pre_pad, brain.vol( 1, 1, 1), 'pre');
pad_brain = padarray( pad_brain, post_pad, brain.vol( 1, 1, 1), 'post');

if ~isempty( options.atlas)
    atlas = padarray( atlas.vol, pre_pad, brain.vol( 1, 1, 1), 'pre');
    atlas = padarray( atlas, post_pad, brain.vol( 1, 1, 1), 'post');
end

if ~isempty( options.refline)
    switch options.refline_view
        case 'x'
            options.refline = options.refline + pre_pad(1);
        case 'y'
            options.refline = options.refline + pre_pad(2);
        case 'z'
            options.refline = options.refline + pre_pad(3);
        otherwise
            options.refline = [];
    end
end

% brain background color
if strcmpi( options.brain_background, 'w')
    mval = max( pad_brain(:));
    pad_brain = ((pad_brain + 1e-8) < 2e-8) * max( [brain_threshold(2), mval])*2 + pad_brain;
end

a = findobj( 'type', 'axes');
if ~isempty( a)
    a = gca;
end

% make colorbar
is_discrete_cm = 'no';
if ~isempty( stats)
    if ~all( brain.volsize == stats.volsize)
        error( 'stats and background images are not in the same standard MNI space.');
    end
    
    pad_stats = padarray( stats.vol, pre_pad, 0, 'pre');
    pad_stats = padarray( pad_stats, post_pad, 0, 'post');
    
    if nb_cm == 1 && size( cm, 1) == (length( unique( stats.vol(:))) - 1 )
        is_discrete_cm = 'yes';
    end
    
    if strcmpi( is_discrete_cm, 'yes')
        stats_val = unique( stats.vol(:));
        stats_val = stats_val( stats_val ~= 0);

        if strcmpi( options.show_colorbar, 'yes')
            figure( 'color', 'w', 'units', 'normalized');
            w = 0.03;
            SetPrintProp( gcf, w*nb_cm, 0.12); 
            
            imagesc( permute( cm, [1 3 2]));
            ctick = 1 : length( stats_val);
            cticklabel = num2cell( stats_val);
            cticklabel = cellfun( @(x) num2str(x), cticklabel, 'UniformOutput', false);
            set( gca, 'box', 'on', 'XTick', [], 'ydir', 'normal', 'YTick', ctick,...
                'YTickLabel', cticklabel, 'LineWidth', 0.5,...
                'FontSize', options.fontsize, 'Units', 'normalized',...
                'TickLength', [1, 1]*w*nb_cm*0.25);
        end
        
    else
        if strcmpi( options.show_colorbar, 'yes')
            figure( 'color', 'w', 'units', 'normalized');
            w = 0.03;
            SetPrintProp( gcf, w*nb_cm, 0.12); 
            for k = 1 : nb_cm  
                subplot( 1, nb_cm, k);
                imagesc( permute( cm( :, :, k), [1 3 2]));
                ctick = linspace( 1, cm_len, nbticks);
                cticklabel = linspace( threshold( k, 1), threshold( k, 2), nbticks);
                cticklabel = round( cticklabel*( 10 ^tick_precision))/ (10 ^tick_precision);
                cticklabel = num2cell( cticklabel);
                cticklabel = cellfun( @(x) num2str(x), cticklabel, 'UniformOutput', false);
                set( gca, 'box', 'on', 'XTick', [], 'ydir', 'normal', 'YTick', ctick,...
                    'YTickLabel', cticklabel, 'LineWidth', 0.5,...
                    'FontSize', options.fontsize, 'Units', 'normalized',...
                    'TickLength', [1, 1]*w*nb_cm*0.25);
            end
        end        
    end

    if strcmpi( options.show_colorbar, 'yes') && ~isempty( options.figdir)
        try
            print( gcf, '-fillpage', '-dpdf', fullfile( options.figdir, 'colorbar.pdf'));
        catch
            print( gcf, '-dpdf', fullfile( options.figdir, 'colorbar.pdf'));
        end
    end
   
else
    pad_stats =[];
end

if ~isempty( a)
    axes( a);
end

% update sz
sz = size( pad_brain);

% build color lookup table
stats_rgb = [];
% slice index to plot
mx = [];
my = [];
mz = [];
if ~isempty( stats)    
    mask_rgb = [];
    mask_ind = [];
    stats_mask = zeros( size( pad_stats)) > 1;
    
    if strcmpi( is_discrete_cm, 'no')
        for k = 1 : nb_threshold    
            cur_threshold = threshold( k, :);

            if cur_threshold(1) <= 0 
                cur_stats_mask = pad_stats <= cur_threshold( 1);
            else
                cur_stats_mask = pad_stats >= cur_threshold( 1);
            end

            tmp_mask_ind = find( cur_stats_mask(:));        
            if ~isempty( tmp_mask_ind)
                stats2rend = pad_stats( tmp_mask_ind);   

                if cur_threshold(1) <= 0
                    max_val = max( stats2rend);  % threshold(1)  
                    if min( stats2rend) < cur_threshold(2)
                        stats2rend( stats2rend < cur_threshold( 2)) = cur_threshold(2);
                    end
                    % rescale to [1, size( cm, 1)]
                    rescale_range = max_val - cur_threshold(2);
                    if abs( rescale_range) < eps
                        stats2rend_ind = cm_len * ones( length( tmp_mask_ind), 1);
                    else
                        stats2rend_ind = round( (( max_val - stats2rend) * (cm_len - 1)) / (max_val - cur_threshold(2)) + 1);
                    end

                else
                    min_val = min( stats2rend);  % threshold(1)  
                    if max( stats2rend) > cur_threshold(2)
                        stats2rend( stats2rend > cur_threshold(2)) = cur_threshold(2);
                    end
                    % rescale to [1, size( cm, 1)]
                    rescale_range = cur_threshold(2) - min_val;
                    if abs( rescale_range) < eps
                        stats2rend_ind = cm_len * ones( length( tmp_mask_ind), 1);
                    else
                        stats2rend_ind = round( (( stats2rend - min_val) * (cm_len - 1)) / (cur_threshold(2) - min_val) + 1);
                    end
                end

                mask_rgb = cat( 1, mask_rgb, cm( stats2rend_ind, :, k));
                mask_ind = cat( 1, mask_ind, tmp_mask_ind);
                stats_mask = stats_mask | cur_stats_mask;
            end
        end
        
    else
        stats_val = unique( stats.vol(:));
        stats_val = stats_val( stats_val ~= 0);
        for k = 1 : length( stats_val)
            ind = find( abs( pad_stats - stats_val(k)) < eps);
            mask_ind = cat( 1, mask_ind, ind);            
            mask_rgb = cat( 1, mask_rgb, repmat( cm( k, :), [length(ind), 1]));
        end
        
        stats_mask = pad_stats ~= 0;
    end
    
    % color look up table for each voxel within mask
    stats_rgb = zeros( sz(1) * sz(2) * sz(3), 3);
    stats_rgb( mask_ind, :) = mask_rgb;
    stats_rgb = reshape( stats_rgb, sz(1), sz(2), sz(3), 3);
    
    % range to plot
    [mx, my, mz] = ind2sub( sz, mask_ind);
    mx = unique( mx);
    my = unique( my);
    mz = unique( mz);
end

% default: plot all slices if coord2plot was not set
if ~isempty( coord2plot)
    % coordinates are specified
    mx = unique( coord2plot( :, 1));
    mx( mx <= 0 | mx > max_sz) = [];
    
    my = unique( coord2plot( :, 2));
    my( my <= 0 | my > max_sz) = [];
    
    mz = unique( coord2plot( :, 3)); 
    mz( mz <= 0 | mz > max_sz) = [];
    
else
    % all axial, sagittal, and cronal slices
    % check the number of figures to plot
    if sum( [length(mx), length( my), length( mz)]) > 100 && isempty( options.figdir)
        val = questdlg( 'Too many figures to be plotted.\n Are you sure to continue?');
        if any( ismember( val, {'no', 'cancel'}))
            exit;
        end
    end
end

% show marker at (x, y, z)
plot_coord = options.show_marker;
if length( mx) == length( my) && length( mx) == length( mz) && length( my) == length( mz)
    %plot_coord = 'yes';
    % do nothing
else
    plot_coord = 'no';
end

% plot images
view_name = {'X', 'Y', 'Z'};
if ~isempty( options.cross)
    if isvector( options.cross)
        options.cross = transpose( options.cross(:));
    end
    
    if size( options.cross, 2) ~= 3
        error( 'Cross locations must be given as a Nx3 matrix.');
    end
    
    options.cross = bsxfun( @plus, options.cross, pre_pad(:)');
end

if ~isempty( options.rect)
    orig_rect = options.rect;
    options.rect = options.rect + repmat( pre_pad(:), [1, 2]);
end

% plot sagittal, cronal and axial slices
for k = 1 : 3
    if k == 1
        % sagittal
        cur_vox = mx;
        if strcmpi( plot_coord, 'yes')
            if ~isempty( coord2plot)
                % user-input coordinates, which should be padded
                vox_coord = [my(:)+pre_pad(2), mz(:)+pre_pad(3)];
            else
                vox_coord = [my(:), mz(:)];
            end
        end
        
        % plot cross
        if ~isempty( options.cross)
            cross_coord = options.cross( :, [2, 3, 1]);
        end
        
    elseif k == 2
        % cronal
        cur_vox = my;
        if strcmpi( plot_coord, 'yes')
            if ~isempty( coord2plot)
                vox_coord = [mx(:)+pre_pad(1), mz(:)+pre_pad(3)];
            else
                vox_coord = [mx(:), mz(:)];
            end
        end
        
        if ~isempty( options.cross)
            cross_coord = options.cross( :, [1, 3, 2]);
        end
        
    else
        % axial
        cur_vox = mz;
        if strcmpi( plot_coord, 'yes')
            if ~isempty( coord2plot)
                vox_coord = [mx(:)+pre_pad(1), my(:)+pre_pad(2)];
            else
                vox_coord = [mx(:), my(:)];
            end
        end
        
        if ~isempty( options.cross)
            cross_coord = options.cross;
        end        
    end

    if ~isempty( cur_vox)
        if ~isempty( coord2plot)
            cur_vox = cur_vox + pre_pad( k);
        end

        for vox_idx = 1 : length( cur_vox)
            slice_idx = cur_vox( vox_idx);

            if strcmpi( options.hold, 'off')
                figure( 'color', 'w', 'units', 'normalized');
                SetPrintProp( gcf, options.img_sz*0.7, options.img_sz);
            end

            if strcmpi( plot_coord, 'yes')
                current_coord = vox_coord( vox_idx, :);
            else
                current_coord = [];
            end
            
            % background brain image
            if k == 1
                brain_img = squeeze( pad_brain( slice_idx, :, :))';
                if ~isempty( stats_rgb)
                    cur_rgb = squeeze( stats_rgb( slice_idx, :, :, :));
                    cur_rgb = permute( cur_rgb, [2 1 3]);
                    cm_mask = squeeze( stats_mask( slice_idx, :, :))';
                    cur_stats = squeeze( pad_stats( slice_idx, :, :))';
                end

                img_coord = [mx( vox_idx); 1; 1];
                if ~isempty( my) && length( my) >= vox_idx
                    img_coord(2) = my( vox_idx);
                end
                if ~isempty( mz) && length( mz) >= vox_idx
                    img_coord(3) = mz( vox_idx);
                end
                
                if ~isempty( options.atlas)
                    atlas_img = squeeze( atlas( slice_idx, :, :))';
                end

            elseif k == 2
                brain_img = squeeze( pad_brain( :, slice_idx, :))';
                if ~isempty( stats_rgb)
                    cur_rgb = squeeze( stats_rgb( :, slice_idx, :, :));
                    cur_rgb = permute( cur_rgb, [2 1 3]);
                    cm_mask = squeeze( stats_mask( :, slice_idx, :))';
                    cur_stats = squeeze( pad_stats( :, slice_idx, :))';
                end

                img_coord = [1; my( vox_idx); 1];
                if ~isempty( mx) && length( mx) >= vox_idx
                    img_coord(1) = mx( vox_idx);
                end
                if ~isempty( mz) && length( mz) >= vox_idx
                    img_coord(3) = mz( vox_idx);
                end
                
                if ~isempty( options.atlas)
                    atlas_img = squeeze( atlas( :, slice_idx, :))';
                end

            else
                brain_img = squeeze( pad_brain( :, :, slice_idx))';
                if ~isempty( stats_rgb)
                    cur_rgb = squeeze( stats_rgb( :, :, slice_idx, :));
                    cur_rgb = permute( cur_rgb, [2 1 3]);
                    cm_mask = squeeze( stats_mask( :, :, slice_idx))';
                    cur_stats = squeeze( pad_stats( :, :, slice_idx))';
                end

                img_coord = [1; 1; mz( vox_idx)];
                if ~isempty( my) && length( my) >= vox_idx
                    img_coord(1) = my( vox_idx);
                end
                if ~isempty( my) && length( my) >= vox_idx
                    img_coord(2) = my( vox_idx);
                end
                
                if ~isempty( options.atlas)
                    atlas_img = squeeze( atlas( :, :, slice_idx))';
                end
            end
                
            % plot brain
            cm_h = imagesc( brain_img);
            colormap( gray);
            caxis( brain_threshold);
            
            % overlay
            if ~isempty( stats_rgb)
                hold on;
                if strcmpi( options.outline_only, 'no')
                    overlay_cm_h = imagesc( cur_rgb);
                    set( overlay_cm_h, 'AlphaData', options.overlay_transparency * cm_mask);
                end
                
                % outline overlay
                if strcmpi( options.outline_overlay, 'yes')
                    if ~iscell( options.outline_color)
                        options.outline_color = repmat( {options.outline_color}, [1, 2]);
                    else
                        if length( options.outline_color) == 1
                            options.outline_color = repmat( options.outline_color, [1, 2]);
                        end
                    end
                    
                    for thresh_idx = 1 : nb_threshold
                        cur_threshold = threshold( thresh_idx, :);

                        if cur_threshold(1) <= 0 
                            G_Outline_Pixels( cm_h, cur_stats <= cur_threshold( 1),...
                                'linecolor', options.outline_color{1},...
                                'linewidth', options.outline_linewidth); 
                        else
                            G_Outline_Pixels( cm_h, cur_stats >= cur_threshold( 1),...
                                'linecolor', options.outline_color{2},...
                                'linewidth', options.outline_linewidth);
                        end
                    end
                end
            end

            % plot dot
            if ~isempty( current_coord)
                hold on;                
                if any( ismember( lower( options.marker_type),  {'dot+cross', 'cross'}))
                    tmp_shift = 0.5;
                    plot( [-1, current_coord(1) - tmp_shift],...
                        [1,1]* current_coord(2), 'Color', [1, 1, 1] * 0.3);
                    plot( [current_coord(1) + tmp_shift, max_sz+1],...
                        [1,1]* current_coord(2), 'Color', [1, 1, 1] * 0.3);
                    plot( current_coord(1) *[1, 1],...
                        [1, current_coord(2) - tmp_shift],  'Color',[1, 1, 1] * 0.3);
                    plot( current_coord(1) *[1, 1],...
                        [ current_coord(2) + tmp_shift, max_sz+1], 'Color', [1, 1, 1] * 0.3);
                end
                
                if any( ismember( lower( options.marker_type),  {'dot+cross', 'dot'}))
                    plot( current_coord(1), current_coord(2), 'ro',...
                        'MarkerFaceColor', 'r',...
                        'MarkerEdgeColor', 'r',...
                        'markersize', options.marker_size);    
                end
            end
            
            % plot cross
            if ~isempty( options.cross)
                hold on;
                tmp_shift = 0.5;
                for coord_idx = 1 : size( cross_coord, 1)
                    if any( ismember( lower( options.cross_type),  {'dot+cross', 'cross'}))
                        plot( [-1, cross_coord(coord_idx, 1) - tmp_shift],...
                            [1,1]* cross_coord( coord_idx, 2),...
                            'Color', [1, 1, 1] * 0.3);
                        plot( [cross_coord( coord_idx, 1) + tmp_shift, max_sz+1],...
                            [1,1]* cross_coord( coord_idx, 2),...
                            'Color', [1, 1, 1] * 0.3);
                        plot( cross_coord( coord_idx, 1) *[1, 1],...
                            [1, cross_coord( coord_idx, 2) - tmp_shift],...
                            'Color',[1, 1, 1] * 0.3);
                        plot( cross_coord( coord_idx, 1) *[1, 1],...
                            [ cross_coord( coord_idx, 2) + tmp_shift, max_sz+1],...
                            'Color', [1, 1, 1] * 0.3);
                    end
                    
                    if any( ismember( lower( options.cross_type),  {'dot+cross', 'dot'}))
                        plot( cross_coord( coord_idx, 1), cross_coord( coord_idx, 2), 'ro',...
                            'MarkerFaceColor', 'r',...
                            'MarkerEdgeColor', 'r',...
                            'markersize', options.marker_size);    
                    end
                end
            end
            
            % atlas outline
            atlas_xlim = [];
            atlas_ylim = [];
            if ~isempty( options.atlas)
                atlas_idx = unique( atlas_img( :));
                atlas_idx( atlas_idx == 0) = [];
                if ~isempty( atlas_idx)
                    hold on;
                    atlas_indx = zeros( length( atlas_idx), 2);
                    atlas_indy = zeros( length( atlas_idx), 2);
                    for tmp_ind = 1 : length( atlas_idx)
                        tmp_atlas = zeros( size( atlas_img));
                        tmp_atlas( atlas_img == atlas_idx( tmp_ind)) = 1;
                        color_loc = unique_atlas == atlas_idx( tmp_ind);
                        G_Outline_Pixels( cm_h, tmp_atlas,...
                            'linecolor', options.atlas_color( color_loc, :),...
                            'linewidth', options.atlas_linewidth); 
                        
                        [tmp_indx, tmp_indy] = ind2sub( size( tmp_atlas), find( tmp_atlas(:) > 0));
                        atlas_indx( tmp_ind, :) = [min( tmp_indx), max( tmp_indx)];
                        atlas_indy( tmp_ind, :) = [min( tmp_indy), max( tmp_indy)];
                    end
                    
                    atlas_xlim = [min( atlas_indx), max( atlas_indx)];
                    atlas_ylim = [min( atlas_indy), max( atlas_indy)];
                end
            end
                       
            axis equal
            axis off
            set( gca, 'YDir', 'normal', 'xlim', [-1, max_sz+1], 'ylim', [-1, max_sz+1],...
                'xtick', [], 'ytick', [], 'ztick', [], 'Box', 'off');                 
                       
            % remove extra white space
            if brain_img( 1, 1) > 0
                [indx, indy] = ind2sub( size( brain_img), find( brain_img(:) < brain_img( 1, 1)));
            else
                [indx, indy] = ind2sub( size( brain_img), find( brain_img(:) > 0));
            end
            
            if ~isempty( indx)    
                % indx_limits = [min( indx)-1, max( indx)+1];
                % indy_limits = [min( indy)-1, max( indy)+1];
                % set( gca, 'xlim', indy_limits);
                % plot( indy_limits( [1, 2, 2, 1 1]), indx_limits( [1 1 2 2 1]), options.text_color);
                indx_limits = [ min( [atlas_xlim, min( indx) - 1.5]),...
                    max( [atlas_xlim, max( indx) + 1.5]) ];
                indy_limits = [ min( [atlas_ylim, min( indy) - 1.5] ),...
                    max( [atlas_ylim, max( indy) + 1.5]) ];               
                ax_pos = get( gca, 'Position');
                raw_ylim = diff( get( gca, 'YLim'));
                raw_xlim = diff( get( gca, 'XLim'));
                vertical_c = ax_pos(2) + ax_pos(4)/2;
                new_height = ax_pos(4) * (diff( indx_limits))/raw_ylim;
                hor_c = ax_pos(1) + ax_pos(3)/2;
                new_width = ax_pos(3)* (diff( indy_limits))/raw_xlim;
                new_pos = [hor_c - new_width/2, vertical_c - new_height/2, new_width, new_height];
                set( gca, 'Position', new_pos);
                set( gca, 'YLim', indx_limits, 'XLim', indy_limits);
                % plot( indy_limits( [1, 2, 2, 1 1]), indx_limits( [1 1 2 2 1]), options.text_color);
            end
            
            % zoom-in
            if ~isempty( options.rect)
                hold on;
                
                rect_mni = round( 10 * brain.vox2ras1 * [orig_rect; ones( 1, 2)]) / 10;
                rect_mni_str = num2cell( rect_mni);
                rect_mni_str = cellfun( @(x) num2str(x), rect_mni_str, 'UniformOutput', false);
                if k == 1
                    plot( options.rect( 2, [1 2 2 1 1]) + 0.5*[-1, 1, 1, -1, -1],...
                        options.rect( 3, [1 1 2 2 1]) + 0.5 * [-1 -1 1 1 -1],...
                        'linewidth', options.rect_linewidth,...
                        'color', options.rect_color);
                    if strcmpi( options.rect_trim, 'yes')
                        set( gca,...
                            'XLim', options.rect( 2, :) + 0.5*[-1, 1],...
                            'YLim', options.rect( 3, :) + 0.5*[-1, 1],...
                            'xtick', options.rect(2, :), 'xticklabel', rect_mni_str( 2, :),...
                            'ytick', options.rect( 3, :), 'yticklabel', rect_mni_str( 3, :),...
                            'tickdir', 'in', 'TickLength', [0, 0] * 0.03);
                        axis on
                        box on
                    end
                elseif k == 2
                    plot( options.rect( 1, [1 2 2 1 1]) + 0.5*[-1 1 1 -1 -1],...
                        options.rect( 3, [1 1 2 2 1]) + 0.5*[-1 -1 1 1 -1],...
                        'linewidth', options.rect_linewidth, 'color',options.rect_color);
                     if strcmpi( options.rect_trim, 'yes')
                         set( gca,...
                             'XLim', options.rect( 1, :) + 0.5*[-1, 1],...
                             'YLim', options.rect( 3, :) + 0.5*[-1, 1] ,...
                            'xtick', options.rect(1, :), 'xticklabel', rect_mni_str( 1, :),...
                            'ytick', options.rect( 3, :), 'yticklabel', rect_mni_str( 3, :),...
                            'tickdir', 'in', 'TickLength', [0, 0] * 0.03);
                        axis on
                        box on
                     end
                else
                    plot( options.rect( 1, [1 2 2 1 1]) + 0.5 * [-1 1 1 -1 -1],...
                        options.rect( 2, [1 1 2 2 1]) + 0.5 * [-1 -1 1 1 -1],...
                        'linewidth', options.rect_linewidth, 'color', options.rect_color);
                     if strcmpi( options.rect_trim, 'yes')
                         set( gca,...
                             'XLim', options.rect( 1, :) + 0.5*[-1, 1],...
                             'YLim', options.rect( 2, :) + 0.5*[-1, 1] ,...
                            'xtick', options.rect(1, :), 'xticklabel', rect_mni_str( 1, :),...
                            'ytick', options.rect( 2, :), 'yticklabel', rect_mni_str( 2, :),...
                            'tickdir', 'in', 'TickLength', [0, 0] * 0.03);
                        axis on
                        box on
                     end
                end
            end                        
                        
            % reference lines
            if ~isempty( options.refline)
                hold on;
                xlim = get( gca, 'xlim');
                ylim = get( gca, 'ylim');
                if k == 1
                    switch options.refline_view
                        case 'y'
                            for tmp_ind = 1 : length( options.refline)
                                plot( options.refline( tmp_ind) * [1, 1], ylim,...
                                    'linestyle', options.refline_style,...
                                    'color', options.refline_color,...
                                    'linewidth', options.refline_width);
                            end
                        case 'z'
                            for tmp_ind = 1 : length( options.refline)
                                plot( xlim, options.refline( tmp_ind) * [1, 1],...
                                    'linestyle', options.refline_style,...
                                    'color', options.refline_color,...
                                    'linewidth', options.refline_width);
                            end
                        otherwise
                            % nothing to plot
                    end
                            
                elseif k == 2
                    switch options.refline_view
                        case 'x'
                            for tmp_ind = 1 : length( options.refline)
                                plot( options.refline( tmp_ind) * [1, 1], ylim,...
                                    'linestyle', options.refline_style,...
                                    'color', options.refline_color,...
                                    'linewidth', options.refline_width);
                            end
                        case 'z'
                            for tmp_ind = 1 : length( options.refline)
                                plot( xlim, options.refline( tmp_ind) * [1, 1],...
                                    'linestyle', options.refline_style,...
                                    'color', options.refline_color,...
                                    'linewidth', options.refline_width);
                            end
                        otherwise
                            % nothing to plot
                    end
                    
                else
                    switch options.refline_view
                        case 'x'
                            for tmp_ind = 1 : length( options.refline)
                                plot( options.refline( tmp_ind) * [1, 1], ylim,...
                                    'linestyle', options.refline_style,...
                                    'color', options.refline_color,...
                                    'linewidth', options.refline_width);
                            end
                        case 'y'
                            for tmp_ind = 1 : length( options.refline)
                                plot( xlim, options.refline( tmp_ind) * [1, 1],...
                                    'linestyle', options.refline_style,...
                                    'color', options.refline_color,...
                                    'linewidth', options.refline_width);
                            end
                        otherwise
                            % nothing to plot
                    end
                end                
            end
            
            % voxel coordinate to MNI coordinate
            if isempty( coord2plot)
                % image coordinates are in padded matrix
                img_coord = img_coord - pre_pad(:);
            end
            mni_coord = round( 10 * brain.vox2ras1 * [img_coord; 1]) / 10;
            
            % display labels
            xlim = get( gca, 'xlim');
            ylim = get( gca, 'YLim');
            
            left_margin = xlim(1) - 4;
            right_margin = xlim(2) + 2;
            top_margin = ylim(2) + 4;
            text_center = diff( ylim)/2 + ylim(1);            
            % % Uncomment the four lines below to use fixed text position
            % left_margin = xlim-5;
            % right_margin = max_sz + 1;
            % top_margin = max_sz - 10;
            % text_center = max_sz/2;
            % show coordinate
            if k == 1
                t = text( left_margin, top_margin, ['x = ', num2str( mni_coord(k)), 'mm']);
            elseif k == 2
                t = text( left_margin, top_margin, ['y = ', num2str( mni_coord(k)), 'mm']);
            else
                t = text( left_margin, top_margin, ['z = ', num2str( mni_coord(k)), 'mm']);
            end            
            set( t, 'Color', options.text_color,...
                'HorizontalAlignment', 'left',...
                'fontsize', options.fontsize);
       
            
            % show orientation label
            if strcmpi( options.show_orientlabel, 'yes')
                if k == 1
                    tlabel = text( left_margin, text_center, 'P');
                    t2label = text( right_margin, text_center, 'A');
                elseif k == 2                    
                    tlabel = text( left_margin, text_center, 'R');
                    t2label = text( right_margin, text_center, 'L');
                else
                    tlabel = text( left_margin, text_center, 'R');
                    t2label = text( right_margin, text_center, 'L');                    
                end

                set( tlabel,...
                    'Color', options.text_color,...
                    'HorizontalAlignment', 'right',...
                    'fontsize', options.fontsize);
                set( t2label,...
                    'Color', options.text_color,...
                    'fontsize', options.fontsize);
            end
            
            if pad_brain( 1, 1, 1) > brain_threshold(2)
                tmp = ones( size( brain_img));
                tmp( brain_img > pad_brain( 1, 1, 1) * 0.7) = 0;
                cm_h.AlphaData = tmp;
            end
            
            % save figure
            if strcmpi( options.hold, 'off') && ~isempty( options.figdir)                 
                fig_name2save = [view_name{k}, '_',...
                    num2str( mni_coord( 1)), '_',...
                    num2str( mni_coord( 2)), '_',...
                    num2str( mni_coord( 3)), '.pdf'];
                try
                    print( gcf, '-dpdf', '-fillpage', fullfile( options.figdir, fig_name2save));
                catch
                    print( gcf, '-dpdf', fullfile( options.figdir, fig_name2save));
                end
                
                if strcmpi( options.close_fig, 'yes')
                    close( gcf);
                end
            end
            
        end % slice
    end
end % view 
