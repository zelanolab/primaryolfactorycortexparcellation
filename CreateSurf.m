function [surf_file, cortex_file, cortex] = CreateSurf( vol, varargin)
% Create surface mask and cortex rendering
% 
% Input
%   vol, full-path-to-original_MR e.g. wkpath/MR.nii 
%        or wkpath/MR.nii.gz
%        or 3d matrix
% 
%   Key-Value pairs
%       'surf_sm', surface smoothing kernel, default 6
%       'surf_threshold', % surface threshold, default 0.1
%       'rend_sm', %smoothing parameter for rendering gray (dim of lattice), default 2
%       'isoparm', % value of desired isosurface after smoothing default 0.65
%                0-1 if render_on_raw is set to 'no' (default), otherwise
%                it's a value between the lowest and highest values of the
%                raw brain image
%       'fold', name of the folder to save the results wkpath/fold
%       're_run', 'yes' | 'no' (run without checking file existence)
%       'render_on_raw', 'yes' | 'no'
%       'transform', option must be set when vol is a 3d matrix
% 
% Output
%    wkpath/MR_surface.nii, surface mask image
%    wkpath/MR_surface.nii.pdf, overlayed surface on raw image
%    wkpath/MR_cortex.mat, cortex rendering
% 
% surfcae reconstruction part was modified from ctmr_guipackage Copyright (C) 2009  D. Hermes & K.J. Miller, Dept of Neurology and Neurosurgery, University Medical Center Utrecht
% 
% 
% See also spm_segment_job.m, 
% 


options = struct( 'surf_sm', 2,... % surface smoothing kernel
    'surf_threshold', 0.01,... % surface threshold
    'rend_sm', 1,... % smoothing parameter for rendering gray (dim of lattice), sigma
    'isoparm', 0.75,...% value of desired isosurface after smoothing
    'fold', 'SurfCortex',...
    're_run', 'no',...
    'render_on_raw', 'no',... % use raw image for rendering instead of tissue probablity imges
    'transform', []); % voxel -> real word transformation matrix

option_names = fieldnames( options);

nbargus = length( varargin);
if mod( nbargus, 2) ~= 0
   error('Key-Value pairs needed.')
end

for pair = reshape( varargin, 2, []) 
   in_name = lower( pair{1}); 
   if any( strcmpi( in_name, option_names))
      options.( in_name) = pair{2}; 
   else
      error('%s is not a recognized parameter name.', in_name)
   end   
end

smooth_kernel = options.surf_sm;
threshold = options.surf_threshold;
isoparm = options.isoparm;
sm_par = options.rend_sm;

if ~isscalar( smooth_kernel) || ~isscalar( threshold) ||  ~isscalar( isoparm) ||  ~isscalar( sm_par)
    error( 'All input arguments must be scalars.');
end

if smooth_kernel < 0 || threshold <= 0 || threshold >= 1 || sm_par < 0
    error( 'Invalid input arguments');
end

%% if input vol is a matrix
if isnumeric( vol)
    mat = options.transform;
    
    if isempty( mat)
        error( 'transform matrix must be set.');
        
    elseif isnumeric( mat)
        if length( size( mat)) ~= 2 || size( mat, 1) ~= 4 || size( mat, 2) ~= 4
            error( 'Input transform was not likely to be a transformation matrix.');
        end

    elseif ischar( mat)
        try
            brain_info = spm_vol( mat);
            mat = brain_info.mat;
            
        catch
            error( ['Failed to read transform matrix from file', options.transform]);
        end

    else
        error( 'Unknown type for transform');
    end
    
    % gaussian smoothing
    if sm_par > 0
        vol = imgaussfilt3( vol, sm_par);
    end
    
    surf_file = '';
    cortex_file = '';
    cortex = InCreateSurf( mat, vol, isoparm);
       
    return;
end


%% create brain surface
% unzip volume if it's in .nii.gz format
home_dir = getenv( 'HOME');
if strcmp( vol( 1:2), '~/')
    vol = [home_dir, filesep, vol(3:end)];
end

[p, nam, ext] = GetFileExtension( vol);
% if strcmpi( ext, 'nii.gz')
if ~isempty( strfind( ext, 'nii.gz'))
    % remove only .nii.gz from extension
    ext = strrep( ext, '.gz', '');
    nii_vol = fullfile( p, [nam, '.', ext]);
    if ~exist( nii_vol, 'file')
        fprintf( 'Unzipping %s to .nii file \n', vol);
        tmp = gunzip( vol);
        vol = tmp{1};
        
    else
        fprintf( 'Uncompressed version of %s was found as %s\n', vol, nii_vol);
        vol = nii_vol;
    end
    
% elseif strcmpi( ext, 'nii')
elseif ~isempty( strfind( ext, 'nii'))
    % do nothing
    
else
    error( '%s data format is not supported. Work with .nii and .nii.gz only.\n');
end

[wkpath, fname, ext] = fileparts( vol);
volname = [fname, ext];
gm_name = ['c1', volname];
wm_name = ['c2', volname];

% directory to save surface and cortex rendering results
savedir = fullfile( wkpath, options.fold);
if ~exist( savedir, 'dir')
    mkdir( savedir);
end

surf_name = [fname, '_surface'];
surf_file = fullfile( savedir, [surf_name, ext]);
surf_bain_file = fullfile( savedir, [surf_name, '_brain', ext]);
cortex_name = [fname, '_cortex'];
cortex_file = fullfile( savedir, [cortex_name, '.mat']);

usm_gm_file = fullfile( wkpath, gm_name);
usm_wm_file = fullfile( wkpath, wm_name);

if ~exist( usm_gm_file, 'file') || ~exist( usm_wm_file, 'file')    
    usm_gm_file = fullfile( wkpath, 'SPMSegment', gm_name);
    usm_wm_file = fullfile( wkpath, 'SPMSegment', wm_name);

    if ~exist( usm_gm_file, 'file') || ~exist( usm_wm_file, 'file')
        % error( 'No segmentation files for %s were found.', vol);
        % run the segmention on fly
        fprintf( '  No segmentation files were found. Running SPM segmentation ...\n');
        spm_segment_job( vol);
    end
    
else
    if strcmpi( options.re_run, 'yes')
        % this apparantly to be unnecessary
        fprintf( '  Re-Running SPM segmentation ...\n');
        spm_segment_job( vol);
        
    else
        % no re-run
        fprintf( '  Use existing gray matter tissue file %s\n', usm_gm_file);
        fprintf( '  Use existing white matter tissue file %s\n', usm_wm_file);
    end
end

% spatial smoothing
% use spm_select to get image file 
if smooth_kernel > 0
    fprintf( '  Smoothing kernel %d \n', smooth_kernel);
    smooth_kernel = [1, 1, 1] * smooth_kernel;
    fprintf( '  Smoothing gray matter image ... \n');
    spm_smooth( usm_gm_file, fullfile( savedir, ['sm_', gm_name]), smooth_kernel);
    gm_file = fullfile( savedir, ['sm_', gm_name]);
    
    fprintf( '  Smoothing white matter image ... \n');
    spm_smooth( usm_wm_file, fullfile( savedir, ['sm_', wm_name]), smooth_kernel);    
    wm_file = fullfile( savedir, ['sm_', wm_name]);

else
    fprintf( '  Unsmoothed images will be used.\n');
    gm_file = usm_gm_file;
    wm_file = usm_wm_file;
end

% add gray and white matter together
g = spm_read_vols( spm_vol( gm_file));
vol_info = spm_vol( wm_file);
w = spm_read_vols( vol_info);

%identifies "enclosed points" for later removal
a= (g+w) > threshold; 
% a = hollow_a(a);
d1 = size( a, 1);
d2 = size( a, 2);
d3 = size( a, 3);
brain = a( 1 : (d1-2), 2:(d2-1), 2:(d3-1)).*...
    a( 2:(d1-1), 1:(d2-2), 2:(d3-1)).*...
    a( 2:(d1-1), 2:(d2-1), 1:(d3-2)).*...
    a( 3:(d1), 2:(d2-1), 2:(d3-1)).*...
    a( 2:(d1-1), 3:(d2), 2:(d3-1)).*...
    a( 2:(d1-1), 2:(d2-1), 3:(d3));
%fill edges back in
b = cat(1, zeros(1, d2-2, d3-2), brain, zeros(1, d2-2, d3-2));
b = cat(2, zeros(d1, 1, d3-2), b, zeros(d1, 1, d3-2));
b = cat(3, zeros(d1, d2, 1), b, zeros(d1, d2, 1));
%remove enclosed points
a = a-b;
a(a<0) =0 ;
bwa = bwlabeln(a);

asize = length(a(:));
size4surface = [asize/10 asize/500]; %set required size for surface
for k=1:max(max(max(bwa)))
    if length(find(bwa==k))<size4surface(1) && length(find(bwa==k))>size4surface(2)
        a = bwa==k;
        fprintf( 'nice surface found \n');
        break;
    end
end
if max( max( max( a)))>1
    fprintf('no good surface representation found, change size of surface in get_mask\n');
end

% % file exists, use a different name
% if exist( surf_file, 'file')
%     surfix = datestr( clock);
%     if any( strfind( surfix, ':'))
%         surfix = strrep( surfix, ':', '-');
%     end
%     if any( strfind( surfix, ' '))
%         surfix = strrep( surfix, ' ', '-');
%     end
% 
%     surf_file = fullfile( savedir, [surf_name, '_', surfix, ext]);
%     surf_bain_file = fullfile( savedir, [surf_name, '_brain_', surfix, ext]);
% end

% disp( strcat( ['saving ' surf_file]));
vol_info.fname = surf_file;
spm_write_vol( vol_info, a);

img = a + b;
img = imfill( img, 'holes');
vol_info.fname = surf_bain_file;
spm_write_vol( vol_info, img);


%% plot images to check the surface mask
nbimgs = 30;
val = squeeze( mean( mean( b, 2)));
nz_val = find( val);
nbslices = length( nz_val);
step = floor( nbslices / (nbimgs-1));
slice2plot = 1 : step : nbslices;
loc = min( [nbimgs, length( slice2plot)]);
slice2plot = slice2plot( 1 : loc);
voldata = spm_read_vols( spm_vol( vol));
figure;
cnt = 0;
for k = 1 : length( slice2plot)
    cur_slice = slice2plot( k);
    cnt = cnt + 1;
    subplot( 5, 6, cnt)
    imagesc( voldata( :, :, nz_val( cur_slice)));
    colormap( gray);
    hold on; 
    [B, ~, N] = bwboundaries( a( :, :, nz_val( cur_slice)));
    for bound_idx = 1:N
       boundary = B{bound_idx};
       plot( boundary(:,2), boundary(:,1), 'y', 'LineWidth', 1);
    end
    axis equal
    title( ['slice ', num2str( cur_slice)]);
end

h = gcf;
set( h, 'Units', 'normalized', 'position', [0, 0, 0.8, 0.8], 'PaperPositionMode', 'auto', 'Units',' inches');
srcn_pos = get( h, 'Position');
set( h, 'PaperPosition', [0, 0, srcn_pos(3:4)], 'PaperSize', [srcn_pos(3:4)]);
print( h, '-dpdf', [surf_file, '.pdf']);
close


%% cortex rendering
% load grey and white segmentations
g = spm_read_vols( spm_vol( usm_gm_file));
w = spm_read_vols( spm_vol( usm_wm_file));
a = g+w;
if strcmpi( options.render_on_raw, 'yes')    
    brain_info = spm_vol( vol);
    img = spm_read_vols( brain_info);
    a = ( a > 0) .* img;
    
else
    brain_info = spm_vol( usm_wm_file); 
end

% smooth
if sm_par > 0
    a = imgaussfilt3( a, sm_par);
end

cortex = InCreateSurf( brain_info.mat, a, isoparm);

% % file exists, use a different name
% if exist( cortex_file, 'file')
%     if ~exist( 'surfix', 'var')
%         surfix = datestr( clock);
%         if any( strfind( surfix, ':'))
%             surfix = strrep( surfix, ':', '-');
%         end
%         if any( strfind( surfix, ' '))
%             surfix = strrep( surfix, ' ', '-');
%         end
%     end
% 
%     cortex_file = fullfile( savedir, [cortex_name, '_', surfix, '.mat']);
% end

% fprintf( 'saving %s\n', cortex_file);
save( cortex_file, 'cortex');

%% check results
fprintf( ['Go to ', surf_file, '.pdf to check surface mask.\n',...
    'Load ', cortex_file, ' to check the cortex rendering.\n'])


end % function

%% sub-functions
function cortex = InCreateSurf( mat, a, isoparm)
% 
st = sum( mat( :, [2 1 3]).^2) .^ 0.5;
sz = size( a);
fv = isosurface([1:sz(2)].* st(1), [1:sz(1)].* st(2), [1:sz(3)].* st(3), a, isoparm); %generate surface that is properly expanded for tesselation later fix indices to be in proper coordinates
vert=fv.vertices;
tri=fv.faces;

%reordering, etc
vert = bsxfun( @rdivide, vert, st);

cx=vert(:,1);
vert(:,1) = vert(:,2);
vert(:,2)=cx;
vert = vert*mat(1:3,1:3)'+repmat(mat(1:3,4)',size(vert,1),1);

cortex = [];
cortex.vert=vert;
cortex.tri=tri; %familiar nomenclature

end