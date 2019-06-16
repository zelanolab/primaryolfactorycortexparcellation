function [redblue, red, blue] = CustomColormap( varargin)
% generate blue and red colormaps
% 
% Syntax
%   [redblue, red, blue] = CustomColormap();
%   [redblue, red, blue] = CustomColormap( varargin);
% 
% Input
%   Optional Key-value pairs
%       'N', length of colormap, default 200
%       'r1', a value between 0 and 1
%       'r2', a value between 0 and 1
%       'r3', a value between 0 and 1, 0 < r1 < r2 < r3 < 1
%             r1, r2 and r3 together controls the proportion of yellow and red (blue and light blue)
%       'clim', [lo, hi] colormap limits
%       'threshold', [] (default), used only when clim is set
%       'mapidx', 1 (default) | 2 | 3
% 
%   white-(white, when threshold is set)-r1(yellow)-r2(red-yellow)-r3(light red)-deep red
% 
% Output
%   red, red colormap
%   blue, blue colormap
%   redblue, red-blue colormap
% 
% Example
%   img = randn( 6, 10);
%   clim = [-2, 2];
%   
%   figure;
%   imagesc( img); colorbar; caxis( clim);
%   [red, blue, redblue] = CustomColormap();
%   colormap( redblue);
% 
%   figure;
%   imagesc( img); colorbar; caxis( clim);
%   [red, blue, redblue] = CustomColormap( 'clim', clim, 'threshold', 0.6, 'r1', 0.6, 'r2', 0.7, 'r3', 0.9);
%   colormap( redblue);
% 

if nargout < 1
    return;
end

% default values
options = struct( 'r1', [],...
    'r2', [],...
    'r3', [],...
    'clim', [],...
    'threshold', [],...
    'n', 200,...
    'mapidx', 1);

% sparse input arguments
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

clim = options.clim;
threshold = options.threshold;
r1 = options.r1;
r2 = options.r2;
r3 = options.r3;
N = options.n;

if ~isscalar( N)
    error( 'The length of colormap must be an integar number');
end
N = floor( N);

if ~isempty( r1) && ~isempty( r2) && ~isempty( r3)
    % r1, r2 and r3 must be all set
    if ~isscalar( r1) || ~isscalar( r2) || ~isscalar( r3)
        error( 'r1, r2, and r3 must be scalars.');

    elseif ~( r1 > 0 && r1 < 1) || ~( r2 > 0 && r2 < 1) || ~(r3 > 0 && r3 < 1) || r1 > r2 || r1 > r2 || r2 > r3
        error( 'r1, r2, and r3 must be numbers between 0 and 1, and r1 < r2 < r3');
    end
    
else
    % use default control points
    r1 = [];
    r2 = [];
    r3 = [];
end

if isempty( threshold)
    threshold_loc = 0;
    
else    
    if isempty( clim)
        error( '"clim" must be set when "threshold" is used.');
        
    elseif ~isvector( threshold) || ~isvector( clim)
        error( '"threshold" and "clim" must be vectors.');
        
    elseif length( threshold) ~= 1 || length( clim) ~= 2
        error( '"threshold" must be a scalar and clim must be a vector with a length of 2.');
    end
    
    threshold = abs( threshold);    
    if threshold >= min( abs( clim))
        error( '"threshold" must be smaller than minimal absolute value of "clim".');
    end
    
    if clim(1)*clim(2) < 0
        threshold_loc_blue = 1 + (threshold + clim(1)) / (-1*clim(1));
        threshold_loc_red = threshold / clim(2);
        
    else        
        if clim(1) < 0
            threshold_loc = 1 + (threshold + clim(1)) / diff( clim);
        else
            threshold_loc = (threshold - clim(1)) / diff( clim);
        end
    end
    
end

% change this part for different colormaps
if options.mapidx == 1
    rval = [1 1 0.6 0.2 0.02 0];
    gval = [1 1 1 0.35 0.1 0];
    bval = [1 1 0.99 0.97 0.9 0.55];

elseif options.mapidx == 2
    rval = [1 1 0.45 0 0 0];
    gval = [1 1 1 0.45 0 0];
    bval = [1 1 0.99 0.98 0.75 0.3];

elseif options.mapidx == 3
    rval = [1 1 0.45 0 0 0];
    gval = [1 1 1 0.45 0 0];
    bval = [1 1 0.99 0.98 0.5 0];
    
else
    error( 'Unknown map index.');
end

if exist( 'threshold_loc_blue', 'var')   
    n_blue = floor( N* abs( clim(1))/diff(clim));
    [r, g, b] = GetRGB( n_blue, threshold_loc_blue, rval, gval, bval, r1, r2, r3);
    blue = [r, g, b];
    
    n_red = ceil( N* clim(2)/diff(clim));
    [r, g, b] = GetRGB( n_red, threshold_loc_red, rval, gval, bval, r1, r2, r3);
    red = [b, g, r];
    
else    
    [r, g, b] = GetRGB( N, threshold_loc, rval, gval, bval, r1, r2, r3);    
    blue = [r, g, b];
    red = [b, g, r];
end

redblue = [flipud( blue( 2:end, :)); red];

end % main function




%% sub-function
function [r, g, b] = GetRGB( N, threshold, rval, gval, bval, r1, r2, r3)
    pnt = zeros( 1, 6);
    pnt( 6) = 1;
    pnt( 2) = threshold;
    step = 1 - threshold;
    if isempty( r1) || isempty( r2) || isempty( r3)
        pnt( 3:5) = pnt(2) + [1, 2, 3]*step/4;
    else       
        pnt( 3:5) = pnt(2) + [r1, r2, r3]*step;
    end

    loc = linspace( 0, 1, N);
    loc = loc';
    if threshold == 0
        r = interp1( pnt( :, 2:end), rval( :, 2:end), loc);
        g = interp1( pnt( :, 2:end), gval( :, 2:end), loc);
        b = interp1( pnt( :, 2:end), bval( :, 2:end), loc);
    else
        r = interp1( pnt, rval, loc);
        g = interp1( pnt, gval, loc);
        b = interp1( pnt, bval, loc);
    end
    
end
