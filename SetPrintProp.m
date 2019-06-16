function h = SetPrintProp( h, width, height, orig, varargin)
% Set figure properties for pdf printing
% 
% Usage
%   SetPrintProp( [], width, height);
%   SetPrintProp( gcf, width, height);
% 
% Input
%   h, handle to figure, e.g. gcf (default), or its children object
%   width, (0, 1], width of the figure (in normalized unit), default 0.5.
%   height, (0, 1], height of the figure (in normalized unit), default 0.8.
%   Optional
%       orig, left-bottom position (in normalized unit) 
% 
% Output
%   h, handle to figure
% 
% print( h, '-dpdf', '-fillpage', ['AM1-Zmap.pdf']);


% default current figure handle
if isempty( h)
    h = gcf;
end

try 
    if ~strcmpi( get( h, 'Type'), 'figure')
        try 
            h = get( h, 'parent');

        catch
            error( 'h must be or handle to figure or its children object.');
        end
    end
    
catch
    error( 'h must be or handle to figure or its children object.');
end

% set default values
if nargin == 1
    width = 0.8;
    height = 0.5;
    
elseif nargin == 2
    height = 0.8;
    
end

if ~isscalar( width) || ~isscalar( height) || width <= 0 || width > 1 || height <= 0 || height > 1
    error( 'Width and Height must be scalars between 0 and 1.');
end

% left-bottom point
if ~exist( 'orig', 'var') || isempty( orig)
    orig = 0.5 * (1 - [width, height]);
    
else
    if ~isvector( orig) || length( orig) ~= 2 || any( orig < 0 ) || any( orig >= 1)
        error( 'invalid left-bottom origs.')
    end
   
    if orig( 1) > width
        orig( 1) = 0;
    end
    if orig( 2) > height
        orig( 2) = 0;
    end
end

set( h, 'Units', 'normalized', 'position', [orig(1), orig(2), width, height]);
set( h, 'PaperPositionMode', 'auto', 'Units',' inches');
srcn_pos = get( h, 'Position');
set( h, 'PaperPosition', srcn_pos, 'PaperSize', [srcn_pos( 3), srcn_pos( 4)]);

% save pdf file.
if ~isempty( varargin) > 0
    if isdir( varargin{1})
        error( 'a file name must be provided.');
    end
    
    file2save = varargin{1};
    [p, nam] = fileparts( file2save);
    if ~exist( p, 'dir')
        mkdir( p);
    end
        
    if exist( file2save, 'file')
        nam = [nam, GetClockSurfix, '.pdf'];
    else
        nam = [nam, '.pdf'];
    end
    file2save = fullfile( p, nam);
    print( h, '-dpdf', '-fillpage',  file2save);
end


end % function