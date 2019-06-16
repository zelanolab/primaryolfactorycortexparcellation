function [h, ch] = PlotSurf3( ax, vertex, faces, varargin)
% 
% see http://stackoverflow.com/questions/30921003/matlab-how-to-make-camera-light-follow-3d-rotation
% 
% vertex_coords, vertex coordinates, one vertex per row
% faces, polygons, one row per polygon
% label_color, label_vertex, cell array of vertex index
% label_vertex, label color, {Mx3, Nx3, ...}
% 

options = struct( 'fc', [1, 1, 1] * 0.5,... % [0.8, 0.8, 0.8],... % background face color
    'label_color', [0.3, 0.5, 0.7],... 
    'label_vertex', [],...
    'view_angle', [-90, 15],...
    'coords', [],... % Nx3 or Nx4
    'marker_size', 8,...
    'marker_color', [],...
    'sfile', '',...
    'hold', 'off'); 

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

if isempty( ax) || ~ishandle( ax)
    ax = gca;
end

set( get( ax, 'Parent'), 'color', 'k');

% face color
fc = options.fc;

label_vertex = options.label_vertex;
if ~isempty( label_vertex)    
    if ~iscell( label_vertex)
        label_vertex = {label_vertex};
    end

    nblabels = length( label_vertex);    
    label_color = options.label_color;
    if ~iscell( label_color)
        label_color = repmat( {label_color}, [nblabels, 1]);
        
    elseif nblabels ~= length( label_color)
        error( 'data length dismatch');
    end
    
    for k = 1 : nblabels
        cur_c = label_color{ k};
        cur_nb_vertex = length( label_vertex{k});
        nb_c = size( cur_c, 1);
        if nb_c ~= cur_nb_vertex && nb_c == 1
            cur_c = repmat( cur_c( 1, :), [cur_nb_vertex, 1]);
            
        elseif nb_c == cur_nb_vertex
            % do nothing
            
        else
            error( 'dimension dismatch');
        end
        fvdata( label_vertex{k}, :)  = cur_c;
    end
end

if iscell( vertex)
    % do nothing
else
    vertex = {vertex};
    faces = {faces};
end

if size( fc, 1) == 1
    fc = repmat( fc, [length( vertex), 1]);
end


for vidx = 1 : length( vertex)
    nb_vertex = size( vertex{ vidx}, 1);
    fvdata = repmat( fc( vidx, :), [nb_vertex, 1]);
%     h( vidx) = patch( ax, 'Vertices', vertex{ vidx}, 'Faces', faces{ vidx},...
%         'FaceVertexCData', fvdata, 'EdgeColor', 'none',...
%         'AmbientStrength', 0.05, 'facealpha', 1, 'CDataMapping', 'direct', 'facecolor', 'interp', 'FaceLighting','flat');
%     material dull;
      h( vidx) = patch( ax, 'Vertices', vertex{ vidx}, 'Faces', faces{ vidx},...
        'FaceVertexCData', fvdata, 'EdgeColor', 'none',...
        'AmbientStrength', 0.05, 'FaceVertexAlphaData', rand( length( vertex{ vidx}), 1),...
        'CDataMapping', 'direct', 'facecolor', 'interp', 'FaceLighting','flat');
        material dull;
end


% 00plot markers at specified coordinates
if ~isempty( options.coords)
    coords = options.coords;
    if size( coords, 2) == 3
        coords( :, 4) = 1;
    end
    
    grp = unique( coords( :, 4));
    nbgroups = length( grp);
    if isempty( options.marker_color)
        cks = {'b', 'r',  'm', 'c', 'k', 'g',  'y'};
    else
        cks = options.marker_color;
    end    
    
    hold( ax, 'on');
    for k = 1 : nbgroups 
        g = coords( :, 4) == grp( k);
        if k > length( cks)
            cks{k} = rand( 1, 3);
        end
        
        plot3( ax, coords( g, 1), coords( g, 2), coords( g, 3),...
            'o', 'color', cks{k}, 'markerfacecolor', cks{k}, 'markersize', options.marker_size);
    end
    hold( ax, 'off');
end

if strcmpi( options.hold, 'off')    
    axis( ax, 'off');
    axis( ax, 'equal');
    axis( ax,  'vis3d');
    lighting( ax, 'gouraud');   
    view( ax, options.view_angle);
    set( ax, 'LooseInset', get(gca, 'TightInset'), 'box', 'off');
    c = camlight( 'headlight');
    set( c, 'style', 'infinite'); 
    ch = rotate3d( ax);
    ch.ActionPostCallback = @RotationCallback;
    ch.Enable = 'on';  
end

% Sub function for callback
    function RotationCallback( ~, ~)
        c = camlight( c, 'headlight');
    end

end % function
