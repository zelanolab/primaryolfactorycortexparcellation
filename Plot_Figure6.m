%% make piriform 3d surface model

% parcellation result
% a 3d model will be created for each subregion
pcfile = '~/Documents/RestfMRI/SourceData/Figure_6_sourcedata/Figure_2_Parcellation_2mm.nii.gz';

% FSL brain surface 
fsl_surf = load( '~/Desktop/Data/tmp/fsl_surf.mat');
fsl_surf = fsl_surf.fsl_surf;


pc = MRIread( pcfile);
pc.vol = permute( pc.vol, [2, 1, 3]);
trasform = pc.vox2ras1;

% affects the appearance of the 3d models
sm_sigma = 1;
isoparm = 0.2;

cortex = cell( 4, 1);
clust_center = zeros( 4, 3);
for struct_ind = 1 : 4
    % get specified brain structure
    img = abs( pc.vol - struct_ind) < eps;
    img = double( img);
    img = pc.vol .* img;
    img = img ./ max( img(:));

    [~, ~, cortex{ struct_ind, 1}] = CreateSurf( img, 'rend_sm', sm_sigma, 'isoparm', isoparm, 'transform', trasform);
 
    N = floor( size( img, 1)/2);
    tmp = img;
    tmp( N+1:end, :, :) = 0;
    [~, ~, cortex{ struct_ind, 2}] = CreateSurf( tmp, 'rend_sm', sm_sigma, 'isoparm', isoparm, 'transform', trasform);
    
    % cluster center
    clust_center( struct_ind, :) = mean( cortex{ struct_ind, 2}.vert, 1);
    
    tmp = img;
    tmp( 1:N, :, :) = 0;
    [~, ~, cortex{ struct_ind, 3}] = CreateSurf( img, 'rend_sm', sm_sigma, 'isoparm', isoparm, 'transform', trasform);
end


% get specified brain structure
img = double( pc.vol > 0);
img = pc.vol .* img;
img = img ./ max( img(:));
tmp = img;
tmp( N+1:end, :, :) = 0;
[~, ~, merge_pc_cortex] = CreateSurf( tmp, 'rend_sm', sm_sigma, 'isoparm', isoparm, 'transform', trasform);
    
roi_center = mean( clust_center, 1);

atlas_file = '/usr/local/fsl/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr50-2mm.nii.gz';
atlas = MRIread( atlas_file);
atlas.vol = permute( atlas.vol, [2, 1, 3]);

tmp = double( atlas.vol == 19);
tmp( N+1:end, :, :) = 0;
[~, ~, hipp] = CreateSurf( tmp, 'rend_sm', sm_sigma, 'isoparm', isoparm, 'transform', trasform);

% central amygala, val = 12
amy_atlas_file = '/usr/local/fsl/data/atlases/Juelich/Juelich-maxprob-thr50-2mm.nii.gz';
amy_atlas = MRIread( amy_atlas_file);
amy_atlas.vol = permute( amy_atlas.vol, [2, 1, 3]);

tmp = double( amy_atlas.vol == 12);
tmp( N+1:end, :, :) = 0;
[~, ~, amy] = CreateSurf( tmp, 'rend_sm', sm_sigma, 'isoparm', isoparm, 'transform', trasform);


d = bsxfun( @minus, fsl_surf.brain.vert.rh, roi_center);
d = sqrt( sum( d .^ 2, 2));


sphere_vertex = find( d < 15);
loc = any( ismember( fsl_surf.brain.face.rh, sphere_vertex), 2);

dd = fsl_surf.brain.vert.rh( sphere_vertex, :);
dd_max = max( dd, [], 1);
dd_min = min( dd, [], 1);
vv = fsl_surf.brain.vert.rh( :, 3) < dd_max(3) & fsl_surf.brain.vert.rh( :, 3) > dd_min(3) & ...
    fsl_surf.brain.vert.rh( :, 2) < dd_max(2) & fsl_surf.brain.vert.rh( :, 2) > dd_min(2);

loc2 = any( ismember( fsl_surf.brain.face.rh, find(vv)), 2);
loc = loc | loc2;


% adjust the view angle, tansparency etc manually

figure;
SetPrintProp( gcf, 0.2, 0.3);
[h] = PlotSurf3( [], fsl_surf.brain.vert.rh, fsl_surf.brain.face.rh( ~loc, :), 'fc', [1, 1, 1] * 0.7);
[h2] = PlotSurf3( [], fsl_surf.brain.vert.rh, fsl_surf.brain.face.rh( loc, :), 'hold', 'on');
alpha( h, 0.5);
alpha( h2, 0.15);

marker_color = {[1 0 0];...
    [204, 51, 255]/255;...
    [0 0 1];...
    [0 1 0];...
    [0 1 1]};
for k = 1 : 4
    PlotSurf3( [], cortex{k, 2}.vert, cortex{k, 2}.tri, 'hold', 'on', 'fc', marker_color{k});
end
print( gcf, '-dpdf', '-fillpage', '~/Documents/RestfMRI/Figures/3dBrain_1.pdf');


for k = 1 : 4
    PlotSurf3( [], cortex{k, 2}.vert, cortex{k, 2}.tri, 'hold', 'on', 'fc', [.3 .5 .6]);
end
print( gcf, '-dpdf', '-fillpage', '~/Documents/RestfMRI/Figures/3dBrain_2.pdf');



% plot hippocampus and amygdala
figure;
SetPrintProp( gcf, 0.2, 0.3);
[h] = PlotSurf3( [], fsl_surf.brain.vert.rh, fsl_surf.brain.face.rh( ~loc, :));
[h2] = PlotSurf3( [], fsl_surf.brain.vert.rh, fsl_surf.brain.face.rh( loc, :), 'hold', 'on');
alpha( h, 0.3);
alpha( h2, 0.3);
for k = 1 : 4
    PlotSurf3( [], cortex{k, 2}.vert, cortex{k, 2}.tri, 'hold', 'on', 'fc', marker_color{k});
end

% hippocampus and maygdala 3d model
PlotSurf3( [], amy.vert, amy.tri, 'hold', 'on', 'fc', [103 197 211]/255);
PlotSurf3( [], hipp.vert, hipp.tri, 'hold', 'on', 'fc', [103 197 211]/255);

print( gcf, '-dpdf', '-fillpage', '~/Documents/RestfMRI/Figures/3dBrain_4.pdf');




