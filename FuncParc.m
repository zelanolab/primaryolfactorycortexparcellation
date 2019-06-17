function FuncParc( wkpath, subjs, savedir, param)
% functional connectivity based parcellation

%% parameters setup
% % % data directory
% % wkpath = '/Volumes/Untitled/RestOlfAnaly/SubjectData';
% % savedir = '/Volumes/Untitled/RestOlfAnaly/Results/Test_Parc25Subjs_Olf_CZ';
% % 
% % % subjects
% % subjs = {'ro_subj001', 'ro_subj002'};
% % 
% % param = [];
% % % relative full path to 4d fmri volume in each subject's folder
% % param.volname = 'Nasal10Min/Nasal10Min.feat/filtered_func_data_rlt_nui_nuicorr_volsorm_filt2std_s3';
% % 
% % % whole-brain gray mask
% % param.brainmask = '/Volumes/Untitled/RestOlfAnaly/SubjectDataScript/ROIs/mask_avg152T1_thr100.nii.gz';
% % 
% % % Region-of-interests
% % param.roidir = '/Volumes/Untitled/RestOlfAnaly/SubjectDataScript/ROIs';
% % param.roinames = {'L_Olf_CZ', 'R_Olf_CZ'};
% % 
% % % correlation type
% % param.corr_type = 'pearson';
% % 
% % % number of clusters
% % param.clust_num = 2 : 10; % number of clusters for k-means clustering
% % 
% % % when doing parcellation for one ROI, whether to exclude other in the list
% % % or not
% % param.exclude_other = 'yes';
% % 
% % % 1 x N cell array of time windows, in volumes, use [] to indicate
% % % whole-time window, e.g. {[], [1, 540], [541, 1075], [1076, 1611]};
% % param.timewindow = {[]};
% % % time window names
% % param.win_names = {'All'};
% % 
% % % save individual functional connectivity matrix
% % param.save_indiv_mat = 'no';
% % % save group-average functional connectivity matrix
% % param.save_group_mat = 'yes';
% % 
% % % Matlab clustering distance method
% % param.dist_method = 'correlation';
% % 


%% TODO. Validate input arguments
if ~exist( 'MRIread', 'file')
    error( 'Freesurfer matlab toolbox is not installed.');
end


%%  
param.wkpath = wkpath;
param.subjs = subjs;
param.savedir = savedir;

% read whole-brain mask
brainmask = MRIread( param.brainmask);

% read all ROIs
nbrois = length( param.roinames);
roimasks = cell( nbrois, 1);
for roi_idx = 1 : length( param.roinames)
    roimask_file = fullfile( param.roidir, param.roinames{ roi_idx});
    roimasks{ roi_idx} = MRIread( roimask_file);
end
        
switch param.exclude_other
    case 'yes'
        % read all ROIs and exclude them from the whole-brain mask
        fprintf( 'Excluding all ROIs from brain mask.\n');
        all_roimask = cellfun( @(x) x.vol, roimasks, 'UniformOutput', false);
        all_roimask = cat( 4, all_roimask{:});
        brainmask.vol( sum( all_roimask, 4) > 0) = 0;        
    case 'no'
        % do nothing
    otherwise
        error('Unknown option.');
end

for roi_idx = 1 : length( param.roinames)
    param.roiname = param.roinames{ roi_idx};

    % exclude current roi
    current_brainmask = brainmask.vol;
    current_brainmask( roimasks{ roi_idx}.vol == 1) = 0;
    nb_brainvox = sum( current_brainmask( :) == 1);
       
    nb_roivox = sum( roimasks{ roi_idx}.vol( :) == 1);
    for win_idx = 1 : length( param.timewindow)
        if isempty( param.timewindow{ win_idx})
            vol_ind = [];
        else
            vol_ind = param.timewindow{ win_idx}( 1) : param.timewindow{ win_idx}( 2);
        end
        
        % resulting correlation matrix
        rval_roi2brain = zeros( nb_roivox, nb_brainvox);
        rval_roi2roi = zeros( nb_roivox, nb_roivox);
        for subj_idx = 1 : length( subjs)
            subj = subjs{ subj_idx};
            fprintf( 'Calculating connectivity matrix for ROI: %s, subject: %s\n', param.roinames{ roi_idx}, subj);
            fmri = fullfile( wkpath, subj, param.volname);

            img = MRIread( fmri);
            nbvols = size( img.vol, 4);
            tmp_mask = repmat( roimasks{ roi_idx}.vol, [1, 1, 1, nbvols]);
            roi_ts = img.vol( tmp_mask == 1);
            roi_ts = reshape( roi_ts, [nb_roivox, nbvols]);

            tmp_mask = repmat( current_brainmask, [1, 1, 1, nbvols]);
            brain_ts = img.vol( tmp_mask == 1);
            brain_ts = reshape( brain_ts, [nb_brainvox, nbvols]);

            if isempty( vol_ind) 
                vol_ind = 1 : size( img.vol, 4);
            end
            
            rval_mat = corr( roi_ts( :, vol_ind)', brain_ts( :, vol_ind)', 'type', param.corr_type);
            % save individual functional connectivity matrix
            if strcmpi( param.save_indiv_mat, 'yes')
                subj_conn_savedir = fullfile( savedir, param.win_names{ win_idx}, 'IndividualConnMatrix', roinames{ roi_idx});
                if ~exist( subj_conn_savedir, 'dir')
                    mkdir( subj_conn_savedir);
                end
                save( fullfile( subj_conn_savedir, subjs{ subj_idx}), 'rval_mat', 'param', '-v7.3');
            end
            
            rval_roi2brain = rval_roi2brain + G_r2z( rval_mat);
            rval_roi2roi = rval_roi2roi + G_r2z( corr( roi_ts( :, vol_ind)', roi_ts( :, vol_ind)', 'type', param.corr_type));
        end

        % group-average
        rval_roi2brain = rval_roi2brain / length( subjs);
        rval_roi2roi = rval_roi2roi / length( subjs);
        
        % convert back to r from Fisher's z
        rval_roi2brain = tanh( rval_roi2brain);        
        rval_roi2roi = tanh( rval_roi2roi);

        % save group average connectivity matrix
        if strcmpi( param.save_group_mat, 'yes')
            tmp_savedir = fullfile( savedir, param.win_names{ win_idx}, 'GroupConnMatrix');
            if ~exist( tmp_savedir, 'dir')
                mkdir( tmp_savedir);
            end
            save( fullfile( tmp_savedir, param.roinames{ roi_idx}), 'rval_roi2brain', 'rval_roi2roi', 'param', '-v7.3');
        end
        
        % different time windows
        rval_roi2brain( :, any( isnan( rval_roi2brain), 1)) = [];

        tmp_savedir = fullfile( savedir, param.win_names{ win_idx}, 'Parc');
        if ~exist( tmp_savedir, 'dir')
            mkdir( tmp_savedir);
        end
        
        for clust_idx = param.clust_num
            % Do kmeans with k clusters
            idx = kmeans( rval_roi2brain, clust_idx, 'distance', param.dist_method);
            tmp = brainmask;
            tmp.vol = 0 * tmp.vol;
            tmp.vol( roimasks{ roi_idx}.vol ~= 0) = idx;
            MRIwrite( tmp, fullfile( tmp_savedir, ['CC_', num2str( clust_idx), '_', param.roinames{ roi_idx}, '.nii.gz']), 'float');
        end
        
    end % time window    
end % roi loop

end % function
