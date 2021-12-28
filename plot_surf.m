path = '/home/eyre/project101/analysis02/seed_groupFC/tbs/ctbs/pre';
dir_sub = dir(fullfile([path,'/net7*']));
list_sub = {dir_sub(:).name};
%% 

for nsub = 1:numel(list_sub)
    sub = char(list_sub(nsub))
    roipath = dir(fullfile([path,'/',sub]))
    path_data = roipath.folder%[path,'/',sub,'/seed*/roi_map'];
%     cmd = ['/bin/csh /home/eyre/pic-generate/snapshot_fcmap_betch.csh ', path_data,' fsaverage6 inflated 0.2 0.6'];
    % system(cmd)
    snapshot(path_data,'rmap')
end