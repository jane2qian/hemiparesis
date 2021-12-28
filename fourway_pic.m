path = '/home/eyre/project101/urge/mydata/seed_pre_post';
dir_sub = dir(fullfile([path,'/t*']));
list_sub = {dir_sub(:).name};
for nsub = 1:numel(list_sub)
    sub = char(list_sub(nsub));
    path_data = [path,'/',sub];
    cmd = ['/home/eyre/project101/snapshot_fcmap_betch.csh ', path_data,' fsaverage6 inflated 1.96 5'];
    system(cmd)
    snapshot(path_data,'tmap')
end