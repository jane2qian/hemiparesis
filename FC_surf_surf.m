clear all;clc

dir_data = '/home/izaac/Documents/NG/code_club/FC/Aphasia';
dir_sub = dir(fullfile(dir_data,'BA*'));
list_subname = {dir_sub(:).name};


seed_lh = [26890];
seed_rh = [];

if isempty(seed_rh)
    seed_ind = seed_lh + 1;
else
    seed_ind = [seed_lh,seed_rh + 40962] + 1;
end

% fcmap_path = fullfile(dir_data,dir_ses.name,'fcmap');
% mkdir(fcmap_path);

for nsub = 1:numel(list_subname)
% for nsub = 1:2
    subname = list_subname(nsub);
    dir_ses = dir(fullfile([dir_data,'/',char(subname)],'*rest*'));
    list_sesname = {dir_ses(:).name};
    %     for nses = 1:numel(list_sesname)
    for nses = 1:2
        sesname = list_sesname(nses);
        dir_surf = fullfile(dir_data,char(subname),char(sesname),'Preprocess','surf');
        
        %% read surf data
        name_runs = dir(fullfile(dir_surf,'*sm6.nii.gz'));
        name_runs = {name_runs(:).name}';
        for runi = 1:numel(name_runs)
            name_runs{runi} = name_runs{runi}(4:end);
        end
        name_runs = unique(name_runs);
        
        run_all = [];
        for runi = 1:numel(name_runs)
            name_subrun = dir(fullfile(dir_surf,['*h*',name_runs{runi}]));
            name_subrun = {name_subrun(:).name}';
            run_tmp = [];
            for hemi = 1:numel(name_subrun)
                name_subrun{hemi}
                hdr = MRIread(fullfile(dir_surf,name_subrun{hemi}));
                subrun_tmp = single(hdr.vol);
                subrun_tmp = reshape(subrun_tmp,[hdr.nvoxels,hdr.nframes]);
                run_tmp = cat(1,run_tmp,subrun_tmp);
                clear subrun_tmp hdr
            end
            run_all = cat(2,run_all,run_tmp);
        end
        
        %% compute correlation
        for seedi = 1:numel(seed_ind)
            data_seed = run_all(seed_ind(seedi),:);
            fcmap_all(:,seedi)=corr(data_seed',run_all');
            fcmap_all(isnan(fcmap_all)==1) = 0;
        end
        
        num_vertex = size(fcmap_all,1)/2;
        for seedi = 1:numel(seed_ind)
            seed = seed_ind(seedi);
            if seed<40962
                fname = ['seed_lh',num2str(seed-1)];
            else
                fname = ['seed_rh',num2str(seed-40962-1)];
            end
            fcmap_path = [dir_data,'/',char(subname),'/',char(sesname)];
            mkdir(fullfile(fcmap_path,'/seeds/',fname));
            save_mgh(fcmap_all(1:num_vertex,seedi),fullfile(fcmap_path,'/seeds/',fname,'/lh_fcmap.mgh'),eye(4));
            save_mgh(fcmap_all(num_vertex+1:num_vertex*2,seedi),fullfile(fcmap_path,'/seeds/',fname,'/rh_fcmap.mgh'),eye(4));
            
            %截图
            data_path = fullfile(fcmap_path,'/seeds/',fname);
            cmd = ['~/Documents/code/snapshot/snapshot_fcmap_betch.csh ', data_path,' fsaverage6 inflated 0.2 0.6'];
            system(cmd)
            snapshot(data_path,'fc_map')

        end
    end
end
