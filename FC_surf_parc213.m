clear all;clc
%% 

dir_data = '/home/izaac/Documents/Aphasia';
dir_sub = dir(fullfile(dir_data,'BA*'));
list_subname = {dir_sub(:).name};

load('/home/izaac/Documents/Aphasia_result/WAB_percent_improve/pre.mat')
load('/home/izaac/Documents/Aphasia_result/WAB_percent_improve/delta.mat')
load('/home/izaac/Documents/Aphasia_result/seed_index_fs6.mat')
load('/home/izaac/Documents/Aphasia_result/lh.mat')
load('/home/izaac/Documents/Aphasia_result/WAB_percent_improve/proportional.mat')

% seed_lh = [2213];
% seed_rh = [];
% 
% if isempty(seed_rh)
%     seed_ind = seed_lh + 1;
% else
%     seed_ind = [seed_lh,seed_rh + 40962] + 1;
% end

% seed_ind = [];
seed_ind = seed_ind + 1;


% fcmap_path = fullfile(dir_data,dir_ses.name,'fcmap');
% mkdir(fcmap_path);
% 
% fcmapW_path = fullfile(dir_data,dir_ses.name,'fcmap_W');
% mkdir(fcmapW_path);

% dir_par = fullfile(dir_data,char(subname),char(sesname),'Parcellation-92');
[lh_vertex,lh_label,lh_ct] = read_annotation('/home/izaac/Documents/code/template/3.6/213/lh_parc_result.annot');
[rh_vertex,rh_label,rh_ct] = read_annotation('/home/izaac/Documents/code/template/3.6/213/rh_parc_result.annot');
% rh_label = rh_label + 1.5;
% rh_ct.table(:,5) = rh_ct.table(:,5) + 1.5;
% label = [lh_label;rh_label];
% table = [lh_ct.table(:,5);rh_ct.table(:,5)];

lh_label = load_mgh('/home/izaac/Documents/code/template/3.6/213/lh.Clustering_108_fs6.mgh');
rh_label = load_mgh('/home/izaac/Documents/code/template/3.6/213/rh.Clustering_108_fs6.mgh');
rh_label = rh_label + 108;
rh_label(find(rh_label == 108)) = 0;
label = [lh_label;rh_label];
table = 1:216;
table = table';
%% 

for nsub = 1:numel(list_subname)
% for nsub = lh'
% for nsub = 30
    subname = list_subname(nsub);
    dir_ses = dir(fullfile([dir_data,'/',char(subname)],'*rest0*'));
    list_sesname = {dir_ses(:).name};
    for nses = 1:numel(list_sesname)
%     for nses = 2
        sesname = list_sesname(nses);
        
        dir_par = fullfile(dir_data,char(subname),char(sesname),'Parcellation-92');
        [lh_vertex,lh_label,lh_ct] = read_annotation([dir_par,'/lh_parc92_fs6_surf.annot']);
        [rh_vertex,rh_label,rh_ct] = read_annotation([dir_par,'/rh_parc92_fs6_surf.annot']);
        rh_label = rh_label + 1;
        rh_ct.table(:,5) = rh_ct.table(:,5) + 1;
        label = [lh_label;rh_label];
        table = [lh_ct.table(:,5);rh_ct.table(:,5)];
        
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
        
        parc218_run = [];
        
        for nnet = 1:numel(table)
            parc218_run(nnet,:) = mean(run_all(find(label == table(nnet)),:));
        end
        
        parc_run = parc218_run;
        mask213 = 1:226;
        mask213 = mask213';
        mask213(isnan(parc_run(:,1)),:) = [];
        parc_run(isnan(parc_run(:,1)),:) = [];
        
        %% compute correlation
        for seedi = 1:numel(seed_ind)
            data_seed = run_all(seed_ind(nsub),:);
            fcmap_all(:,1)=corr(data_seed',parc_run');     %   这部分变量需检查
            fcmap_all(isnan(fcmap_all)==1) = 0;
%         end
        
        fcmap_all = Fisherz(fcmap_all);
%         fcmap_all(find(fcmap_all>0)) = 0;     %取正或者负的部分
%         pdata_Nmn = (pdata-min(pdata))./(max(pdata)-min(pdata));     %极值归一
%         delta_z = atan(delta)*(2/pi);      %正切归一
%         delta_z = zscore(delta);     %z变换归一
%         pre_z = zscore(pre);     %z变换归一
%         W = delta./(100-pre);
%         W_z = zscore(W);
%         fcmap_all_W = fcmap_all .* W_z(nsub);
%         
%         fcmap_parc218 = zeros(218,1);
%         fcmap_parc218(mask213) = fcmap_all;
%         fcmap_213 = [];
%         for nv = 1:81924
%             if label(nv) == 0 || label(nv) == 1
%                 fcmap_213(nv,1) = 0;
%             else
%                 fcmap_213(nv,1) = fcmap_parc218(find(table == label(nv)),1);
%             end
%         end
            
        
        num_vertex = size(fcmap_all,1)/2;
%         for seedi = 1:numel(seed_ind)
            seed = seed_ind(nsub);
            if seed<40962
                fname = ['seed_lh',num2str(seed-1)];
            else
                fname = ['seed_rh',num2str(seed-40962-1)];
            end
            fcmap_path = [dir_data,'/',char(subname),'/',char(sesname)];
            mkdir(fullfile(fcmap_path,fname,'/fcmap_213_template'));
            
            save([fcmap_path,'/',fname,'/fcmap_213_template/fcmap213.mat'],'fcmap_all')
%             save_mgh(fcmap_all(1:118,1),fullfile(fcmap_path,fname,'/fcmap_213/lh_fcmap.mgh'),eye(4));
%             save_mgh(fcmap_all(119:end,1),fullfile(fcmap_path,fname,'/fcmap_213/rh_fcmap.mgh'),eye(4));
              
%             fcmapW_path = [dir_data,'/',char(subname),'/',char(sesname)];
%             mkdir(fullfile(fcmapW_path,fname,'/fcmap_potentialZW'));
%             
%             save_mgh(fcmap_all_W(1:num_vertex,1),fullfile(fcmapW_path,fname,'/fcmap_potentialZW/lh_fcmap.mgh'),eye(4));
%             save_mgh(fcmap_all_W(num_vertex+1:num_vertex*2,1),fullfile(fcmapW_path,fname,'/fcmap_potentialZW/rh_fcmap.mgh'),eye(4));
%             save_mgh(fcmap_213(1:40962,1),fullfile(fcmap_path,fname,'/fcmap_213/lh_fcmap.mgh'),eye(4));
%             save_mgh(fcmap_213(40963:81924,1),fullfile(fcmap_path,fname,'/fcmap_213/rh_fcmap.mgh'),eye(4));
%         end
        end
    end
end
