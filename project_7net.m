clear all;clc

%% 

dir_par = fullfile('/home/eyre/project101/Yeo_FreeSurfer/fsaverage6/label/');

[lh_vertex,lh_label,lh_ct] = read_annotation([dir_par, 'lh.Yeo2011_7Networks_N1000.annot']);
[rh_vertex,rh_label,rh_ct] = read_annotation([dir_par, 'rh.Yeo2011_7Networks_N1000.annot']);

%% 

lh_surf = load_mgh(fullfile('/home/eyre/project101/analysis02/seed_groupFC/seedfc_pre/lh.avg_corr_across_fs6.mgh')); 
rh_surf = load_mgh(fullfile('/home/eyre/project101/analysis02/seed_groupFC/seedfc_post/rh.avg_corr_across_fs6.mgh')); 
%% 
lh_table = lh_ct.table(:,5);
lh_net = [];
       
   for nnet = 1:numel(lh_table)
       lh_net(nnet,:) = mean(lh_surf(find(lh_label == lh_table(nnet)),:));
   end
  

rh_table = rh_ct.table(:,5);
rh_net = [];
       
   for nnet = 1:numel(rh_table)
       rh_net(nnet,:) = mean(rh_surf(find(rh_label == rh_table(nnet)),:));
   end
%%
lh_surf_7net = [];
for nv = 1:40962
   if lh_label(nv) == 0 
        lh_surf_7net(nv,1) = 0;
   else
        lh_surf_7net(nv,1) = lh_net(find(lh_table == lh_label(nv)),1);
   end
end

rh_surf_7net = [];
for nv = 1:40962
   if rh_label(nv) == 0 
        rh_surf_7net(nv,1) = 0;
   else
        rh_surf_7net(nv,1) = rh_net(find(rh_table == rh_label(nv)),1);
   end
end
%% 

Netpath = '/home/eyre/project101/collosum';
Net1_lh = squeeze(load_mgh([Netpath '/lh_network_1_fs6.mgh']));
Net1_rh = squeeze(load_mgh([Netpath '/rh_network_1_fs6.mgh']));
ind_lh = find(Net1_lh==0);
ind_rh = find(Net1_rh==0);
 
lh_ind = find(Net1_lh==0);
lh = 0*Net1_lh;
lh(lh_ind) = lh_surf_7net(lh_ind);
rh_ind = find(Net1_rh==0);
rh = 0*Net1_rh;
rh(rh_ind) = rh_surf_7net(rh_ind);
%% 
opath = ['/home/eyre/project101/analysis02/seed_groupFC/net7/'] 
mkdir(opath);
save_mgh(lh,fullfile(opath,'lh_net7.mgh'),eye(4));   
save_mgh(rh,fullfile(opath,'rh_net7.mgh'),eye(4));   

  