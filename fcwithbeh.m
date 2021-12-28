%------------------------ roi FC 批量计算 ---------------------------------%
clear all
clc
%% 

Netpath = '/home/eyre/project101/collosum';
Net1_lh = squeeze(load_mgh([Netpath '/lh_network_1_fs6.mgh']));
Net1_rh = squeeze(load_mgh([Netpath '/rh_network_1_fs6.mgh']));
ind_lh = find(Net1_lh==0);
ind_rh = find(Net1_rh==0);
Net1_all = [Net1_lh;Net1_rh];
ind_all = find(Net1_all==1);
Ind_all = find(Net1_all==0); % 有效范围
%% 

path = '/home/eyre/project101';
[ndata, text] = xlsread('/home/eyre/project101/seedfs_info.xlsx',2);
text(1,:) = [];
%% 

subs = text(:,1);
hemi = text(:,3);
%% 
groups = text(:,10);
itbs = subs(find(strcmp(groups,'iTBS')));
ctbs = subs(find(strcmp(groups,'cTBS')));

%% 
jiansub = subs(find(strcmp(hemi,'jian')));
illsub = subs(find(strcmp(hemi,'ill')));
%% 

lhsubs = subs(find(strcmp(hemi,'lh')));
rhsubs = subs(find(strcmp(hemi,'rh')));
%% 

lh_indice = find(strcmp(hemi,'lh'));
rh_indice = find(strcmp(hemi,'rh'));

%% 
%     cond = ['rest',num2str(i)];
    path1 = [path,'/seedfc_pre'];
    path2 = [path,'/seedfc_post'];
    % 所有被试ID
    group_fclist = [];
    group_fclist1 = [];
    ind_avg_pre = [];
    ind_avg_post =[];
    allsub_avg_pre = [];
    allsub_avg_post=[];
    lh =[];lh1=[];
    rh=[];rh1=[];allsub_delta=[];
    %% 
    
    for j = 1:length(itbs)

        sub = itbs{j};
        spp = [path1,'/',sub];
        [~,tpath]=system(['ls -f ', [spp,'/*/fcmap/lh_corr.mgh']]);
        lh_path=strsplit(tpath)';
        lh_path(end)=[];
        lh_path = lh_path{1};
        rh_path = strrep(lh_path,'lh_','rh_');

        lh_d = load_mgh(lh_path);
        rh_d = load_mgh(rh_path);
        
        if size(lh_d,1)==40962
            all_d = [lh_d;rh_d];
        else
            all_d = [lh_d,rh_d]';
        end

      
        group_fclist = [group_fclist,all_d];


        spp1 = [path2,'/',sub];
        [~,tpath1]=system(['ls -f ', [spp1,'/*/fcmap/lh_corr.mgh']]);
        lh_path1=strsplit(tpath1)';
        lh_path1(end)=[];
        lh_path1 = lh_path1{1};
        rh_path1 = strrep(lh_path1,'lh_','rh_');

        lh_d1 = load_mgh(lh_path1);
        rh_d1 = load_mgh(rh_path1);
        
        if size(lh_d1,1)==40962
            all_d1 = [lh_d1;rh_d1];
        else
            all_d1 = [lh_d1,rh_d1]';
        end
%         lh_pre = lh_d;
%         lh_pre(lh_pre<0.2)=0;
%         lh_post = lh_d1;
%          lh_post(lh_post<0.2)=0;
%         lh_delta = abs(lh_post-lh_pre);
%         lh = [lh,mean(lh_delta(lh_post~=0))];
%         rh_pre = rh_d;
%         rh_pre(rh_pre<0.2)=0;
%         rh_post = rh_d1;
%         rh_post(rh_post<0.2)=0;
%         rh_delta = abs(rh_post-rh_pre);
%         rh = [rh,mean(rh_delta(rh_post~=0))];
        group_fclist1 = [group_fclist1,all_d1];
%         ind_avg_pre = mean(all_d(all_d>0.2));
%         ind_avg_post = mean(all_d1(all_d1>0.2));
%         allsub_avg_pre = [allsub_avg_pre,ind_avg_pre];
%         allsub_avg_post=[allsub_avg_post,ind_avg_post];
% %         all_d1(all_d1<0.2)=0;all_d(all_d<0.2)=0;
%         deltafc = abs(all_d1-all_d);
%         allsub_delta = [allsub_delta,mean(deltafc(deltafc~=0))];
    
    end

%% 
   mean_group_fc = mean(group_fclist,2);

    len = length(mean_group_fc);
    outpath = [path,'/analysis02/seed_groupFC/tbs/ctbs/pre'];
    mkdir(outpath);
    save_mgh(mean_group_fc(1:len/2), [outpath '/lh.avg_corr_across_fs6.mgh'],eye(4))
    save_mgh(mean_group_fc(len/2+1:len), [outpath '/rh.avg_corr_across_fs6.mgh'],eye(4))
    %% 
    
    
[h,p,CI,STATS] = ttest(group_fclist1',group_fclist','Alpha',0.05);
tval = STATS.tstat;
group_tval = tval;
group_tval(h~=1) = 0;
% group_tval(p>0.05) = 0;
dline = find(group_tval==0 | isnan(group_tval)==1);
%group_tval= -log10(group_tval);
group_tval(dline) = 0;
len = length(group_tval);
% group_ttest_plist = p;
% group_ttest_plist(group_ttest_plist>0.05) = 0;
% dline = find(group_ttest_plist==0 | isnan(group_ttest_plist)==1);
% group_ttest_plist = -log10(group_ttest_plist);
% group_ttest_plist(dline) = 0;
% len = length(group_ttest_plist);
outpath = [path,'/analysis02/seed_groupFC/tbs/ctbs/seed_pre_post'];
 mkdir(outpath);

save_mgh(group_tval(1:len/2), [outpath '/lh.pre_post_tmap_p0.05_fs6.mgh'],eye(4))
save_mgh(group_tval(len/2+1:len), [outpath '/rh.pre_post_tmap_p0.05_fs6.mgh'],eye(4))
    
%% 
FMA = ndata(:,3);
FMA_upper = ndata(:,4);
FMA_lower = ndata(:,5);

%% 
lh_FMA = ndata(lh_indice,9);
rh_FMA = ndata(rh_indice,9);

%% 

deltaFC = (group_fclist1'-group_fclist');%./group_fclist';
pos_deltaFC = deltaFC;
pos_deltaFC(pos_deltaFC<0)=0;

%% 
z_deltaFC = deltaFC;
z_deltaFC = normalize(z_deltaFC);
%% 
FMA_tbs = FMA_upper(find(strcmp(groups,'iTBS')));
[rho,pval] = corr(pos_deltaFC,FMA_tbs);
rho(pval<0.05)=0;
[h, crit_p, adj_ci_cvrg, adj_p]= fdr_bh(pval', 0.01);
len = length(rho);
behpath = [path '/analysis02/seed_groupFC/tbs/itbs/FCwithFMA_upper/uncorp0.05/'];
mkdir(behpath);
save_mgh(rho(1:len/2), [behpath 'lh.abs_deltafc_FMA.mgh'],eye(4));
save_mgh(rho(len/2+1:len), [behpath 'rh.abs_deltafc_FMA.mgh'],eye(4));
%% 

pos_rho=rho;
pos_rho(pos_rho<0)=0;
save_mgh(pos_rho(1:len/2), [path '/FCwithBeh/lh.pos_abs_deltafc_FMA_lower_0.001.mgh'],eye(4));
save_mgh(pos_rho(len/2+1:len), [path '/FCwithBeh/rh.pos_abs_deltafc_FMA_lower_0.001.mgh'],eye(4));
%% 

neg_rho=rho;
neg_rho(neg_rho>0)=0;
save_mgh(pos_rho(1:len/2), [path '/FCwithBeh/lh.neg_abs_deltafc_FMA_lower_0.001.mgh'],eye(4));
save_mgh(pos_rho(len/2+1:len), [path '/FCwithBeh/rh.neg_abs_deltafc_FMA_lower_0.001.mgh'],eye(4));
