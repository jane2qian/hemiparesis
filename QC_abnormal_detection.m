
%-------------------------- baseline abnorm ------------------------------%
clear
clc


% cerebrum V6
load('/home/ruiqing/abnormal-detection/script_v2/utility/FSNon_VolSeedVoxel_uninclude_wm_ven.mat');
% LeftSeedVox
% RightSeedVox
basepath = '/home/ruiqing/abnormal-detection/script_v2/cerebrum_uninclude_wm_ven';

len_lh = length(LeftSeedVox);
len_rh = length(RightSeedVox);

load([basepath,'/mean_vector_z.mat']);
load([basepath,'/LH_mean_corr_z.mat']);
load([basepath,'/RH_mean_corr_z.mat']);

% APP 下载文件夹一级文件路径 
lmu_path = '/home/ruiqing/QC/dset/POINT_HENAN_MDD/';
inpath = lmu_path;
rdir = dir(inpath);
dnames = {rdir.name}';
dnames(1:2) = [];

pathlist = cell(length(dnames),1);
sublist = dnames;
for i = 1:length(dnames)
    dir1 = dir([lmu_path,'/',dnames{i}]);
    snames = {dir1.name}';
    snames(1:2) = [];
    for j = 1:size(snames,1)
        path = [lmu_path,dnames{i},'/',snames{j},'/Preprocess/'];
            pathlist{i} = path;
          
    end
end

hx_patient_opath = '/home/ruiqing/QC/dset/FSNon_VolSeedVoxel_uninclude_wm_ven/';
outpath_1 = hx_patient_opath;
mkdir(outpath_1);
% template
hdr=MRIread('/home/ruiqing/abnormal-detection/script_v2/utility/Brain.1.100.1.download_subjects_sorted1002.txt.2mm.nii.gz');

for pp = 1:length(pathlist)
    
    % all brain
    inpath = pathlist{pp};
    subs = {sublist{pp}};
    
    outpath = [outpath_1,'Corr'];
    mkdir(outpath);
    
    sub = subs{1};
    pat = '\_|\s';
    DataPath=fullfile(inpath,'vol',strcat(sub,'*rest_reorient_skip_faln_mc_g1000000000_bpss_resid_FS1mm_FS2mm_sm6*.nii.gz'));
    func_data_list = dir(DataPath);
    %ruiqing
    % extract baseline mean signal 
    % compute correlation within and between left&right hemisphere of all
    % bold in abnormal subjects 
    %
    %
    if size(func_data_list,1)>0
        func_data = [];
        func_data_3d=[];
        lh_seed = [];
        rh_seed = [];
        for bld = 1:length(func_data_list)
            bld_mri = MRIread(fullfile(inpath,'vol',func_data_list(bld,1).name));
            bld_mri_vol = bld_mri.vol;

            [d1,d2,d3,t4] = size(bld_mri_vol);
            d = d1*d2*d3;
            volume = reshape(bld_mri_vol, d,t4);
            clear bld_mri_vol bld_mri
            
            lh_seed1 = zeros(len_lh,t4);
            rh_seed1 = zeros(len_rh,t4);
            
            for lh = 1:len_lh
                lh_seed1(lh,:) = mean(volume(LeftSeedVox{lh},:),1);
            end
            for rh = 1:len_rh
                rh_seed1(rh,:) = mean(volume(RightSeedVox{rh},:),1);
            end
            
            lh_seed = [lh_seed,lh_seed1];
            rh_seed = [rh_seed,rh_seed1];
            clear volume
        end
        
        %%%%%%% LH-LH
        [rho,~] = corr(lh_seed',lh_seed');
        rho(isnan(rho))=0;
        rho_z=0.5*(log2(1+rho)-log2(1-rho));
        for m=1:size(rho_z,1)
            rho_z(m,m)=1;
        end
        save([outpath '/' sub '_lh_to_lh.mat'],'rho','rho_z');
        
        %%%%%%% LH-RH
        [rho,~] = corr(lh_seed',rh_seed');
        rho(isnan(rho))=0;
        rho_z=0.5*(log2(1+rho)-log2(1-rho));
        
        save([outpath '/' sub '_lh_to_rh.mat'],'rho','rho_z');
        
        %%%%%%% RH-RH
        [rho,~] = corr(rh_seed',rh_seed');
        rho(isnan(rho))=0;
        rho_z=0.5*(log2(1+rho)-log2(1-rho));
        for m=1:size(rho_z,1)
            rho_z(m,m)=1;
        end
        save([outpath '/' sub '_rh_to_rh.mat'],'rho','rho_z');
        
        %%%%%%% RH-LH
        [rho,~] = corr(rh_seed',lh_seed');
        rho(isnan(rho))=0;
        rho_z=0.5*(log2(1+rho)-log2(1-rho));
        
        save([outpath '/' sub '_rh_to_lh.mat'],'rho','rho_z');
    else
        continue

    end
    
    % S2
    
    path =outpath;
    SIDs = subs;
    corr_LH=zeros(length(SIDs),len_lh);
    outpath2 = outpath_1;
    mkdir([outpath2,'abnor_score']);
    mkdir([outpath2,'abnor_score_final']);
    Abnormality_score=zeros(length(SIDs),len_lh);
    
    for s=1:length(SIDs)
        sid = SIDs{s};
        
        [LL, LR] = R2Z_CECE_LH(path,sid);
        L = [LL LR];
        for v =1:len_lh
            corr1=corrcoef(LH_Big_mean(v,:),L(v,:));
            corr_LH(s,v)=corr1(2,1);
            Abnormality_score(s,v)=abs(corr_LH(s,v)-LH_corr_mean(:,v))./LH_corr_std(:,v); %%% z-score
        end
    end
    save([outpath2,'abnor_score/Abnormality_Score_',sid,'_LH_z'], 'Abnormality_score');
    
    for s=1:length(SIDs)
        sid = SIDs{s};
        
        mergevol=0*hdr.vol;
        Ind_Abnormality_score=squeeze(Abnormality_score(s,:));
        
        for k=1:length(LeftSeedVox)
            mergevol(LeftSeedVox{k})=Ind_Abnormality_score(k);
        end
        
        hdr.vol=mergevol;
        MRIwrite(hdr,[outpath2 'abnor_score_final/' sid '_abnor_score_LH_z.nii.gz']);
        
    end
    
    
    %%%%%%%%
    %% RH %%
    %%%%%%%%
    
    corr_RH=zeros(length(SIDs),len_rh);
    
    Abnormality_score=zeros(length(SIDs),len_rh);
    for s=1:length(SIDs)
        sid = SIDs{s};
        
        [RR,RL] = R2Z_CECE_RH(path,sid);
        R = [RR RL];
        for v =1:len_rh
            corr1=corrcoef(RH_Big_mean(v,:),R(v,:));
            corr_RH(s,v)=corr1(2,1);
            Abnormality_score(s,v)=abs(corr_RH(s,v)-RH_corr_mean(:,v))./RH_corr_std(:,v); %%% z-score
        end
    end
    save([outpath2,'abnor_score/Abnormality_Score_',sid,'_RH_z'], 'Abnormality_score');
    
    for s=1:length(SIDs)
        sid = SIDs{s};
        
        mergevol=0*hdr.vol;
        Ind_Abnormality_score=squeeze(Abnormality_score(s,:));
        
        for k=1:length(RightSeedVox)
            mergevol(RightSeedVox{k})=Ind_Abnormality_score(k);
        end
        
        
        hdr.vol=mergevol;
        MRIwrite(hdr,[outpath2 'abnor_score_final/' sid '_abnor_score_RH_z.nii.gz']);
        
    end
    
    % S3
    
    inpath = outpath2;
    outpath3 = [inpath,'sum/'];
    mkdir(outpath3);
    outpath3 = [inpath,'sum/',sub,'/'];
    mkdir(outpath3);
    
    for i=1:length(SIDs)
        
        [t1,t2,t3,t4]=size(hdr.vol);
        sub = SIDs{i};
        a = MRIread([outpath2 'abnor_score_final/' sub '_abnor_score_LH_z.nii.gz']);
        b = MRIread([outpath2 'abnor_score_final/' sub '_abnor_score_RH_z.nii.gz']);
        
        lh = reshape(a.vol,t1*t2*t3,1);
        rh = reshape(b.vol,t1*t2*t3,1);
        result = zeros(t1*t2*t3,1);
        result = [result,lh,rh];
        
        abnor_result = zeros(t1*t2*t3,1);
        for j=1:size(result,1)
            abnor_result(j,1) = max(abs(result(j,:)));
        end
        
        hdr.vol = reshape(abnor_result,t1,t2,t3,1);
        MRIwrite(hdr,[outpath3 sub '_ALT_Abnormality_Score_all_z40_z.nii.gz']);
    end
    
    
end







