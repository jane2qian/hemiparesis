%% parcellation QC
%  necessary files: 
%  lh_parc18_native_surf.annot,rh_parc18_native_surf.annot
%%
clear
clc

% input dir
path = '/home/ruiqing/QC/dset/POINT_HENAN_MDD';
[~,tpath]=system(['ls -d ', [path,'/*/*']]);
tsf=strsplit(tpath)';
tsf(end)=[];

% recon dir
[~,tpath1]=system(['ls -d ', [path,'/*/*/Recon/*reconall']]);
rsf=strsplit(tpath1)';
rsf(end)=[];

% % recon dir
% [~,tpath1]=system(['ls -d ', [path,'/*/*/Recon']]);
% rsf=strsplit(tpath1)';
% rsf(end)=[];
%% 18 parcellation projection from native space to fs6 space
% % parc18 native--> fs6
% for i = 1:length(tsf)
%     
%     % parc18
%     pp2 = [tsf{i},'/Parcellation-18'];
% 
%     % recon path
%     spath = rsf{i};
%     [~,recon_id] = system(['echo $(basename ',spath,')']);
%     recon_id = strtrim(recon_id);
%     
% %     % recon path
% %     spath = rsf{i};
% %     recon_id = 'Recon';
% 
%     rpath = strrep(spath,['/',recon_id],'');
%     
%     % parc18 surf2surf
%     cmd = ['cd ',rpath,';ln -s /usr/local/freesurfer/subjects/fsaverage6 fsaverage6'];
%     system(cmd);
% 
%     cmd = ['SUBJECTS_DIR=',rpath,';mri_surf2surf --srcsubject ',recon_id,32,'--sval-annot ',pp2,'/lh_parc18_native_surf.annot',...
%         32,'--trgsubject fsaverage6 --tval ',pp2,'/lh_parc18_fs6_surf.annot --hemi lh;',...
%         'mri_surf2surf --srcsubject ',recon_id,32,'--sval-annot ',pp2,'/rh_parc18_native_surf.annot',...
%         32,'--trgsubject fsaverage6 --tval ',pp2,'/rh_parc18_fs6_surf.annot --hemi rh'];
%     system(cmd);
% 
% end
%% 
% 
% % snap and 4way (indi)
for i = 1:length(tsf)
    
    % recon path
    spath = rsf{i};
    [~,recon_id] = system(['echo $(basename ',spath,')']);
    recon_id = strtrim(recon_id);
    
    recon_path = strrep(spath,['/',recon_id],'');

%     %recon path
       pp1 = [tsf{i},'/Parcellation-18'];
        lh_path = [pp1,'/lh_parc18_native_surf.annot'];

%      pp2 = [tsf{i},'/Parcellation-213'];
%    lh_path = [pp2,'/lh_parc213_native_surf.annot'];
    rh_path = strrep(lh_path,'lh','rh');
    
    
    
    % inflated
    % lh
    cmd = ['SUBJECTS_DIR=',recon_path,';tksurfer ',recon_id,' lh inflated -annot ',lh_path,...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_lat.tcl;',...
        'mv 1.tiff ',lh_path,'_indi_inflated.tiff'];
    system(cmd);
    
    cmd = ['SUBJECTS_DIR=',recon_path,';tksurfer ',recon_id,' lh inflated -annot ',lh_path,...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_med.tcl;',...
        'mv 1.tiff ',lh_path,'_indi_inflated_med.tiff'];
    system(cmd);
    
    % rh
    cmd = ['SUBJECTS_DIR=',recon_path,';tksurfer ',recon_id,' rh inflated -annot ',rh_path,...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_lat.tcl;',...
        'mv 1.tiff ',rh_path,'_indi_inflated.tiff'];
    system(cmd);
    
    cmd = ['SUBJECTS_DIR=',recon_path,';tksurfer ',recon_id,' rh inflated -annot ',rh_path,...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_med.tcl;',...
        'mv 1.tiff ',rh_path,'_indi_inflated_med.tiff'];
    system(cmd);

    
    %% pial
    % lh
    cmd = ['SUBJECTS_DIR=',recon_path,';tksurfer ',recon_id,' lh pial -annot ',lh_path,...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_lat.tcl;',...
        'mv 1.tiff ',lh_path,'_indi_pial.tiff'];
    system(cmd);
    
    cmd = ['SUBJECTS_DIR=',recon_path,';tksurfer ',recon_id,' lh pial -annot ',lh_path,...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_med.tcl;',...
        'mv 1.tiff ',lh_path,'_indi_pial_med.tiff'];
    system(cmd);
    
    % rh
    cmd = ['SUBJECTS_DIR=',recon_path,';tksurfer ',recon_id,' rh pial -annot ',rh_path,...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_lat.tcl;',...
        'mv 1.tiff ',rh_path,'_indi_pial.tiff'];
    system(cmd);
    
    cmd = ['SUBJECTS_DIR=',recon_path,';tksurfer ',recon_id,' rh pial -annot ',rh_path,...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_med.tcl;',...
        'mv 1.tiff ',rh_path,'_indi_pial_med.tiff'];
    system(cmd);
    
    % 4-way pic
    [~,lh_name] = system(['echo $(basename ',lh_path,')']);
    lh_name = strtrim(lh_name);
    sname1 = strrep(lh_name,'lh_','');
    sname = strrep(sname1,'.annot','');
    
    OutPath = pp1; %
    InPath = OutPath;
    
    % V2
    %     range1 = 100:500;
    %     range2 = 60:505;
    %     inflated
    range1 = 40:500;
    range2 = 20:600;

    Model = cell(1,1);
    Model{1,1} = sname;
    
    for kk = 1:size(Model,1)
        model = Model{kk,1}
        if exist([OutPath '/' model '_inflated_M22.tiff'],'file')
            continue
        end
        files = dir([InPath  '/*' model '*inflated*.tiff']);
        CombinedMap = cell(1,4);
        
        if length(files)==0
            continue
        end
        for f = 1:length(files)
            
            fileName = files(f).name;
            tmp =imread([InPath '/' fileName]);
            tmp = tmp(range1,range2,:);
            tmp(tmp==0) = 255;
            CombinedMap{f} = tmp;
        end
        
        Combined_Lat = cat(2,CombinedMap{1},CombinedMap{3});
        Combined_Med = cat(2,CombinedMap{2},CombinedMap{4});
        
        Combined = cat(1,Combined_Lat,Combined_Med);
        
        imwrite(Combined,[OutPath '/' model '_inflated_M22.tiff']);
    end
    
    %V1
    range1 = 150:450;
    range2 = 105:505;
    
    Model = cell(1,1);
    Model{1,1} = sname;
    
    for kk = 1:size(Model,1)
        model = Model{kk,1}
        if exist([OutPath '/' model '_pial_M22.tiff'],'file')
            continue
        end
        files = dir([InPath  '/*' model '*pial*.tiff']);
        CombinedMap = cell(1,4);
        
        if length(files)==0
            continue
        end
        for f = 1:length(files)
            
            fileName = files(f).name;
            tmp =imread([InPath '/' fileName]);
            tmp = tmp(range1,range2,:);
            tmp(tmp==0) = 255;
            CombinedMap{f} = tmp;
        end
        
        Combined_Lat = cat(2,CombinedMap{1},CombinedMap{3});
        Combined_Med = cat(2,CombinedMap{2},CombinedMap{4});
        
        Combined = cat(1,Combined_Lat,Combined_Med);
        
        imwrite(Combined,[OutPath '/' model '_pial_M22.tiff']);
    end
    
end


% snap and 4way (fs6)
for i = 1:length(tsf)
    
    setenv('SUBJECTS_DIR', '/usr/local/freesurfer/subjects');
     pp1 = [tsf{i},'/Parcellation-18'];
     lh_path = [pp1,  '/lh_parc18_fs6_surf.annot'];
%     
%        pp2 = [tsf{i},'/Parcellation-213'];
%    lh_path = [pp2,'/lh_parc213_fs6_surf.annot'];
%     lh_path = [pp1,'/lh_parc_result_indi_s02_fs6.annot'];
    rh_path = strrep(lh_path,'lh','rh');
    
    % inflated
    % lh
    cmd = ['tksurfer fsaverage6 lh inflated -annot ',lh_path,...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_lat.tcl;',...
        'mv 1.tiff ',lh_path,'_fs6_inflated.tiff'];
    system(cmd);
    
    cmd = ['tksurfer fsaverage6 lh inflated -annot ',lh_path,...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_med.tcl;',...
        'mv 1.tiff ',lh_path,'_fs6_inflated_med.tiff'];
    system(cmd);
    
    % rh
    cmd = ['tksurfer fsaverage6 rh inflated -annot ',rh_path,...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_lat.tcl;',...
        'mv 1.tiff ',rh_path,'_fs6_inflated.tiff'];
    system(cmd);
    
    cmd = ['tksurfer fsaverage6 rh inflated -annot ',rh_path,...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_med.tcl;',...
        'mv 1.tiff ',rh_path,'_fs6_inflated_med.tiff'];
    system(cmd);
    
%     pial
%     lh
    cmd = ['tksurfer fsaverage6 lh pial -annot ',lh_path,...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_lat.tcl;',...
        'mv 1.tiff ',lh_path,'_fs6_pial.tiff'];
    system(cmd);
    
    cmd = ['tksurfer fsaverage6 lh pial -annot ',lh_path,...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_med.tcl;',...
        'mv 1.tiff ',lh_path,'_fs6_pial_med.tiff'];
    system(cmd);
    
%     rh
    cmd = ['tksurfer fsaverage6 rh pial -annot ',rh_path,...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_lat.tcl;',...
        'mv 1.tiff ',rh_path,'_fs6_pial.tiff'];
    system(cmd);
    
    cmd = ['tksurfer fsaverage6 rh pial -annot ',rh_path,...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_med.tcl;',...
        'mv 1.tiff ',rh_path,'_fs6_pial_med.tiff'];
    system(cmd);
    
    % 4-way pic
    [~,lh_name] = system(['echo $(basename ',lh_path,')']);
    lh_name = strtrim(lh_name);
    sname1 = strrep(lh_name,'lh_','');
    sname = sname1;
%     sname = strrep(sname1,'.annot','');
    
    OutPath = pp1; %
    InPath = OutPath;
    
    % V1
    range1 = 120:470;
    range2 = 65:535;
    
    Model = cell(1,1);
    Model{1,1} = sname;
    
    for kk = 1:size(Model,1)
        model = Model{kk,1}
        if exist([OutPath '/' model '_inflated_M22.tiff'],'file')
            continue
        end
        files = dir([InPath  '/*' model '*inflated*.tiff']);
        CombinedMap = cell(1,4);
        
        if length(files)==0
            continue
        end
        for f = 1:length(files)
            
            fileName = files(f).name;
            tmp =imread([InPath '/' fileName]);
            tmp = tmp(range1,range2,:);
            tmp(tmp==0) = 255;
            CombinedMap{f} = tmp;
        end
        
        Combined_Lat = cat(2,CombinedMap{1},CombinedMap{3});
        Combined_Med = cat(2,CombinedMap{2},CombinedMap{4});
        
        Combined = cat(1,Combined_Lat,Combined_Med);
        
        imwrite(Combined,[OutPath '/' model '_inflated_M22.tiff']);
    end
    
    %V1
    range1 = 150:450;
    range2 = 105:505;
    
    for kk = 1:size(Model,1)
        model = Model{kk,1}
        if exist([OutPath '/' model '_pial_M22.tiff'],'file')
            continue
        end
        files = dir([InPath  '/*' model '*pial*.tiff']);
        CombinedMap = cell(1,4);
        
        if length(files)==0
            continue
        end
        for f = 1:length(files)
            
            fileName = files(f).name;
            tmp =imread([InPath '/' fileName]);
            tmp = tmp(range1,range2,:);
            tmp(tmp==0) = 255;
            CombinedMap{f} = tmp;
        end
        
        Combined_Lat = cat(2,CombinedMap{1},CombinedMap{3});
        Combined_Med = cat(2,CombinedMap{2},CombinedMap{4});
        
        Combined = cat(1,Combined_Lat,Combined_Med);
        
        imwrite(Combined,[OutPath '/' model '_pial_M22.tiff']);
    end
    
end








