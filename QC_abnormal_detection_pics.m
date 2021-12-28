% abnorm vol2surf and snap batch

%% indi
clear
clc

% input dir
path = '/home/ruiqing/QC/sum';
[~,tpath]=system(['ls -d ', [path,'/*']]);
tsf=strsplit(tpath)';
tsf(end)=[];


% recon dir
path1 = '/home/ruiqing/QC/dset/POINT_HENAN_MDD/';
[~,tpath1]=system(['ls -d ', [path1,'/*/*/Recon/*_reconall']]);
rsf=strsplit(tpath1)';
rsf(end)=[];
%% 
% 
% 
% for i = 1:length(tsf)
%     
% %     % recon path
% %     spath = rsf{i};
% %     recon_id = 'Recon';
% %     recon_path = strrep(spath,['/',recon_id],'');
% 
%     % recon path
%     spath = rsf{i};
%     [~,recon_id] = system(['echo $(basename ',spath,')']);
%     recon_id = strtrim(recon_id);
%     
%     recon_path = strrep(spath,['/',recon_id],'');
%     
%     spath = tsf{i};
% 
%     % vol2surf
%     % input file
%     [~,path_a]=system(['ls -f ', [spath,'/*.nii.gz']]);
%     afile=strsplit(path_a)';
%     afile(end)=[];
%     afile = afile{1};
%     
%     [~,fn] = system(['echo $(basename ',afile,')']);
%     fn = strtrim(fn);
%     
%     %opath
%     indi_path = [spath,'/indi_surf'];
%     mkdir(indi_path);
%     
%     cmd = ['SUBJECTS_DIR=',recon_path,';cd $SUBJECTS_DIR;ln -s /home/ruiqing/abnormal-detection/dset/FS2mm FS2mm;',...
%         'mri_vol2surf --mov ',afile,' --trgsubject ',recon_id,' --hemi lh --regheader FS2mm',...
%         ' --o ',indi_path,'/lh.',fn,'_indi.mgh --projfrac 0.5 --reshape --interp trilinear;',...
%         'mri_vol2surf --mov ',afile,' --trgsubject ',recon_id,' --hemi rh --regheader FS2mm',...
%         ' --o ',indi_path,'/rh.',fn,'_indi.mgh --projfrac 0.5 --reshape --interp trilinear'];
% 
%     system(cmd);
%     
%     % inflated
%     % lh
%     
%     [~,fpath]=system(['ls -f ', [spath,'/indi_surf/lh*.mgh']]);
%     lh_file = strsplit(fpath)';
%     lh_file(end)=[];
%     lh_path = lh_file{1};
%     rh_path = strrep(lh_path,'lh','rh');
%     
%     cmd = ['SUBJECTS_DIR=',recon_path,';tksurfer ',recon_id,' lh inflated -overlay ',lh_path,' -fminmax  2 4',...
%         '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/dset/test_lat.tcl;',...
%         'mv 1.tiff ',lh_path,'_indi_inflated_2_4.tiff'];
%     system(cmd);
%     
%     cmd = ['SUBJECTS_DIR=',recon_path,';tksurfer ',recon_id,' lh inflated -overlay ',lh_path,' -fminmax  2 4',...
%         '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/dset/test_med.tcl;',...
%         'mv 1.tiff ',lh_path,'_indi_inflated_med_2_4.tiff'];
%     system(cmd);
%     
%     % rh
%     cmd = ['SUBJECTS_DIR=',recon_path,';tksurfer ',recon_id,' rh inflated -overlay ',rh_path,' -fminmax  2 4',...
%         '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/dset/test_lat.tcl;',...
%         'mv 1.tiff ',rh_path,'_indi_inflated_2_4.tiff'];
%     system(cmd);
%     
%     cmd = ['SUBJECTS_DIR=',recon_path,';tksurfer ',recon_id,' rh inflated -overlay ',rh_path,' -fminmax  2 4',...
%         '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/dset/test_med.tcl;',...
%         'mv 1.tiff ',rh_path,'_indi_inflated_med_2_4.tiff'];
%     system(cmd);
%     
%     % pial
%     % lh
%     cmd = ['SUBJECTS_DIR=',recon_path,';tksurfer ',recon_id,' lh pial -overlay ',lh_path,' -fminmax  2 4',...
%         '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/dset/test_lat.tcl;',...
%         'mv 1.tiff ',lh_path,'_indi_pial_2_4.tiff'];
%     system(cmd);
%     
%     cmd = ['SUBJECTS_DIR=',recon_path,';tksurfer ',recon_id,' lh pial -overlay ',lh_path,' -fminmax  2 4',...
%         '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/dset/test_med.tcl;',...
%         'mv 1.tiff ',lh_path,'_indi_pial_med_2_4.tiff'];
%     system(cmd);
%     
%     % rh
%     cmd = ['SUBJECTS_DIR=',recon_path,';tksurfer ',recon_id,' rh pial -overlay ',rh_path,' -fminmax  2 4',...
%         '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/dset/test_lat.tcl;',...
%         'mv 1.tiff ',rh_path,'_indi_pial_2_4.tiff'];
%     system(cmd);
%     
%     cmd = ['SUBJECTS_DIR=',recon_path,';tksurfer ',recon_id,' rh pial -overlay ',rh_path,' -fminmax  2 4',...
%         '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/dset/test_med.tcl;',...
%         'mv 1.tiff ',rh_path,'_indi_pial_med_2_4.tiff'];
%     system(cmd);
%     
%     % 4-way pic
%     [~,lh_name] = system(['echo $(basename ',lh_path,')']);
%     lh_name = strtrim(lh_name);
%     sname1 = strrep(lh_name,'lh.','');
%     sname = strrep(sname1,'.mgh','');
%     
%     OutPath = [spath,'/indi_surf']; %
%     InPath = OutPath;
%     
%     % V2
%     %     range1 = 100:500;
%     %     range2 = 60:505;
%     %     inflated
%     range1 = 40:500;
%     range2 = 20:600;
%     
%     Model = cell(1,1);
%     Model{1,1} = sname;
%     
%     for kk = 1:size(Model,1)
%         model = Model{kk,1};
%         if exist([OutPath '/' model '_inflated_M22.tiff'],'file')
%             continue
%         end
%         files = dir([InPath  '/*' model '*inflated*.tiff']);
%         CombinedMap = cell(1,4);
%         
%         if length(files)==0
%             continue
%         end
%         for f = 1:length(files)
%             
%             fileName = files(f).name;
%             tmp =imread([InPath '/' fileName]);
%             tmp = tmp(range1,range2,:);
%             tmp(tmp==0) = 255;
%             CombinedMap{f} = tmp;
%         end
%         
%         Combined_Lat = cat(2,CombinedMap{1},CombinedMap{3});
%         Combined_Med = cat(2,CombinedMap{2},CombinedMap{4});
%         
%         Combined = cat(1,Combined_Lat,Combined_Med);
%         
%         imwrite(Combined,[OutPath '/' model '_inflated_M22.tiff']);
%     end
%     
%     %V1
%     range1 = 150:450;
%     range2 = 105:505;
%     
%     for kk = 1:size(Model,1)
%         model = Model{kk,1}
%         if exist([OutPath '/' model '_pial_M22.tiff'],'file')
%             continue
%         end
%         files = dir([InPath  '/*' model '*pial*.tiff']);
%         CombinedMap = cell(1,4);
%         
%         if length(files)==0
%             continue
%         end
%         for f = 1:length(files)
%             
%             fileName = files(f).name;
%             tmp =imread([InPath '/' fileName]);
%             tmp = tmp(range1,range2,:);
%             tmp(tmp==0) = 255;
%             CombinedMap{f} = tmp;
%         end
%         
%         Combined_Lat = cat(2,CombinedMap{1},CombinedMap{3});
%         Combined_Med = cat(2,CombinedMap{2},CombinedMap{4});
%         
%         Combined = cat(1,Combined_Lat,Combined_Med);
%         
%         imwrite(Combined,[OutPath '/' model '_pial_M22.tiff']);
%     end
%     
%     
% end


%% fs6

for i = 1:length(tsf)
    
    spath = tsf{i};
    
    [~,subname] = system(['echo $(basename ',spath,')']);
    subname = strtrim(subname);
    
    % surf2vol
    % input file
    [~,path_a]=system(['ls -f ', [spath,'/*.nii.gz']]);
    afile=strsplit(path_a)';
    afile(end)=[];
    afile = afile{1};
    
    [~,fn] = system(['echo $(basename ',afile,')']);
    fn = strtrim(fn);
    
    %opath
    indi_path = [spath,'/fs6_surf'];
    mkdir(indi_path);
%     
%     cmd = ['SUBJECTS_DIR=/home/ruiqing/abnormal-detection/dset/;',...
%         'mri_vol2surf --mov ',afile,' --trgsubject fsaverage6 --hemi lh --regheader FS2mm',...
%         ' --o ',indi_path,'/lh.',fn,'_fs6.mgh --projfrac 0.5 --reshape --interp trilinear;',...
%         'mri_vol2surf --mov ',afile,' --trgsubject fsaverage6 --hemi rh --regheader FS2mm',...
%         ' --o ',indi_path,'/rh.',fn,'_fs6.mgh --projfrac 0.5 --reshape --interp trilinear'];
%       
    cmd = ['SUBJECTS_DIR=/home/ruiqing/abnormal-detection/dset;',...
        'mri_vol2surf --mov ',afile,' --trgsubject fsaverage6 --hemi lh --regheader FS2mm',...
        ' --o ',indi_path,'/lh.',fn,'_fs6.mgh --projfrac 0.5 --reshape --interp trilinear;',...
        'mri_vol2surf --mov ',afile,' --trgsubject fsaverage6 --hemi rh --regheader FS2mm',...
        ' --o ',indi_path,'/rh.',fn,'_fs6.mgh --projfrac 0.5 --reshape --interp trilinear'];
    system(cmd);
    
    [~,fpath]=system(['ls -f ', [spath,'/fs6_surf/lh*.mgh']]);
    lh_file = strsplit(fpath)';
    lh_file(end)=[];
    lh_path = lh_file{1};
    rh_path = strrep(lh_path,'lh','rh');
    
    % inflated
    % lh
    cmd = ['tksurfer fsaverage6 lh inflated -overlay ',lh_path,' -fminmax  2 4',...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/dset/test_lat.tcl;',...
        'mv 1.tiff ',lh_path,'_fs6_inflated_2_4.tiff'];
    system(cmd);
    
    cmd = ['tksurfer fsaverage6 lh inflated -overlay ',lh_path,' -fminmax  2 4',...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/dset/test_med.tcl;',...
        'mv 1.tiff ',lh_path,'_fs6_inflated_med_2_4.tiff'];
    system(cmd);
    
    % rh
    cmd = ['tksurfer fsaverage6 rh inflated -overlay ',rh_path,' -fminmax  2 4',...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/dset/test_lat.tcl;',...
        'mv 1.tiff ',rh_path,'_fs6_inflated_2_4.tiff'];
    system(cmd);
    
    cmd = ['tksurfer fsaverage6 rh inflated -overlay ',rh_path,' -fminmax  2 4',...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/dset/test_med.tcl;',...
        'mv 1.tiff ',rh_path,'_fs6_inflated_med_2_4.tiff'];
    system(cmd);
    
    % pial
    % lh
    cmd = ['tksurfer fsaverage6 lh pial -overlay ',lh_path,' -fminmax  2 4',...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/dset/test_lat.tcl;',...
        'mv 1.tiff ',lh_path,'_fs6_pial_2_4.tiff'];
    system(cmd);
    
    cmd = ['tksurfer fsaverage6 lh pial -overlay ',lh_path,' -fminmax  2 4',...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/dset/test_med.tcl;',...
        'mv 1.tiff ',lh_path,'_fs6_pial_med_2_4.tiff'];
    system(cmd);
    
    % rh
    cmd = ['tksurfer fsaverage6 rh pial -overlay ',rh_path,' -fminmax  2 4',...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/dset/test_lat.tcl;',...
        'mv 1.tiff ',rh_path,'_fs6_pial_2_4.tiff'];
    system(cmd);
    
    cmd = ['tksurfer fsaverage6 rh pial -overlay ',rh_path,' -fminmax  2 4',...
        '  -colscalebarflag 0 -invphaseflag 0 -tcl /home/ruiqing/abnormal-detection/dset/test_med.tcl;',...
        'mv 1.tiff ',rh_path,'_fs6_pial_med_2_4.tiff'];
    system(cmd);
    
    % 4-way pic
    [~,lh_name] = system(['echo $(basename ',lh_path,')']);
    lh_name = strtrim(lh_name);
    sname1 = strrep(lh_name,'lh.','');
    sname = strrep(sname1,'.mgh','');
    
    OutPath = [spath,'/fs6_surf']; %
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




