function snapshot_4way_combine_inflated_fs6(input_path,input_dir,fmin,fmax)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
[~,lh_name] = system(['echo $(basename ',input_path,')']);
lh_name = strtrim(lh_name);
sname1 = strrep(lh_name,'lh.','');
% sname = strrep(sname1,'.mgh','');
sname = sname1;

OutPath = input_dir; %
InPath = OutPath;

% V1
range1 = 120:470;
range2 = 65:535;

Model = cell(1,1);
Model{1,1} = sname;

for kk = 1:size(Model,1)
    model = Model{kk,1}
    if exist([OutPath '/' model '_inflated_' num2str(fmin) '_' num2str(fmax) '_M22.tiff'],'file')
        continue
    end
    files = dir([InPath  '/*' model '*inflated*' num2str(fmin) '*' num2str(fmax) '*.tiff']);
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
    
    imwrite(Combined,[OutPath '/' model '_inflated_' num2str(fmin) '_' num2str(fmax) '_M22.tiff']);
end
end

