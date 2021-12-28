function snapshot(inpath,outname)

%% tmp
% range1 = 95:505;
% range2 = 10:600;
%         
% inflated_fs6 
range1 = 135:470; 
range2 = 50:540;
% 
% % inflated _indi
% range1 = 115:485; 
% range2 = 45:565;
% 
% % pial_fs6
% range1 = 143:457;
% range2 = 90:520; 

% % pial_indi
% range1 = 160:435;
% range2 = 110:505; 

%Mode1 = ['C' num2str(i)];
files = dir(fullfile(inpath,'*h*tiff'));
CombinedMap = cell(1,4);

bar_tmp=imread([inpath '/' files(1).name]);
bar=bar_tmp(421:600,521:600,:);
bar1=bar;
bar(bar1==0)=255;
bar(bar1==255)=1;

for f = 1:length(files)
    fileName = files(f).name;
    tmp = imread([inpath '/' fileName ]);
    tmp = tmp(range1,range2,:); 
    tmp(tmp==0) = 255;
    CombinedMap{f} = tmp;
    system(['rm ',inpath,'/',fileName]);
end

Combined_Lat = cat(2,CombinedMap{1},CombinedMap{3});
Combined_Med = cat(2,CombinedMap{2},CombinedMap{4});
Combined = cat(1,Combined_Lat,Combined_Med);

sizeC=size(Combined);
bar_posx=(sizeC(1)-180)/2+6;
bar_posy=(sizeC(2)-80)/2-9;
Combined_bar=Combined;
Combined_bar(bar_posx:bar_posx+179,bar_posy:bar_posy+79,:)=bar;
if strcmp(outname(1:2),'fc')
    imwrite(Combined_bar,[inpath '/' outname '.tiff']);
else
    imwrite(Combined,[inpath '/' outname '.tiff']);
end
% imwrite(Combined,[inpath '/' outname '.tiff']);

end
