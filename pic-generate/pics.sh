homepath=/home/ruiqing/language-task
#ls ${homepath}/semantic > file1.txt
input=${homepath}/file1.txt
while IFS= read -r line
do
for sess in reading1 reading2 ;do
        
          tksurfer fsaverage6 lh inflated -fminmax 1.96 3 -colscalebarflag 0 -invphaseflag 0 \
        -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_lat.tcl \
         mv 1.tiff rh_inflated.tiff

        tksurfer fsaverage6 lh inflated -fminmax 1.96 3 -colscalebarflag 0 -invphaseflag 0 \
        -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_med.tcl \
         mv 1.tiff rh_inflated_med.tiff
        
           tksurfer fsaverage6 rh inflated -fminmax 1.96 3 -colscalebarflag 0 -invphaseflag 0 \
        -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_lat.tcl \
         mv 1.tiff rh_inflated.tiff

         tksurfer fsaverage6 rh inflated -fminmax 1.96 3 -colscalebarflag 0 -invphaseflag 0 \
        -tcl /home/ruiqing/abnormal-detection/script_v2/utility/test_med.tcl \
         mv 1.tiff rh_inflated_med.tiff 
       
    
done
done