set inpath = /home/mengyuan/ASD/data/henan/Sub012/2021-05-20-21-43-37-012/fcmap_abnorm/seed_lh383

set files = `ls $inpath/lh*.mgh`
set template = (fsaverage6)
set model = (inflated)
set fmin = 0.2
set fmax = 0.6


foreach file ($files)
echo $file
tksurfer $template lh $model -overlay $file -fminmax  $fmin $fmax -invphaseflag 0 -tcl /home/mengyuan/Snapshots/test_lat.tcl
mv 1.tiff ${file}.tiff

tksurfer $template lh $model -overlay $file -fminmax $fmin $fmax  -invphaseflag 0 -tcl /home/mengyuan/Snapshots/test_med.tcl
mv 1.tiff ${file}_med.tiff

end

set files = `ls $inpath/rh*.mgh`

foreach file ($files)
echo $file
tksurfer $template rh $model -overlay $file -fminmax $fmin $fmax -colscalebarflag 0 -invphaseflag 0 -tcl /home/mengyuan/Snapshots/test_lat.tcl
mv 1.tiff ${file}.tiff

tksurfer $template rh $model -overlay $file -fminmax $fmin $fmax -colscalebarflag 0 -invphaseflag 0 -tcl /home/mengyuan/Snapshots/test_med.tcl
mv 1.tiff ${file}_med.tiff

end

end
