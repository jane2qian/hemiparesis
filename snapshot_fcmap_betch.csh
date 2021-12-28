#!/bin/tcsh 


set inpath = $1

set files = `ls $inpath/lh*`
set template = ($2)
set model = ($3)
set fmin = $4
set fmax = $5


foreach file ($files)
echo $file
tksurfer $template lh $model -overlay $file -fminmax  $fmin $fmax -invphaseflag 0 -tcl /home/cajal/Documents/code/FC/FCmap/script/test_lat.tcl
mv 1.tiff ${file}.tiff

tksurfer $template lh $model -overlay $file -fminmax $fmin $fmax  -invphaseflag 0 -tcl /home/cajal/Documents/code/FC/FCmap/script/test_med.tcl
mv 1.tiff ${file}_med.tiff

end

set files = `ls $inpath/rh*`

foreach file ($files)
echo $file
tksurfer $template rh $model -overlay $file -fminmax $fmin $fmax -colscalebarflag 0 -invphaseflag 0 -tcl /home/cajal/Documents/code/FC/FCmap/script/test_lat.tcl
mv 1.tiff ${file}.tiff

tksurfer $template rh $model -overlay $file -fminmax $fmin $fmax -colscalebarflag 0 -invphaseflag 0 -tcl /home/cajal/Documents/code/FC/FCmap/script/test_med.tcl
mv 1.tiff ${file}_med.tiff

end

end
