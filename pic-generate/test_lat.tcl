set colscalebarflag 0
redraw;
setfile rgb "1.tiff";
rotate_brain_x 0
puts $rgb;
redraw;
save_tiff 1.tiff;
redraw;
exit
