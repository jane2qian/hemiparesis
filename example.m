% Those functions are from Freesurfer
% nifti, mgh, annot

%--------------1. Volumetric---------------%
hdr = MRIread('/home/zoey/ngwork/BoAi/2020-08-26-22-54-13-wangyang0824/Recon/wangyang_reconall/surf/lh.area.nii.gz');
vol = reshape(hdr.vol, [hdr.nvoxels, hdr.nframes]); % reshape to 2D
hdr.vol = reshape(vol, [hdr.volsize, hdr.nframes]); % reshape back to 4D
MRIwrite(hdr, ['output_example_preop.nii.gz']);

%--------------2. mgh file ---------------%
metric = load_mgh('lh.example.mgh');
save_mgh(metric, 'lh.output.mgh', eye(4))

%--------------3. annot file -------------%
[vertex, label, ct] = read_annotation('lh_parc_result.annot');
write_annotation('lh_parc_result_out.annot', vertex, label, ct)


[curv, fnum] = read_curv('/media/zoey/ElementsSE_4/download/ABIDE_TR30/A00032677/2020-06-14-11-29-12-A00032677/Recon/A00032677_reconall/surf/lh.volume');
[vertex_coords, faces] = read_surf('/media/zoey/ElementsSE_4/download/ABIDE_TR30/A00032677/2020-06-14-11-29-12-A00032677/Recon/A00032677_reconall/surf/lh.pial');
[vertex, label, ct] = read_annotation('/media/zoey/ElementsSE_4/download/ABIDE_TR30/A00032677/2020-06-14-11-29-12-A00032677/Parcellation-18/lh_parc18_native_surf.annot');


[vertex_coords, faces] = read_surf('/usr/local/freesurfer/subjects/fsaverage/surf/lh.orig');

