def mri_surf2surf(srcsubject, sval, trgsubject, tval, hemi):
    sh.mri_surf2surf(
        '--srcsubject', srcsubject,
        '--sval', sval,
        '--trgsubject', trgsubject,
        '--tval', tval,
        '--hemi', hemi)

pial_file = subject_dir / subj / 'surf' / 'lh.pial'
coords, faces = nib.freesurfer.read_geometry(pial_file)
target_mask_np = np.zeros(coords.shape[0], np.float32)
target_mask_np[target_idx] = 1
target_mask = ants.from_numpy(target_mask_np.reshape((-1, 1, 1)))
target_mask_file = workdir / 'target_mask.mgh'
ants.image_write(target_mask, str(target_mask_file))
target_mask_fs_file = workdir / 'target_mask_fs.mgh'
mri_surf2surf(subj, target_mask_file, 'fsaverage6', target_mask_fs_file, 'rh')
target_mask_fs = ants.image_read(str(target_mask_fs_file))
target_mask_fs_np = target_mask_fs.numpy().flatten()
target_idx_fs = np.argmax(target_mask_fs_np)
