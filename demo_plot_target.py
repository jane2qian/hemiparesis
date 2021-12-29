import os
import warnings
import zipfile
from pathlib import Path
import nibabel as nib
import numpy as np
import pandas as pd
import ants
import sh
import pyvista as pv


def surf2mesh(surf_file):
    coords, faces = nib.freesurfer.read_geometry(surf_file)
    pv_vertices = coords
    dims = np.full(shape=(faces.shape[0], 1), fill_value=3)
    pv_faces = np.hstack((dims, faces))
    mesh = pv.PolyData(pv_vertices, pv_faces)
    return mesh


if __name__ == '__main__':

    df = pd.read_excel('/home/eyre/project101/seedinfs_behavior.xlsx')
    target_idx_fs_list = df['fs'].to_list()
    FMA_recovery_list = df['FMA_recovery'].to_list()

    lh_pial_fs_file = '/usr/local/freesurfer/subjects/fsaverage/surf/lh.pial'
    lh_pial_fs_mesh = surf2mesh(lh_pial_fs_file)
    lh_pial_fs_file = '/usr/local/freesurfer/subjects/fsaverage/surf/lh.pial'
    lh_pial_fs_mesh = surf2mesh(lh_pial_fs_file)

    plotter = pv.Plotter()
    plotter.add_mesh(lh_pial_fs_mesh)
    for idx, target_idx_fs in enumerate(target_idx_fs_list):
        center = lh_pial_fs_mesh.points[target_idx_fs, :]
        sphere = pv.Sphere(radius=1, center=center)
        sphere.point_data['FMA_recovery'] = np.full(sphere.n_points, FMA_recovery_list[idx])
        # plotter.add_mesh(sphere, color='red')
        plotter.add_mesh(sphere)

    cpos = [(-269.68481076167365, 231.92558789366757, 122.0485812367097),
            (-33.50578773021698, -17.931861877441406, 15.487567901611328),
            (0.08385772328621494, -0.32280152727754197, 0.9427444278448635)]
    cpos = plotter.show(cpos=cpos, return_cpos=True)
    print(cpos)
