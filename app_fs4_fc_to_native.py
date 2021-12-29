# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import nibabel as nib
import numpy as np
import sh 
import sys
import os
import glob

def generate_vertex_map_between_native_and_fs4(subject, viz_recon_path):
    for hemi in ['lh', 'rh']:
        # construct artificial fs4 mgh file with content from 0 to 2561
        artificial_fs4 = np.arange(2562, dtype=np.float32)
        artificial_fs4_hemi = os.path.join(viz_recon_path, 'fs4_%s_artificial.mgz' % hemi)
        img = nib.MGHImage(artificial_fs4, np.eye(4))
        nib.save(img, artificial_fs4_hemi)
        
        os.system(os.path.join('cd ', viz_recon_path))
        os.system('ln -s /usr/local/freesurfer/subjects/fsaverage4 {}'.format(os.path.join(viz_recon_path,'fsaverage4')))
        os.system(os.path.join('SUBJECTS_DIR=',viz_recon_path))

        # map the artificial fs4 file to native via mri_surf2surf
        native_artificial_name = 'native_%s_artificial.mgz' % hemi
        artificial_native_hemi = os.path.join(viz_recon_path, native_artificial_name)
        args = [
            "--srcsubject", "fsaverage4",
            "--sval", artificial_fs4_hemi,
            "--trgsubject", subject,
            "--tval", artificial_native_hemi,
            "--hemi", hemi,
            "--mapmethod", "nnf"
        ]
        sh.mri_surf2surf(*args, _out=sys.stdout)

        artificial_native_img = nib.load(artificial_native_hemi)
        artificial_native = artificial_native_img.get_fdata(dtype=np.float32).flatten()

        # Write .txt files for interactive connectivity from native to fs4 and zip them
        fname = os.path.join(viz_recon_path, 'fs4_map_%s_native.txt' % hemi)
        print(fname)
        with open(fname, 'w') as f:
            f.writelines(['%s\n' % item for item in map(int, artificial_native)])
#        (fname + '.gz').remove_p()
#        sh.gzip(fname, _out=sys.stdout)

def generate_vertex_map_between_native_and_fs6(subject, viz_recon_path):
    for hemi in ['lh', 'rh']:
        # construct artificial fs6 mgh file with content from 0 to 2561
        artificial_fs6 = np.arange(40962, dtype=np.float32)
        artificial_fs6_hemi = os.path.join(viz_recon_path, 'fs6_%s_artificial.mgz' % hemi)
        img = nib.MGHImage(artificial_fs6, np.eye(4))
        nib.save(img, artificial_fs6_hemi)
        os.system(os.path.join('cd ', viz_recon_path))
        os.system('sudo ln -s /usr/local/freesurfer/subjects/fsaverage6 {}'.format(os.path.join(viz_recon_path,'fsaverage6')))
        os.system(os.path.join('SUBJECTS_DIR=',viz_recon_path))
        # map the artificial fs6 file to native via mri_surf2surf
        native_artificial_name = 'native_%s_artificial.mgz' % hemi
        artificial_native_hemi = os.path.join(viz_recon_path, native_artificial_name)
        args = [
            "--srcsubject", "fsaverage6",
            "--sval", artificial_fs6_hemi,
            "--trgsubject", subject,
            "--tval", artificial_native_hemi,
            "--hemi", hemi,
            "--mapmethod", "nnf"
        ]
        sh.mri_surf2surf(*args, _out=sys.stdout)
        artificial_native_img = nib.load(artificial_native_hemi)
        artificial_native = artificial_native_img.get_fdata(dtype=np.float32).flatten()
        # Write .txt files for interactive connectivity from native to fs4 and zip them
        fname = os.path.join(viz_recon_path, 'fs6_map_%s_native.txt' % hemi)
        print(fname)
        with open(fname, 'w') as f:
            f.writelines(['%s\n' % item for item in map(int, artificial_native)])
#        (fname + '.gz').remove_p()
#        sh.gzip(fname, _out=sys.stdout)

#generate_vertex_map_between_native_and_fs4(subject, viz_recon_path)

def generate_vertex_map_between_native_and_fs(subject, viz_recon_path):
    for hemi in ['lh', 'rh']:
        # construct artificial fs6 mgh file with content from 0 to 2561
        artificial_fs = np.arange(163842, dtype=np.float32)
        artificial_fs_hemi = os.path.join(viz_recon_path, 'fs_%s_artificial.mgz' % hemi)
        img = nib.MGHImage(artificial_fs, np.eye(4))
        nib.save(img, artificial_fs_hemi)   
        os.system(os.path.join('cd ', viz_recon_path))
        os.system('sudo ln -s /usr/local/freesurfer/subjects/fsaverage {}'.format(os.path.join(viz_recon_path,'fsaverage')))
        os.system(os.path.join('SUBJECTS_DIR=',viz_recon_path))
        # map the artificial fs file to native via mri_surf2surf
        native_artificial_name = 'native_%s_artificial.mgz' % hemi
        artificial_native_hemi = os.path.join(viz_recon_path, native_artificial_name)
        args = [
            "--srcsubject", "fsaverage",
            "--sval", artificial_fs_hemi,
            "--trgsubject", subject,
            "--tval", artificial_native_hemi,
            "--hemi", hemi,
            "--mapmethod", "nnf"
        ]
        sh.mri_surf2surf(*args, _out=sys.stdout)
        artificial_native_img = nib.load(artificial_native_hemi)
        artificial_native = artificial_native_img.get_fdata(dtype=np.float32).flatten()
        # Write .txt files for interactive connectivity from native to fs4 and zip them
        fname = os.path.join(viz_recon_path, 'fs_map_%s_native.txt' % hemi)
        print(fname)
        with open(fname, 'w') as f:
            f.writelines(['%s\n' % item for item in map(int, artificial_native)])
#        (fname + '.gz').remove_p()
#        sh.gzip(fname, _out=sys.stdout)

####################################################
path = '/home/eyre/project101/dset/hemiparasis/'
slist = os.listdir(path)


for i in range(len(slist)):
sub = slist[i]
# BM010_LiuHuaiyu
BM014_YueWenping
BM015_CaoGuiling
BM016_ZhenLi
BM018_MaJunqing
BM022_ZhangChangch
BM025_LiuXin
BM028_zhaojianyong
BM030_Wusiyu
BM036_qiyunlong


x=['BM008_LiYanbin','BM009_QiYan']
#['BM010_LiuHuaiyu','BM014_YueWenping','BM015_CaoGuiling','BM016_ZhenLi','BM018_MaJunqing',
#'BM022_ZhangChangch'
#,'BM025_LiuXin'
#,'BM028_zhaojianyong'
#,'BM030_Wusiyu'
#,'BM036_qiyunlong']
#['BM008_LiYanbin','BM009_QiYan']
#,'BM010_LiuHuaiyu',
#            'BM014_YueWenping','BM015_CaoGuiling','BM016_ZhenLi']
y=[55352,41554]#[51015,60528]#,52487,48692,50931,48854,56833,52450,55064,35597]
z=['rh','rh']#['lh','rh','lh','lh','rh','lh','lh','lh','rh','rh']

ot=[]
for s,nfs,h in zip(x,y,z):
    sub = s
    plist = glob.glob(os.path.join(path+sub+'/*0/Recon'))
    viz_recon_path = plist[0]
    os.environ['SUBJECTS_DIR']=viz_recon_path
    subject = (sub + '_reconall')
    generate_vertex_map_between_native_and_fs6(subject, viz_recon_path)
    a= nib.load(os.path.join(viz_recon_path,'native_{}_artificial.mgz'.format(h))).get_fdata()
    print(a[nfs])
    ot.append(a[nfs])


