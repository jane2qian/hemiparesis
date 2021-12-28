import os
from os.path import abspath, dirname, join
import h5py
import nibabel as nib
import numpy as np
from nilearn import masking
import pandas as pd
import glob
from scipy.stats import ttest_rel
import time
from nilearn.connectome import ConnectivityMeasure
from nilearn.regions import RegionExtractor
import gc

root_dir = dirname(abspath("_file_"))
res_dir = join(root_dir,'analysis01')
pdir = join(root_dir,'dset/hemiparasis')

#fdir = [join(x,y,'Preprocess/surf') for x in subdir for y in glob.glob(join(x,'*rest[0:1]'))]
#print('subjects number :',len(subdir),'/n'
#     ,'functional dir number:',len(fdir))

#fsdset = [join(x,y) for x in fdir for y in glob.glob(join(x,'*fsaverage6_sm6.nii.gz'))]
subjects = os.listdir(pdir)

#def Single_Sess_Pre_Post(subjects):
for sub in subjects:
subdir = join(root_dir,'dset/hemiparasis',sub) 
print('subject path:',subdir)
ot_dir = join(res_dir,sub[0:5])
if not os.path.exists(ot_dir):
    os.mkdir(ot_dir)
#for run in np.arange(3,6):
            #matrice =[]
            fp = join(ot_dir, 'run{0:03d}corr_lr'.format(run) +'.h5')
            if not os.path.exists(fp):
                print('Create new .h5 file at', fp)
                read_type = 'w'
            else:
                print('Read existed .h5 file', fp)
                read_type = 'r+'
            f = h5py.File(fp, read_type)
            #f.close()
            for sess in np.arange(0,2):
                sessdir = glob.glob(join(subdir,'*_rest{0:01d}'.format(sess)))
                print('session dir path:',sessdir[0])
                fdir = join(subdir,sessdir[0],'Preprocess/surf') 
                #print('fdir : ',fdir)
                files = [glob.glob(join(fdir,h + 'h_allrun.nii.gz'.format(run))) for h in ['r','l'] ]
                print('file[0] :',files[0])
                imgr = nib.load(files[0][0]).get_fdata()
                imgl = nib.load(files[1][0]).get_fdata()
                print('@@@@@@@@ HERE IS OK!!!!!!')
                nvox = imgr.shape[0]*imgr.shape[1]*imgr.shape[2]
                nTRs = imgr.shape[3]
                ts_r = imgr.reshape(nvox,nTRs).astype('single')
                ts_l = imgl.reshape(nvox,nTRs).astype('single')
                #ts = np.append(ts_l,ts_r,axis=0)
                print('######### HERE IS OK!!!!!!')
                corr_lr = np.zeros((nvox,nvox))
                for iv in np.arange(nvox):
                   # print('ts_l[iv,:]: ',ts_l[iv,:].shape)
                   # stime=time.time()
                    for ivv in np.arange(nvox):
                        corr = np.corrcoef(ts_l[iv,:],ts_r[ivv,:])
                    #print('corr',corr.shape)
                        corr_lr[ivv,iv] = corr[np.tril_indices(corr.shape[0],k=-1)][0]
                        del corr
                        gc.collect()
                    print('&&&&&&&&&& HERE IS OK!!!!!!')
                    #matrice.append(corr)
                    #etime=time.time()
               # dur=etime-stime
                #print('one session processed time:',dur)
                #print('one run correlation done!')
                #print('run{} in both hemisphere correlation matrice computed!'.format(run))
                #print("Directory '%s' created" %ot_dir)
                #f = h5py.File(join(ot_dir, 'run{0:03d}corr'.format(run) +'.h5'), 'r')
                f.create_dataset('lrh_rest{}_corr'.format(run,sess), data=corr_lr)
                print('dataset created')
                del corr_lr
                gc.collect()
          #  stats,p = group_compare(matrice[0],matrice[1])
          #  return stats,p
    #f.close()