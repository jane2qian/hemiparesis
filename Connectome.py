from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import accuracy_score
import numpy as np
import os
from os.path import abspath, dirname, join
import nibabel as nib
from nilearn.input_data import NiftiLabelsMasker
from nilearn.connectome import ConnectivityMeasure
from nilearn.regions import RegionExtractor
import glob
import pandas as pd
# Connectomes per measure
from sklearn.covariance import LedoitWolf
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import cross_val_score,cross_val_predict
import my_estimator
from my_estimator import sklearn_regressor

### BM026 datapath:join(root_dir,'SpecialCase/withbeh-nopost')
root_dir = dirname(abspath("_file_"))
res_dir = join(root_dir,'analysis02/Conn-FMA')
pdir = join(root_dir,'dset/hemiparasis/TBS')
subjects = os.listdir(pdir)
subjects.append('BM026_ZhaoZhen')
subjects = sorted(subjects)
#pre_funcs = [ sorted(glob.glob(join(pdir,sub,'*_rest0/Preprocess/vol/allrun_MNI.nii.gz'))) for sub in subjects ]
#pre_funcs.append(glob.glob(join(root_dir,'SpecialCase/withbeh-nopost/BM026_ZhaoZhen','*_rest0/Preprocess/vol/allrun_MNI.nii.gz')))
#pre_funcs = sorted(pre_funcs)
post_funcs = [ sorted(glob.glob(join(pdir,sub,'*_rest1/Preprocess/vol/allrun_MNI.nii.gz'))) for sub in subjects ]
post_funcs.append(glob.glob(join(root_dir,'SpecialCase/withbeh-nopost/BM026_ZhaoZhen','*_rest1/Preprocess/vol/allrun_MNI.nii.gz')))
post_funcs = sorted(post_funcs)
kinds = ['correlation', 'partial correlation', 'tangent']
beh_performance = pd.read_csv(join(root_dir,'beh-26sub/beh.csv'),header=0)
fma = beh_performance['FMA_recovery'].to_list()
fma_upper = beh_performance['upper_recovery'].to_list()
fma_lower = beh_performance['lower_recovery'].to_list()
#cv = LeaveOneOut()
from sklearn.model_selection import ShuffleSplit
cv = ShuffleSplit(n_splits=10, test_size=0.10,
                            random_state=0)


aal170_atlas = join(join(root_dir,'AAL3/AAL3v1.nii.gz'))
# Loading atlas data stored in 'labels'
labels = pd.read_csv(join(root_dir,'AAL3/AAL3labels.txt'),names=['regions'])

masker =NiftiLabelsMasker(labels_img=aal170_atlas, standardize=True,
                         memory='nilearn_cache', verbose=5)
timeseries =[]
for idx in np.arange(len(post_funcs)):
    masked_funcs = masker.fit_transform(post_funcs[idx][0])
    subdir=join(root_dir,'analysis01',subjects[idx][0:5])
    if not os.path.exists(subdir):
        os.makedirs(subdir)
    np.save(join(subdir,'AAL170_timeseries_post'),masked_funcs)
    timeseries.append(masked_funcs)


for idx in subjects:
    subdir=join(root_dir,'analysis01/AAL/pre',idx[0:5])
    masked_funcs = np.load(join(subdir,'AAL170_timeseries.npy'))
    timeseries.append(masked_funcs)

iter_for_prediction = cv.split(timeseries, fma)
columns = ['measure', 'classifier', 'scores', 'iter_shuffle_split',
           'covariance_estimator']
results = dict()
for column_name in columns:
    results.setdefault(column_name, [])

## motor 
index = np.array([0,1,4,5,14,15,60,61,66,67,72,73])
motor_timeseries = []
for i in range(len(timeseries)):    
    motor_timeseries.append(timeseries[i-1][:,index])

connections = ConnectivityMeasure(
  #  cov_estimator=LedoitWolf(assume_centered=True),
    kind=kinds[0])#,vectorize=True)
conn_coefs = connections.fit_transform(timeseries)
# save conn_coef
for i in np.arange(len(subjects)):
    subdir=join(root_dir,'analysis01/AAL/pre',subjects[i][0:5])
    np.savetxt(join(subdir,'AAL170_timeseries.txt'),conn_coefs[i])


for index, (train_index, test_index) in enumerate(iter_for_prediction):
       # for est_key in sklearn_regressor.keys():
       #         print('Supervised learning: classification {0}'.format(est_key))
        #        estimator = sklearn_regressor[est_key]
                score = cross_val_score(estimator, conn_coefs,
                                        fma, scoring='r2',
                                        cv=[(train_index, test_index)])
                print('@@@@@ score = ',score)
                results['iter_shuffle_split'].append(index)
                results['measure'].append(kind)
                results['classifier'].append(est_key)
                results['scores'].append(score)
                results['covariance_estimator'].append('LedoitWolf')
    res = pd.DataFrame(results)
    res.to_csv(join(res_dir, 'scores.csv'))
all_results = pd.DataFrame(results)
all_results.to_csv(join(res_dir,'predictions.csv'))

score = cross_val_score(estimator, conn_coefs,
                                        fma, scoring='r2',
                                        cv=cv)        


