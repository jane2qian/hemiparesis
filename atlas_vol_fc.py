import os
from os.path import abspath, dirname, join
import h5py
import nibabel as nib
import numpy as np
from nilearn import masking,plotting,NiftiLabelsMasker
import pandas as pd
import glob
from scipy.stats import ttest_rel
from nilearn.connectome import ConnectivityMeasure
from nilearn.regions import RegionExtractor
from scipy.stats import ttest_rel
import seaborn as sns
import matplotlib.pyplot as plt

root_dir = dirname(abspath("_file_"))
res_dir = join(root_dir,'analysis01')
pdir = join(root_dir,'dset/hemiparasis/TBS')

#fdir = [join(x,y,'Preprocess/surf') for x in subdir for y in glob.glob(join(x,'*rest[0:1]'))]
#print('subjects number :',len(subdir),'/n'
#     ,'functional dir number:',len(fdir))

#fsdset = [join(x,y) for x in fdir for y in glob.glob(join(x,'*fsaverage6_sm6.nii.gz'))]
subjects = os.listdir(pdir)

#def Single_Sess_Pre_Post(subjects):
for sub in subjects:
    ot_dir = join(res_dir,sub[0:5])
    if not os.path.exists(ot_dir):
        os.mkdir(ot_dir)

pre_funcs = [ glob.glob(join(pdir,sub,'*_rest0/Preprocess/vol/allrun_MNI.nii.gz')) for sub in subjects ]
post_funcs = [ glob.glob(join(pdir,sub,'*_rest1/Preprocess/vol/allrun_MNI.nii.gz')) for sub in subjects ]

aal170_atlas = join(join(root_dir,'AAL3/AAL3v1.nii.gz'))
# Loading atlas data stored in 'labels'
labels = pd.read_csv(join(root_dir,'AAL3/AAL3labels.txt'),names=['regions'])

from nilearn.input_data import NiftiLabelsMasker
from nilearn.connectome import ConnectivityMeasure
masker =NiftiLabelsMasker(labels_img=aal170_atlas, standardize=True,
                         memory='nilearn_cache', verbose=5)
time_series =[]
for f in pre_funcs:
    time_series.append(masker.fit_transform(f[0]))

post_time_series = []
for ff in post_funcs:
    post_time_series.append(masker.fit_transform(ff[0]))

correlation_measure = ConnectivityMeasure(kind='correlation')
post_correlation_matrix = correlation_measure.fit_transform(post_time_series)

# ttest 
res_dir = join(root_dir,'analysis01')
subjects = os.listdir(join(res_dir,'AAL/pre'))
time_series =[]
post_time_series = []
for idx in np.arange(len(subjects)):
    predir=join(res_dir,'AAL/pre',subjects[idx])
    time_series.append(np.load(join(predir,'AAL170_timeseries.npy')))
    postdir=join(res_dir,'AAL/post',subjects[idx])
    post_time_series.append(np.load(join(postdir,'AAL170_timeseries_post.npy')))



pre_corr = []
post_corr=[]
for i in np.arange(31):
    #pre_corr.append(correlation_matrix[i][np.tril_indices(correlation_matrix[i].shape[0],-1)])
    post_corr.append(post_correlation_matrix[i][np.tril_indices(post_correlation_matrix[i].shape[0],-1)])

pre_corr = np.asarray(pre_corr)
post_corr = np.asarray(post_corr)

tval,pval=ttest_rel(correlation_matrix, post_correlation_matrix,alternative='less') 
from statsmodels.stats.multitest import fdrcorrection
np.fill_diagonal(pval,1)

lh_label = labels[0:108:2]
rh_label =labels[1:108:2]
lh_matrix=pval[0:108:2,0:108:2]
rh_matrix =pval[1:108:2,1:108:2]
heatmap = sns.heatmap(crosshemi_matrix,  annot=False,cbar=True, 
         cmap='ocean',xticklabels=True, yticklabels=True,vmax=0.01, vmin=0.001) 
#plt.setp(heatmap.get_xticklabels(), rotation=-15, ha="right",
#             rotation_mode="anchor")
heatmap.set_xticklabels(rh_label) 
heatmap.set_yticklabels(lh_label)  

plt.show()


# cluster map 
thresp = pval[:]
thresp[thresp>0.05]
sns.set_theme()
used_networks = crosshemi_labels
network_pal = sns.husl_palette(54, s=.45)
network_lut = dict(zip(map(str, used_networks), network_pal))
network_colors = pd.Series(used_networks).map(network_lut)
g = sns.clustermap(test, center=0, cmap="vlag",
                   row_colors=network_colors, col_colors=network_colors,
                   dendrogram_ratio=(.1, .2,.3,.4),
                   cbar_pos=(.02, .32, .03, .2),
                   linewidths=.75, figsize=(12, 13))

g.ax_row_dendrogram.remove()



crosshemi_matrix = pval[0:108:2,1:108:2]
crosshemi_labels = labels[0:108:2]


plotting.plot_connectome(pval, coordinates,colorbar=True, edge_threshold="85%")
plotting.show()

# left and right 
# 0:112,121:166 in pandas left:even number right:odd numer 
lh_labels = labels['regions'][0:108:2]
rh_labels =labels['regions'][1:108:2]
lh_matrix=pval[0:108:2,0:108:2]
rh_matrix =pval[1:108:2,1:108:2]
plotting.plot_matrix(rh_matrix,colorbar=True,labels=rh_labels.to_list(), cmap='ocean',
                     vmax=0.05, vmin=0.001,reorder=True)
plotting.show()



lh_coordinates = coordinates[0:108:2,:]
rh_coordinates = coordinates[1:108:2,:]
indice = np.where(rh_matrix<0.01)
plotting.plot_connectome(rh_matrix[indice][0], rh_coordinates[indice,:][0][0]#,colorbar=True
            , edge_vmin =0.95,edge_vmax = 0.99,edge_threshold='90%')
plotting.show()
plotting.plot_connectome(lh_matrix, lh_coordinates,colorbar=True, edge_vmin =0.95,edge_vmax = 0.99,edge_threshold='90%')
plotting.show()

crosshemi_matrix = pval[0:108:2,1:108:2]
crosshemi_labels = labels['regions'][0:108:2]
heatmap = sns.heatmap(crosshemi_matrix,  annot=False,cbar=True,  cmap='YlGnBu',xticklabels=True, yticklabels=True,vmax=0.05, vmin=0.001) 
heatmap.set_xticklabels(rh_labels) 
heatmap.set_yticklabels(lh_labels)  
plt.show()

plotting.plot_matrix(rh_matrix,colorbar=True,labels=rh_labels.to_list(), cmap='ocean',
                     vmax=0.05, vmin=0.001,reorder=True)
plotting.show()
######
post_group_matrix = post_correlation_matrix.mean(axis=0)
# Get number of labels that we have
# 35,36:Cingulate_Ant 81,82:Thalamus missing 
NUM_LABELS =170
cleaned_time_series = []
for i in np.arange(len(subjects)):
    signal = time_series[i]
    print(signal.shape)
    num_timepoints = signal.shape[0]
    # Create an array of zeros that has the correct size
    final_signal = np.zeros((num_timepoints, NUM_LABELS))
    # Get regions that are kept
    regions_kept = np.array(masker.labels_)
    # Fill columns matching labels with signal values
    cleaned_time_series.append(final_signal[:, regions_kept])

from nilearn.connectome import ConnectivityMeasure
correlation_measure = ConnectivityMeasure(kind='correlation')
correlation_matrix = correlation_measure.fit_transform(time_series)
group_correlation_matrix = correlation_matrix.mean(axis=0)
# Display the correlation matrix
# Mask out the major diagonal
np.fill_diagonal(group_correlation_matrix, 0)
labels.drop(labels.iloc[[34,35,80,81]].index,inplace=True)
plotting.plot_matrix(group_correlation_matrix, labels=labels['regions'].to_list(), colorbar=True,
                     vmax=0.8, vmin=-0.8)
# plot brain part regions


fronto_matrix = group_correlation_matrix[0:38,0:38]
np.fill_diagonal(fronto_matrix,0)
fronto_labels = labels['regions'][0:38]

##### 
fronto_matrix = post_group_matrix[0:38,0:38]
np.fill_diagonal(fronto_matrix,0)
fronto_labels = labels['regions'][0:38]

plotting.plot_matrix(fronto_matrix, labels=fronto_labels.to_list(), colorbar=True,
                     vmax=0.8, vmin=-0.8)

#####
tempo_matrix = group_correlation_matrix[38:90,38:90]
tempo_labels = labels['regions'][38:90]
np.fill_diagonal(tempo_matrix,0)
plotting.plot_matrix(tempo_matrix, labels=tempo_labels.to_list(), colorbar=True,
                     vmax=0.8, vmin=-0.8)
#####

tempo_matrix = post_group_matrix[38:90,38:90]
tempo_labels = labels['regions'][38:90]
np.fill_diagonal(tempo_matrix,0)
plotting.plot_matrix(tempo_matrix, labels=tempo_labels.to_list(), colorbar=True,
                     vmax=0.8, vmin=-0.8)


#####
subco_matrix = group_correlation_matrix[90:170,90:170]
subco_labels = labels['regions'][90:170]
np.fill_diagonal(subco_matrix,0)
plotting.plot_matrix(subco_matrix, labels=subco_labels.to_list(), colorbar=True,
                     vmax=0.8, vmin=-0.8)


###
subco_matrix = post_group_matrix[90:170,90:170]
subco_labels = labels['regions'][90:170]
np.fill_diagonal(subco_matrix,0)
plotting.plot_matrix(subco_matrix, labels=subco_labels.to_list(), colorbar=True,
                     vmax=0.8, vmin=-0.8)

## plot connectome 
coordinates = plotting.find_parcellation_cut_coords(labels_img=aal170_atlas)
fronto_coordinates = coordinates[0:38,:]
plotting.plot_connectome(fronto_matrix, fronto_coordinates,colorbar=True, edge_threshold="85%")
plotting.show()

tempo_coordinates = coordinates[38:90,:]
plotting.plot_connectome(tempo_matrix, tempo_coordinates,colorbar=True, edge_threshold="85%")
plotting.show()

subco_coordinates = coordinates[90:170,:]
plotting.plot_connectome(subco_matrix, subco_coordinates,colorbar=True, edge_threshold="85%")
plotting.show()

#####3
index = np.array([0,1,4,5,14,15,60,61,66,67,72,73])
motor_labels = labels['regions'][index]
motor_matrix = post_group_matrix[index[:, np.newaxis],index]
motor_coordinates = coordinates[index,:]
plotting.plot_connectome(motor_matrix, motor_coordinates,colorbar=True, edge_threshold="85%")
plotting.show()
motor_matrix = group_correlation_matrix[index[:, np.newaxis],index]
np.fill_diagonal(motor_matrix,0)
plotting.plot_matrix(motor_matrix, labels=motor_labels.to_list(), colorbar=True,
                     vmax=0.8, vmin=-0.8)
# precentral : 1,2 poscentral:61,62 SMA:15,16
# supermarginal:67,68,Paracentral_Lobule_L:73,74,mid frontal:5,6
# 75,76,77,78 :Caudate_L,Putamen_L;79:82 Pallidum_L,Thalamus_R
# 95:112 cerebellum 113:120 vermis 121:150：thalamus
index = np.array([74,75
        ,76,77,78,79,80,81,94,95,96,97,98,99,100,101,102,103,104,105,106
        ,107,108,109,110,111]
        ,112,113,114,115,115,116,117,118,119,120,121,122,123,124,
        125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149])

index = np.array([0,1,4,5,14,15,60,61,66,67,72,73])
motor_labels = labels['regions'][index]
motor_matrix = group_correlation_matrix[index[:, np.newaxis],index]
motor_coordinates = coordinates[index,:]
plotting.plot_connectome(motor_matrix, motor_coordinates,colorbar=True, edge_threshold="85%")
plotting.show()
# single subject plot 
_, axes = plt.subplots(1, 3, figsize=(15, 5))
for i, (matrix, ax) in enumerate(zip(correlation_matrices, axes)):
    plotting.plot_matrix(matrix, tri='lower', colorbar=False, axes=ax,
                         title='correlation, child {}'.format(i))