
import os
import sys
import numpy as np
import pandas as pd
import pylab as pl
from sklearn.cluster import k_means
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from scipy.spatial.distance import cdist
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from pdblib.num import *


df = pd.read_csv("DNA_step_library_frame_refine.csv")
steps = df.step_name.unique()

library_all = Pdb('DNA_step_library.pdb')
library_clust = Pdb()

colors = ['red','blue','yellow','green','orange','cyan','skyblue','lime','gold','gray','violet','firebrick','purple','teal']

nclust_dic = {'CG/CG':5,'GC/GC':6,'CG/AT':5,'TA/GC':5,
              'AT/AT':3,'TA/TA':3,'AT/GC':6,'CG/TA':4,
              'GC/CG':4,'TA/CG':4,'GC/AT':5,'CG/GC':4,
              'AT/TA':4,'TA/AT':3,'GC/TA':4,'AT/CG':4,
             }

fig = pl.figure(figsize=(16,16))

subdf_list = []

for idx, step in enumerate(steps):

    print("--- Working on %s [%d/%d] ---"%(step, idx, len(steps)))
    df_step = df.loc[df['step_name'] == step]
    step_num = len(df_step)
    X = np.array([map(float,a.split(',')) for a in df_step.origin+','+df_step.M_axis])
    X2 = np.hstack((X[0:,0:3],X[0:,3:]))
    pca = PCA(n_components=2)
    pca.fit(X2)
    y = pca.fit_transform(X2)
   
    nclust = nclust_dic[step]

    model = KMeans(n_clusters=nclust)
    clusassign = model.fit_predict(X2)
    
    min_dist = np.min(cdist(X2,model.cluster_centers_,'euclidean'),axis=1)
    X2_df = pd.DataFrame(X2)
    Y = pd.DataFrame(min_dist,index=X2_df.index,columns=['Center_euclidean'])
    Z = pd.DataFrame(clusassign,index=X2_df.index,columns=['cluster_ID'])
    PAP = pd.concat([Y,Z], axis=1)

    grouped = PAP.groupby(['cluster_ID'])['Center_euclidean'].idxmin()
    count = PAP.groupby(['cluster_ID']).count()
    energy = count/len(PAP)
    clust = pd.concat([grouped,count,energy],axis=1)
    clust.columns = ['argmin','count','pop']
    clust['old_id'] = df_step.index.values[grouped.values]
    clust['step_idx'] = df_step.step_idx.values[grouped.values]
    clust['origin'] = df_step.origin.values[grouped.values]
    clust['M_axis'] = df_step.M_axis.values[grouped.values]
    clust['step_name'] = step
    clust = clust.reset_index()
    clust = clust[['step_name','cluster_ID','old_id','step_idx','argmin','count','pop','origin','M_axis']]
    subdf_list.append(clust)

    for i in range(nclust):
        c = clust['old_id'].ix[i]
        md = library_all.mds[c]
        library_clust.mds.append(md)

    for i in range(nclust):
        c = np.where(model.labels_ == i)[0]
        ax = fig.add_subplot(4,4,idx+1)
        ax.plot(y[c,0],y[c,1],'o',color=colors[i],markersize=16)
        ax.set_title(step)
        ax.set_ylim((-3,3))
        ax.set_xlim((-3,3))

    for i in range(nclust):
        c = clust['argmin'].ix[i]
        ax = fig.add_subplot(4,4,idx+1)
        ax.plot(y[c,0],y[c,1],'*',color='black',markersize=24)
        ax.set_title(step+" "+"("+str(step_num)+")")
        ax.set_ylim((-3,3))
        ax.set_xlim((-3,3))

fig.tight_layout()

pl.savefig("all_step_nclust.pdf")
pl.show()

df_final = pd.concat(subdf_list,axis=0)
df_final.reset_index(drop=True)
df_final.to_csv("DNA_step_cluster_frame.csv",index=False)

library_clust.write('DNA_step_cluster.pdb')

