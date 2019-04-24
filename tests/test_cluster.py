
import os
import sys
import numpy as np
import pandas as pd
import pylab as pl
from scipy.cluster.hierarchy import *
from scipy.cluster.vq import *
from scipy.spatial.distance import pdist
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans
from sklearn.cluster import k_means
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def mydist(p1,p2):
    
    dscore = np.linalg.norm(p1[0:3]-p2[0:3])
    rscore = 2*np.sum(np.abs(p1[3:]-p2[3:]))
    
    return dscore+rscore

df = pd.read_csv("DNA_step_library_frame.csv")

steps = df.step_name.unique()

colors = ['red','blue','yellow','green','orange','cyan','skyblue','lime','gold','gray','violet','firebrick','purple','teal']
nclust = 4

nclusts = range(1,11)
for j, nclust in enumerate(nclusts):
    fig = pl.figure(figsize=(10,10))
    
    for idx, step in enumerate(steps):
    
        aa = df.loc[df['step_name'] == step]
        X = np.array([map(float,a.split(',')) for a in aa.origin+','+aa.M_axis])
        X2 = np.hstack((X[0:,0:3],X[0:,3:]))
        pca = PCA(n_components=2)
        pca.fit(X2)
        y = pca.fit_transform(X2)
        
        kmeans = k_means(X2,n_clusters=nclust,precompute_distances=True)
    
        for i in range(nclust):
            c = np.where(kmeans[1]==i)[0]
            ax = fig.add_subplot(4,4,idx+1)
            ax.plot(y[c,0],y[c,1],'o',color=colors[i],markersize=12)
            ax.set_title(step)
    
    fig.tight_layout()
    
    pl.savefig("all_step_nclust_%s.pdf"%nclust)


sys.exit()

aa = df.loc[df['step_name'] == 'GC/CG']
X = np.array([map(float,a.split(',')) for a in aa.origin+','+aa.M_axis])
#M = np.array([[mydist(u,v) for u in X] for v in X])
#print(M)
#print(M.shape)
print(X.shape)

X2 = np.hstack((X[0:,0:3],X[0:,3:]))
#print X2[0]
#fclust = fclusterdata(X, t=0.5, criterion='distance', metric=mydist)
#kmeans = KMeans(n_clusters=10,random_state=0).fit(X)
#print(kmeans.labels_)
#kmeans = KMeans(n_clusters=10,random_state=0).fit(X2)
#print(kmeans.labels_)
#kmeans = k_means(X2,n_clusters=10,precompute_distances=True)
#print(kmeans[1])
#print(kmeans[2])

#x = range(1,30)
#y = []
#for i in range(1,30):
#    kmeans = k_means(X2,n_clusters=i,precompute_distances=True)
#    y.append(kmeans[2])

#pl.plot(x,y,'-or')
#pl.show()


pca = PCA(n_components=2)
pca.fit(X2)
#print(pca.explained_variance_ratio_)
#print(pca.singular_values_)
#print(pca.components_)
y = pca.fit_transform(X2)
#print(y)

#tsne = TSNE(n_components=3)
#tsne.fit(X2)
#print(tsne.embedding_)
#print(tsne.kl_divergence_)
#y = tsne.fit_transform(X2)
#print(y)

colors = ['red','blue','yellow','green','orange','cyan','skyblue','lime','gold','gray','violet','firebrick','purple','teal']
nclust = 4
kmeans = k_means(X2,n_clusters=nclust,precompute_distances=True)

for i in range(nclust):
    c = np.where(kmeans[1]==i)[0]
    pl.plot(y[c,0],y[c,1],'o',color=colors[i],markersize=12)

pl.show()

#fig = pl.figure()
#ax = fig.add_subplot(111,projection='3d')
#plt.scatter(x=y[:,0],y=y[:,1],zs=y[:,2],zdir='z',c='red')

#pl.show()

#fig = pl.figure()
#ax = fig.add_subplot(111,projection='3d')
#plt.scatter(x=X[:,0],y=X[:,1],zs=X[:,2],zdir='z',c='red')

#pl.show()




#print(X[0])


#pl.plot(range(1,11),pca.singular_values_,'-or')
#pl.show()

#print(kmeans.cluster_centers_)


#
#print(X)
#print(fclust)
#aa['cluster'] = fclust
#
##print(aa)
#print(aa.groupby(['cluster'])['step_id'].count()/len(aa))

#clustering = DBSCAN(eps=1.0, min_samples=5, metric='precomputed').fit_predict(M)
#print(clustering)
#print(clustering.labels_)
#print(clustering.components_)
#print(clustering.core_sample_indices_)
#print(len(clustering.labels_))

