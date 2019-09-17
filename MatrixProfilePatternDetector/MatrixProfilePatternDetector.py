import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
import stumpy
import pickle
#from MatrixProfilePatternDetector import MatrixProfilePatternDetector as mppd
from scipy import signal
from matrixprofile import *
import seaborn as sns
import scipy as scipy
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import pyplot, transforms
from scipy import ndimage
from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import KMeans

def butt_filter(T, order, freq, btype = "low"): #Butterworth filter
    """ 
    T - timeseries,
    order - the order of the filter.
    freq - cut-off frequency of the filter
    btype - the type of filter
    """
    fs = len(T)
    fc = freq  
    w = fc / (fs / 2)
    b, a = signal.butter(order, w, btype)
    output = signal.filtfilt(b, a, T)
    return output

def plot_motifs(mtfs, labels, ax):

    colori = 0
    colors = 'rgbcm'
    for ms,l in zip(mtfs,labels):
        c =colors[colori % len(colors)]
        starts = list(ms)
        ends = [min(s + m,len(pattern)-1) for s in starts]
        ax.plot(starts, pattern[starts],  c +'o',  label=l)
        ax.plot(ends, pattern[ends],  c +'o', markerfacecolor='none')
        for nn in ms:
            ax.plot(range(nn,nn+m),pattern[nn:nn+m], c , linewidth=2)
        colori += 1

    ax.plot(pattern, 'k', linewidth=1, label="data")
  



scaler = MinMaxScaler()
data = pd.read_csv("gas_turbine.csv")
t1 = data['AP'].values

matrix_pofile = pd.DataFrame()
ft1 = t1 - butt_filter(t1, order = 5, freq = 3500, btype = 'low')


print("Length of timeseries -", len(ft1))

pattern = ft1[:250]
x = []
for point in range(len(pattern)):
    x.append(point)

m = 10
print('Start')


for i in range(len(pattern)-m+1):
    query = pattern[i:(i+m)]
    distanceProfile = []
    for j in range(len(pattern)-m+1):
        
        distanceProfile.append(scipy.spatial.distance.euclidean(query,pattern[j:(j+m)]))

    dp = np.array(distanceProfile)
   
    
    
    matrix_pofile = pd.concat([matrix_pofile,pd.DataFrame(dp)], axis=1)











mp_t = pd.DataFrame(np.transpose(matrix_pofile))
mp_t = pd.DataFrame(scaler.fit_transform(mp_t))
tuples = tuple(mp_t.values[0:])
one_tuple = tuple(tuples)
print(one_tuple[0])
kmeans = KMeans(n_clusters=5, random_state=0)
mp = matrixProfile.stomp(pattern,m)
mtfs ,motif_d  = motifs.motifs(pattern, tuples, max_motifs=10)




matrix_of_dist = pd.DataFrame()
dist_of_dist=[]
for i in range(len(pattern)-m+1):
    dist_of_dist=[]
    for j in range(len(pattern)-m+1):
         diff = abs(mp_t.iloc[i,:].values - mp_t.iloc[j,:].values)
         dist_of_dist.append(sum(diff)/(len(pattern)-m))
    dp = np.array(dist_of_dist)
    matrix_of_dist = pd.concat([matrix_of_dist,pd.DataFrame(dp)], axis=1)

matrix_of_dist_t =  pd.DataFrame(np.transpose(matrix_of_dist))
sns.heatmap(matrix_of_dist_t, cmap="YlGnBu",xticklabels = 10, yticklabels=10, square=True)

    
plt.show()


index_of_common = np.where((matrix_of_dist_t.values < 0.08) & (matrix_of_dist_t.values != 0))
starts_y = index_of_common[0]
print(starts_y)
#starts_x = index_of_common[1]
#print(starts_x)
#kmeans.fit(mp_t)
#labels = kmeans.labels_
#print(len(labels))

#new_matrix = matrix_of_dist_t.values



fig, axScatter = plt.subplots(figsize=(5.5, 5.5))
axScatter = sns.heatmap(mp_t, cmap="YlGnBu",xticklabels = 10, yticklabels=10, square=True)
divider = make_axes_locatable(axScatter)
axHistx = divider.append_axes("top", 1.2, pad=0.1, sharex=axScatter)
axHisty = divider.append_axes("left", 1.2, pad=0.1, sharey=axScatter)
axHistx.plot(x, pattern)
axHisty.plot(-1*pattern, x)
print(mtfs)
plot_motifs(mtfs, [f"{md:.3f}" for md in motif_d], axHistx)
#ax3.set_ylabel('Motifs', size=22)
#colors = 'rgbcmk'
#for i in range(len(pattern)-m+1):
 #   axHistx.plot(x[i:i+m], pattern[i:i+m],c=colors[labels[i]])

plt.show()



colors = 'rgbcm'
c = []
i = 0
prev = -1
for e in starts_y:
    if e == prev:
        c.append(colors[i % len(colors)])
    else:
        i += 1
        c.append(colors[i % len(colors)])


for j in range(len(starts_y)):
    col = c[j]
    e = starts_y[j]
    #ex = starts_x[j]
    axHistx.plot(x[e:e+m], pattern[e:e+m],c='r')
    #axHistx.plot(x[ex:ex+m], pattern[ex:ex+m],c=col)
        
    

   
 
plt.show()





print(len(dp[0]))


def plot_motifs(mtfs, labels, ax):

    colori = 0
    colors = 'rgbcm'
    for ms,l in zip(mtfs,labels):
        c =colors[colori % len(colors)]
        starts = list(ms)
        ends = [min(s + m,len(pattern)-1) for s in starts]
        ax.plot(starts, pattern[starts],  c +'o',  label=l)
        ax.plot(ends, pattern[ends],  c +'o', markerfacecolor='none')
        for nn in ms:
            ax.plot(range(nn,nn+m),pattern[nn:nn+m], c , linewidth=2)
        colori += 1

    ax.plot(pattern, 'k', linewidth=1, label="data")
    ax.legend()


mtfs ,motif_d  = motifs.motifs(pattern, mp, max_motifs=10)
#Append np.nan to Matrix profile to enable plotting against raw data
mp_adj = np.append(mp[0],np.zeros(m-1)+np.nan)

#Plot the signal data
fig, (ax1, ax2, ax3) = plt.subplots(3,1,sharex=True,figsize=(20,10))
ax1.plot(np.arange(len(pattern)),pattern, label="Synthetic Data")
ax1.set_ylabel('Signal', size=22)

#Plot the Matrix Profile
ax2.plot(np.arange(len(mp_adj)),mp_adj, label="Matrix Profile", color='red')
ax2.set_ylabel('Matrix Profile', size=22)

#Plot the Motifs
plot_motifs(mtfs, [f"{md:.3f}" for md in motif_d], ax3)
ax3.set_ylabel('Motifs', size=22)
#plt.xlim((0,100))
plt.show()



fig, (ax1, ax2, ax3) = plt.subplots(3,1,sharex=True,figsize=(20,10))


mtfs ,motif_d  = motifs.motifs(pattern, mp, max_motifs=5, n_neighbors=4)
plot_motifs(mtfs, [f"{md:.3f}" for md in motif_d], ax1)
ax1.set_ylabel('4 Neigbhors', size=22)

mtfs ,motif_d  = motifs.motifs(pattern, mp, max_motifs=5, radius=10)
plot_motifs(mtfs, [f"{md:.3f}" for md in motif_d], ax2)
ax2.set_ylabel('Radius = 10', size=22)

mtfs ,motif_d  = motifs.motifs(pattern, mp, max_motifs=5, ex_zone=2*m)
plot_motifs(mtfs, [f"{md:.3f}" for md in motif_d], ax3)
ax3.set_ylabel('Exclude 2*m', size=22)
plt.show()

#cac = fluss.fluss(mp[1], m)
#mp_adj = np.append(mp[0],np.zeros(m-1)+np.nan)

#Plot the signal data
#fig, (ax1, ax2, ax3) = plt.subplots(3,1,sharex=True,figsize=(20,10))
#ax1.plot(np.arange(len(pattern)),pattern, label="Synthetic Data")
#ax1.set_ylabel('Signal', size=22)

#Plot the Matrix Profile
#ax2.plot(np.arange(len(mp_adj)),mp_adj, label="Matrix Profile", color='red')
#ax2.set_ylabel('Matrix Profile', size=22)
#ax2.set_xlabel('Sample', size=22)

#Plot the CAC
#ax3.plot(np.arange(len(cac)),cac, label="CAC", color='green')
#ax3.set_ylabel('CAC', size=22)
#ax3.set_xlabel('Sample', size=22)




plt.show()