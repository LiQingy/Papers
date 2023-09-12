import numpy as np 
import h5py
import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split

import sys
sys.path.append('/home/qyli/oPDF2/')

fileCM = '../proxy/cluster/G3X-CM-masscomplete-ALL.txt'
filefit = '../proxy/cluster/GXsub_DM_fitmc_rcin200_100000_TMP.txt'
dataCM = np.loadtxt(fileCM)
datafit = np.loadtxt(filefit)
import oPDFplus as opd
mmtrue = datafit[:,0] / dataCM[:,2] * 1e15
cctrue = datafit[:,1] / dataCM[:,6]
loc = opd.sigma3sel(np.log10(mmtrue),np.log10(cctrue))

target = np.loadtxt('../proxy/cluster/GXsub_DM_InL_rcin200_100000_TMP.txt')[loc]

maind = np.loadtxt('../proxy/GadgetX_DS_main_v2')
maind = maind[loc,1:]

#features for Xray
Xrayd = np.loadtxt('../proxy/GadgetX_DS_Xray')
Xrayd = Xrayd[loc,1:]

#features for SZ
SZd = np.loadtxt('../proxy/GadgetX_DS_SZ')
SZd = SZd[loc,1:]

#features for Offset
Offsetd = np.loadtxt('../proxy/GadgetX_DS_Offset')
Offsetd = Offsetd[loc,1:]


def cal_onefeature(fea0):
    k = 0
    
    oob_leaf = np.zeros(25)
    for i in range(25):
        x_train,x_test,y_train,y_test = train_test_split(fea0.reshape(310,1), target, test_size=0.3, random_state=k)
        forest = RandomForestRegressor(n_estimators=400, min_samples_leaf = i*2+2, random_state=k,  oob_score = True)
        forest.fit(x_train, y_train)
        oob_leaf[i] = forest.oob_score_
        
    oob_split = np.zeros(25)
    for i in range(25):
        x_train,x_test,y_train,y_test = train_test_split(fea0.reshape(310,1), target, test_size=0.3, random_state=k)
        forest = RandomForestRegressor(n_estimators=400, min_samples_split = i*2+2, random_state=k,  oob_score = True)
        forest.fit(x_train, y_train)
        oob_split[i] = forest.oob_score_
        
    oob_depth = np.zeros(20)  
    for i in range(20):
        x_train,x_test,y_train,y_test = train_test_split(fea0.reshape(310,1), target, test_size=0.3, random_state=k)
        forest = RandomForestRegressor(n_estimators=400, max_depth = i+1, random_state=k,  oob_score = True)
        forest.fit(x_train, y_train)
        oob_depth[i] = forest.oob_score_
        
    oob_ntrees = np.zeros(20)  
    ntrees = np.arange(20,420,20,dtype = np.int64)
    for i in range(20):
        x_train,x_test,y_train,y_test = train_test_split(fea0.reshape(310,1), target, test_size=0.3, random_state=k)
        forest = RandomForestRegressor(n_estimators=ntrees[i], random_state=k,  oob_score = True)
        forest.fit(x_train, y_train)
        oob_ntrees[i] = forest.oob_score_
           
    return oob_leaf, oob_split, oob_depth, oob_ntrees

def plot_figures(dpy, nfea, oob_leaf, oob_split, oob_depth, oob_ntrees):
	plt.figure(figsize = (8,7))

	#---min_samples_leaf---
	plt.subplot(221)
	xx = np.arange(50,step = 2)+2
	plt.plot(xx, oob_leaf, 'o')
	plt.xlabel('min_samples_leaf')
	plt.ylabel('OOB score')
	plt.grid(c = 'grey', alpha = 0.6, ls = '--')

	#---min_samples_split---
	plt.subplot(222)
	xx = np.arange(50, step = 2)+2
	plt.plot(xx, oob_split, 'o')
	plt.xlabel('min_samples_split')
	plt.ylabel('OOB score')
	plt.grid(c = 'grey', alpha = 0.6, ls = '--')

	#---max_depth---
	plt.subplot(223)
	xx = np.arange(20)+1
	plt.plot(xx, oob_depth, 'o')
	plt.xlabel('max_depth')
	plt.ylabel('OOB score')
	plt.grid(c = 'grey', alpha = 0.6, ls = '--')

	#---min_samples_leaf---
	plt.subplot(224)
	xx = np.arange(20,420,20,dtype = np.int64)
	plt.plot(xx, oob_ntrees, 'o')
	plt.xlabel('n_estimators')
	plt.ylabel('OOB score')
	plt.grid(c = 'grey', alpha = 0.6, ls = '--')

	plt.tight_layout()

	plt.savefig('./%s/%s.png' %(dpy,nfea+1))

def main(dpy):

	if dpy == 'D3':
		data = maind
	elif dpy == 'Xray':
		data = Xrayd
	elif dpy == 'SZ':
		data = SZd
	else:
		data = Offsetd 


	for i in range(data.shape[1]):
		oob_leaf, oob_split, oob_depth, oob_ntrees = cal_onefeature(data[:,i])
		plot_figures(dpy, i, oob_leaf, oob_split, oob_depth, oob_ntrees)
		print(i)

main('D3')
main('Xray')
main('SZ')
main('Offset')













