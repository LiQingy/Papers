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


def cal_onefeature(fea0,para_leaf):

	oob = np.zeros(50)
	for i in range(50):
		x_train,x_test,y_train,y_test = train_test_split(fea0.reshape(310,1), target, test_size=0.3, random_state=i)
		forest = RandomForestRegressor(n_estimators=400, min_samples_leaf = para_leaf, random_state=i,  oob_score = True)
		forest.fit(x_train, y_train)
		oob[i] = forest.oob_score_
	       
	return np.mean(oob),np.std(oob)

def main(dpy):

	if dpy == 'D3':
		data = maind
		# para_leaf = [14,60,40,45,60,30,24,40,60,68,68,68,75,75,20,50,68,68,40,60,75,32,68,68,68,68]
		para_leaf = [14,60,62,45,60, 30,24,55,45,55, 68,50,75,70,70, 50,68,68,40,60, 75,60,68,68,68, 68] #To be updated
	elif dpy == 'Xray':
		data = Xrayd
		# para_leaf = [68,30,40,68,68,40,42]
		para_leaf = [50,30,40,60,65,40,42] #To be updated
	elif dpy == 'SZ':
		data = SZd
		# para_leaf = [40,56,68,68,50,68,68]
		para_leaf = [58,56,65,40,50,40,64] #To be updated
	else:
		data = Offsetd
		# para_leaf = [22,6,40,10,12,8,36,30,34,68,24,12,44,42,68]
		para_leaf = [68,20,40,25,40, 30,36,30,34,68, 24,12,44,42,68] #To be updated


	score_one = np.zeros((data.shape[1],2))
	for i in range(data.shape[1]):
		score_one[i,0], score_one[i,1] = cal_onefeature(data[:,i], para_leaf[i])
		print(i)
	np.savetxt('./data_nRandom/%s_one_OOB' %dpy, score_one)

main('D3')
main('Xray')
main('SZ')
main('Offset')













