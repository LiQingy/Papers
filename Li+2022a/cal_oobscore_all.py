import numpy as np 
import h5py
import matplotlib.pyplot as plt

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split

from scipy.special import comb
from itertools import combinations

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def unify_data(inputdata, untype = 'MM'):
	dmin = np.min(inputdata, axis = 0)
	dmax = np.max(inputdata, axis = 0)

	if untype == 'MM':
		return (inputdata - dmin) / (dmax - dmin)


def cal_nodes(ntotal,rank0,size0):
    if size0 >= ntotal:
        if rank0 >= ntotal:
            return 0,0
        else:
            return rank0,rank0+1

    remainder = ntotal%size0
    neach = int(ntotal/size0)

    
    if rank0 < remainder:
        ibegin = (neach+1)*rank0
        iend = ibegin + neach + 1
    else:
        ibegin = (neach+1)*remainder+(rank0-remainder)*neach
        iend = ibegin + neach
    # print(int(ibegin),int(iend))
    return int(ibegin),int(iend)

#find more combination 
#nfea: choose n features
def oob_n(feature, target, totfea, nfea, ibegin, iend, nerr):
      
	nper = int(comb(totfea,nfea)) #permutation number
	score = np.zeros(shape = (iend - ibegin, nerr), dtype = np.float64)
	per = np.array(list(combinations(np.arange(0,totfea,1), nfea))) #permutation

	for i in range(ibegin,iend):
		feature0 = feature[:,per[i]]
		for j in range(nerr):
			x_train,x_test,y_train,y_test = train_test_split(feature0, target, test_size=0.3, random_state=j)

			forest = RandomForestRegressor(n_estimators=400, min_samples_leaf = 5, random_state=j, oob_score = True)
			forest.fit(x_train, y_train)
			score[i-ibegin, j] = forest.oob_score_

	return score, per

#=================
#nfea: selected feature number
#model: oPDF or Jeans
#totfea: total feature choose
#error: subsample times
#=================

def main(feature, target, model, nfea, totfea, nerr):
	ntotal = int(comb(totfea,nfea))
	ibegin, iend = cal_nodes(ntotal,rank,size)

	score, per = oob_n(feature, target, totfea, nfea, ibegin, iend, nerr)
	caln = iend - ibegin

	calnall = comm.gather(caln, root = 0)
	if rank != 0:
		comm.Send([score,MPI.FLOAT], dest = 0, tag = rank)
		
	else:
		scoreall = np.zeros((ntotal,nerr),dtype = np.float64)
		i0 = 0
		for i in range(size):
			if i != 0:
				caln = calnall[i]
				score = np.empty(shape = (caln,nerr), dtype = np.float64)
				comm.Recv([score, MPI.FLOAT],source = i, tag = i)

			scoreall[i0:(i0+caln)] = score
			i0 += caln
		# save data with oob sorted
		print(np.where(scoreall == 0))
		np.savetxt('./data/oob_%s/%s_oob_score%s' %(model, model, nfea),scoreall, fmt = '%f')
		np.savetxt('./data/oob_%s/%s_oob_per%s' %(model, model, nfea),per, fmt = '%d') #save for No. proxy

import sys
sys.path.append('/home/qyli/oPDF2/')

fileCM = './proxy/cluster/G3X-CM-masscomplete-ALL.txt'
filefit = './proxy/cluster/GXsub_DM_fitmc_rcin200_100000_TMP.txt'
dataCM = np.loadtxt(fileCM)
datafit = np.loadtxt(filefit)
import oPDFplus as opd
mmtrue = datafit[:,0] / dataCM[:,2] * 1e15
cctrue = datafit[:,1] / dataCM[:,6]
loc = opd.sigma3sel(np.log10(mmtrue),np.log10(cctrue))

target = np.loadtxt('./proxy/cluster/GXsub_DM_InL_rcin200_100000_TMP.txt')[loc]
#target = unify_data(np.log10(target))
#target = unify_data(target)
maind = np.loadtxt('./proxy/GadgetX_DS_main_v2')

#feathers for DS
#mar = maind[:,6]
#maind[mar<0,6] = np.min(mar[mar>0])
#maind[:,6] = np.log10(maind[:,6]) #make mar to be log
#maind = unify_data(maind[loc,1:])
maind = maind[loc,1:]

#features for Xray
Xrayd = np.loadtxt('./proxy/GadgetX_DS_Xray')
#Xrayd = unify_data(Xrayd[loc,1:])
Xrayd = Xrayd[loc,1:]

#features for SZ
SZd = np.loadtxt('./proxy/GadgetX_DS_SZ')
#SZd = unify_data(SZd[loc,1:])
SZd = SZd[loc,1:]

#features for Offset
Offsetd = np.loadtxt('./proxy/GadgetX_DS_Offset')
#Offsetd = unify_data(Offsetd[loc,1:])
Offsetd = Offsetd[loc,1:]
print('begin calculation')
for k in range(1,7):	
	main(Xrayd, target, model = 'Xray',  nfea = k, totfea = Xrayd.shape[1], nerr = 50)
	main(SZd, target, model = 'SZ',  nfea = k, totfea = SZd.shape[1], nerr = 50)
	main(Offsetd, target, model = 'Offset',  nfea = k, totfea = Offsetd.shape[1], nerr = 50)
	main(maind, target, model = 'oPDF',  nfea = k, totfea = maind.shape[1], nerr = 50)



