import numpy as np 
from scipy.optimize import curve_fit
from iminuit import Minuit
from numpy.linalg import inv,pinv
from oPDF.oPDF import *
import h5py
#import scipy
from numpy import matrix

# the model of TMP
def tmp(xx,M,c,nclu):
	# halofit = Halo(halotype=HaloTypes.NFWMC)
	halofit = Halo(halotype=HaloTypes.TMPMC, TMPid = nclu)
	m0 = 10**M/1e10
	c0 = 10**c
	halofit.set_param([m0,c0])
	return halofit.mass(xx)*1e10

#calculate covariance matrix from bootstrap sample
def calcov(i,ncut,dn,ptype):
	#import bootstrap sample. 
	#Here I drop out the inner 7 bins because the bins are within softning length. 
	#so the shape of matrix is 18 * 18

	fboots = '/home/qyli/oPDFnew/data/VirM/subsample_%s_l200n5000/GXsub_%s_VMboot%s.txt' %(ptype,ptype,i+1)
	boots = np.loadtxt(fboots)[:,ncut:]

	#compute average value
	ave_boots = np.mean(boots, axis = 0)

	#calculate covariance matrix
	matrix_c0 = np.zeros(shape = (dn,dn))
	for j in range(dn):
		for k in range(dn):
			matrix_c0[j,k] = np.sum((boots[:,j] - ave_boots[j]) * (boots[:,k] - ave_boots[k])) / boots.shape[0]
	#get inverse matrix
	matrix_c = inv(matrix_c0)
	print(np.allclose(np.dot(matrix_c,matrix_c0), np.eye(dn)))
	return ave_boots,matrix_c

#fit M and c
def mainfit(dataCM,nbin,rcut,ptype,npp):

	# the parent sample
	filemain = '/home/qyli/oPDFnew/data/VirM/VirM_GXsub_%s_Mpro_l200n%s.txt' %(ptype,npp)
	#filemain = '/home/qyli/VirM_GXsub_%s_Mpro_%s_test.txt' %(ptype,npp)
	maindata = np.loadtxt(filemain)
	fileCM = '/home/qyli/oPDFnew/data/cluster/G3X-CM-masscomplete-ALL.txt'
	dataCM = np.loadtxt(fileCM)

	mainMC = np.zeros(shape = (324,2))
	fvalall = np.zeros(shape = (100,100))
	for i in range(324): #compute for 324 clusters
		# i = 55
		rvir = dataCM[i][3]
		rbin = np.linspace(np.log10(200),np.log10(rvir),nbin+1)
		nr = np.log10(rvir-200) / nbin
		rlin = 10**(rbin[:-1] + nr / 2)
		
		#locrcut = np.where(rlin > rcut)[0]
		#ncut = int(locrcut[0])
		ncut = 0
		#dn = locrcut.shape[0]
		dn = nbin	
		Mtrue = np.log10(dataCM[i][2])
		ctrue = np.log10(dataCM[i,6])
		if i == 293:
			continue
	
		ave_boots, matrix_c = calcov(i,ncut,dn,ptype)
		pdata = maindata[i][ncut:]
		#pdata = ave_boots

		#the Chi square function 
		def minichi2(M,c):
			chi2 = 0
			modelM = tmp(rlin[ncut:],M,c,i) #TMP model value
			
			#compute Chi2
			for ii in range(dn):
				for jj in range(dn):
					chi2 += (modelM[ii] - pdata[ii]) * matrix_c[ii,jj] * (modelM[jj] - pdata[jj])
			return chi2

		fitmc = Minuit(minichi2, M = Mtrue, c=ctrue)
		fitmc.limits["M"]=(14,16)
		fitmc.limits["c"]=(-1.,2.)
		fitmc.errors["M"]=1
		fitmc.errors["c"]=0.1
		fitmc.errordef = Minuit.LEAST_SQUARES

		fitmc.migrad()
		fitmc.hesse()
		if i != 46 and i != 231:
			fitmc.minos()

		mainMC[i][0] = 10**fitmc.values['M']
		mainMC[i][1] = 10**fitmc.values['c']
		# print(mainMC[i][0]-10**Mtrue,mainMC[i][1]-10**ctrue)
		print(mainMC[i])
		print(10**Mtrue,10**ctrue)

		Mmin = fitmc.values['M'] - 0.02
		Mmax = fitmc.values['M'] + 0.02
		Cmin = fitmc.values['c'] - 0.02
		Cmax = fitmc.values['c'] + 0.02

		fx,fy,fval0 = fitmc.contour('M','c',size = 100,bound = [[Mmin,Mmax],[Cmin,Cmax]])
		if (not i in [ 4,  13,  46,  67, 207, 234] ): #3 sigma cut 
			fvalall += fval0.T - np.min(fval0)
		#print(fval0.T-np.min(fval0))
		# if i == 55 or i == 120 or i == 240 or i == 110 or i == 314 or i == 315:
		# 	con1 = fitmc.mncontour('M', 'c', size = 100)
		# 	np.savetxt('/home/qyli/JeansE_con%s.txt' %(i+1),con1)

		# 	fx,fy,fval0 = fitmc.contour('M','c',size = 100,bound = [[Mmin,Mmax],[Cmin,Cmax]])
		# 	conf = h5py.File('/home/qyli/JeansE_conv%s.hdf5' %(i+1),'w')
		# 	conf['x'] = fx
		# 	conf['y'] = fy
		# 	conf['fval'] = fval0
		# 	conf.close()
		print(i)

	np.savetxt('/home/qyli/oPDFnew/data/VirM/VirM_Chi2_l200n5000.txt', fvalall/317)
	return mainMC

def MC(nbin, ptype, rcut, npp):

	fileCM = '/home/qyli/oPDFnew/data/cluster/G3X-CM-masscomplete-ALL.txt'
	dataCM = np.loadtxt(fileCM) #The true cluster information
	mainMC = mainfit(dataCM,nbin,rcut,ptype,npp) #get M and c
	
	np.savetxt('/home/qyli/oPDFnew/data/VirM/VirM_Chi2fit_MC%s_rcin200_TMP_l200n5000.txt' %ptype, mainMC)

def gather_data(nbin):

	mainMpro = np.zeros(shape = (324,nbin))
	file1 = '/home/qyli/oPDFnew/data/VirM/VirM_GXsub_DM_Mpro_5000_l200n5000_1.txt'
	file2 = '/home/qyli/oPDFnew/data/VirM/VirM_GXsub_DM_Mpro_5000_l200n5000_2.txt'
	file3 = '/home/qyli/oPDFnew/data/VirM/VirM_GXsub_DM_Mpro_5000_l200n5000_3.txt'
	d1 = np.loadtxt(file1)
	d2 = np.loadtxt(file2)
	d3 = np.loadtxt(file3)

	mainMpro[:101] = d1[:101]
	mainMpro[101:201] = d2[101:201]
	mainMpro[201:] = d3[201:]

	np.savetxt('/home/qyli/oPDFnew/data/VirM/VirM_GXsub_DM_Mpro_l200n5000.txt', mainMpro)

gather_data(nbin = 20)
MC(nbin = 20, ptype = 'DM', rcut = 200, npp = 5000)

