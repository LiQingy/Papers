import numpy as np 
import h5py
# from oPDF.oPDF import *
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import oPDFplus as opd
import random


def main(model,nx):
	fileCM = '/Users/liqy/Documents/data/oPDFnew/cluster/G3X-CM-masscomplete-ALL.txt'
	if model == 'JeansE':
		filefit = '/Users/liqy/Documents/data/oPDF1/JeansE/JeansE_Chi2fit_MCDM_rcin200_TMP_n100000_bin20_boot200.txt'
		ntracer = float(1e5)
	elif model == 'HE':
		filefit = '/Users/liqy/Documents/data/oPDF1/HE/HE_Chi2fit_MChotgas_rcin200_TMP_n100000_bin20_boot200.txt'
		ntracer = float(1e5)
	elif model == 'VirM':
		filefit = '/Users/liqy/Documents/data/oPDF1/VirM/VirM_Chi2fit_MCDM_rcin200_TMP_n1000_bin20_boot200.txt'
		ntracer = float(1000)

	dataCM = np.loadtxt(fileCM)
	datafit = np.loadtxt(filefit)
	
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	plt.tick_params(top = 'on', right = 'on', which='both')
	
	mmtrue = datafit[:,0] / dataCM[:,2]
	cctrue = datafit[:,1] / dataCM[:,6]
	loc = opd.sigma3sel(np.log10(mmtrue),np.log10(cctrue))
	print(loc.shape)

	m0sig = np.std(np.log10(mmtrue[loc]),ddof = 1)
	c0sig = np.std(np.log10(cctrue[loc]),ddof = 1)
	covall = np.cov(np.log10(mmtrue[loc]),np.log10(cctrue[loc]))
	# print("The total scatter is ", covall[0][0]**0.5)

	sigm_sta = m0sig / nx
	sigc_sta = c0sig / nx 
	Neff = ntracer / (nx**2 - 1) 

	#use DM statistic to estimate satellite error
	# sigm_sta = 0.00536 * (1e5/200)**0.5
	# Neff = ntracer / ((m0sig/sigm_sta)**2 - 1)

	print("The total scatter is ", m0sig, c0sig)
	print("The statistic1 error for m and c are", sigm_sta,sigc_sta)
	print("The statistic2 error for m and c are", m0sig/(ntracer/Neff-1)**0.5,c0sig/(ntracer/Neff-1)**0.5)
	print("The systematic error for m and c are", (ntracer/Neff)**0.5*m0sig/nx,(ntracer/Neff)**0.5*c0sig/nx)

	msig = np.zeros(400)
	csig = np.zeros(400)
	for i in range(400):
		locin = random.choices(loc,k=loc.shape[0])
		m_m0 = np.log10(mmtrue[locin])
		c_c0 = np.log10(cctrue[locin])
		cov = np.cov(m_m0,c_c0)
		msig[i] = cov[0,0]**0.5
		csig[i] = cov[1,1]**0.5
	sigmm = np.std(msig)
	sigcc = np.std(csig)

	Neff0 = ntracer / (((m0sig + sigmm)/sigm_sta)**2 - 1)
	Neff1 = ntracer/ (((m0sig - sigmm)/sigm_sta)**2 - 1)
	print("The total scatter variance are ", sigmm, sigcc)

	return Neff,Neff0,Neff1

def Main(model,nx):
	Neff,Neff0,Neff1 = main(model,nx)
	print("The Neff is ", Neff)
	print("The Neff adding variance are ", Neff0, Neff1)
	print("The Neff variance are ", Neff-Neff0, Neff1-Neff)
	

'''
nx0: initial value to the expansion multiple of statistic error 
'''

# Main(model = 'JeansE', nx = 42)
Main(model = 'HE',  nx = 77)
# Main(model = 'VirM',  nx = 4)





