import numpy as np 
import h5py
# from oPDF.oPDF import *
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import oPDFplus as opd
import random


def main(pptype,nx):
	if pptype == 'DM':
		datadir = '/Users/liqy/Documents/data/oPDF1/WentingFits/'
		MRII=np.loadtxt(datadir+'MRIIfits.dat') #isolated
		MRIIb=np.loadtxt(datadir+'MRIIfitsbinary.dat') #binary
		MRall = np.append(MRII,MRIIb[:,:8],axis = 0)
		MR=np.vstack([MRII[:,2:4], MRIIb[:,2:4]]) #combine
		# MR=MRII[:,2:4]
		print(MR.shape)
		sel=(MR[:,0]>0.1)&(MR[:,1]>0.1)&(MR[:,0]<10)&(MR[:,1]<10)

		mmtrue = MR[:,0]
		cctrue = MR[:,1]
		loc = sel
		ntracer = 1e5

	elif pptype == 'star':
		filefit = '/Users/liqy/Documents/data/oPDF1/WentingFits/SVfof_mcfit_4_Wang.txt'
		datafit = np.loadtxt(filefit)
		mmtrue = 10**datafit[:,0]
		cctrue = 10**datafit[:,1]
		loc = np.arange(0,24,1)
		ntracer = 4.5e4

	else:
		f1 = '/Users/liqy/Documents/data/oPDF1/ForQiangyang/MRII-HBTSat-20-400kpc-Rbin15.np.npy'
		d1 = np.load(f1)
		loc = np.where((10 > d1[:,0]) & (d1[:,0] > 1e-1) & (10 > d1[:,1]) & (d1[:,1] > 1e-1) & (d1[:,-1] > 0))[0]
		mmtrue = d1[:,0]
		cctrue = d1[:,1]
		ntracer = 500

	m0sig = np.std(np.log10(mmtrue[loc]),ddof = 1)
	c0sig = np.std(np.log10(cctrue[loc]),ddof = 1)
	covall = np.cov(np.log10(mmtrue[loc]),np.log10(cctrue[loc]))

	sigm_sta = m0sig / nx
	sigc_sta = c0sig / nx 
	Neff = ntracer / (nx**2 - 1)
	print(c0sig)

	if pptype != 'ste':
		print("The total scatter is ", m0sig, c0sig)
		print("The statistic error for m and c are", sigm_sta,sigc_sta)
		print("The systematic error for m and c are", (ntracer/Neff)**0.5*m0sig/nx,(ntracer/Neff)**0.5*c0sig/nx)

	if pptype == 'ste':
		sigm_sta = m0sig / nx * (10000/500)**0.5  
		sigc_sta = c0sig / nx * (10000/500)**0.5  
		Neff = ntracer / (nx**2/20 - 1)
		
		print("The total scatter is ", m0sig, c0sig)
		print("The statistic error for m and c are", sigm_sta,sigc_sta)
		print("The systematic error for m and c are", (ntracer/Neff)**0.5*m0sig/nx*20**0.5,
			(ntracer/Neff)**0.5*c0sig/nx*20**0.5)

	#use DM data


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

	print("The total scatter variance are ", sigmm, sigcc)

	Neff0 = ntracer / (((m0sig + sigmm)/sigm_sta)**2 - 1)
	Neff1 = ntracer/ (((m0sig - sigmm)/sigm_sta)**2 - 1)

	return Neff,Neff0,Neff1

def Main(pptype,nx):
	Neff,Neff0,Neff1 = main(pptype,nx)
	print("The Neff is ", Neff)
	print("The Neff adding variance are ", Neff0, Neff1)
	print("The Neff variance are ", Neff-Neff0, Neff1-Neff)
	

'''
nx0: initial value to the expansion multiple of statistic error 
'''
print("MW halo")
Main(pptype = 'DM', nx = 10)
Main(pptype = 'star', nx = 1/0.03)
Main(pptype = 'ste', nx = 5.7)







