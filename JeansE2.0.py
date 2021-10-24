from astropy import constants as const
import numpy as np 
from iminuit import Minuit
from numpy.linalg import inv
import sys
sys.path.append('/home/qyli/oPDFnew')
from oPDF.oPDF import *
import h5py
import random  

def Mcal(rvir,nbin,xx,vv):
    # rbin = np.linspace(np.log10(rvir)/nbin, np.log10(rvir), nbin+1)[:-1] #the bin number is n
	rbin = np.linspace(np.log10(200),np.log10(rvir),nbin+1)
	rr = np.log10(np.sqrt(np.sum(xx**2,axis = 1)))

	betav = np.zeros(nbin)
	vsigma2 = np.zeros(nbin)
	sigma_r2 = np.zeros(nbin)
	numden = np.zeros(nbin)

	for kk in range(nbin):
		rbin0 = rbin[kk+1]
		rbin1 = rbin[kk]
		loc = np.where((rr > rbin1) & (rr < rbin0))[0]

		vtrue = loc.shape[0] / (4 / 3 * np.pi * ((10**rbin0)**3 - (10**rbin1)**3)) #number density

		rabs = np.sqrt(np.sum(xx[loc]**2, axis = 1))
		v_r = np.sum(xx[loc] * vv[loc], axis = 1) / rabs
		v_phi = (xx[loc,0] * vv[loc,1] - xx[loc,1] * vv[loc,0]) / np.sqrt(xx[loc,0]**2 + xx[loc,1]**2) 
		v_theta = (xx[loc,0]*xx[loc,2]*vv[loc,0] + xx[loc,1]*xx[loc,2]*vv[loc,1] - xx[loc,0]**2*vv[loc,2] - 
					xx[loc,1]**2*vv[loc,2]) / rabs / np.sqrt(xx[loc,0]**2 + xx[loc,1]**2) 

		sigma_vr = np.std(v_r)
		sigma_vphi = np.std(v_phi)
		sigma_vtheta = np.std(v_theta)

		betav[kk] = 1 - (sigma_vphi**2 + sigma_vtheta**2) / 2 / sigma_vr**2
		vsigma2[kk] = vtrue * sigma_vr**2
		sigma_r2[kk] = sigma_vr**2
		numden[kk] = vtrue

	nr = np.log10(rvir-200) / nbin
	rlin = 10**(rbin[:-1] + nr / 2)
	G = const.G.to('kpc3 / (Msun s2)').value  #kpc^3/Msun/s^2
	unit = const.kpc.value / 1e3

	Mpro = np.zeros(nbin)
	slopev = np.gradient(vsigma2[:],np.log10(rlin[:]),edge_order=1)
	slopev = slopev / rlin[:] / np.log(10)
	Mpro[:] = - rlin[:]**2 / G * (1/numden[:]*slopev[:] + 2 * betav[:] * sigma_r2[:] / rlin[:]) / unit**2
	#Mpro[:7] = 0
	return Mpro

def bootstrap(nbin,Nboot,npp,rvir,xx,vv):
	bootdata = np.zeros(shape = (Nboot,nbin))
	for k in range(Nboot):
		Npploc = np.array(random.choices(range(0,npp),k=npp))
		xx0 = xx[Npploc]
		vv0 = vv[Npploc]
		bootdata[k] = Mcal(rvir,nbin,xx0,vv0)
	return bootdata


def cal_matrix(maindata,bootdata,rvir,nbin,Mtrue,ctrue,nclu):
    # the parent sample

	rbin = np.linspace(np.log10(200),np.log10(rvir),nbin+1)
	nr = np.log10(rvir-200) / nbin
	rlin = 10**(rbin[:-1] + nr / 2)

	pdata = maindata / 1e15
	boots = bootdata / 1e15
	ave_boots = np.mean(boots, axis = 0)

	matrix_c0 = np.zeros(shape = (nbin,nbin))
	for j in range(nbin):
		for k in range(nbin):
			matrix_c0[j,k] = np.sum((boots[:,j] - ave_boots[j]) * (boots[:,k] - ave_boots[k])) / boots.shape[0]
	matrix_c = inv(matrix_c0)

	#the Chi square function 
	def minichi2(M,c):
		chi2 = 0
		halofit = Halo(halotype = HaloTypes.TMPMC, TMPid = nclu)
		m0 = 10**M/1e10
		c0 = 10**c
		halofit.set_param([m0,c0])
		modelM = halofit.mass(rlin)/1e5  
		for ii in range(nbin):
			for jj in range(nbin):
				chi2 += (modelM[ii] - pdata[ii]) * matrix_c[ii,jj] * (modelM[jj] - pdata[jj])
		return chi2

	#----------------------------------------------------
	fitmc = Minuit(minichi2, M = Mtrue, c = ctrue)
	fitmc.limits["M"]=(14,17)
	fitmc.limits["c"]=(-1.,2.)
	fitmc.errors["M"]=1
	fitmc.errors["c"]=1
	fitmc.errordef = Minuit.LEAST_SQUARES

	fitmc.migrad()
	fitmc.hesse()

	if fitmc.valid == False:
		return 0,0,0

	fit_M = 10**fitmc.values['M']
	fit_C = 10**fitmc.values['c']

	Mmin = fitmc.values['M'] - 0.02
	Mmax = fitmc.values['M'] + 0.02
	Cmin = fitmc.values['c'] - 0.02
	Cmax = fitmc.values['c'] + 0.02

	fx,fy,fval0 = fitmc.contour('M','c',size = 100,bound = [[Mmin,Mmax],[Cmin,Cmax]])

	return fit_M,fit_C,fval0

def subsample(nbin,data,rvir,npp,Nboot,Mtrue,ctrue,nclu):

	xx = np.array(data['x'][:])
	vv = np.array(data['v'][:])
	data.close()

	rr = np.sqrt(np.sum(xx**2, axis = 1))
	loc200 = np.where(rr > 200)[0]
	N200 = loc200.shape[0]
	print("The total particle number larger than 200 kpc/h is",N200)

	Npploc = np.array(random.sample(list(loc200),npp))
	xx0 = xx[Npploc]
	vv0 = vv[Npploc]

	maindata = Mcal(rvir,nbin,xx0,vv0)
	fit_M = 0

	while fit_M == 0:
		bootdata = bootstrap(nbin,Nboot,npp,rvir,xx0,vv0)
		fit_M,fit_C,fval0 = cal_matrix(maindata,bootdata,rvir,nbin,Mtrue,ctrue,nclu)

	return fit_M,fit_C,fval0,maindata,bootdata


def main(nbin,Nboot,npp):
	filecen = '/home/qyli/oPDFnew/data/cluster/G3X_Mass_snap_128-center-cluster.txt'
	datacen = np.loadtxt(filecen)
	fileCM = '/home/qyli/oPDFnew/data/cluster/G3X-CM-masscomplete-ALL.txt'
	dataCM = np.loadtxt(fileCM)

	mcfit = np.zeros(shape = (324,2))
	mainall = np.zeros(shape = (324,nbin))
	fchi2 = 0

	for i in range(324):
		rvir = datacen[i][6]
		Mtrue = np.log10(dataCM[i,2])
		ctrue = np.log10(dataCM[i,6])

		dataf = '/home/qyli/oPDFnew/data/inidata/GXsub_DM/GXsub_DM_vx%s.h5' %(i+1)
		data = h5py.File(dataf,'r')

		mcfit[i,0], mcfit[i,1], fval0, mainall[i], bootdata = subsample(nbin, data, rvir, npp, Nboot, Mtrue,ctrue ,i ) 
		fchi2 += fval0.T - np.min(fval0)

		np.savetxt('/home/qyli/oPDFnew/data/JeansE/bootstrap_DM_n%s_bin%s_boot%s/GXsub_DM_JEboot%s.txt' %(npp,nbin,Nboot,i+1), bootdata)

		print(i)
	fchi2 = fchi2 / 324

	np.savetxt('/home/qyli/oPDFnew/data/JeansE/JeansE_GXsub_DM_Mpro_n%s_bin%s_boot%s.txt' %(npp,nbin,Nboot),mainall)
	np.savetxt('/home/qyli/oPDFnew/data/JeansE/JeansE_Chi2_rcin200_TMP_n%s_bin%s_boot%s.txt' %(npp,nbin,Nboot), fchi2)
	np.savetxt('/home/qyli/oPDFnew/data/JeansE/JeansE_Chi2fit_MCDM_rcin200_TMP_n%s_bin%s_boot%s.txt' %(npp,nbin,Nboot), mcfit)
main(nbin = 20, Nboot = 200,npp = 100000)













