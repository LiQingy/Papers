import numpy as np 
import h5py
import os
from oPDF.oPDF import *

def oPDFfit(filevx,ppn,fitmodel,dx0,rin_cut,rout_cut,sigcont,i,MCtrue):
	Globals.cosmology.OmegaM0=0.307
	Globals.cosmology.OmegaL0=0.693
	filevx = filevx.encode('utf8')

	if rin_cut == 1 or rin_cut == 0:
		Sample = Tracer(filevx, rmax = rout_cut)
		Sample = Sample.copy(0,ppn)
	else:
		Sample = Tracer(filevx, rmin = rin_cut, rmax = rout_cut)
		Sample = Sample.copy(0,ppn)

	if fitmodel == 'NFW':
		Sample.halo.set_type(HaloTypes.NFWMC,scales = [1e5, 1])
	elif fitmodel == 'TMP':
		Sample.halo.set_type(HaloTypes.TMPMC,scales = [1e5, 1], TMPid = i)
	
	Estimators.RBinLike.nbin=10
	mc,fval,status = Sample.dyn_fit(Estimators.RBinLike, x0 = MCtrue)
	f0=Sample.likelihood(MCtrue, Estimators.RBinLike)
	InLi=(fval-f0)
	print(fval-f0)

	if sigcont:
		mline,cline,sig,like=Sample.scan_confidence(Estimators.RBinLike, mc, ngrids=[20,20], dx=dx0, logscale=True, maxlike=fval)
		return mc,sig,like,fval,InLi
	else:
		return mc,InLi


def main(model,pp,ppn,fitmodel,rin,rout,dx0,sigcont,InL):

	if model == 'GXsub' or model == 'GX':
		fileCM = '/home/qyli/oPDF/data/cluster/G3X-CM-masscomplete-ALL.txt'
		filecen = '/home/qyli/oPDF/data/cluster/G3X_Mass_snap_128-center-cluster.txt'
	
	if model == 'MDPL2sub':
		fileCM = '/home/qyli/oPDF/data/cluster/MDPL2-CM-masscomplete-ALL.txt'
		filecen = '/home/qyli/oPDF/data/cluster/MDPL2_Mass_snap_128-center-cluster.txt'

	if model == 'GXall': # all galaxy
		filengal = '/home/qyli/oPDF/data/cluster/GX_galall_ngal.txt'
	else: # galaxy with stellar mass cut and within r200
		filengal = '/home/qyli/oPDF/data/cluster/GX_ste_truen.txt'
	
	datangal = np.loadtxt(filengal)
	dataCM = np.loadtxt(fileCM)
	datacen = np.loadtxt(filecen)
	
	Chi2 = 0

	mcfit = np.zeros(shape = (324,2))
	InLvalue = np.zeros(324)
	ppntemp = 0

	for i in range(324):
		#read snapshots data

		if model == 'MDPL2sub':
			filevx = '/home/qyli/oPDFnew/data/inidata/%s/%s_DM_vx%s.h5' %(model,model,i+1)
		else:
			filevx = '/home/qyli/oPDFnew/data/inidata/%s_%s/%s_%s_vx%s.h5' %(model,pp,model,pp,i+1)

		if rin < 1.5:
			rin_cut = rin * dataCM[i][3]
		else:
			rin_cut = rin

		rout_cut = dataCM[i][3] * rout
		MCtrue = [dataCM[i, 2] / 1e15, dataCM[i, 6]]
		#MCtrue = [1, 2]

		if ppn == 'nsame':
			ppntemp = datangal[i]
		else:
			ppntemp = int(ppn * 1)

		#if i == 46:
		#	mcfit[i] = 0
		#	continue

		if sigcont == True:
			mc,sig,like,fval,InLi=oPDFfit(filevx,ppntemp,fitmodel,dx0,rin_cut,rout_cut,sigcont,i,MCtrue)
			Chi2 += 2*(fval - like)
		else:	
			mc,InLi = oPDFfit(filevx,ppntemp,fitmodel,dx0,rin_cut,rout_cut,sigcont,i,MCtrue)

		mcfit[i] = mc
		InLvalue[i] = InLi
		print(mc)
		np.savetxt('/home/qyli/oPDFnew/data/oPDF_stat/%s_%s/%s_%s_stat%s' %(pp,ppn,model,pp,i+1), Chi2)

	# average Chi2
	if sigcont == True:
		Chi2 = Chi2 / 324
		np.savetxt('/home/qyli/oPDFnew/data/Trueini/%s_%s_Chi2_rcin%s_%s_%s.txt' %(model,pp,rin,ppn,fitmodel),Chi2)
	if InL == True:
		np.savetxt('/home/qyli/oPDFnew/data/Trueini/%s_%s_InL_rcin%s_%s_%s.txt' %(model,pp,rin,ppn,fitmodel),InLvalue)
	# radius
	if rout == 1:
		np.savetxt('/home/qyli/oPDFnew/data/Trueini/%s_%s_fitmc_rcin%s_%s_%s.txt' %(model,pp,rin,ppn,fitmodel),mcfit)
	else:
		np.savetxt('/home/qyli/oPDFnew/data/Trueini/%s_%s_fitmc_rcout%s_%s_%s.txt' %(model,pp,rout,ppn,fitmodel),mcfit)

	print('The tracer distribution is %s' %model)
	print('The tracer type is %s' %pp)
	print('Tracer number for each cluster is %s' %ppn)
	print('The tracer selection region is %s and %s' %(rin,rout))

#=================================================================
'''
model: the region in simulations
pp: tracer type
ppn: tracer number
rin: inner radius cut (two type: 1. direct cut; 2. relative radius)
rout: outer radius cut (if set 1 taking rvir)

dx0: the intergal for likelihood
sigcont: calculate the statistic error
InL: calculate the likelihood difference between best-fit and true

dx0: 0.02(DM); 0.2(star); 0.6(ste);  
'''
main(model = 'GXsub', pp = 'DM', ppn = '100000', fitmodel = 'TMP', rin = 200, rout = 1, dx0 = 0.02, sigcont = False, InL = False)











