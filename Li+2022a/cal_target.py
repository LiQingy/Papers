import numpy as np 
import h5py
import os
import sys
sys.path.append('/home/qyli/oPDFnew')
from oPDF.oPDF import *

def oPDFfit(filevx,pp,ppn,fitmodel,rin_cut,rout_cut,i,MCtrue):
	Globals.cosmology.OmegaM0=0.307
	Globals.cosmology.OmegaL0=0.693
	filevx = filevx.encode('utf8')

	Sample = Tracer(filevx, rmin = rin_cut, rmax = rout_cut)
	Sample = Sample.copy(0,ppn)

	Sample.halo.set_type(HaloTypes.TMPMC,scales = [1e5, 1], TMPid = i)
	
	Estimators.RBinLike.nbin=10
	mc,fval,status = Sample.dyn_fit(Estimators.RBinLike, x0 = MCtrue)
	
	f0 = Sample.likelihood(MCtrue, Estimators.RBinLike)
	AD = Sample.likelihood(mc, Estimators.AD)
	MP = Sample.likelihood(mc, Estimators.MeanPhase)
	MPR = Sample.likelihood(mc, Estimators.MeanPhase)
	lnLike=(fval-f0)

	return mc,lnLike,AD,MP,MPR


def main(model,pp,ppn,fitmodel,rin,rout):

	fileCM = '/home/qyli/oPDF/data/cluster/G3X-CM-masscomplete-ALL.txt'
	filecen = '/home/qyli/oPDF/data/cluster/G3X_Mass_snap_128-center-cluster.txt'
	filengal = '/home/qyli/oPDF/data/cluster/GX_ste_truen.txt'

	dataCM = np.loadtxt(fileCM)
	datacen = np.loadtxt(filecen)
	datangal = np.loadtxt(filengal)
	
	ppntemp = 0

	mcfit = np.zeros(shape = (324,2))
	lnLvalue = np.zeros(324)
	ADvalue = np.zeros(324)
	MPvalue = np.zeros(324)
	MPRvalue = np.zeros(324)

	for i in range(324):
		#read snapshots data

		filevx = '/home/qyli/oPDFnew/data/inidata/%s_%s/%s_%s_vx%s.h5' %(model,pp,model,pp,i+1)
		rin_cut = rin
		rout_cut = dataCM[i][3]

		MCtrue = [dataCM[i, 2] / 1e15, dataCM[i, 6]]
		
		if ppn == 'nsame':
			ppntemp = datangal[i]		
		else:
			ppntemp = int(ppn * 1)
	
		mc,lnLike,AD,MP,MPR = oPDFfit(filevx,pp,ppntemp,fitmodel,rin_cut,rout_cut,i,MCtrue)

		mcfit[i] = mc
		print(mc)

		lnLvalue[i] = lnLike
		ADvalue[i] = AD
		MPvalue[i] = MP
		MPRvalue[i] = MPR

		print(lnLike,AD,MP,MPR)

	np.savetxt('/home/qyli/oPDF2/data/target/%s_%s_fitmc_rcin%s_%s_TMP.txt' %(model,pp,rin,ppn),mcfit)
	np.savetxt('/home/qyli/oPDF2/data/target/%s_%s_lnL_rcin%s_%s_TMP.txt' %(model,pp,rin,ppn),lnLvalue)
	# np.savetxt('/home/qyli/oPDF2/data/target/%s_%s_AD_rcin%s_%s_TMP.txt' %(model,pp,rin,ppn),ADvalue)
	# np.savetxt('/home/qyli/oPDF2/data/target/%s_%s_MP_rcin%s_%s_TMP.txt' %(model,pp,rin,ppn),MPvalue)
	# np.savetxt('/home/qyli/oPDF2/data/target/%s_%s_MPR_rcin%s_%s_TMP.txt' %(model,pp,rin,ppn),MPRvalue)

	print('The tracer distribution is %s' %model)
	print('The tracer type is %s' %pp)
	print('Tracer number for each cluster is %s' %ppn)
	print('The tracer selection region is %s and %s' %(rin,rout))

#------------------------------------------------------------------
'''
model: the region in simulations
pp: tracer type
ppn: tracer number
rin: inner radius cut (two type: 1. direct cut; 2. relative radius)
rout: outer radius cut (if set 1 taking rvir)

dx: the intergal for likelihood
sigcont: calculate the statistic error
'''

main(model = 'GX', pp = 'ste', ppn = 'nsame', fitmodel = 'TMP', rin = 200, rout = 1)











