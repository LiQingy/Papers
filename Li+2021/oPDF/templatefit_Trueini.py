import numpy as np 
from oPDF.oPDF import *

def templatefit(fitmodel, model, pp, rin,rout, nubtype, nbin):
	if model == 'GXsub':
		filecen = '/home/qyli/oPDFnew/data/cluster/G3X_Mass_snap_128-center-cluster.txt'
	else:
		filecen = '/home/qyli/oPDFnew/data/cluster/MDPL2_Mass_snap_128-center-cluster.txt'
	datacen = np.loadtxt(filecen)
	if fitmodel == 'TMP':
		filepp = '/home/qyli/oPDFnew/data/Trueini/%s_%s_fitmc_rcin%s_%s_TMP.txt' %(model,pp,rin,nubtype)
	else:
		filepp = '/home/qyli/oPDFnew/data/Trueini/%s_%s_fitmc_rcin%s_%s_NFW.txt' %(model,pp,rin,nubtype)
	datapp = np.loadtxt(filepp)

	template_M = np.zeros(shape = (324,nbin))
	for i in range(324):
		#rv = np.log10(datacen[i][6])
		if fitmodel == 'TMP':
			halo = Halo(halotype=HaloTypes.TMPMC, TMPid = i)
		else:
			halo = Halo(halotype=HaloTypes.NFWMC)
		m = datapp[i][0] * 1e5
		c = datapp[i][1]
		halo.set_param([m,c])
		# rs = (m * 1e10 / ( 4/ 3 * np.pi * 2.7755e2 * 200))**(1/3) / c
		rv = datacen[i][6]
		r = 10**(np.linspace(0, np.log10(rv), nbin+1)[1:])  #postion is r/rv, where the rv is true value
		template_M[i] = halo.mass(r)
	np.savetxt('/home/qyli/oPDFnew/data/Trueini/%s_%s_rin%s_rout%s_%s_%s_M.txt' %(model,pp,rin,rout,nubtype,fitmodel),template_M)

templatefit(fitmodel = 'TMP',model = 'GXsub', pp = 'DM', rin = 0, rout = 1, nubtype = '100000', nbin = 20)
	
