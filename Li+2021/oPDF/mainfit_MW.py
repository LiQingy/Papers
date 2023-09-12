import numpy as np 
import h5py
import os
from oPDF.oPDF import *

def oPDFfit(filevx,rin_cut,rout_cut,ppn,i, Mtrue, Cture, hub):
	filevx = filevx.encode('utf8')
	Globals.set_units(1e10*1*(0.704/hub), 1*(0.704/hub), 1) #set to 1e10Msun, kpc, km/s
	Globals.get_units()
	

	Sample = Tracer(filevx, rmin = rin_cut, rmax = rout_cut)
	Sample = Sample.select(((Sample.data["subid"]==0) + (Sample.data["subid"]>1000000)))
	Sample = Sample.copy(0,ppn)

	Sample.halo.set_type(HaloTypes.TMPMC,scales = [1, 1], TMPid =i)
	Estimators.RBinLike.nbin=10
	mc,fval,status = Sample.dyn_fit(Estimators.RBinLike, x0 = [Mtrue, Cture])
	return mc

def main(model,pptype,ppn,hub):
	if hub == 0.704:
		fdir = '/home/qyli/oPDFnew/data/apostle_MR_h/'
	else:
		fdir = '/home/qyli/oPDFnew/data/apostle_MR/'

	if model == 'S':
		rvfof = np.array([224.1709,172.7433,236.8058,206.0765,192.1194,173.7896,233.5474,194.9914,199.178,189.7775,197.9979,170.9422])
		Ctruefof = np.array([5.8976, 8.8382, 8.8435, 6.9903, 12.776, 8.7318, 10.0735, 7.18, 9.8334, 8.2011, 12.346, 11.1352])
		Mtruefof = np.array([129.82, 59.4, 153.03, 100.85, 81.72, 60.49, 146.8, 85.44, 91.06, 78.77, 89.45, 57.56])
	else:
		rvfof = np.array([242.2135,206.9113,187.9074,186.0567,235.1194,218.5439,221.8365,221.4039,197.6648,193.7516,265.3684,214.9532])
		Ctruefof = np.array([12.4923,4.1164,11.2907,9.8764,8.7675,4.2887,13.519,3.701,8.1664,11.4837,8.0264,8.6198])
		Mtruefof = np.array([163.76, 102.08, 76.46, 74.22, 149.79, 120.29, 125.81, 125.07, 89, 83.82, 215.35, 114.46])
	
	
	mcfit = np.zeros(shape = (12,2))
	for i in range(11,12):
		rin_cut = 20 * hub
		rout_cut = rvfof[i] * hub

		Mtrue = Mtruefof[i] * hub
		Cture = Ctruefof[i]
		filevx = fdir + '%s%sfof%spart%s.hdf5' %(model,int(i/2)+1,i%2+1,pptype)
		
		mc = oPDFfit(filevx,rin_cut,rout_cut,ppn,i,Mtrue,Cture,hub)
		mcfit[i] = mc
		print(i,mc)
	#np.savetxt('/home/qyli/oPDFnew/data/MW/%sfof_mcfit_%s.txt' %(model,pptype), mcfit)

main(model = 'S', pptype = '4', ppn = 0, hub = 1.)
