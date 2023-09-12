""" This file creates the potential and cumulative density profile templates for a given DM file"""
from oPDF.oPDF import *
from oPDF.myutils import *
import h5py,os,sys
from scipy.stats import chi2
#plt.ion()
def calPotential(model,ptype,filevx,rv,massDM,hub):
	#set units if needed
	Globals.set_units(1e10*1 * (0.704 / hub), 1 * (0.704/hub), 1) #set to 1e10Msun, kpc, km/s
	Globals.get_units()

	DMfile = filevx
	DMfile = DMfile.encode('utf8')
	#real parameters, for comparison with analytical profile:

	nbin=100 #do not change this, unless you change TemplateData.h as well.
	npart=0 #number of particles to use. 0 means FullSample.

	FullSample=Tracer(DMfile)
	Sample=FullSample.copy(0,npart)
	# if model == 'S':
	# 	if ptype == 'DM':
	# 		Sample_M = np.tile(massDM, h5py.File(filevx,'r')['x'][:].shape[0])
	# 	elif ptype == 'gas':
	# 		Sample_M = np.array(h5py.File(filevx,'r')['PartMass'][:]) * hub
	# 	elif ptype == 'star':
	# 		Sample_M = np.array(h5py.File(filevx,'r')['PartMass'][:]) * hub
	# 	else:
	# 		Sample_M = np.array(h5py.File(filevx,'r')['PartMass'][:]) * hub
	# else:
	# 	if ptype == 'DM':
	# 		Sample_M = np.tile(massDM, h5py.File(filevx,'r')['x'][:].shape[0])
	# 	elif ptype == 'gas':
	# 		Sample_M = np.array(h5py.File(filevx,'r')['PartMass'][:]) * hub
	# 	elif ptype == 'star':
	# 		Sample_M = np.array(h5py.File(filevx,'r')['PartMass'][:]) * hub
	# 	else:
	# 		Sample_M = np.array(h5py.File(filevx,'r')['PartMass'][:]) * hub
	
	FullSample.clean()

	
	xbin=np.logspace(np.log10(0.5 * hub), np.log10(rv), nbin)
	if ptype == 'DM':
		Sample.mP = massDM * 1.
		Sample_M = Sample.mP * Sample.data['w']
		countM,tmp=np.histogram(Sample.data['r'] * hub, np.hstack([0.,xbin]))#dM
		countM=countM*Sample.mP
		countR,tmp=np.histogram(Sample.data['r'] * hub, np.hstack([0.,xbin]), weights=Sample.mP/Sample.data['r']/hub)#dM/R
	else:
		Sample_M = Sample.mP * Sample.data['w'] * hub
		countM,tmp=np.histogram(Sample.data['r'] * hub, np.hstack([0., xbin]), weights = Sample_M)#dM
		countR,tmp=np.histogram(Sample.data['r'] * hub, np.hstack([0., xbin]), weights=Sample_M/Sample.data['r']/hub)#dM/R
	pot=countM.cumsum()/xbin+countR.sum()-countR.cumsum() #M(<r)/r+\int_r^rmax dM/r
	# G = 4.3007*10.**4
	# pot*=G
	pot*=Globals.units.Const.G
	density_cum=countM.cumsum()/xbin**3/(4*np.pi/3) 
	#pad with bin 0
	xbin=np.hstack([0., xbin])
	# pot=np.hstack([countR.sum()*G, pot])
	pot=np.hstack([countR.sum()*Globals.units.Const.G, pot])
	density_cum=np.hstack([density_cum[0], density_cum])
	Sample.clean()
	return xbin,pot,density_cum

def main(model,hub):
		if model == 'S':
			masstable=[3.868920866941207E-5,3.868920866941207E-5,3.960878906555616E-5,3.960878906555616E-5,3.8168845221394696E-5,3.8168845221394696E-5,3.8646087907235435E-5,3.8646087907235435E-5,3.876494719513764E-5,3.876494719513764E-5,3.8168845221394696E-5,3.8168845221394696E-5]
			rvfof = np.array([224.1709,172.7433,236.8058,206.0765,192.1194,173.7896,233.5474,194.9914,199.178,189.7775,197.9979,170.9422]) 
		else:
			masstable=[4.216002105295047E-5,4.216002105295047E-5,4.362418207941065E-5,4.362418207941065E-5,4.374081474372147E-5,4.374081474372147E-5,4.284229790380992E-5,4.284229790380992E-5,4.362416693231139E-5,4.362416693231139E-5,4.448360728554904E-5,4.448360728554904E-5]
			rvfof = np.array([242.2135,206.9113,187.9074,186.0567,235.1194,218.5439,221.8365,221.4039,197.6648,193.7516,265.3684,214.9532]) 
		pp = np.array(['DM','star','gas','BH'])
		pnub = ['1','4long','0','5']
		fdir = '/home/qyli/oPDFnew/data/apostle_MR/'
		
		for i in range(12):
			rv = 1.5 * rvfof[i] * hub
			template = np.zeros(shape = (3,101))
			massDM = masstable[i] / (0.704 / hub) * 10

			for pn in range(4):
				filevx = fdir + '%s%sfof%spart%s.hdf5' %(model,int(i/2)+1,i%2+1,pnub[pn])
				rbin, pot, dens_cul = calPotential(model,pp[pn],filevx,rv,massDM,hub)

				template[0] = rbin
				template[1] += pot
				template[2] += dens_cul
			np.savetxt('/home/qyli/oPDFnew/data/template/MW_%s/rin0_rout1_rpd%s%s.txt' %(model,int(i/2)+1,i%2+1), template)

#calculate radius, potential, culmulative density
#hub 1: not include h; hub 0.704: include h
main(model = 'S',hub = 1.)
main(model = 'V',hub = 1.)


