""" This file creates the potential and cumulative density profile templates for a given DM file"""
from oPDF.oPDF import *
from oPDF.myutils import *
import h5py,os,sys
from scipy.stats import chi2
#plt.ion()
def calPotential(model,ptype,filevx,rv):
	#set units if needed
	hubble=0.678
	Globals.set_units(1e10, 1, 1) #set to 1e10Msun/h, kpc/h, km/s with the current h

	DMfile = filevx
	DMfile = DMfile.encode('utf8')
	#real parameters, for comparison with analytical profile:

	nbin=100 #do not change this, unless you change TemplateData.h as well.
	npart=0 #number of particles to use. 0 means FullSample.

	FullSample=Tracer(DMfile)
	Sample=FullSample.copy(0,npart)
	if model == 'MDPL2':
		Sample_M = np.tile(0.150541, h5py.File(filevx,'r')['x'][:].shape[0])
	elif model == 'GX' and ptype == 'DM':
		Sample_M = np.tile(0.127, h5py.File(filevx,'r')['x'][:].shape[0])
	elif model == 'GX' and ptype == 'gas':
		Sample_M = np.tile(0.0236, h5py.File(filevx,'r')['x'][:].shape[0])
	elif model == 'GX' and ptype == 'star':
		Sample_M = np.array(h5py.File(filevx,'r')['mass'][:])
	print(Sample_M[0])
	FullSample.clean()
	
	xbin=np.logspace(np.log10(0.1), np.log10(rv), nbin)

	countM,tmp=np.histogram(Sample.data['r'], np.hstack([0., xbin]), weights = Sample_M)#dM
	countR,tmp=np.histogram(Sample.data['r'], np.hstack([0., xbin]), weights=Sample_M/Sample.data['r'])#dM/R
	pot=countM.cumsum()/xbin+countR.sum()-countR.cumsum() #M(<r)/r+\int_r^rmax dM/r
	pot*=Globals.units.Const.G
	density_cum=countM.cumsum()/xbin**3/(4*np.pi/3) 
	#pad with bin 0
	xbin=np.hstack([0., xbin])
	pot=np.hstack([countR.sum()*Globals.units.Const.G, pot])
	density_cum=np.hstack([density_cum[0], density_cum])
	return xbin,pot,density_cum

def main(model):
		if model == 'MDPL2':
			filecen = '/home/qyli/oPDFnew/data/cluster/MDPL2_Mass_snap_128-center-cluster.txt'
			npar = 1
		elif model == 'GX':
			filecen = '/home/qyli/oPDFnew/data/cluster/G3X_Mass_snap_128-center-cluster.txt'
			npar = 3
		datacen = np.loadtxt(filecen)
		pp = np.array(['DM','star','gas'])

		for i in range(324):
			rv = datacen[i][6]
			template = np.zeros(shape = (3,101))
			for pn in range(npar):
				if model == 'MDPL2':
					filevx = '/home/qyli/oPDFnew/data/inidata/%s/%s_%s_vx%s.h5'  %(model, model, pp[pn], i+1)
				else:
					filevx = '/home/qyli/oPDFnew/data/inidata/%s_%s/%s_%s_vx%s.h5'  %(model, pp[pn], model, pp[pn], i+1)
				rbin, pot, dens_cul = calPotential(model,pp[pn],filevx,rv)
				template[0] = rbin
				template[1] += pot
				template[2] += dens_cul
			np.savetxt('/home/qyli/oPDFnew/data/template/%s/rin0_rout1_rpd%s.txt' %(model,i+1), template)
#calculate radius, potential, culmulative density
main(model = 'GX')


