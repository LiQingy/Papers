import numpy as np 
import h5py
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

dgroup = np.loadtxt('/home/qyli/clu_finder/CLAUDS_HSC/catalogue/odata3/CLAUDS_HSC_iband_group')
d2 = np.loadtxt('/home/qyli/clu_finder/CLAUDS_HSC/catalogue/odata3/iCLAUDS_HSC_iband_2')
digal = np.loadtxt('/home/qyli/clu_finder/CLAUDS_HSC/catalogue/odata3/CLAUDS_HSC_iband_igal')
digal0 = np.loadtxt('/home/qyli/clu_finder/CLAUDS_HSC/catalogue/odata3/CLAUDS_HSC_iband_igal0')

#compute the median value and percentile error [16,84]
def cal_err(lumall): #percentile error
	nbin = lumall.shape[1]
	med = np.zeros(nbin)
	errup = np.zeros(nbin)
	errdown = np.zeros(nbin)
	for i in range(nbin):
		selpos = np.where(lumall[:,i] >= 0)[0]
		med[i] = np.median(lumall[selpos,i])
		if len(selpos):
			err = np.percentile(lumall[selpos,i], [16,84])
			errup[i] = err[1] - med[i]
			errdown[i] = med[i] - err[0]
	return med,errup,errdown


#calculate the physical distance
def cal_dr(dz,sel_clu):
	from astropy.coordinates import SkyCoord
	from astropy.cosmology import FlatLambdaCDM
	from astropy import units as u
	cosmo = FlatLambdaCDM(H0=67.4, Om0=0.315)

	clu_ra = dgroup[sel_clu,2]
	clu_dec = dgroup[sel_clu,3]
	clu_redz = dgroup[sel_clu,4]
	d_A = cosmo.angular_diameter_distance(z=clu_redz)

	dtheta = (11*u.Mpc / d_A).to(u.degree, u.dimensionless_angles()).value
	dredz = dz

	idx_region = np.where((np.abs(clu_ra - digal[:,1]) < dtheta) & 
			((clu_dec - digal[:,2]) <  dtheta)  & 
			(np.abs(clu_redz - digal[:,3]) < dredz))[0] 
	gal_lum = digal0[idx_region, -1]
	gal_ra = digal[idx_region, 1]
	gal_dec = digal[idx_region, 2]

	cgal = SkyCoord(gal_ra, gal_dec, unit="deg")
	cclu = SkyCoord(clu_ra, clu_dec, unit="deg")
	sep = cgal.separation(cclu)
	d_r = (sep * d_A).to(u.Mpc, u.dimensionless_angles())

	return gal_lum,d_r

#calculate the total luminosity profile
def cal_lumall(xbin, nbin, zbin, dz, nm):
	sel_redz = np.where((dgroup[:,4] > zbin[0]) & (dgroup[:,4] <= zbin[1]) & (dgroup[:,1] > nm))[0]
	Ngroup = sel_redz.shape[0]
	print('redshift: %s < z < %s, Ngroup is %s at Ngal > %s' %(zbin[0], zbin[1], Ngroup, nm))

	lumall = np.zeros((Ngroup,nbin))
	numall = np.zeros((Ngroup,nbin))
	for i in range(Ngroup):
		gal_lum, d_r = cal_dr(dz, sel_redz[i])
		rr = d_r.value #Mpc
		sumlum,xxb = np.histogram(rr, bins=xbin, weights = gal_lum)
		ngal,xxb = np.histogram(rr, bins=xbin)
		#print(ngal)
		lumall[i] = sumlum
		numall[i] = ngal

	return lumall,numall

def main():

	nbin = 15
	xbin = np.linspace(0, 10, nbin+1) #from 0 Mpc to 10 Mpc
	rbin = (xbin[1:] + xbin[:-1])/2
	nz = 11 #bin number of redshift 
	dz0 = 0.2 #redshift interval

	fh5 =  h5py.File('./edges.hdf5','w')
	fh5['rbin'] = rbin

	nmall = [100,100,50,20,15,10,5,3,2,2,2] #limit for the number member galaxies
	for i in range(nz):
		zbin = np.array([i*0.4,(i+1)*0.4]).round(1)
		lumall,numall = cal_lumall(xbin, nbin, zbin, dz = dz0, nm = nmall[i])
		med,errup,errdown = cal_err(lumall)

		zspace = fh5.create_group('z%s-%s' %(zbin[0],zbin[1]))
		zspace['median'] = med
		zspace['percerr_16'] = errdown
		zspace['percerr_84'] = errup

		print('Finish redshift range (%s, %s]' %(zbin[0],zbin[1]))

	fh5.close()

#calculat the total luminosity profile within 11 Mpc physcial distance 
main()













