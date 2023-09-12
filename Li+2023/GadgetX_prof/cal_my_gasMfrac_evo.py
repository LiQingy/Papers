import numpy as np
import h5py
# from readsnapsgl import readsnapsgl
from astropy import constants as const
from astropy import units
from scipy.spatial import cKDTree
import os, sys
import glob
import readsnapsgl as readsgl



#global units=======================================================
mean_mol_weight = 0.588
prtn = 1.67373522381e-24  # (proton mass in g)
bk = 1.3806488e-16        # (Boltzman constant in CGS)

# Cosmologies critical density in 10^10 M_sun/h/(kpc/h)**3
from astropy.cosmology import FlatLambdaCDM
cosm = FlatLambdaCDM(H0=67.77, Om0=0.307115)
# rhoc500=500*cosm.critical_density0.to("M_sun/kpc**3")/1.0e10/0.6777**2
#=======================================================
Xspath = '/home/nifty2014/TheThreeHundred/simulation/GadgetX/'

'''
head info:
M200: {:.3e} [Msun/h], 
R200 {:.3f} [kpc/h], 
M500: {:.3e} [Msun/h], 
R500 {:.3f} [kpc/h], 
center: {:.3f} {:.3f} {:.3f} [kpc/h] 250 ratial bins from 0.001 - 1.5 R200 \n gas density [Msun/h/(kpc/h)^3]; 
MW temp [k]; 
Pressure []; 
electron number density []; 
entropy []; 
metal [];

format( M200, r200, cc[0], cc[1], cc[2])
'''

#produce data
def produce_prof(sngrp, red, cn, snapn, prop, hostid):

	#read snapshot data 
	filename = Xspath+cn+'/'+snapn

	j = 0 #particle type: 0: gas; 1: DM; 4: stars
	pos = readsgl.readsnapsgl(filename, block = 'POS ',ptype = j)
	rr = np.sqrt(np.sum((pos-prop['cc'])**2, axis=1))
	mas = readsgl.readsnapsgl(filename, block = 'MASS',ptype = j)
	# sgden=readsgl.readsnapsgl(filename, block = 'RHO ',ptype = j)* 1.0e10 /0.6777/(1/0.6777)**3  #in Msun/kpc^3
	# temp= f['PartType'+str(j)+'/InternalEnergy'][:] * (5. / 3 - 1) * v_unit**2 * prtn * mean_mol_weight / bk # in k
	temp = readsgl.readsnapsgl(filename, block = 'TEMP',ptype = j)

	#log space bins
	ids0_hot = np.where((temp>=1.0e7)&(rr<1.0*prop['r500']))[0] #select satisfied particles
	ids0_warm = np.where((temp>=1.0e5)&(temp<1e7)&(rr<1.0*prop['r500']))[0] #select satisfied particles
	ids0_warm0 = np.where((temp>=1.0e5)&(temp<1e6)&(rr<1.0*prop['r500']))[0] #select satisfied particles
	ids0_warm1 = np.where((temp>=1.0e6)&(temp<1e7)&(rr<1.0*prop['r500']))[0] #select satisfied particles
	ids0_cold = np.where((temp<1.0e5)&(rr<1.0*prop['r500']))[0] #select satisfied particles

	hot_Mfrac = np.sum(mas[ids0_hot]) * 1e10 / prop['M500']
	warm_Mfrac = np.sum(mas[ids0_warm]) * 1e10 / prop['M500']
	warm_Mfrac0 = np.sum(mas[ids0_warm0]) * 1e10 / prop['M500']
	warm_Mfrac1 = np.sum(mas[ids0_warm1]) * 1e10 / prop['M500']
	cold_Mfrac = np.sum(mas[ids0_cold]) * 1e10 / prop['M500']
	#--------------------------------------------------------------------------------------

	np.seterr(divide='ignore', invalid='ignore')
	if 'Hotfrac' in sngrp.keys():  #Rbins
		sngrp['Hotfrac'][:]=hot_Mfrac
	else:
		sngrp.create_dataset('Hotfrac', data = hot_Mfrac) 

	if 'Warmfrac' in sngrp.keys():  #Rbins
		sngrp['Warmfrac'][:]=warm_Mfrac
	else:
		sngrp.create_dataset('Warmfrac', data = warm_Mfrac) 

	if 'Warmfrac_cold' in sngrp.keys():  #Rbins
		sngrp['Warmfrac_cold'][:]=warm_Mfrac0
	else:
		sngrp.create_dataset('Warmfrac_cold', data = warm_Mfrac0) 

	if 'Warmfrac_hot' in sngrp.keys():  #Rbins
		sngrp['Warmfrac_hot'][:]=warm_Mfrac1
	else:
		sngrp.create_dataset('Warmfrac_hot', data = warm_Mfrac1) 

	if 'Coldfrac' in sngrp.keys():  #Rbins
		sngrp['Coldfrac'][:]=cold_Mfrac
	else:
		sngrp.create_dataset('Coldfrac', data = cold_Mfrac) 


#sel which snap in which region to compute and prepare basic host halo properties
def sel_sn(h5fname, progenIDs, sn_sa_red):

	for sn in range(128,24,-1): # from snap_128 to snap_25
		snapn='snap_%03d'%sn
		red=sn_sa_red[sn,2]

		for crn in np.arange(1,325):
			cn = 'NewMDCLUSTER_%04d'%crn

			#h5 write files
			fh5w = h5py.File(h5fname, "r+")
			if cn in fh5w.keys(): #which region: 1-324
				cgrp = fh5w[cn]
			else:
				cgrp = fh5w.create_group(cn)

			if snapn in cgrp.keys(): #which snapshot
				sngrp=cgrp[snapn]
			else:
				sngrp = cgrp.create_group(snapn)

			#Determine the basic property of host halo at z = 0
			if (progenIDs[crn-1, sn]<1):
				print('No progenitors at ',cn, snapn)
				continue
			else:
				tmpd = np.load('/home/nifty2014/TheThreeHundred/playground/qingyang/GadgetX_prof/data/GadgetX/G3X_Mass_snap_%03dinfo.npy' %sn)
				#ReginIDs HIDs  HosthaloID Mvir(4) Xc(5)   Yc(6)   Zc(7)  Rvir(8) fMhires(38) cNFW (42) Mgas200 M*200 M500(13)  R500(14) fgas500 f*500
				ids = np.where((np.int32(tmpd[:,0])==crn) & (np.int64(tmpd[:,1]) == progenIDs[crn-1, sn]))[0]
				if len(ids) != 1:
					print("can not find central galaxy for ", crn, snapn, "len(ids):", len(ids))
					continue

			cc = tmpd[ids,4:7][0]
			r200 = tmpd[ids,7][0]
			r500 = tmpd[ids, 13][0]
			M500 = tmpd[ids, 12][0] * 1e10
			M200 = tmpd[ids, 3][0]
			print(M200,r200,cc)

			sngrp.attrs['redshift'] = red; sngrp.attrs['snapshot'] = sn;
			sngrp.attrs['center_x'] = cc[0]; sngrp.attrs['center_y'] = cc[1]; sngrp.attrs['center_z'] = cc[2];
			sngrp.attrs['r200'] = r200; sngrp.attrs['M200'] = M200;
			sngrp.attrs['r500'] = r500; sngrp.attrs['M500'] = M500;


			prop = {'cc':cc, 'r200':r200, 'r500':r500, 'M500':M500}
			
			produce_prof(sngrp, red, cn, snapn, prop, progenIDs[crn-1, sn])

			fh5w.close()

			print('redshift: ', sn, 'region: ', crn)

#calculate profiles from Gadget-X profiles in nifty2014@castor.ft.uam.es
#consider evolution profiles (2021.4.13)
def main():

	progenIDs=np.loadtxt('/home/nifty2014/TheThreeHundred/playground/qingyang/GadgetX_prof/data/GadgetX/Progenitor-IDs-for-center-cluster.txt',dtype=np.int64)
	sn_sa_red=np.loadtxt('/home/nifty2014/TheThreeHundred/playground/weiguang/redshifts.txt')

	##Open h5 file for write
	h5fname="/home/nifty2014/TheThreeHundred/playground/qingyang/GadgetX_prof/data/G3X-profiles-progenitors_myevoz_gasMfrac_plus.hdf5"
	if not os.path.isfile(h5fname):
		fh5w = h5py.File(h5fname, "w")
		fh5w.close()

	sel_sn(h5fname, progenIDs, sn_sa_red)
main()