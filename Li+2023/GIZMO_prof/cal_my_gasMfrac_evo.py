import numpy as np
import h5py
# from readsnapsgl import readsnapsgl
from astropy import constants as const
from astropy import units
from scipy.spatial import cKDTree
import os, sys
import glob
import readsnapsgl



#global units=======================================================
mean_mol_weight = 0.588
prtn = 1.67373522381e-24  # (proton mass in g)
bk = 1.3806488e-16        # (Boltzman constant in CGS)

from astropy.cosmology import FlatLambdaCDM
cosm = FlatLambdaCDM(H0=67.77, Om0=0.307115)

#=======================================================
Xspath = "/home2/weiguang/data7/Gizmo-Simba/"

#produce data
def produce_prof(sngrp, red, cn, snapn, prop, hostid):

	#load data 
	filename = Xspath+cn+'/'+snapn+'.hdf5'
	f=h5py.File(filename, 'r')

	#read snapshot data: particle type: 0: gas; 1: DM; 4: stars
	#--------------------------------------------------------------------------------------
	#--------------------------------------------------------------------------------------
	j = 0
	pos = f['PartType'+str(j)+'/Coordinates'][:]
	rr = np.sqrt(np.sum((pos-prop['cc'])**2, axis=1))
	mas = f['PartType'+str(j)+'/Masses'][:]
	sgden=f['PartType'+str(j)+'/Density'][:] * 1.0e10 /0.6777/(1/0.6777)**3  #in Msun/kpc^3
	# temp= f['PartType'+str(j)+'/InternalEnergy'][:] * (5. / 3 - 1) * v_unit**2 * prtn * mean_mol_weight / bk # in k
	temp = readsnapsgl.readhdf5data(filename, block = 'temperature', quiet=False, ptype=j)
	delaytime = f['PartType'+str(j)+'/DelayTime'][:]

	# #--------------------------------------------------------------------------------------
	# #--------------------------------------------------------------------------------------

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
	if 'Hotfrac' in sngrp.keys():
		sngrp['Hotfrac'][:]=hot_Mfrac
	else:
		sngrp.create_dataset('Hotfrac', data = hot_Mfrac) 

	if 'Warmfrac' in sngrp.keys():
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

	if 'Coldfrac' in sngrp.keys(): 
		sngrp['Coldfrac'][:]=cold_Mfrac
	else:
		sngrp.create_dataset('Coldfrac', data = cold_Mfrac) 

	f.close()

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

			#Now we use the cluster withID
			if (progenIDs[crn-1, sn]<1):
				print('No progenitors at ',cn, snapn)
				continue
			else:
				tmpd=np.load('/home2/weiguang/Project-300-Clusters/Halo_mass_function_mass-difference/GIZMO/GS_Mass_snap_%03dinfo.npy' %sn)
				#ReginIDs HIDs  HosthaloID Mvir(4) Xc(5)   Yc(6)   Zc(7)  Rvir(8) fMhires(9) cNFW (10) Mgas200 M*200 M500(13)  R500(14) fgas500 f*500
				ids = np.where((np.int32(tmpd[:,0])==crn) & (np.int64(tmpd[:,1]) == progenIDs[crn-1, sn]))[0]
				if len(ids) != 1:
					print("can not find central galaxy for ", crn, snapn, "len(ids):", len(ids))
					continue

			cc = tmpd[ids, 4:7][0]
			r200 = tmpd[ids, 7][0]
			r500 = tmpd[ids, 13][0]
			M500 = tmpd[ids, 12][0]
			M200 = tmpd[ids,3][0]
			print(M200, r200, cc)

			sngrp.attrs['redshift'] = red; sngrp.attrs['snapshot'] = sn;
			sngrp.attrs['center_x'] = cc[0]; sngrp.attrs['center_y'] = cc[1]; sngrp.attrs['center_z'] = cc[2];
			sngrp.attrs['r200'] = r200; sngrp.attrs['r500'] = r500;
			sngrp.attrs['M200'] = M200; sngrp.attrs['M500'] = M500;

			prop = { 'cc':cc,'r200':r200, 'r500':r500, 'M500': M500}

			produce_prof(sngrp, red, cn, snapn, prop, progenIDs[crn-1, sn])

			fh5w.close()

			print('redshift: ', sn, 'region: ', crn)


#calculate profiles from GIZMO profiles in qingyang@brutus.ft.uam.es
def main():

	progenIDs=np.loadtxt("/home2/weiguang/Project-300-Clusters/Halo_mass_function_mass-difference/GIZMO/Progenitor-IDs-for-center-cluster.txt",dtype=np.int64)
	sn_sa_red=np.loadtxt('/home2/weiguang/Project-300-Clusters/redshifts.txt')
	with open('/home2/qingyang/GIZMO_prof/readme.txt') as freadme:
  		readme = freadme.read()

	##Open h5 file for write
	h5fname="/home2/qingyang/GIZMO_prof/data/GIZ-profiles-progenitors_myevoz_gasMfrac_plus.hdf5"
	if not os.path.isfile(h5fname):
		fh5w = h5py.File(h5fname, "w")
		fh5w.close()

	sel_sn(h5fname, progenIDs, sn_sa_red)

main()




