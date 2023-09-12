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
	#--------------------------------------------------------------------------------------
	#--------------------------------------------------------------------------------------

	fahf = glob.glob('/data4/niftydata/TheThreeHundred/data/catalogues/AHF/GIZMO/%s/GIZMO-%s.%s.z*.AHF_halos' %(cn,cn,snapn))
	ahf = np.loadtxt(fahf[0])
	idv = np.where(ahf[:,1] == hostid)[0]
	subhalo_x = ahf[idv,5:8]
	subhalo_stellarmass = ahf[idv,64]
	rrsubhalo = np.sqrt(np.sum((subhalo_x-prop['cc'])**2, axis=1))
	rrsubhalo2D = np.sqrt(np.sum((subhalo_x[:,:2]-prop['cc'][:2])**2, axis=1)) #kpc/h
	#--------------------------------------------------------------------------------------
	#--------------------------------------------------------------------------------------
	#===save data===
	np.seterr(divide='ignore', invalid='ignore')
	if 'Subhalo_rr2D' in sngrp.keys():  #subhalo radius
		sngrp['Subhalo_rr2D'][:]=rrsubhalo2D
	else:
		sngrp.create_dataset('Subhalo_rr2D', data = rrsubhalo2D) 

	if 'Subhalo_rr3D' in sngrp.keys():  #subhalo radius
		sngrp['Subhalo_rr3D'][:]=rrsubhalo
	else:
		sngrp.create_dataset('Subhalo_rr3D', data = rrsubhalo) 

	if 'Subhalo_stellarmass' in sngrp.keys():  #subhalo radius
		sngrp['Subhalo_stellarmass'][:]=subhalo_stellarmass
	else:
		sngrp.create_dataset('Subhalo_stellarmass', data = subhalo_stellarmass) 


#sel which snap in which region to compute and prepare basic host halo properties
def sel_sn(sn, h5fname, progenIDs, sn_sa_red):
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

		prop = { 'cc':cc,'r200':r200}

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
	h5fname="/home2/qingyang/GIZMO_prof/data/GIZ-profiles-progenitors_galrr.hdf5"
	if not os.path.isfile(h5fname):
		fh5w = h5py.File(h5fname, "w")
		fh5w.close()

	fh5w = h5py.File(h5fname, "r+")
	cn = 'README'
	if cn in fh5w.keys(): #which region: 1-324
		head = fh5w[cn]
	else:
		head = fh5w.create_group(cn)
	head.attrs['README'] = readme
	fh5w.close()

	sel_sn(128,h5fname, progenIDs, sn_sa_red)
	sel_sn(96,h5fname, progenIDs, sn_sa_red)
	sel_sn(78,h5fname, progenIDs, sn_sa_red)
	sel_sn(65,h5fname, progenIDs, sn_sa_red)
	sel_sn(55,h5fname, progenIDs, sn_sa_red)
	sel_sn(47,h5fname, progenIDs, sn_sa_red)
	sel_sn(40,h5fname, progenIDs, sn_sa_red)

main()




