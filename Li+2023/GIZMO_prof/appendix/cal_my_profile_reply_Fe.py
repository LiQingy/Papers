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
	metl= f['PartType'+str(j)+'/Metallicity'][:]

	metFe = metl[:,-1] #Fe
	###particle selection conditions
	ids0 = np.where((temp>1.0e6)&(rr<1.5*prop['r200'])&(delaytime<=0)&(sgden<2.88e6))[0] #select satisfied particles
	rr=rr[ids0]
	mas=mas[ids0]*1.0e10
	metFe=metFe[ids0]

	# #--------------------------------------------------------------------------------------
	# #--------------------------------------------------------------------------------------
	# f.close()
	# #log space bins
	nbin = 50
	rbins=np.logspace(np.log10(0.001*prop['r200']), np.log10(1.5*prop['r200']), num=nbin+1)

	#--------------------------------------------------------------------------------------
	#--------------------------------------------------------------------------------------
	# data=np.zeros((dens.size,4),dtype=np.float64)
	# data[:,0] = (rad[1:]+rad[:-1])/2
	num0,xn0 = np.histogram(rr, bins=rbins)
	met,xe = np.histogram(rr, bins=rbins, weights=mas*metFe)
	tmass,xe = np.histogram(rr, bins=rbins, weights=mas)
	
	#--------------------------------------------------------------------------------------
	#--------------------------------------------------------------------------------------
	#===save data===
	np.seterr(divide='ignore', invalid='ignore')
	if 'Rbin' in sngrp.keys():  #Rbins
		sngrp['Rbin'][:]=(rbins[1:] + rbins[:-1])/2
	else:
		sngrp.create_dataset('Rbin', data = (rbins[1:] + rbins[:-1])/2) 

	#particle number in bins
	#shape:(3,nbin); 1: gas;
	if 'NuminBin' in sngrp.keys():  #Rbins
		sngrp['NuminBin'][:]=num0
	else:
		sngrp.create_dataset('NuminBin', data = num0)

	# Mass-Weighted Fe metallicity
	if 'MWFe' in sngrp.keys():
		sngrp['MWFe'][:]=met/tmass
	else:
		sngrp.create_dataset('MWFe', data = met/tmass) 


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
	h5fname = "/home2/qingyang/GIZMO_prof/data/GIZ-profiles-progenitors_reply_Fe.hdf5"
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

	sel_sn(128, h5fname, progenIDs, sn_sa_red)
	sel_sn(96, h5fname, progenIDs, sn_sa_red)
	sel_sn(78, h5fname, progenIDs, sn_sa_red)
	sel_sn(65, h5fname, progenIDs, sn_sa_red)
	sel_sn(55, h5fname, progenIDs, sn_sa_red)
main()




