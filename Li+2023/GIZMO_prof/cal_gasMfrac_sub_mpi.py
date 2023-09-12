import numpy as np
import h5py
# from readsnapsgl import readsnapsgl
from astropy import constants as const
from astropy import units
from scipy.spatial import cKDTree
import os, sys
import glob
import readsnapsgl
import readsubpp
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

#produce data
Xspath = "/home2/weiguang/data7/Gizmo-Simba/"
def produce_prof(sngrp, red, crn, cn, snapn, prop, hostid):
	#particle type: 0: gas; 1: DM; 4: stars

	filename = Xspath+cn+'/'+snapn+'.hdf5'
	f=h5py.File(filename, 'r')

	sn = sngrp.attrs['snapshot']

	#read snapshot data: particle type: 0: gas; 1: DM; 4: stars
	#--------------------------------------------------------------------------------------
	#--------------------------------------------------------------------------------------
	j = 0
	subhaloid, outsubid, locsubid, locoutsubid = readsubpp.readsubid(sn, crn, pt = j, ptype = 'SUB')
	print('finish reading particles ID')
	pos = f['PartType'+str(j)+'/Coordinates'][:]
	pos = pos[locoutsubid]
	rr = np.sqrt(np.sum((pos-prop['cc'])**2, axis=1))
	mas = f['PartType'+str(j)+'/Masses'][:]
	mas = mas[locoutsubid]
	temp = readsnapsgl.readhdf5data(filename, block = 'temperature', quiet=False, ptype=j)
	temp = temp[locoutsubid]

	#log space bins
	bcgratio = 0.05
	ids0_hot = np.where((temp>=1.0e7)&(rr<1.0*prop['r500'])&(rr>bcgratio*prop['r500']))[0] #select satisfied particles
	ids0_warm = np.where((temp>=1.0e5)&(temp<1e7)&(rr<1.0*prop['r500'])&(rr>bcgratio*prop['r500']))[0] #select satisfied particles
	ids0_warm0 = np.where((temp>=1.0e5)&(temp<1e6)&(rr<1.0*prop['r500'])&(rr>bcgratio*prop['r500']))[0] #select satisfied particles
	ids0_warm1 = np.where((temp>=1.0e6)&(temp<1e7)&(rr<1.0*prop['r500'])&(rr>bcgratio*prop['r500']))[0] #select satisfied particles
	ids0_cold = np.where((temp<1.0e5)&(rr<1.0*prop['r500'])&(rr>bcgratio*prop['r500']))[0] #select satisfied particles

	hot_Mfrac = np.sum(mas[ids0_hot]) * 1e10 / prop['M500']
	warm_Mfrac = np.sum(mas[ids0_warm]) * 1e10 / prop['M500']
	warm_Mfrac0 = np.sum(mas[ids0_warm0]) * 1e10 / prop['M500']
	warm_Mfrac1 = np.sum(mas[ids0_warm1]) * 1e10 / prop['M500']
	cold_Mfrac = np.sum(mas[ids0_cold]) * 1e10 / prop['M500']

	#--------------------------------------------------------------------------------------
	#--------------------------------------------------------------------------------------
	j = 4
	subhaloid, outsubid, locsubid, locoutsubid = readsubpp.readsubid(sn, crn, pt = j, ptype = 'SUB')
	pos0 = f['PartType'+str(j)+'/Coordinates'][:]
	rr0 = np.sqrt(np.sum((pos0-prop['cc'])**2, axis=1))
	mas0 = f['PartType'+str(j)+'/Masses'][:]
	pos = pos0[locoutsubid]
	rr = rr0[locoutsubid]
	mas = mas0[locoutsubid]

	#log space bins
	bcgratio = 0.05
	ids_ICL = np.where((rr<1.0*prop['r500'])&(rr>bcgratio*prop['r500']))[0] #select satisfied particles
	ICL_Mfrac = np.sum(mas[ids_ICL]) * 1e10 / prop['M500']
	ids_tot = np.where((rr0<1.0*prop['r500']))[0] #select satisfied particles
	tot_Mfrac = np.sum(mas0[ids_tot]) * 1e10 / prop['M500']

	#===save data===
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

	if 'ICLfrac' in sngrp.keys():  #Rbins
		sngrp['ICLfrac'][:]=ICL_Mfrac
	else:
		sngrp.create_dataset('ICLfrac', data = ICL_Mfrac)

	if 'Stefrac' in sngrp.keys():  #Rbins
		sngrp['Stefrac'][:]=tot_Mfrac
	else:
		sngrp.create_dataset('Stefrac', data = tot_Mfrac)
	
	f.close()

#sel which snap in which region to compute and prepare basic host halo properties
def sel_sn(progenIDs, sn_sa_red, snapseq):
	'''
	progenIDs: progenitor ID
	sn_sa_red: redshift information
	snapseq: splited snapshot sequence
	'''
	for sn in snapseq: # from snap_128 to snap_39
		snapn='snap_%03d'%sn
		red=sn_sa_red[sn,2]

		h5fname="/home2/qingyang/GIZMO_prof/data/subgasMfrac_allsn/GIZ-profiles-progenitors_subgasMfrac_sn%s.hdf5" %sn
		if not os.path.isfile(h5fname):
			fh5w = h5py.File(h5fname, "w")
			fh5w.close()

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
				fh5w.close()
				print('No progenitors at ',cn, snapn)
				continue

			else:
				tmpd = np.load('/home2/weiguang/Project-300-Clusters/Halo_mass_function_mass-difference/GIZMO/GS_Mass_snap_%03dinfo.npy' %sn)
				#ReginIDs HIDs  HosthaloID Mvir(4) Xc(5)   Yc(6)   Zc(7)  Rvir(8) fMhires(38) cNFW (42) Mgas200 M*200 M500(13)  R500(14) fgas500 f*500
				ids = np.where((np.int32(tmpd[:,0])==crn) & (np.int64(tmpd[:,1]) == progenIDs[crn-1, sn]))[0]
				if len(ids) != 1:
					print("can not find central galaxy for ", crn, snapn, "len(ids):", len(ids))
					fh5w.close()
					continue

			cc = tmpd[ids,4:7][0]
			r200 = tmpd[ids,7][0]
			r500 = tmpd[ids, 13][0]
			M500 = tmpd[ids, 12][0] * 1e10
			M200 = tmpd[ids, 3][0]

			sngrp.attrs['redshift'] = red; sngrp.attrs['snapshot'] = sn;
			sngrp.attrs['center_x'] = cc[0]; sngrp.attrs['center_y'] = cc[1]; sngrp.attrs['center_z'] = cc[2];
			sngrp.attrs['r200'] = r200; sngrp.attrs['M200'] = M200;
			sngrp.attrs['r500'] = r500; sngrp.attrs['M500'] = M500;

			prop = {'cc':cc, 'r200':r200, 'r500':r500, 'M500':M500}
			
			produce_prof(sngrp, red, crn, cn, snapn, prop, progenIDs[crn-1, sn])

			fh5w.close()

#calculate profiles from Gadget-X profiles in nifty2014@castor.ft.uam.es
#consider evolution profiles (2021.4.13)

def main():
	progenIDs=np.loadtxt("/home2/weiguang/Project-300-Clusters/Halo_mass_function_mass-difference/GIZMO/Progenitor-IDs-for-center-cluster.txt",dtype=np.int64)
	sn_sa_red=np.loadtxt('/home2/weiguang/Project-300-Clusters/redshifts.txt')

	seq = np.arange(128,39-1,-1)
	snapseq = np.array_split(seq, size)
	sel_sn(progenIDs, sn_sa_red, snapseq[rank])
main()
