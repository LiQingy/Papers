'''
script to calculate the data which needs to reply referee
'''

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
def produce_prof(sngrp, cn, snapn, prop):

	#load data 
	filename = Xspath+cn+'/'+snapn+'.hdf5'
	f=h5py.File(filename, 'r')

	#set unit
	v_unit = 1.0e5 * np.sqrt(f['Header'].attrs['Time'])  # (e.g. 1.0 km/sec)
	m_ptoMsun = const.m_p.to('Msun').value #proton mass in unit Msun
	KtokeV = const.k_B.value / const.e.value * 1e-3 # [K] to [keV] == k_B*T
	kinunit = 1/ units.keV.to('Msun km2 s-2') #from Msun*(km/s)^2 to keV

	#read snapshot data: particle type: 0: gas; 1: DM; 4: stars
	#--------------------------------------------------------------------------------------
	#--------------------------------------------------------------------------------------
	j = 0
	vel = f['PartType'+str(j)+'/Velocities'][:]
	pos = f['PartType'+str(j)+'/Coordinates'][:]
	rr = np.sqrt(np.sum((pos-prop['cc'])**2, axis=1))
	mas = f['PartType'+str(j)+'/Masses'][:]
	sgden=f['PartType'+str(j)+'/Density'][:] * 1.0e10 /0.6777/(1/0.6777)**3  #in Msun/kpc^3
	# temp= f['PartType'+str(j)+'/InternalEnergy'][:] * (5. / 3 - 1) * v_unit**2 * prtn * mean_mol_weight / bk # in k
	temp = readsnapsgl.readhdf5data(filename, block = 'temperature', quiet=False, ptype=j)
	eleab= f['PartType'+str(j)+'/ElectronAbundance'][:]
	delaytime = f['PartType'+str(j)+'/DelayTime'][:]
	sfr = f['PartType'+str(j)+'/StarFormationRate'][:]
	metl= f['PartType'+str(j)+'/Metallicity'][:]
	pot = f['PartType'+str(j)+'/Potential'][:]

	metHe = metl[:,1] #He 
	metM = metl[:,0] #metallicity
	nH = (1 - metHe - metM) * mas * 1.0e10 / 0.6777  #proton number 
	nele = nH * eleab #electron number: proton number multiply relative value

	###particle selection conditions
	ids0 = np.where((temp>1.0e6)&(rr<2*prop['r500'])&(delaytime<=0)&(sgden<2.88e6))[0] #select satisfied particles
	
	rr=rr[ids0]
	mas=mas[ids0]*1.0e10
	temp=temp[ids0]
	nele=nele[ids0]
	metM=metM[ids0]
	pot = pot[ids0]
	vel = vel[ids0]
	hubblev = pos[ids0] / 0.6777 * 67.77 / 1e3 # Hubble flow vH = xH0; unit: km/s
	vv = vel + hubblev

	# kin_energy_gas=np.sum(vel**2,axis=1,dtype=np.float64)*mas/2./0.6777 #unit: Msun (km/s)^2
	# them_energy_gas=3.*KtokeV*mas / 0.6777 * temp / 2. #unit: keV
	# #--------------------------------------------------------------------------------------
	# #--------------------------------------------------------------------------------------
	#select star mass
	j = 4
	vel4 = f['PartType'+str(j)+'/Velocities'][:]
	pos4 = f['PartType'+str(j)+'/Coordinates'][:]
	pot4 = f['PartType'+str(j)+'/Potential'][:]
	rr4 = np.sqrt(np.sum((pos4-prop['cc'])**2, axis=1))
	mas4 = f['PartType'+str(j)+'/Masses'][:]
	ids4 = np.where(rr4<2*prop['r500'])[0]
	mas4 = mas4[ids4]*1.0e10
	rr4 = rr4[ids4]
	pot4 = pot4[ids4]
	vel4 = vel4[ids4]

	hubblev4 = pos4[ids4] / 0.6777 * 67.77 / 1e3 #Hubble flow; unit: km/s
	vv4 = vel4 + hubblev4
	# kin_energy_star=np.sum(vel4**2,axis=1,dtype=np.float64)*mas4/2.
	#--------------------------------------------------------------------------------------
	#--------------------------------------------------------------------------------------
	#select DM mass
	j = 1
	vel1 = f['PartType'+str(j)+'/Velocities'][:]
	pos1 = f['PartType'+str(j)+'/Coordinates'][:]
	pot1 = f['PartType'+str(j)+'/Potential'][:]
	rr1 = np.sqrt(np.sum((pos1-prop['cc'])**2, axis=1))
	mas1 = f['PartType'+str(j)+'/Masses'][:]
	ids1 = np.where(rr1<2*prop['r500'])[0]
	mas1 = mas1[ids1]*1.0e10
	rr1 = rr1[ids1]
	pot1 = pot1[ids1]
	vel1 = vel1[ids1]
	hubblev1 = pos1[ids1] / 0.6777 * 67.77 / 1e3 #Hubble flow; unit: km/s
	vv1 = vel1 + hubblev1
	# kin_energy_dm=np.sum(vel1**2,axis=1,dtype=np.float64)*mas1/2.

	#--------------------------------------------------------------------------------------
	#--------------------------------------------------------------------------------------
	#select BH mass
	j = 5
	if 'PartType5' in f.keys():
		vel5 = f['PartType'+str(j)+'/Velocities'][:]
		pos5 = f['PartType'+str(j)+'/Coordinates'][:]
		pot5 = f['PartType'+str(j)+'/Potential'][:]
		rr5 = np.sqrt(np.sum((pos5-prop['cc'])**2, axis=1))
		mas5 = f['PartType'+str(j)+'/BH_Mass'][:]
		ids5 = np.where(rr5<2*prop['r500'])[0]
		mas5 = mas5[ids5]*1.0e10
		rr5 = rr5[ids5]
		pot5 = pot5[ids5]
		vel5 = vel5[ids5]
		hubblev5 = pos5[ids5] / 0.6777 * 67.77 / 1e3 #Hubble flow; unit: km/s
		vv5 = vel5 + hubblev5

	# #--------------------------------------------------------------------------------------
	# #--------------------------------------------------------------------------------------
	# f.close()
	# #log space bins
	nbin = 50
	rbins=np.logspace(np.log10(0.001*prop['r500']), np.log10(2*prop['r500']), num=nbin+1)
	vol=4*np.pi*rbins**3/3
	vol=vol[1:]-vol[:-1] #unit: (kpc/h)^3
	volcm=vol*const.kpc.to('cm').value**3 #unit: (cm/h)^3


	#--------------------------------------------------------------------------------------
	#--------------------------------------------------------------------------------------
	#satellite galaxies
	#!!! here we select r200
	# fahf = glob.glob('/data4/niftydata/TheThreeHundred/data/catalogues/AHF/GIZMO/%s/GIZMO-%s.%s.z*.AHF_halos' %(cn,cn,snapn))
	# ahf = np.loadtxt(fahf[0])
	# idv = np.where(ahf[:,1] == hostid)[0]
	# subhalo_x = ahf[idv,5:8]
	# rrsubhalo = np.sqrt(np.sum((subhalo_x-prop['cc'])**2, axis=1))
	# rrsubhalo2D = np.sqrt(np.sum((subhalo_x[:,:2]-prop['cc'][:2])**2, axis=1))

	# rbins20=np.logspace(np.log10(0.001*prop['r200']), np.log10(1.5*prop['r200']), num=20+1)
	# vol20=4*np.pi*rbins20**3/3
	# vol20=vol20[1:]-vol20[:-1] #unit: (kpc/h)^3
	# vol20_2D = np.pi*rbins20**2
	# vol20_2D=vol20_2D[1:]-vol20_2D[:-1] #unit: (kpc/h)^2
	# numgal,xngal = np.histogram(rrsubhalo, bins=rbins20)
	# numgal20,xngal20 = np.histogram(rrsubhalo2D, bins=rbins20)

	#--------------------------------------------------------------------------------------
	#--------------------------------------------------------------------------------------
	# data=np.zeros((dens.size,4),dtype=np.float64)
	# data[:,0] = (rad[1:]+rad[:-1])/2
	num0,xn0 = np.histogram(rr, bins=rbins)
	num1,xn1 = np.histogram(rr1, bins=rbins)
	num4,xn4 = np.histogram(rr4, bins=rbins)
	if 'PartType5' in f.keys():
		num5,xn5 = np.histogram(rr5, bins=rbins)
		numdis = np.vstack((num0,num4,num0+num1+num4+num5))
		tmass5,xe5 = np.histogram(rr5, bins=rbins, weights=mas5)
	else:
		numdis = np.vstack((num0,num4,num0+num1+num4))

	temp = temp *  KtokeV #
	tmass,xe = np.histogram(rr, bins=rbins, weights=mas)
	tmass4,xe4 = np.histogram(rr4, bins=rbins, weights=mas4)
	tmass1,xe1 = np.histogram(rr1, bins=rbins, weights=mas1)
	avetem,xe = np.histogram(rr, bins=rbins, weights=temp) 
	tem,xe = np.histogram(rr, bins=rbins, weights=mas*temp) #mass-weighted==MW
	met,xe = np.histogram(rr, bins=rbins, weights=mas*metM)
	tele,xe = np.histogram(rr, bins=rbins, weights=nele)
	
	#--------------------------------------------------------------------------------------
	#--------------------------------------------------------------------------------------
	#===save data===
	np.seterr(divide='ignore', invalid='ignore')
	if 'Rbin' in sngrp.keys():  #Rbins
		sngrp['Rbin'][:]=(rbins[1:] + rbins[:-1])/2
	else:
		sngrp.create_dataset('Rbin', data = (rbins[1:] + rbins[:-1])/2) 

	#particle number in bins
	#shape:(3,nbin); 1: gas; 2: star; 3: total(dm+gas+star+BH)
	if 'NuminBin' in sngrp.keys():  #Rbins
		sngrp['NuminBin'][:]=numdis
	else:
		sngrp.create_dataset('NuminBin', data = numdis)

	# #satellite galaxies number density
	# #[0]: 3D; [1]: 2D
	# # del sngrp['Numgal']
	# dnumgal = np.vstack((numgal/vol20, numgal20/vol20_2D))
	# if 'Numgal' in sngrp.keys():  #Rbins
	# 	sngrp['Numgal'][:]=dnumgal
	# else:
	# 	sngrp.create_dataset('Numgal', data = dnumgal, dtype = 'f')  

	# gas density
	#unit: Msun/h/(kpc/h)^3
	if 'Gdens' in sngrp.keys():
		sngrp['Gdens'][:]=tmass/vol
	else:
		sngrp.create_dataset('Gdens', data = tmass/vol)   

	# Mass-Weighted temperature
	#unit: keV
	if 'MWTemp' in sngrp.keys():
		sngrp['MWTemp'][:]=tem/tmass
	else:
		sngrp.create_dataset('MWTemp', data = tem/tmass)   

	# pressure: P = \sum (k_B/\mu/m_p)\rho_i*T_i (Planelles+2017)
	#unit: keV * (cm/h)^(-3)
	#mass-weighted
	if 'Pressure' in sngrp.keys():
		sngrp['Pressure'][:]=tem/volcm/0.6125/m_ptoMsun/0.6777
	else:
		sngrp.create_dataset('Pressure', data = tem/volcm/0.6125/m_ptoMsun/0.6777)

	# Mass-Weighted metallicity
	#unit: Zsun
	if 'MWMetal' in sngrp.keys():
		sngrp['MWMetal'][:]=met/tmass/0.0134
	else:
		sngrp.create_dataset('MWMetal', data = met/tmass/0.0134) 

	# electron number density 
	# unit: {cm/h}^{-3}
	if 'Eledens' in sngrp.keys():  
		sngrp['Eledens'][:]=tele /volcm / m_ptoMsun
	else:
		sngrp.create_dataset('Eledens', data = tele/volcm/ m_ptoMsun)

	# Entropy: S = k_b * T / n_e^{2/3}
	#unit: [keV (cm/h)^2]
	if 'Entropy' in sngrp.keys():
		sngrp['Entropy'][:]= tem/tmass / (tele/volcm / m_ptoMsun)**(2/3)
	else:
		sngrp.create_dataset('Entropy', data = tem/tmass / (tele/volcm / m_ptoMsun)**(2/3)) 

	# stellar density in Msun/h/(kpc/h)^3
	if 'Stellardens' in sngrp.keys(): 
		sngrp['Stellardens'][:]=tmass4/vol
	else:
		sngrp.create_dataset('Stellardens', data = tmass4/vol)

	# total density in Msun/h/(kpc/h)^3
	if 'Totdens' in sngrp.keys(): 
		if 'PartType5' in f.keys():
			sngrp['Totdens'][:]=(tmass4+tmass+tmass1+tmass5)/vol
		else:
			sngrp['Totdens'][:]=(tmass4+tmass+tmass1)/vol
	else:
		if 'PartType5' in f.keys():
			sngrp.create_dataset('Totdens', data = (tmass4+tmass1+tmass+tmass5)/vol)
		else:
			sngrp.create_dataset('Totdens', data = (tmass4+tmass1+tmass)/vol)

#sel which snap in which region to compute and prepare basic host halo properties
def sel_sn(sn, h5fname, progenIDs, sn_sa_red):
	snapn='snap_%03d'%sn #which redshift...
	red=sn_sa_red[sn,2]

	for crn in range(1,325):
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

		#Now we use the cluster withID M500 
		fahf = glob.glob('/data6/users/aknebe/TMP/GIZMO_R500c/%s/GIZMO-%s.%s.z*.AHF_halos' %(cn,cn,snapn))
		ahf = np.loadtxt(fahf[0])
		idv = np.where((np.log10(ahf[:,3] / 0.678) >= 14) & (np.log10(ahf[:,3] / 0.678) <= 14.3) & (ahf[:,37] > 0.98))[0]

		if len(idv) < 1:
			print("Not find any satisfied clusters for ", crn, snapn)
			continue
			
		print('snapshot: %s, region: %s, matched clusters: %s' %(sn, crn, len(idv)))

		cc_all = ahf[idv, 5:8]
		r500_all = ahf[idv, 11]
		M500_all = ahf[idv, 3]
		for ni,seqgal in enumerate(idv):
			cc = cc_all[ni]
			r500 = r500_all[ni]
			M500 = M500_all[ni]
			print(M500, r500, cc)

			gal_key = 'gal_%s' %(seqgal+1)
			if gal_key in sngrp.keys(): #which snapshot
				galid = sngrp[gal_key]
			else:
				galid = sngrp.create_group(gal_key)

			galid.attrs['redshift'] = red; galid.attrs['snapshot'] = sn;
			galid.attrs['center_x'] = cc[0]; galid.attrs['center_y'] = cc[1]; galid.attrs['center_z'] = cc[2];
			galid.attrs['r500'] = r500
			galid.attrs['M500'] = M500

			prop = {'cc': cc, 'r500': r500}

			produce_prof(galid, cn, snapn, prop)

		fh5w.close()

		print('redshift: ', sn, 'region: ', crn)


#calculate profiles from GIZMO profiles in qingyang@brutus.ft.uam.es
def main():

	progenIDs=np.loadtxt("/home2/weiguang/Project-300-Clusters/Halo_mass_function_mass-difference/GIZMO/Progenitor-IDs-for-center-cluster.txt",dtype=np.int64)
	sn_sa_red=np.loadtxt('/home2/weiguang/Project-300-Clusters/redshifts.txt')

	##Open h5 file for write
	#check the data for reply of referee  
	h5fname="/home2/qingyang/GIZMO_prof/data/GIZ-profiles-progenitors_reply_massive.hdf5"
	if not os.path.isfile(h5fname):
		fh5w = h5py.File(h5fname, "w")
		fh5w.close()

	sel_sn(128, h5fname, progenIDs, sn_sa_red)
	sel_sn(96, h5fname, progenIDs, sn_sa_red)
	sel_sn(78, h5fname, progenIDs, sn_sa_red)
	# sel_sn(65, h5fname, progenIDs, sn_sa_red)
	# sel_sn(55, h5fname, progenIDs, sn_sa_red)

main()




