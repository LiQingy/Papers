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
	m_ptoMsun = const.m_p.to('Msun').value #proton mass in unit Msun
	KtokeV = const.k_B.value / const.e.value * 1e-3 # [K] to [keV] == k_B*T

	j = 0 #particle type: 0: gas; 1: DM; 4: stars
	vel = readsgl.readsnapsgl(filename, block = 'VEL ',ptype = j)
	pos = readsgl.readsnapsgl(filename, block = 'POS ',ptype = j)
	rr = np.sqrt(np.sum((pos-prop['cc'])**2, axis=1))
	mas = readsgl.readsnapsgl(filename, block = 'MASS',ptype = j)
	sgden=readsgl.readsnapsgl(filename, block = 'RHO ',ptype = j)* 1.0e10 /0.6777/(1/0.6777)**3  #in Msun/kpc^3
	# temp= f['PartType'+str(j)+'/InternalEnergy'][:] * (5. / 3 - 1) * v_unit**2 * prtn * mean_mol_weight / bk # in k
	temp = readsgl.readsnapsgl(filename, block = 'TEMP',ptype = j)
	eleab= readsgl.readsnapsgl(filename, block = 'NE  ',ptype = j)
	sfr = readsgl.readsnapsgl(filename, block = 'SFR ',ptype = j)
	metl= readsgl.readsnapsgl(filename, block = 'Z   ',ptype = j)
	metHe = 0.24 # fixed He abundance 
	metM = metl #metallicity
	nH = (1 - metHe - metM)*mas *1.0e10 / 0.6777  #proton number 
	nele = nH * eleab #electron number: proton number multiply relative value

	#particle conditions
	ids = np.where((temp>1.0e6)&(rr<1.5*prop['r200'])&(sgden<2.88e6))[0] #select satisfied particles
	
	rr=rr[ids]
	mas=mas[ids]*1.0e10
	temp=temp[ids]
	nele=nele[ids]
	metM=metM[ids]
	vel=vel[ids]
	hubblev = pos[ids] / 0.6777 * 67.77 / 1e3 #Hubble flow; unit: km/s
	vv = vel + hubblev

	#select star mass
	j = 4
	pos4 = readsgl.readsnapsgl(filename, block = 'POS ',ptype = j)
	vel4 = readsgl.readsnapsgl(filename, block = 'VEL ',ptype = j)
	rr4 = np.sqrt(np.sum((pos4-prop['cc'])**2, axis=1))
	mas4 = readsgl.readsnapsgl(filename, block = 'MASS',ptype = j)
	metl4= readsgl.readsnapsgl(filename, block = 'Z   ',ptype = j)
	ids4 = np.where(rr4<1.5*prop['r200'])[0]
	mas4 = mas4[ids4]*1.0e10
	rr4 = rr4[ids4]
	vel4 = vel4[ids4]
	metM4 = metl4[ids4]

	hubblev4 = pos4[ids4] / 0.6777 * 67.77 / 1e3 #Hubble flow; unit: km/s
	vv4 = vel4 + hubblev4

	#select DM mass
	j = 1
	pos1 = readsgl.readsnapsgl(filename, block = 'POS ',ptype = j)
	vel1 = readsgl.readsnapsgl(filename, block = 'VEL ',ptype = j)
	rr1 = np.sqrt(np.sum((pos1-prop['cc'])**2, axis=1))
	mas1 = np.tile(0.12691148, pos1.shape[0])
	ids1 = np.where(rr1<1.5*prop['r200'])[0]
	mas1 = mas1[ids1]*1.0e10
	rr1 = rr1[ids1]
	vel1=vel1[ids1]

	hubblev1 = pos1[ids1] / 0.6777 * 67.77 / 1e3 #Hubble flow; unit: km/s
	vv1 = vel1 + hubblev1

	#select BH mass
	j = 5
	pos5 = readsgl.readsnapsgl(filename, block = 'POS ',ptype = j)
	vel5 = readsgl.readsnapsgl(filename, block = 'VEL ',ptype = j)
	rr5 = np.sqrt(np.sum((pos5-prop['cc'])**2, axis=1))
	mas5 = readsgl.readsnapsgl(filename, block = 'MASS',ptype = j)
	ids5 = np.where(rr5<1.5*prop['r200'])[0]
	if len(ids5):
		mas5 = mas5[ids5]*1.0e10
		rr5 = rr5[ids5]
		vel5 = vel5[ids5]

		hubblev5 = pos5[ids5] / 0.6777 * 67.77 / 1e3 #Hubble flow; unit: km/s
		vv5 = vel5 + hubblev5
	else:
		print(cn,snapn,mas5,ids5)

	#log space bins
	nbin = 50
	rbins=np.logspace(np.log10(0.001*prop['r200']), np.log10(1.5*prop['r200']), num=nbin+1)
	vol=4*np.pi*rbins**3/3
	vol=vol[1:]-vol[:-1] #unit: (kpc/h)^3
	volcm=vol*const.kpc.to('cm').value**3 #unit: (cm/h)^3

	#--------------------------------------------------------------------
	#--------------------------------------------------------------------
	#satellite galaxies
	fahf = glob.glob('/home/nifty2014/TheThreeHundred/catalogues/AHF/GadgetX/%s/GadgetX-%s.%s.z*.AHF_halos' %(cn,cn,snapn))
	ahf = np.loadtxt(fahf[0])
	idv = np.where(ahf[:,1] == hostid)[0]
	subhalo_x = ahf[idv,5:8]
	rrsubhalo = np.sqrt(np.sum((subhalo_x-prop['cc'])**2, axis=1))
	numgal,xngal = np.histogram(rrsubhalo, bins=rbins)

	# data=np.zeros((dens.size,4),dtype=np.float64)
	# data[:,0] = (rad[1:]+rad[:-1])/2
	num0,xn0 = np.histogram(rr, bins=rbins)
	num1,xn1 = np.histogram(rr1, bins=rbins)
	num4,xn4 = np.histogram(rr4, bins=rbins)
	if len(ids5):
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
	met4,xe = np.histogram(rr4, bins=rbins, weights=mas4*metM4)
	tele,xe = np.histogram(rr, bins=rbins, weights=nele)


	np.seterr(divide='ignore', invalid='ignore')
	if 'Rbin' in sngrp.keys():  #Rbins
		sngrp['Rbin'][:]=(rbins[1:] + rbins[:-1])/2
	else:
		sngrp.create_dataset('Rbin', data = (rbins[1:] + rbins[:-1])/2)  

	#particle number in bins
	#1: gas; 2: star; 3: total(dm+gas+star+BH)
	if 'NuminBin' in sngrp.keys():  #Rbins
		sngrp['NuminBin'][:]=numdis
	else:
		sngrp.create_dataset('NuminBin', data = numdis) 
        
        #satellite galaxies number density
	# del sngrp['Numgal']
	if 'Numgal' in sngrp.keys():  #Rbins
		sngrp['Numgal'][:]=numgal/vol
	else:
		sngrp.create_dataset('Numgal', data = numgal/vol, dtype = 'f')

	# gas density
	#unit: Msun/h/(kpc/h)^3
	if 'Gdens' in sngrp.keys():
		sngrp['Gdens'][:]=tmass/vol
	else:
		sngrp.create_dataset('Gdens', data = tmass/vol)   

	# MW temp
	#unit: keV
	if 'MWTemp' in sngrp.keys():
		sngrp['MWTemp'][:]=tem/tmass
	else:
		sngrp.create_dataset('MWTemp', data = tem/tmass)   

	# pressure: P = \sum (k_B/\mu/m_p)\rho_i*T_i
	#unit: keV * (cm/h)^(-3)
	#mass-weighted
	if 'Pressure' in sngrp.keys():
		sngrp['Pressure'][:]=tem/volcm/0.6125/m_ptoMsun/0.6777
	else:
		sngrp.create_dataset('Pressure', data = tem/volcm/0.6125/m_ptoMsun/0.6777)

	# MW metallicity
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
		sngrp.create_dataset('Eledens', data = tele/volcm/ m_ptoMsun )

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
		if len(ids5):
			sngrp['Totdens'][:]=(tmass4+tmass+tmass1+tmass5)/vol
		else:
			sngrp['Totdens'][:]=(tmass4+tmass+tmass1)/vol
	else:
		if len(ids5):
			sngrp.create_dataset('Totdens', data = (tmass4+tmass1+tmass+tmass5)/vol)
		else:
			sngrp.create_dataset('Totdens', data = (tmass4+tmass1+tmass)/vol)	



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


		prop = {'cc':cc, 'r200':r200}
		
		produce_prof(sngrp, red, cn, snapn, prop, progenIDs[crn-1, sn])

		fh5w.close()

		print('redshift: ', sn, 'region: ', crn)

#calculate profiles from Gadget-X profiles in nifty2014@castor.ft.uam.es
#consider evolution profiles (2021.4.13)
def main():

	progenIDs=np.loadtxt('/home/nifty2014/TheThreeHundred/playground/qingyang/GadgetX_prof/data/GadgetX/Progenitor-IDs-for-center-cluster.txt',dtype=np.int64)
	sn_sa_red=np.loadtxt('/home/nifty2014/TheThreeHundred/playground/weiguang/redshifts.txt')
	
	with open('/home/nifty2014/TheThreeHundred/playground/qingyang/GadgetX_prof/readme.txt') as freadme:
  		readme = freadme.read()

	##Open h5 file for write
	h5fname="/home/nifty2014/TheThreeHundred/playground/qingyang/GadgetX_prof/data/G3X-profiles-progenitors_myevoz.hdf5"
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

# Velocity dispersion !!add hubble flow !! minus halo bulk Velocity has no effects on sigma_v!!
#sigam = std(v)= <(v-<v>)^2> = (given <v> = sum(m_i*v_i)/M) ... prove omit = <v^2> - <v>^2
# here <v^2> = sum(m_i*v_i^2)/M, Must in 3D sigma = np.sqrt(sigma_x^2+sigma_y^2+sigma_z^2)
        # hvel    =svel+scfa*Hubblez*spos/HubbleParam
        # hbv     =np.sum(hvel[:idr200]*np.tile(smas[:idr200],(3,1)).T, axis=0)/m200
        # sgm_200=np.sqrt(np.sum((hvel[:idr200]-hbv)**2*np.tile(smas[:idr200],(3,1)).T)/m200)
