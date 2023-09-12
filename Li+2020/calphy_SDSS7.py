import numpy as np
import math
from astropy.cosmology import Planck15 as cosmo

def position_xyz(xx2,xx3,z):
	degs = math.pi / 180
	dist = (cosmo.comoving_distance(z)).value * 1000 #kpc
	x1 = xx2 * degs 
	x2 = xx3 * degs
	unit = 1
	rx = (dist * math.cos(x2) * math.cos(x1) * unit)
	ry = (dist * math.cos(x2) * math.sin(x1) * unit)
	rz = (dist * math.sin(x2) * unit)

	return rx,ry,rz

def phyprofile(n,groupc,sdss7am,sdss7sfr,galtype,zmax):
	groupsize = groupc.shape[0]
	starn4 = sdss7_totaln4 = np.zeros(shape = (groupsize,n),dtype = 'f8')
	sdss7_totaln4 = np.zeros(shape = (groupsize,n))
	sdss7_density4 = np.zeros(shape = (groupsize,n))
	sdss7_rr4 = np.zeros(shape = (groupsize,n))
	mass4 = np.zeros(shape = (groupsize,n))

	sdss7_starz = np.zeros(shape = (groupsize,n))
	sdss7_age = np.zeros(shape = (groupsize,n))
	sdss7_rr = np.zeros(shape = (groupsize,n))
	mass4_age_z = np.zeros(shape = (groupsize,n))
	
	nbin = 60
	sdss7_sfr = np.zeros(shape = (nbin,nbin))
	sdss7_ssfr = np.zeros(shape = (nbin,nbin))
	sdss7_age2D = np.zeros(shape = (nbin,nbin))
	sdss7_starz2D = np.zeros(shape = (nbin,nbin))
	sdss7_age2Done = []
	sdss7_starz2Done = []
	sdss7_ageBCG = []
	sdss7_starzBCG = []
	sdss7_sfr2Done = []
	sdss7_ssfr2Done = []
	sdss7_sfrBCG = []
	sdss7_ssfrBCG = []



	ntotal200 = 0
	ntotalM = 0
	nubage = 0
	nextra = 0
	nmass = 0
	nBCG = 0
	zall = 0
	nsatellite = 0
	ninage = 0

	j = 0
	while j < sdss7am.shape[0]:
		seqgroup = 0
			# ntotalM += 1
		groupid = sdss7am[j][0]
		for zz in range(groupc.shape[0]):
			if groupc[zz][0] == groupid:
				seqgroup = zz *1 
				zz = -1
				break
		locn = np.where(sdss7am[j][1] == galtype[:,0])[0]
		startype = galtype[locn,2]

		groupcx = groupc[seqgroup][1]
		groupcy = groupc[seqgroup][2]
		groupcz = groupc[seqgroup][3]
		xc,yc,zc = position_xyz(groupcx,groupcy,groupcz)
		rx,ry,rz = position_xyz(sdss7am[j][4],sdss7am[j][5],sdss7am[j][6])
		if rx == xc and ry == yc and rz ==zc and sdss7am[j][18] < 5 *10**10:
			nBCG += 1
			print(sdss7am[j][18])
		if sdss7am[j][18] > 5 *10**10 and startype == 1:
				
			if zz == -1:
				nmass += 1
				groupcx = groupc[seqgroup][1]
				groupcy = groupc[seqgroup][2]
				groupcz = groupc[seqgroup][3]
				rv = groupc[seqgroup][4]
				xc,yc,zc = position_xyz(groupcx,groupcy,groupcz)
				rx,ry,rz = position_xyz(sdss7am[j][4],sdss7am[j][5],sdss7am[j][6])

				rcen = math.sqrt(xc**2 + yc**2 + zc**2)
				rgal = math.sqrt(rx**2 + ry**2 + rz**2)
				r1 = abs(rgal - rcen)
				r0 = math.sqrt((xc - rx)**2 + (yc - ry)**2 + (zc - rz)**2)
				rr = math.sqrt(r0**2 - r1**2)
				

				if rr == 0:
					if sdss7am[j][8] * 1e-9 < 13.798:
						sdss7_ageBCG.extend([sdss7am[j][8]*1e-9])
					sdss7_starzBCG.extend([sdss7am[j][13]* 0.02 / 0.0134])
					ntotalM += 1
					# ntotal200 += sdss7am[j][18]
					# zall += sdss7am[j][6]* sdss7am[j][18]

					sfr = 10**sdss7sfr[galaxyid - 1][3]
					
					if sfr < 0.001:
						sfr = np.log10(0.001)
						ssfr = -15
					else:
						sfr = np.log10(sfr)
						ssfr = sdss7sfr[galaxyid - 1][6]
					sdss7_sfrBCG.extend([sfr])
					sdss7_ssfrBCG.extend([ssfr])

				if 0 < rr < rv:
					h = rv / n
					rn = int(rr / h)
					nsatellite += 1
					ntotal200 += sdss7am[j][18]
					zall += sdss7am[j][6] * sdss7am[j][18]
					galaxyid = int(sdss7am[j][1])
					groupid = int(groupc[seqgroup][0])
					#------------------------------
					zgalaxy = zmax[galaxyid-1][2]
					nclu = 0
					for k in range(groupc.shape[0]):
						if groupc[k][3] <= zgalaxy:
							nclu += 1
						if groupc[k][0] == groupid:
							sdss7_totaln4[k][rn] += 1
							mass4[k][rn] += sdss7am[j][18]
							sdss7_rr4[k][rn] += rr / rv * sdss7am[j][18]

							sdss7_age[k][rn] += sdss7am[j][18] * sdss7am[j][8]
							sdss7_starz[k][rn] += sdss7am[j][18] * sdss7am[j][13]
							mass4_age_z[k][rn] += sdss7am[j][18]
							sdss7_rr[k][rn] += rr / rv * sdss7am[j][18]
					
					for k in range(groupc.shape[0]):
						if groupc[k][3] > zgalaxy and groupc[k][0] != groupid:
							sdss7_totaln4[k][rn] += 1/nclu
							mass4[k][rn] += sdss7am[j][18] / nclu
							sdss7_rr4[k][rn] += rr / rv * sdss7am[j][18] / nclu
							nextra += 1/nclu
					#--------------------------------------------------------------
					#age2D and metallicity2D
					nr = int(rr / ( rv / nbin))
					if sdss7am[j][8] * 1e-9 < 13.789:
						ninage += 1
						nage = int((sdss7am[j][8] * 1e-9 - 4) / 10 * nbin)
						if nage >= 60:
							nubage += 1
						if 60 > nage >= 0 :	
							sdss7_age2D[nage][nr] += 1
							sdss7_age2Done.extend([sdss7am[j][8]*1e-9])					
					nstarz = int(sdss7am[j][13] * 0.02 / 0.0134 / 3 * nbin)
					if nstarz < nbin:
						sdss7_starz2D[nstarz][nr] += 1
						sdss7_starz2Done.extend([sdss7am[j][13]* 0.02 / 0.0134])
						#SFR
					sfr = 10**sdss7sfr[galaxyid - 1][3]
				
					if sfr < 0.001:
						sfr = np.log10(0.001)
						ssfr = -15
					else:
						sfr = np.log10(sfr)
						ssfr = sdss7sfr[galaxyid - 1][6]
					nr = int(rr / ( rv / nbin))
					nsfr = int((sfr + 3) / 6 * nbin)
					nssfr = int((ssfr + 15) / 6 * nbin)
					sdss7_sfr[nsfr][nr] += 1
					sdss7_ssfr[nssfr][nr] += 1	
					
					sdss7_sfr2Done.extend([sfr])
					sdss7_ssfr2Done.extend([ssfr])
		j += 1

	starn4 = sdss7_totaln4 * 1
	sdss7_totaln4ac = np.zeros(shape = (92,20))
	sdss7_mass4 = np.zeros(shape = (92,20))
	sdss7_density4ac = np.zeros(shape = (92,20))
	for i in range(groupc.shape[0]):
		h = groupc[i][4] / n
		for j in range(n):
			# dv = 4. / 3 * math.pi *(pow((j+1)*h,3) - pow((j*h),3))
			dv = math.pi *(pow((j+1)*h,2))
			sdss7_totaln4ac[i][j] = np.sum(sdss7_totaln4[i][:j+1]) / dv
			sdss7_density4ac[i][j] = np.sum(mass4[i][:j+1]) / dv
			sdss7_mass4[i][j] = np.sum(mass4[i][:j+1])
	sdss7_rr4 = sdss7_rr4 / mass4

	sdss7_age = sdss7_age / mass4_age_z
	sdss7_starz = sdss7_starz / mass4_age_z
	sdss7_rr = sdss7_rr / mass4_age_z
	sdss7_age2Done = np.array(sdss7_age2Done)
	sdss7_starz2Done = np.array(sdss7_starz2Done)
	sdss7_ageBCG = np.array(sdss7_ageBCG)
	sdss7_starzBCG = np.array(sdss7_starzBCG)
	sdss7_sfr2Done = np.array(sdss7_sfr2Done)
	sdss7_ssfr2Done = np.array(sdss7_ssfr2Done)
	sdss7_sfrBCG = np.array(sdss7_sfrBCG)
	sdss7_ssfrBCG = np.array(sdss7_ssfrBCG)
	mass200gal = np.sum(mass4,axis = 1) * 1

	print(nmass,ntotal200,nBCG,ntotalM,nubage,nextra)
	print('satellite galaxies within rv is %s' %nsatellite)
	print('mass weighted galaxies average redshift is %s' %(zall / ntotal200))
	print('the galaxies number which age are smaller than 13.789 is %s' %ninage)
	return sdss7_mass4,mass200gal,sdss7_density4ac,sdss7_age,sdss7_starz,sdss7_sfr,sdss7_ssfr,sdss7_rr4,sdss7_totaln4ac,sdss7_rr,sdss7_age2D,sdss7_starz2D,sdss7_age2Done,sdss7_starz2Done,sdss7_ageBCG,sdss7_starzBCG,sdss7_sfr2Done,sdss7_ssfr2Done,sdss7_sfrBCG,sdss7_ssfrBCG,starn4

def main(savedata):

	filegroupc = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss_group_centerPlanckBCG.txt' 
	fileam = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7am_Planck.txt' 
	filesfr = '/Users/liqy/Documents/data/SDSS7/SDSS7/SDSS7_SFR'
	filezmax = '/Users/liqy/Documents/data/SDSS7/Zmax/SDSS7_model_zmax.dat'
	filegaltype = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_firefly_matchgal.txt'
	groupc = np.loadtxt(filegroupc)
	print(groupc.shape)
	sdss7am = np.loadtxt(fileam)
	print(sdss7am.shape)
	sdss7sfr = np.loadtxt(filesfr)
	galtype = np.loadtxt(filegaltype)
	zmax = np.loadtxt(filezmax)

	sdss7_mass4,mass200gal,sdss7_density4ac,sdss7_age,sdss7_starz,sdss7_sfr,sdss7_ssfr,sdss7_rr4,sdss7_totaln4ac,sdss7_rr,sdss7_age2D,sdss7_starz2D,sdss7_age2Done,sdss7_starz2Done,sdss7_ageBCG,sdss7_starzBCG,sdss7_sfr2Done,sdss7_ssfr2Done,sdss7_sfrBCG,sdss7_ssfrBCG,starn4 = phyprofile(20,groupc,sdss7am,sdss7sfr,galtype,zmax)
	if savedata == True:
		print('#')
		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_starn4PlanckBCG.txt' ,starn4)
		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_mass4PlanckBCG.txt' ,sdss7_mass4)
		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_density4PlanckBCG.txt' ,sdss7_density4ac) #Msun / kpc^-3
		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_mass200gal.txt' ,mass200gal) #Msun / kpc^-3
		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_agePlanckBCG.txt' ,sdss7_age) #yr
		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_rrPlanckBCG.txt' ,sdss7_rr) #yr

		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_starzPlanckBCG.txt' ,sdss7_starz) #Zsun
		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_sfrPlanckBCG.txt' ,sdss7_sfr) #Msun * yr^-1
		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_ssfrPlanckBCG.txt' ,sdss7_ssfr) #yr^-1
		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_age2DPlanckBCG.txt' ,sdss7_age2D)
		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_starz2DPlanckBCG.txt' ,sdss7_starz2D)
		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_age2DonePlanckBCG.txt' ,sdss7_age2Done)
		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_starz2DonePlanckBCG.txt' ,sdss7_starz2Done)
		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_age2DoneBCG.txt' ,sdss7_ageBCG)
		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_starz2DoneBCG.txt' ,sdss7_starzBCG)

		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_sfr2Done.txt' ,sdss7_sfr2Done)
		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_ssfr2Done.txt' ,sdss7_ssfr2Done)
		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_sfrBCG.txt' ,sdss7_sfrBCG)
		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_ssfrBCG.txt' ,sdss7_ssfrBCG)

		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_rr4PlanckBCG.txt' ,sdss7_rr4) #R/R200
		# np.savetxt('/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_totaln4PlanckBCG.txt' ,sdss7_totaln4ac) 

main(savedata = True)