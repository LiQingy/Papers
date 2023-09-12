import numpy as np
import math
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15 as cosmo
from mpl_toolkits.mplot3d import Axes3D
from scipy import integrate

def comvdis(z):
	f = lambda x: 1 / (0.238*(1 + x)**3 + 0.762)**0.5
	return integrate.quad(f, 0, z)[0] * 3000 / 0.73

def position_xyz(xx2,xx3,z):
	degs = math.pi / 180
	dist = (cosmo.comoving_distance(z)).value * 1000 #kpc
	# dist = comvdis(z) * 1000 #kpc
	x1 = xx2 * degs 
	x2 = xx3 * degs
	unit = 1
	rx = (dist * math.cos(x2) * math.cos(x1) * unit)
	ry = (dist * math.cos(x2) * math.sin(x1) * unit)
	rz = (dist * math.sin(x2) * unit)

	return rx,ry,rz

def findgalaxy(filepC,fileip1,fileSDSS7,filest):
	petroC = np.loadtxt(filepC)
	ipetroC1 = np.loadtxt(fileip1)
	sdss7 = np.loadtxt(fileSDSS7)
	st = np.loadtxt(filest)

	sizepC = petroC.shape
	sizeip1 = ipetroC1.shape

	i = 0
	group = []
	while i < sizepC[0]:
		if 10**petroC[i][6] / 0.678 > 9.483e14:
			group.append(petroC[i])
		i += 1
	group = np.array(group) #all group information
	print(group.shape[0])
	
	j = 0
	gg =[]
	while j < group.shape[0]:
		k = 0
		while k < sizeip1[0]:
			if int(group[j][0]) == int(ipetroC1[k][2]):
				gg.append([group[j][0],ipetroC1[k][0],ipetroC1[k][1],ipetroC1[k][4]])
			k += 1		
		j += 1
	gg = np.array(gg)#1.groupid;2.galaxyid;3.galaxyid(NYU);4.whether most massive galaxy 
	ngalaxy = gg.shape[0]
	print('The number of galaxy is %s' %ngalaxy)

	#add galaxy mass and position 
	i1 = 0
	galaxy = []
	while i1 < ngalaxy:
		i2 = 0
		while i2 < sdss7.shape[0]:
			if gg[i1][1] == sdss7[i2][0]:
				galaxy.append([gg[i1][0],gg[i1][1],gg[i1][2],gg[i1][3],sdss7[i2][2],sdss7[i2][3],sdss7[i2][4],st[i2][3]])
				break
			i2 += 1
		i1 += 1
	galaxy = np.array(galaxy)#galaxy:1.groupid;2.galaxyid;3.galaxyid(NYU);4 whether most massive galaxy;5.RA;6.dec;7,z;8.M*


	print('galaxy is ok')
	return group,galaxy
#----------------------------------------------------------------------------

def allphy(fileam,filegalaxy):
	dataam = np.loadtxt(fileam)
	galaxy = np.loadtxt(filegalaxy)
	sdss7am = np.zeros(shape = (galaxy.shape[0],19))

	nzero = 0
	i = 0
	while i < galaxy.shape[0]:
		j = 0
		while j < dataam.shape[0]:
			if int(galaxy[i][1]) == int(dataam[j][0]):
				sdss7am[i][:8] = galaxy[i]
				sdss7am[i][8:] = dataam[j][5:16]
				break
			j += 1

		if j == dataam.shape[0]: #no age and metallicity
			sdss7am[i][:8] = galaxy[i]
			sdss7am[i][8:] = 0
		i += 1
	sdss7am = np.array(sdss7am)
	print('Adding age and metal is ok')
	return sdss7am


def findcen(group,galaxy,info,d):
	# fig = plt.figure()
	# ax = fig.add_subplot(111, projection='3d')

	sdss_group_center = []
	sizegroup = group.shape[0]
	sizegalaxy = galaxy.shape[0]
	
	pcr = 2.775 * 10**11 
	nn0 = 0 # the number group which exist no galaxy in 500kpc
	nn1 = 0 # the number group which exist galaxy in 500kpc
	nm  = 0
	nrall = 0

	rc500 = 0
	rcv = 0
	for j1 in range(sizegroup):
		groupid = group[j1][0]
		groupmass = group[j1][6] #halo mass based on group_L

		r200 = pow(10**groupmass * 3 / math.pi / 800/pcr, 1 / 3) * 1000 / 0.678#kpc 	

		rx,ry,rz = position_xyz(group[j1][1],group[j1][2],group[j1][3])#luminosity weghted group centre 
		rxg = rx
		ryg = ry
		rzg = rz
		rcen = math.sqrt(rxg**2 + ryg**2 + rzg**2)

		massmax = 0

		ju_BCG = -1
		ju_500 = -1
		mcx = 0
		mcy = 0
		mcz = 0
		mgmass = 0
		mcg = 0
		nr500 = 0
		nr200 = 0

		bigmass = 0
		rcr = 0

		k1 = 0
		while k1 < sizegalaxy:
			if galaxy[k1][0] == groupid and galaxy[k1][18] > 5* 10**10:
				nrall += 1
				galaxyid = int(galaxy[k1][1])
				rx,ry,rz = position_xyz(galaxy[k1][4],galaxy[k1][5],galaxy[k1][6])
				rgal = math.sqrt(rx**2 + ry**2 + rz**2)
				r1 = abs(rgal - rcen)
				r0 = math.sqrt((rxg - rx)**2 + (ryg - ry)**2 + (rzg - rz)**2)
				r0r = math.sqrt(abs(r0**2 - r1**2))
				# if groupid == 30:
					# print(r0r)
					# ax.scatter3D(rx-rxg,ry-ryg,rz-rzg,color = 'red')
				if galaxy[k1][18] > bigmass:
					bigmass = galaxy[k1][18] * 1
					rcr = r0r * 1
				if r0r < 500:
					nr500+=1
				if r0r > 0:
					#mass-weighted centoid
					mcx += galaxy[k1][4] * galaxy[k1][18]
					mcy += galaxy[k1][5] * galaxy[k1][18]
					mcz += galaxy[k1][6] * galaxy[k1][18]
					mgmass += galaxy[k1][18]
					nr200 += 1
					#BCGs
					if info[galaxyid - 1][22] < mcg:
						xcg = galaxy[k1][4]
						ycg = galaxy[k1][5]
						zcg = galaxy[k1][6]
						mcg = info[galaxyid - 1][22]
						ju_BCG = -2
					if galaxy[k1][18] > 5*10**10:
						nm += 1	
				#Massive galaxy	within 500kpc		
				if r0r < 100 and galaxy[k1][18] > massmax:
					massmax = galaxy[k1][18]
					xc = galaxy[k1][4]
					yc = galaxy[k1][5]
					zc = galaxy[k1][6]
					ju_500 = -2
			k1 += 1	
		# print(rcr)
		if rcr < 500:
			rc500 += 1
		if rcr > r200:
			rcv += 1
		if ju_500 == -2:
			sdss_group_center.append([groupid,group[j1][1],group[j1][2],group[j1][3],r200])	
			nn1 += 1
		elif ju_BCG == -2:
			# mcx = mcx / mgmass
			# mcy = mcy / mgmass
			# mcz = mcz / mgmass
			
			# sdss_group_center.append([groupid,mcx,mcy,mcz,r200])
			sdss_group_center.append([groupid,group[j1][1],group[j1][2],group[j1][3],r200])
			
			print(groupid,r200,nr200,nr500)
			nn0 += 1
	print('#',rc500,rcv)
	# u = np.linspace(0,2*np.pi,1000)
	# v = np.linspace(0,2*np.pi,1000)
	# x=500*np.outer(np.cos(u),np.sin(v))
	# y=500*np.outer(np.sin(u),np.sin(v))
	# z=500*np.outer(np.ones(np.size(u)),np.cos(v))
	# ax.plot_surface(x,y,z,color = 'b',alpha = 0.1)
	# ax.scatter3D(0,0,0,color = 'black')
	# # plt.savefig('/home/qyli/Desktop/group112_c.png')
	# plt.show()

	print('Already Finding all group center',nn0,nn1,nrall,nm)
	sdss_group_center = np.array(sdss_group_center)#rx,ry,rz is kpc
	return sdss_group_center



#-----------------------------------------------------------------------------

def phyprofile(n,groupc,sdss7am,sdss7sfr):
	# fig = plt.figure()
	# ax = fig.add_subplot(111, projection='3d')

	groupsize = groupc.shape[0]
	sdss7_totaln4 = np.zeros(shape = (groupsize,n))
	sdss7_starz = np.zeros(shape = (groupsize,n))
	sdss7_age = np.zeros(shape = (groupsize,n))
	nbin = 60
	sdss7_sfr = np.zeros(shape = (nbin,nbin))
	sdss7_ssfr = np.zeros(shape = (nbin,nbin))
	sdss7_density4 = np.zeros(shape = (groupsize,n))
	sdss7_rr4 = np.zeros(shape = (groupsize,n))

	totalm = 0
	total200 = 0
	for i in range(groupsize):
		xc,yc,zc = position_xyz(groupc[i][1],groupc[i][2],groupc[i][3])
		rcen = math.sqrt(xc**2 + yc**2 + zc**2)
		rv = groupc[i][4]
		h = rv / n

		mass4 = np.zeros(n)
		metal = np.zeros(n)
		age = np.zeros(n)
		nub = np.zeros(n)
		rr4 = np.zeros(n)

		j = 0
		# 	ax.scatter(0, 0, 0, marker = '^', color = 'b')
		while j < sdss7am.shape[0]:
			if sdss7am[j][0] == groupc[i][0] and sdss7am[j][18] > 5 * 10**10:
				rx,ry,rz = position_xyz(sdss7am[j][4],sdss7am[j][5],sdss7am[j][6])
				rgal = math.sqrt(rx**2 + ry**2 + rz**2)
				r1 = abs(rgal - rcen)
				r0 = math.sqrt((xc - rx)**2 + (yc - ry)**2 + (zc - rz)**2)
				rr = math.sqrt(r0**2 - r1**2)  
				totalm += 1
				# 	ax.scatter(rx - xc, ry - yc, rz - zc, color = 'r')
				# if groupc[i][0] == 112:
				# 	ax.scatter3D(rx-xc,ry-yc,rz-zc,color = 'r')
				if rr < rv:

					total200 += 1
					rn = int(rr / h)
					mass4[rn] = sdss7am[j][18] + mass4[rn]
					metal[rn] = sdss7am[j][18] * sdss7am[j][13] + metal[rn]
					age[rn] = sdss7am[j][18] * sdss7am[j][8] + age[rn]

					rr4[rn] = rr / rv * sdss7am[j][18] + rr4[rn]
					nub[rn] += 1

					#SFR
					gid = int(sdss7am[j][1])
					sfr = 10**sdss7sfr[gid - 1][3]
					

					if sfr < 0.001:
						sfr = np.log10(0.001)
						ssfr = -15
					else:
						sfr = np.log10(sfr)
						ssfr = sdss7sfr[gid - 1][6]
					nr = int(rr / ( rv / nbin))
					nsfr = int((sfr + 3) / 6 * nbin)
					nssfr = int((ssfr + 15) / 6 * nbin)
					sdss7_sfr[nsfr][nr] += 1
					sdss7_ssfr[nssfr][nr] += 1				
			j += 1

		# h = np.median(groupc[:,4]) / n
		for j in range(n):
			dv = 4. / 3 * math.pi *(pow((j+1)*h,3) - pow((j*h),3))
			sdss7_density4[i][j] = mass4[j] / dv
			sdss7_totaln4[i][j] = nub[j] / dv

		sdss7_starz[i] = metal / mass4
		sdss7_age[i] = age / mass4
		sdss7_rr4[i] = rr4 / mass4

	# u = np.linspace(0,2*np.pi,1000)
	# v = np.linspace(0,2*np.pi,1000)
	# x=2531*np.outer(np.cos(u),np.sin(v))
	# y=2531*np.outer(np.sin(u),np.sin(v))
	# z=2531*np.outer(np.ones(np.size(u)),np.cos(v))
	# ax.plot_surface(x,y,z,color = 'b',alpha = 0.1)
	# ax.scatter(0,0,0,color = 'blue')
	# plt.savefig('/home/qyli/Desktop/group112_MC.png')
	# plt.show()

	print('The number of groups is %s' %groupsize)
	print('The number of galaxy satisfing stellar mass condition is %s' %totalm)
	print('The number of galaxy in R200 is %s' %total200)	

	return sdss7_density4,sdss7_age,sdss7_starz,sdss7_sfr,sdss7_ssfr,sdss7_rr4,sdss7_totaln4
#------------------------------------------------------------------------------
def main(n,d,savedata):
	#get the center of group and information about galaxy
	if n == 1:
				
		filepC = '/home/qyli/Desktop/SDSS7/SDSS7_realspace/reReal_Plack15/modelC_L_m'
		fileip1 = '/home/qyli/Desktop/SDSS7/SDSS7_realspace/reReal_Plack15/imodelC_1'
		fileSDSS7 = '/home/qyli/Desktop/SDSS7/SDSS_DR7/SDSS7'
		filest = '/home/qyli/Desktop/SDSS7/SDSS_DR7/SDSS7_ST'
		group,galaxy = findgalaxy(filepC,fileip1,fileSDSS7,filest)

		np.savetxt('/home/qyli/Desktop/300data/sdss7/Planck/galaxy_Planck.txt',galaxy)
		# galaxy:1.groupid;2.galaxyid;3.galaxyid(NYU);4 whether most massive galaxy;5.RA;6.dec;7,z;8.M*
		np.savetxt('/home/qyli/Desktop/300data/sdss7/Planck/group_Planck.txt',group)
		#1.all information with petroC_group
	#------------------------------------------
	# add the age and metallicity into the galaxy
	if n == 2:
	
		fileam = '/home/qyli/Desktop/300data/sdss7/sdss_firefly-26_SDSS7_sort.txt'
		filegalaxy = '/home/qyli/Desktop/300data/sdss7/Planck/galaxy_Planck.txt' 
		sdss7am = allphy(fileam,filegalaxy)
		if savedata == True:
			np.savetxt('/home/qyli/Desktop/300data/sdss7/Planck/sdss7am_Planck.txt',sdss7am)
		# 1.groupid;2.galaxyid;3.galaxyid(NYU);4 whether most massive galaxy;5.RA;6.dec;7,z;8.M*;9-13.age;14-18.metallicity;19.stellar Mass
	#-------------------------------------------
	#Find cluster center
	if n == 2.1:
		filegalaxy = '/home/qyli/Desktop/300data/sdss7/Planck/sdss7am_Planck.txt' 
		filegroup = '/home/qyli/Desktop/300data/sdss7/Planck/group_Planck.txt'
		fileinfo = '/home/qyli/Desktop/SDSS7/SDSS_DR7/SDSS7_INFO'
		group = np.loadtxt(filegroup)
		galaxy = np.loadtxt(filegalaxy)
		info = np.loadtxt(fileinfo)
				
		sdss_group_center = findcen(group,galaxy,info,d)
		print(sdss_group_center.shape)

		if savedata == True:
			np.savetxt('/home/qyli/Desktop/300data/sdss7/Planck/sdss_group_centerPlanckBCG.txt',sdss_group_center)
		# 1.groupid;2.xc;3.yc;4.zc;5.r200	
	#-----------------------------------------------
	# get the starz of SDSS7
	if n == 3:
		filegroupc = '/home/qyli/Desktop/300data/sdss7/Planck/sdss_group_centerPlanckBCG.txt' 
		fileam = '/home/qyli/Desktop/300data/sdss7/Planck/sdss7am_Planck.txt' 
		filesfr = '/home/qyli/Desktop/SDSS7/SDSS7/SDSS7_SFR'
		groupc = np.loadtxt(filegroupc)
		print(groupc.shape)
		sdss7am = np.loadtxt(fileam)
		sdss7sfr = np.loadtxt(filesfr)

		sdss7_density4,sdss7_age,sdss7_starz,sdss7_sfr,sdss7_ssfr,sdss7_rr4,sdss7_totaln4 = phyprofile(20,groupc,sdss7am,sdss7sfr)
		if savedata == True:
			# np.savetxt('/home/qyli/Desktop/300data/sdss7/Planck/sdss7_density4.txt' ,sdss7_density4) #Msun / kpc^-3
			# np.savetxt('/home/qyli/Desktop/300data/sdss7/Planck/sdss7_age.txt' ,sdss7_age) #yr
			# np.savetxt('/home/qyli/Desktop/300data/sdss7/Planck/sdss7_starz.txt' ,sdss7_starz) #Zsun
			# np.savetxt('/home/qyli/Desktop/300data/sdss7/Planck/sdss7_sfr.txt' ,sdss7_sfr) #Msun * yr^-1
			# np.savetxt('/home/qyli/Desktop/300data/sdss7/Planck/sdss7_ssfr.txt' ,sdss7_ssfr) #yr^-1
			np.savetxt('/home/qyli/Desktop/300data/sdss7/Planck/sdss7_rr4PlanckBCG.txt' ,sdss7_rr4) #R/R200
			np.savetxt('/home/qyli/Desktop/300data/sdss7/Planck/sdss7_totaln4PlanckBCG.txt' ,sdss7_totaln4) 

# main(1,d = 500,savedata = True)
# main(2,d = 500,savedata = True)
main(2.1,d = 500,savedata = True)
# main(3,d = 500,savedata = True)




				







