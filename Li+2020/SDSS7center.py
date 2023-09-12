import numpy as np
import math
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15 as cosmo
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d import Axes3D

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

def BCGloc(groupid,galaxy,info,d):
	k2 = 0
	sizegalaxy = galaxy.shape[0]
	while k2 < sizegalaxy:
		galaxyid = int(galaxy[k2][1])
		if galaxy[k2][0] == groupid and galaxy[k2][3] == 1:
			massBCG = galaxy[k2][18]
			break
		k2 += 1
	return k2,massBCG

def findcen(group,galaxy,info,d):
	
	sdss_group_center = []
	sizegroup = group.shape[0]
	sizegalaxy = galaxy.shape[0]
	
	pcr = 2.775 * 10**11 

	nBCG = 0
	n500 = 0
	nonmass = 0
	for j1 in range(sizegroup):

		groupid = group[j1][0]
		groupmass = group[j1][6] #halo mass based on group_L

		if groupid == 7:
			fig = plt.figure()
			ax = fig.add_subplot(111, projection='3d')

		r200 = pow(10**groupmass * 3 / math.pi / 800/pcr, 1 / 3) * 1000 / 0.678#kpc 	
		lrxg,lryg,lrzg = position_xyz(group[j1][1],group[j1][2],group[j1][3])#luminosity weghted group centre 
		
		kk,massBCG = BCGloc(groupid,galaxy,info,d)
		rxg,ryg,rzg = position_xyz(galaxy[kk][4],galaxy[kk][5],galaxy[kk][6])

		rcen = math.sqrt(rxg**2 + ryg**2 + rzg**2)

		ju_BCG = -1
		ju_500 = -1
		massmax = massBCG
		if massBCG <= 5*10**10:
			# print(massBCG)
			nonmass += 1
		mcg = 0

		k1 = 0
		while k1 < sizegalaxy:
			if galaxy[k1][0] == groupid and galaxy[k1][18] > 5* 10**10:
				galaxyid = int(galaxy[k1][1])

				rx,ry,rz = position_xyz(galaxy[k1][4],galaxy[k1][5],galaxy[k1][6])
				rgal = math.sqrt(rx**2 + ry**2 + rz**2)
				r1 = abs(rgal - rcen)
				r0 = math.sqrt((rxg - rx)**2 + (ryg - ry)**2 + (rzg - rz)**2)
				r0r = math.sqrt(r0**2 - r1**2)

				if groupid == 7:
					ax.scatter3D(rx-rxg,ry-ryg,rz-rzg,color = 'b')

				# #Massive galaxy	within 500kpc		
				if r0r < 150 and galaxy[k1][18] > massmax:
					massmax = galaxy[k1][18]
					xc = galaxy[k1][4]
					yc = galaxy[k1][5]
					zc = galaxy[k1][6]
					xxc = rx
					yyc = ry
					zzc = rz
					ju_500 = -2
			k1 += 1	
		
		if massBCG > 5 * 10**10:
			sdss_group_center.append([groupid,galaxy[kk][4],galaxy[kk][5],galaxy[kk][6],r200])
			nBCG += 1
		else:
			if ju_500 == -2:
				sdss_group_center.append([groupid,xc,yc,zc,r200])
				print(massBCG,groupid)
				n500 += 1

		# u = np.linspace(0,2*np.pi,1000)
		# v = np.linspace(0,2*np.pi,1000)
		# x=50*np.outer(np.cos(u),np.sin(v))
		# y=50*np.outer(np.sin(u),np.sin(v))
		# z=50*np.outer(np.ones(np.size(u)),np.cos(v))
		# ax.plot_surface(x,y,z,color = 'b',alpha = 0.1)
		if groupid == 7:
			if massBCG > 5 * 10**10:
				ax.scatter3D(0,0,0,color = 'red',label = 'BCG',s =30,marker = '^')
			else:
				ax.scatter3D(0,0,0,color = 'green',label = 'BCG',s = 30,marker = '^')
			ax.scatter3D(lrxg-rxg,lryg - ryg,lrzg - rzg,s = 30,marker = 'd',color = 'black',label = 'luminosity weighted C')
			ax.scatter3D(xxc-rxg,yyc - ryg,zzc - rzg,s = 35,marker = 'o',color = 'orange',label = 'New centre')
			ax.legend()
			plt.savefig('/home/qyli/Desktop/test/300cluster/centerplot/%s.png' %int(groupid))
	
	# plt.show()

	print(n500,nBCG,nonmass)
	sdss_group_center = np.array(sdss_group_center)#rx,ry,rz is kpc
	return sdss_group_center

def main(d,savedata):
	filegalaxy = '/home/qyli/Desktop/300data/sdss7/Planck/sdss7am_Planck.txt' 
	filegroup = '/home/qyli/Desktop/300data/sdss7/Planck/group_Planck.txt'
	fileinfo = '/home/qyli/Desktop/SDSS7/SDSS_DR7/SDSS7_INFO'
	group = np.loadtxt(filegroup)
	galaxy = np.loadtxt(filegalaxy)
	info = np.loadtxt(fileinfo)

	sdss_group_center = findcen(group,galaxy,info,d)
	if savedata == True:
			np.savetxt('/home/qyli/Desktop/300data/sdss7/Planck/sdss_group_centerPlanckBCG.txt',sdss_group_center)
main(d = 500,savedata = True)

