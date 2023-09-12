import numpy as np 
from readsnapsgl import *
from astropy import constants as const
from astropy.cosmology import Planck15 as cosmo

def BCGphy(filename,m500,xc,yc,zc,rv):
	head = readsnapsgl(filename,block = 'HEAD',rhb = False)
	temp = readsnapsgl(filename,block = 'TEMP', ptype = 0,rhb = False)
	pos = readsnapsgl(filename,block = 'POS ', ptype = 0,rhb = False)
	mass = readsnapsgl(filename,block = 'MASS', ptype = 0,rhb = False)

	rr = np.log10(np.sqrt((pos[:,0] - xc)**2 + (pos[:,1] - yc)**2 + (pos[:,2] - zc)**2))
	nrv = np.log10(rv) / 20

	densityone0 = np.zeros(20)
	densityone1 = np.zeros(20)
	densityone2 = np.zeros(20)
	gasxone0 = np.zeros(20)
	gasxone1 = np.zeros(20)
	gasxone2 = np.zeros(20)

	for i in range(20):
		r0 = nrv * i
		r1 = nrv * (i+1)
		loc0 = np.where((rr < r1) & (rr > r0) & (temp > 1e7))[0] #hot
		loc1 = np.where((rr < r1) & (rr > r0) & (temp < 1e7) & (temp > 1e5))[0] #warm
		loc2 = np.where((rr < r1) & (rr > r0) & (temp < 1e5))[0] #cold

		densityone0[i] = np.sum(mass[loc0]) * 1e10 / (4/3*np.pi*((10**r1)**3 - (10**r0)**3))
		densityone1[i] = np.sum(mass[loc1]) * 1e10 / (4/3*np.pi*((10**r1)**3 - (10**r0)**3))  
		densityone2[i] = np.sum(mass[loc2]) * 1e10 / (4/3*np.pi*((10**r1)**3 - (10**r0)**3))   
		gasxone0[i] = 10**(np.sum(mass[loc0] * rr[loc0]) / np.sum(mass[loc0])) / rv
		gasxone1[i] = 10**(np.sum(mass[loc1] * rr[loc1]) / np.sum(mass[loc1])) / rv
		gasxone2[i] = 10**(np.sum(mass[loc2] * rr[loc2]) / np.sum(mass[loc2])) / rv

	return densityone0,densityone1,densityone2,gasxone0,gasxone1,gasxone2

def main(model):
	if model == 'Gadget-X':
		speciname = '/home/nifty2014/TheThreeHundred/catalogues/AHF/complete_sample/G3X_Mass_snap_128-center-cluster.txt'
		lab = 'GX'
	else:
		speciname = '/home/nifty2014/TheThreeHundred/catalogues/AHF/complete_sample/Music_Mass_snap_017-center-cluster.txt'
		lab = 'GM'

	datacen = np.loadtxt(speciname)
	size = datacen.shape[0]

	densityhot = np.zeros(shape = (324,20))
	densitywarm = np.zeros(shape = (324,20))
	densitycold = np.zeros(shape = (324,20))
	gasxhot = np.zeros(shape = (324,20))
	gasxwarm = np.zeros(shape = (324,20))
	gasxcold = np.zeros(shape = (324,20))

	ii = 0
	nnan = 0
	while  ii < size:
		
		region = int(datacen[ii][0])
		xc = datacen[ii][3]
		yc = datacen[ii][4]
		zc = datacen[ii][5]
		rv = datacen[ii][10]
		m500 = datacen[ii][9] *1e10 / 0.678

		if model == 'Gadget-X':
			filename = '/home/nifty2014/TheThreeHundred/simulation/GadgetX/NewMDCLUSTER_0%03d/snap_128' %region
		else:
			filename = '/home/nifty2014/TheThreeHundred/simulation/GadgetMUSIC/The300_MUSIC/NewMDCLUSTER_0%03d/snap_017' %region
		
		densityone0,densityone1,densityone2,gasxone0,gasxone1,gasxone2 = BCGphy(filename,m500,xc,yc,zc,rv)
		densityhot[ii] = densityone0
		densitywarm[ii] = densityone1
		densitycold[ii] = densityone2
		gasxhot[ii] = gasxone0
		gasxwarm[ii] = gasxone1
		gasxcold[ii] = gasxone2
		
		ii += 1
		print(ii)

	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/%s_snaplog_gasxhot.txt' %lab, gasxhot)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/%s_snaplog_gasxwarm.txt' %lab, gasxwarm)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/%s_snaplog_gasxcold.txt' %lab, gasxcold)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/%s_snaplog_densityhot.txt' %lab, densityhot)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/%s_snaplog_densitywarm.txt' %lab, densitywarm)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/%s_snaplog_densitycold.txt' %lab, densitycold)

main('Gadget-X')
main('Gadget-MUSIC')


