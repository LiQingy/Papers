import numpy as np 
from readsnapsgl import *
from astropy import constants as const
from astropy.cosmology import Planck15 as cosmo

def BCGphy(filename,m500,xc,yc,zc,rv):
	head = readsnapsgl(filename,block = 'HEAD',rhb = False)
	temp0 = readsnapsgl(filename,block = 'TEMP', ptype = 0,rhb = False)
	temp = temp0 * const.k_B.value / const.e.value * 1e-3
	metal = readsnapsgl(filename, block = 'Z   ', ptype = 0,rhb = False)
	pos = readsnapsgl(filename,block = 'POS ', ptype = 0,rhb = False)
	mass = readsnapsgl(filename,block = 'MASS', ptype = 0,rhb = False)
	rho = readsnapsgl(filename, block = 'RHO ', ptype = 0,rhb = False)

	rr = np.log10(np.sqrt((pos[:,0] - xc)**2 + (pos[:,1] - yc)**2 + (pos[:,2] - zc)**2))
	nrv = np.log10(rv) / 20

	metalone = np.zeros(20)
	densityone = np.zeros(20)
	tempone = np.zeros(20)
	gasxone = np.zeros(20)
	gasxden = np.zeros(20)
	for i in range(20):
		r0 = nrv * i
		r1 = nrv * (i+1)
		loc0 = np.where((rr< r1)& (rr > r0))[0]
		loc = np.where((rr < r1) & (rr > r0) & (temp > 0.3) & (rho < 0.000627091))[0]
		metalone[i] = np.sum(mass[loc] *rho[loc]*np.sqrt(temp[loc]) * metal[loc]) / np.sum(mass[loc] *rho[loc]*np.sqrt(temp[loc]))
		densityone[i] = np.sum(mass[loc0]) * 1e10 / (4/3*np.pi*((10**r1)**3 - (10**r0)**3))  
		tempone[i] = np.sum(mass[loc] *rho[loc]*np.sqrt(temp[loc]) * temp[loc]) / np.sum(mass[loc] *rho[loc]*np.sqrt(temp[loc]))
		
		gasxone[i] = 10**(np.sum(mass[loc] *rho[loc]*np.sqrt(temp[loc]) * rr[loc]) / np.sum(mass[loc] *rho[loc]*np.sqrt(temp[loc]))) / rv
		gasxden[i] = 10**(np.sum(mass[loc0] *rr[loc0]) / np.sum(mass[loc0])) / rv
	
	t500 = 8.85 * (m500 / 1e15) **(2/3) * (0.588/0.6) * (0.307 + 0.693)
	tempone = tempone / t500

	return metalone,densityone,tempone,gasxone,gasxden

def main(model):
	if model == 'Gadget-X':
		speciname = '/home/nifty2014/TheThreeHundred/catalogues/AHF/complete_sample/G3X_Mass_snap_128-center-cluster.txt'
		lab = 'GX'
	else:
		speciname = '/home/nifty2014/TheThreeHundred/catalogues/AHF/complete_sample/Music_Mass_snap_017-center-cluster.txt'
		lab = 'GM'

	datacen = np.loadtxt(speciname)
	size = datacen.shape[0]

	densityall = np.zeros(shape = (324,20))
	metalall = np.zeros(shape = (324,20))
	tempall = np.zeros(shape = (324,20))
	gasxall = np.zeros(shape = (324,20))
	gasxdenall = np.zeros(shape = (324,20))

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
		
		metalone,densityone,tempone,gasxone,gasxden = BCGphy(filename,m500,xc,yc,zc,rv)
		densityall[ii] = densityone
		metalall[ii] = metalone
		tempall[ii] = tempone
		gasxall[ii] = gasxone
		gasxdenall[ii] = gasxden
		
		ii += 1
		print(ii)

	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/%s_snaplog_gasx.txt' %lab, gasxall)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/%s_snaplog_gasdenx.txt' %lab, gasxdenall)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/%s_snaplog_metal.txt' %lab, metalall)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/%s_snaplog_temp.txt' %lab, tempall)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/%s_snaplog_density.txt' %lab, densityall)

main('Gadget-X')
main('Gadget-MUSIC')


