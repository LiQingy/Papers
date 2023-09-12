import numpy as np 
import readsnapsgl

def BCGphy(filename,rscale,xc,yc,zc,rv):
	age = readsnapsgl(filename,block = 'AGE ', ptype = 4)
	metal = readsnapsgl(filename, block = 'Z   ', ptype = 4)
	pos = readsnapsgl(filename,block = 'POS ', ptype = 4)
	mass = readsnapsgl(filename,block = 'MASS', ptype = 4)
	print(age.shape,metal.shape,pos.shape,mass.shape)

	rr = np.sqrt((pos[:,0] - xc)**2 + (pos[:,1] - yc)**2 + (pos[:,2] - zc)**2)
	loc = np.where(rr < rscale * rv)
	print(rr.shape)
	ageBCG = np.sum(mass[loc] * age[loc]) / np.sum(mass[loc])
	metalBCG = np.sum(mass[loc] * metal[loc]) / np.sum(mass[loc])

	return ageBCG,metalBCG

def main(model,rscale):
	if model == 'Gadget-X':
		speciname = '/home/nifty2014/TheThreeHundred/catalogues/AHF/complete_sample/G3X_Mass_snap_128-center-cluster.txt'
		lab = 'GX'
	else:
		speciname = '/home/nifty2014/TheThreeHundred/catalogues/AHF/complete_sample/Music_Mass_snap_017-center-cluster.txt'
		lab = 'GM'

	data = np.loadtxt(speciname)
	size = data.shape[0]

	ageall = np.zeros(324)
	metalall = np.zeros(324)

	ii = 0
	while  ii < size:
		
		region = int(data[ii][0])
		xc = data[ii][3]
		yc = data[ii][4]
		zc = data[ii][5]
		rv = data[ii][10]

		if model == 'Gadget-X':
			filename = '/home/nifty2014/TheThreeHundred/simulation/GadgetX/NewMDCLUSTER_0%03d/snap_128' %region
		else:
			filename = '/home/nifty2014/TheThreeHundred/simulation/GadgetMUSIC/The300_MUSIC/NewMDCLUSTER_0%03d/snap_017' %region
		
		ageBCG,metalBCG = BCGphy(filename,rscale,xc,yc,zc,rv)
		ageall[ii] = ageBCG
		metalall[ii] = metalBCG

		np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/GX500_snap_ageBCG.txt', ageall)
		np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/GX500_snap_metalBCG.txt', metalall)
main('Gadget-X',0.01)