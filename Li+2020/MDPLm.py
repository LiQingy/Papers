import numpy as np
import math
import readsnapsgl as readsgl

		

def SAM(n,xc,yc,zc,rv,region):
	filename = '/TheThreeHundred/simulation/MDPL2/NewMDCLUSTER_0%03d/snap_128' %region
	data = readsgl.readsnapsgl(filename, block = 'POS ')
	
	h = rv / n

	massall = np.zeros(n)
	rrall = np.zeros(n)

	size = data.shape[0]	
	for j in range(size):
		rx = data[j][0]  / 0.678
		ry = data[j][1]  / 0.678
		rz = data[j][2]  / 0.678
		rr = math.sqrt((rx - xc)**2 + (ry - yc)**2 + (rz - zc)**2)		
		if rr < rv:
			nn = int(rr / h)
			massall[nn] +=  1.5e9 / 0.678#mass
			rrall[nn] += rr / rv * (1.5e9 / 0.678) 
#---------------------------------------------------------------------------------------------
	ndv = np.zeros(n)
	for i in range(n):
		dv = 4. / 3. * math.pi * (pow((i + 1) * h, 3) - pow(i * h, 3))
		ndv[i] = dv
	
	densityall = massall / ndv
	rrall = rrall / massall

	return densityall,rrall


def main(n,savedata):
	filecen = '/home/nifty2014/TheThreeHundred/playground/weiguang/MDPL2_Mass_snap_128-center-cluster.txt'
	datacen = np.loadtxt(filecen)

	densityall_MDPL = np.zeros(shape = (324,n))
	rrall_MDPL = np.zeros(shape = (324,n))

	region = 1
	while region <= 324:

		rv = datacen[region - 1][6] / 0.678
		xc = datacen[region - 1][3]  / 0.678
		yc = datacen[region - 1][4]  / 0.678
		zc = datacen[region - 1][5]  / 0.678
		
		densityall,rrall = SAM(n,xc,yc,zc,rv,region)
		densityall_MDPL[region - 1] = densityall
		rrall_MDPL[region - 1] = rrall

		print(region)
		region += 1
		

	if savedata == True:
		np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/SAMdata/SAM_MDPL_densityall.txt' , densityall_MDPL)
		np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/SAMdata/SAM_MDPL_rrall.txt' , rrall_MDPL)

main(20,savedata = True)
