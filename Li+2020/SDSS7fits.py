import numpy as np 
import time
import math
from astropy.io import fits 


def main():
	filephy = '/home/nifty2014/TheThreeHundred/playground/qingyang/300code/sdss_firefly-26.fits'
	dataphy0 = fits.open(filephy)
	dataphy = dataphy0[1]

	size = dataphy.data.shape[0]
	i = 0
	dphy = []
	print('loading is ok')
	while i < size:
		data1 = dataphy.data[i][61:65]
		data2 = dataphy.data[i][143:154]
		data = np.insert(data2,0,data1)
		dphy.append(data)
		i += 1
	dphy = np.array(dphy)
	print('extracting data is ok')
	dataphy0.close()
#--------------------------------------------------------------------------------------------------
	ready = np.sort(dphy.view('f8,' * 15),axis = 0,order = ['f2','f0','f3']).view(np.float)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/300code/sdss_firefly_ready.txt',ready)
	fileid = '/home/nifty2014/TheThreeHundred/playground/qingyang/Center/SDSS7_ID_new.txt'
	dataid = np.loadtxt(fileid)
	sizeid = dataid.shape[0]
	sizephy = ready.shape[0]
	j = 41202
	loc = 0
	galaxy = []

	print('sorting is ok')
	while j < sizeid:
		nid = dataid[j][0]
		mjd = dataid[j][2]
		plate = dataid[j][3]
		fiberid = dataid[j][4]
		k = loc
		while  k < sizephy:
			if mjd == int(ready[k][2]):
				loc = k
				if plate < int(ready[k][0]):
					break
				if plate == int(ready[k][0]):
					loc = k
					if fiberid < int(ready[k][3]):
						break
					if fiberid == int(ready[k][3]):
						loc = k + 1
						dd = np.insert(ready[k],0,nid)
						galaxy.append(dd)
						break
			k += 1
		j += 1

	print('matching is ok')
	galaxy = np.array(galaxy)
	print(galaxy.shape)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/300code/sdss_firefly-26_SDSS7.txt',galaxy,delimiter = ' ',fmt = '%i ' * 5 + '%f ' * 11 )
	print('saving data is ok')
main()