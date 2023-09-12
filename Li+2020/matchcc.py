import numpy as np 

def main(model):
	if model == 'G3X': 
		file180 = '/Users/liqy/Documents/data/300data/center/G3X_M100_M180_snap_128info.txt'
		filecc = '/Users/liqy/Documents/data/300data/center/pseudoentropy_g3x.txt'
		filecen = '/Users/liqy/Documents/data/300data/center/G3X_Mass_snap_128-center-cluster.txt'
	else:
		file180 = '/Users/liqy/Documents/data/300data/center/Music_M100_M180_snap_017info.txt'
		filecc = '/Users/liqy/Documents/data/300data/center/pseudoentropy_music.txt'
		filecen = '/Users/liqy/Documents/data/300data/center/Music_Mass_snap_017-center-cluster.txt'


	data180 = np.loadtxt(file180)
	datacc = np.loadtxt(filecc)
	datacen = np.loadtxt(filecen)

	cc = np.zeros(324)

	for i in range(datacen.shape[0]):
		regionID = datacen[i][0]
		AHID = datacen[i][1]
		loc = np.where((regionID == data180[:,0]) & (AHID == data180[:,1]))[0]
		print(loc)
		if len(loc):
			cc[i] = datacc[loc,3]
		else:
			print('#')
			cc[i] = -1
	np.savetxt('/Users/liqy/Documents/data/300data/center/%s_core.txt' %model,cc)
main('Music')

