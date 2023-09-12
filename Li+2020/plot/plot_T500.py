import numpy as np 
import matplotlib.pyplot as plt

def plot(model,col):
	fileslT500 = '/Users/liqy/Documents/data/300data/simulation/%s_snaplog_slT500.txt' %model
	if model == 'GX':
		filecen = '/Users/liqy/Documents/data/oPDF/cluster/G3X_Mass_snap_128-center-cluster.txt' 
		filecste = '/Users/liqy/Documents/data/300data/center/DS-G3X_Mass_snap_128-center-cluster.txt'
		lab = 'Gadget-X'
		mak = 'o'
	elif model =='GM':
		filecen = '/Users/liqy/Documents/data/oPDF/cluster/Music_Mass_snap_017-center-cluster.txt'
		filecste = '/Users/liqy/Documents/data/300data/center/DS-Music_Mass_snap_017-center-cluster.txt' 
		lab = 'Gadget-MUSIC'
		mak = '^'

	slT500 = np.loadtxt(fileslT500)
	datacen = np.loadtxt(filecen)
	datacste = np.loadtxt(filecste)

	Gt500 = 8.85 * (datacen[:, 9] * 1e10 / 0.7 / 1e15) **(2/3) * (0.6125/0.6)
	print(np.min(slT500),np.max(slT500),np.median(slT500))

	#relaxed clusters
	loc1 = np.where(datacste[:,2] == 1)[0]
	plt.scatter( datacen[loc1, 9] * 1e10/ 0.678 / 1e15, Gt500[loc1] / slT500[loc1], 
		marker = mak,color = col,label = '%s (relaxed)' %lab, s= 25)

	#un-relaxed clusters
	loc2 = np.where(datacste[:,2] == 0)[0]
	plt.scatter( datacen[loc2, 9] * 1e10/ 0.678 / 1e15, Gt500[loc2] / slT500[loc2],  
		marker = mak,facecolor = 'none',edgecolor = col,label = '%s (un-relaxed)' %lab, s= 25)

	# T500all = np.hstack((slT500.reshape(324,1),Gt500.reshape(324,1)))
	# np.savetxt('/Users/liqy/Desktop/%s_T500unit.txt' %model, T500all)
	# print(T500all.shape,slT500[:5])

def main():
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	plt.tick_params(top = 'on', right = 'on',which='both', labelsize = 11)

	plot(model = 'GX', col = 'red')
	plot(model = 'GM', col = 'blue')

	plt.xlabel(r'$M_{500}[10^{15} M_\odot]$', fontsize = 14)
	plt.ylabel(r'$T_{500,G+19}/T_{500,sl}$', fontsize = 14)
	plt.axhline(1, color = 'k', linestyle = '--')

	plt.legend(fontsize = 'small')
	plt.savefig('/Users/liqy/Desktop/T500.pdf')
	plt.show()
main()

