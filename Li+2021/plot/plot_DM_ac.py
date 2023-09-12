import numpy as np 
import h5py
# from oPDF.oPDF import *
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import oPDFplus as opd
import itertools

def mc_all(loc,dataCM,datafit,col,lab,linet):
	rvir = np.mean(dataCM[loc,3])

	mx = datafit[loc,0] / dataCM[loc,2] * 1e15
	cy = datafit[loc,1] / dataCM[loc,6]

	m_m0 = np.log10(mx)
	c_c0 = np.log10(cy)
	xs,ys,zs,siglevel = opd.sigma2dis(m_m0,c_c0,prisig = False)

	plt.contour(xs,ys,zs,levels = siglevel,linestyles = [linet],colors = col)

	p1 = plt.scatter(mx,cy,c = col,s = 3,label = lab)

def main():
	fileCM = '/Users/liqy/Documents/data/oPDFnew/cluster/G3X-CM-masscomplete-ALL.txt'
	filefit = '/Users/liqy/Desktop/Trueini/GXsub_DM_fitmc_rcin200_100000_TMP.txt'
	fileac = '/Users/liqy/Documents/data/oPDFnew/cluster/GX_cluphy_c.txt'

	dataCM = np.loadtxt(fileCM)
	datafit = np.loadtxt(filefit)
	dataac = np.loadtxt(fileac)
	
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	plt.tick_params(top = 'on', right = 'on', which='both')

	#--------------------------------------------------------------------------------------	
	plt.text(0.23,0.25,'DM', fontdict = {'fontsize':9,'fontweight': 400})
	# plt.text(0.135,0.018,'200 kpc/h cut', fontdict = {'fontsize':10,'fontweight': 300})
	# plt.text(0.14,0.04,'$10^5$ particles', fontdict = {'fontsize':10,'fontweight': 500})
	# loc = np.array(dataCM[:,0] - 1, dtype = 'i4')

	mmtrue = datafit[:,0] / dataCM[:,2] * 1e15
	cctrue = datafit[:,1] / dataCM[:,6]
	Trueliers = opd.sigma3sel(np.log10(mmtrue),np.log10(cctrue))
	print(Trueliers.shape)
	dataCM = dataCM[Trueliers]
	datafit = datafit[Trueliers]
	dataac = dataac[Trueliers]

	nac = 0.7
	loc = np.where(dataac < nac)[0]
	mc_all(loc,dataCM,datafit,col = 'cyan',lab = 'c/a < 0.7',linet = '-')
	loc = np.where(dataac > nac)[0]
	mc_all(loc,dataCM,datafit,col = 'blue',lab = 'c/a > 0.7',linet = '-')

	plt.tick_params(labelsize = 10)
	plt.xlim(1e-1,10)
	plt.xlabel(r'$M/M_{\rm true}$',fontsize = 12)
	plt.ylabel(r'$c/c_{\rm true}$',fontsize = 12)
	plt.loglog()
	plt.legend(fontsize = 'small')

	bwith = 1.2 #边框宽度设置为2
	ax = plt.gca()#获取边框
	ax.spines['bottom'].set_linewidth(bwith)
	ax.spines['left'].set_linewidth(bwith)
	ax.spines['top'].set_linewidth(bwith)
	ax.spines['right'].set_linewidth(bwith)
	
def Main():
	plt.figure(figsize = (3.5,3.4))
	main()

	plt.text(0.23,0.25,'DM', fontdict = {'fontsize':9,'fontweight': 400})
	plt.text(0.14,0.180,'200 kpc/$h$ cut', fontdict = {'fontsize':9,'fontweight': 400})
	plt.plot(plt.xlim(1e-1,10),[1,1], 'k--', [1,1], plt.ylim(1e-1,10),  'k--', lw = 1.2)
	plt.tight_layout()
	plt.savefig('/Users/liqy/Desktop/figs_Trueini/axis.pdf')
	plt.show()
Main()


	


