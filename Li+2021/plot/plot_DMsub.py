import numpy as np 
import h5py
# from oPDF.oPDF import *
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import oPDFplus as opd
import itertools

def mc_all(loc,dataCM,datafit,col,mak):
	rvir = np.mean(dataCM[loc,3])

	mx = datafit[loc,0] / dataCM[loc,2] * 1e15
	cy = datafit[loc,1] / dataCM[loc,6]

	m_m0 = np.log10(mx)
	c_c0 = np.log10(cy)
	xs,ys,zs,siglevel = opd.sigma2dis(m_m0,c_c0,prisig = True)

	plt.contour(xs,ys,zs,levels = siglevel,linestyles = ['-'],colors = col)
	p1 = plt.scatter(mx,cy,c = col, s = 3, alpha = 0.55, marker = mak)
	plt.plot(10**np.mean(m_m0),10**np.mean(c_c0),color = 'green',marker = mak,zorder = 10,markersize = 3.5)
	print("The mean M/M0 and c/c0 in log space is",np.mean(m_m0),np.mean(c_c0))

def main(model,pptype,nubtype,rcut,col,mak):
	fileCM = '/Users/liqy/Documents/data/oPDFnew/cluster/G3X-CM-masscomplete-ALL.txt'
	filefit = '/Users/liqy/Documents/oPDF1/data/%s_%s_fitmc_rcin%s_%s_TMP.txt' %(model,pptype,rcut,nubtype)

	dataCM = np.loadtxt(fileCM)
	datafit = np.loadtxt(filefit)
	
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	plt.tick_params(top = 'on', right = 'on', which='both')

	#--------------------------------------------------------------------------------------	

	mmtrue = datafit[:,0] / dataCM[:,2] * 1e15
	cctrue = datafit[:,1] / dataCM[:,6]

	loc = opd.sigma3sel(np.log10(mmtrue),np.log10(cctrue))
	print(loc.shape)
	mc_all(loc,dataCM,datafit,col,mak)


	
def Main():
	plt.figure(figsize = (3.5,3.4))

	main(model = 'GXsub', pptype = 'DM', nubtype = '100000', rcut = '200',col = 'grey',mak = 'o')
	main(model = 'GX', pptype = 'DM', nubtype = '100000', rcut = '200',col = 'red',mak = 's')

	plt.text(0.23,0.25,'DM', fontdict = {'fontsize':9,'fontweight': 400})
	plt.text(0.14,0.180,'$10^5$ particles', fontdict = {'fontsize':9,'fontweight': 400})
	plt.text(0.135,0.13,r'$200\ \mathrm{kpc}/h\ \rm cut$', fontdict = {'fontsize':9,'fontweight': 400})
	
	clu_patch = Patch(edgecolor='red', label='All particles',facecolor='white')
	clu0_patch = Patch(edgecolor='grey',label='Excluding subhalo particles',facecolor='white')
	plt.legend(handles=[clu_patch,clu0_patch],loc = 1,fontsize = 'x-small')
	plt.plot(plt.xlim(1e-1,10),[1,1], 'k--', [1,1], plt.ylim(1e-1,10),  'k--', lw = 1)

	plt.tick_params(labelsize = 10)
	plt.xlim(1e-1,10)
	plt.ylim(1e-1,10)
	plt.xlabel(r'$M/M_{\rm true}$',fontsize = 12)
	plt.ylabel(r'$c/c_{\rm true}$',fontsize = 12)
	plt.loglog()
	
	bwith = 1.2 #边框宽度设置为2
	ax = plt.gca()#获取边框
	ax.spines['bottom'].set_linewidth(bwith)
	ax.spines['left'].set_linewidth(bwith)
	ax.spines['top'].set_linewidth(bwith)
	ax.spines['right'].set_linewidth(bwith)

	plt.tight_layout()
	# plt.savefig('/Users/liqy/Documents/oPDF1/fig/DMsub.pdf')
	plt.show()
Main()

