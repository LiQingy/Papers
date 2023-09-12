import numpy as np 
import h5py
# from oPDF.oPDF import *
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import oPDFplus as opd
import itertools

def onlycontour(model, pptype, nubtype, rcut,col,lin):
	fileCM = '/Users/liqy/Documents/data/oPDFnew/cluster/G3X-CM-masscomplete-ALL.txt'
	filefit = '/Users/liqy/Documents/oPDF1/data/%s_%s_fitmc_rcin%s_%s_TMP.txt' %(model,pptype,rcut,nubtype)
	dataCM = np.loadtxt(fileCM)
	datafit = np.loadtxt(filefit)

	mmtrue = datafit[:,0] / dataCM[:,2] * 1e15
	cctrue = datafit[:,1] / dataCM[:,6]
	loc = opd.sigma3sel(np.log10(mmtrue),np.log10(cctrue))
	print(loc.shape)

	m_m0 = np.log10(mmtrue[loc])
	c_c0 = np.log10(cctrue[loc])
	xs,ys,zs,siglevel = opd.sigma2dis(m_m0,c_c0,prisig = True)

	plt.contour(xs,ys,zs,levels = siglevel,linestyles = [lin],colors = col,label = pptype)
	if lin == '-':
		plt.plot(10**np.mean(m_m0),10**np.mean(c_c0),c = col, marker = 'o',zorder = 5,markersize= 4)


def Main():
	plt.figure(figsize = (3.5,3.4),dpi = 500)
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	plt.tick_params(top = 'on', right = 'on', which='both')

	onlycontour(model = 'GXsub', pptype = 'DM', nubtype = 'nsame', rcut = '200',col = 'grey', lin = '-')
	# onlycontour(model = 'GXsub', pptype = 'DM', nubtype = '100000', rcut = '200',col = 'grey', lin = '--')
	onlycontour(model = 'GXsub', pptype = 'star', nubtype = 'nsame', rcut = '200',col = 'orange', lin = '-')
	# onlycontour(model = 'GXsub', pptype = 'star', nubtype = '10000', rcut = '200',col = 'orange', lin = '--')
	onlycontour(model = 'GX', pptype = 'ste', nubtype = 'nsame', rcut = '200', col = 'blue', lin = '-')

	gal_patch = Patch(edgecolor='blue', label='Satellite galaxies',facecolor='white')
	DM_patch = Patch(edgecolor='grey', label='DM',facecolor='white')
	star_patch = Patch(edgecolor='orange', label='Halo stars',facecolor='white')
	plt.legend(handles=[DM_patch,star_patch,gal_patch], fontsize = 'x-small',loc = 1)

	plt.plot(plt.xlim(1e-1,10),[1,1], 'k--', [1,1], plt.ylim(1e-1,10),  'k--', lw = 1)
	plt.tick_params(labelsize = 10)
	plt.loglog()
	plt.xlim(1e-1,10)
	plt.ylim(1e-1,10)

	plt.xlabel(r'$M/M_{\rm true}$',fontsize = 12)
	plt.ylabel(r'$c/c_{\rm true}$',fontsize = 12)

	bwith = 1.2 #边框宽度设置为2
	ax = plt.gca()#获取边框
	ax.spines['bottom'].set_linewidth(bwith)
	ax.spines['left'].set_linewidth(bwith)
	ax.spines['top'].set_linewidth(bwith)
	ax.spines['right'].set_linewidth(bwith)

	# plt.text(0.13,0.22,r'$M_{\bigstar}>10^9\ \rm M_{\odot}$', fontdict = {'fontsize':10,'fontweight': 500})
	# plt.text(0.13,0.12,r'$200\ \mathrm{kpc}/h\ \rm cut$', fontdict = {'fontsize':9,'fontweight': 400})
	# plt.text(0.108,0.17,'Satellite galaxies', fontdict = {'fontsize':9,'fontweight': 400})

	plt.tight_layout()
	plt.savefig('/Users/liqy/Documents/oPDF1/fig/tracers.png' )
	plt.show()
Main()

