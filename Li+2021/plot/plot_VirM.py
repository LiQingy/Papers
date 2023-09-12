import numpy as np 
import h5py
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import oPDFplus as opd


def onlycontour(model, pptype, nubtype, rcut):

	fileCM = '/Users/liqy/Documents/data/oPDFnew/cluster/G3X-CM-masscomplete-ALL.txt'
	filefit = '/Users/liqy/Desktop/Trueini/%s_%s_fitmc_rcin%s_%s_TMP.txt' %(model,pptype,rcut,nubtype)
	dataCM = np.loadtxt(fileCM)
	datafit = np.loadtxt(filefit)

	mmtrue = datafit[:,0] / dataCM[:,2] * 1e15
	cctrue = datafit[:,1] / dataCM[:,6]
	loc = opd.sigma3sel(np.log10(mmtrue),np.log10(cctrue))
	print(loc.shape)
	
	print('The cluster halo number with oPDF is %s' %loc.shape[0])
	m_m0 = np.log10(mmtrue[loc])
	c_c0 = np.log10(cctrue[loc])
	xs,ys,zs,siglevel = opd.sigma2dis(m_m0,c_c0,prisig = True)
	plt.contour(xs,ys,zs,levels = siglevel,linestyles = ['-'],colors = 'grey')
	plt.plot(10**np.mean(m_m0),10**np.mean(c_c0), marker = 'o', color = 'grey', zorder = 14,markersize= 3)

def mc_all(loc,dataCM,datafit,col):

	mx = datafit[loc,0]  / dataCM[loc,2]
	cy = datafit[loc,1] / dataCM[loc,6]

	m_m0 = np.log10(mx)
	c_c0 = np.log10(cy)
	xs,ys,zs,siglevel = opd.sigma2dis(m_m0,c_c0,prisig = True)
	
	print('The cluster halo number with SJE is %s' %mx.shape[0])
	print("The mean M/Mtrue and C/Ctrue ",np.mean(mx),np.mean(cy))

	plt.contour(xs,ys,zs,levels = siglevel,linestyles = ['-'],colors = 'green')
	p2 = plt.plot(10**np.mean(m_m0),10**np.mean(c_c0),'go',zorder = 14,markersize= 3)
	p1 = plt.scatter(mx,cy,c = col, s = 5, marker = 'o')

	print("The mean log M/Mtrue and log C/Ctrue ",np.mean(m_m0),np.mean(c_c0))

def main():
	fileCM = '/Users/liqy/Documents/data/oPDFnew/cluster/G3X-CM-masscomplete-ALL.txt'
	filefit = '/Users/liqy/Documents/data/oPDF1/VirM/VirM_Chi2fit_MCDM_rcin200_TMP_n1000_bin20_boot200.txt'
	# filefit = '/Users/liqy/Desktop/VirMArt_Chi2fit_MCDM_rcin200_TMP_n1000_bin20_boot200.txt'

	dataCM = np.loadtxt(fileCM)
	datafit = np.loadtxt(filefit)
	
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	plt.tick_params(top = 'on', right = 'on', which='both')

	#--------------------------------------------------------------------------------------	
	mmtrue = datafit[:,0]  / dataCM[:,2]
	cctrue = datafit[:,1] / dataCM[:,6]
	loc = opd.sigma3sel(np.log10(mmtrue),np.log10(cctrue))

	print("The outliers 3 sigma cut loc",np.delete(range(0,324),loc))

	mc_all(loc,dataCM,datafit,col = 'magenta')

def Main():
	plt.figure(figsize = (3.5,3.4))
	main()
	onlycontour(model = 'GXsub', pptype = 'DM', nubtype = '100000', rcut = 200)

	gas_patch = Patch(edgecolor='brown', label='VT',facecolor='white')
	DM_patch = Patch(edgecolor='grey', label='oPDF',facecolor='white')
	plt.legend(handles=[gas_patch,DM_patch], fontsize = 'small', loc = 1)

	plt.tick_params(labelsize = 10)
	plt.xlim(1e-1,10)
	plt.ylim(1e-1,10)
	plt.xlabel(r'$M/M_{\rm true}$',fontsize = 12)
	plt.ylabel(r'$c/c_{\rm true}$',fontsize = 12)
	plt.loglog()

	plt.text(0.23,0.25,'DM', fontdict = {'fontsize':9,'fontweight': 400})
	plt.text(0.14,0.180,'200 kpc/$h$ cut', fontdict = {'fontsize':9,'fontweight': 400})
	plt.plot(plt.xlim(1e-1,10),[1,1], 'k--', [1,1], plt.ylim(1e-1,10),  'k--', lw = 1.)

	bwith = 1.2 #边框宽度设置为2
	ax = plt.gca()#获取边框
	ax.spines['bottom'].set_linewidth(bwith)
	ax.spines['left'].set_linewidth(bwith)
	ax.spines['top'].set_linewidth(bwith)
	ax.spines['right'].set_linewidth(bwith)

	plt.tight_layout()
	plt.savefig('/Users/liqy/Desktop/figs_Trueini/overall_VirM.pdf')
	plt.show()
Main()
