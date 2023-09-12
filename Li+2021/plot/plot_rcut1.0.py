import numpy as np 
import h5py
# from oPDF.oPDF import *
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import oPDFplus as opd
import itertools

def mc_all(loc,dataCM,datafit,col,axn):

	mx = datafit[loc,0] / dataCM[loc,2] * 1e15
	cy = datafit[loc,1] / dataCM[loc,6]

	m_m0 = np.log10(mx)
	c_c0 = np.log10(cy)
	print("the mean log value are",np.mean(m_m0),np.mean(c_c0))
	xs,ys,zs,siglevel = opd.sigma2dis(m_m0,c_c0,prisig = False)

	axn.contour(xs,ys,zs,levels = siglevel,linestyles = ['-'],colors = col, linewidths = 1.5)
	axn.scatter(mx,cy,c = col, s = 2, alpha = 0.5)

	bwith = 1.2 #边框宽度设置为2
	axn.spines['bottom'].set_linewidth(bwith)
	axn.spines['left'].set_linewidth(bwith)
	axn.spines['top'].set_linewidth(bwith)
	axn.spines['right'].set_linewidth(bwith)

def main(model,pptype,nubtype,rcut,col,wellfit,axn):
	fileCM = '/Users/liqy/Documents/data/oPDFnew/cluster/G3X-CM-masscomplete-ALL.txt'
	filefit = '/Users/liqy/Desktop/Trueini/%s_%s_fitmc_%s_%s_TMP.txt' %(model,pptype,rcut,nubtype)

	dataCM = np.loadtxt(fileCM)
	datafit = np.loadtxt(filefit)
	#-----------------------------------------------------------------------------

	mmtrue = datafit[:,0] / dataCM[:,2] * 1e15
	cctrue = datafit[:,1] / dataCM[:,6]

	if pptype == 'ste':
		mmtrue[46] = 0.0001
		cctrue[46] = 0.0001

	loc = opd.sigma3sel(np.log10(mmtrue),np.log10(cctrue))
	print(loc.shape)

	print(pptype,rcut)
	mc_all(loc,dataCM,datafit,col,axn)
	
def Main(mod, pp, npp, rr, ax):
	if rr == 'inner':
		rcutn = np.array(['rcin0.01','rcin0.1'])
	elif rr == 'outer':
		# rcutn = np.array(['rcout1.2','rcout0.8','rcout0.7'])
		rcutn = np.array(['rcin200','rcout0.6'])

	#--------------------------------------------------------------------------------------
	main(model = mod, pptype = pp, nubtype = npp, rcut = rcutn[1], col = 'r',wellfit = True,axn = ax)
	main(model = mod, pptype = pp, nubtype = npp, rcut = rcutn[0], col = 'k',wellfit = True,axn = ax)
	#--------------------------------------------------------------------------------------

	if pp == 'star':
		ax.text(0.15,0.25,'Halo stars', fontdict = {'fontsize':10,'fontweight': 300})
	elif pp == 'ste':
		ax.text(0.15,0.25,'Satellites', fontdict = {'fontsize':10,'fontweight': 300})
	else:
		ax.text(0.22,0.25,'DM', fontdict = {'fontsize':10,'fontweight': 300})

	# if npp == '10000':
	# 	ax.text(0.14,0.215,'$10^4$ particles', fontdict = {'fontsize':10,'fontweight': 300})
	# elif npp == '100000':
	# 	ax.text(0.14,0.215,'$10^5$ particles', fontdict = {'fontsize':10,'fontweight': 500})	

	if rr == 'inner':
		ax.text(0.11,0.15,'Inner radius cut', fontdict = {'fontsize':10,'fontweight': 500})
	else:

		ax.text(0.11,0.15,'Outer radius cut', fontdict = {'fontsize':10,'fontweight': 500})
	

def plotall():
	fig = plt.figure(figsize=(8.2, 5))
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	plt.tick_params(top = 'on', right = 'on', which='both')

	ax1 = plt.subplot2grid((6,9),(0,0),colspan = 3,rowspan = 3)
	ax2 = plt.subplot2grid((6,9),(3,0),colspan = 3,rowspan = 3)
	ax3 = plt.subplot2grid((6,9),(0,3),colspan = 3,rowspan = 3)
	ax4 = plt.subplot2grid((6,9),(3,3),colspan = 3,rowspan = 3)
	ax5 = plt.subplot2grid((6,9),(0,6),colspan = 3,rowspan = 3)
	ax6 = plt.subplot2grid((6,9),(3,6),colspan = 3,rowspan = 3)
	#-----------------------------------------------------------------------------
	Main(mod = 'GXsub', pp = 'DM', npp = '10000', rr = 'inner', ax = ax1)
	Main(mod = 'GXsub', pp = 'DM', npp = '10000', rr = 'outer', ax = ax2)
	Main(mod = 'GXsub', pp = 'star', npp = '8000', rr = 'inner', ax = ax3)
	Main(mod = 'GXsub', pp = 'star', npp = '8000', rr = 'outer', ax = ax4)
	Main(mod = 'GX', pp = 'ste', npp = '50', rr = 'inner', ax = ax5)
	Main(mod = 'GX', pp = 'ste', npp = '50', rr = 'outer', ax = ax6)
	#-----------------------------------------------------------------------------
	subplots = np.array([ax1,ax2,ax3,ax4,ax5,ax6])
	for subid in subplots:
		subid.tick_params(labelsize = 11)
		subid.plot(plt.xlim(1e-1,10),[1,1], 'k--', [1,1], plt.ylim(1e-1,10),  'k--', lw = 1)
		subid.set_xlim(1e-1,10)
		subid.set_ylim(1e-1,10)
		subid.loglog()
		subid.tick_params(top = 'on', right = 'on', which = 'both', direction = 'in')

		if subid == ax1 or subid ==ax2:
			subid.set_ylabel(r'$c/c_{\rm true}$',fontsize = 12)
		else:
			subid.axes.yaxis.set_ticklabels([])

	ax1.axes.xaxis.set_ticklabels([])
	ax3.axes.xaxis.set_ticklabels([])
	ax5.axes.xaxis.set_ticklabels([])
	ax1.axes.set_yticks([1,10])
	ax2.axes.set_xticks([1e-1,1,10])
	ax4.axes.set_xticks([1,10])
	ax6.axes.set_xticks([1,10])
	ax2.set_xlabel(r'$M/M_{\rm true}$',fontsize = 12)
	ax4.set_xlabel(r'$M/M_{\rm true}$',fontsize = 12)
	ax6.set_xlabel(r'$M/M_{\rm true}$',fontsize = 12)


	rcut0_patch = Patch(edgecolor='k', label=r'$0.01r_{200}$',facecolor='white',linestyle = '-')
	rcut1_patch = Patch(edgecolor='r', label=r'$0.10r_{200}$',facecolor='white',linestyle = '-')
	ax1.legend(handles=[rcut0_patch,rcut1_patch],loc = 1, fontsize = 'small')

	rcut0_patch = Patch(edgecolor='k', label=r'$1.00r_{200}$',facecolor='white',linestyle = '-')
	rcut1_patch = Patch(edgecolor='r', label=r'$0.60r_{200}$',facecolor='white',linestyle = '-')
	ax2.legend(handles=[rcut0_patch,rcut1_patch],loc = 1, fontsize = 'small')

	plt.tight_layout()
	plt.subplots_adjust(wspace = 0.1, hspace = 0.1)

	plt.savefig('/Users/liqy/Desktop/figs_Trueini/rcut.pdf')
	plt.show()

plotall()





