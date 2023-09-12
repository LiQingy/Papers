import numpy as np 
import h5py
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import oPDFplus as opd
import itertools

def mc_all(loc,mmtrue,cctrue,state,axn):

	m_m0 = np.log10(mmtrue[loc])
	c_c0 = np.log10(cctrue[loc])
	xs,ys,zs,siglevel = opd.sigma2dis(m_m0,c_c0,prisig = True)

	axn.contour(xs,ys,zs,levels = siglevel,linestyles = '-',colors = state['col'])
	# axn.scatter(10**np.mean(m_m0),10**np.mean(c_c0),color = state['col'], marker = 'o', facecolor = 'none',label = state['lab1'],zorder = 10,s= 40)
	
	axn.plot(plt.xlim(1e-1,10),[1,1], 'k--', [1,1], plt.ylim(1e-1,10),  'k--', lw = 1)
	axn.scatter(mmtrue[loc],cctrue[loc],c = state['col'], marker = state['sym'],s = 3)


#------------------------------------------------

def DS(model,dataCM,datastate,datafit,axn):

	if model == 'oPDF':
		mmtrue = datafit[:,0] / dataCM[:,2] * 1e15
	else:
		mmtrue = datafit[:,0] / dataCM[:,2]
	cctrue = datafit[:,1] / dataCM[:,6]

	unrelax = {'sym':'o','lab':'un-relaxed clusters','col':'cyan','linet':'--','lab1':'average center (un-relaxed)'}
	loc = opd.sigma3sel(np.log10(mmtrue),np.log10(cctrue))
	print(loc.shape)

	loc_unrex = np.where((datastate[loc,2] == 0))[0]
	print(loc_unrex.shape[0])
	mc_all(loc[loc_unrex],mmtrue,cctrue,unrelax,axn)
	
	relax = {'sym':'o','lab':'relaxed clusters','col':'blue','linet':'-','lab1':'average center (relaxed)'}
	loc_rex = np.where((datastate[loc,2] == 1))[0]
	print(loc_rex.shape[0])
	mc_all(loc[loc_rex],mmtrue,cctrue,relax,axn)

#--------------------------
def main(model,axn):

	fileCM = '/Users/liqy/Documents/data/oPDFnew/cluster/G3X-CM-masscomplete-ALL.txt'
	filestate = '/Users/liqy/Documents/data/oPDFnew/cluster/DS-G3X_Mass_snap_128-center-cluster.txt'
	if model == 'oPDF':
		filefit = '/Users/liqy/Documents/oPDF1/data/GXsub_DM_fitmc_rcin200_100000_TMP.txt'
		axn.text(0.15,0.2,'oPDF', fontdict = {'fontsize':9,'fontweight': 400})
	else:
		filefit = '/Users/liqy/Documents/data/oPDF1/JeansE/JeansE_Chi2fit_MCDM_rcin200_TMP.txt'
		axn.text(0.15,0.2,'SJE', fontdict = {'fontsize':9,'fontweight': 400})

	dataCM = np.loadtxt(fileCM)
	datastate = np.loadtxt(filestate)
	datafit = np.loadtxt(filefit)
	DS(model,dataCM,datastate,datafit,axn = axn)
	#---------------------------------------------------------------------------------------------------------
		

	axn.set_xlim(1e-1,10)
	axn.set_ylim(1e-1,10)
	axn.tick_params(labelsize = 10)

	axn.set_xlabel(r'$M/M_{\rm true}$',fontsize = 12)
	axn.set_ylabel(r'$c/c_{\rm true}$', fontsize = 12)
	axn.loglog()

	bwith = 1.2 #边框宽度设置为2
	ax = plt.gca()#获取边框
	axn.spines['bottom'].set_linewidth(bwith)
	axn.spines['left'].set_linewidth(bwith)
	axn.spines['top'].set_linewidth(bwith)
	axn.spines['right'].set_linewidth(bwith)
	axn.tick_params(top = 'on', right = 'on', which = 'both', direction = 'in')

def plotall():
	fig = plt.figure(figsize=(7.2, 3.5), dpi = 500)
	plt.rcParams['figure.constrained_layout.use'] = True

	ax1 = plt.subplot2grid((4,8),(0,0),colspan = 4,rowspan = 4)
	ax2 = plt.subplot2grid((4,8),(0,4),colspan = 4,rowspan = 4)

	main(model = 'oPDF',axn = ax1)
	main(model = 'JSE',axn = ax2)

	relax_patch = Patch(edgecolor='blue', label='relaxed',facecolor='white')
	unrelax_patch1 = Patch(edgecolor='cyan', label='un-relaxed',facecolor='white')
	legend2 = plt.legend(handles=[relax_patch,unrelax_patch1],loc = 1)

	plt.subplots_adjust(wspace = 1.4 ,hspace = 0)
	plt.tight_layout()
	plt.savefig('/Users/liqy/Documents/oPDF1/fig/DS.png')
	plt.show()
plotall()
#for gal
# main(model = 'GX', pptype = 'ste', rcut = '0',nubtype = 'nsame')







