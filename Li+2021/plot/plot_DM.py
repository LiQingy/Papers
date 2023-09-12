import numpy as np 
import h5py
# from oPDF.oPDF import *
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import oPDFplus as opd
import myutils as mys


def MW(pptype):
	def cont2(cov,av,rr):
		msigma = np.sqrt(cov[0][0])
		csigma = np.sqrt(cov[1][1])
		pcov = cov[0][1] / msigma / csigma
		mave = av[0]
		cave = av[1]
		def gaussian(xx,yy):
			A = 1 / 2 / np.pi / msigma / csigma / np.sqrt(1 - pcov**2)
			B = (xx - mave)**2 / msigma**2 + (yy - cave)**2 / csigma**2 - 2*pcov*(xx - mave)*(yy - cave) / msigma / csigma
			return A * np.exp(-1 / 2 / (1 - pcov**2) * B)
		xx = np.linspace(-rr,rr,500)
		yy = np.linspace(-rr,rr,500)
		X,Y = np.meshgrid(xx,yy)
		Z = gaussian(X,Y)

		k = 2.3
		k1 = np.exp(-k / 2) / 2 / np.pi / msigma / csigma / np.sqrt(1 - pcov**2)
		sigmalevel = [10**k1]
		return 10**X,10**Y,10**Z,sigmalevel
	
	datadir = '/Users/liqy/Documents/data/oPDF1/WentingFits/'
	MRII=np.loadtxt(datadir+'MRIIfits.dat') #isolated
	MRIIb=np.loadtxt(datadir+'MRIIfitsbinary.dat') #binary
	print(MRIIb.shape, MRII.shape)
	MRall = np.append(MRII,MRIIb[:,:8],axis = 0)
	MRall = MRII
	MR=np.vstack([MRII[:,2:4], MRIIb[:,2:4]]) #combine
	# MR=MRII[:,2:4]
	print(MR.shape)
	sel=(MR[:,0]>0.1)&(MR[:,1]>0.1)&(MR[:,0]<10)&(MR[:,1]<10)
	# loc=np.arange(0,MR.shape[0],1)

	covMR=np.cov(np.log10(MR[sel,:]).T)
	avMR=np.log10(MR[sel,:]).mean(axis=0)
	MW_mc = MR[sel,:]
	print('MW halo mean log M/Mtrue and log C/Ctrue', avMR)
	print('MW halo mean M/Mtrue and C/Ctrue', MR[sel,:].mean(axis=0))
	print('MW halo covariance matrix is %s' %covMR)
	print('sigma M and sigma c is %s, %s' %(covMR[0][0]**0.5, covMR[1][1]**0.5))
	print('correlation coefficience is %s' %(covMR[0][1]/covMR[0][0]**0.5/covMR[1][1]**0.5))
	print('MW halo number is %s' %MW_mc.shape[0])
	# plt.scatter(MW_mc[:,0], MW_mc[:,1],c = 'r', s = 3)
	plt.scatter(MW_mc[:,0], MW_mc[:,1], color = 'k',s = 2.2,zorder = 1, alpha = 0.85, marker = '^',label='Galactic halos')
	plt.plot(10**avMR[0],10**avMR[1], 'g^',zorder = 10,markersize= 4)
	xs,ys,zs,siglevel = cont2(covMR,avMR,rr=0.5)
	plt.contour(xs,ys,zs,levels = siglevel,linestyles = ['-'],colors = 'k', zorder = 9)

def mc_all(loc,dataCM,datafit,col):
	rvir = np.mean(dataCM[loc,3])

	mx = datafit[loc,0] / dataCM[loc,2] * 1e15
	cy = datafit[loc,1] / dataCM[loc,6]

	print('##',np.sqrt(np.sum((mx - 1)**2) / mx.shape[0]))
	print('The cluster halo number is %s' %mx.shape[0])
	print("The cluster fit mean M/Mtrue and C/Ctrue",np.mean(mx),np.mean(cy))

	xs,ys,zs,siglevel = opd.sigma2dis(mx,cy,prisig = True)

	m_m0 = np.log10(mx)
	c_c0 = np.log10(cy)
	xs,ys,zs,siglevel = opd.sigma2dis(m_m0,c_c0,prisig = True)

	plt.contour(xs,ys,zs,levels = siglevel,linestyles = ['-'],colors = 'grey')
	p2 = plt.plot(10**np.mean(m_m0),10**np.mean(c_c0),'go',zorder = 10,markersize= 4)
	
	plt.plot(plt.xlim(1e-1,10),[1,1], 'k--', [1,1], plt.ylim(1e-1,10),  'k--', lw = 1)
	p1 = plt.plot(mx,cy,'.',c = col, ms = 1.8, label='Clusters', alpha = 1, zorder = 5)

	# print("The cluster fit mean M/Mtrue and C/Ctrue",np.mean(m_m0),np.mean(c_c0))

	print("The cluster fit mean log M/Mtrue and log C/Ctrue",np.mean(m_m0),np.mean(c_c0))

def allpoint(dataChi2,dx0): # add all the log-likelihood radio for single cluster and get 1 sigma conficence 
	xl=np.log10(3)-dx0
	xu=np.log10(3)+dx0
	x=np.logspace(xl, xu, 20) 
	y=np.logspace(xl, xu, 20)

	# print(dataChi2)
	# dataChi2 = mys.Chi2Sig(dataChi2, dof=2)
	'''convert chi-square value to significance level, for dof degrees of freedom
		in order to guage the times
	'''

	plt.contour(x,y,dataChi2,levels=[2.3],colors = 'red',linewidths = [0.5],zorder = 100)


def main(model,pptype,nubtype,rcut):
	fileCM = '/Users/liqy/Documents/data/oPDFnew/cluster/G3X-CM-masscomplete-ALL.txt'
	filefit = '/Users/liqy/Documents/oPDF1/data/%s_%s_fitmc_rcin%s_%s_TMP.txt' %(model,pptype,rcut,nubtype)
	fileChi2 = '/Users/liqy/Documents/oPDF1/data/%s_%s_Chi2_rcin%s_%s_TMP.txt' %(model,pptype,rcut,nubtype)

	dataCM = np.loadtxt(fileCM)
	datafit = np.loadtxt(filefit)
	dataChi2 = np.loadtxt(fileChi2)
	
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	plt.tick_params(top = 'on', right = 'on', which='both')

	#--------------------------------------------------------------------------------------	
	plt.text(0.23,0.25,'DM', fontdict = {'fontsize':9,'fontweight': 400})
	plt.text(0.14,0.180,'$10^5$ particles', fontdict = {'fontsize':9,'fontweight': 400})
	plt.text(0.135,0.13,r'$200\ \mathrm{kpc}/h\ \rm cut$', fontdict = {'fontsize':9,'fontweight': 400})
	
	mmtrue = datafit[:,0] / dataCM[:,2] * 1e15
	cctrue = datafit[:,1] / dataCM[:,6]
	loc = opd.sigma3sel(np.log10(mmtrue),np.log10(cctrue))
	print(loc.shape)
	print(np.delete(np.arange(0,324), loc))

	col = 'grey'
	mc_all(loc,dataCM,datafit,col)
	MW(pptype)
	dx0 = 0.02
	allpoint(dataChi2,dx0)

	plt.legend(loc = 1,fontsize = 'x-small')

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
	
def Main(model,pptype,nubtype,rcut):
	plt.figure(figsize = (3.5,3.4), dpi = 500)

	main(model,pptype,nubtype,rcut)

	plt.tight_layout()
	plt.savefig('/Users/liqy/Documents/oPDF1/fig/overall_DM.png')
	plt.show()
Main(model = 'GXsub', pptype = 'DM', nubtype = '100000', rcut = '200')

