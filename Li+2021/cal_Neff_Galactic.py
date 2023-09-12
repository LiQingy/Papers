import numpy as np 
import h5py
# from oPDF.oPDF import *
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.optimize import minimize
import oPDFplus as opd
from iminuit import Minuit

def quantity_size(cov,xcon1,ycon1,nx0):
	sig2m = cov[0,0]
	sig2c = cov[1,1]
	sigmc = cov[0,1]
	n1 = 1/2*(sig2m+sig2c+((sig2m-sig2c)**2+4*sigmc**2)**0.5) #to large eig
	n2 = 1/2*(sig2m+sig2c-((sig2m-sig2c)**2+4*sigmc**2)**0.5) #to samll eig
	vals, vecs = np.linalg.eig(cov)
	invT = np.linalg.inv(vecs)
	a = invT[0,0]
	b = invT[0,1]
	c = invT[1,0]
	d = invT[1,1]

	global No
	No = xcon1.shape[0]

	def minichi2(nx):
		Kchi=2.3
		resi = 0
		xcon11 = xcon1*nx
		ycon11 = ycon1*nx
		for i in range(No):
			k = ycon11[i]/xcon11[i]
			xall = Kchi / ((a+b*k)**2/n2 + (c+d*k)**2/n1)
			
			if xcon1[i] > 0:
				x1 = xall**0.5
			else:
				x1 = -xall**0.5
			resi += (x1-xcon11[i])**2+(k*x1-ycon11[i])**2

		return resi

	fitnx = Minuit(minichi2, nx = nx0, error_nx=0.01,print_level=0,errordef=2.3)
	fitnx.migrad()
	nxbest = fitnx.values['nx']
	return nxbest


def mc_all(locin,datafit,dataChi2,col,nx0):

	m_m0 = np.log10(datafit[locin,0])
	c_c0 = np.log10(datafit[locin,1])
	avx = 10**np.mean(m_m0)
	avy = 10**np.mean(c_c0)

	xs,ys,zs,siglevel = opd.sigma2dis(m_m0,c_c0,prisig = False)
	cen = plt.plot(np.log10(avx),np.log10(avy),'go',zorder = 10,markersize=1)
	cs = plt.contour(np.log10(xs),np.log10(ys),zs,levels = siglevel,linestyles = ['-'],colors = 'grey',zorder=100)
	pc=cs.collections[0].get_paths()[0]
	con=pc.vertices

	#move center to (0,0)
	cov = np.cov(m_m0,c_c0)
	xcon=con[:,0] - np.log10(avx)
	ycon=con[:,1] - np.log10(avy)
	plt.scatter(xcon,ycon,c = 'g', s = 5)

	xcon1=dataChi2[:,0]
	ycon1=dataChi2[:,1]
	# plt.scatter(xcon1,ycon1,c = 'k', s = 5)
	nxbest = quantity_size(cov,xcon1,ycon1,nx0)

	plt.scatter(xcon1*nxbest,ycon1*nxbest,c = 'cyan', s = 5)
	plt.plot(xcon1,ycon1,c = 'red')
	
	return nxbest

def main(model,nx0):
	filefit = '/Users/liqy/Documents/data/oPDF1/ForQiangyang/MRII-HBTSat-20-400kpc-Rbin15.np.npy'
	fileChi2 = '/Users/liqy/Documents/data/oPDF1/ForQiangyang/1sigma-10000particles.20-400kpc.dat'

	datafit = np.load(filefit)
	dataChi2 = np.loadtxt(fileChi2)
	
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	plt.tick_params(top = 'on', right = 'on', which='both')
	
	mmtrue = datafit[:,0]
	cctrue = datafit[:,1]
	loc = np.where((10 > mmtrue) & (mmtrue > 1e-1) & (10 > cctrue) & (cctrue > 1e-1))[0]
	col = 'grey'
	ntracer = 500
	nxbest = mc_all(loc,datafit,dataChi2,col,nx0)
	print("The best nx is ",nxbest)
	print("The Neff is ", ntracer / (nxbest**2/20-1))

	
def Main(model,nx0):
	main(model,nx0)


	plt.legend(loc = 1,fontsize = 'x-small')

	plt.tick_params(labelsize = 10)
	plt.xlabel(r'$M/M_{\rm true}$',fontsize = 12)
	plt.ylabel(r'$c/c_{\rm true}$',fontsize = 12)

	bwith = 1.2 #边框宽度设置为2
	ax = plt.gca()#获取边框
	ax.spines['bottom'].set_linewidth(bwith)
	ax.spines['left'].set_linewidth(bwith)
	ax.spines['top'].set_linewidth(bwith)
	ax.spines['right'].set_linewidth(bwith)
	plt.axhline(0, ls = '--', c = 'k', lw = 1)
	plt.axvline(0, ls = '--', c = 'k', lw =1 )
	plt.xlim(-1,1)
	plt.ylim(-1,1)
	plt.title('%s' %model)
	plt.show()

	

'''
nx0: initial value to the expansion multiple of statistic error 
'''
Main(model = 'MW-Galactic', nx0 = 6)






