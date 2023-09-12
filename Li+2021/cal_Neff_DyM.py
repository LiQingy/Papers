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

	fitnx = Minuit(minichi2, nx = nx0, error_nx=0.1,print_level=0,errordef=1)
	fitnx.migrad()
	fitnx.hesse()
	nxbest = fitnx.values['nx']
	err = fitnx.covariance['nx','nx']

	print("The nx error is ", err**0.5)
	return nxbest


def mc_all(locin,dataCM,datafit,dataChi2,col,nx0,nd):

	mx = datafit[locin,0] / dataCM[locin,2]
	cy = datafit[locin,1] / dataCM[locin,6]
	m_m0 = np.log10(mx)
	c_c0 = np.log10(cy)
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
	plt.plot(xcon,ycon,c = 'g')
	

	x = np.linspace(-nd,nd,100)  
	y = np.linspace(-nd,nd,100)
	cs1 = plt.contour(x,y,dataChi2,levels=[2.3],colors = 'red',lw = 0.8)
	pc1=cs1.collections[0].get_paths()[0]
	con1=pc1.vertices
	xcon1=con1[:,0]
	ycon1=con1[:,1]
	# plt.scatter(xcon1,ycon1,c = 'k', s = 5)
	nxbest = quantity_size(cov,xcon1,ycon1,nx0)
	# nxbest = nx0

	# plt.scatter(xcon1*nxbest,ycon1*nxbest,c = 'cyan', s = 5)
	plt.plot(xcon1*nxbest,ycon1*nxbest,c = 'red', ls = '--', label = 'multiple size')
	plt.legend()
	
	return nxbest

def main(model,nx0):
	fileCM = '/Users/liqy/Documents/data/oPDFnew/cluster/G3X-CM-masscomplete-ALL.txt'
	if model == 'JeansE':
		filefit = '/Users/liqy/Documents/data/oPDF1/JeansE/JeansE_Chi2fit_MCDM_rcin200_TMP_n100000_bin20_boot200.txt'
		fileChi2 = '/Users/liqy/Documents/data/oPDF1/JeansE/JeansE_Chi2_rcin200_TMP_n100000_bin20_boot200.txt'
		# filefit = '/Users/liqy/Desktop/JeansE_Chi2fit_MCDM_rcin200_TMP_n1000_bin20_subs200.txt'
		# fileChi2 = '/Users/liqy/Desktop/JeansE_Chi2_rcin200_TMP_n1000_bin20_subs200.txt'

		ntracer = float(1e5)
		nd = 0.02 #the grid 
	elif model == 'HE':
		filefit = '/Users/liqy/Documents/data/oPDF1/HE/HE_Chi2fit_MChotgas_rcin200_TMP_n100000_bin20_boot200.txt'
		fileChi2 = '/Users/liqy/Documents/data/oPDF1/HE/HE_Chi2_rcin200_TMP_n100000_bin20_boot200.txt'
		ntracer = float(1e5)
		nd = 0.01
		# fileChi2 = '/Users/liqy/Desktop/HE_Chi2.txt'
	elif model == 'VirM':
		filefit = '/Users/liqy/Documents/data/oPDF1/VirM/VirM_Chi2fit_MCDM_rcin200_TMP_n1000_bin20_boot200.txt'
		fileChi2 = '/Users/liqy/Documents/data/oPDF1/VirM/VirM_Chi2_rcin200_TMP_n1000_bin20_boot200.txt'		
		ntracer = float(1000)
		nd = 0.08

	dataCM = np.loadtxt(fileCM)
	datafit = np.loadtxt(filefit)
	dataChi2 = np.loadtxt(fileChi2)
	
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	plt.tick_params(top = 'on', right = 'on', which='both')
	
	mmtrue = datafit[:,0] / dataCM[:,2]
	cctrue = datafit[:,1] / dataCM[:,6]
	loc = opd.sigma3sel(np.log10(mmtrue),np.log10(cctrue))
	print(loc.shape)
	col = 'grey'
	nxbest = mc_all(loc,dataCM,datafit,dataChi2,col,nx0,nd)

	print("The best nx is ",nxbest)
	print("The Neff is ", ntracer / (nxbest**2-1))

	
def Main(model,nx0):
	plt.figure(figsize = (4.5,4))

	main(model,nx0)


	plt.legend(loc = 1,fontsize = 'x-small')

	plt.tick_params(labelsize = 10)
	plt.xlabel(r'$lg(M/M_{\rm true})$',fontsize = 12)
	plt.ylabel(r'$lg(c/c_{\rm true})$',fontsize = 12)

	bwith = 1.2 #边框宽度设置为2
	ax = plt.gca()#获取边框
	ax.spines['bottom'].set_linewidth(bwith)
	ax.spines['left'].set_linewidth(bwith)
	ax.spines['top'].set_linewidth(bwith)
	ax.spines['right'].set_linewidth(bwith)
	plt.axhline(0, ls = '--', c = 'k', lw = 1)
	plt.axvline(0, ls = '--', c = 'k', lw =1 )
	plt.xlim(-0.5,0.5)
	plt.ylim(-0.5,0.5)
	plt.title('%s' %model)
	plt.tight_layout()
	plt.show()


'''
nx0: initial value to the expansion multiple of statistic error 
'''
# Main(model = 'JeansE', nx0 = 42) 
Main(model = 'HE', nx0 = 23)
# Main(model = 'VirM', nx0 = 4)






