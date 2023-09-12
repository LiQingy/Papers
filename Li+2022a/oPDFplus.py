import numpy as np 
import h5py
# import matplotlib.pyplot as plt

# def Nto1sigma(nclu,datasigma,mcfit,avecen,nsigma,ntour):
	# Amplify combined clusters sigma to 1 cluster simga  
	# x0 = mcfit #fitdata
	# x1 = avecen #average true data

	# ngrids=[20,20]
	# dx=[0.1,0.1]
	# xl=np.log10(x0)-dx
	# xu=np.log10(x0)+dx
	# mline=np.logspace(xl[0], xu[0], ngrids[0])
	# cline=np.logspace(xl[1], xu[1], ngrids[1])

	# '''Get the contour value'''
	# cs = plt.contour(mline/x1[0],cline/x1[1],datasigma,levels = [nsigma])
	# p = cs.collections[0].get_paths()[ntour]
	# v = p.vertices
	# x = v[:,0]
	# y = v[:,1]

	# amp = pow(nclu,1/2)
	# c01 = np.log10(x0 / x1)

	# xnew = c01[0] - (c01[0] - np.log10(x)) * amp
	# ynew = c01[1] - (c01[1] - np.log10(y)) * amp
	# return 10**xnew,10**ynew

def sigma3sel(m_m0,c_c0):	
	mc = np.array([m_m0,c_c0]).T
	cluid = np.arange(0,324,1, dtype = 'i4')
	Nsel = 324
	
	for i in range(10):
		pcov = np.cov(m_m0[cluid],c_c0[cluid])
		locsel = []
		for i in range(Nsel):
			nclu = cluid[i] 
			x = mc[nclu]
			if (x-mc[cluid].mean(0)).dot(np.linalg.inv(pcov)).dot(x-mc[cluid].mean(0))<11.8:
				locsel.append(nclu)

		locsel = np.array(locsel)
		Nsel = locsel.shape[0]
		#print(Nsel)
		if np.sum(locsel) == np.sum(cluid):
			break
		cluid = locsel
	return locsel

def sigma2dis(m_m0,c_c0,prisig):	
	mave = np.mean(m_m0)
	cave = np.mean(c_c0)
	msigma = np.std(m_m0,ddof = 1)
	csigma = np.std(c_c0,ddof = 1)
	pcov = np.cov(m_m0,c_c0)[0][1] / msigma / csigma
	if prisig == True:
		print(np.cov(m_m0,c_c0),msigma,csigma)
		print('The correlation coefficience is %s' %pcov)
	# pos = np.zeros(shape = (m_m0.shape[0],2))
	# pos[:,0] = datafit[loc,0] / dataCM[:,2] * 1e15
	# pos[:,1] = datafit[loc,1] / dataCM[:,6]
	# cov=np.cov(np.log10(pos).T)
	# l1=plot_cov_ellipse(cov, np.log10(pos).mean(axis=0), color='green', fill=0)
	# print(np.log10(pos).mean(axis=0),mave,cave)

	def gaussian2(xx,yy):
		A = 1 / 2 / np.pi / msigma / csigma / np.sqrt(1 - pcov**2)
		B = (xx - mave)**2 / msigma**2 + (yy - cave)**2 / csigma**2 - 2*pcov*(xx - mave)*(yy - cave) / msigma / csigma
		return A * np.exp(-1 / 2 / (1 - pcov**2) * B)
	
	xx = np.linspace(-3,2,400)
	yy = np.linspace(-3,2,400)
	X,Y = np.meshgrid(xx,yy)
	Z = gaussian2(X,Y)
	
	k = 2.3
	k1 = np.exp(-k / 2) / 2 / np.pi / msigma / csigma / np.sqrt(1 - pcov**2)
	sigmalevel = [10**k1]
	return 10**X,10**Y,10**Z,sigmalevel
	






	

#-------------------------------------------------------
