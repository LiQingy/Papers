import numpy as np
import matplotlib.pyplot as plt
import itertools as itert
from scipy.optimize import curve_fit

def massdensity(ax1,ax10,d,col,mark,lab1,lab2,linet):
	if d == 'GX' or d =='GM':
		file = '/home/qyli/Desktop/300data/AHF/AHFdata/%s_AHFxy_density4.txt' %d
		rr = '/home/qyli/Desktop/300data/AHF/AHFdata/%s_AHFxy_rr4.txt' %d
		if d == 'GX':
			filecen = '/home/qyli/Desktop/300data/center/G3X_Mass_snap_128-center-cluster.txt'
		else:
			filecen = '/home/qyli/Desktop/300data/center/Music_Mass_snap_017-center-cluster.txt'
		filegal200 = '/home/qyli/Desktop/300data/AHF/AHFdata/%s_AHFxy_mass200gal.txt' %d
	elif d == 'SDSS7':
		file = '/home/qyli/Desktop/300data/sdss7/Planck/sdss7_density4PlanckBCG.txt'
		rr = '/home/qyli/Desktop/300data/sdss7/Planck/sdss7_rr4PlanckBCG.txt'
		filecen = '/home/qyli/Desktop/300data/sdss7/Planck/sdss_group_centerPlanckBCG.txt'
		filegal200 = '/home/qyli/Desktop/300data/sdss7/Planck/sdss7_mass200gal.txt'
	else:
		file = '/home/qyli/Desktop/300data/SAM/SAMdata/SAM_%sxy_density4.txt' %d
		rr = '/home/qyli/Desktop/300data/SAM/SAMdata/SAM_%sxy_rr4.txt' %d
		filecen = '/home/qyli/Desktop/300data/center/MDPL2_Mass_snap_128-center-cluster.txt'
		filegal200 = '/home/qyli/Desktop/300data/SAM/SAMdata/SAM_%sxy_mass200gal.txt' %d
	data = np.loadtxt(file)
	datacen = np.loadtxt(filecen)
	gal200 = np.loadtxt(filegal200)
	rr = np.loadtxt(rr)

	cmedian = np.zeros(datacen.shape[0])
	for i in range(datacen.shape[0]):
	
		loc = np.where(np.isnan(rr[i,:]) == False)[0]
		x = rr[i,loc]
		y = np.log10(data[i,loc])

		#-----------------------------
		#fitting
		x0 = x[1:]
		y0 = y[1:]

		if d == 'SDSS7':
			r200 = datacen[i,4]
		else:
			r200 = datacen[i,6] / 0.678
		M200 = gal200[i]
	
		def func(x,b):
			return np.log10(M200 * (np.log(1 + x*b) - (x*b)/(1 + x*b))/ (np.log(1 + b) - (b)/(1 + b)) )

		coeffs, matcov = curve_fit(func, x0, y0, p0 = 5,bounds = (1,20))
		yfit = func(x, coeffs[0])
		cmedian[i] = coeffs[0]

		nv = np.pi * (r200 / 20 * (loc + 1))**2
		yfit = np.log10(10**yfit / nv)
		ytrue = np.log10(10**y / nv)
		resd = 10**ytrue/10**yfit - 1

		#residuals
		if lab1 == 'SDSS7':
			ax1.plot(x, ytrue, color = col, linestyle = linet,zorder = 2)
			ax10.plot(x, resd, color = col, linestyle = linet)
		else:
			ax1.plot(x, ytrue, color = col, linestyle = linet,zorder = 2)
			ax10.plot(x, resd, color = col, linestyle = linet)
	print(np.median(cmedian),np.std(cmedian,ddof = 1),np.max(cmedian))
		# if d == 'SDSS7':
		# 	ax1.errorbar(x, y, yerr = [y - yerrdown, yerrup - y], fmt = mark, markeredgewidth = 1.7, markersize = 7,color = col,zorder = 3,capsize=4, elinewidth = 1.7,label = lab1)
		# else:
		# 	ax1.fill_between(x,yerrdown,yerrup,alpha = 0.25,color = col,zorder = 1)
		# 	ax1.scatter(x,y,marker = mark,s = 6**2,color = col,zorder = 3,label = lab1)

def main():
	#---------------------------------------------------------------------------------
	sym = ['GX','GM','Galacticus','SAG','SAGE','SDSS7']
	# sym = ['GX']
	lab11 = itert.cycle(['Gadget-X','Gadget-MUSIC','Galacticus','SAG','SAGE','SDSS7'])
	lab22 = itert.cycle(['fit to Gadget-X','fit to Gadget-MUSIC','fit to Galacticus','fit to SAG','fit to SAGE','fit to SDSS7'])
	linett = itert.cycle(['-','-','--','--','--',':'])
	coll = itert.cycle(['red','blue','g','orange','k','c'])
	markk = itert.cycle(['o','^','d','s','*','x'])

	for d in sym:
		plt.rcParams['xtick.direction'] = 'in'
		plt.rcParams['ytick.direction'] = 'in'
		fig = plt.figure(dpi =100,figsize=(9, 8))

		ax1 = plt.subplot2grid((8,8),(0,0),colspan = 8,rowspan = 6)
		ax10 = plt.subplot2grid((8,8),(6,0),colspan = 8,rowspan = 2)
		col = next(coll)
		mark = next(markk)
		lab1 = next(lab11)
		lab2 = next(lab22)
		linet = next(linett)
		massdensity(ax1,ax10,d,col,mark,lab1,lab2,linet)

	#---------------------------------------------------------------------------------
		ax1.tick_params(top = 'on', right = 'on',which='both',labelsize = 11)
		# ax1.set_yscale('log')
		ax1.set_ylabel(r'$Log_{10}\ \rho_{\bigstar}(<r)\ [M_⊙kpc^{-2}]$',fontsize = 14)
		ax1.set_xlim(0,1)
		ax1.set_ylim(5,8)
		ax1.set_title('%s' %d)
		ax1.axes.xaxis.set_ticklabels([])

		ax10.tick_params(top = 'off', right = 'on',which='both',labelsize = 11)
		ax10.set_xlabel('$r_p$/$r_{200}$',fontsize = 14)
		ax10.set_ylabel(r'$\rho_{\bigstar}/\rho_{fit}-1$',fontsize = 14)
		ax10.set_xlim(0,1)
		ax10.set_ylim(-1,1)
		ax10.axhline(0,c='grey',zorder = 30)
		ax10.set_yticks([-0.8,-0.4,0,0.4,0.8])

		plt.tight_layout()
		plt.subplots_adjust(wspace =0 ,hspace = 0)

		ax1.legend(loc = 'best', ncol = 2,fontsize = 'small')
		# plt.savefig('/home/qyli/Desktop/test/300cluster/stellar_density_%s.png' %d)
		# plt.show()
main()