import numpy as np
import matplotlib.pyplot as plt
import itertools as itert
from scipy.optimize import curve_fit

def massdensity(ax1,ax10,ax20,resid):

	sym = ['GX','GM','Galacticus','SAG','SAGE','SDSS7']
	lab11 = np.array(['Gadget-X','Gadget-MUSIC','Galacticus','SAG','SAGE','SDSS7'])
	lab22 = np.array(['fit to Gadget-X','fit to Gadget-MUSIC','fit to Galacticus','fit to SAG','fit to SAGE','fit to SDSS7'])
	linett = np.array(['-','-','--','--','--',':'])
	coll = np.array(['red','blue','g','orange','k','c'])
	markk = np.array(['o','^','d','s','*','x'])
	cc = np.array([6.2524286,9.7396373,3.1278967,3.23723932,2.3286938948,6.13155616])

	datax = np.zeros(shape = (6,20))
	datay = np.zeros(shape = (6,20))
	nd = 0
	for d in sym:
		if d == 'GX' or d =='GM':
			file = '//Users/liqy/Documents/data/300data/AHF/AHFdata/%s_AHFxy_density4.txt' %d
			rr = '/Users/liqy/Documents/data/300data/AHF/AHFdata/%s_AHFxy_rr4.txt' %d
			if d == 'GX':
				filecen = '/Users/liqy/Documents/data/300data/center/G3X_Mass_snap_128-center-cluster.txt'
			else:
				filecen = '/Users/liqy/Documents/data/300data/center/Music_Mass_snap_017-center-cluster.txt'
			filegal200 = '/Users/liqy/Documents/data/300data/AHF/AHFdata/%s_AHFxy_mass200gal.txt' %d
		elif d == 'SDSS7':
			file = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_density4PlanckBCG.txt'
			rr = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_rr4PlanckBCG.txt'
			filecen = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss_group_centerPlanckBCG.txt'
			filegal200 = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_mass200gal.txt'
		else:
			file = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_density4.txt' %d
			rr = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_rr4.txt' %d
			filecen = '/Users/liqy/Documents/data/300data/center/MDPL2_Mass_snap_128-center-cluster.txt'
			filegal200 = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_mass200gal.txt' %d
		data = np.loadtxt(file)
		datacen = np.loadtxt(filecen)
		gal200 = np.loadtxt(filegal200)
		rr = np.loadtxt(rr)

		y = np.zeros(20)
		x = np.zeros(20)
		yerrdown = np.zeros(20)
		yerrup = np.zeros(20)
		for k in range(20):
			loc = np.where(np.isnan(rr[:,k]) == False)[0]
			x[k] = np.median(rr[loc,k])
			y[k] = np.median(np.log10(data[loc,k]))
			yerr = np.percentile(np.log10(data[loc,k]), [16, 84])
			yerrdown[k] = yerr[0]
			yerrup[k] = yerr[1]

		col = coll[nd]
		mark = markk[nd]
		lab1 = lab11[nd]
		lab2 = lab22[nd]
		linet = linett[nd]
		ec = cc[nd]
		print(x)
		#-----------------------------
		#fitting
		x0 = x[1:]
		y0 = y[1:]
		
		pcr = 2.775 * 10**11 * 0.678**2 * 1e-9
		# r200 =  np.power(3 * M200 / 800/ np.pi / pcr,1/3)
		if d == 'SDSS7':
			r200 = np.median(datacen[:,4])
			M200 = np.median(gal200)
			# print(np.std(gal200,ddof = 1)/1e12,np.std(datacen[:,4],ddof = 1))
		else:
			r200 = np.median(datacen[:,6]) / 0.678
			M200 = np.median(gal200)
			# print(np.std(gal200,ddof = 1)/1e12,np.std(datacen[:,6] / 0.678,ddof = 1))
		xline = np.linspace(0,1,100)		
		def func(x,b):
			return np.log10(M200 / (np.pi * (x * r200)**2) * (np.log(1 + x*b) - (x*b)/(1 + x*b))/ (np.log(1 + b) - (b)/(1 + b)) )
		
		# def func(x,a,b):
		# 	return np.log10(a *  (np.log(1 + x*b) - (x*b)/(1 + x*b))/ (np.log(1 + b) - (b)/(1 + b)) / x **2)

		coeffs, matcov = curve_fit(func, x0, y0, p0 = 20,bounds = (1,20))
		yline = func(xline, coeffs[0])
		print(coeffs[0],r200/1e2,M200/1e12)

		res = 10**y/10**func(x, coeffs[0])- 1

		#residuals
		if d == 'SDSS7':
			ax1.errorbar(x,y,yerr = [y - yerrdown, yerrup - y],fmt = mark,markeredgewidth = 1.7,markersize = 7,color = col,zorder = 3,capsize=4, elinewidth = 1.7,label = lab1)
			ax1.plot(xline, yline, color = col, linestyle = linet,linewidth = 2.2,zorder = 2,label = lab2)
		else:
			ax1.scatter(x,y,marker = mark,s = 6**2,color = col,zorder = 3,label = lab1)
			ax1.fill_between(x,yerrdown,yerrup,alpha = 0.25,color = col,zorder = 1)
			ax1.plot(xline, yline, color = col, linestyle = linet,linewidth = 1.5,zorder = 2,label = lab2)
		
		if resid:
			ax10.plot(x,res,color = col, linestyle = linet,linewidth = 1.5)
		else:
			ax10.plot(x,yerrup - y,color = col, linestyle = linet)
			ax10.plot(x,yerrdown - y,color = col, linestyle = linet)

		datax[nd] = x
		datay[nd] = y
		nd += 1

	#ax20
	for i in range(6):
		ax20.plot(datax[i],10**datay[i]/10**datay[5]-1,color = coll[i],linestyle = linett[i])


def main(resid):
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	fig = plt.figure(dpi =100,figsize=(7.4,7.5))

	ax1 = plt.subplot2grid((12,8),(0,0),colspan = 8,rowspan = 8)
	ax10 = plt.subplot2grid((12,8),(8,0),colspan = 8,rowspan = 2)
	ax20 = plt.subplot2grid((12,8),(10,0),colspan = 8,rowspan = 2)
	
	massdensity(ax1,ax10,ax20,resid)

	ax1.tick_params(top = 'on', right = 'on',which='both',labelsize = 11)
	# ax1.set_yscale('log')
	ax1.set_ylabel(r'$Log_{10}\ \rho_{\bigstar}(<r)\ [M_⊙kpc^{-2}]$',fontsize = 14)
	ax1.set_xlim(0,1)
	ax1.set_ylim(5,8)
	ax1.axes.xaxis.set_ticklabels([])

	ax10.tick_params(top = 'off', right = 'on',which='both',labelsize = 11)
	ax10.set_xlabel('$r_p$/$r_{200}$',fontsize = 14)
	ax10.set_xlim(0,1)
	ax10.set_ylim(-0.5,0.5)
	ax10.axhline(0,c='grey',zorder = 30)
	ax10.set_yticks([-0.3,0,0.3])
	ax10.axes.xaxis.set_ticklabels([])

	ax20.tick_params(top = 'off', right = 'on',which='both',labelsize = 11)
	ax20.set_xlabel('$r_p$/$r_{200}$',fontsize = 14)
	ax20.set_ylabel(r'$\rho/\rho_{SDSS7} - 1$',fontsize = 14)
	ax20.set_xlim(0,1)
	ax20.set_ylim(-1,10)
	ax20.set_yticks([0,5,10])

	if resid:
		ax10.set_ylabel(r'$\rho_{\bigstar}/\rho_{fit}-1$',fontsize = 14)
	else:
		ax10.set_ylabel('scatter',fontsize = 14)

	plt.tight_layout()
	plt.subplots_adjust(wspace =0 ,hspace = 0)

	ax1.legend(loc = 'best', ncol = 2)
	plt.savefig('/Users/liqy/Desktop/stellar_density.pdf')
	plt.show()

main(resid = False)

