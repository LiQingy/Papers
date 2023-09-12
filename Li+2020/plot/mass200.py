import numpy as np
import matplotlib.pyplot as plt
import itertools as itert
from scipy.optimize import curve_fit

def massdensity(ax1,ax10):

	sym = ['GX','GM','Galacticus','SAG','SAGE','SDSS7']
	lab11 = itert.cycle(['Gadget-X','Gadget-MUSIC','Galacticus','SAG','SAGE','SDSS7'])
	lab22 = itert.cycle(['fit to Gadget-X','fit to Gadget-MUSIC','fit to Galacticus','fit to SAG','fit to SAGE','fit to SDSS7'])
	linett = itert.cycle(['-','-','--','--','--',':'])
	coll = itert.cycle(['red','blue','g','orange','k','c'])
	markk = itert.cycle(['o','^','d','s','*','x'])

	for d in sym:
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

		col = next(coll)
		mark = next(markk)
		lab1 = next(lab11)
		lab2 = next(lab22)
		linet = next(linett)

		#-----------------------------
		#fitting
		x0 = x[1:]
		y0 = y[1:]
		

		M200 = np.median(gal200)
		xline = np.linspace(0,1,100)		
		def func(x,b):
			return np.log10(M200 * (np.log(1 + x*b) - (x*b)/(1 + x*b))/ (np.log(1 + b) - (b)/(1 + b)) )

		coeffs, matcov = curve_fit(func, x0, y0, p0 = 20,bounds = (1,20))
		yline = func(xline, coeffs[0])
		print(coeffs,M200 / 1e12)

		#residuals
		if lab1 == 'SDSS7':
			ax1.plot(xline, yline, color = col, linestyle = linet,linewidth = 2.2,zorder = 2,label = lab2)
			ax10.plot(x,yerrup - y,color = col, linestyle = linet)
			ax10.plot(x,yerrdown - y,color = col, linestyle = linet)
		else:
			ax1.plot(xline, yline, color = col, linestyle = linet,linewidth = 1.5,zorder = 2,label = lab2)
			ax10.plot(x,yerrup - y,color = col, linestyle = linet)
			ax10.plot(x,yerrdown - y,color = col, linestyle = linet)

		if d == 'SDSS7':
			ax1.errorbar(x,y,yerr = [y - yerrdown, yerrup - y],fmt = mark,markeredgewidth = 1.7,markersize = 7,color = col,zorder = 3,capsize=4, elinewidth = 1.7,label = lab1)
		else:
			ax1.fill_between(x,yerrdown,yerrup,alpha = 0.25,color = col,zorder = 1)
			ax1.scatter(x,y,marker = mark,s = 6**2,color = col,zorder = 3,label = lab1)


def main():
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	# fig = plt.figure(dpi =100,figsize=(12, 6))

	ax1 = plt.subplot2grid((8,8),(0,0),colspan = 8,rowspan = 6)
	ax10 = plt.subplot2grid((8,8),(6,0),colspan = 8,rowspan = 2)

	massdensity(ax1,ax10)

	ax1.tick_params(top = 'on', right = 'on',which='both',labelsize = 11)
	# ax1.set_yscale('log')
	ax1.set_ylabel(r'$Log_{10}\ M_{\bigstar}(<r)\ [M_⊙]$',fontsize = 14)
	ax1.set_xlim(0,1)
	ax1.set_ylim(10,14)
	ax1.axes.xaxis.set_ticklabels([])

	ax10.tick_params(top = 'off', right = 'on',which='both',labelsize = 11)
	ax10.set_xlabel('$r_p$/$r_{200}$',fontsize = 14)
	ax10.set_ylabel('scatter',fontsize = 14)
	ax10.set_xlim(0,1)
	ax10.set_ylim(-0.5,0.5)
	ax10.axhline(0,c='grey',zorder = 30)
	ax10.set_yticks([-0.3,0,0.3])

	plt.tight_layout()
	plt.subplots_adjust(wspace =0 ,hspace = 0)

	ax1.legend(loc = 'best', ncol = 2,fontsize = 'small')
	plt.savefig('/home/qyli/Desktop/stellar_mass.pdf')
	plt.show()
main()