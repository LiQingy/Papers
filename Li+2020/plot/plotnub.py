import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.tick_params(top = 'on', right = 'on', which='both')

def plotnn(datax,data,line,colorn,labeln,sym,lab2):
	
	x = []
	y = []
	yerrdown = []
	yerrup = []
	for i in range(20):
		xx = []
		yy = []
		nn = 0
		for j in range(data.shape[0]):
			if data[j][i] != 0: # and math.isnan(datax[j][i]) == False
				yy.extend([data[j][i]])
				xx.extend([datax[j][i]])
				nn += 1
		yy = np.array(yy)
		xx = np.array(xx)
		if nn > 1:
			if labeln == 'SDSS7':
				xn = np.median(xx)
				yn = np.median(yy)
			else:
				xn = np.median(xx)
				yn = np.median(yy)
			x.extend([xn])
			y.extend([yn])
			yerr = np.percentile(yy, [16, 84])
			yerrdown.extend([yerr[0]])
			yerrup.extend([yerr[1]])
	x = np.array(x)
	y = np.array(y)
	yerrdown = np.array(yerrdown)
	yerrup = np.array(yerrup)
	if labeln == 'SDSS7':
		plt.plot(x, y, linestyle = line, color = colorn, label = labeln,linewidth = 2.5 )
	else:
		plt.plot(x, y, linestyle = line, color = colorn, label = labeln)
	plt.fill_between(x, yerrdown, yerrup , color = colorn, alpha = 0.25)

def plotnumber_all():

	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'


	file4 = '/home/qyli/Desktop/300data/AHF/AHFdata/GX_AHFxy_totaln4.txt'
	filex = '/home/qyli/Desktop/300data/AHF/AHFdata/GX_AHFxy_rr4.txt'
	data = np.loadtxt(file4)
	datax = np.loadtxt(filex)
	plotnn(datax,data,line = '-',colorn = 'red',labeln = 'Gadget-X',sym = 'o',lab2 = 'fit to Gadget-X')


	file5 = '/home/qyli/Desktop/300data/AHF/AHFdata/GM_AHFxy_totaln4.txt'
	filex = '/home/qyli/Desktop/300data/AHF/AHFdata/GM_AHFxy_rr4.txt'
	data = np.loadtxt(file5)
	datax = np.loadtxt(filex)
	plotnn(datax,data,line = '-',colorn = 'blue',labeln = 'Gadget-MUSIC',sym = '^',lab2 = 'fit to Gadget-MUSIC')
		
	file1 = '/home/qyli/Desktop/300data/SAM/SAMdata/SAM_Galacticusxy_totaln4.txt'
	filex = '/home/qyli/Desktop/300data/SAM/SAMdata/SAM_Galacticusxy_rr4.txt'
	data = np.loadtxt(file1)
	datax = np.loadtxt(filex)
	plotnn(datax,data,line = '--',colorn = 'green',labeln = 'Galacticus',sym = 'd',lab2 = 'fit to Galacticus')


	file2 = '/home/qyli/Desktop/300data/SAM/SAMdata/SAM_SAGxy_totaln4.txt'
	filex = '/home/qyli/Desktop/300data/SAM/SAMdata/SAM_SAGxy_rr4.txt'
	data = np.loadtxt(file2)
	datax = np.loadtxt(filex)
	plotnn(datax,data,line = '--',colorn = 'orange',labeln = 'SAG',sym = 's',lab2 = 'fit to SAG')


	file3 = '/home/qyli/Desktop/300data/SAM/SAMdata/SAM_SAGExy_totaln4.txt'
	filex = '/home/qyli/Desktop/300data/SAM/SAMdata/SAM_SAGExy_rr4.txt'
	data = np.loadtxt(file3)
	datax = np.loadtxt(filex)
	plotnn(datax,data,line = '--',colorn = 'black',labeln = 'SAGE',sym = '*',lab2 = 'fit to SAGE')


	file6 = '/home/qyli/Desktop/300data/sdss7/Planck/sdss7_totaln4PlanckBCG.txt'
	filex = '/home/qyli/Desktop/300data/sdss7/Planck/sdss7_rr4PlanckBCG.txt'
	data = np.loadtxt(file6)
	datax = np.loadtxt(filex)
	plotnn(datax,data,line = ':',colorn = 'c',labeln = 'SDSS7',sym = 'x',lab2 = 'fit to SDSS7')


	plt.xlim(0,1)
	plt.yscale('log')
	plt.ylim(3e-7,1e-4)
	plt.tick_params(top = 'on', right = 'on', which='both',labelsize = 11)
	plt.xlabel('$r_p$/$r_{200}$',fontsize = 14, weight = 'heavy')
	plt.ylabel(r'$galaxy\ number\ density\ [kpc^{-2}]$', fontsize = 14,fontweight = 800)
	plt.text(0.2, 10**-4.4, r'$M_{\bigstar} > 5 \times 10^{10}M_⊙$', fontsize = 16, fontweight = 1000)
	# plt.text(0.23, 10**-7.5, r'$SDSS7:\ 0.01\ <\ z\ <\ 0.15$', fontsize = 12, fontweight = 'heavy')
	plt.legend(loc = 'best',ncol = 1)
	plt.savefig('/home/qyli/Desktop/stellar_n.pdf')
	plt.show()

plotnumber_all()


