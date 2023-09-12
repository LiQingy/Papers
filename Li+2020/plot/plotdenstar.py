import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit,minimize
from scipy.optimize import fmin as simplex
# import lnr
global x
global y
from numba import vectorize
def median(phy,ax1,ax2,ax10,ax20,n,data,datax,mark,col,linet,lab1,lab2):

	x = []
	y = []
	ydown = []
	yup = []

	for k1 in range(n):
		yy = []
		xx = []
		nx = 0
		for k2 in range(data.shape[0]):
			if data[k2][k1] != 0 and math.isnan(data[k2][k1]) == False:
				xx.extend([datax[k2][k1]])
				yy.extend([math.log10(data[k2][k1])])
				nx += 1
		
		if nx > 1:
			yy = np.array(yy)
			xx = np.array(xx)
			
			ymedian = np.median(yy)
			xmedian = np.median(xx)
			yerr = np.percentile(yy, [16,84])

			x.extend([xmedian])
			y.extend([ymedian])
			ydown.extend([yerr[0]])
			yup.extend([yerr[1]])

	
	x = np.array(x)
	y = np.array(y)
	ydown = np.array(ydown)
	yup = np.array(yup)



	if lab2 == 'fit to SDSS7':
		ax1.errorbar(x,y,yerr = [y - ydown, yup - y],fmt = mark,markeredgewidth = 1.7,markersize = 7,color = col,zorder = 3,capsize=4, elinewidth = 1.7,label = lab1)
		ax2.errorbar(x,y,yerr = [y - ydown, yup - y],fmt = mark,markeredgewidth = 1.7,markersize = 7,color = col,zorder = 3,capsize=4, elinewidth = 1.7,label = lab1)
	else:
		ax1.fill_between(x,ydown,yup,alpha = 0.25,color = col,zorder = 1)
		ax1.scatter(x,y,marker = mark,s = 6**2,color = col,zorder = 3,label = lab1)

		ax2.fill_between(x,ydown,yup,alpha = 0.25,color = col,zorder = 1)
		ax2.scatter(x,y,marker = mark,s = 6**2,color = col,zorder = 3,label = lab1)

	#---------------------------------------------------------------------------------------------
	#---------------------------------------------------------------------------------------------
	#fitting

	x0 = np.array(x)[1:]
	y0 = np.array(y)[1:]
	ydown0 = ydown[1:]
	yup0 = yup[1:]
	sigma = (ydown0 + yup0) / 2
	#'SEC'
	xline = np.linspace(0,1,100)
	#--------------------------
	# def func1(params, X, Y, Err):
	# # extract current values of fit parameters from input array
	# 	a = params[0]
	# 	b = params[1]
	# 	c = params[2]
	# 	# compute chi-square
	# 	chi2 = 0.0
	# 	for n in range(len(X)):
	# 		x = X[n]
	# 		# The function y(x)=a+b*x+c*x^2 is a polynomial
	# 		# in this example.
	# 		y = np.log10(a * np.exp(-b * (pow(x / c, 1/4) - 1)))

	# 		chi2 = chi2 + (Y[n] - y)*(Y[n] - y)/(Err[n]*Err[n])
	# 	return chi2
	# coeffs = simplex(func1, x0 = [5,15,40], args=(x0, y0, sigma), full_output=0)

	def func11(x,a,b,c):
		return np.log10(a * np.exp(-b * (pow(x / c, 1/4) - 1)))
	coeffs, matcov = curve_fit(func11, x0, y0)
	yline = func11(xline, coeffs[0],coeffs[1],coeffs[2])
	print(lab1,np.round(coeffs,decimals = 3))
	
	#residuals
	res1 = 10**y/10**func11(x, coeffs[0],coeffs[1],coeffs[2]) - 1
	if lab1 == 'SDSS7':
		ax1.plot(xline, yline, color = col, linestyle = linet,linewidth = 2.2,zorder = 2,label = lab2)
		ax10.plot(x,res1,color = col, linestyle = linet,linewidth = 1.9)
	else:
		ax1.plot(xline, yline, color = col, linestyle = linet,linewidth = 1.5,zorder = 2,label = lab2)
		ax10.plot(x,res1,color = col, linestyle = linet)

	#-------------------------------------------------
	#'NFW'
	def func2(x,a,b):
		return np.log10(a / (x/b * (1 + x/b)**2))

	# if lab1 != 'SAGE' and lab1 != 'Galacticus':

	# coeffs, matcov = curve_fit(func2, x0, y0,sigma = (ydown[1:] + yup[1:]) / 2,maxfev= 5000)
	# a = 1e4
	# b = 3
	# def y1(a,b):
	# 	return np.sum(np.log10(a / (x0/b * (1 + x0/b)**2)))
	# print(lnr.bces(x0,y1(a,b)))

#------------------------------------------------------------------------------
	def func(params, X, Y, Err):
	# extract current values of fit parameters from input array
		a = params[0]
		b = params[1]
		# compute chi-squaresss
		chi2 = 0.0
		for n in range(len(X)):
			x = X[n]
			# The function y(x)=a+b*x+c*x^2 is a polynomial
			# in this example.
			y = np.log10(a / (x/b * (1 + x/b)**2))

			chi2 = chi2 + (Y[n] - y)*(Y[n] - y)/(Err[n]*Err[n])
		return chi2
	coeffs = simplex(func, x0 = [2.92e5,20], args=(x0, y0, sigma), full_output=0, maxiter = 150)
	
	print(coeffs)
	yline = func2(xline, coeffs[0],coeffs[1])
	#residuals
	res2 = 10**y/10**func2(x, coeffs[0],coeffs[1]) - 1
	if lab1 == 'SDSS7':
		ax2.plot(xline, yline, color = col, linestyle = linet,linewidth = 2.2,zorder = 2,label = lab2)
		ax20.plot(x,res2,color = col, linestyle = linet,linewidth = 1.9)
	else:
		ax2.plot(xline, yline, color = col, linestyle = linet,linewidth = 1.5,zorder = 2,label = lab2)
		ax20.plot(x,res2,color = col, linestyle = linet)
#-----------------------------------------------------------------------------

	# else:
	# 	ydown0 = np.array(ydown)[1:]
	# 	yup0 = np.array(yup)[1:]
	# 	xline = np.linspace(0,1,100)
	# 	if lab1 == 'Galacticus':
	# 		x00 = [3e4, 10]
	# 		mod = 'Powell'
	# 	else:
	# 		x00 = [1e3,20]
	# 		mod = 'Powell'
	# 	def func3(para):
	# 		a = para[0]
	# 		b = para[1]
	# 		return np.sum((y0 - np.log10(a / (x0/b * (1 + x0/b)**2)))**2 / ((ydown0 + yup0 - 2 * y0) / 2)**2)
	# 	coeff3 = minimize(func3,x0 = x00,method = mod)
	# 	print(lab2,coeff3.x)
	# 	yline = func2(xline,coeff3.x[0],coeff3.x[1])
	# 	ax2.plot(xline, yline, color = col, linestyle = linet,linewidth = 1.5,zorder = 2,label = lab2)
	# 	res2 = 10**y/10**func2(x, coeff3.x[0],coeff3.x[1]) - 1
	# 	ax20.plot(x,res2,color = col, linestyle = linet)

	#--------------------------------------------------------------------------------------------------
	# @vectorize
	# def func2DNFW(x,rc,c):
	# 	pcr = 2.775 * 10**11 / 0.678**2 * 1e-9
	# 	# rc = (3 * M / 800 / np.pi/ pcr) **(1/3)	/ c
	# 	M = np.pi * 4 /3 * 200 * pcr * (rc * c)**3
	# 	xx = x * c
	# 	Ic = 1/c**3 * (np.log(1 + c) + 1/(1+c) -1)

	# 	if xx < 1:
	# 		fx = 1/(xx**2 -1) *(1 - np.log((1 + (1 - xx**2)**0.5)/xx)/(1 - xx**2)**0.5)
	# 		gx = 2/xx**2 *(np.log(xx/2) + np.log((1 + (1 - xx**2)**0.5)/xx)/(1 - xx**2)**0.5)
	# 	elif xx == 1:
	# 		fx = 1/3
	# 		gx = 2 + 2 * np.log(1/2)
	# 	else:
	# 		fx = 1/(xx**2 -1) *(1 - np.arctan((xx**2 - 1)**0.5)/(xx**2 - 1)**0.5)
	# 		gx = 2/xx**2 *(np.log(xx/2) + np.arctan((xx**2 - 1)**0.5)/(xx**2 - 1)**0.5)
	# 	return np.log10(M / (np.pi * 2 * rc**2) / Ic * (gx-fx))
	# coeffs,matcov = curve_fit(func2DNFW,x0,y0,p0 = [1,1])
	# print(coeffs)
	# xline = np.linspace(0.001,1,100)
	# # coeffs = [1e4,10]
	# yline = func2DNFW(xline,coeffs[0],coeffs[1])
	# res2 = 10**y/10**func2DNFW(x, coeffs[0],coeffs[1]) - 1
	# if lab1 == 'SDSS7':
	# 	ax2.plot(xline, yline, color = col, linestyle = linet,linewidth = 2.2,zorder = 2,label = lab2)
	# 	ax20.plot(x,res2,color = col, linestyle = linet,linewidth = 1.9)
	# else:
	# 	ax2.plot(xline, yline, color = col, linestyle = linet,linewidth = 1.5,zorder = 2,label = lab2)
	# 	ax20.plot(x,res2,color = col, linestyle = linet)

#------------------------------------------------
#------------------------------------------------
#------------------------------------------------
def Simulation(phy,SIMmodel,ax1,ax2,ax10,ax20):

		
	file = '/home/qyli/Desktop/300data/simulation/%s500_snap_density0.txt' %(SIMmodel)
	data = np.loadtxt(file)

	if SIMmodel == 'GX':
		filex = '/home/qyli/Desktop/300data/simulation/GX500_snap_rr0.txt'
		col = 'red'
		mark = 'o'
		linet = '-'
		lab1 = 'Gadget-X'
		lab2 = 'fit to Gadget-X'
	
	if SIMmodel == 'GM':
		filex = '/home/qyli/Desktop/300data/simulation/GM500_snap_rr0.txt'
		col = 'blue'
		mark = '^'
		linet = '-'
		lab1 = 'Gadget-MUSIC'
		lab2 = 'fit to Gadget-MUSIC'	

	data = data * 10**10 * 0.678**2
	datax = np.loadtxt(filex)

	median(phy,ax1,ax2,ax10,ax20,20,data,datax,mark,col,linet,lab1,lab2)

def AHF(phy,SIMmodel,ax1,ax2,ax10,ax20):
	if SIMmodel == 'GX':
		col = 'red'
		mark = 'o'
		linet = '-'
		lab1 = 'Gadget3-X'
		lab2 = 'fit to Gadget-X'

	if SIMmodel == 'GM':
		col = 'blue'
		mark = '^'
		linet = '-'
		lab1 = 'Gadget3-MUSIC'
		lab2 = 'fit to Gadget-MUSIC'

	if phy == 'total':
		file = '/home/qyli/Desktop/300data/AHF/AHFdata/%s500_AHF_densityall.txt' %(SIMmodel)
		filex = '/home/qyli/Desktop/300data/AHF/AHFdata/%s500_AHF_rrall.txt' %SIMmodel
		data = np.loadtxt(file) * 10**10 / 0.678
	else:
		file = '/home/qyli/Desktop/300data/AHF/AHFdata/%s_AHFxy_density4.txt' %(SIMmodel)
		filex = '/home/qyli/Desktop/300data/AHF/AHFdata/%s_AHFxy_rr4.txt' %SIMmodel
		data = np.loadtxt(file)
	datax = np.loadtxt(filex)
	

	median(phy,ax1,ax2,ax10,ax20,20,data,datax,mark,col,linet,lab1,lab2)

def SAM(phy,SAMmodel,ax1,ax2,ax10,ax20):

	if SAMmodel == 'Galacticus':
		mark = 'd'
		linet = '--'
		col = 'green'
		lab1 = 'Galacticus'
		lab2 = 'fit to Galacticus'
	elif SAMmodel == 'SAG':
		mark = 's'
		linet = '--'
		col = 'orange'
		lab1 = 'SAG'
		lab2 = 'fit to SAG'
	elif SAMmodel == 'SAGE':
		mark = '*'
		linet = '--'
		col = 'black'
		lab1 = 'SAGE'
		lab2 = 'fit to SAGE'
	elif SAMmodel == 'MDPL':
		mark = '*'
		linet = '--'
		col = 'black'
		lab1 = 'MDPL2'
		lab2 = 'fit to MDPL2'

	if phy == 'total':
		file = '/home/qyli/Desktop/300data/SAM/SAMdata/SAM500_MDPL_densityall.txt' 
		filex = '/home/qyli/Desktop/300data/SAM/SAMdata/SAM500_MDPL_rrall.txt' 
	else:
		file = '/home/qyli/Desktop/300data/SAM/SAMdata/SAM_%sxy_density4.txt' %(SAMmodel)
		filex = '/home/qyli/Desktop/300data/SAM/SAMdata/SAM_%sxy_rr4.txt' %SAMmodel
	data = np.loadtxt(file)
	datax = np.loadtxt(filex)
	median(phy,ax1,ax2,ax10,ax20,20,data,datax,mark,col,linet,lab1,lab2)

def SDSS7(phy,ax1,ax2,ax10,ax20):

	
	filedata = '/home/qyli/Desktop/300data/sdss7/Planck/sdss7_density4PlanckBCG.txt' 
	filex = '/home/qyli/Desktop/300data/sdss7/Planck/sdss7_rr4PlanckBCG.txt'

	data = np.loadtxt(filedata)
	datax = np.loadtxt(filex)

	mark = 'x'
	linet = ':'
	col = 'c'
	lab1 = 'SDSS7'
	lab2 = 'fit to SDSS7'
	median(phy,ax1,ax2,ax10,ax20,20,data,datax,mark,col,linet,lab1,lab2)

def main(phy):
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	fig = plt.figure(dpi =100,figsize=(12, 6))

	ax1 = plt.subplot2grid((7,14),(0,0),colspan = 7,rowspan = 5)
	ax2 = plt.subplot2grid((7,14),(0,7),colspan = 7,rowspan = 5)
	ax10 = plt.subplot2grid((7,14),(5,0),colspan = 7,rowspan = 2)
	ax20 = plt.subplot2grid((7,14),(5,7),colspan = 7,rowspan = 2)

	AHF(phy,'GX',ax1,ax2,ax10,ax20)
	AHF(phy,'GM',ax1,ax2,ax10,ax20)
	SAM(phy,'Galacticus',ax1,ax2,ax10,ax20)
	SAM(phy,'SAG',ax1,ax2,ax10,ax20)
	SAM(phy,'SAGE',ax1,ax2,ax10,ax20)
	SDSS7(phy,ax1,ax2,ax10,ax20)
	ax1.get_legend_handles_labels()
	ax1.legend(loc = 'best',ncol = 2)
	ax2.get_legend_handles_labels()
	ax2.legend(loc = 'best',ncol = 2)
	ax10.set_ylim(-0.5,0.5)
	ax20.set_ylim(-0.5,0.5)
	ax10.axes.set_yticks([-0.4,-0.2,0,0.2,0.4])
	ax20.axes.set_yticks([-0.4,-0.2,0,0.2,0.4])
	ax1.set_ylim(4,8)
	ax2.set_ylim(4,8)
	ax1.set_ylabel(r'$Log_{10}\ \rho_{\bigstar}\ [M_⊙kpc^{-2}]$',fontsize = 14)
	ax1.text(0.06,4.4,r'$S\acute{e}rsic\ profile\ fit$', fontsize = 15,fontweight = 30)
	ax2.text(0.06,4.4,r'$NFW\ profile\ fit$', fontsize = 15,fontweight = 30)
	ax10.set_ylabel(r'$\rho_{\bigstar}/\rho_{fit}-1$',fontsize = 14)
	ax10.set_xlabel('$r_p$/$r_{200}$',fontsize = 14)
	ax20.set_xlabel('$r_p$/$r_{200}$',fontsize = 14)
		# ax10.axes.set_yticks([-0.1,0,0.1])


	ax1.tick_params(top = 'on', right = 'on',which='both',labelsize = 11)
	ax2.tick_params(top = 'on',right = 'on', left = 'off',which='both',labelsize = 11)

	ax10.tick_params(top = 'off', right = 'on',which='both',labelsize = 11)
	ax20.tick_params(top = 'off',right = 'on', left = 'off',which='both',labelsize = 11)

	ax1.set_xlim(0,1)
	ax2.set_xlim(0,1)

	ax1.axes.xaxis.set_ticklabels([])
	ax2.axes.xaxis.set_ticklabels([])
	ax2.axes.yaxis.set_ticklabels([])
	ax20.axes.yaxis.set_ticklabels([])


	ax10.axhline(0,c='grey',zorder = 30)
	ax10.set_xlim(0,1)
	ax10.axes.set_xticks([0,0.2,0.4,0.6,0.8])
	

	# ax20.plot(plt.xlim(),'k--',[0,0])
	ax20.axhline(0,c='grey',zorder = 30)
	ax20.set_xlim(0,1)
	
	plt.tight_layout()
	plt.subplots_adjust(wspace =0 ,hspace = 0)
	
	plt.savefig('/home/qyli/Desktop/stellar_density.pdf')
	plt.show()
main(phy = 'star')
