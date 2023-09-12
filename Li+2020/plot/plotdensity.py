import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize
global x
global y
def pltob(ax,name,npoint,colorr,marker,dataM):
	
	
	sizeM = dataM.shape
	if npoint == 5:

		Mn = int(sizeM[0]/5)
		data1 = np.zeros(Mn)
		data2 = np.zeros(Mn)
		data3 = np.zeros(shape = (Mn,2))
		data4 = np.zeros(Mn)
		data5 = np.zeros(Mn)

		i = 0
		j = 0
		while i < sizeM[0]:
			datavol = np.sort(dataM[i+1:i+4], axis = 0)
			data3[j][0] = datavol[1][0]
			data3[j][1] = datavol[1][1]

			data1[j] = data3[j][0] - dataM[i][0]
			data2[j] = datavol[0][1] - data3[j][1]
			data4[j] = data3[j][1] - datavol[2][1]
			data5[j] = dataM[i+4][0] - data3[j][0]
			i += 5
			j += 1
		# plt.ylim(0,0.9)

		ax.errorbar(data3[:,0], data3[:,1], yerr=[data4, data2], xerr=[data1, data5], fmt = marker,markerfacecolor = 'none',zorder = 26, markersize = 6, elinewidth = 1.3,label = name, color = colorr)

	if npoint == 2:
		ax.plot(dataM[:,0],dataM[:,1], color = colorr, linestyle = '--', linewidth = 2,label = name,zorder = 25)

	if npoint == 3:
		Mn = int(sizeM[0]/3)
		datax = np.zeros(Mn)
		datay = np.zeros(Mn)
		xerr0 = np.zeros(Mn)
		xerr1 = np.zeros(Mn)
		i = 0
		j = 0
		while i < sizeM[0]:
			datavol = np.sort(dataM[i:i+3], axis = 0)
			datax[j] = datavol[1][0]
			datay[j] = datavol[1][1]
			xerr0[j] = datavol[1][0] - datavol[0][0]
			xerr1[j] = datavol[2][0] - datavol[1][0]
			i += 3
			j += 1
		ax.errorbar(datax,datay,xerr = [xerr0,xerr1],label = name,markerfacecolor = 'none',color = colorr,fmt = marker,markersize = 6)

def density0(ax):
	#A. Vikhlinin
	marker = ','
	A478n = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A478n.csv'
	dataM = np.loadtxt(open(A478n, 'r'),delimiter = ',')
	for i in range(dataM.shape[0]):
		dataM[i,0] = dataM[i,0] / 1359 #x: R/R500
		dataM[i,1] = math.log10(dataM[i,1] / (10**3)**3 ) #y: MsunMpc-3 -> Msunkpc-3	
	pltob(ax,'A478',2,'purple',marker,dataM)

	marker = ','
	A1795n = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A1795n.csv'
	dataM = np.loadtxt(open(A1795n, 'r'),delimiter = ',')
	for i in range(dataM.shape[0]):
		dataM[i,0] = dataM[i,0] / 1283
		dataM[i,1] = math.log10(dataM[i,1] / (10**3)**3 ) #MsunMpc-3 -> Msunkpc-3
	pltob(ax,'A1795',2,'orange',marker,dataM)

	marker = ','
	A2029n = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A2029/A2029n.csv'
	dataM = np.loadtxt(open(A2029n, 'r'),delimiter = ',')
	for i in range(dataM.shape[0]):
		dataM[i,0] = dataM[i,0] / 1380
		dataM[i,1] = math.log10(dataM[i,1] / (10**3)**3 ) #MsunMpc-3 -> Msunkpc-3
	pltob(ax,'A2029',2,'green',marker,dataM)
	#---
	pcri = math.log10(2.775 * 1e11 / 1e9 * 0.7**2) #Msunkpc-3
	x = np.array([0.01,0.02,0.04,0.07,0.11,0.16,0.22,0.29,0.37,0.45,0.54,0.64,0.74,0.84,0.94,1.04])
	y = np.array([3.57,3.25,3.13,2.96,2.80,2.66,2.52,2.37,2.21,2.06,1.92,1.79,1.67,1.55,1.43,1.32])
	y = y + pcri

	yer = np.array([0.14,0.11,0.09,0.08,0.06,0.04,0.03,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.02])
	ax.plot(x,y,color = 'black',label = 'McDonald + 2017',linestyle = '-',linewidth = 2,zorder = 30)
	ax.fill_between(x,y-yer,y+yer,color = 'black',alpha = 0.3,zorder = 29)

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

	x0 = np.array(x)[:]
	y0 = np.array(y)[:]
	#---------------------------------------------------------------------------------
	#'SEC'
	# xline = np.linspace(0,1,100)
	# def func1(x,a,b,c):
	# 	return np.log10(a * np.exp(-b * (pow(x / c, 1/4) - 1)))
	# coeffs, matcov = curve_fit(func1, x0, y0)
	# yline = func1(xline, coeffs[0],coeffs[1],coeffs[2])
	# # print(np.round(coeffs,decimals = 3))
	
	# #residuals
	# res1 = 10**y/10**func1(x, coeffs[0],coeffs[1],coeffs[2]) - 1
	# if lab1 == 'SDSS7':
	# 	ax1.plot(xline, yline, color = col, linestyle = linet,linewidth = 2.2,zorder = 2,label = lab2)
	# 	ax10.plot(x,res1,color = col, linestyle = linet,linewidth = 1.9)
	# else:
	# 	ax1.plot(xline, yline, color = col, linestyle = linet,linewidth = 1.5,zorder = 2,label = lab2)
	# 	ax10.plot(x,res1,color = col, linestyle = linet)

	# #-------------------------------------------------
	# # #'NFW'
	# def func2(x,a,b):
	# 	return np.log10(a / (x/b * (1 + x/b)**2))
	# coeffs, matcov = curve_fit(func2, x0, y0)
	# yline = func2(xline, coeffs[0],coeffs[1])
	# print(coeffs)
	# #residuals
	# res2 = 10**y/10**func2(x, coeffs[0],coeffs[1]) - 1
	# if lab1 == 'SDSS7':
	# 	ax2.plot(xline, yline, color = col, linestyle = linet,linewidth = 2.2,zorder = 2,label = lab2)
	# 	ax20.plot(x,res2,color = col, linestyle = linet,linewidth = 1.9)
	# else:
	# 	ax2.plot(xline, yline, color = col, linestyle = linet,linewidth = 1.5,zorder = 2,label = lab2)
	# 	ax20.plot(x,res2,color = col, linestyle = linet)
	#----------------------------------------------------------------------------------
	#Gas beta model
	xline = np.linspace(0,1,100)
	def func1(x,a,b,c):
		return np.log10(a / (1 + (x/b)**2)**(3*c/2))
	coeffs, matcov = curve_fit(func1, x0, y0,p0 = (1e6,1,1))
	# coeffs = [1e7,1,1]
	yline = func1(xline, coeffs[0],coeffs[1],coeffs[2])
	# print(np.round(coeffs,decimals = 3))
	print(coeffs)
	
	#residuals
	res1 = 10**y/10**func1(x, coeffs[0],coeffs[1],coeffs[2]) - 1
	if lab1 == 'SDSS7':
		ax1.plot(xline, yline, color = col, linestyle = linet,linewidth = 2.2,zorder = 2,label = lab2)
		ax10.plot(x,res1,color = col, linestyle = linet,linewidth = 1.9)
	else:
		ax1.plot(xline, yline, color = col, linestyle = linet,linewidth = 1.5,zorder = 2,label = lab2)
		ax10.plot(x,res1,color = col, linestyle = linet)

	#-------------------------------------------------
	# #Gas double beta model
	def func2(x,a1,b1,a2,b2,c):
		return np.log10(a1 / (1 + (x/b1)**2)**(3*c/2) + a2 / (1 + (x/b2)**2)**(3*c/2))
	coeffs, matcov = curve_fit(func2, x0, y0,p0 = (1e6,1,1e6,1,1))
	yline = func2(xline, coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4])
	print(coeffs)
	#residuals
	res2 = 10**y/10**func2(x, coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4]) - 1
	if lab1 == 'SDSS7':
		ax2.plot(xline, yline, color = col, linestyle = linet,linewidth = 2.2,zorder = 2,label = lab2)
		ax20.plot(x,res2,color = col, linestyle = linet,linewidth = 1.9)
	else:
		ax2.plot(xline, yline, color = col, linestyle = linet,linewidth = 1.5,zorder = 2,label = lab2)
		ax20.plot(x,res2,color = col, linestyle = linet)


#------------------------------------------------
#------------------------------------------------
#------------------------------------------------
def Simulation(phy,SIMmodel,ax1,ax2,ax10,ax20):

		
	file = '/Users/liqy/Documents/data/300data/simulation/%s500_snap_density0.txt' %(SIMmodel)
	data = np.loadtxt(file)

	if SIMmodel == 'GX':
		filex = '/Users/liqy/Documents/data/300data/simulation/GX500_snap_rr0.txt'
		col = 'red'
		mark = 'o'
		linet = '-'
		lab1 = 'Gadget-X'
		lab2 = 'fit to Gadget-X'
	
	if SIMmodel == 'GM':
		filex = '/Users/liqy/Documents/data/300data/simulation/GM500_snap_rr0.txt'
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
		file = '/Users/liqy/Documents/data/300data/AHF/AHFdata/%s500_AHF_densityall.txt' %(SIMmodel)
		filex = '/Users/liqy/Documents/data/300data/AHF/AHFdata/%s500_AHF_rrall.txt' %SIMmodel
		data = np.loadtxt(file) * 10**10 / 0.678
	else:
		file = '/Users/liqy/Documents/data/300data/AHF/AHFdata/%s_AHFxy_density4.txt' %(SIMmodel)
		filex = '/Users/liqy/Documents/data/300data/AHF/AHFdata/%s_AHFxy_rr4.txt' %SIMmodel
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
		file = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM500_MDPL_densityall.txt' 
		filex = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM500_MDPL_rrall.txt' 
	else:
		file = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_density4.txt' %(SAMmodel)
		filex = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_rr4.txt' %SAMmodel
	data = np.loadtxt(file)
	datax = np.loadtxt(filex)
	median(phy,ax1,ax2,ax10,ax20,20,data,datax,mark,col,linet,lab1,lab2)

def SDSS7(phy,ax1,ax2,ax10,ax20):

	
	filedata = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_density4PlanckBCG.txt' 
	filex = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_rr4PlanckBCG.txt'

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

	if phy == 'star' or phy == 'total':
		AHF(phy,'GX',ax1,ax2,ax10,ax20)
		AHF(phy,'GM',ax1,ax2,ax10,ax20)

		if phy == 'total':
			SAM(phy,'MDPL',ax1,ax2,ax10,ax20)
			ax1.get_legend_handles_labels()
			ax1.legend(loc = 'best',ncol = 2)
			ax2.get_legend_handles_labels()
			ax2.legend(loc = 'best',ncol = 2)
			ax1.set_ylim(3.5,7)
			ax2.set_ylim(3.5,7)
			ax10.set_ylim(-0.3,0.3)
			ax20.set_ylim(-0.3,0.3)
			ax1.set_ylabel(r'$Log_{10}\ \rho\ [M_⊙kpc^{-3}]$')
			ax10.set_ylabel(r'$\rho/\rho_{fit}-1$')
			ax1.text(0.06,4,r'$S\acute{e}rsic\ profile\ fit$', fontsize = 14,fontweight = 30)
			ax2.text(0.06,4,'NFW profile fit', fontsize = 14,fontweight = 30)
			ax1.text(0.6,6,'total density', fontsize = 14,fontweight = 30)
			ax2.text(0.6,6,'total density', fontsize = 14,fontweight = 30)
		else:
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

	if phy == 'gas':
		Simulation(phy,'GX',ax1,ax2,ax10,ax20)
		Simulation(phy,'GM',ax1,ax2,ax10,ax20)
		density0(ax1)
		density0(ax2)
		ax1.set_ylim(3,6)
		ax2.set_ylim(3,6)
		ax1.set_ylabel(r'$Log_{10}\ \rho_{gas}\ [M_⊙kpc^{-3}]$',fontsize = 14)
		ax1.text(0.06,3.4,r'$Beta\ model\ fit$', fontsize = 15,fontweight = 30)
		ax2.text(0.06,3.4,r'$Double\ beta\ model\ fit$', fontsize = 15,fontweight = 30)
		# ax1.text(0.06,3.4,r'$S\acute{e}rsic\ profile\ fit$', fontsize = 15,fontweight = 30)
		# ax2.text(0.06,3.4,r'$NFW\ profile\ fit$', fontsize = 15,fontweight = 30)
		ax1.get_legend_handles_labels()
		ax1.legend(loc = 'best')
		ax2.get_legend_handles_labels()
		ax2.legend(loc = 'best')
		ax10.set_ylabel(r'$\rho_{gas}/\rho_{fit}-1$',fontsize = 14)

		ax10.set_ylim(-0.3,0.3)
		ax20.set_ylim(-0.3,0.3)
		# ax10.axes.set_yticks([-0.08,-0.04,0,0.04,0.08])
		ax10.set_xlabel('R/$R_{500}$',fontsize = 14)
		ax20.set_xlabel('R/$R_{500}$',fontsize = 14)


	
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
	
	plt.savefig('/Users/liqy/Desktop/%s_density.pdf' %phy)
	plt.show()
main(phy = 'gas')
