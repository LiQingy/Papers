import numpy as np 
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

global pi
pi = 3.1415926

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

		plt.errorbar(data3[:,0], data3[:,1], yerr=[data4, data2], xerr=[data1, data5], fmt = marker,markerfacecolor = 'none',markeredgewidth = 1.7,zorder = 26, markersize = 7, elinewidth = 2,label = name, color = colorr)

	if npoint == 2:
		ax.plot(dataM[:,0],dataM[:,1], color = colorr, linestyle = '--', label = name,zorder = 25)

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
		plt.errorbar(datax,datay,xerr = [xerr0,xerr1],label = name,markerfacecolor = 'none',color = colorr,fmt = marker,elinewidth = 2.5,markersize = 6)

def density0(ax1):
	#A. Vikhlinin
	marker = ','
	A478n = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A478n.csv'
	dataM = np.loadtxt(open(A478n, 'r'),delimiter = ',')
	for i in range(dataM.shape[0]):
		dataM[i,0] = dataM[i,0] / 1359 #x: R/R500
		dataM[i,1] = math.log10(dataM[i,1] / (10**3)**3 ) #y: MsunMpc-3 -> Msunkpc-3h2		
	pltob(ax1,'A478',2,'purple',marker,dataM)

	marker = ','
	A1795n = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A1795n.csv'
	dataM = np.loadtxt(open(A1795n, 'r'),delimiter = ',')
	for i in range(dataM.shape[0]):
		dataM[i,0] = dataM[i,0] / 1283
		dataM[i,1] = math.log10(dataM[i,1] / (10**3)**3 ) #MsunMpc-3 -> Msunkpc-3h2
	pltob(ax1,'A1795',2,'orange',marker,dataM)

	marker = ','
	A2029n = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A2029/A2029n.csv'
	dataM = np.loadtxt(open(A2029n, 'r'),delimiter = ',')
	for i in range(dataM.shape[0]):
		dataM[i,0] = dataM[i,0] / 1380
		dataM[i,1] = math.log10(dataM[i,1] / (10**3)**3 ) #MsunMpc-3 -> Msunkpc-3h2
	pltob(ax1,'A2029',2,'green',marker,dataM)
	#---
	pcri = math.log10(2.775 * 1e11 / 1e9 * 0.7**2) #Msunkpc-3h2
	x = np.array([0.01,0.02,0.04,0.07,0.11,0.16,0.22,0.29,0.37,0.45,0.54,0.64,0.74,0.84,0.94,1.04])
	y = np.array([3.57,3.25,3.13,2.96,2.80,2.66,2.52,2.37,2.21,2.06,1.92,1.79,1.67,1.55,1.43,1.32])
	y = y + pcri

	yer = np.array([0.14,0.11,0.09,0.08,0.06,0.04,0.03,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.02])
	ax1.plot(x,y,color = 'black',label = 'McDonald + 2017',linestyle = '-',linewidth = 2,zorder = 30)
	ax1.fill_between(x,y-yer,y+yer,color = 'black',alpha = 0.3,zorder = 29)


def median(n,phy,ax1,ax10,data,datax,mark,col,linet,lab1,lab2):
	if mark == 'o':
		filecen = '/Users/liqy/Documents/data/300data/center/G3X_Mass_snap_128-center-cluster.txt'
	if mark == '^':
		filecen = '/Users/liqy/Documents/data/300data/center/Music_Mass_snap_017-center-cluster.txt'
	datacen = np.loadtxt(filecen)

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
				if phy == 'density':
					yy.extend([np.log10(data[k2][k1] * 0.678**2)])
				else:
					yy.extend([data[k2][k1]])
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
	sigmay = ydown + yup - 2 * y


	# if mark == '^':
	# 	ax1.errorbar(x,y,yerr = [y - ydown, yup - y],fmt = mark,markersize = 6,color = col,zorder = 3,capsize=4,label = lab1)
	# if mark == 'o':
	ax1.fill_between(x,ydown,yup,alpha = 0.25,color = col,zorder = 1)
	ax1.scatter(x,y,marker = mark,s = 6**2,color = col,zorder = 3,label = lab1)
	
	x0 = np.array(x)[7:]
	y0 = np.array(y)[7:]
	sigmay = sigmay[7:]
	xline = np.logspace(-2,0,100)


	def func1(x,a1,b1,a2,b2,c):
		return np.log10(a1 / (1 + (x/b1)**2)**(3*c/2) + a2 / (1 + (x/b2)**2)**(3*c/2))
	if col == 'red':
		# coeffs, matcov = curve_fit(func1, x0, y0,p0 = (5e4,0.4,5e4,0.4,0.08))
		coeffs, matcov = curve_fit(func1, x0, y0,p0 = (1e5,1,1e5,1,1)) # normal one
	elif col == 'blue':
		# coeffs, matcov = curve_fit(func1, x0, y0,p0 = (1e4,0.3,1e4,0.3,0.8))
		coeffs, matcov = curve_fit(func1, x0, y0,p0 = (1e5,0.3,1e5,0.3,0.8)) # normal one
	yline = func1(xline, coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4])
	yline1 = func1(x, coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4])
	print(coeffs)

	def func2(x,a,b):
		return np.log10(a / (x/b * (1 + x/b)**2))
	coeffs, matcov = curve_fit(func2, x0, y0)
	yline2 = func2(x, coeffs[0],coeffs[1])

	def func3(x,a,b,c):
		return np.log10(a * np.exp(-b * (pow(x / c, 1/4) - 1)))
	coeffs, matcov = curve_fit(func3, x0, y0)
	yline3 = func3(x, coeffs[0],coeffs[1],coeffs[2])

	def func4(x,a,b,c):
		return np.log10(a / (1 + (x/b)**2)**(3*c/2))
	coeffs, matcov = curve_fit(func4, x0, y0,p0 = (1e6,1,1))
	yline4 = func4(x, coeffs[0],coeffs[1],coeffs[2])

	ax1.plot(xline,yline,color = col,label = lab2)
	ax10.plot(x,10**yline1/10**y - 1,color = col,linestyle = '-',label = r'double $\rm \beta$-model for %s' %lab1)
	ax10.plot(x,10**yline4/10**y - 1,color = col,linestyle = '--',label = r'$\rm \beta$-model for %s' %lab1)
	ax10.plot(x,10**yline3/10**y - 1,color = col,linestyle = '-.',label = r'S$\rm \acute{e}$rsic profile for %s' %lab1)
	ax10.plot(x,10**yline2/10**y - 1,color = col,linestyle = ':',label = 'NFW profile for %s' %lab1)

##----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

def Simulation(SIMmodel,phy,ax1,ax10,lab1,lab2,col,linet):

		
	file = '/Users/liqy/Documents/data/300data/simulation/%s_snaplog_%s.txt' %(SIMmodel,phy)

			
	if SIMmodel == 'GX':
		filex = '/Users/liqy/Documents/data/300data/simulation/GX_snaplog_gasdenx.txt'
		mark = 'o'

	if SIMmodel == 'GM':
		filex = '/Users/liqy/Documents/data/300data/simulation/GM_snaplog_gasdenx.txt'
		mark = '^'
	data = np.loadtxt(file)
	if phy == 'metal':
		data = data / 0.0134
	datax = np.loadtxt(filex)

	median(20,phy,ax1,ax10,data,datax,mark,col,linet,lab1,lab2)

def main(phy):
	fig = plt.figure(figsize=(6, 7.6))

	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	plt.tick_params(top = 'on', right = 'on',which='both')
	

	ax1 = plt.subplot2grid((8,5),(0,0),colspan = 5,rowspan = 5)
	ax10 = plt.subplot2grid((8,5),(5,0),colspan = 5,rowspan = 3)

#---------------------------------
	
	Simulation('GX',phy,ax1,ax10,lab1 = 'Gadget-X',lab2 = 'fit to Gadget-X',col = 'red',linet = '-')
	Simulation('GM',phy,ax1,ax10,lab1 = 'Gadget-MUSIC',lab2 = 'fit to Gadget-MUSIC',col = 'blue',linet = '-')

#-----------------------------------

	density0(ax1)
	ax1.set_ylim(3.5,6.0)
	ax1.set_xlim(0.03,1)
	ax1.set_ylabel(r'$Log_{10}\ \rho_{gas}\ [M_⊙kpc^{-3}]$',fontsize = 14)
	ax1.set_xscale('log')
	ax1.axes.xaxis.set_ticklabels([])
	ax1.legend(loc = 1, fontsize = 'small')


	ax10.set_xlabel('R/$R_{500}$',fontsize = 14)
	ax10.tick_params(labelsize = 11)
	ax10.set_xlim(0.03, 1)
	ax10.set_xscale('log')
	ax10.axhline(0,c='grey',zorder = 0)
	ax10.set_ylabel(r'$\rho_{fit}/\rho_{gas}-1$',fontsize = 14)
	ax10.set_ylim(-0.5,0.5)
	ax10.axes.set_yticks([-0.4,-0.2,0,0.2,0.4])
	ax10.legend(loc = 'best', fontsize = 'x-small', ncol = 2)
	# # plt.yticks([0.03,0.04,0.05,0.06,0.07])

	# plt.axvline(0.2,c = 'grey',linestyle = '--')
	# plt.legend(fontsize = 'small',ncol = 2,loc = 'best')
	ax1.text(0.035,3.8,r'Double $\rm \beta$-model fit', fontsize = 15,fontweight = 30)
	plt.tight_layout()
	plt.subplots_adjust(wspace =0 ,hspace = 0)
	plt.savefig('/Users/liqy/Desktop/gas_density.pdf' )
	plt.show()

main(phy = 'density')
