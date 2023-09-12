import numpy as np 
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.optimize import curve_fit

from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import LogLocator
from matplotlib.ticker import NullFormatter
import matplotlib.ticker as mtick

global pi
pi = 3.1415926

def pltob(name,npoint,colorr,marker,dataM):
	
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

		if name == 'Lovisari + 2019':
			# for kk in range(data3.shape[0]):
			# 	xy = (data3[kk,0] - data1[kk],data3[kk,1]+data4[kk])
			# 	wd = data5[kk] + data1[kk]
			# 	ht = abs(data2[kk] + data4[kk])
			# 	if kk == 0:
			# 		rect = Rectangle(xy, width = wd, height = ht,color = 'cyan',edgecolor = '',alpha = 0.3)
			# 	else:
			# 		rect = Rectangle(xy, width = wd, height = ht,color = 'cyan',edgecolor = '',alpha = 0.3)
			# 	plt.gca().add_patch(rect)
			plt.fill_between(data3[:,0],data3[:,1] - data4, data3[:,1] + data2, color = 'cyan', alpha = 0.3)
			plt.plot(data3[:,0], data3[:,1], color = 'magenta', marker = 's', label = 'Lovisari + 2019', zorder = 100)
				
			# plt.errorbar(data3[:,0], data3[:,1], yerr=[data4, data2], xerr=[data1, data5], fmt = marker,markerfacecolor = 'none',markeredgewidth = 1.7,zorder = 26, markersize = 7, elinewidth = 2,label = name, color = colorr)
		else:
			plt.errorbar(data3[:,0], data3[:,1], yerr=[data4, data2], xerr=[data1, data5], fmt = marker,markerfacecolor = 'none',markeredgewidth = 1.7,zorder = 26, markersize = 7, elinewidth = 2,label = name, color = colorr)

	if npoint == 2:
		plt.plot(dataM[:,0],dataM[:,1], color = colorr, linestyle = '--', label = name,zorder = 25)

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


def GasZ():

	A478Zc = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A478Zc.csv'
	A478Zl = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A478Zl.csv'
	A478Zt = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A478Zt.csv'
	A478Zb = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A478Zb.csv'
	A478Zr = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A478Zr.csv'
	A478Zc = np.loadtxt(open(A478Zc, 'r'),delimiter = ',')
	A478Zl = np.loadtxt(open(A478Zl, 'r'),delimiter = ',')	
	A478Zt = np.loadtxt(open(A478Zt, 'r'),delimiter = ',')			
	A478Zb = np.loadtxt(open(A478Zb, 'r'),delimiter = ',')	
	A478Zr = np.loadtxt(open(A478Zr, 'r'),delimiter = ',')	

	xl = (A478Zc[:,0] - A478Zl[:,0]) / 1359
	xr = (A478Zr[:,0] - A478Zc[:,0]) / 1359
	yt = A478Zt[:,1] - A478Zc[:,1]
	yb = A478Zc[:,1] - A478Zb[:,1]
	plt.errorbar(A478Zc[:,0]/1359,A478Zc[:,1]* 0.01886 / 0.0134,xerr = [xl,xr],yerr = [yb* 0.01886 / 0.0134,yt* 0.01886 / 0.0134],elinewidth = 2,color = 'purple',label = 'A478',fmt = ',')

	A478Z00 = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A478Z00.csv'
	A478Z00 = np.loadtxt(open(A478Z00, 'r'),delimiter = ',')
	plt.plot(A478Z00[:,0]/1359,A478Z00[:,1],marker = ',',color = 'purple',linewidth = 2)
	#-----------------------------------------------------------------------
	A1795Zc = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A1795Zc.csv'
	A1795Zl = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A1795Zl.csv'
	A1795Zt = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A1795Zt.csv'
	A1795Zb = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A1795Zb.csv'
	A1795Zr = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A1795Zr.csv'
	A1795Zc = np.loadtxt(open(A1795Zc, 'r'),delimiter = ',')
	A1795Zl = np.loadtxt(open(A1795Zl, 'r'),delimiter = ',')	
	A1795Zt = np.loadtxt(open(A1795Zt, 'r'),delimiter = ',')			
	A1795Zb = np.loadtxt(open(A1795Zb, 'r'),delimiter = ',')	
	A1795Zr = np.loadtxt(open(A1795Zr, 'r'),delimiter = ',')	

	xl = (A1795Zc[:,0] - A1795Zl[:,0]) / 1283
	xr = (A1795Zr[:,0] - A1795Zc[:,0]) / 1283
	yt = A1795Zt[:,1] - A1795Zc[:,1]
	yb = A1795Zc[:,1] - A1795Zb[:,1]
	plt.errorbar(A1795Zc[:,0]/1283,A1795Zc[:,1]* 0.01886 / 0.0134,xerr = [xl,xr],yerr = [yb* 0.01886 / 0.0134,yt* 0.01886 / 0.0134],color = 'orange',label = 'A1795',elinewidth = 2,fmt = ',')

	A1795Z00 = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A1795Z00.csv'
	A1795Z00 = np.loadtxt(open(A1795Z00, 'r'),delimiter = ',')
	plt.plot(A1795Z00[:,0]/1283,A1795Z00[:,1],marker = ',',color = 'orange',linewidth = 2)
	#-----------------------------------------------------------------------
	A2029Zc = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A2029/A2029Zc.csv'
	A2029Zl = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A2029/A2029Zl.csv'
	A2029Zt = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A2029/A2029Zt.csv'
	A2029Zb = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A2029/A2029Zb.csv'
	A2029Zr = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A2029/A2029Zr.csv'
	A2029Zc = np.loadtxt(open(A2029Zc, 'r'),delimiter = ',')
	A2029Zl = np.loadtxt(open(A2029Zl, 'r'),delimiter = ',')	
	A2029Zt = np.loadtxt(open(A2029Zt, 'r'),delimiter = ',')			
	A2029Zb = np.loadtxt(open(A2029Zb, 'r'),delimiter = ',')	
	A2029Zr = np.loadtxt(open(A2029Zr, 'r'),delimiter = ',')	
	xl = (A2029Zc[:,0] - A2029Zl[:,0]) / 1380
	xr = (A2029Zr[:,0] - A2029Zc[:,0]) / 1380
	yt = A2029Zt[:,1] - A2029Zc[:,1]
	yb = A2029Zc[:,1] - A2029Zb[:,1]
	plt.errorbar(A2029Zc[:,0] / 1380, A2029Zc[:,1] * 0.01886 / 0.0134, xerr = [xl,xr],yerr = [yb* 0.01886 / 0.0134,yt* 0.01886 / 0.0134],color = 'green',label = 'A2029',elinewidth = 2,fmt = ',')

	A2029Z00 = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A2029/A2029Z00.csv'
	A2029Z00 = np.loadtxt(open(A2029Z00, 'r'),delimiter = ',')
	plt.plot(A2029Z00[:,0]/1380, A2029Z00[:,1], marker = ',',color = 'green',linewidth = 2)

	#S. Majerowicz
	marker = 'o'
	Lov = '/Users/liqy/Documents/data/300data/observation/Lovisari/1811.05987.csv'
	Lov = np.loadtxt(open(Lov, 'r'),delimiter = ',')
	pltob('Lovisari + 2019',5,'cyan',marker,Lov)

def median(n,phy,data,datax,mark,col,linet,lab1,lab2):
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

	# if mark == '^':
	# 	plt.errorbar(x,y,yerr = [y - ydown, yup - y],fmt = mark,markersize = 6,color = col,zorder = 3,capsize=4,label = lab1)
	# if mark == 'o':
	plt.fill_between(x,ydown,yup,alpha = 0.25,color = col,zorder = 1)
	plt.scatter(x,y,marker = mark,s = 6**2,color = col,zorder = 3,label = lab1)
	

	xline = np.logspace(-2,1,100)
	x0 = x[7:]
	y0 = y[7:]
	def func1(x,a,b,c):
		return a * np.exp(-b * (pow(x / c, 1/4) - 1))
	coeffs, matcov = curve_fit(func1, x0, y0)
	print(coeffs)
	yline = func1(xline, coeffs[0],coeffs[1],coeffs[2])
	plt.plot(xline,yline,color = col,label = lab2)

	#---------------------------------------------------------------------------------------------
	#---------------------------------------------------------------------------------------------
	#fitting
	
##----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

def Simulation(SIMmodel,phy,lab1,lab2,col,linet):

		
	file = '/Users/liqy/Documents/data/300data/simulation/%s_snaplog_%s.txt' %(SIMmodel,phy)

			
	if SIMmodel == 'GX':
		filex = '/Users/liqy/Documents/data/300data/simulation/GX_snaplog_gasx.txt'
		mark = 'o'

	if SIMmodel == 'GM':
		filex = '/Users/liqy/Documents/data/300data/simulation/GM_snaplog_gasx.txt'
		mark = '^'
	data = np.loadtxt(file)
	if phy == 'metal':
		data = data / 0.0134

	datax = np.loadtxt(filex)

	median(20,phy,data,datax,mark,col,linet,lab1,lab2)

def main(phy):

	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	plt.tick_params(top = 'on', right = 'on',which='both')

#---------------------------------
	
	Simulation('GX',phy,lab1 = 'Gadget-X',lab2 = 'fit to Gadget-X',col = 'red',linet = '-')
	Simulation('GM',phy,lab1 = 'Gadget-MUSIC',lab2 = 'fit to Gadget-MUSIC',col = 'blue',linet = '-')

#-----------------------------------
	if phy == 'metal':
		GasZ()
		plt.ylim(0,1)
		plt.ylabel(r'$Z_{gas}/Z_{\odot}$',fontsize = 14)

	plt.xlabel('R/$R_{500}$',fontsize = 14)
	# plt.tick_params(labelsize = 11)
	plt.xlim(0.03,1)

	plt.xscale('log')
	plt.legend(ncol = 2)
	

	# plt.axvline(0.2,c = 'grey',linestyle = '--')
	# plt.legend(fontsize = 'small',ncol = 2,loc = 'best')
	plt.tight_layout()
	plt.text(0.06,3.4,r'$Sersic$', fontsize = 15,fontweight = 30)
	plt.savefig('/Users/liqy/Desktop/gas_%s.pdf' %phy)
	plt.show()

main(phy = 'metal')
