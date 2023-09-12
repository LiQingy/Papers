import numpy as np 
import math
import matplotlib.pyplot as plt
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

		if name == 'Lovisari + 2018':
			plt.errorbar(data3[:,0], data3[:,1] * 0.0134/0.01886, yerr=[data4* 0.0134/0.01886, data2* 0.0134/0.01886], xerr=[data1, data5], fmt = marker,markerfacecolor = 'none',markeredgewidth = 1.7,zorder = 26, markersize = 7, elinewidth = 2,label = name, color = colorr)
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

def Temp():
	A478Tc = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A478Tc.csv'
	A478Tl = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A478Tl.csv'
	A478Tt = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A478Tt.csv'
	A478Tb = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A478Tb.csv'
	A478Tr = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A478Tr.csv'
	A478Tc = np.loadtxt(open(A478Tc, 'r'),delimiter = ',')
	A478Tl = np.loadtxt(open(A478Tl, 'r'),delimiter = ',')	
	A478Tt = np.loadtxt(open(A478Tt, 'r'),delimiter = ',')			
	A478Tb = np.loadtxt(open(A478Tb, 'r'),delimiter = ',')	
	A478Tr = np.loadtxt(open(A478Tr, 'r'),delimiter = ',')	

	t500 = 8.85 * (7.83e14 / (1e15*0.7/0.71) ) **(2/3) * (0.3*(1 + 0.0881)**3 + 0.7)**(1/3) * (0.6125/0.6)
	print(t500)

	xl = (A478Tc[:,0] - A478Tl[:,0]) / 1359
	xr = (A478Tr[:,0] - A478Tc[:,0]) / 1359
	yt = (A478Tt[:,1] - A478Tc[:,1]) / t500
	yb = (A478Tc[:,1] - A478Tb[:,1]) / t500
	plt.errorbar(A478Tc[:,0]/1359,A478Tc[:,1]/t500,xerr = [xl,xr],yerr = [yb,yt],color = 'purple',label = 'A478',fmt = ',',zorder = 40,elinewidth = 2)
	# temp_obfit(x = A478Tc[1:12,0]/1359,y=A478Tc[1:12,1]/t500,col = 'purple',lab2 = 'fit to A478',pp0 = (1.4,0.65,3.6e-3,15,0.5,0.6))
	#-----------------------------------------------------------------------
	A1795Tc = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A1795Tc.csv'
	A1795Tl = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A1795Tl.csv'
	A1795Tt = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A1795Tt.csv'
	A1795Tb = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A1795Tb.csv'
	A1795Tr = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A1795Tr.csv'
	A1795Tc = np.loadtxt(open(A1795Tc, 'r'),delimiter = ',')
	A1795Tl = np.loadtxt(open(A1795Tl, 'r'),delimiter = ',')	
	A1795Tt = np.loadtxt(open(A1795Tt, 'r'),delimiter = ',')			
	A1795Tb = np.loadtxt(open(A1795Tb, 'r'),delimiter = ',')	
	A1795Tr = np.loadtxt(open(A1795Tr, 'r'),delimiter = ',')	

	t500 = 8.85 * (6.57e14 / (1e15*0.7/0.71)) **(2/3) * (0.3*(1 + 0.0622)**3 + 0.7)**(1/3) * (0.6125/0.6)
	print(t500)

	xl = (A1795Tc[:,0] - A1795Tl[:,0]) / 1283
	xr = (A1795Tr[:,0] - A1795Tc[:,0]) / 1283
	yt = (A1795Tt[:,1] - A1795Tc[:,1]) / t500
	yb = (A1795Tc[:,1] - A1795Tb[:,1]) / t500
	plt.errorbar(A1795Tc[:,0]/1283,A1795Tc[:,1]/t500,xerr = [xl,xr],yerr = [yb,yt],color = 'orange',label = 'A1795',fmt = ',',zorder = 41,elinewidth = 2)
	# temp_obfit(A1795Tc[:,0]/1359,A1795Tc[:,1]/t500,col = 'orange',lab2 = 'fit to A1795',pp0 = (1.4,0.65,3.6e-3,15,0.5,0.6))
	#-----------------------------------------------------------------------

	A2029Tc = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A2029/A2029Tc.csv'
	A2029Tl = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A2029/A2029Tl.csv'
	A2029Tt = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A2029/A2029Tt.csv'
	A2029Tb = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A2029/A2029Tb.csv'
	A2029Tr = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A2029/A2029Tr.csv'
	A2029Tc = np.loadtxt(open(A2029Tc, 'r'),delimiter = ',')
	A2029Tl = np.loadtxt(open(A2029Tl, 'r'),delimiter = ',')	
	A2029Tt = np.loadtxt(open(A2029Tt, 'r'),delimiter = ',')			
	A2029Tb = np.loadtxt(open(A2029Tb, 'r'),delimiter = ',')	
	A2029Tr = np.loadtxt(open(A2029Tr, 'r'),delimiter = ',')	

	t500 = 8.85 * (8.29e14 / (1e15*0.7/0.71)) **(2/3) * (0.3*(1 + 0.0779)**3 + 0.7)**(1/3) * (0.6125/0.6)
	print(t500)
	print(A2029Tc.shape)
	xl = (A2029Tc[:,0] - A2029Tl[:,0]) / 1380
	xr = (A2029Tr[:,0] - A2029Tc[:,0]) / 1380
	yt = (A2029Tt[:,1] - A2029Tc[:,1]) / t500
	yb = (A2029Tc[:,1] - A2029Tb[:,1]) / t500
	plt.errorbar(A2029Tc[:,0] / 1380,A2029Tc[:,1]/t500,xerr = [xl,xr],yerr = [yb,yt],color = 'green',label = 'A2029',fmt = ',',zorder = 42,elinewidth = 2)
	# temp_obfit(A2029Tc[:,0]/1359,A2029Tc[:,1]/t500,col = 'green',lab2 = 'fit to A478',pp0 = (1.4,0.65,3.6e-3,15,0.5,0.6))

	#------------------------------------------------------------------------------------------------
	# x = np.array([2.283e-2,2.909e-2,3.706e-2,4.723e-2,6.018e-2,7.669e-2,9.772e-2,1.245e-1,1.587e-1,2.022e-1,2.577e-1,3.283e-1,4.184e-1,5.331e-1,6.794e-1,8.657e-1])
	# y = np.array([7.606e-1,7.776e-1,8.259e-1,8.485e-1,8.759e-1,9.102e-1,9.460e-1,9.669e-1,9.717e-1,9.816e-1,9.751e-1,9.540e-1,9.317e-1,8.714e-1,7.529e-1,6.540e-1])


	# ylow = np.array([7.049e-1,7.338e-1,7.814e-1,8.086e-1,8.434e-1,8.705e-1,9.028e-1,9.314e-1,9.425e-1,9.354e-1,9.518e-1,9.347e-1,8.946e-1,8.409e-1,7.093e-1,6.061e-1])
	# yhigh = np.array([8.741e-1,8.823e-1,9.534e-1,9.598e-1,9.510e-1,9.650e-1,9.913e-1,1.011,1.021,1.024,1.005,9.794e-1,9.662e-1,8.970e-1,7.937e-1,6.917e-1])
	# plt.plot(x,y,color = 'black',label = 'Ghirardini + 2019',linestyle = '-',linewidth = 2,zorder = 0)
	# plt.fill_between(x,ylow,yhigh,color = 'black',alpha = 0.3,zorder = 29)
	#------------------------------------------------------------------------------------------------
	def fvikh(x,A0,A1,A2,A3,A4,A5):
		return A0*(A3+(x/A1)**A4)/(1.+(x/A1)**A4)/(1.+(x/A2)**2.)**A5

	data_multi2=np.loadtxt('/Users/liqy/Desktop/T_bfit.dat')
	x=np.exp(np.linspace(np.log(0.01),np.log(2),200))
	t_fitv=np.zeros((len(x),3))
	for i in range(len(x)):
		t_fitv[i,:]=np.percentile(fvikh(x[i],data_multi2[:,0],np.exp(data_multi2[:,1]),data_multi2[:,2],data_multi2[:,3],data_multi2[:,4],data_multi2[:,5]),(16,50,84))

	plt.plot(x,t_fitv[:,1],color='k',label = 'Ghirardini + 2019')
	plt.fill_between(x,t_fitv[:,0],t_fitv[:,2],color='k',alpha=0.3)

	# def ft(x,t0,tmint0,rcool,acool,rt,c2):
	# 	return t0 * (tmint0 + (x/rcool)**acool) / (1 + (x/rcool)**acool) / (1 + (x/rt)**2)**c2

	# t0 = 1.21
	# tmint0 = 0.5
	# rcool = np.e**(-2.8)
	# acool = 1.03
	# rt = 0.34
	# c2 = 0.27
	# xft = np.logspace(-2,0,100)
	# yft = ft(xft,t0,tmint0,rcool,acool,rt,c2)
	# plt.plot(xft,yft,color = 'black',label = 'Ghirardini + 2019',zorder = 0,linewidth = 2)

	# t0 = 1.21 - 0.23
	# tmint0 = 0.5 - 0.24
	# rcool = np.e**(-2.8 - 1.1)
	# acool = 1.03 - 0.78
	# rt = 0.34 - 0.1
	# c2 = 0.27 - 0.04
	# xlow = np.logspace(-2,0,100)
	# ylow = ft(xft,t0,tmint0,rcool,acool,rt,c2)

	# t0 = 1.21 + 0.23
	# tmint0 = 0.5 + 0.24
	# rcool = np.e**(-2.8 + 1.1)
	# acool = 1.03 + 0.78
	# rt = 0.34 + 0.1
	# c2 = 0.27 + 0.04
	# xup = np.logspace(-2,0,100)
	# yup = ft(xft,t0,tmint0,rcool,acool,rt,c2)

	# plt.fill_between(xft,ylow,yup,color = 'black',alpha = 0.3,zorder = 1)

def median(n,phy,data,datax,mark,col,linet,lab1,lab2,maindata):
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

	if phy == 'retemp':
		revalue = 1.14
	else:
		revalue = 1
	x = np.array(x) 
	y = np.array(y) * 1.14 * 0.6125/0.588
	ydown = np.array(ydown) * 1.14 * 0.6125/0.588
	yup = np.array(yup) * 1.14 * 0.6125/0.588

	if maindata == True:
		plt.fill_between(x,ydown,yup,alpha = 0.25,color = col,zorder = 1)
		plt.scatter(x,y,marker = mark,s = 6**2,color = col,zorder = 3,label = lab1)
	else:
		plt.plot(x,y * 1.14 * 1.14 * 0.6125/0.588,ls = '--',c = col, linewidth = 1.5,zorder = 5, label = r'%s on $\rm T_{500,G+19}$' %lab1)

	#---------------------------------------------------------------------------------------------
	#---------------------------------------------------------------------------------------------
	#fitting
	
	if maindata == True:
		x0 = np.array(x)[9:]
		y0 = np.array(y)[9:]
		xline = np.logspace(-2,0,100)

		def func(x,t0,tmin,rcool,acool,rt,c):
			a0 = (tmin / t0 + (x / rcool)**acool)
			a1 = (1 + (x/rcool)**acool)
			b1 = (1 + (x/rt)**2)**(c/2)
			return t0 * a0 / a1 / b1
			# return t0 * (tmin + (x / rcool)**acool) / (1 + (x/rcool)**acool) / (1 + (x/rt)**2)**(c)
		if lab1 == 'Gadget-X':
			if phy == 'retemp':
				coeffs, matcov = curve_fit(func, x0, y0, p0 = (0.86,0.64,0.24,2,0.34,0.55)) #rescale 
			else:
				coeffs, matcov = curve_fit(func, x0, y0, p0 = (0.7,0.66,0.02,4,0.34,0.55)) # not rescale 
		if lab1 == 'Gadget-MUSIC':
			if phy == 'retemp':
				coeffs, matcov = curve_fit(func, x0, y0, p0 = (0.86,0.44,0.02,4,0.34,0.55)) # rescale
			else:
				coeffs, matcov = curve_fit(func, x0, y0, p0 = (0.7,0.66,0.02,4,0.34,0.55)) # not rescale
			# coeffs, matcov = curve_fit(func, x0, y0, p0 =  (0.86,0.66,0.04,4,0.44,0.54))
		# coeffs, matcov = curve_fit(func, x0, y0, p0 = (1.2982171, 0.75643626, 0.01101951, 5.95338365, 0.13094973, 0.24326953))
		yline = func(xline, coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4],coeffs[5])
		print(coeffs)
		plt.plot(xline, yline, color = col, linestyle = linet,linewidth = 1.5,zorder = 2,label = lab2)

##----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

def Simulation(SIMmodel,phy,lab1,lab2,col,linet):

	file = '/Users/liqy/Documents/data/300data/simulation/%s_snaplog_spectemp.txt' %(SIMmodel)			
	if SIMmodel == 'GX':
		fileslT500 = '/Users/liqy/Documents/data/300data/simulation/GX_snaplog_slT500.txt'
		filex = '/Users/liqy/Documents/data/300data/simulation/GX_snaplog_gasx.txt'
		filecen = '/Users/liqy/Documents/data/oPDF/cluster/G3X_Mass_snap_128-center-cluster.txt'
		mark = 'o'

	if SIMmodel == 'GM':
		fileslT500 = '/Users/liqy/Documents/data/300data/simulation/GM_snaplog_slT500.txt'
		filex = '/Users/liqy/Documents/data/300data/simulation/GM_snaplog_gasx.txt'
		filecen = '/Users/liqy/Documents/data/oPDF/cluster/Music_Mass_snap_017-center-cluster.txt'
		mark = '^'
	slT500 = np.loadtxt(fileslT500)
	datacen = np.loadtxt(filecen)
	t500 = 8.85 * (datacen[:, 9] * 1e10 / 0.7 / 1e15) **(2/3) * (0.6125/0.6)
	datax = np.loadtxt(filex)

	data = np.loadtxt(file) / t500.reshape(324,1)
	median(20,phy,data,datax,mark,col,linet,lab1,lab2, maindata = True)
	# data = np.loadtxt(file) / t500.reshape(324,1)
	# median(20,phy,data,datax,mark,col,linet,lab1,lab2, maindata = False)

def main(phy):

	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	plt.tick_params(top = 'on', right = 'on',which='both')

#---------------------------------
	
	Simulation('GX',phy,lab1 = 'Gadget-X',lab2 = 'fit to Gadget-X',col = 'red',linet = '-')
	Simulation('GM',phy,lab1 = 'Gadget-MUSIC',lab2 = 'fit to Gadget-MUSIC',col = 'blue',linet = '-')

#-----------------------------------

	Temp()
	plt.ylim(0.4,1.5)
	plt.yscale('log')
	plt.ylabel(r'$T/T_{500}$',fontsize = 14)
	y_major = LogLocator(base = 10,subs = [1])
	y_minor = LogLocator(base = 10,subs = [0.4,0.5,0.6,0.7,0.8,0.9,2])
	ax = plt.gca()
	ax.yaxis.set_major_locator(y_major)
	ax.yaxis.set_minor_locator(y_minor)

	ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
	ax.yaxis.set_minor_formatter(ScalarFormatter())
	# ax.yaxis.set_minor_formatter(NullFormatter())
	# plt.yticks([0.4,0.5,0.6,0.7,0.8])

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

main(phy = 'retemp')
