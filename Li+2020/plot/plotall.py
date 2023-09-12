import numpy as np 
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

global pi
pi = 3.1415926

def nfw(xr,rv,mvir,cnfw):
	rs = rv / cnfw
	r = 0
	ps0 = rs**3 *(rs / (rs + r) + np.log(r + rs))
	r = rv
	psrvir = rs**3 *(rs / (rs + r) + np.log(r + rs))
	ps = mvir / 4 / pi / (psrvir - ps0) * 0.678**3 #Msun * kpc^-3
	print(rs,ps)
	y = ps / (xr / rs * (1 + xr/ rs)**2)
	return y

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



def density0():
	#A. Vikhlinin
	marker = ','
	A478n = '/home/qyli/Desktop/300data/observation/Vikhlinin/A478n.csv'
	dataM = np.loadtxt(open(A478n, 'r'),delimiter = ',')
	for i in range(dataM.shape[0]):
		dataM[i,0] = dataM[i,0] / 1359 #x: R/R500
		dataM[i,1] = math.log10(dataM[i,1] / (10**3)**3 ) #y: MsunMpc-3 -> Msunkpc-3h2		
	pltob(ax,'A478',2,'purple',marker,dataM)

	marker = ','
	A1795n = '/home/qyli/Desktop/300data/observation/Vikhlinin/A1795n.csv'
	dataM = np.loadtxt(open(A1795n, 'r'),delimiter = ',')
	for i in range(dataM.shape[0]):
		dataM[i,0] = dataM[i,0] / 1283
		dataM[i,1] = math.log10(dataM[i,1] / (10**3)**3 ) #MsunMpc-3 -> Msunkpc-3h2
	pltob(ax,'A1795',2,'yellow',marker,dataM)

	marker = ','
	A2029n = '/home/qyli/Desktop/300data/observation/Vikhlinin/A2029n.csv'
	dataM = np.loadtxt(open(A2029n, 'r'),delimiter = ',')
	for i in range(dataM.shape[0]):
		dataM[i,0] = dataM[i,0] / 1380
		dataM[i,1] = math.log10(dataM[i,1] / (10**3)**3 ) #MsunMpc-3 -> Msunkpc-3h2
	pltob(ax,'A2029',2,'green',marker,dataM)
	#---
	pcri = math.log10(2.775 * 1e11 / 1e9 * 0.7**2) #Msunkpc-3h2
	x = np.array([0.01,0.02,0.04,0.07,0.11,0.16,0.22,0.29,0.37,0.45,0.54,0.64,0.74,0.84,0.94,1.04])
	y = np.array([3.57,3.25,3.13,2.96,2.80,2.66,2.52,2.37,2.21,2.06,1.92,1.79,1.67,1.55,1.43,1.32])
	y = y + pcri

	yer = np.array([0.14,0.11,0.09,0.08,0.06,0.04,0.03,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.02])
	plt.plot(x,y,color = 'black',label = 'McDonald + 2017',linestyle = '-',zorder = 30)
	plt.fill_between(x,y-yer,y+yer,color = 'black',alpha = 0.3,zorder = 29)
	#---NFW
	xr= np.linspace(0,2242.26,100)
	rv = 2242.26
	mvir = 2.62134e15
	cnfw = 3.521
	y = nfw(xr,rv,mvir,cnfw)
	y = np.log10(y) - 0.75
	r500 = 1431.0203
	plt.plot(xr / r500, y, c = 'yellow', label = r'$10^{-0.75}\ NFW\ profile$',linewidth = 2, zorder = 30) 


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
	plt.errorbar(A478Zc[:,0]/1359,A478Zc[:,1],xerr = [xl,xr],yerr = [yb,yt],elinewidth = 2,color = 'purple',label = 'A478',fmt = ',')

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
	plt.errorbar(A1795Zc[:,0]/1283,A1795Zc[:,1],xerr = [xl,xr],yerr = [yb,yt],color = 'orange',label = 'A1795',elinewidth = 2,fmt = ',')

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
	plt.errorbar(A2029Zc[:,0] / 1380,A2029Zc[:,1],xerr = [xl,xr],yerr = [yb,yt],color = 'green',label = 'A2029',elinewidth = 2,fmt = ',')

	A2029Z00 = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A2029/A2029Z00.csv'
	A2029Z00 = np.loadtxt(open(A2029Z00, 'r'),delimiter = ',')
	plt.plot(A2029Z00[:,0]/1380,A2029Z00[:,1],marker = ',',color = 'green',linewidth = 2)

	#S. Majerowicz
	marker = 'o'
	Lov = '/Users/liqy/Documents/data/300data/observation/Lovisari/1811.05987.csv'
	Lov = np.loadtxt(open(Lov, 'r'),delimiter = ',')
	pltob('Lovisari + 2018',5,'cyan',marker,Lov)

	# x = np.array([0.02,0.06,0.1,0.15,0.21,0.28,0.385])
	# y = np.array([0.45,0.352,0.294,0.262,0.23,0.245,0.248])
	# yer = np.array([0.01,0.009,0.01,0.01,0.013,0.017,0.023])
	# xer = np.array([0.02,0.02,0.02,0.03,0.03,0.04,0.065])
	# plt.errorbar(x / 0.6, y, yerr = yer, xerr = xer / 0.6,fmt = 's',markerfacecolor = 'none',zorder = 26, markersize = 6, elinewidth = 1.3,label = 'Leccardi + 2008', color = 'black')
	# xup = '/home/qyli/Desktop/300data/observation/Leccardi/xup.csv'
	# xdown = '/home/qyli/Desktop/300data/observation/Leccardi/xdown.csv'
	# xup = np.loadtxt(open(xup, 'r'),delimiter = ',')
	# xdown = np.loadtxt(open(xdown, 'r'),delimiter = ',')
	# xave = (xup[:,0] + xdown[:,0]) / 2 / 0.6
	# plt.fill_between(xave,xdown[:,1],xup[:,1],alpha = 0.3,color = 'black',zorder = 1)

	# x = np.linspace(0.4/0.6,1)
	# plt.fill_between(x,0.06,0.246,alpha = 0.3,color = 'green',label = 'Molendi + 2016')
def temp_obfit(x,y,col,lab2,pp0):
	
	xline = np.linspace(0,1,100)
	# func = lambda cc: 

	def func2(x,t0,tmin,rcool,acool,rt,c):
		a0 = (tmin / t0 + (x / rcool)**acool)
		a1 = (1 + (x/rcool)**acool)
		b1 = (1 + (x/rt)**2)**(c/2)
		return t0 * a0 / a1 / b1
		# return t0 * (tmin + (x / rcool)**acool) / (1 + (x/rcool)**acool) / (1 + (x/rt)**2)**(c)

	coeffs, matcov = curve_fit(func2, x, y, p0 = pp0)
	yline = func2(xline, coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4],coeffs[5])
	print(coeffs)
	
	plt.plot(xline, yline, color = col,linestyle = '--',linewidth = 1.5,zorder = 2,label = lab2)
	
def Temp():
	A478Tc = '/home/qyli/Desktop/300data/observation/Vikhlinin/A478Tc.csv'
	A478Tl = '/home/qyli/Desktop/300data/observation/Vikhlinin/A478Tl.csv'
	A478Tt = '/home/qyli/Desktop/300data/observation/Vikhlinin/A478Tt.csv'
	A478Tb = '/home/qyli/Desktop/300data/observation/Vikhlinin/A478Tb.csv'
	A478Tr = '/home/qyli/Desktop/300data/observation/Vikhlinin/A478Tr.csv'
	A478Tc = np.loadtxt(open(A478Tc, 'r'),delimiter = ',')
	A478Tl = np.loadtxt(open(A478Tl, 'r'),delimiter = ',')	
	A478Tt = np.loadtxt(open(A478Tt, 'r'),delimiter = ',')			
	A478Tb = np.loadtxt(open(A478Tb, 'r'),delimiter = ',')	
	A478Tr = np.loadtxt(open(A478Tr, 'r'),delimiter = ',')	

	t500 = 8.85 * (7.83e14 / 1e15 * 0.7) **(2/3) * (0.3*(1 + 0.0881)**3 + 0.7) * (0.588/0.6)	

	xl = (A478Tc[:,0] - A478Tl[:,0]) / 1359
	xr = (A478Tr[:,0] - A478Tc[:,0]) / 1359
	yt = (A478Tt[:,1] - A478Tc[:,1]) / t500
	yb = (A478Tc[:,1] - A478Tb[:,1]) / t500
	plt.errorbar(A478Tc[:,0]/1359,A478Tc[:,1]/t500,xerr = [xl,xr],yerr = [yb,yt],color = 'purple',label = 'A478',fmt = ',',zorder = 40,elinewidth = 2)
	# temp_obfit(x = A478Tc[1:12,0]/1359,y=A478Tc[1:12,1]/t500,col = 'purple',lab2 = 'fit to A478',pp0 = (1.4,0.65,3.6e-3,15,0.5,0.6))
	#-----------------------------------------------------------------------
	A1795Tc = '/home/qyli/Desktop/300data/observation/Vikhlinin/A1795Tc.csv'
	A1795Tl = '/home/qyli/Desktop/300data/observation/Vikhlinin/A1795Tl.csv'
	A1795Tt = '/home/qyli/Desktop/300data/observation/Vikhlinin/A1795Tt.csv'
	A1795Tb = '/home/qyli/Desktop/300data/observation/Vikhlinin/A1795Tb.csv'
	A1795Tr = '/home/qyli/Desktop/300data/observation/Vikhlinin/A1795Tr.csv'
	A1795Tc = np.loadtxt(open(A1795Tc, 'r'),delimiter = ',')
	A1795Tl = np.loadtxt(open(A1795Tl, 'r'),delimiter = ',')	
	A1795Tt = np.loadtxt(open(A1795Tt, 'r'),delimiter = ',')			
	A1795Tb = np.loadtxt(open(A1795Tb, 'r'),delimiter = ',')	
	A1795Tr = np.loadtxt(open(A1795Tr, 'r'),delimiter = ',')	

	t500 = 8.85 * (6.57e14 / 1e15 * 0.7) **(2/3) * (0.3*(1 + 0.0622)**3 + 0.7) * (0.588/0.6)

	xl = (A1795Tc[:,0] - A1795Tl[:,0]) / 1283
	xr = (A1795Tr[:,0] - A1795Tc[:,0]) / 1283
	yt = (A1795Tt[:,1] - A1795Tc[:,1]) / t500
	yb = (A1795Tc[:,1] - A1795Tb[:,1]) / t500
	plt.errorbar(A1795Tc[:,0]/1283,A1795Tc[:,1]/t500,xerr = [xl,xr],yerr = [yb,yt],color = 'orange',label = 'A1795',fmt = ',',zorder = 41,elinewidth = 2)
	# temp_obfit(A1795Tc[:,0]/1359,A1795Tc[:,1]/t500,col = 'orange',lab2 = 'fit to A1795',pp0 = (1.4,0.65,3.6e-3,15,0.5,0.6))
	#-----------------------------------------------------------------------

	A2029Tc = '/Users/liqy/Documents/data/300data/observation/Vikhlinin/A2029/A2029Tc.csv'
	A2029Tl = '/home/qyli/Desktop/300data/observation/Vikhlinin/A2029/A2029Tl.csv'
	A2029Tt = '/home/qyli/Desktop/300data/observation/Vikhlinin/A2029/A2029Tt.csv'
	A2029Tb = '/home/qyli/Desktop/300data/observation/Vikhlinin/A2029/A2029Tb.csv'
	A2029Tr = '/home/qyli/Desktop/300data/observation/Vikhlinin/A2029/A2029Tr.csv'
	A2029Tc = np.loadtxt(open(A2029Tc, 'r'),delimiter = ',')
	A2029Tl = np.loadtxt(open(A2029Tl, 'r'),delimiter = ',')	
	A2029Tt = np.loadtxt(open(A2029Tt, 'r'),delimiter = ',')			
	A2029Tb = np.loadtxt(open(A2029Tb, 'r'),delimiter = ',')	
	A2029Tr = np.loadtxt(open(A2029Tr, 'r'),delimiter = ',')	

	t500 = 8.85 * (8.29e14 / 1e15 * 0.7) **(2/3) * (0.3*(1 + 0.0779)**3 + 0.7) * (0.588/0.6)
	print(A2029Tc.shape)
	xl = (A2029Tc[:,0] - A2029Tl[:,0]) / 1380
	xr = (A2029Tr[:,0] - A2029Tc[:,0]) / 1380
	yt = (A2029Tt[:,1] - A2029Tc[:,1]) / t500
	yb = (A2029Tc[:,1] - A2029Tb[:,1]) / t500
	plt.errorbar(A2029Tc[:,0] / 1380,A2029Tc[:,1]/t500,xerr = [xl,xr],yerr = [yb,yt],color = 'green',label = 'A2029',fmt = ',',zorder = 42,elinewidth = 2)
	# temp_obfit(A2029Tc[:,0]/1359,A2029Tc[:,1]/t500,col = 'green',lab2 = 'fit to A478',pp0 = (1.4,0.65,3.6e-3,15,0.5,0.6))

	def ft(x,t0,tmint0,rcool,acool,rt,c2):
		return t0 * (tmint0 + (x/rcool)**acool) / (1 + (x/rcool)**acool) / (1 + (x/rt)**2)**c2

	t0 = 1.21
	tmint0 = 0.5
	rcool = 10**-2.8
	acool = 1.03
	rt = 0.34
	c2 = 0.27
	xft = np.linspace(0,1,100)
	yft = ft(xft,t0,tmint0,rcool,acool,rt,c2)
	plt.plot(xft,yft,color = 'black',label = 'Ghirardini + 2019',zorder = 0,linewidth = 2)

	t0 = 1.21 - 0.23
	tmint0 = 0.5 - 0.24
	rcool = 10**(-2.8 - 1.1)
	acool = 1.03 - 0.78
	rt = 0.34 - 0.1
	c2 = 0.27 - 0.04
	xlow = np.linspace(0,1,100)
	ylow = ft(xft,t0,tmint0,rcool,acool,rt,c2)

	t0 = 1.21 + 0.23
	tmint0 = 0.5 + 0.24
	rcool = 10**(-2.8 + 1.1)
	acool = 1.03 + 0.78
	rt = 0.34 + 0.1
	c2 = 0.27 + 0.04
	xup = np.linspace(0,1,100)
	yup = ft(xft,t0,tmint0,rcool,acool,rt,c2)

	plt.fill_between(xft,ylow,yup,color = 'black',alpha = 0.3,zorder = 1)





	
#------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------

def plot(n,phy,data,datarv,loglog,lab,volumnrv):

	for i in range(data.shape[0]):
		x = []
		y = []
		for j in range(n):
			if math.isnan(data[i][j]) == False and data[i][j] != 0:
				rx = math.log10(datarv[i][volumnrv]) / n * (j + 1)
				x.extend([rx])
				if loglog == True:
					y.extend([math.log10(data[i][j])])
				else:
					y.extend([data[i][j]])
		x = np.array(x)
		y = np.array(y)
		if phy == 'gasz':
			coolr = 'skyblue'
		if phy == 'starz':
			coolr = 'cadetblue'
		if phy == 'age':
			coolr = 'lightsalmon'
		if phy == 'sfr':
			coolr = 'plum'
		if phy == 'ssfr':
			coolr = 'palegreen'
		if phy == 'density0':
			coolr = 'darkgrey'
		if phy == 'density4':
			coolr = 'tan'
		if phy == 'temp':
			coolr = 'gold'
		if i == 0:
			plt.plot(x,y,color = coolr,label = lab,alpha = 0.3, zorder = 1)
		else:
			plt.plot(x,y,color = coolr,alpha = 0.15, zorder = 1)


def median(n,phy,data,datax,loglog,allfit,mark,col,col1,linet,lab1,lab2):
	if lab1 == 'Gadget-X':
		filecen = '/Users/liqy/Documents/data/300data/center/G3X_Mass_snap_128-center-cluster.txt'
		datacen = np.loadtxt(filecen)
	if lab1 == 'Gadget-MUSIC':
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
				if phy == 'temp':
					m500 = datacen[k2][9] *1e10 / 0.678
					t500 = 8.85 * (m500 / 1e15 * 0.7) **(2/3) * (0.588/0.6) * (0.307 + 0.693)			
					ttem = data[k2][k1] / t500
					yy.extend([ttem])
				elif loglog == True:
					yy.extend([math.log10(data[k2][k1])])
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


	if allfit == True:
		if lab2 == 'fit to SDSS7':
			plt.errorbar(x,y,yerr = [y - ydown, yup - y],fmt = mark,markeredgewidth = 1.7,markersize = 7,color = col,zorder = 3,capsize=4,elinewidth = 1.7,label = lab1)
		else:
			plt.fill_between(x,ydown,yup,alpha = 0.25,color = col,zorder = 1)
			plt.scatter(x,y,marker = mark,s = 6**2,color = col,zorder = 3,label = lab1)

	#---------------------------------------------------------------------------------------------
	#---------------------------------------------------------------------------------------------
	#fitting
	if phy == 'starz' or phy == 'age' or phy == 'sfr' or phy == 'ssfr':

		x = x[1:]
		y = y[1:]

		xi2 = 0
		xiyi = 0
		xsum = 0
		ysum = 0
		nub = 0

		for i in range(x.shape[0]):
			xj = x[i]
			yy = y[i]
					
			xi2 = xj**2 + xi2
			xiyi = yy*xj + xiyi
			xsum = xj + xsum
			ysum = yy + ysum
			nub += 1
		
		xave = xsum / nub
		yave = ysum / nub
		b = (xiyi - nub*xave*yave) / (xi2 - nub*xave**2)
		a = yave - b * xave
		print(b,a)
		xline = np.linspace(0,3.5)
		yline = xline * b + a
		if lab1 == 'SDSS7':
			plt.plot(xline, yline, color = col, linestyle = linet,linewidth = 2.2,zorder = 2,label = lab2)
		else:	
			plt.plot(xline,yline,color = col1,linestyle = linet,linewidth = 1.5,zorder = 2,label = lab2)
	else:
		xline = np.linspace(0,1,500)
		if phy == 'gasz':
			x = np.array(x)[1:]
			y = np.array(y)[1:]
			def func(x,a,b,c):
				return a * np.exp(-b * (pow(x / c, 1/4) - 1))
			coeffs, matcov = curve_fit(func, x, y)
			yline = func(xline, coeffs[0],coeffs[1],coeffs[2])
			print(coeffs)

		if phy == 'temp':
			x = np.array(x)[:]
			y = np.array(y)[:]
			def func(x,t0,tmin,rcool,acool,rt,c):
				a0 = (tmin / t0 + (x / rcool)**acool)
				a1 = (1 + (x/rcool)**acool)
				b1 = (1 + (x/rt)**2)**(c/2)
				return t0 * a0 / a1 / b1
				# return t0 * (tmin + (x / rcool)**acool) / (1 + (x/rcool)**acool) / (1 + (x/rt)**2)**(c)
			coeffs, matcov = curve_fit(func, x, y, p0 = (1.21,1.19,1.6e-3,1.03,0.34,0.54))
			yline = func(xline, coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4],coeffs[5])
		print(coeffs)
		
		plt.plot(xline, yline, color = col1,linestyle = linet,linewidth = 1.5,zorder = 2,label = lab2)


##----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

def Simulation(SIMmodel,phy,loglog,allfit,lines,lab1,lab2):

		
	file = '/Users/liqy/Documents/data/300data/simulation/%s500_snap_%s.txt' %(SIMmodel,phy)

			
	if SIMmodel == 'GX':
		filex = '/Users/liqy/Documents/data/300data/simulation/GX500_snap_rr0.txt'
		filecste = '/Users/liqy/Documents/data/300data/center/DS-G3X_Mass_snap_128-center-cluster.txt'
		col = 'red'
		col1 = 'red'
		mark = 'o'
		linet = '-'
		filecen = '/Users/liqy/Documents/data/300data/center/G3X_Mass_snap_128-center-cluster.txt'
		
	
	if SIMmodel == 'GM':
		filex = '/Users/liqy/Documents/data/300data/simulation/GM500_snap_rr0.txt'
		filecste = '/Users/liqy/Documents/data/300data/center/DS-Music_Mass_snap_017-center-cluster.txt'
		col = 'blue'
		col1 = 'blue'
		mark = '^'
		linet = '-'
		filecen = '/Users/liqy/Documents/data/300data/center/Music_Mass_snap_017-center-cluster.txt'
	datacen = np.loadtxt(filecen)
	datacste = np.loadtxt(filecste)	

	data = np.loadtxt(file)
	if phy == 'gasz':
		data = data/0.01886
	if phy == 'density0':
		data = data * 10**10 * 0.678**2

	datacste = np.loadtxt(filecste)
	datax = np.loadtxt(filex)

	median(20,phy,data,datax,loglog,allfit,mark,col,col1,linet,lab1,lab2)

	if lines == True:
		plot(20,phy,data,datarv,loglog,lab,volumnrv = 6)


def AHF(SIMmodel,phy,loglog,allfit,lines,lab1,lab2):
	if SIMmodel == 'GX':
		col = 'red'
		col1 = 'red'
		mark = 'o'
		linet = '-'
		lab = 'AHF for Gadget3-X'
		filecste = '/Users/liqy/Documents/data/300data/center/DS-G3X_Mass_snap_128-center-cluster.txt'

	if SIMmodel == 'GM':
		col = 'blue'
		col1 = 'blue'
		mark = '^'
		linet = '-'
		lab = 'AHF for Gadget3-MUSIC'
		filecste = '/Users/liqy/Documents/data/300data/center/DS-Music_Mass_snap_017-center-cluster.txt'

	file = '/Users/liqy/Documents/data/300data/AHF/AHFdata/%s_AHFxy_%s.txt' %(SIMmodel,phy)
	data = np.loadtxt(file)

			
	if phy == 'starz' or phy == 'density4':
		filex = '/Users/liqy/Documents/data/300data/AHF/AHFdata/%s_AHFxy_rr4.txt' %SIMmodel
	if phy == 'age':
		filex = '/Users/liqy/Documents/data/300data/AHF/AHFdata/%s_AHFxy_rrage.txt' %SIMmodel
	datax = np.loadtxt(filex)


	if phy == 'starz':
		data = data / 0.01886
	
	median(20,phy,data,datax,loglog,allfit,mark,col,col1,linet,lab1,lab2)

	if lines == True:
		plot(20,phy,data,datarv,loglog,lab,volumnrv = 6)


def SAM(SAMmodel,phy,loglog,allfit,lines,lab1,lab2):

	file = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_%s.txt' %(SAMmodel,phy)
	data = np.loadtxt(file)


	if phy == 'starz':
		data = data / 0.01886

	if SAMmodel == 'Galacticus':
		mark = 'd'
		linet = '--'
		col = 'green'
		col1 = 'green'
		lab = 'Galacticus'
	if SAMmodel == 'SAG':
		mark = 's'
		linet = '--'
		col = 'orange'
		col1 = 'orange'
		lab = 'SAG'
	if SAMmodel == 'SAGE':
		mark = '*'
		linet = '--'
		col = 'black'
		col1 = 'black'
		lab = 'SAGE'

	filex = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_rr4.txt' %SAMmodel
	datax = np.loadtxt(filex)


	median(20,phy,data,datax,loglog,allfit,mark,col,col1,linet,lab1,lab2)

	if lines == True:
		plot(20,phy,data,datarv,loglog,lab,volumnrv = 6)


def SDSS7(phy,loglog,allfit,lines,lab1,lab2):

	
	filedata = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_%sPlanckBCG.txt' %phy
	if phy == 'age' or phy == 'starz':
		filex = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_rrPlanckBCG.txt'
	else:
		filex = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_rr4PlanckBCG.txt'

	data = np.loadtxt(filedata)
	datax = np.loadtxt(filex)

	if phy == 'age':
		data = 1e-9 * data
		# print(data)

	mark = 'x'
	linet = ':'
	col = 'c'
	col1 = 'c'
	lab = 'SDSS7'



	median(20,phy,data,datax,loglog,allfit,mark,col,col1,linet,lab1,lab2)
	if lines == True:
		plot(20,phy,data,datarv,loglog,lab,volumnrv = 5)

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

def main(phy,ob,allfit,lines,savepic):
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	if  phy == 'density0' or phy == 'density4':
		loglog = True
	else:
		loglog = False

#----------------------------------
#allfit: fillbetween for  median
#lines: all the data
	if phy == 'density0' or phy == 'temp' or phy == 'gasz':
		Simulation('GX',phy,loglog,allfit,lines,lab1 = 'Gadget-X',lab2 = 'fit to Gadget-X')
		Simulation('GM',phy,loglog,allfit,lines,lab1 = 'Gadget-MUSIC',lab2 = 'fit to Gadget-MUSIC')
	else:
		AHF('GX',phy,loglog,allfit,lines,lab1 = 'Gadget-X',lab2 = 'fit to Gadget-X')
		AHF('GM',phy,loglog,allfit,lines,lab1 = 'Gadget-MUSIC', lab2 = 'fit to Gadget-MUSIC')
		if phy != 'age':
			SAM('Galacticus',phy,loglog,allfit,lines,lab1 = 'Galacticus',lab2 = 'fit to Galacticus')
		SAM('SAG',phy,loglog,allfit,lines,lab1 = 'SAG',lab2 = 'fit to SAG')
		SAM('SAGE',phy,loglog,allfit,lines,lab1 = 'SAGE',lab2 = 'fit to SAGE')
		if phy == 'sfr' or phy == 'ssfr':
			loglog = False
		SDSS7(phy,loglog,allfit,lines,lab1 = 'SDSS7',lab2 = 'fit to SDSS7')

#-----------------------------------
	
	if phy == 'density4':
		
		plt.ylabel(r'$Log_{10}\ \rho_{\bigstar}\ [M_⊙kpc^{-3}]$')
		plt.xlabel('R/$R_{200}$')
		plt.ylim(0,7)
		plt.text(0.06,6,'Median Profile', fontsize = 14,fontweight = 30)
		first_legend = plt.legend(loc = 'best',ncol = 2,fontsize = 'small')
		
	if phy == 'age':
		plt.ylabel('AGE [Gyr]',fontsize = 14)
		plt.xlabel('$r_p$/$r_{200}$',fontsize = 14)
		plt.ylim(4,14)
		plt.legend(loc = 'best', fontsize = 'small', ncol = 2)
		plt.text(0.03,13.2,'Median Profile', fontsize = 14,fontweight = 30)
		plt.tick_params(labelsize = 11)
		# plt.text(0.015,4.8,'Relaxed clusters', style='italic',fontsize = 14,fontweight = 30,bbox=dict(facecolor='none',edgecolor = 'black'))

	if phy == 'starz':
		plt.ylabel(r'$Z_{\bigstar}$/$Z_⊙$',fontsize = 14)
		plt.xlabel('$r_p$/$r_{200}$', fontsize = 14)
		plt.ylim(0,4)
		plt.legend(loc = 'upper right', fontsize = 'x-small', ncol = 3)
		plt.text(0.02,3.5,'Median Profile', fontsize = 14,fontweight = 30)
		plt.tick_params(labelsize = 11)

	if phy == 'density0':
		plt.ylabel(r'$Log_{10}\ \rho_{gas}\ [M_⊙kpc^{-3}]$')
		plt.xlabel('R/$R_{500}$')
		plt.ylim(3,6)
		if ob == True:
			density0()
		plt.text(0.2,5.6,'Median Profile', fontsize = 14,fontweight = 30)
		plt.legend(loc = 'best')

	if phy == 'temp':
		plt.ylabel(r'$T/T_{500}$',fontsize = 14)
		plt.xlabel('R/$R_{500}$',fontsize = 14)
		plt.ylim(0.4,1.5)
		# plt.yticks([0.03,0.04,0.05,0.06,0.07])
		if ob == True:
			Temp()
		# plt.axvline(0.2,c = 'grey',linestyle = '--')

		# plt.text(0.03,3.5,'Unrelaxed clusters', style='italic',fontsize = 14,fontweight = 30,bbox=dict(facecolor='none',edgecolor = 'black'))
		# handles,labels = plt.get_legend_handles_labels()
		# handles = [handles[0], handles[2], handles[1]]
		# labels = [labels[0], labels[2], labels[1]]
		# plt.legend(handles,labels,loc=1)
		plt.legend(loc = 1,fontsize = 'x-small',ncol = 2)
		plt.tick_params(labelsize = 11)
	if phy == 'gasz':
		plt.tick_params(labelsize = 11)
		plt.ylabel(r'$Z_{gas}$/$Z_⊙$',fontsize = 14)
		plt.xlabel('R/$R_{500}$',fontsize = 14)
		plt.ylim(0,0.7)
		if ob == True:
			GasZ()
		plt.text(0.2,0.6,'Median Profile', fontsize = 15,fontweight = 30)
		plt.legend(loc = 'best',fontsize = 'small')

	plt.tick_params(top = 'on', right = 'on', which='both')	
	plt.xlim(0,1)

	# axes = plt.axes()
	# axes.set_yticks([0.35,0.4,0.45,0.5,0.55,0.6])
	#plt.legend(ncol = 2,ncol = 2,edgecolor = 'black',fontsize = 'xx-small',bbox_to_anchor=(0.01, 0.55,))
	#--------------------------------------------------------------------------------------------
	if savepic == True:
		if phy == 'density0' or phy == 'temp' or phy == 'gasz':
			plt.savefig('/Users/liqy/Desktop/gas_%s.pdf' %phy)
		else:
			plt.savefig('/Users/liqy/Desktop/stellar_%s.pdf' %phy)
	plt.show()

#phy: physics property
#ob:observation
#allfit: linear fitting
#lines: all datas
print('Physical property is:')
phy = input()
main(phy,ob = True,allfit = True,lines = False,savepic = True)
