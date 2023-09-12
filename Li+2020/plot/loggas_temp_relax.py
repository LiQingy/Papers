import numpy as np 
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize

from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import LogLocator
from matplotlib.ticker import NullFormatter
import matplotlib.ticker as mtick

global pi
pi = 3.1415926


def median(n,phy,ax1,ax10,data,datax,mark,col,linet,lab1,lab2,raxn,rel):
	if col == 'red':
		filecen = '/Users/liqy/Documents/data/300data/center/G3X_Mass_snap_128-center-cluster.txt'
		datacen = np.loadtxt(filecen)
	if col == 'blue':
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
			if datax[k2][k1] != 0 and math.isnan(datax[k2][k1]) == False:
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

	#---------------------------------------------------------------------------------------------
	#---------------------------------------------------------------------------------------------
	#fitting

	if phy == 'temp':
		x0 = x[9:]
		y0 = y[9:]
		xline = np.logspace(-2,0,100)
		def func2(x,t0,tmin,rcool,acool,rt,c):
			a0 = (tmin / t0 + (x / rcool)**acool)
			a1 = (1 + (x/rcool)**acool)
			b1 = (1 + (x/rt)**2)**(c/2)
			return t0 * a0 / a1 / b1
			# return t0 * (tmin + (x / rcool)**acool) / (1 + (x/rcool)**acool) / (1 + (x/rt)**2)**(c)
		if rel == 'state':
			if raxn == 1:
				if col == 'red':
					pp0 = (0.86,0.88,0.23,3.40,0.374,0.345)#best one
					coeffs, matcov = curve_fit(func2, x0, y0, p0 = pp0)
				else:
					pp0 = (0.88,0.76,0.06,4,0.348,0.54) #better one
					# pp0 = (0.68,0.788,0.003,1.24,0.98,0.9) #tests
					coeffs, matcov = curve_fit(func2, x0, y0, p0 = pp0)
					# coeffs = [0.6,1.18,0.8,0.6,0.44,0.44]
			else:
				if col == 'red':
					pp0 = (0.7,0.66,0.02,4,0.34,0.55)#best one
					coeffs, matcov = curve_fit(func2, x0, y0, p0 = pp0)
				else:
					pp0 = (0.88,0.66,0.02,4,0.34,0.55)#better one
					# pp0 = (0.7,0.66,0.027,4,0.44,0.54)#tests
					coeffs, matcov = curve_fit(func2, x0, y0, p0 = pp0)
		else:
			if raxn == 1:
				if col == 'red':
					pp0 = (0.7,0.66,0.03,2,0.42,0.55)
					coeffs, matcov = curve_fit(func2, x0, y0, p0 = pp0)
				else:
					pp0 = (0.78,0.65,0.03,2,0.38,0.55) #done、
					coeffs, matcov = curve_fit(func2, x0, y0, p0 = pp0)
					# coeffs = [0.6,1.18,0.8,0.6,0.44,0.44]
			else:
				if col == 'red':
					pp0 = (0.7,0.66,0.28,4,0.4,0.54) #done
					coeffs, matcov = curve_fit(func2, x0, y0, p0 = pp0)
				else:
					pp0 = (0.86,0.46,0.04,4.8,0.48,0.54)#done
					coeffs, matcov = curve_fit(func2, x0, y0, p0 = pp0)

		
		yline = func2(xline, coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4],coeffs[5])
		print(coeffs)
		
		ax1.plot(xline, yline, color = col, linestyle = linet,linewidth = 1.5,zorder = 2,label = lab2)

	if raxn == 0:
		ax1.errorbar(x,y,yerr = [y - ydown, yup - y],fmt = ',',markersize = 4,color = col,zorder = 3,capsize=4,markerfacecolor = 'none',markeredgewidth = 1.5,capthick = 1 )
		s2 = 7**2
		sedge = 1.2
		ax1.scatter(x,y,marker = mark,s = s2,color = col,zorder = 10,label = lab1,facecolor = 'none',linewidths = sedge)
		col = col
	else:
		ax1.fill_between(x,ydown,yup,alpha = 0.25,color = col,zorder = 1)
		ax1.scatter(x,y,marker = mark,s = 7**2,color = col,zorder = 9,label = lab1)
	#residuals
	
	ax10.plot(x,yup - y,color = col, linestyle = linet)
	ax10.plot(x,ydown - y,color = col, linestyle = linet)
##----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

def Simulation(rel,SIMmodel,phy,ax1,ax10,lab1,lab2,col,linet,raxn):

		
	file = '/Users/liqy/Documents/data/300data/simulation/%s_snaplog_spec%s.txt' %(SIMmodel,phy)

			
	if SIMmodel == 'GX':
		filex = '/Users/liqy/Documents/data/300data/simulation/GX_snaplog_gasx.txt'
		fileslT500 = '/Users/liqy/Documents/data/300data/simulation/GX_snaplog_slT500.txt'
		filecste = '/Users/liqy/Documents/data/300data/center/DS-G3X_Mass_snap_128-center-cluster.txt'
		filecc = '/Users/liqy/Documents/data/300data/center/G3X_core.txt'
		filecen = '/Users/liqy/Documents/data/300data/center/G3X_Mass_snap_128-center-cluster.txt'
		mark = 'o'

	if SIMmodel == 'GM':
		filex = '/Users/liqy/Documents/data/300data/simulation/GM_snaplog_gasx.txt'
		fileslT500 = '/Users/liqy/Documents/data/300data/simulation/GM_snaplog_slT500.txt'
		filecste = '/Users/liqy/Documents/data/300data/center/DS-Music_Mass_snap_017-center-cluster.txt'
		filecc = '/Users/liqy/Documents/data/300data/center/Music_core.txt'
		filecen = '/Users/liqy/Documents/data/300data/center/Music_Mass_snap_017-center-cluster.txt'
		mark = '^'
	datacen = np.loadtxt(filecen)	
	slT500 = np.loadtxt(fileslT500)
	datacste = np.loadtxt(filecste)
	datacc = np.loadtxt(filecc)
	if rel == 'cc':
		if raxn == 1:
			loc = np.where(datacc < 0.55)[0] #CC
		else:
			loc = np.where(datacc >= 0.55)[0] #NCC
	else:
		loc = np.where(datacste[:,2] == raxn)[0]

	t500 = 8.85 * (datacen[:, 9] * 1e10 / 0.7 / 1e15) **(2/3) * (0.6125/0.6)
	data = np.loadtxt(file)[loc] 
	nnloc = data.shape[0]
	data = np.loadtxt(file)[loc] / t500.reshape(324,1)[loc]  * 1.14 * (0.6125 / 0.588)
	datax = np.loadtxt(filex)[loc]

	median(20,phy,ax1,ax10,data,datax,mark,col,linet,lab1,lab2,raxn,rel)


def main(phy,rel):
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'

	ax1 = plt.subplot2grid((7,7),(0,0),colspan = 7,rowspan = 5)
	ax10 = plt.subplot2grid((7,7),(5,0),colspan = 7,rowspan = 2)

#----------------------------------
#allfit: fillbetween for  median
#lines: all the data
	if rel == 'cc':
		Simulation(rel,'GX',phy,ax1,ax10,lab1 = 'Gadget-X (CC)',lab2 = 'fit to Gadget-X (CC)',col = 'red',linet = '-',raxn = 1)
		Simulation(rel,'GX',phy,ax1,ax10,lab1 = 'Gadget-X (NCC)',lab2 = 'fit to Gadget-X (NCC)',col = 'red',linet = '--',raxn = 0)
		Simulation(rel,'GM',phy,ax1,ax10,lab1 = 'Gadget-MUSIC (CC)',lab2 = 'fit to Gadget-MUSIC (CC)',col = 'blue',linet = '-',raxn = 1)
		Simulation(rel,'GM',phy,ax1,ax10,lab1 = 'Gadget-MUSIC (NCC)',lab2 = 'fit to Gadget-MUSIC (NCC)',col = 'blue',linet = '--',raxn = 0)
	else:
		Simulation(rel,'GX',phy,ax1,ax10,lab1 = 'Gadget-X (relaxed)',lab2 = 'fit to Gadget-X (relaxed)',col = 'red',linet = '-',raxn = 1)
		Simulation(rel,'GX',phy,ax1,ax10,lab1 = 'Gadget-X (unrelaxed)',lab2 = 'fit to Gadget-X (unrelaxed)',col = 'red',linet = '--',raxn = 0)
		Simulation(rel,'GM',phy,ax1,ax10,lab1 = 'Gadget-MUSIC (relaxed)',lab2 = 'fit to Gadget-MUSIC (relaxed)',col = 'blue',linet = '-',raxn = 1)
		Simulation(rel,'GM',phy,ax1,ax10,lab1 = 'Gadget-MUSIC (unrelaxed)',lab2 = 'fit to Gadget-MUSIC (unrelaxed)',col = 'blue',linet = '--',raxn = 0)


#-----------------------------------

	ax1.set_xlabel('R/$R_{500}$',fontsize = 14)
	ax1.set_ylabel(r'$T/T_{500}$',fontsize = 14)
	ax10.set_ylabel(r'$scatter$',fontsize = 14)
	ax10.set_xlabel('R/$R_{500}$',fontsize = 14)

	if rel == 'state':
		ax1.set_ylim(0.4,1.6)
		y_major = LogLocator(base = 10,subs = [1])
		y_minor = LogLocator(base = 10,subs = [0.4,0.5,0.6,0.7,0.8,0.9,2])
	if rel == 'cc':
		ax1.set_ylim(0.4,1.6)
		y_major = LogLocator(base = 10,subs = [1])
		y_minor = LogLocator(base = 10,subs = [0.4,0.5,0.6,0.7,0.8,0.9,2,3])

	ax10.set_ylim(-0.5,0.5)
	ax10.set_yticks([-0.3,0,0.3])
	ax1.set_xlim(0.03,1)
	ax10.set_xlim(0.03,1)
	ax1.set_xscale('log')
	ax10.set_xscale('log')
	ax1.set_yscale('log')


	ax1.yaxis.set_major_locator(y_major)
	ax1.yaxis.set_minor_locator(y_minor)
	ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
	ax1.yaxis.set_minor_formatter(ScalarFormatter())

	ax1.tick_params(top = 'on', right = 'on',which='both')
	ax10.tick_params(right = 'on',which='both')
	# plt.text(0.3,9,'Median Profile', fontsize = 14,fontweight = 30)
	# ax1.axvline(0.2,c = 'grey',linestyle = '--')
	# ax10.axvline(0.2,c = 'grey',linestyle = '--')
	ax10.axhline(0,c='grey',zorder = 30)
	ax1.legend(loc = 'best',fontsize = 'x-small',ncol = 2)
	ax1.axes.tick_params(labelsize = 11)
	ax10.axes.tick_params(labelsize = 11)
	ax1.axes.xaxis.set_ticklabels([])

	plt.tight_layout()
	plt.subplots_adjust(wspace =0 ,hspace = 0)
	#--------------------------------------------------------------------------------------------
	plt.savefig('/Users/liqy/Desktop/gas_%s_%s.pdf' %(phy,rel))
	plt.show()

main(phy = 'temp', rel = 'state')
