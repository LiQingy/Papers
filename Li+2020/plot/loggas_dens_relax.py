import numpy as np 
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

global pi
pi = 3.1415926

def median(n,phy,ax1,ax10,data,datax,mark,col,linet,lab1,lab2,raxn,rel):
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
	# print(x,y)

	x0 = np.array(x)[7:]
	y0 = np.array(y)[7:]

	xline = np.logspace(-2,0,100)

	def func1(x,a1,b1,a2,b2,c):
		return np.log10(a1 / (1 + (x/b1)**2)**(3*c/2) + a2 / (1 + (x/b2)**2)**(3*c/2))
	if rel == 'cc':
		if col == 'red' and linet == '-':
			coeffs, matcov = curve_fit(func1, x0, y0,p0 = (1e5,0.2,1e5,0.2,0.8))
		else:
			coeffs, matcov = curve_fit(func1, x0, y0,p0 = (1e5,0.2,1e5,0.2,1))
	else:
		coeffs, matcov = curve_fit(func1, x0, y0,p0 = (1e5,0.2,1e5,0.2,0.8))
	yline = func1(xline, coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4])
	print(coeffs)
	ax1.plot(xline,yline,color = col,label = lab2,linestyle = linet)



		# x = x[8:]
		# y = y[8:]
	if raxn == 0:
		ax1.errorbar(x,y,yerr = [y - ydown, yup - y],fmt = ',',markersize = 4,color = col,zorder = 3,capsize=4,markerfacecolor = 'none',markeredgewidth = 1.5,capthick = 1 )
		if col == 'black':
			s2 = 8.5**2
			sedge = 0.8
		else:
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

		
	file = '/Users/liqy/Documents/data/300data/simulation/%s_snaplog_%s.txt' %(SIMmodel,phy)

			
	if SIMmodel == 'GX':
		filex = '/Users/liqy/Documents/data/300data/simulation/GX_snaplog_gasdenx.txt'
		filecste = '/Users/liqy/Documents/data/300data/center/DS-G3X_Mass_snap_128-center-cluster.txt'
		filecc = '/Users/liqy/Documents/data/300data/center/G3X_core.txt'
		mark = 'o'

	if SIMmodel == 'GM':
		filex = '/Users/liqy/Documents/data/300data/simulation/GM_snaplog_gasdenx.txt'
		filecste = '/Users/liqy/Documents/data/300data/center/DS-Music_Mass_snap_017-center-cluster.txt'
		filecc = '/Users/liqy/Documents/data/300data/center/Music_core.txt'
		mark = '^'

	datacste = np.loadtxt(filecste)
	datacc = np.loadtxt(filecc)
	if rel == 'cc':
		if raxn == 1:
			loc = np.where(datacc < 0.55)[0]
		else:
			loc = np.where(datacc >= 0.55)[0]
	else:
		loc = np.where(datacste[:,2] == raxn)[0]

	if phy == 'metal':
		data = data / 0.0134
	data = np.loadtxt(file)[loc]
	datax = np.loadtxt(filex)[loc]
	median(20,phy,ax1,ax10,data,datax,mark,col,linet,lab1,lab2,raxn,rel)

def main(phy,rel):
	fig = plt.figure(figsize=(6, 7.1))

	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	plt.tick_params(top = 'on', right = 'on',which='both')
	

	ax1 = plt.subplot2grid((7,7),(0,0),colspan = 7,rowspan = 5)
	ax10 = plt.subplot2grid((7,7),(5,0),colspan = 7,rowspan = 2)

#---------------------------------
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

	ax1.set_ylim(3.5,6)
	ax1.set_xlim(0.03,1)
	ax1.set_ylabel(r'$Log_{10}\ \rho_{gas}\ [M_⊙kpc^{-3}]$',fontsize = 14)
	ax1.set_xscale('log')
	ax1.axes.xaxis.set_ticklabels([])
	ax1.legend(loc = 1, fontsize = 'small')


	ax10.set_xlabel('R/$R_{500}$',fontsize = 14)
	ax10.tick_params(labelsize = 11)
	ax10.set_xlim(0.03,1)
	ax10.set_xscale('log')
	ax10.axhline(0,c='grey',zorder = 30)
	ax10.set_ylabel('scatter',fontsize = 14)
	ax10.set_ylim(-0.7,0.7)
	ax10.axes.set_yticks([-0.5,0,0.5])
	# # plt.yticks([0.03,0.04,0.05,0.06,0.07])

	# plt.axvline(0.2,c = 'grey',linestyle = '--')
	# plt.legend(fontsize = 'small',ncol = 2,loc = 'best')
	ax1.text(0.035,3.8,r'Double $\rm \beta$-model fit', fontsize = 15,fontweight = 30)
	plt.tight_layout()
	plt.subplots_adjust(wspace =0 ,hspace = 0)
	plt.savefig('/Users/liqy/Desktop/gas_%s_%s.pdf' %(phy,rel))
	plt.show()

main(phy = 'density', rel = 'cc')
