import numpy as np 
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

global pi
pi = 3.1415926

def median(n,phy,data,datax,mark,col,linet,lab1,lab2):
	if mark == 'o':
		filecen = '/Users/liqy/Documents/data/300data/center/G3X_Mass_snap_128-center-cluster.txt'
		datacen = np.loadtxt(filecen)
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



	if mark == '^':
		plt.errorbar(x,y,yerr = [y - ydown, yup - y],fmt = mark,markersize = 6,color = col,zorder = 3,capsize=4,label = lab1)
	if mark == 'o':
		plt.fill_between(x,ydown,yup,alpha = 0.25,color = col,zorder = 1)
		plt.scatter(x,y,marker = mark,s = 6**2,color = col,zorder = 3,label = lab1)
	# if phy == 'cold':
		# plt.errorbar(x,y,yerr = [y - ydown, yup - y],fmt = mark,markersize = 6,color = col,zorder = 3,capsize=4,label = lab1)

	#---------------------------------------------------------------------------------------------
	#---------------------------------------------------------------------------------------------

##----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

def Simulation(SIMmodel,phy,lab1,lab2,col,linet):

		
	file = '/Users/liqy/Documents/data/300data/simulation/%s_snaplog_density%s.txt' %(SIMmodel,phy)

			
	if SIMmodel == 'GX':
		filex = '/Users/liqy/Documents/data/300data/simulation/GX_snaplog_gasx%s.txt' %phy
		mark = 'o'

	if SIMmodel == 'GM':
		filex = '/Users/liqy/Documents/data/300data/simulation/GM_snaplog_gasx%s.txt' %phy
		mark = '^'
	
	data = np.loadtxt(file) * 0.678**2
	datax = np.loadtxt(filex)

	median(20,phy,data,datax,mark,col,linet,lab1,lab2)




def main():
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	plt.tick_params(top = 'on', right = 'on',which='both')

#----------------------------------
	phy = 'hot'
	Simulation('GX',phy,lab1 = 'Gadget-X (hot)',lab2 = 'fit to Gadget-X',col = 'red',linet = '-')
	Simulation('GM',phy,lab1 = 'Gadget-MUSIC (hot)',lab2 = 'fit to Gadget-MUSIC',col = 'orange',linet = '-')
	phy = 'warm'
	Simulation('GX',phy,lab1 = 'Gadget-X (warm)',lab2 = 'fit to Gadget-X',col = 'green',linet = '--')
	Simulation('GM',phy,lab1 = 'Gadget-MUSIC (warm)',lab2 = 'fit to Gadget-MUSIC',col = 'lime',linet = '--')
	phy = 'cold'
	Simulation('GX',phy,lab1 = 'Gadget-X (cold)',lab2 = 'fit to Gadget-X',col = 'blue',linet = ':')
	Simulation('GM',phy,lab1 = 'Gadget-MUSIC (cold)',lab2 = 'fit to Gadget-MUSIC',col = 'c',linet = ':')
#-----------------------------------

	plt.ylabel(r'$Log_{10}\ \rho_{gas}\ [M_⊙kpc^{-3}]$',fontsize = 14)
	plt.xlabel('R/$R_{500}$',fontsize = 14)
	plt.tick_params(labelsize = 11)
	plt.ylim(-0.5,7)
	plt.xlim(0.03,1)
	plt.xscale('log')
	# plt.yscale('log')
	# plt.yticks([0.03,0.04,0.05,0.06,0.07])

	# plt.axvline(0.2,c = 'grey',linestyle = '--')
	plt.legend(fontsize = 'small',ncol = 2,loc = 'best')

	plt.savefig('/Users/liqy/Desktop/gas_tempdensity.pdf' )
	plt.show()


main()
