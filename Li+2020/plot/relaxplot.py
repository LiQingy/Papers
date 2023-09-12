import numpy as np 
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize

global pi
pi = 3.1415926


def median(n,phy,data,datax,r200,M200,loglog,mark,col,col1,linet,lab1,lab2,raxn,ax0,ax00):
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

	#---------------------------------------------------------------------------------------------
	#---------------------------------------------------------------------------------------------
	#fitting

	if phy == 'age' or phy == 'temp':
		x0 = x[:]
		y0 = y[:]

		xline = np.linspace(0,1,100)
		def func2(x,t0,tmin,rcool,acool,rt,c):
			a0 = (tmin / t0 + (x / rcool)**acool)
			a1 = (1 + (x/rcool)**acool)
			b1 = (1 + (x/rt)**2)**(c/2)
			return t0 * a0 / a1 / b1
			# return t0 * (tmin + (x / rcool)**acool) / (1 + (x/rcool)**acool) / (1 + (x/rt)**2)**(c)
		if raxn == 1:
			if col == 'red':
				pp0 = (1.22,1.14,4.4e-3,1.04,0.33,0.54)
				coeffs, matcov = curve_fit(func2, x0, y0, p0 = pp0, bounds = ([0,0.01,0,0,0,0],[2,5,2,30,1,1]))
			else:
				pp0 = (1.21,0.86,4.6e-3,1.03,0.34,0.54) #done、
				coeffs, matcov = curve_fit(func2, x0, y0, p0 = pp0)
		else:
			if col == 'red':
				pp0 = (1.24,1.14,1.6e-3,1.08,0.35,0.56)
				coeffs, matcov = curve_fit(func2, x0, y0, p0 = pp0, bounds = ([0,0.01,0,0,0,0],[2,5,2,25,1,1]))
			else:
				pp0 = (1.22,1.2,1.7e-3,1.03,0.34,0.54)#done
				coeffs, matcov = curve_fit(func2, x0, y0, p0 = pp0)
		
		yline = func2(xline, coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4],coeffs[5])
		print(coeffs)
		
		ax0.plot(xline, yline, color = col, linestyle = linet,linewidth = 1.5,zorder = 2,label = lab2)

	else:
		if phy == 'gasz' or phy == 'density0':

			xline = np.linspace(0,1,100)
			def func2(x,a,b):
				return np.log10(a / (x/b * (1 + x/b)**2))
			coeffs, matcov = curve_fit(func2, x0, y0)
			yline = func2(xline, coeffs[0],coeffs[1])
			print(coeffs)
			ax0.plot(xline, yline, color = col, linestyle = linet,linewidth = 1.5,zorder = 2,label = lab2)
		
		if phy == 'density4' or phy == 'mass4':
			
			x0 = x[1:]
			y0 = y[1:]
			
			xline = np.linspace(0,1,100)		
			def func(x,b):
				return np.log10(M200 / (4/3*np.pi*(r200*x)**3)* (np.log(1 + x*b) - (x*b)/(1 + x*b))/ (np.log(1 + b) - (b)/(1 + b)) )

			coeffs, matcov = curve_fit(func, x0, y0, p0 = 5,bounds = (0.01,30))
			yline = func(xline, coeffs[0])
			print(coeffs,r200/1e3,M200/1e12)

			# nv = np.zeros(20)
			# for i in range(20):
			# 	nv[i] = np.pi * (r200 * x[i])**3 * 4 /3
			# y = np.log10(10**y / nv)
			# yup = np.log10(10**yup / nv)
			# ydown = np.log10(10**ydown / nv)
			# nvv = np.zeros(100)
			# for j in range(100):
			# 	nvv[j] = np.pi * (r200 * xline[j])**3 * 4 /3
			# yline = np.log10(10**yline / nvv)

			ax0.plot(xline, yline, color = col, linestyle = linet,linewidth = 1.5,zorder = 2,label = lab2)
	if raxn == 0:
		ax0.errorbar(x,y,yerr = [y - ydown, yup - y],fmt = ',',markersize = 4,color = col1,zorder = 3,capsize=4,markerfacecolor = 'none',markeredgewidth = 1.5,capthick = 1 )
		if col1 == 'black':
			s2 = 8.5**2
			sedge = 0.8
		else:
			s2 = 7**2
			sedge = 1.2
		ax0.scatter(x,y,marker = mark,s = s2,color = col1,zorder = 10,label = lab1,facecolor = 'none',linewidths = sedge)
		col = col1
	else:
		ax0.fill_between(x,ydown,yup,alpha = 0.25,color = col,zorder = 1)
		ax0.scatter(x,y,marker = mark,s = 7**2,color = col,zorder = 9,label = lab1)
	#residuals
	
	ax00.plot(x,yup - y,color = col, linestyle = linet)
	ax00.plot(x,ydown - y,color = col, linestyle = linet)
##----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

def Simulation(SIMmodel,phy,loglog,ax1,ax2,ax10,ax20,raxn):
	if raxn == 1:
		linet = '--'
		state = 'relaxed'
	else:
		linet = '-'
		state = 'unrelaxed'

	file = '/Users/liqy/Documents/data/300data/simulation/%s500_snap_%s.txt' %(SIMmodel,phy)

			
	if SIMmodel == 'GX':
		filex = '/Users/liqy/Documents/data/300data/simulation/GX500_snap_rr0.txt'
		filecste = '/Users/liqy/Documents/data/300data/center/DS-G3X_Mass_snap_128-center-cluster.txt'
		col = 'red'
		col1 = 'red'
		mark = 'o'
		lab1 = 'Gadget-X (%s)' %state
		lab2 = 'fit to Gadget-X (%s)' %state
		
		filecen = '/Users/liqy/Documents/data/300data/center/G3X_Mass_snap_128-center-cluster.txt'
		
	
	if SIMmodel == 'GM':
		filex = '/Users/liqy/Documents/data/300data/simulation/GM500_snap_rr0.txt'
		filecste = '/Users/liqy/Documents/data/300data/center/DS-Music_Mass_snap_017-center-cluster.txt'
		col = 'blue'
		col1 = 'blue'
		mark = '^'
		lab1 = 'Gadget-MUSIC (%s)' %state
		lab2 = 'fit to Gadget-MUSIC (%s)' %state

		filecen = '/Users/liqy/Documents/data/300data/center/Music_Mass_snap_017-center-cluster.txt'
	datacen = np.loadtxt(filecen)
	datacste = np.loadtxt(filecste)
	loc = np.where(datacste[:,2] == raxn)

	data = np.loadtxt(file)[loc]
	if phy == 'gasz':
		data = data/0.0194
	if phy == 'density0':
		data = data * 10**10 * 0.678**2

	datax = np.loadtxt(filex)[loc]

	r200 = 1
	M200 = 1
	median(20,phy,data,datax,r200,M200,loglog,mark,col,col1,linet,lab1,lab2,raxn,ax1,ax10)



def AHF(SIMmodel,phy,loglog,ax1,ax2,ax10,ax20,raxn):
	if raxn == 1:
		linet = '-'
		state = 'relaxed'
	else:
		linet = '--'
		state = 'unrelaxed'

	if SIMmodel == 'GX':
		col = 'red'
		col1 = 'red'
		mark = 'o'
		lab1 = 'Gadget-X (%s)' %state
		lab2 = 'fit to Gadget-X (%s)' %state

		filecen = '/Users/liqy/Documents/data/300data/center/G3X_Mass_snap_128-center-cluster.txt'
		filecste = '/Users/liqy/Documents/data/300data/center/DS-G3X_Mass_snap_128-center-cluster.txt'

	if SIMmodel == 'GM':
		col = 'blue'
		col1 = 'blue'
		mark = '^'
		lab1 = 'Gadget-MUSIC (%s)' %state 
		lab2 = 'fit to Gadget-MUSIC (%s)' %state

		filecen = '/Users/liqy/Documents/data/300data/center/Music_Mass_snap_017-center-cluster.txt'
		filecste = '/Users/liqy/Documents/data/300data/center/DS-Music_Mass_snap_017-center-cluster.txt'
	filegal200 = '/Users/liqy/Documents/data/300data/AHF/AHFdata/%s_AHF_mass200gal.txt' %SIMmodel
	gal200 = np.loadtxt(filegal200)

	datacste = np.loadtxt(filecste)
	loc = np.where(datacste[:,2] == raxn)
	datacen = np.loadtxt(filecen)

	file = '/Users/liqy/Documents/data/300data/AHF/AHFdata/%s_AHF_%s.txt' %(SIMmodel,phy)
	data = np.loadtxt(file)
	data = data[loc]
			
	if phy == 'starz' or phy == 'density4' or phy == 'mass4':
		filex = '/Users/liqy/Documents/data/300data/AHF/AHFdata/%s_AHF_rr4.txt' %SIMmodel
	if phy == 'age':
		filex = '/Users/liqy/Documents/data/300data/AHF/AHFdata/%s_AHF_rrage.txt' %SIMmodel
	datax = np.loadtxt(filex)
	datax = datax[loc]
	r200 = np.median(datacen[loc[0],6]) / 0.678
	M200 = np.median(gal200[loc[0]])
		
	median(20,phy,data,datax,r200,M200,loglog,mark,col,col1,linet,lab1,lab2,raxn,ax1,ax10)

def SAM(SAMmodel,phy,loglog,ax1,ax2,ax10,ax20,raxn):
	filecen = '/Users/liqy/Documents/data/300data/center/MDPL2_Mass_snap_128-center-cluster.txt'
	filecste = '/Users/liqy/Documents/data/300data/center/DS-MDPL2_Mass_snap_128-center-cluster.txt'
	datacste = np.loadtxt(filecste)
	datacen = np.loadtxt(filecen)
	loc = np.where(datacste[:,2] == raxn)

	file = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%s_%s.txt' %(SAMmodel,phy)
	data = np.loadtxt(file)
	data = data[loc]

	if raxn == 1:
		linet = '-'
		state = 'relaxed'
	else:
		linet = '--'
		state = 'unrelaxed'

	if SAMmodel == 'Galacticus':
		mark = 'd'
		col = 'green'
		col1 = 'green'
		lab1 = 'Galacticus (%s)' %state
		lab2 = 'fit to Galacticus (%s)' %state

	if SAMmodel == 'SAG':
		mark = 's'
		col = 'orange'
		col1 = 'orange'
		lab1 = 'SAG (%s)' %state
		lab2 = 'fit to SAG (%s)' %state

	if SAMmodel == 'SAGE':
		mark = '*'
		col = 'black'
		col1 = 'black'
		lab1 = 'SAGE (%s)' %state
		lab2 = 'fit to SAGE (%s)' %state

	filex = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%s_rr4.txt' %SAMmodel
	datax = np.loadtxt(filex)
	datax = datax[loc]
	filegal200 = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%s_mass200gal.txt' %SAMmodel
	gal200 = np.loadtxt(filegal200)
	r200 = np.median(datacen[loc[0],6]) / 0.678
	M200 = np.median(gal200[loc[0]])

	median(20,phy,data,datax,r200,M200,loglog,mark,col,col1,linet,lab1,lab2,raxn,ax2,ax20)


def main(phy,ob,savepic):
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	if phy == 'density4' or phy == 'age' or phy == 'mass4':
		fig = plt.figure(dpi =100,figsize=(12, 6))

		ax1 = plt.subplot2grid((7,14),(0,0),colspan = 7,rowspan = 5)
		ax2 = plt.subplot2grid((7,14),(0,7),colspan = 7,rowspan = 5)
		ax10 = plt.subplot2grid((7,14),(5,0),colspan = 7,rowspan = 2)
		ax20 = plt.subplot2grid((7,14),(5,7),colspan = 7,rowspan = 2)
	if phy == 'density4' or phy == 'mass4':
		loglog = True
	else:
		loglog = False
	if phy == 'temp':
		ax1 = plt.subplot2grid((7,7),(0,0),colspan = 7,rowspan = 5)
		ax10 = plt.subplot2grid((7,7),(5,0),colspan = 7,rowspan = 2)
		ax2 = 0
		ax20 =0
#----------------------------------
#allfit: fillbetween for  median
#lines: all the data
	if phy == 'density0' or phy == 'temp' or phy == 'gasz':
		Simulation('GX',phy,loglog,ax1,ax2,ax10,ax20,raxn = 1)
		Simulation('GM',phy,loglog,ax1,ax2,ax10,ax20,raxn = 1)
	else:
		AHF('GX',phy,loglog,ax1,ax2,ax10,ax20,raxn = 1)
		AHF('GM',phy,loglog,ax1,ax2,ax10,ax20,raxn = 1)
		if phy != 'age':
			SAM('Galacticus',phy,loglog,ax1,ax2,ax10,ax20,raxn  =1)
		SAM('SAG',phy,loglog,ax1,ax2,ax10,ax20,raxn =1)
		SAM('SAGE',phy,loglog,ax1,ax2,ax10,ax20,raxn = 1)

	if phy == 'density0' or phy == 'temp' or phy == 'gasz':
		Simulation('GX',phy,loglog,ax1,ax2,ax10,ax20,raxn = 0)
		Simulation('GM',phy,loglog,ax1,ax2,ax10,ax20,raxn = 0)
	else:
		AHF('GX',phy,loglog,ax1,ax2,ax10,ax20,raxn = 0)
		AHF('GM',phy,loglog,ax1,ax2,ax10,ax20,raxn = 0)
		if phy != 'age':
			SAM('Galacticus',phy,loglog,ax1,ax2,ax10,ax20,raxn  =0)
		SAM('SAG',phy,loglog,ax1,ax2,ax10,ax20,raxn =0)
		SAM('SAGE',phy,loglog,ax1,ax2,ax10,ax20,raxn = 0)

#-----------------------------------
	
	if phy == 'density4' or  phy == 'mass4':
		
		ax1.set_ylabel(r'$Log_{10}\ \rho_{\bigstar}(<R)\ [M_⊙kpc^{-3}]$',fontsize = 14)
		ax10.set_ylabel(r'$scatter$',fontsize = 14)
		ax10.set_xlabel('R/$R_{200}$',fontsize = 14)
		ax20.set_xlabel('R/$R_{200}$',fontsize = 14)
		ax1.set_ylim(1,6)
		ax2.set_ylim(1,6)
		ax10.set_ylim(-1,1)
		ax20.set_ylim(-1,1)
		ax10.set_yticks([-1,-0.5,0,0.5])
		# ax1.text(0.06,6,'Median Profile', fontsize = 14,fontweight = 30)
		# ax2.text(0.06,6,'Median Profile', fontsize = 14,fontweight = 30)
	
	if phy == 'age':
		ax1.set_ylabel('AGE [Gyr]',fontsize = 14)
		ax10.set_ylabel(r'$scatter$',fontsize = 14)
		
		ax10.set_xlabel('R/$R_{200}$',fontsize = 14)
		ax20.set_xlabel('R/$R_{200}$',fontsize = 14)
		ax1.set_ylim(6,13)
		ax2.set_ylim(6,13)
		ax10.set_ylim(-3,3)
		ax20.set_ylim(-3,3)


		# plt.text(0.2,13,'Median Profile', fontsize = 14,fontweight = 30)
		# plt.text(0.015,4.8,'Relaxed clusters', style='italic',fontsize = 14,fontweight = 30,bbox=dict(facecolor='none',edgecolor = 'black'))


	if phy == 'temp':
		# ax1.axes.xaxis.set_ticklabels([])
		ax1.set_xlabel('R/$R_{500}$',fontsize = 14)
		ax1.set_ylabel(r'$T/T_{500}$',fontsize = 14)
		ax10.set_ylabel(r'$scatter$',fontsize = 14)
		ax10.set_xlabel('R/$R_{500}$',fontsize = 14)
		ax1.axes.xaxis.set_ticklabels([])
		ax1.set_ylim(0.4,1.5)
		# ax10.set_ylim(-0.05,0.05)
		ax1.set_xlim(0,1)
		ax10.set_xlim(0,1)
		ax1.tick_params(top = 'on', right = 'on',which='both')
		ax10.tick_params(right = 'on',which='both')
		# ax10.set_yticks([-0.03,0,0.03])
		# plt.text(0.3,9,'Median Profile', fontsize = 14,fontweight = 30)
		ax1.axvline(0.2,c = 'grey',linestyle = '--')
		ax10.axvline(0.2,c = 'grey',linestyle = '--')
		ax10.axhline(0,c='grey',zorder = 30)
		ax1.legend(loc = 1,fontsize = 'small',ncol = 2)
		ax1.axes.tick_params(labelsize = 11)
		ax10.axes.tick_params(labelsize = 11)

	if phy == 'density4' or phy == 'age' or phy == 'mass4':
		ax1.axes.xaxis.set_ticklabels([])
		ax2.axes.xaxis.set_ticklabels([])
		ax2.axes.yaxis.set_ticklabels([])
		ax20.axes.yaxis.set_ticklabels([])
		ax1.tick_params(top = 'on', right = 'on',which='both')
		ax2.tick_params(top = 'on',right = 'on', left = 'off',which='both')

		ax10.tick_params(top = 'off', right = 'on',which='both')
		ax20.tick_params(top = 'off',right = 'on', left = 'off',which='both')
		# ax10.set_xlim(0,1)
		# ax20.set_xlim(0,1)
		ax1.set_xlim(0,1)
		ax2.set_xlim(0,1)
		ax10.axhline(0,c='grey',zorder = 30)
		ax10.set_xlim(0,1)
		ax20.set_xlim(0,1)
		ax10.axes.set_xticks([0,0.2,0.4,0.6,0.8])

		ax1.get_legend_handles_labels()
		ax1.legend(loc = 'best',ncol = 2,fontsize = 'small')
		ax2.get_legend_handles_labels()
		ax2.legend(loc = 'best',ncol = 2,fontsize = 'small')
		ax20.axhline(0,c='grey',zorder = 30)
	if phy != 'temp':
		axlist = [ax1,ax10,ax2,ax20]
		for sym in axlist:
			sym.axes.tick_params(labelsize = 11)
		# ax20.set_xlim(0,1)
	plt.tight_layout()
	plt.subplots_adjust(wspace =0 ,hspace = 0)
	#--------------------------------------------------------------------------------------------
	if savepic == True:
		if phy == 'temp':
			plt.savefig('/Users/liqy/Desktop/gas_%s_state.pdf' %phy)
		else:
			plt.savefig('/Users/liqy/Desktop/stellar_%s_state.pdf' %phy)
	plt.show()

#phy: physics property
#ob:observation
#allfit: linear fitting
#lines: all datas
print('Physical property is:')
phy = input()
main(phy,ob = True,savepic = True)
