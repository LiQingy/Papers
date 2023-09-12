import numpy as np
import matplotlib.pyplot as plt
import itertools as itert

def numb(data,rr,color,size,m,d):
	med = np.zeros(20)
	rx = np.zeros(20)
	for i in range(20):
		loc = np.where(data[:,i] != 0)[0]
		med[i] = np.median(data[loc,i], axis = 0)
		rx[i] = np.median(rr[loc,i])
		
	# print(med)
	plt.plot(rx,med,marker = m,markersize = size,label = d,color = color)
	plt.xlabel(r'$r/r_{200}$',fontsize = 14)
	plt.ylabel('median number',fontsize = 14)

def numdensity(ax1,ax10,data,datacen,rr,color,m,d):

	if d == 'SDSS7':
		nclu = data.shape[0]
		nrv = 4
	else:
		nclu = 324
		nrv = 6

	y = np.zeros(20)
	x = np.zeros(20)
	yerrdown = np.zeros(20)
	yerrup = np.zeros(20)
	for k in range(20):
		loc = np.where(np.isnan(rr[:,k]) == False)[0]
		x[k] = np.median(rr[loc,k])
		y[k] = np.median(data[loc,k])
		yerr = np.percentile(data[loc,k], [16, 84])
		yerrdown[k] = yerr[0]
		yerrup[k] = yerr[1]
	# print(y)

	if d == 'SDSS7':
		ax1.plot(x, y, linestyle = ':', color = color, label = d,linewidth = 2.5 )
		print(x)
	elif d == 'GX':
		ax1.plot(x, y, linestyle = '-', color = color, label = 'Gadget-X')
	elif d == 'GM':
		ax1.plot(x, y, linestyle = '-', color = color, label = 'Gadget-MUSIC')
	else:
		ax1.plot(x, y, linestyle = '--', color = color, label = d)
	ax1.fill_between(x, yerrdown, yerrup , color = color, alpha = 0.25)
	
	return x,y

def main(): 	
	sym = ['GX','GM','Galacticus','SAG','SAGE','SDSS7']
	col = np.array(['red','blue','g','orange','k','c'])
	ss = np.array([7,8,9,8,8,8])
	marker = np.array(['o','^','d','s','*','x'])
	ll = np.array(['-','-','--','--','--',':'])

	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	plt.tick_params(top = 'on', right = 'on', which='both')

	ax1 = plt.subplot2grid((9,8),(0,0),colspan = 8,rowspan = 7)
	ax10 = plt.subplot2grid((9,8),(7,0),colspan = 8,rowspan = 2)

	datax = np.zeros(shape = (6,20))
	datay = np.zeros(shape = (6,20))
	nd = 0
	for d in sym:
		if d == 'GX' or d =='GM':
			file = '/Users/liqy/Documents/data/300data/AHF/AHFdata/%s_AHFxy_totaln4.txt' %d
			rr = '/Users/liqy/Documents/data/300data/AHF/AHFdata/%s_AHFxy_rr4.txt' %d
			if d == 'GX':
				filecen = '/Users/liqy/Documents/data/300data/center/G3X_Mass_snap_128-center-cluster.txt'
			else:
				filecen = '/Users/liqy/Documents/data/300data/center/Music_Mass_snap_017-center-cluster.txt'
		elif d == 'SDSS7':
			file = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_totaln4PlanckBCG.txt'
			rr = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_rr4PlanckBCG.txt'
			filecen = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss_group_centerPlanckBCG.txt' 
		else:
			file = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_totaln4.txt' %d
			rr = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_rr4.txt' %d
			filecen = '/Users/liqy/Documents/data/300data/center/MDPL2_Mass_snap_128-center-cluster.txt'
		data = np.loadtxt(file)
		print(np.sum(data))
		rr = np.loadtxt(rr)
		datacen = np.loadtxt(filecen)


		# numb(data,rr,color,size,m,d)
		x,y = numdensity(ax1,ax10,data,datacen,rr,col[nd],marker[nd],d)
		datax[nd] = x
		datay[nd] = y
		nd += 1

	for i in range(6):
		ax10.plot(datax[i],datay[i]/datay[5] -1,color = col[i],linestyle = ll[i])

	ax10.tick_params(top = 'off', right = 'on',which='both',labelsize = 11)
	ax10.set_xlabel('$r_p$/$r_{200}$',fontsize = 14)
	ax10.set_ylabel(r'$\phi/\phi_{SDSS7} - 1$',fontsize = 9)
	ax10.set_xlim(0,1)
	ax10.set_ylim(-1,2.5)
	# ax10.set_ylim(-0.5,0.5)
	# ax10.axhline(0,c='cyan',linestyle = ':',zorder = 30)
	ax10.set_yticks([-0.8,0,1,2])

	ax1.set_xlim(0,1)
	ax1.set_yscale('log')
	ax1.tick_params(top = 'on', right = 'on', which='both',labelsize = 11)
	ax1.set_xlabel('$r_p$/$r_{200}$',fontsize = 14, weight = 'heavy')
	ax1.set_ylabel(r'$\phi (<r) \ [kpc^{-2}]$', fontsize = 14,fontweight = 'heavy')

	ax1.legend()
	ax1.axes.xaxis.set_ticklabels([])
	ax1.text(0.2, 4.5 * 10**-5, r'$M_{\bigstar} > 5 \times 10^{10}M_⊙$', fontsize = 16)
	ax1.text(0.22, 8 * 10**-5, 'Median profile', fontsize = 16) 

	plt.subplots_adjust(wspace =0 ,hspace = 0)
	plt.savefig('/Users/liqy/Desktop/stellar_n.pdf')
	plt.show()

main()
