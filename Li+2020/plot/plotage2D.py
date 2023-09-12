import numpy as np 
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from matplotlib import transforms
from matplotlib import ticker

def main(phy,rscale,savepic):
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	fig = plt.figure(figsize=(20, 14))

	# gs = gridspec.GridSpec(15,16)
	gs0 = gridspec.GridSpec(3,2)
	gs0.update(wspace = 0.15,hspace = 0.2,left = 0.04,right = 1.01, top = 0.96,bottom = 0.05)
	gs00 = gridspec.GridSpecFromSubplotSpec(5, 8, wspace = 0,subplot_spec=gs0[0,0])
	gs01 = gridspec.GridSpecFromSubplotSpec(5, 8, wspace = 0,subplot_spec=gs0[0,1])
	gs10 = gridspec.GridSpecFromSubplotSpec(5, 8, wspace = 0,subplot_spec=gs0[1,0])
	gs11 = gridspec.GridSpecFromSubplotSpec(5, 8, wspace = 0,subplot_spec=gs0[1,1])
	gs20 = gridspec.GridSpecFromSubplotSpec(5, 16, wspace = 0,subplot_spec=gs0[2,:])

	
	ax1 = plt.subplot(gs00[:,:6])
	ax10 = plt.subplot(gs00[:,6:8])

	ax2 = plt.subplot(gs01[:,:6])
	ax20 = plt.subplot(gs01[:,6:8])

	ax3 = plt.subplot(gs10[:,:6])
	ax30 = plt.subplot(gs10[:,6:8])

	ax4 = plt.subplot(gs11[:,:6])
	ax40 = plt.subplot(gs11[:,6:8])

	ax5 = plt.subplot(gs20[:,4:10])
	ax50 = plt.subplot(gs20[:,10:12])

	# ax6 = plt.subplot(gs21[:,:6])
	# ax60 = plt.subplot(gs21[:,6:8])
	#--------------------------------------------------------------------
	if phy == 'age2D':
		dy = 10 / 60
		ylab = 'AGE [Gyr]'
		high = 13
		xlevel = 0.09
		xy = [0,1,4,14]

	binare = dy * 1 / 60 

	model = 'SAG'
	file = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_%s.txt' %(model,phy)
	fileone = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_%sone.txt' %(model,phy)
	fileBCG = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_%soneBCG.txt' %(model,phy)
	SAG = np.loadtxt(file)
	SAGone = np.loadtxt(fileone)
	SAGBCG = np.loadtxt(fileBCG)
	SAG = SAG / 324 / binare 
	
	model = 'SAGE'
	file = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_%s.txt' %(model,phy)
	fileone = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_%sone.txt' %(model,phy)
	fileBCG = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_%soneBCG.txt' %(model,phy)
	SAGE = np.loadtxt(file)
	SAGEone = np.loadtxt(fileone)
	SAGEBCG = np.loadtxt(fileBCG)
	SAGE = SAGE / 324 / binare
	
	model = 'SDSS7'
	file = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_%sPlanckBCG.txt' %phy
	fileone = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_%sonePlanckBCG.txt' %phy
	fileBCG = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_%soneBCG.txt' %phy
	SDSS7 = np.loadtxt(file)
	SDSS7one = np.loadtxt(fileone)
	SDSS7BCG = np.loadtxt(fileBCG)
	SDSS7 = SDSS7 / 100 / binare 
	
	model = 'Gadget-X'
	file = '/Users/liqy/Documents/data/300data/AHF/AHFdata/GX_AHFxy_%s.txt' %phy
	fileone = '/Users/liqy/Documents/data/300data/AHF/AHFdata/GX_AHFxy_%sone.txt' %(phy)
	fileBCG = '/Users/liqy/Documents/data/300data/AHF/AHFdata/GX_snap_ageBCG_r%s.txt' %rscale
	GX = np.loadtxt(file)
	GXone = np.loadtxt(fileone)
	GXBCG = np.loadtxt(fileBCG)
	GX = GX / 324 / binare 

	model = 'Gadget-MUSIC'
	file = '/Users/liqy/Documents/data/300data/AHF/AHFdata/GM_AHFxy_%s.txt' %phy
	fileone = '/Users/liqy/Documents/data/300data/AHF/AHFdata/GM_AHFxy_%sone.txt' %(phy)
	fileBCG = '/Users/liqy/Documents/data/300data/AHF/AHFdata/GM_snap_ageBCG_r%s.txt' %rscale
	GM = np.loadtxt(file)
	GMone = np.loadtxt(fileone)
	GMBCG = np.loadtxt(fileBCG)
	GM = GM / 324 / binare 

	print(np.max(SDSS7one))
	# print(S0)
	csigma = 0.5
	corder = 0

	SAG = ndimage.gaussian_filter(SAG, sigma = csigma)
	SAGE = ndimage.gaussian_filter(SAGE, sigma = csigma)  
	SDSS7 = ndimage.gaussian_filter(SDSS7, sigma = csigma)
	GX = ndimage.gaussian_filter(GX, sigma = csigma) 
	GM = ndimage.gaussian_filter(GM, sigma = csigma)

	SDSS7[-1] = 0

	cSAG = 'orange'
	cSAGE = 'black'
	cSDSS7 = 'c'
	cGX = 'red'
	cGM = 'blue'

	lwd = 2
	bin0 = 30
	bin1 = 30
	binr = (4,14)
	lvf = np.logspace(np.log10(1e-4),np.log10(100),100)
	mapp = 'jet'

	SAGp = ax3.contourf(SAG, levels = lvf,cmap = mapp,extent = xy,norm = LogNorm())
	ax3.text(xlevel,high,'SAG', color = cSAG,fontsize = 19,fontweight = 500)
	ax30.hist(SAGone,bins = bin0,range = binr,histtype = 'step',orientation='horizontal',linewidth = 2,color = cSAG,label = 'Satellite')
	ax30.hist(SAGBCG,bins = bin1,range = binr,histtype = 'stepfilled', alpha = 0.5,
		orientation='horizontal',linewidth = lwd,color = cSAG,linestyle = '--',label = 'BCG')
	# ax30.axes.set_xticks([200,400,600])

	SAGEp = ax4.contourf(SAGE, levels = lvf,cmap = mapp, extent = xy,norm = LogNorm())
	ax4.text(xlevel,high,'SAGE', color = cSAGE,fontsize = 19,fontweight = 500)
	ax40.hist(SAGEone,bins = bin0,range = binr,histtype = 'step',orientation='horizontal',linewidth = 2,color = cSAGE,label = 'Satellite')
	ax40.hist(SAGEBCG,bins = bin1,range = binr,histtype = 'stepfilled', alpha = 0.5,
		orientation='horizontal',linewidth = lwd,color = cSAGE,linestyle = '--',label = 'BCG')
	# ax40.axes.set_xticks([200,400,600])

	SDSS7p = ax5.contourf(SDSS7, levels = lvf,cmap = mapp, extent = xy, norm = LogNorm(),order = 1)
	ax5.text(xlevel,high,'SDSS7', color = cSDSS7,fontsize = 19,fontweight = 500)
	ax50.hist(SDSS7one,bins = bin0,range = binr,histtype = 'step',orientation='horizontal',linewidth = 2,color = cSDSS7,label = 'Satellite')
	ax50.hist(SDSS7BCG,bins = bin1,range = binr,histtype = 'stepfilled', alpha = 0.5,
		orientation='horizontal',linewidth = lwd,color = cSDSS7,linestyle = '--',label = 'BCG')
	# ax50.axes.set_xticks([10,20,30,40,50])

	GXp = ax1.contourf(GX, levels = lvf,cmap = mapp, extent = xy,norm = LogNorm())
	ax1.text(xlevel,high,'Gadget-X', color = cGX,fontsize = 19,fontweight = 500)
	ax10.hist(GXone,bins = bin0,range = binr,histtype = 'step',orientation='horizontal',linewidth = 2,color = cGX,label = 'Satellite')
	ax10.hist(GXBCG,bins = bin1,range = binr,histtype = 'stepfilled', alpha = 0.5,
		orientation='horizontal',linewidth = lwd,color = cGX,linestyle = '--',label = 'BCG')
	# ax10.axes.set_xticks([300,600,900])

	GMp = ax2.contourf(GM, levels = lvf,cmap = mapp, extent = xy,norm = LogNorm())
	ax2.text(xlevel,high,'Gadget-MUSIC', color = cGM,fontsize = 19,fontweight = 500)
	ax20.hist(GMone,bins = bin0,range = binr,histtype = 'step',orientation='horizontal',linewidth = 2,color = cGM,label = 'Satellite')
	ax20.hist(GMBCG,bins = bin1,range = binr,histtype = 'stepfilled', alpha = 0.5,
		orientation='horizontal',linewidth = lwd,color = cGM,linestyle = '--',label = 'BCG')
	# ax20.axes.set_xticks([300,600,900])

	axlist = [ax1,ax2,ax3,ax4,ax5,ax10,ax20,ax30,ax40,ax50]
	axlist1 = [ax1,ax2,ax3,ax4,ax5]
	axlist10 = [ax10,ax20,ax30,ax40,ax50]
	for sym in axlist:
		sym.axes.tick_params(labelsize = 11)
	for sym in axlist1:
		sym.set_ylabel(ylab,fontsize = 15, fontweight = 500)
		sym.set_xlabel('$r_p$/$r_{200}$',fontsize = 16, fontweight = 800)
	for sym in axlist10:
		sym.set_ylim(4,14)
		sym.set_xlabel('N',fontsize = 16, fontweight = 200)
		sym.set_xscale('log')
		sym.set_xlim(0,4000)
		sym.axes.set_xticks([10,100,1000])
		sym.legend(loc = 1)

		sym.axes.yaxis.set_ticklabels([])

	fig.tight_layout()
	cbar = fig.colorbar(SDSS7p, ax = axlist,ticks = [1e-4,0.001,0.01,0.1,1,10,100],shrink = 0.8)
	cbar.ax.tick_params(labelsize = 12)
	fig.text(0.93,0.66,'normalised  galaxy  number',rotation=90, fontsize = 20, fontweight = 100)
	plt.savefig('/Users/liqy/Desktop/%s%s.pdf' %(phy,rscale))
	plt.show()

main(phy = 'age2D',rscale = 0.015,savepic = True)