import numpy as np 
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib import ticker

def main(phy,savepic,rscale):
	fig = plt.figure(figsize=(20, 14))
	# gs = gridspec.GridSpec(15,16)
	gs0 = gridspec.GridSpec(3,2)
	gs0.update(wspace = 0.15,hspace = 0.2,left = 0.04,right = 1.01, top = 0.96,bottom = 0.05)
	gs00 = gridspec.GridSpecFromSubplotSpec(5, 8, wspace = 0,subplot_spec=gs0[0,0])
	gs01 = gridspec.GridSpecFromSubplotSpec(5, 8, wspace = 0,subplot_spec=gs0[0,1])
	gs10 = gridspec.GridSpecFromSubplotSpec(5, 8, wspace = 0,subplot_spec=gs0[1,0])
	gs11 = gridspec.GridSpecFromSubplotSpec(5, 8, wspace = 0,subplot_spec=gs0[1,1])
	gs20 = gridspec.GridSpecFromSubplotSpec(5, 8, wspace = 0,subplot_spec=gs0[2,0])
	gs21 = gridspec.GridSpecFromSubplotSpec(5, 8, wspace = 0,subplot_spec=gs0[2,1])

	ax1 = plt.subplot(gs00[:,:6])
	ax10 = plt.subplot(gs00[:,6:8])

	ax2 = plt.subplot(gs01[:,:6])
	ax20 = plt.subplot(gs01[:,6:8])

	ax3 = plt.subplot(gs10[:,:6])
	ax30 = plt.subplot(gs10[:,6:8])

	ax4 = plt.subplot(gs11[:,:6])
	ax40 = plt.subplot(gs11[:,6:8])

	ax5 = plt.subplot(gs20[:,:6])
	ax50 = plt.subplot(gs20[:,6:8])

	ax6 = plt.subplot(gs21[:,:6])
	ax60 = plt.subplot(gs21[:,6:8])

	#-------------------------------------------
	bcg = 0
	binare = 0.1 * 1 / 60 
	model = 'SAG'
	file = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_%s.txt' %(model,phy)
	# file = '/Users/liqy/Desktop/SAM_%sxy_%s.txt' %(model,phy)
	fileone = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_%s2Done.txt' %(model,phy)
	fileBCG = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_%s2DoneBCG.txt' %(model,phy)
	sfr_SAG = np.loadtxt(file)
	SAGone = np.loadtxt(fileone)
	SAGBCG = np.loadtxt(fileBCG)
	sfr_SAG = sfr_SAG / 324 / binare + bcg
	
	model = 'SAGE'
	file = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_%s.txt' %(model,phy)
	# file = '/Users/liqy/Desktop/SAM_%sxy_%s.txt' %(model,phy)
	fileone = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_%s2Done.txt' %(model,phy)
	fileBCG = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_%s2DoneBCG.txt' %(model,phy)
	sfr_SAGE = np.loadtxt(file)
	SAGEone = np.loadtxt(fileone)
	SAGEBCG = np.loadtxt(fileBCG)
	sfr_SAGE = sfr_SAGE / 324 / binare + bcg
	
	model = 'Galacticus'
	file = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_%s.txt' %(model,phy)
	# file = '/Users/liqy/Desktop/SAM_%sxy_%s.txt' %(model,phy)
	fileone = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_%s2Done.txt' %(model,phy)
	fileBCG = '/Users/liqy/Documents/data/300data/SAM/SAMdata/SAM_%sxy_%s2DoneBCG.txt' %(model,phy)
	sfr_Gala = np.loadtxt(file)
	Galaone = np.loadtxt(fileone)
	GalaBCG = np.loadtxt(fileBCG)
	sfr_Gala = sfr_Gala / 324 / binare + bcg
	print(GalaBCG.shape, np.where((GalaBCG >= -3) & (GalaBCG < 3))[0].shape[0])
	
	model = 'SDSS7'
	file = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_%sPlanckBCG.txt' %phy
	fileone = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_%s2Done.txt' %phy
	fileBCG = '/Users/liqy/Documents/data/300data/sdss7/Planck/sdss7_%sBCG.txt' %phy
	sfr_SDSS7 = np.loadtxt(file)
	SDSS7one = np.loadtxt(fileone)
	SDSS7BCG = np.loadtxt(fileBCG)
	# print(np.min(sfr_SDSS7[sfr_SDSS7 > 0]), np.max(sfr_SDSS7), np.sum(sfr_SDSS7[sfr_SDSS7 > 0]))
	sfr_SDSS7 = sfr_SDSS7 / 100 / binare + bcg
	
	model = 'Gadget-X'
	file = '/Users/liqy/Documents/data/300data/AHF/AHFdata/GX_AHFxy_%s.txt' %phy
	fileone = '/Users/liqy/Documents/data/300data/AHF/AHFdata/GX_AHFxy_%s2Done.txt' %(phy)
	fileBCG = '/Users/liqy/Documents/data/300data/AHF/AHFdata/GX_snap_%sBCG_r%s.txt' %(phy,rscale)
	sfr_GX = np.loadtxt(file)
	GXone = np.loadtxt(fileone)
	GXBCG = np.loadtxt(fileBCG)
	sfr_GX = sfr_GX / 324 / binare + bcg

	model = 'Gadget-MUSIC'
	file = '/Users/liqy/Documents/data/300data/AHF/AHFdata/GM_AHFxy_%s.txt' %phy
	fileone = '/Users/liqy/Documents/data/300data/AHF/AHFdata/GM_AHFxy_%s2Done.txt' %(phy)
	fileBCG = '/Users/liqy/Documents/data/300data/AHF/AHFdata/GM_snap_%sBCG_r%s.txt' %(phy,rscale)
	sfr_GM = np.loadtxt(file)
	GMone = np.loadtxt(fileone)
	GMBCG = np.loadtxt(fileBCG)
	sfr_GM = sfr_GM / 324 / binare + bcg
	# print('#',np.min(GMBCG),np.max(GMBCG),GMBCG.shape,GMBCG)

	# print(np.min(sfr_SAG[sfr_SAG > 0]), np.max(sfr_SAG))
	# print(np.min(sfr_SAGE[sfr_SAGE > 0]), np.max(sfr_SAGE))
	# print(np.min(sfr_Gala[sfr_Gala > 0]), np.max(sfr_Gala))
	# print(np.min(sfr_SDSS7[sfr_SDSS7 > 0]),np.max(sfr_SDSS7))
	# print(np.min(sfr_GX[sfr_GX > 0]),np.max(sfr_GX))
	# print(np.min(sfr_GM[sfr_GM > 0]),np.max(sfr_GM))


	# print(S0)
	csigma = 0.5
	corder = 0

	sfr_SAG = ndimage.gaussian_filter(sfr_SAG, sigma = csigma) + bcg
	sfr_SAGE = ndimage.gaussian_filter(sfr_SAGE, sigma = csigma) + bcg
	sfr_Gala = ndimage.gaussian_filter(sfr_Gala, sigma = csigma) + bcg
	sfr_SDSS7 = ndimage.gaussian_filter(sfr_SDSS7, sigma = csigma) + bcg
	sfr_GX = ndimage.gaussian_filter(sfr_GX, sigma = csigma) + bcg
	sfr_GM = ndimage.gaussian_filter(sfr_GM, sigma = csigma) + bcg

	# print('#')
	# print(np.min(sfr_SAG[sfr_SAG > bcg ]), np.max(sfr_SAG))
	# print(np.min(sfr_SAGE[sfr_SAGE > bcg]), np.max(sfr_SAGE))
	# print(np.min(sfr_Gala[sfr_Gala > bcg]), np.max(sfr_Gala))
	# print(np.min(sfr_SDSS7[sfr_SDSS7 > bcg]),np.max(sfr_SDSS7))
	# print(np.min(sfr_GX[sfr_GX > bcg]),np.max(sfr_GX))
	# print(np.min(sfr_GM[sfr_GM > bcg]),np.max(sfr_GM))

	S0 = sfr_SDSS7[sfr_SDSS7 > 0]
	nlv = np.array(S0.shape[0] * np.array([0.16,0.5,0.84]),dtype = 'int')
	lv = np.sort(S0)[nlv]
	# lv = [np.min(sfr_SDSS7),0.01438571]
	# print(lv)
	# print(SAGBCG.shape,SAGEBCG.shape,GalaBCG.shape,SDSS7BCG.shape,GXBCG.shape,GMBCG.shape)

	cSAG = 'orange'
	cSAGE = 'black'
	cGala = 'green'
	cSDSS7 = 'c'
	cGX = 'red'
	cGM = 'blue'
	clv = ['black','black','black']
	llv = ['-','--',':']
	lwd = 2
	bin0 = 30
	bin1 = 30
	#------------------------------------------------------------------------------------------------------------------

	if phy == 'sfr':
		high = 2.5
		xlevel = 0.04
		xy = [0,1,-3,3]
		binr = (-3,3)
	else:
		high = -9.5
		xlevel = 0.09
		xy = [0,1,-15,-9]
		binr = (-15,-9)
	# lv = np.max(sfr_SDSS7) * np.array([0.16,0.5,0.84])

	lvf = np.logspace(np.log10(1e-4),np.log10(100),100)
	# print(lvf)
	mapp = 'jet'

	ax3.contour(sfr_SDSS7,levels = lv, colors = clv, extent = xy,linestyles = llv)
	SAG = ax3.contourf(sfr_SAG, levels = lvf,cmap = mapp,extent = xy,norm = LogNorm())
	ax3.text(xlevel,high,'SAG', color = cSAG,fontsize = 19,fontweight = 70)
	ax30.hist(SAGone,bins = bin0,range = binr,histtype = 'step',orientation='horizontal',linewidth = lwd,color = cSAG,label = 'Satellite')
	ax30.hist(SAGBCG,bins = bin1,range = binr,histtype = 'stepfilled', alpha = 0.5,
		orientation='horizontal',linewidth = lwd,color = cSAG,linestyle = '--',label = 'BCG')

	ax4.contour(sfr_SDSS7,levels = lv, colors = clv,extent = xy,linestyles = llv)
	SAGE = ax4.contourf(sfr_SAGE, levels = lvf,cmap = mapp, extent = xy,norm = LogNorm())
	ax4.text(xlevel,high,'SAGE', color = cSAGE,fontsize = 19,fontweight = 50)
	ax40.hist(SAGEone,bins = bin0,range = binr,histtype = 'step',orientation='horizontal',linewidth = lwd,color = cSAGE,label = 'Satellite')
	ax40.hist(SAGEBCG,bins = bin1,range = binr,histtype = 'stepfilled', alpha = 0.5,
		orientation='horizontal',linewidth = lwd,color = cSAGE,linestyle = '--',label = 'BCG')

	ax5.contour(sfr_SDSS7,levels = lv, colors = clv,extent = xy,linestyles = llv)
	# ax3.contour(sfr_Gala,levels = lv, colors = cGala,extent = xy)
	Gala = ax5.contourf(sfr_Gala, levels = lvf,cmap = mapp, extent = xy,norm = LogNorm())
	ax5.text(xlevel,high,'Galacticus', color = cGala,fontsize = 19,fontweight = 50)
	ax50.hist(Galaone,bins = bin0,range = binr,histtype = 'step',orientation='horizontal',linewidth = lwd,color = cGala,label = 'Satellite')
	Gala0 = ax50.hist(GalaBCG,bins = bin1,range = binr, histtype = 'stepfilled', alpha = 0.5,
		orientation='horizontal',linewidth = lwd,color = cGala,linestyle = '--',label = 'BCG')
	print(np.sum(Gala0[0]))

	ax6.contour(sfr_SDSS7,levels = lv, colors = clv,extent = xy,linestyles = llv)
	SDSS7 = ax6.contourf(sfr_SDSS7, levels = lvf,cmap = mapp, extent = xy, norm = LogNorm(), order = 1)
	ax6.text(xlevel,high,'SDSS7', color = cSDSS7,fontsize = 19,fontweight = 50)
	ax60.hist(SDSS7one,bins = bin0,range = binr,histtype = 'step',orientation='horizontal',linewidth = lwd,color = cSDSS7,label = 'Satellite')
	ax60.hist(SDSS7BCG,bins = bin1,range = binr,histtype = 'stepfilled', alpha = 0.5,
		orientation='horizontal',linewidth = lwd,color = cSDSS7,linestyle = '--',label = 'BCG')

	ax1.contour(sfr_SDSS7,levels = lv, colors = clv,extent = xy,linestyles = llv)
	GX = ax1.contourf(sfr_GX, levels = lvf,cmap = mapp, extent = xy, norm = LogNorm())
	ax1.text(xlevel,high,'Gadget-X', color = cGX,fontsize = 19,fontweight = 50)
	ax10.hist(GXone,bins = bin0,range = binr,histtype = 'step',orientation='horizontal',linewidth = lwd,color = cGX,label = 'Satellite')
	ax10.hist(GXBCG,bins = bin1,range = binr,histtype = 'stepfilled', alpha = 0.5,
		orientation='horizontal',linewidth = lwd,color = cGX,linestyle = '--',label = 'BCG')
	
	ax2.text(xlevel,high,'Gadget-MUSIC', color = cGM,fontsize = 19,fontweight = 50)
	ax2.contour(sfr_SDSS7,levels = lv, colors = clv,extent = xy,linestyles = llv)
	GM = ax2.contourf(sfr_GM, levels = lvf,cmap = mapp, extent = xy,norm = LogNorm())
	ax20.hist(GMone,bins = bin0,range = binr,histtype = 'step',orientation='horizontal',linewidth = lwd,color = cGM,label = 'Satellite')	
	ax20.hist(GMBCG,bins = bin1,range = binr,histtype = 'stepfilled', alpha = 0.5,
		orientation='horizontal',linewidth = lwd,color = cGM,linestyle = '--',label = 'BCG',zorder = 100)

	axlist = [ax1,ax2,ax3,ax4,ax5,ax6,ax10,ax20,ax30,ax40,ax50,ax60]
	axlist1 = [ax1,ax2,ax3,ax4,ax5,ax6]
	axlist10 = [ax10,ax20,ax30,ax40,ax50,ax60]

	for sym in axlist:
		sym.axes.tick_params(labelsize = 12)
	for sym in axlist1:
		if phy == 'sfr':
			sym.set_ylabel(r'$log_{10}\ SFR\ [M_⊙yr^{-1}]$',fontsize = 15, fontweight = 800)
		else:
			sym.set_ylabel(r'$log_{10}\ sSFR\ [yr^{-1}]$',fontsize = 15, fontweight = 800)
		sym.set_xlabel('$r_p$/$r_{200}$',fontsize = 16, fontweight = 800)
	for sym in axlist10:
		if phy == 'sfr':
			sym.set_ylim(-3,3)
		else:
			sym.set_ylim(-15,-9)
		sym.set_xlabel('N',fontsize = 16, fontweight = 200)
		sym.set_xscale('log')
		sym.set_xlim(0,2e4)
		sym.legend(loc = 1)
		sym.axes.yaxis.set_ticklabels([])

	ax10.axes.set_xticks([10,100,1000,1e4])
	ax20.axes.set_xticks([10,100,1000,1e4])
	ax30.axes.set_xticks([10,100,1000,1e4])
	ax40.axes.set_xticks([10,100,1000,1e4])
	ax50.axes.set_xticks([10,100,1000,1e4])
	ax60.axes.set_xticks([10,100,1000,1e4])

	# print(np.min(GXBCG),np.max(GXBCG))
	fig.tight_layout()
	# fig.subplots_adjust(wspace=0.1)
	# fig.subplots_adjust(left = 0.06,bottom=0.06, right=0.999, top=0.96, wspace=0.2, hspace=0.25)
	cbar = fig.colorbar(SDSS7, ax = axlist,pad = 0.055,ticks = [1e-4,0.001,0.01,0.1,1,10,100],shrink = 0.8)
	cbar.ax.tick_params(labelsize = 12)
	fig.text(0.93,0.66,'normalised  galaxy  number',rotation=90, fontsize = 20, fontweight = 100)
	
	# plt.tight_layout()

	if savepic == True:
		plt.savefig('/Users/liqy/Desktop/%s_%s.pdf' %(phy,rscale))
	plt.show()

main(phy = 'ssfr',savepic = True,rscale = 0.015)
