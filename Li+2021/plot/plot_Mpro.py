import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import LogLocator
from matplotlib.ticker import NullFormatter
import matplotlib.ticker as mtick


def pro(nbin, model, pptype, fitm, nubtype, rcin, linet, col, mak):
	if model == 'GX' or model == 'GXsub':
		filecen = '/Users/liqy/Documents/data/oPDFnew/cluster/G3X_Mass_snap_128-center-cluster.txt'
		fileMtrue = '/Users/liqy/Documents/data/oPDF1/Masspro/GX_MassAllprofile_log.txt'
		
	else:
		filecen = '/Users/liqy/Documents/data/oPDFnew/cluster/MDPL2_Mass_snap_128-center-cluster.txt'
		fileMtrue = '/Users/liqy/Documents/data/oPDF1/Masspro/MDPL2_MassAllprofile_log.txt'
	
	fileMpro = '/Users/liqy/Desktop/Trueini/%ssub_%s_rin%s_rout1_%s_%s_M.txt' %(model, pptype, rcin,nubtype,fitm)
	datacen = np.loadtxt(filecen)
	Mtrue = np.loadtxt(fileMtrue)
	dataMpro = np.loadtxt(fileMpro) * 1e10

	Mbias = dataMpro / Mtrue
	Mpro = np.median(Mbias, axis = 0)
	sigma = np.std(Mbias,axis = 0)
	perc = np.percentile(Mbias, [14,86],axis = 0)
	# print(Mpro_TMP,perc_TMP)
	
	rv = datacen[:,6]
	xx = 10**np.linspace(np.log10(rv)/nbin,np.log10(rv),nbin)/rv
	xx = np.median(xx, axis = 1)

	yup = Mpro + sigma
	ydown = Mpro - sigma

	TMPs = plt.scatter(xx, Mpro, label = '%s(%s)' %(model,fitm), color = col, marker=mak, zorder = 10)
	TMPl, = plt.plot(xx, perc[0],  color = col, ls = linet, label = 'percentile error: %s(%s)'%(model,fitm),zorder = 5)
	plt.plot(xx, perc[1],  color = col, ls = linet, zorder = 5)
		# plt.errorbar(xx, Mpro_TMP, yerr = sigma_TMP,label = '%s(TMP)'%model, elinewidth = 1, capsize = 1.5
	# 	,color = col, linestyle = linet, fmt = 'o',zorder = 10,ms = 4)
	return TMPs,TMPl


def main(nbin):
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['ytick.direction'] = 'in'
	plt.tick_params(top = 'on', right = 'on', which='both')

	#--------------
	GXNFWs,GXNFWl = pro(nbin, model = 'GX', pptype = 'DM', fitm = 'NFW', nubtype = '100000', rcin = '0', linet = '--', col = 'r', mak = 'o')
	MNFWs,MNFWl = pro(nbin, model = 'MDPL2', pptype = 'DM', fitm = 'NFW', nubtype = '100000', rcin = '0', linet = '--', col = 'blue', mak = 'o')

	GXTMPs,GXTMPl = pro(nbin, model = 'GX', pptype = 'DM', fitm = 'TMP', nubtype = '100000', rcin = '0', linet = '-', col = 'r', mak = '^')
	MTMPs, MTMPl = pro(nbin, model = 'MDPL2', pptype = 'DM', fitm = 'TMP', nubtype = '100000', rcin = '0', linet = '-', col = 'blue', mak ='^')
	#---------------
	leg0 = plt.legend(handles = [GXNFWl,MNFWl,GXTMPl,MTMPl],loc = 8)
	plt.gca().add_artist(leg0)
	plt.legend(handles = [GXNFWs,MNFWs,GXTMPs,MTMPs],loc = 1)
	#---------------
	plt.xscale('log')
	plt.xlabel(r'$r/r_{200}$',fontsize = 14)
	plt.ylim(4e-1,2)
	plt.ylabel(r'$M(<r)/M_{\mathrm{true}}(<r)$', fontsize = 14)
	plt.yscale('log')

	plt.tick_params(labelsize = 11)
	plt.xlim(1e-2,1)
	plt.axhline(1,color = 'k',linestyle='-', lw = 1.2)

	bwith = 1.2 #边框宽度设置为2
	ax = plt.gca()
	ax.spines['bottom'].set_linewidth(bwith)
	ax.spines['left'].set_linewidth(bwith)
	ax.spines['top'].set_linewidth(bwith)
	ax.spines['right'].set_linewidth(bwith)
	
	y_major = LogLocator(base = 10,subs = [1])
	y_minor = LogLocator(base = 10,subs = [0.4,0.5,0.6,0.7,0.8,0.9,2])
	ax = plt.gca()
	ax.yaxis.set_major_locator(y_major)
	ax.yaxis.set_minor_locator(y_minor)

	ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
	ax.yaxis.set_minor_formatter(ScalarFormatter())

	plt.savefig('/Users/liqy/Desktop/figs_Trueini/DMO_SPH.pdf')
	plt.show()
main(nbin = 20)






