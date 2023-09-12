import numpy as np 
import matplotlib.pyplot as plt

def single(SIMmodel, nclu):
	datan = np.zeros(shape = (20,4))
	filecen = '/Users/liqy/Documents/data/oPDF/cluster/G3X_Mass_snap_128-center-cluster.txt'
	datacen = np.loadtxt(filecen)

	fileden = '/Users/liqy/Documents/data/300data/simulation/%s_snaplog_density.txt' %(SIMmodel)
	filetemp = '/Users/liqy/Documents/data/300data/simulation/%s_snaplog_temp.txt' %(SIMmodel)
	filespectemp = '/Users/liqy/Documents/data/300data/simulation/%s_snaplog_spectemp.txt' %(SIMmodel)
	filemetal = '/Users/liqy/Documents/data/300data/simulation/%s_snaplog_metal.txt' %(SIMmodel)
	filex = '/Users/liqy/Documents/data/300data/simulation/GX_snaplog_gasdenx.txt'

	t500 = 8.85 * (datacen[nclu-1, 9] * 1e10 / 0.678 / 1e15) **(2/3) * (0.6125/0.6) * (0.307 + 0.693)
	den = np.loadtxt(fileden)[nclu - 1] * 0.678**2
	temp = np.loadtxt(filetemp)[nclu - 1] / 6.783
	spectemp = np.loadtxt(filespectemp)[nclu - 1]  / 6.783
	print(t500/11.491)
	metal = np.loadtxt(filemetal)[nclu - 1] / 0.01886
	datax = np.loadtxt(filex)[nclu - 1]

	# datan[:,0] = datax
	# datan[:,1] = den
	# datan[:,2] = temp
	# datan[:,3] = metal
	# np.savetxt('/Users/liqy/Desktop/GX_Cluster%03d.txt' %nclu, datan)


	plt.plot(datax,temp, label = 'EW',color = 'b',linestyle = '-')
	plt.plot(datax,spectemp, label = 'SL',color = 'b',linestyle = '--')

	plt.title('same T500 from Elena')
	plt.xscale('log')
	plt.ylim(0.5,2)
	plt.yscale('log')
	plt.xlim(1e-2,1)
	plt.xlabel('R/R500')
	plt.ylabel('T/T500')
	plt.legend()
	plt.tight_layout()
	# plt.savefig('/Users/liqy/Desktop/CL059_sameT500.pdf')
	plt.show()

#2, 12, 59, 160
# single(SIMmodel = 'GX', nclu = 2)

def alltempplot(SIMmodel):

	filecen = '/Users/liqy/Documents/data/oPDF/cluster/G3X_Mass_snap_128-center-cluster.txt'
	datacen = np.loadtxt(filecen)

	nid = np.array([2,12,59,160])
	col = np.array(['k','r','b','m'])
	tcl500 = np.array([11.207, 6.381, 6.755, 7.586])

	
	for i in range(4):
		nclu = nid[i]
		filetemp = '/Users/liqy/Documents/data/300data/simulation/%s_snaplog_temp.txt' %(SIMmodel)
		filespectemp = '/Users/liqy/Documents/data/300data/simulation/%s_snaplog_spectemp.txt' %(SIMmodel)
		filex = '/Users/liqy/Documents/data/300data/simulation/GX_snaplog_gasx.txt'

		t500 = 8.85 * (datacen[nclu-1, 9] * 1e10 / 0.678 / 1e15) **(2/3) * (0.6125/0.6) * (0.307 + 0.693)

		temp = np.loadtxt(filetemp)[nclu - 1] / tcl500[i]
		spectemp = np.loadtxt(filespectemp)[nclu - 1]  / tcl500[i]
		datax = np.loadtxt(filex)[nclu - 1]

		plt.plot(datax,temp, label = 'CL%s-EW'%nclu,color = col[i],linestyle = '-')
		plt.plot(datax,spectemp, label = 'CL%s-SL'%nclu,color = col[i],linestyle = '--')

	plt.xscale('log')
	plt.ylim(0.5,3)
	plt.yscale('log')
	plt.xlim(1e-2,1)
	plt.xlabel('R/R500')
	plt.ylabel('T/T500')
	plt.legend(ncol = 2)
	plt.tight_layout()
	plt.savefig('/Users/liqy/Desktop/temp.pdf')
	plt.show()
alltempplot(SIMmodel = 'GX')


