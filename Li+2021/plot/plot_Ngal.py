import numpy as np 
import matplotlib.pyplot as plt

plt.figure(figsize = (3.5,3.4))

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.tick_params(top = 'on', right = 'on',which='both')

filestenub = '/Users/liqy/Documents/data/oPDFnew/cluster/GX_ste_Ngal.txt'
stenub = np.loadtxt(filestenub)
fileMtrue = '/Users/liqy/Documents/data/oPDFnew/cluster/G3X_Mass_snap_128-center-cluster.txt'
Mtrue = np.loadtxt(fileMtrue)

print(np.min(stenub))
print(np.mean(stenub),np.median(stenub))

plt.scatter(np.log10(Mtrue[:,2]), stenub, s = 3)

plt.ylabel(r'$N_{\rm gal}$', fontsize = 12)
plt.xlabel(r'$\mathrm{log_{10}}\ M_{\mathrm{true}}[\mathrm{M_{\odot}}/h]$', fontsize = 12)
plt.tick_params(labelsize = 10)

plt.text(14.65,520, r'$200\ \mathrm{kpc}/h\ \rm cut$', fontdict = {'fontsize':9,'fontweight': 400})
plt.text(14.67,580, r'$M_{\bigstar}>10^9\ \rm M_{\odot}$', fontdict = {'fontsize':9,'fontweight': 400})

bwith = 1.2 #边框宽度设置为2
ax = plt.gca()#获取边框
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)

plt.tight_layout()

plt.savefig('/Users/liqy/Desktop/oPDF1/paper/Ngal.pdf')
plt.show()






