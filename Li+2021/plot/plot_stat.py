import numpy as np 
import h5py
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import oPDFplus as opd

plt.figure(figsize = (7,3.5),dpi = 500)
plt.subplot(121)

plt.tick_params(top = 'on', right = 'on', which='both',direction = 'in')
# plt.rcParams['xtick.direction'] = 'in'
# plt.rcParams['ytick.direction'] = 'in'

fileChi2 = '/Users/liqy/Documents/oPDF1/data/GXsub_DM_Chi2_rcin200_100000_TMP.txt'
dataChi2 = np.loadtxt(fileChi2)
dx0 = 0.02
xl=np.log10(1)-dx0
xu=np.log10(1)+dx0
x=np.logspace(xl, xu, 20) 
y=np.logspace(xl, xu, 20)

# plt.figure(figsize = (3.6,3.4))


plt.contour(np.log10(x),np.log10(y),dataChi2,levels=[2.3],colors = 'grey',zorder = 100)

fval0= np.loadtxt('/Users/liqy/Documents/data/oPDF1/HE/HE_Chi2_rcin200_TMP_n100000_bin20_boot200.txt')
xx = np.linspace(-0.01,0.01,100)  
yy = np.linspace(-0.01,0.01,100)
plt.contour(xx, yy, fval0, levels = [2.3], colors = ['purple'])

fval0= np.loadtxt('/Users/liqy/Documents/data/oPDF1/JeansE/JeansE_Chi2_rcin200_TMP_n100000_bin20_boot200.txt')
xx = np.linspace(-0.02,0.02,100)  
yy = np.linspace(-0.02,0.02,100)
cs1 = plt.contour(xx, yy, fval0, levels = [2.3], colors = ['brown'])

plt.xlim(-0.02,0.02)
plt.ylim(-0.02,0.02)
plt.ylabel(r'$\mathrm{log}(c/c_{\rm best})$',fontsize = 12)
plt.xlabel(r'$\mathrm{log}(M/M_{\rm best})$',fontsize = 12)
plt.axvline(0,color = 'k',linestyle = '--',lw = 1)
plt.axhline(0,color = 'k',linestyle = '--',lw = 1)
plt.tick_params(labelsize = 10)
plt.xticks([-0.02,-0.01,0,0.01,0.02])
plt.yticks([-0.02,-0.01,0,0.01,0.02])
plt.text(-0.019,-0.015,'Statistical error', fontdict = {'fontsize':10,'fontweight': 400})
plt.text(-0.017,0.015,'Efficiency', fontsize = 10, weight = 'bold')

patch1 = Patch(edgecolor='grey', label='oPDF',facecolor='white')
patch2 = Patch(edgecolor='brown', label='SJE',facecolor='white')
patch3 = Patch(edgecolor='purple', label='GHE',facecolor='white')
# patch4 = Patch(edgecolor='magenta', label='VT',facecolor='white')

plt.legend(handles=[patch1,patch2,patch3], loc = 1)

bwith = 1.2 #边框宽度设置为2
ax = plt.gca()#获取边框
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)

plt.subplot(122)

labels = ['oPDF', 'SJE', 'GHE']
Nint = [200,60,20]
col = ['grey', 'brown', 'purple']
plt.bar(labels, Nint, color = col, width = 0.3, alpha = 0.7)
plt.ylabel(r'$N_{\rm int}$', fontsize = 12)
plt.tick_params(top = 'on', right = 'on', which='both',direction = 'in',labelsize = 10)
plt.tick_params(axis = 'x', labelsize = 12)

plt.text(0.9,185,'Robustness', fontsize = 10, weight = 'bold')

bwith = 1.2 #边框宽度设置为2
ax = plt.gca()#获取边框
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)


plt.tight_layout()
plt.savefig('/Users/liqy/Documents/oPDF1/fig/stat.png')
plt.show()

#test for Galactic halo
# fval = np.loadtxt(open('/Users/liqy/Desktop/Gala_DM_staerr.csv', 'r'),delimiter = ',')
# Gala = plt.scatter(fval[:,0]-0.2,fval[:,1]-0.2,c = 'k',label = 'Galactic-oPDF')

#virial statistic error
# fval0= np.loadtxt('/Users/liqy/Documents/data/oPDF1/VirM/VirM_Chi2_rcin200_TMP_n1000_bin20_boot200.txt')
# xx = np.linspace(-0.08,0.08,100)  
# yy = np.linspace(-0.08,0.08,100)
# cs1 = plt.contour(xx, yy, fval0, levels = [2.3], colors = ['magenta'])
# pc1=cs1.collections[0].get_paths()[0]
# con1=pc1.vertices
# xcon1=con1[:,0]
# ycon1=con1[:,1]
# np.savetxt('/Users/liqy/Documents/data/oPDF1/VirM/VirM_sta_linedata',np.vstack((xcon1,ycon1)))

# VirMline = np.loadtxt('/Users/liqy/Documents/data/oPDF1/VirM/VirM_sta_linedata')
# print(VirMline.shape)
# plt.plot(VirMline[0]/10,VirMline[1]/10,c = 'magenta')
