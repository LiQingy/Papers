from oPDF.oPDF import *
from oPDF.myutils import *
import h5py,os,sys
from scipy.stats import chi2
import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec

#pl.ion()

cenpos1=np.array([[6.60408545,12.36700916,59.14881516],
[20.28102112,46.4294548,12.07352352],
[36.32865524,20.41593742,8.41809368],
[44.44973373,14.18295765,51.0191803],
[30.0664444,61.80264282,65.65539551],
[25.31764984,42.23421097,30.82751656],
[3.25924897,51.21772766,60.2011261],
[29.28999329,30.9513073,45.64191437],
[40.15666962,27.09303284,52.6264534],
[43.61115646,17.05785561,34.00644684],
[44.04451752,33.84285736,24.80952072],
[33.47823334,27.86206627,15.55616951]])

cenpos2=np.array([[6.32896137,11.9668541,59.49020386],
[19.97343063,46.44976425,11.58646107],
[35.97160721,19.88262939,8.46210289],
[44.8338356,13.74028492,50.96887589],
[30.17848778,62.28554535,65.94735718],
[25.45763397,41.6999855,30.46043205],
[3.24130416,51.75031662,59.98556137],
[29.44311333,30.48748016,45.67856598],
[40.47468567,26.67650032,52.41774368],
[43.55581284,17.55819321,33.88335037],
[43.9861908,33.63922882,25.31379318],
[33.31719208,28.26669693,15.78494644]]
)


sep=np.zeros(24)
for i in range(12):
  sep[i*2]=np.sqrt((cenpos1[i,0]-cenpos2[i,0])**2+(cenpos1[i,1]-cenpos2[i,1])**2+(cenpos1[i,2]-cenpos2[i,2])**2)*1000./0.704
  sep[i*2+1]=np.sqrt((cenpos1[i,0]-cenpos2[i,0])**2+(cenpos1[i,1]-cenpos2[i,1])**2+(cenpos1[i,2]-cenpos2[i,2])**2)*1000./0.704

print(sep)

nbin=100
G=4.3007*10.**4 #10**10Msun
datadir='/home/qyli/oPDFnew/data/apostle_MR/'
Globals.set_units(0.704*10.**10,0.704,1.)
Globals.get_units()

masstable=[3.868920866941207E-5,3.868920866941207E-5,3.960878906555616E-5,3.960878906555616E-5,3.8168845221394696E-5,3.8168845221394696E-5,3.8646087907235435E-5,3.8646087907235435E-5,3.876494719513764E-5,3.876494719513764E-5,3.8168845221394696E-5,3.8168845221394696E-5]

f=open('template.txt','w')

icount=35

filelist=[('S1fof1part',224.1709,129.82,5.8976),('S1fof2part',172.7433,59.4,8.8382),('S2fof1part',236.8058,153.03,8.8435),('S2fof2part',206.0765,100.85,6.9903),('S3fof1part',192.1194,81.72,12.776),('S3fof2part',173.7896,60.49,8.7318),('S4fof1part',233.5474,146.8,10.0735),('S4fof2part',194.9914,85.44,7.18),('S5fof1part',199.178,91.06,9.8334),('S5fof2part',189.7775,78.77,8.2011),('S6fof1part',197.9979,89.45,12.346),('S6fof2part',170.9422,57.56,11.1352)]

for base_file,rvir,M0,C0 in filelist:

  icount=icount+1
  if rvir*1.5<sep[icount-36+12]:
    up=rvir*1.5
  else:
    up=sep[icount-36+12]/2.
  print(up)
  xbin=np.logspace(np.log10(0.5),np.log10(up),nbin)
  xcen=(xbin[1:]+xbin[:-1])/2.
  vol=np.diff(np.hstack([0., xbin])**3)*np.pi*4./3.

  datafile=datadir+base_file+"1.hdf5"
  datafile = datafile.encode('utf8')
  Sample=Tracer(datafile)
  Sample.mP=masstable[icount-36]*10./0.704########
  countM1,tmp=np.histogram(Sample.data['r'], np.hstack([0.,xbin]))#dM
  countM1=countM1*Sample.mP
  countR1,tmp=np.histogram(Sample.data['r'], np.hstack([0.,xbin]), weights=Sample.mP/Sample.data['r'])#dM/R
  Sample.clean()

  datafile=datadir+base_file+'4.hdf5'
  datafile = datafile.encode('utf8')
  Sample=Tracer(datafile)
  countM4,tmp=np.histogram(Sample.data['r'], np.hstack([0.,xbin]), weights=Sample.mP*Sample.data['w'])#dM
  countR4,tmp=np.histogram(Sample.data['r'], np.hstack([0.,xbin]), weights=Sample.mP*Sample.data['w']/Sample.data['r'])#dM/R
  Sample.clean()

  datafile=datadir+base_file+'0.hdf5'
  datafile = datafile.encode('utf8')
  Sample=Tracer(datafile)
  countM0,tmp=np.histogram(Sample.data['r'], np.hstack([0.,xbin]), weights=Sample.mP*Sample.data['w'])#dM
  countR0,tmp=np.histogram(Sample.data['r'], np.hstack([0.,xbin]), weights=Sample.mP*Sample.data['w']/Sample.data['r'])#dM/R
  Sample.clean()

  datafile=datadir+base_file+'5.hdf5'
  datafile = datafile.encode('utf8')
  Sample=Tracer(datafile)
  countM5,tmp=np.histogram(Sample.data['r'], np.hstack([0.,xbin]), weights=Sample.mP*Sample.data['w'])#dM
  countR5,tmp=np.histogram(Sample.data['r'], np.hstack([0.,xbin]), weights=Sample.mP*Sample.data['w']/Sample.data['r'])#dM/R
  Sample.clean()

  countM=countM0+countM1+countM4+countM5
  countR=countR0+countR1+countR4+countR5

  density=countM/vol
  pot=countM.cumsum()/xbin+countR.sum()-countR.cumsum()
  pot*=G
  density_cum=countM.cumsum()/xbin**3/(4*np.pi/3)
  #pad with bin 0
  xbin=np.hstack([0., xbin])
  pot=np.hstack([countR.sum()*G, pot])
  density_cum=np.hstack([density_cum[0], density_cum])#?????
  mass=np.hstack([0.,countM.cumsum()])

  halo=Halo()
  halo.set_param([M0,C0])
  potNFW=-halo.pot(xbin)
  massNFW=halo.mass(xbin)

  print("x bin are",xbin)
  print("accumulative density is",density_cum)
  print("potential is ",pot)

  break

  # iref=-1 
  # fig = pl.figure(figsize=(7, 7))
  # gs = gridspec.GridSpec(2, 1, wspace=0.2, hspace=0.2)
  # ax1 = pl.subplot(gs[1])
  # ax0 = pl.subplot(gs[0])
  # pl.sca(ax0)
  # iref=np.abs(xbin-halo.Rs).argmin()
  # pl.plot(xbin, pot-pot[iref]+potNFW[iref], 'gx')
  # pl.plot(xbin, potNFW, 'k')
  # pl.loglog()
  # pl.xlabel('R')
  # pl.ylabel(r'$\psi$')
  # pl.legend(('Data','NFWfit'),loc='upper right')
  # pl.sca(ax1)
  # pl.plot(xbin, mass, 'gx')
  # pl.plot(xbin, massNFW, 'k')
  # pl.loglog()
  # pl.xlabel('R')
  # pl.ylabel(r'$M(<R)$')
  # pl.savefig('tmp_'+base_file+'.ps') #rasterize=True, dpi=300
  # pl.show()
  
  # print 'Profile template to be added to C/TemplateData.h:'
  # f.write('//'+base_file + str(icount)+'\n')
  # f.write('{{'+','.join(['{:f}'.format(i) for i in xbin])+'},'+'\n')
  # f.write('{'+','.join(['{:f}'.format(i) for i in pot])+'},'+'\n')
  # f.write('{'+','.join(['{:g}'.format(i) for i in density_cum])+'}},'+'\n')

  ## Now recompile and try the newly added template
  # tmpid=icount #id of the newly added template
  # xnew=np.logspace(-1,3,50)
  # tmphalo=Halo(halotype=HaloTypes.TMPMC, TMPid=tmpid)
  # tmphalo.set_param([M0,C0])
  # potNew=-tmphalo.pot(xnew)
  # tmphalo2=Halo(halotype=HaloTypes.TMPMC, TMPid=tmpid)
  # tmphalo2.set_param([2*M0,C0])
  # potNew2=-tmphalo2.pot(xnew)
  # pl.plot(xnew, potNew, 'ro')
  # pl.plot(xnew, potNew2, 'gs')
  # pl.plot(xbin, pot, 'k-')
  # pl.loglog()
  # pl.show()
  



