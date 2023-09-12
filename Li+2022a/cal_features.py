import numpy as np 
import matplotlib
import glob


'''
Age-weighted mass accretion rate
'''
def masshistroy(hostid, nn, k):
	snapall = np.linspace(128,129-nn,nn,dtype = 'i8')

	mass_accrete = np.zeros(nn)
	redshift = np.zeros(nn)
	snapn = 128
	while snapn >= 129 - nn:
		filegal0 = glob.glob('/home/nifty2014/TheThreeHundred/catalogues/AHF/GadgetX/NewMDCLUSTER_0%03d/GadgetX-NewMDCLUSTER_0%03d.snap_%03d.z*.AHF_halos' %(k+1,k+1,snapn))
		datagal0 = np.loadtxt(filegal0[0])
		filemaintree = glob.glob('/home/nifty2014/TheThreeHundred/catalogues/AHF/GadgetX/NewMDCLUSTER_0%03d/GadgetX-NewMDCLUSTER_0%03d.snap_%03d.*.AHF_mtree_idx'%(k+1,k+1,snapn))
		mtree = np.loadtxt(filemaintree[0])
		filegal1 = glob.glob('/home/nifty2014/TheThreeHundred/catalogues/AHF/GadgetX/NewMDCLUSTER_0%03d/GadgetX-NewMDCLUSTER_0%03d.snap_%03d.z*.AHF_halos' %(k+1,k+1,snapn-1))
		
		if filegal1 == []:
			filegal1 = glob.glob('/home/nifty2014/TheThreeHundred/catalogues/AHF/GadgetX/NewMDCLUSTER_0%03d/GadgetX-NewMDCLUSTER_0%03d.snap_%03d.*.AHF_halos'%(k+1,k+1,snapn-2))
			skip = 2
		else:
			skip = 1
		datagal1 = np.loadtxt(filegal1[0])
		redshift[128-snapn] = filegal0[0].split('z')[1].split('.A')[0]

		loclast = np.where(mtree[:,0] == hostid)[0]
		hostlast = mtree[loclast,1]

		loc0 = np.where(datagal0[:,0] == hostid)[0]
		loc1 = np.where(datagal1[:,0] == hostlast)[0]
		mass_accrete[128-snapn] = datagal0[loc0,3] - datagal1[loc1,3]

		hostid = hostlast
		if skip == 2:
			redshift[128-snapn+1] = -1
		snapn = int(snapn - skip)

	return mass_accrete,redshift

def main(nn):
	filecen = '/home/nifty2014/TheThreeHundred/catalogues/AHF/complete_sample/G3X_Mass_snap_128-center-cluster.txt'
	datacen = np.loadtxt(filecen)

	massall = np.zeros(shape = (324, nn))
	zall = np.zeros(shape = (324, nn))
	for k in range(324):
		hostid = np.array([datacen[k][1]], dtype = 'i8')
		massall[k],zall[k] = masshistroy(hostid, nn, k)
		print(k)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/oPDFnew2/data/Mhistory/Mhistory.txt',massall)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/oPDFnew2/data/Mhistory/zall.txt',zall)
# nn: the number of snapshot
main(nn = 30)

def masshistroy():
	filemass = './Mhistory.txt'
	filez = './zall.txt'
	amass = np.loadtxt(filemass)[:,:15]
	redz = np.loadtxt(filez)[:,:16]

	loc0 = np.where(redz == -1)
	redz[loc0] = 0
	ageall = cosmo.age(redz).value
	ageall[loc0] = 0

	deltaage = np.diff(-ageall)
	age = ageall[:,:15]

	amassfunc = np.sum(age * amass / deltaage, axis = 1) / np.sum(age, axis = 1)
	np.savetxt('./mar.txt',amassfunc)
	return amassfunc
#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================
'''
ICL and BCG stellar mass fraction
'''

def indexed(x, y, missing='raise', return_missing=False):
    x, y = np.asarray(x), np.asarray(y)

    x_index = np.argsort(x)
    y_index_sorted = np.searchsorted(x[x_index], y, side='left')
    index = np.take(x_index, y_index_sorted, mode="clip")

    if missing != 'ignore' or return_missing:
        invalid = x[index] != y

    if missing != 'ignore':
        if missing == 'raise':
            if np.any(invalid):
                raise ValueError('Not all elements in `y` are present in `x`')
        elif missing == 'mask':
            index = np.ma.array(index, mask=invalid)
        else:
            index[invalid] = missing

    if return_missing:
        return index, invalid
    else:
        return index

def subpp(model,subid,hostid,datagal,datacen,i,pt):
	size = datagal.shape[0]
	rc = datacen[i][3:6]
	rvir = datacen[i][6]

	subhaloid = []	
	nskip = 2
	for j in range(size):
		rgal = datagal[j][5:8]
		rvgal = datagal[j][11]
		hid = datagal[j][1]
		npp = int(datagal[j][4])

		rr = (np.sum((rc - rgal)**2))**0.5
		if 0 < rr < (rvir + rvgal) and hid == hostid: #exclude satellite structure particles
			ppdata = np.genfromtxt(subid, skip_header=nskip, max_rows=npp, dtype='i8')
			loc = np.where(ppdata[:,1] == pt)[0]
			ppdata = ppdata[loc,0]
			subhaloid.extend(ppdata)
			nskip = 1
		else:
			nskip += (npp + 1)

	subhaloid = np.array(subhaloid)
	print(subhaloid.shape)
	return subhaloid
		
def readvx(model, pt, ptype):
	filecen = '/home/nifty2014/TheThreeHundred/playground/weiguang/G3X_Mass_snap_128-center-cluster.txt'
	datacen = np.loadtxt(filecen)
	
	ICLfrac = np.zeros(324)
	for i in range(324):
		file = '/home/nifty2014/TheThreeHundred/simulation/GadgetX/NewMDCLUSTER_0%03d/snap_128' %(i + 1)
		filegal = '/home/nifty2014/TheThreeHundred/catalogues/AHF/GadgetX/NewMDCLUSTER_0%03d/GadgetX-NewMDCLUSTER_0%03d.snap_128.z0.000.AHF_halos' %(i + 1,i + 1)
		if i == 0:
			filesubid = '/home/nifty2014/TheThreeHundred/catalogues/AHF/GadgetX/NewMDCLUSTER_0%03d/GadgetX-NewMDCLUSTER_0%03d.snap_128.z0.000.AHF_particles' %(i+1,i+1)
			subid = open(filesubid, 'r')
		else:
			filesubid = '/home/nifty2014/TheThreeHundred/catalogues/AHF/GadgetX/NewMDCLUSTER_0%03d/GadgetX-NewMDCLUSTER_0%03d.snap_128.z0.000.AHF_particles.gz' %(i+1,i+1)
			subid = gzip.open(filesubid,'r')
		hostid = datacen[i][1]

		datagal = np.loadtxt(filegal)
		pos = readsnap.readsnapsgl(file , block = 'POS ', ptype=pt)
		mass = readsnap.readsnapsgl(file, block = 'MASS', ptype=pt)
		gid = readsnap.readsnapsgl(file , block = 'ID  ', ptype=pt)
		subhaloid = subpp(model,subid,hostid,datagal,datacen,i,pt)

		rc = datacen[i][3:6]
		rv = datacen[i][6]
		loc = np.where(datagal[:,0] == hostid)[0]

		rr = np.sqrt(np.sum((pos - rc)**2, axis = 1))
		locr = np.where(rr < rv)[0]
		pos = pos[locr]
		mass = mass[locr]
		gid = gid[locr]

		idwell = list( set(gid).difference( set(subhaloid) ) )   
		idwell = np.array(idwell)
		locid = indexed(gid, idwell)

		ICLfrac[i] = np.sum(mass[locid]) / np.sum(mass)
		print(i)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/oPDFnew2/data/ICLfmass/ICLfrac.txt', ICLfrac)
readvx(model = 'Gadget-X', pt = 4, ptype = 'star')

#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================
























