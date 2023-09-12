import numpy as np 
import h5py
import readsnapsgl as readsnap
import gzip

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

def subpp(subid,hostid,datagal,datacen,i,pt,ptype):
	size = datagal.shape[0]
	rc = datacen[i][3:6]
	rvir = datacen[i][6]

	subhaloid = []	
	nskip = 2
	for j in range(size):
		rgal = datagal[j][5:8]
		rvgal = datagal[j][11]
		hid = datagal[j][1]
		galid = datagal[j][0]
		npp = int(datagal[j][4])

		rr = (np.sum((rc - rgal)**2))**0.5
		if ptype == 'ICL':
			if (0 < rr < (rvir + rvgal) and hid == hostid) or (galid == hostid): #exclude satellite structure particles
				ppdata = np.genfromtxt(subid, skip_header=nskip, max_rows=npp, dtype='i8')
				loc = np.where(ppdata[:,1] == pt)[0]
				ppdata = ppdata[loc,0]
				subhaloid.extend(ppdata)
				nskip = 1
			else:
				nskip += (npp + 1)
		else:
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
	if model == 'GXsub':
		filecen = '/home/nifty2014/TheThreeHundred/playground/weiguang/G3X_Mass_snap_128-center-cluster.txt'
	elif model == 'MDPL2sub':
		filecen = '/home/nifty2014/TheThreeHundred/playground/weiguang/MDPL2_Mass_snap_128-center-cluster.txt'
	datacen = np.loadtxt(filecen)
	
	
	for i in range(324):
		if model == 'GXsub':
			file = '/home/nifty2014/TheThreeHundred/simulation/GadgetX/NewMDCLUSTER_0%03d/snap_128' %(i + 1)
			filegal = '/home/nifty2014/TheThreeHundred/catalogues/AHF/GadgetX/NewMDCLUSTER_0%03d/GadgetX-NewMDCLUSTER_0%03d.snap_128.z0.000.AHF_halos' %(i + 1,i + 1)
			if i == 0:
				filesubid = '/home/nifty2014/TheThreeHundred/catalogues/AHF/GadgetX/NewMDCLUSTER_0%03d/GadgetX-NewMDCLUSTER_0%03d.snap_128.z0.000.AHF_particles' %(i+1,i+1)
				subid = open(filesubid, 'r')
			else:
				filesubid = '/home/nifty2014/TheThreeHundred/catalogues/AHF/GadgetX/NewMDCLUSTER_0%03d/GadgetX-NewMDCLUSTER_0%03d.snap_128.z0.000.AHF_particles.gz' %(i+1,i+1)
				subid = gzip.open(filesubid,'r')
			
			vx = h5py.File('/home/nifty2014/TheThreeHundred/playground/qingyang/oPDFnew/data/%s_1.2%s/GXsub_1.2%s_vx%s.h5' %(model,ptype,ptype,i + 1),'w')
			hostid = datacen[i][1]

		elif model == 'MDPL2sub':
			file = '/home/nifty2014/TheThreeHundred/simulation/MDPL2/NewMDCLUSTER_0%03d/snap_128' %(i + 1)
			filegal = '/home/nifty2014/TheThreeHundred/catalogues/AHF/MDPL2/NewMDCLUSTER_0%03d/MDPL2-NewMDCLUSTER_0%03d.snap_128.z0.000.AHF_halos' %(i + 1,i + 1)
			filesubid = '/home/nifty2014/TheThreeHundred/catalogues/AHF/MDPL2/NewMDCLUSTER_0%03d/MDPL2-NewMDCLUSTER_0%03d.snap_128.z0.000.AHF_particles' %(i+1,i+1)
			subid = open(filesubid,'r')
			
			vx = h5py.File('/home/nifty2014/TheThreeHundred/playground/qingyang/oPDFnew/data/MDPL2sub/MDPL2sub_%s_vx%s.h5' %(ptype,i + 1),'w')
			hostid = datacen[i][1]

		datagal = np.loadtxt(filegal)
		pos = readsnap.readsnapsgl(file , block = 'POS ', ptype=pt)
		vel = readsnap.readsnapsgl(file , block = 'VEL ', ptype=pt)
		#mass = readsnap.readsnapsgl(file , block = 'MASS', ptype=pt)
		gid = readsnap.readsnapsgl(file , block = 'ID  ', ptype=pt)
            
		print(gid.shape)
		subhaloid = subpp(subid,hostid,datagal,datacen,i,pt,ptype)

		rc = datacen[i][3:6]
		rv = datacen[i][6]
		loc = np.where(datagal[:,0] == hostid)[0]
		vc = datagal[loc,8:11][0]

		rr = np.sqrt(np.sum((pos - rc)**2, axis = 1))
		locr = np.where(rr < rv)[0]
		pos = pos[locr]
		vel = vel[locr]
		gid = gid[locr]
		#mass = mass[locr]

		idwell = list( set(gid).difference( set(subhaloid)) )   
		idwell = np.array(idwell)
		locid = indexed(gid, idwell)

		v = np.array(vel[locid] - vc,dtype = 'f4')
		x = np.array(pos[locid] - rc,dtype = 'f4')
		#m0 = np.array(mass[locid],dtype = 'f4')

		vx['v'] = v
		vx['x'] = x
		#vx['Partmass'] = m0
		vx.close()

		print(i)
readvx(model = 'GXsub',pt = 1,ptype = 'DM')
