import numpy as np 
import h5py
import readsnapsgl as readsnap
import gzip
import glob

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

def readsubid(sn, crn, pt, ptype):
    '''
    read particle id in subhalos for AHF files

    Parameters:
    -----------
    sn: snapshot 
    crn: cluster id
    pt: type of particles
    ptype: 'TOT': select all particles in host halos; 'SUB': select all particles in subhalos

    Return:
    -------
	subhaloid: subhalo particle id 
	outsubid: outsubhalo particle id
	locsubid: index of subhalo particle id
	locoutsubid: index of out subhalo particle id
    
    Note also:
    ----------
    1. a few particles belong to not only one subhalos
    2. a few particles can still exist in the host halo but out of the content of host halo particles due to the unbounding problem
    '''

    progenIDs=np.loadtxt("/home2/weiguang/Project-300-Clusters/Halo_mass_function_mass-difference/GIZMO/Progenitor-IDs-for-center-cluster.txt",dtype=np.int64)
    Xspath = "/home2/weiguang/data7/Gizmo-Simba/"

    cn = 'NewMDCLUSTER_%04d'%crn
    #ReginIDs HIDs  HosthaloID Mvir(4) Xc(5)   Yc(6)   Zc(7)  Rvir(8) fMhires(38) cNFW (42) Mgas200 M*200 M500(13)  R500(14) fgas500 f*500	
    tmpd = np.load('/home2/weiguang/Project-300-Clusters/Halo_mass_function_mass-difference/GIZMO/GS_Mass_snap_%03dinfo.npy' %sn)
    ids = np.where((np.int32(tmpd[:,0])==crn) & (np.int64(tmpd[:,1]) == progenIDs[crn-1, sn]))[0]

    filename = Xspath+cn+'/'+'snap_%03d'%sn+'.hdf5'
    f=h5py.File(filename, 'r')
    gid = f['PartType'+str(pt)+'/ParticleIDs'][:]
    f.close()

    filegal = glob.glob('/data4/niftydata/TheThreeHundred/data/catalogues/AHF/GIZMO/%s/GIZMO-%s.snap_%03d.z*.AHF_halos' %(cn,cn,sn))[0]
    filesubid = glob.glob('/data4/niftydata/TheThreeHundred/data/catalogues/AHF/GIZMO/%s/GIZMO-%s.snap_%03d.z*.AHF_particles.gz' %(cn,cn,sn))[0]

    datagal = np.loadtxt(filegal)
    subid = gzip.open(filesubid,'r')

    hostid = progenIDs[crn-1, sn]
    #--------------------------------------------------------------------
    size = datagal.shape[0]

    subhaloid = []	
    nskip = 2
    idxsub = np.where(datagal[:,1] == hostid)[0]
    for j in range(size):
        hid = datagal[j][1] #host halo id
        galid = datagal[j][0]
        npp = int(datagal[j][4])

        if ptype == 'TOT':
            if galid == hostid: #find particles in host halo
                ppdata = np.genfromtxt(subid, skip_header=nskip, max_rows=npp, dtype='i8')
                loc = np.where(ppdata[:,1] == pt)[0]
                ppdata = ppdata[loc,0]
                subhaloid.extend(ppdata)
                break
        else:
            if hid == hostid: #find particles in subhalo structures 
                ppdata = np.genfromtxt(subid, skip_header=nskip, max_rows=npp, dtype='i8')
                loc = np.where(ppdata[:,1] == pt)[0]
                ppdata = ppdata[loc,0]
                subhaloid.extend(ppdata)
                nskip = 1
            else:
                nskip += (npp + 1)

    subhaloid = np.array(subhaloid, dtype = np.int64)

    #exclude subhalo particle id in total particle id
    outsubid = list(set(gid).difference(set(subhaloid)))
    outsubid = np.array(outsubid, dtype = np.int64)
    locoutsubid = indexed(gid, outsubid)
    locsubid = indexed(gid, subhaloid)

    return subhaloid, outsubid, locsubid, locoutsubid

# subhaloid, outsubid, locsubid, locoutsubid = readsubid(sn = 128, crn = 5, pt = 0, ptype = 'TOT')
