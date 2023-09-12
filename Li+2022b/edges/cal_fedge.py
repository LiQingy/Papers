import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy.coordinates import SkyCoord
import h5py
import os


h5fname="./CLAUDS_HSC_iband_random200_gals_masked.hdf5"
if not os.path.isfile(h5fname):
	fh5w = h5py.File(h5fname, "w")
	fh5w.close()
fh5w = h5py.File(h5fname, "r+")

dmasked = np.loadtxt('./CLAUDS_HSC_iband_random200_gals_masked')

fh5w.create_dataset('groupid', data = dmasked[:,0])
fh5w.create_dataset('ra', data = dmasked[:,1])
fh5w.create_dataset('dec', data = dmasked[:,2])
fh5w.create_dataset('masked_flag', data = dmasked[:,3])
fh5w.close()


# Ngroup = 2232134

# fedge = np.zeros((Ngroup,2))

# iseq = 0
# for i in range(Ngroup):
# 	iend = iseq + 200
# 	masked_flag = dmasked[iseq:iend,3]
# 	Nremain = np.sum(masked_flag == 1)

# 	fedge[i,0] = i+1
# 	fedge[i,1] = Nremain / 200

# 	iseq = iend

# np.savetxt('./CLAUDS_HSC_iband_grpfedge', fedge, header = '1: group id; 2: fedge (Nremain/200)', fmt = '%d %.6f')   
