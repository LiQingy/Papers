The radial bins are given by rbins=np.logspace(np.log10(0.001*r200), np.log10(1.5*r200), num=151)
The profiles are calculated from snapshot_128 (z = 0) to snapshot_25 (z = 8.79).
The mean atomic weight is 0.6125.

The head information in hdf5 files (attributions) for each cluster include:
'redshift', 'snapshot', 'center_x', 'center_y', 'center_z', 'r200', 'r500', 'M200', 'M500'. For these quantities, distance unit in kpc/h, mass unit in M_sun/h, velocity unit in km/s.

The gas particles are selected with T > 10^6 K and single density < 2.88*10^6 M_sun/kpc^3.

Profile descriptions and units:
'Rbin' -- mean position in each bin interval -- kpc/h
'NuminBin' -- particle numbers in each bin (gas, star and total particles in each row respectively) -- None
'Gdens' -- gas mass density -- M_sun*h^2*kpc^-3; 
'MWTemp' -- mass weighted temperature -- keV; 
'Pressure' -- gas pressure -- keV*(cm/h)^-3;
'Potential' -- total gravitational potential -- (kpc/s)^2*h
'MWMetal' -- mass-weighted gas metallicity -- Zsun (Zsun==0.0134 taken from Asplund+2009);
'Eledens' -- electron number density -- (cm/h)^-3;
'Entropy' -- gas entropy -- keV*(cm/h)^2;
'Veldisp' -- velocity dispersion including Hubble flow (gas, subhaloes and total particles in each row respectively) -- km/s; 
'Stellardens' -- stellar mass density -- M_sun*h^2*kpc^-3;
'Totdens' -- total mass density -- M_sun*h^2*kpc^-3;
