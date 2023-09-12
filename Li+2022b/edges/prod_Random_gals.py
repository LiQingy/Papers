import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy.coordinates import SkyCoord
cosmo = FlatLambdaCDM(H0=67.4, Om0=0.315)
import random

def cal_rp0(ag_ra,ag_dec,c_ra,c_dec,c_redz):
    from astropy.cosmology import FlatLambdaCDM
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    cosmo = FlatLambdaCDM(H0=67.4, Om0=0.315)
    
    dr_clu = cosmo.comoving_distance(c_redz).value * 0.674 #Mpc/h

    cgal = SkyCoord(ra=ag_ra*u.degree, dec=ag_dec*u.degree, distance=dr_clu)
    cclu = SkyCoord(ra=c_ra*u.degree, dec=c_dec*u.degree, distance=dr_clu)
    cgal_x = cgal.cartesian.x
    cgal_y = cgal.cartesian.y
    cgal_z = cgal.cartesian.z

    cclu_x = cclu.cartesian.x
    cclu_y = cclu.cartesian.y
    cclu_z = cclu.cartesian.z

    l = np.array([cgal_x+cclu_x, cgal_y+cclu_y, cgal_z+cclu_z]).T / 2
    s = np.array([cclu_x - cgal_x, cclu_y - cgal_y, cclu_z - cgal_z]).T
    r_pi = np.sum(l*s,axis = 1) / np.sqrt(np.sum(l**2, axis = 1)) 
    r_p = np.sqrt(np.sum(s**2,axis = 1) - r_pi**2)
#     r1 = np.sqrt(cgal_x**2 + cgal_y**2 + cgal_z**2)
#     r2 = np.sqrt(cclu_x**2 + cclu_y**2 + cclu_z**2)
#     r_p = (cgal_x - cclu_x)**2 + (cgal_y - cclu_y)**2 + (cgal_z - cclu_z)**2 - (r1 - r2)**2 
    
    return r_p

def cal_sep(clu_r180, clu_redz):
    from astropy.cosmology import FlatLambdaCDM
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    cosmo = FlatLambdaCDM(H0=67.4, Om0=0.315)
    
    d_A = cosmo.comoving_distance(z=clu_redz) 
    r180u = (clu_r180) * u.Mpc / 0.674
    theta_radian = (r180u / d_A).to(u.degree, u.dimensionless_angles()).value # unit is Mpc only now
    
    return theta_radian 

def prod_random_gals():

	table_grp = []
	for cid in range(1,Ngrp+1): 
		clu_ra = dgroup[cid-1, 2]
		clu_dec = dgroup[cid-1, 3]
		clu_redz = dgroup[cid-1, 4]
		clu_Mh = dgroup[cid-1, -2]
		# clu_r180 = read_r180(cid)

		clu_r180 = 0.781 * (10**clu_Mh / 0.315 / 1e14)**(1/3)

		theta_radian = cal_sep(clu_r180, clu_redz) #calculate the maxmium separation

		N_tot = 0
		ra_ag = []
		dec_ag = []
		while N_tot < 200:
		    ra_random = np.random.random(3000) - 0.5
		    dec_random = np.random.random(3000) - 0.5
		    ra_r = ra_random * 2 * (theta_radian+0.2/(1+clu_redz*1.5)) + clu_ra
		    dec_r = dec_random * 2 * (theta_radian+0.15/(1+clu_redz*1.5)) + clu_dec
		    
		    d_r = cal_rp0(ra_r,dec_r,clu_ra,clu_dec,clu_redz)
		    idx_sat = d_r <= clu_r180
		    ra_sat = ra_r[idx_sat]
		    dec_sat = dec_r[idx_sat]
		    N_sat = ra_sat.shape[0]

		    ra_ag.extend(ra_sat)
		    dec_ag.extend(dec_sat)
		    N_tot += N_sat
		    
		idx_200 = np.int64(random.sample(range(0, N_tot), 200))
		ra_ag = np.array(ra_ag)[idx_200]
		dec_ag = np.array(dec_ag)[idx_200]

		cluid = np.tile(cid, ra_ag.shape[0])
		one_coor = np.vstack((cluid,ra_ag,dec_ag)).T
		table_grp.extend(one_coor)

		if cid%100000 == 0:
			print(cid)

	table_grp = np.array(table_grp)
	np.savetxt('./CLAUDS_HSC_iband_random200_gals', table_grp, header = '1: group id; 2. random ra; 3. random dec', fmt = '%d  %.6f %.6f')


dgroup = np.loadtxt('../odata/CLAUDS_HSC_iband_group')
digal = np.loadtxt('../odata/CLAUDS_HSC_iband_igal')
d2 = np.loadtxt('../odata/iCLAUDS_HSC_iband_2')


Ngrp = dgroup.shape[0]

prod_random_gals()

