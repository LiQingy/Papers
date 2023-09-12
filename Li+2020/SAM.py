import numpy as np
import math
from astropy import constants as const

		

def SAM(n,nbin,xc,yc,zc,rv,region,model,hostid):

	if region <= 9:
		file = '/home/nifty2014/TheThreeHundred/SAMregions/%s/MDPL2_%s_z0.00_region000%s.txt' %(model,model,region)
	elif 9 < region <= 99:
		file = '/home/nifty2014/TheThreeHundred/SAMregions/%s/MDPL2_%s_z0.00_region00%s.txt' %(model,model,region)
	else:
		file = '/home/nifty2014/TheThreeHundred/SAMregions/%s/MDPL2_%s_z0.00_region0%s.txt' %(model,model,region)
	data = np.loadtxt(file)
	
	h = rv / n

	mass0 = np.zeros(n)
	mass4 = np.zeros(n)
	density4 = np.zeros(n)
	metal4 = np.zeros(n)
	age = np.zeros(n)
	rr4 = np.zeros(n)
	totaln4 = np.zeros(n)
	sfr_SAM = np.zeros(shape = (nbin,nbin))

	size = data.shape[0]	
	for j in range(size):
		rx = data[j][4] * 1e3 / 0.678
		ry = data[j][5] * 1e3 / 0.678
		rz = data[j][6] * 1e3 / 0.678
		rr = math.sqrt((rx - xc)**2 + (ry - yc)**2 + (rz - zc)**2)
		Mstar = (data[j][10] + data[j][11]) / 0.678
		Mgas = (data[j][12] + data[j][13] + data[j][14]) / 0.678
			
		if rr < rv and Mstar > 5 * 10**10:

			nn = int(rr / h)
			
			mass4[nn] = mass4[nn] + Mstar #mass
			mass0[nn] = mass0[nn] + Mgas

			massmetal = (data[j][25] + data[j][26]) / 0.678 / Mstar #Z = (MZstarDisk+MZstarSpheroid) / (MstarSpheroid+MstarDisk)
			metal4[nn] = metal4[nn] + massmetal * Mstar
			age[nn] = Mstar * data[j][19] + age[nn]

			#SFR
			if model == 'SAGE':
				limit = 1
			else:
				limit = 1e-9
			sfr0 = data[j][16] * limit / 0.678
			if sfr0 < 0.001:
				sfr0 = np.log10(0.001)
			else:
				sfr0 = np.log10(sfr0)
			
			nr = int(rr / ( rv / nbin))
			nsfr = int((sfr0 + 3) / 0.15)			
			sfr_SAM[nsfr][nr] += 1
 
			rr4[nn] = rr / rv * Mstar + rr4[nn]
			totaln4[nn] += 1
	#print(totaln4)
#---------------------------------------------------------------------------------------------
	ndv = np.zeros(n)
	for i in range(n):
		dv = 4. / 3. * math.pi * (pow((i + 1) * h, 3) - pow(i * h, 3))
		ndv[i] = dv
	
	density4 = mass4 / ndv
	totaln4 = totaln4 / ndv
	age = age / mass4
	starz = metal4 / mass4
	rr4 = rr4 / mass4

	return density4,age,starz,rr4,totaln4,sfr_SAM


def main(n,model,savedata):
	filerv = '/home/nifty2014/TheThreeHundred/playground/weiguang/MDPL2_Mass_snap_128-center-cluster.txt'
	filecen = '/home/nifty2014/TheThreeHundred/playground/qingyang/Center/%s_center.txt' %model
	datarv = np.loadtxt(filerv)
	datacen = np.loadtxt(filecen)

	density4_SAM = np.zeros(shape = (324,n))
	age_SAM = np.zeros(shape = (324,n))
	nbin = 40
	sfr_ALL = np.zeros(shape = (nbin,nbin))
	starz_SAM = np.zeros(shape = (324,n))
	rr4_SAM = np.zeros(shape = (324,n))
	totaln4_SAM = np.zeros(shape = (324,n))

	region = 1
	while region <= 324:

		rv = datarv[region - 1][6] / 0.678
		xc = datacen[region - 1][1] * 1e3 / 0.678
		yc = datacen[region - 1][2] * 1e3 / 0.678
		zc = datacen[region - 1][3] * 1e3 / 0.678
		hostid = datacen[region - 1][13]
		
		density4,age,starz,rr4,totaln4,sfr_SAM = SAM(n,nbin,xc,yc,zc,rv,region,model,hostid)


		density4_SAM[region - 1] = density4
		age_SAM[region - 1] = age
		starz_SAM[region - 1] = starz
		rr4_SAM[region - 1] = rr4
		totaln4_SAM[region - 1] = totaln4
		sfr_ALL += sfr_SAM

		print(region)
		region += 1
		

	if savedata == True:
		np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/SAMdata/SAM_%s_density4.txt' %model, density4_SAM)
		np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/SAMdata/SAM_%s_age.txt' %model, age_SAM)
		np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/SAMdata/SAM_%s_sfr.txt' %model, sfr_ALL)
		np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/SAMdata/SAM_%s_starz.txt' %model, starz_SAM)
		np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/SAMdata/SAM_%s_rr4.txt' %model, rr4_SAM)
		np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/SAMdata/SAM_%s_totaln4.txt' %model, totaln4_SAM)
	print(model)

main(20,'SAGE',savedata = True)
main(20,'SAG',savedata = True)
main(20,'Galacticus',savedata = True)
