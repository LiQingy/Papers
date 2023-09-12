
import numpy as np
import math


#-------------------------------------------------------------------------------------------------------------------
global pi
pi = 3.1415926


#----------------------------------------------------------------------------------------------
def AHF1(n,nbin,xc,yc,zc,rv,region,model):
	if model == 'Gadget-X':
		if region <= 9:
			file = '/home/nifty2014/TheThreeHundred/playground/weiguang/GadgetX-new/NewMDCLUSTER_000%s_snap_128.txt' %region
		elif 9 < region <= 99:
			file = '/home/nifty2014/TheThreeHundred/playground/weiguang/GadgetX-new/NewMDCLUSTER_00%s_snap_128.txt' %region
		else:
			file = '/home/nifty2014/TheThreeHundred/playground/weiguang/GadgetX-new/NewMDCLUSTER_0%s_snap_128.txt' %region
	if model == 'Gadget-MUSIC':
		if region <= 9:
			file = '/home/nifty2014/TheThreeHundred/playground/weiguang/GadgetMUSIC-new/GadgetMUSIC/NewMDCLUSTER_000%s_snap_017.txt' %region
		elif 9 < region <= 99:
			file = '/home/nifty2014/TheThreeHundred/playground/weiguang/GadgetMUSIC-new/GadgetMUSIC/NewMDCLUSTER_00%s_snap_017.txt' %region
		else:
			file = '/home/nifty2014/TheThreeHundred/playground/weiguang/GadgetMUSIC-new/GadgetMUSIC/NewMDCLUSTER_0%s_snap_017.txt' %region
	data = np.loadtxt(file)

	h = rv / n

	mass4_age = np.zeros(n)
	age = np.zeros(n)
	sfr_AHF = np.zeros(shape = (nbin,nbin))
	rrage = np.zeros(n)

	size = data.shape[0]
	for k in range(size):
		rx = data[k][5] / 0.678
		ry = data[k][6] / 0.678
		rz = data[k][7] / 0.678
		rr = math.sqrt((rx - xc)**2 + (ry - yc)**2 + (rz - zc)**2)
		if  rr < rv:
			nn = int(rr / h)
			# #AHF HID [1] Mvir[2] Rvir[3] Mgas[4] M*[5] XYZ[678] AGE[9 Gyr] SFR[10 M*/yr]
			if data[k][8] > 0 and data[k][4] / 0.678 > 5 * 10**10:
				mass4_age[nn] = data[k][4] / 0.678 + mass4_age[nn]
				age[nn] = data[k][4] / 0.678 * data[k][8] + age[nn]
				rrage[nn] = rr / rv * data[k][4] / 0.678+ rrage[nn]

			if data[k][4] / 0.678 > 5 * 10**10:
				sfr = data[k][9]
				if sfr < 0.001:
					sfr = np.log10(0.001)
				else:
					sfr = np.log10(sfr)
				nr = int(rr / ( rv / nbin))
				nsfr = int((sfr + 3) / 0.15)
				if nsfr < 40:
					sfr_AHF[nsfr][nr] += 1
				else:
					print('#&')
#--------------------------------------------------------------------------------------
	'''
	for i in range(1,n):
		if nub_age[i] < 5:
			for j in range(i):
				nbin = i - j - 1
				if nub_age[nbin] > 0:
					age[nbin] = age[nbin] + age[i]
					mass4_age[nbin] = mass4_age[nbin] + mass4_age[i]
					rrage[nbin] = rrage[nbin] + rrage[i]
					nub_age[nbin] = nub_age[nbin] + nub_age[i]
					break

			age[i] = 0
			mass4_age[i] = 0
			rrage[i] = 0
			nub_age[i] = 0

		if nub_sfr[i] < 5:
			for k in range(i):
				nbin = i - k - 1
				if nub_sfr[nbin] > 0:
					sfr[nbin] = sfr[nbin] + sfr[i]
					mass4_sfr[nbin] = mass4_sfr[nbin] + mass4_sfr[i]
					rrsfr[nbin] = rrsfr[nbin] + rrsfr[i]
					nub_sfr[nbin] = nub_sfr[nbin] + nub_sfr[i]
					break

			sfr[i] = 0
			mass4_sfr[i] = 0
			rrsfr[i] = 0
			nub_sfr[i] = 0
	'''
#--------------------------------------------------------------------------------------
	age = age / mass4_age
	rrage = rrage / mass4_age

	return age,rrage,sfr_AHF


def AHF2(n,xc,yc,zc,rv,region,model):

	if model == 'Gadget-X':
		if region <= 9:
			file = '/home/nifty2014/TheThreeHundred/catalogues/AHF/GadgetX/NewMDCLUSTER_000%s/GadgetX-NewMDCLUSTER_000%s.snap_128.z0.000.AHF_halos' %(region,region)
		elif 9 < region <= 99:
			file = '/home/nifty2014/TheThreeHundred/catalogues/AHF/GadgetX/NewMDCLUSTER_00%s/GadgetX-NewMDCLUSTER_00%s.snap_128.z0.000.AHF_halos' %(region,region)
		else:
			file = '/home/nifty2014/TheThreeHundred/catalogues/AHF/GadgetX/NewMDCLUSTER_0%s/GadgetX-NewMDCLUSTER_0%s.snap_128.z0.000.AHF_halos' %(region,region)
	if model == 'Gadget-MUSIC':
		if region <= 9:
			file = '/home/nifty2014/TheThreeHundred/catalogues/AHF/GadgetMUSIC/NewMDCLUSTER_000%s/GadgetMUSIC-NewMDCLUSTER_000%s.z0.000.AHF_halos' %(region,region)
		elif 9 < region <= 99:
			file = '/home/nifty2014/TheThreeHundred/catalogues/AHF/GadgetMUSIC/NewMDCLUSTER_00%s/GadgetMUSIC-NewMDCLUSTER_00%s.z0.000.AHF_halos' %(region,region)
		else:
			file = '/home/nifty2014/TheThreeHundred/catalogues/AHF/GadgetMUSIC/NewMDCLUSTER_0%s/GadgetMUSIC-NewMDCLUSTER_0%s.z0.000.AHF_halos' %(region,region)
	data = np.loadtxt(file)

	h = rv / n

	mass4 = np.zeros(n)
	density4 = np.zeros(n)
	metal4 = np.zeros(n)
	totaln4 = np.zeros(n)
	rr4 = np.zeros(n)

	size = data.shape[0]	
	for j in range(size):
		rx = data[j][5] / 0.678
		ry = data[j][6] / 0.678
		rz = data[j][7] / 0.678
		rr = math.sqrt((rx - xc)**2 + (ry - yc)**2 + (rz - zc)**2)
		if rr <= rv and data[j][64] / 0.678 >5 * 10**10:

			nn = int(rr / h)

			mass4[nn] = data[j][64] / 0.678 + mass4[nn] #stellar mass; M_star(65)
			metal4[nn] = data[j][64] / 0.678 * data[j][87] + metal4[nn] #stellar metal; mean_z_star(88)
			rr4[nn] = rr / rv * data[j][64] / 0.678 + rr4[nn] 
			totaln4[nn] += 1
#-----------------------------------------------------------------------------------------------------
	ndv = np.zeros(n)
	for i in range(n):
		dv = 4. / 3 * math.pi * (pow((i + 1) * h,3) - pow(i * h,3))
		ndv[i] = dv
	'''
	for i in range(1,n):
		if totaln4[i] < 5:
			for k in range(i):
				nbin = i - k - 1
				if totaln4[nbin] > 0:
					mass4[nbin] = mass4[nbin] + mass4[i]
					metal4[nbin] = metal4[nbin] + metal4[i]
					rr4[nbin] = rr4[nbin] + rr4[i]
					totaln4[nbin] = totaln4[nbin] + totaln4[i]
					ndv[nbin] = ndv[nbin] + ndv[i]
					break

			mass4[i] = 0
			metal4[i] = 0
			rr4[i] = 0
			totaln4[i] = 0
			ndv[i] = 0
	
#-----------------------------------------------------------------------------------------------------
	'''
	'''				
	for k in range(n):
		dv = 4. / 3 * math.pi *(pow((k + 1) * h,3) - pow(k * h,3))
		density4[k] = mass4[k] / dv
		totaln4[k] = totaln4[k] / dv

	'''
	density4 = mass4 / ndv
	totaln4 = totaln4 / ndv
	starz = metal4 / mass4
	rr4 = rr4 / mass4
	return density4,starz,rr4,totaln4

def main(n,model,savetxt):
	if model == 'Gadget-X':
		filecen = '/home/nifty2014/TheThreeHundred/playground/weiguang/G3X_Mass_snap_128-center-cluster.txt'
		name = 'GX'
	if model == 'Gadget-MUSIC':
		filecen = '/home/nifty2014/TheThreeHundred/playground/weiguang/Music_Mass_snap_017-center-cluster.txt'
		name = 'GM'
	datacen = np.loadtxt(filecen)

	density4_AHF = np.zeros(shape = (324,n))
	age_AHF = np.zeros(shape = (324,n))
	starz_AHF = np.zeros(shape = (324,n))
	rr4_AHF = np.zeros(shape = (324,n))
	rrage_AHF = np.zeros(shape = (324,n))
	totaln4_AHF = np.zeros(shape = (324,n))
	nbin = 40
	sfr_ALL = np.zeros(shape = (nbin,nbin))

	region = 1
	while  region <= 324:


		rv = datacen[region - 1][6] / 0.678
		xc = datacen[region - 1][3] / 0.678
		yc = datacen[region - 1][4] / 0.678
		zc = datacen[region - 1][5] / 0.678
		
		age,rrage,sfr_AHF = AHF1(n,nbin,xc,yc,zc,rv,region,model)
		density4,starz,rr4,totaln4 = AHF2(n,xc,yc,zc,rv,region,model)

		density4_AHF[region - 1] = density4
		age_AHF[region - 1] = age
		starz_AHF[region - 1] = starz
		rr4_AHF[region - 1] = rr4
		rrage_AHF[region - 1] = rrage
		totaln4_AHF[region - 1] = totaln4
		sfr_ALL += sfr_AHF

		print(region)
		region += 1
	if savetxt == True:
		np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/AHFdata/%s_AHF_density4.txt' %name,density4_AHF)
		np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/AHFdata/%s_AHF_age.txt' %name,age_AHF)
		np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/AHFdata/%s_AHF_starz.txt' %name,starz_AHF)
		np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/AHFdata/%s_AHF_sfr.txt' %name,sfr_ALL)
		np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/AHFdata/%s_AHF_rr4.txt' %name,rr4_AHF)
		np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/AHFdata/%s_AHF_rrage.txt' %name,rrage_AHF)
		np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/AHFdata/%s_AHF_totaln4.txt' %name,totaln4_AHF)
	print(model)

main(20, 'Gadget-MUSIC', savetxt = True)
main(20, 'Gadget-X', savetxt = True)
