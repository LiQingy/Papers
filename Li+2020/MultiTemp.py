import numpy as np 
from struct import unpack
from os import fstat
from astropy import constants as const
import math
from astropy.cosmology import Planck15 as cosmo


global pi
pi = 3.1415926

#------------------------------------------------------------------------------------

#read the number of particles
def read_npart(filename):
	file = open(filename, 'rb')	
	space_8 = unpack('<i',file.read(4))
	bname = file.read(4).decode('ascii')
	bsize = unpack('<i',file.read(4))
	space8 = file.read(8)
	npart = unpack('<iiiiii', file.read(4 * 6))
	file.close()
	return npart

#read the datas of block(except 'MASS')
def read_datas(filename,block, dtype, column, dt, file):
	npart = read_npart(filename)
	if block == 'AGE ':
		data = file.read(npart[4] * dt.itemsize * column)
	elif block == 'ZTOT' or block == 'Z   ' or block == 'ZS  ':
		if dtype == 0:
			data = file.read(npart[0] * 4 * column)
		else:
			file.seek(file.tell() + npart[0] * 4 * column)
			data = file.read(npart[4] * 4 * column)
	else:	
		if dtype == 6:
			data = file.read(np.sum(npart[:dtype]) * dt.itemsize * column)
			nsum = np.sum(npart[:])
			arr = np.ndarray(shape = (np.int32(nsum), column), dtype = dt, buffer = data)
		else:
			if dtype == 0:
				data =  file.read(npart[0] * dt.itemsize * column)
			else:
				file.seek(file.tell() + np.sum(npart[:dtype]) * dt.itemsize * column)
				data = file.read(npart[dtype] * dt.itemsize * column)

	if column != 1 and dtype != 6:
		arr = np.ndarray(shape = (np.int32(npart[dtype]), column), dtype = dt, buffer = data)
	elif column == 1 and dtype != 6:
		arr = np.ndarray(shape = (np.int32(npart[dtype])), dtype = dt, buffer = data)
	
	return arr

#find the location of block
def read_all(filename, block, dtype, column ,dt): 
	file = open(filename, 'rb')
	allsize = fstat(file.fileno()).st_size
	
	size0 = 24 + 256
	file.seek(size0)
	while file.tell() < allsize:
		space4 = file.read(4)
		bname = file.read(4).decode('ascii')
		bsize = unpack('<i',file.read(4))[0] 
		space8 = file.read(8)
		if (block != bname):
			size0 =  bsize + size0 +16
			file.seek(size0)
		else:
			return read_datas(filename,block, dtype, column, dt, file)
		pass
	file.close()
	return None

#the mass of particles
def mass_type(filename,ntype):
	npart = read_npart(filename)
	file = open(filename, 'rb')
	file.seek(20 + 24)
	massarr = unpack('<dddddd', file.read(8 * 6))

	if massarr[ntype] != 0:
		mas = massarr[ntype]
		mass = np.ones(npart[ntype])
		mass = mass * mas
	elif npart[ntype] == 0:
		mass = 0 
	else:
		size1 = 280
		for i in range(3):
			file.seek(size1)
			file.read(4)
			bname = file.read(4).decode('ascii')
			bsize = unpack('<i',file.read(4))[0] 
			size1 = size1 + bsize + 16
		locntype = 0
		for j in range(ntype):
			if massarr[j] == 0:
				locntype = locntype + npart[j] * 4		
		file.seek(size1 + 20 + locntype)
		data = file.read(npart[ntype] * 4)
		mass = np.ndarray(shape = np.int32(npart[ntype]), dtype = 'float32', buffer = data)
		
	file.close()
	return mass

#the temperture of gas 
def caltemp(filename):
	file = open(filename, 'rb')
	file.read(20 + 24 + 48)
	time = unpack('<d',file.read(8))
	temp = read_all(filename,'U   ',0,1,np.dtype('float32'))
	xH = 0.76
	yhelium = (1. - xH) / (4 * xH)
	NE = read_all(filename,'NE  ',0,1,np.dtype('float32'))
	mean_mol_weight = (1. + 4. * yhelium) / (1. + yhelium + NE)
	v_unit = 1.0e5 * np.sqrt(time)
	prtn = 1.67373522381e-24
	bk = 1.3806488e-16
	temp = temp * (5. / 3 - 1) * v_unit**2 * prtn * mean_mol_weight / bk

	temp = temp * const.k_B.value / const.e.value * 1e-3
	file.close()
	return temp

#the stars age
def calage(filename):
	a = read_all(filename,'AGE ', 4, 1, np.dtype('float32'))
	z = 1. / a - 1
	age = z
	age = cosmo.age(z).value
	return age




#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#find particles within virial radius
def distribution(xc,yc,zc,rv,filename,nn):# h size of interval

	#read datas
	npart = read_npart(filename)
	pos0 = read_all(filename,'POS ', 0, 3, np.dtype('float32'))
	temp = caltemp(filename)
	mass_0 = mass_type(filename,0)
	

	h = rv / nn

	temp_hot = np.zeros(nn)
	temp_cold = np.zeros(nn)
	temp_warm = np.zeros(nn)

	mass_hot = np.zeros(nn)
	mass_cold = np.zeros(nn)
	mass_warm = np.zeros(nn)

	rr_hot = np.zeros(nn)
	rr_cold = np.zeros(nn)
	rr_warm = np.zeros(nn)

	i = 0
	while npart[0] > i:
		rc = math.sqrt( (pos0[i][0] - xc)**2 + (pos0[i][1] - yc)**2 + (pos0[i][2] - zc)**2) 
		tempk = temp[i] / const.k_B.value * const.e.value * 1e3
		if rc < rv:
			n = int(rc / h)
			if tempk <= 1e5:
				mass_cold[n] += mass_0[i]
				temp_cold[n] += temp[i]*mass_0[i]
				rr_cold[n] += rc / rv * mass_0[i]
			elif 1e5 < tempk < 1e7:
				mass_warm[n] += mass_0[i]
				temp_warm[n] += temp[i]*mass_0[i]
				rr_warm[n] += rc / rv * mass_0[i] 
			else:
				mass_hot[n] += mass_0[i]
				temp_hot[n] += temp[i]*mass_0[i]
				rr_hot[n] = rc / rv * mass_0[i] 
		i += 1

	tempcold = temp_cold / mass_cold
	tempwarm = temp_warm / mass_warm
	temphot = temp_hot / mass_hot
	rrcold = rr_cold / mass_cold
	rrwarm = rr_warm / mass_warm
	rrhot = rr_hot / mass_hot

	return tempcold,tempwarm,temphot,rrcold,rrwarm,rrhot

	
#-------------------------------------------------------------------------------------

def main(model,nn):
	if model == 'Gadget-X':
		speciname = '/home/nifty2014/TheThreeHundred/catalogues/AHF/complete_sample/G3X_Mass_snap_128-center-cluster.txt'
		lab = 'GX'
	if model == 'Gadget-MUSIC':
		speciname = '/home/nifty2014/TheThreeHundred/catalogues/AHF/complete_sample/Music_Mass_snap_017-center-cluster.txt'
		lab = 'GM'
	data = np.loadtxt(speciname)
	size = data.shape[0]

	cold = np.zeros(shape = (324,nn))
	hot = np.zeros(shape = (324,nn))
	warm = np.zeros(shape = (324,nn))
	rcold = np.zeros(shape = (324,nn))
	rwarm = np.zeros(shape = (324,nn))
	rhot = np.zeros(shape = (324,nn))

	ii = 0
	while  ii < size:
		region = int(data[ii][0])

		xc = data[ii][3]
		yc = data[ii][4]
		zc = data[ii][5]
		rv = data[ii][10]

		if model == 'Gadget-X':
			if region <= 9:
				filename = '/home/nifty2014/TheThreeHundred/simulation/GadgetX/NewMDCLUSTER_000%s/snap_128' %region
			elif  9 < region <= 99:
				filename = '/home/nifty2014/TheThreeHundred/simulation/GadgetX/NewMDCLUSTER_00%s/snap_128' %region
			else:
				filename = '/home/nifty2014/TheThreeHundred/simulation/GadgetX/NewMDCLUSTER_0%s/snap_128' %region
		else:
			if region <= 9:
				filename = '/home/nifty2014/TheThreeHundred/simulation/GadgetMUSIC/The300_MUSIC/NewMDCLUSTER_000%s/snap_017' %region
			elif  9 < region <= 99:
				filename = '/home/nifty2014/TheThreeHundred/simulation/GadgetMUSIC/The300_MUSIC/NewMDCLUSTER_00%s/snap_017' %region
			else:
				filename = '/home/nifty2014/TheThreeHundred/simulation/GadgetMUSIC/The300_MUSIC/NewMDCLUSTER_0%s/snap_017' %region

	
		tempcold,tempwarm,temphot,rrcold,rrwarm,rrhot = distribution(xc,yc,zc,rv,filename,nn)

		cold[ii] = tempcold
		warm[ii] = tempwarm
		hot[ii] = temphot
		rcold[ii] = rrcold
		rwarm[ii] = rrwarm
		rhot[ii] = rrhot
			
		ii += 1
		print(ii)

	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/%s500_snap_tempcold.txt' %lab, cold)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/%s500_snap_temphot.txt' %lab, hot)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/%s500_snap_tempwarm.txt' %lab, warm)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/%s500_snap_rrcold.txt' %lab, rcold)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/%s500_snap_rrhot.txt' %lab, rhot)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/%s500_snap_rrwarm.txt' %lab, rwarm)

		
main(model,20)







