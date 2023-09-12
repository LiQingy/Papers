import numpy as np 
import math
from struct import unpack
from os import fstat
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

def AHF(n,xc,yc,zc,rv,region,model,file):
	pos0 = read_all(file,'POS ', 0, 3, np.dtype('float32'))
	pos1 = read_all(file,'POS ', 1, 3, np.dtype('float32'))
	pos4 = read_all(file,'POS ', 4, 3, np.dtype('float32'))
	if model == 'Gadget-X':
		pos5 = read_all(file,'POS ', 5, 3, np.dtype('float32'))
		mass5 = mass_type(file,5)
	mass0 = mass_type(file,0)
	mass1 = mass_type(file,1)
	mass4 = mass_type(file,4)
	

	size0 = pos0.shape[0]
	size1 = pos1.shape[0]
	size4 = pos4.shape[0]
	h = rv / n

	massall = np.zeros(n)
	rrall0 = np.zeros(n)
	
	for j in range(size0):
		rx = pos0[j][0] / 0.678
		ry = pos0[j][1] / 0.678
		rz = pos0[j][2] / 0.678
		rr = math.sqrt((rx - xc)**2 + (ry - yc)**2 + (rz - zc)**2)
		if rr < rv:
			nn = int(rr / h)
			massall[nn] += mass0[j] / 0.678 #all mass
			rrall0[nn] += rr / rv * mass0[j] / 0.678

	for j in range(size1):
		rx = pos1[j][0] / 0.678
		ry = pos1[j][1] / 0.678
		rz = pos1[j][2] / 0.678
		rr = math.sqrt((rx - xc)**2 + (ry - yc)**2 + (rz - zc)**2)
		if rr < rv:
			nn = int(rr / h)
			massall[nn] += mass1[j] / 0.678 #all mass
			rrall0[nn] += rr / rv * mass1[j] / 0.678

	for j in range(size4):
		rx = pos4[j][0] / 0.678
		ry = pos4[j][1] / 0.678
		rz = pos4[j][2] / 0.678
		rr = math.sqrt((rx - xc)**2 + (ry - yc)**2 + (rz - zc)**2)
		if rr < rv:
			nn = int(rr / h)
			massall[nn] += mass4[j] / 0.678 #all mass
			rrall0[nn] += rr / rv * mass4[j] / 0.678

	if model == 'Gadget-X':
		for j in range(size5):
			rx = pos5[j][0] / 0.678
			ry = pos5[j][1] / 0.678
			rz = pos5[j][2] / 0.678
			rr = math.sqrt((rx - xc)**2 + (ry - yc)**2 + (rz - zc)**2)
			if rr < rv:
				nn = int(rr / h)
				massall[nn] += mass5[j] / 0.678 #all mass
				rrall0[nn] += rr / rv * mass5[j] / 0.678						

#-----------------------------------------------------------------------------------------------------
	ndv = np.zeros(n)
	for i in range(n):
		dv = 4. / 3 * math.pi * (pow((i + 1) * h,3) - pow(i * h,3))
		ndv[i] = dv

	densityall0 = massall / ndv
	rrall0 = rrall0 / massall

	return densityall0,rrall0

def main(n,model,savetxt):
	if model == 'Gadget-X':
		filecen = '/home/nifty2014/TheThreeHundred/playground/weiguang/G3X_Mass_snap_128-center-cluster.txt'
		name = 'GX'
	if model == 'Gadget-MUSIC':
		filecen = '/home/nifty2014/TheThreeHundred/playground/weiguang/Music_Mass_snap_017-center-cluster.txt'
		name = 'GM'
	datacen = np.loadtxt(filecen)

	densityall = np.zeros(shape = (324,n))
	rrall = np.zeros(shape = (324,n))

	region = 1
	while  region <= 324:

		rv = datacen[region - 1][6] / 0.678
		xc = datacen[region - 1][3] / 0.678
		yc = datacen[region - 1][4] / 0.678
		zc = datacen[region - 1][5] / 0.678

		if model == 'GadgetX':
			filename = '/home/nifty2014/TheThreeHundred/simulation/GadgetX/NewMDCLUSTER_0%03d/snap_128' %region
		if model == 'Gadget-MUSIC':
			filename = '/home/nifty2014/TheThreeHundred/simulation/GadgetMUSIC/The300_MUSIC/NewMDCLUSTER_0%03d/snap_017' %region


		densityall0,rrall0 = AHF(n,xc,yc,zc,rv,region,model,filename)

		densityall[region - 1] = densityall0
		rrall[region - 1] = rrall0

		print(region)
		region += 1

	if savetxt == True:
		np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/AHFdata/%s200_AHF_densityall.txt' %name,densityall)
		np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/AHFdata/%s200_AHF_rrall.txt' %name,rrall)

	print(model)

main(20, 'Gadget-MUSIC', savetxt = True)
main(20, 'Gadget-X', savetxt = True)