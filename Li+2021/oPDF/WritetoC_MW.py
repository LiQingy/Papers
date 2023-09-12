import numpy as np 

def main(model,hub):
	if model == 'S':
		rvfof = np.array([224.1709,172.7433,236.8058,206.0765,192.1194,173.7896,233.5474,194.9914,199.178,189.7775,197.9979,170.9422])
		Ctruefof = np.array([5.8976, 8.8382, 8.8435, 6.9903, 12.776, 8.7318, 10.0735, 7.18, 9.8334, 8.2011, 12.346, 11.1352])
	else:
		rvfof = np.array([242.2135,206.9113,187.9074,186.0567,235.1194,218.5439,221.8365,221.4039,197.6648,193.7516,265.3684,214.9532])
		Ctruefof = np.array([12.4923,4.1164,11.2907,9.8764,8.7675,4.2887,13.519,3.701,8.1664,11.4837,8.0264,8.6198])
	rs = rvfof * hub / Ctruefof

	with open('/home/qyli/oPDFnew/oPDF/C/TemplateData.h','w') as obj:
		obj.write('#define LEN_PROF 101\n')
		obj.write('double PotentialTemplate[][3][LEN_PROF]={\n')

	for i in range(12):
		file = '/home/qyli/oPDFnew/data/template/MW_%s/rin0_rout1_rpd%s%s.txt' %(model,int(i/2)+1,i%2+1)
		data = np.loadtxt(file)

		with open('/home/qyli/oPDFnew/oPDF/C/TemplateData.h','a') as obj:
			rbin = str(list(data[0,:])).replace('[','').replace(']','')
			pot = str(list(data[1,:])).replace('[','').replace(']','')
			dens_cul = str(list(data[2,:])).replace('[','').replace(']','')
			if i != 11:
				obj.write('{ {'+ rbin + '},\n' + '{' + pot  + '},\n' + '{' + dens_cul + '} }, \n ')
			else:
				obj.write('{ {'+ rbin + '},\n' + '{' + pot  + '},\n' + '{' + dens_cul + '} } \n ')
	with open('/home/qyli/oPDFnew/oPDF/C/TemplateData.h','a') as obj:
		obj.write('};\n')
		rs0 = str(list(rs)).replace('[','').replace(']','')
		obj.write('double TemplateScale[]={'+ rs0 + '\n};')

#'GX' or 'MDPL2'
main(model = 'S',hub = 1.)
