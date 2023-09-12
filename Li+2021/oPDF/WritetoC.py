import numpy as np 

def main(model,stack):
	if model == 'GX':
		fileCM = '/home/qyli/oPDFnew/data/cluster/G3X-CM-masscomplete-ALL.txt'
	else:
		fileCM = '/home/qyli/oPDFnew/data/cluster/MDPL2-CM-masscomplete-ALL.txt'
	dataCM = np.loadtxt(fileCM)
	rs = dataCM[:,3] / dataCM[:,6]

	with open('/home/qyli/oPDFnew/oPDF/C/TemplateData.h','w') as obj:
		obj.write('#define LEN_PROF 101\n')
		obj.write('double PotentialTemplate[][3][LEN_PROF]={\n')

	if stack == True:
		rbinall = 0
		potall = 0
		denall = 0

	for i in range(324):
		file = '/home/qyli/oPDFnew/data/1.5template/%s/rin0_rout1_rpd%s.txt' %(model, i+1)
		data = np.loadtxt(file)

		with open('/home/qyli/oPDFnew/oPDF/C/TemplateData.h','a') as obj:
			rbin0 = data[0,:]
			pot0 = data[1,:]
			dens0 = data[2,:]

			if stack == False:
				rbin = str(list(rbin0)).replace('[','').replace(']','')
				pot = str(list(pot0)).replace('[','').replace(']','')
				dens_cul = str(list(dens0)).replace('[','').replace(']','')
				if i != 323:
					obj.write('{ {'+ rbin + '},\n' + '{' + pot  + '},\n' + '{' + dens_cul + '} }, \n ')
				else:
					obj.write('{ {'+ rbin + '},\n' + '{' + pot  + '},\n' + '{' + dens_cul + '} } \n ')
			else:
				rbinall+=rbin0
				potall+=pot0
				denall+=dens0


	if stack == True:
		rs = np.mean(rs)
		rbinall = rbinall / 324
		potall = potall / 324
		denall = denall / 324
		with open('/home/qyli/oPDFnew/oPDF/C/TemplateData.h','a') as obj:

			rbin = str(list(rbin0)).replace('[','').replace(']','')
			pot = str(list(pot0)).replace('[','').replace(']','')
			dens_cul = str(list(dens0)).replace('[','').replace(']','')

			obj.write('{ {'+ rbin + '},\n' + '{' + pot  + '},\n' + '{' + dens_cul + '} } \n ')
			obj.write('};\n')
			# rs0 = str(list(rs)).replace('[','').replace(']','')
			obj.write('double TemplateScale[]={'+ str(rs) + '\n};')
	else:
		with open('/home/qyli/oPDFnew/oPDF/C/TemplateData.h','a') as obj:
			obj.write('};\n')
			rs0 = str(list(rs)).replace('[','').replace(']','')
			obj.write('double TemplateScale[]={'+ rs0 + '\n};')
#'GX' or 'MDPL2'
#fix stack is False
main(model = 'GX',stack = False)
