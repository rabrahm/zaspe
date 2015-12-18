import new2
import numpy as np
import pyfits
import time
import os

f = open('zaspe.pars','r')
lines = f.readlines()

for line in lines:
	cos = line.split()
	if len(cos)==2:
		if cos[0] == 'mod':
			mod = cos[1]
		elif cos[0] == 'library':
			library = cos[1]
		elif cos[0] == 'spec':
			spec = cos[1]
		elif cos[0] == 'RV0':
			RV0 = float(cos[1])
		elif cos[0] == 'vsini':
			guess_vsini = float(cos[1])
		elif cos[0] == 'RESI':
			RESI = float(cos[1])
		elif cos[0] == 'ncores':
			ncores = int(cos[1])
		elif cos[0] == 'trunc':
			trunc = int(cos[1])
		elif cos[0] == 'nit':
			nit = int(cos[1])
		elif cos[0] == 'nsim':
			nsim = int(cos[1])
		elif cos[0] == 'fixG':
			fixG = float(cos[1])
		elif cos[0] == 'efixG':
			efixG = float(cos[1])
		elif cos[0] == 'T':
			T = float(cos[1])
		elif cos[0] == 'G':
			G = float(cos[1])
		elif cos[0] == 'Z':
			Z = float(cos[1])
		elif cos[0] == 'R':
			R = float(cos[1])
		elif cos[0] == 'V':
			V = float(cos[1])
isfits = True
try:
	sc = pyfits.getdata(spec)
	ttemp = spec
except:
	isfits = False
	d = np.loadtxt(spec)
	ords = np.unique(d[:,0]).astype('int')
	first = True
	for o in ords:
		I = np.where(d[:,0]==o)[0]
		if first:
			wave = d[:,1][I]
			flux = d[:,2][I]
			first = False
		else:
			wave = np.vstack((wave,d[:,1][I]))
			flux = np.vstack((flux,d[:,2][I]))
	sc = np.zeros((4,len(ords),wave.shape[1]))
	sc[0] = wave
	sc[3] = flux
	hdu = pyfits.PrimaryHDU(sc)
	rtemp = 'temp_zaspe_spectra.fits'
	if os.access(rtemp,os.F_OK):
		os.system('rm '+rtemp)
	hdu.writeto(rtemp)
	ttemp = spec
	spec = rtemp
	
if library == 'R':
	mod_n = 'Brahm et al. (2016)'
elif library == 'P':
	mod_n = 'Husser et al. (2013)'
elif library == 'C':
	mod_n = 'Coelho et al. (2005)'

print '\n\t ZASPE'

print '\tinput spectrum:', ttemp
print '\twill use the', mod_n, 'library of synthetic spectra.\n'

if 'P' in mod:
	print '\tPerforming the search of the optimal parameters ...'
	pars = new2.get_rough_pars(spec,RV0=RV0,guess_vsini=guess_vsini,RESI=RESI,ncores=ncores,\
                              trunc=trunc,printing=True,use_masks=True,fixG=fixG,nit=nit)
	print '\tZAPE parameters:'
	print '\tTeff=', pars[0], 'log(g)=', pars[1], '[Fe/H]=',pars[2], 'vsin(i)=', pars[3], 'RV=', pars[4]
	
else:
	if T!=-1 and G!=-1 and Z!=-1 and R!=-1 and V!=-1:
		pars = np.array([T,G,Z,R,V])
	else:
		raise ValueError('You have not requested the computation of the stellar parameters and you have not provided them')

if 'E' in mod:
	print '\tPerforming the Monte Carlo simulation to obtain the covariance in the parameters ...'
	mc  = new2.get_precise_parameters(spec,pars,RESI=RESI,ncores=ncores,trunc=trunc,fixG=fixG,nsim=nsim,efixG=efixG)
	print '\tSimulation done.'

	et = np.around(np.sqrt(np.var(mc[:,0])))
	eg = np.around(np.sqrt(np.var(mc[:,1])),3)
	ez = np.around(np.sqrt(np.var(mc[:,2])),3)
	er = np.around(np.sqrt(np.var(mc[:,3])),3)
	ev = np.around(np.sqrt(np.var(mc[:,4])),3)
	
	print '\tZAPE parameters:'
	print '\tTeff=',pars[0], '+/-', et
	print '\tlog(g)=',pars[1], '+/-', eg
	print '\t[Fe/H]=',pars[2], '+/-', ez
	print '\tvsin(i)=',pars[3], '+/-', er
	print '\tRV=', pars[4], '+/-', ev

	out = np.vstack((pars,mc))
	hdu = pyfits.PrimaryHDU(out)
	date = str(time.localtime()).split('=')
	odat = '_'+date[1].split(',')[0]+'_'+date[2].split(',')[0]+'_'+date[3].split(',')[0]+\
			'_'+date[4].split(',')[0]+'_'+date[4].split(',')[0]+'_'+date[6].split(',')[0] + '_' 
	rout = ttemp.split('/')[-1][:-5] + odat + 'zaspe_out.fits'
	hdu.header.update('DATE',odat)
	hdu.header.update('SPEC',ttemp)
	hdu.header.update('RES',RESI)
	hdu.header.update('TRUNC',trunc)
	hdu.header.update('FIXG',False)
	if fixG != -1:
	    hdu.header['FIXG'] = True
	hdu.header.update('TEFF',pars[0])
	hdu.header.update('LOGG',pars[1])
	hdu.header.update('FEH',pars[2])
	hdu.header.update('VSINI',pars[3])
	hdu.header.update('RV',pars[4])
	hdu.header.update('ETEFF',et)
	hdu.header.update('ELOGG',eg)
	hdu.header.update('EFEH',ez)
	hdu.header.update('EVSINI',er)
	hdu.header.update('ERV',ev)
	hdu.header.update('LIB',library)
	hdu.header.update('WI',wi)
	hdu.header.update('WF',wf)
	hdu.writeto(rout)

if not isfits:
	os.system('rm ' + spec)
