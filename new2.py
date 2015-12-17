import sys
#sys.path.append('/data/echelle/AUTOMOOGPUC/degradations/')
#sys.path.append('/data/echelle/AUTOMOOGPUC/rot_conv/')
import integration
import rot_conv
import pyfits
#import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import ndimage
from scipy import interpolate
from scipy import optimize
from scipy import signal
from multiprocessing import Pool
import time
import pickle
import os
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages
#from pyevolve import *

def n_Edlen(l):
    sigma = 1e4 / l
    sigma2 = sigma*sigma
    n = 1 + 1e-8 * (8342.13 + 2406030 / (130-sigma2) + 15997/(38.9-sigma2))
    return n

def n_Morton(l):
    sigma = 1e4 / l
    sigma2 = sigma*sigma
    n = 1 + 6.4328e-5 + 2.94981e-2 / (146.-sigma2) + 2.5540e-4/(41.-sigma2)
    return n

def ToAir(l):
    return (l / n_Edlen(l))

def ToVacuum(l):
    cond = 1
    l_prev = l.copy()
    while(cond):
        l_new = n_Edlen(l_prev) * l
        if (max(abs(l_new - l_prev)) < 1e-10): cond = 0
        l_prev = l_new
    
    return l_prev

f = open('zaspe.pars','r')
lines = f.readlines()
for line in lines:
	cos = line.split()
	if len(cos)==2:
		if cos[0] == 'library':
			library = cos[1]
		elif cos[0] == 'wi':
			wi = float(cos[1])
		elif cos[0] == 'wf':
			wf = float(cos[1])
		elif cos[0] == 'wavP':
			wavP = cos[1]
		elif cos[0] == 'modP':
			modP = cos[1]
		elif cos[0] == 'modC':
			modC = cos[1]
		elif cos[0] == 'modR':
			modR = cos[1]
#library = 'C'
#library = 'P'
#library = 'R'
#wi,wf = 5300,6000

linear_interpolation = False
sim_type = 'factors'

lux = 299792.458 # speed of light in km/s

if library == 'P':
	wav_dir = wavP
	mod_dir = modP
	mw = pyfits.getdata(wav_dir)
	ti,tf,dt = 4000,7000,100
	gi,gf,dg = 1.0,5.0,0.5
	zi,zf,dz = -1.0,1.0,0.5
	ts = np.arange(ti,tf+1,dt)
	gs = np.arange(gi,gf+.1,dg)
	zs = np.arange(zi,zf+.1,dz)
	rs = np.arange(0.5,15.5,5.)
elif library == 'C':
	mod_dir = modC
	mhd = pyfits.getheader(mod_dir+'7000_50_p05p00.ms.fits')
	mw = ToVacuum(np.arange(mhd['NAXIS1'])*mhd['CD1_1'] + mhd['CRVAL1'])
	ti,tf,dt = 4000,7000,250
	gi,gf,dg = 1.0,5.0,0.5
	zi,zf,dz = -2.0,0.5,0.5
	ts = np.arange(ti,tf+1,dt)
	gs = np.arange(gi,gf+.1,dg)
	zs = np.arange(zi,zf+0.1,dz)
elif library == 'R':
	mod_dir = modR
	mhd = pyfits.getheader(mod_dir+'4000_50_p000.fits')
	mw = ToVacuum(np.arange(mhd['NAXIS1'])*float(mhd['CD1_1']) + float(mhd['CRVAL1']))
	ti,tf,dt = 4000,7000,200
	gi,gf,dg = 1.0,5.0,0.5
	zi,zf,dz = -1.0,0.5,0.25
	ts = np.arange(ti,tf+1,dt)
	gs = np.arange(gi,gf+.1,dg)
	zs = np.arange(-1.0,0.6,0.25)

Iwav = np.where((mw>wi-100)&(mw<wf+100))[0]
mw = mw[Iwav]

def get_vmac(T):
	"""
	this function returns the macroturbulence velocity in km/s given the
	effective temperature of the star following valenti & fisher 2005
	"""
	m = 0.001529
	n = -4.788
	vmac = m*T + n
	if vmac < 1.5:
		vmac = 1.5
	return vmac

def get_oriname(t,g,z):
	st = str(int(t))
	#print t,st
	sg = str(g)
	sz = str(np.absolute(z))
	zsign = '+'
	zsignc = 'p'
	if z < 0.:
		zsign = '-'
		zsignc = 'm'
	if z==0:
		zsign = '-'
	if library == 'P':
		return mod_dir + 'Z' + zsign + sz + '/lte0' + st + '-' + sg + '0' + zsign + sz + '.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits' 
	elif library == 'C':
		return mod_dir + st + '_' + sg[0] + sg[2] + '_' + zsignc + sz[0] + sz[2] + 'p00.ms.fits'
	elif library == 'R':
		sg = str(int(10*g))
		sign = 'p'
		if z < 0:
			sign = 'm'
		sz = str(int(np.absolute(z*100)))
		if np.absolute(z)<1:
			sz = '0'+sz
		if len(sz)<3:
			sz = '0'+sz
		name = str(int(t))+'_'+sg + '_' + sign + sz + '.fits'
		return mod_dir + name

def get_interpolated(t,g,z):
	dts = np.argmin((ts - t)**2)

	if dts > 1 and dts < len(ts) - 2:
		dts = ts[dts-2:dts+3]
	elif dts == 1:
		dts = ts[dts-1:dts+4]
	elif dts == 0:
		dts = ts[dts-0:dts+5]
	elif dts == len(ts) - 2:
		dts = ts[dts-3:dts+2]
	elif dts == len(ts) - 1:
		dts = ts[dts-4:dts+1]
	
	dgs = np.argmin((gs - g)**2)
	if dgs > 1 and dgs < len(gs) - 2:
		dgs = gs[dgs-2:dgs+3]
	elif dgs == 1:
		dgs = gs[dgs-1:dgs+4]
	elif dgs == 0:
		dgs = gs[dgs-0:dgs+5]
	elif dgs == len(gs) - 2:
		dgs = gs[dgs-3:dgs+2]
	elif dgs == len(gs) - 1:
		dgs = gs[dgs-4:dgs+1]

	dzs = np.argmin((zs - z)**2)
	if dzs > 1 and dzs < len(zs) - 2:
		dzs = zs[dzs-2:dzs+3]
	elif dzs == 1:
		dzs = zs[dzs-1:dzs+4]
	elif dzs == 0:
		dzs = zs[dzs-0:dzs+5]
	elif dzs == len(zs) - 2:
		dzs = zs[dzs-3:dzs+2]
	elif dzs == len(zs) - 1:
		dzs = zs[dzs-4:dzs+1]
	data = np.zeros((len(dts),len(dgs),len(dzs),len(mw)))
	data2 = np.zeros((len(dts),len(dgs),len(dzs)))

	for it in np.arange(len(dts)):
		for ig in np.arange(len(dgs)):
			for iz in np.arange(len(dzs)):
				if library == 'P':
					data[it,ig,iz] = pyfits.getdata(get_oriname(dts[it],dgs[ig],dzs[iz]))[Iwav]
					mmhd = pyfits.getheader(get_oriname(dts[it],dgs[ig],dzs[iz]))
					if 'PHXXI_L' in mmhd.keys():
						data2[it,ig,iz] = mmhd['PHXXI_L']
					else:
						data2[it,ig,iz] = 2.0
				elif library == 'C':
					data[it,ig,iz] = pyfits.getdata(get_oriname(dts[it],dgs[ig],dzs[iz]))[0][Iwav]
				elif library == 'R':
					data[it,ig,iz] = pyfits.getdata(get_oriname(dts[it],dgs[ig],dzs[iz]))[Iwav]

	tckt = interpolate.splrep(dts,np.arange(len(dts)),k=1)
	tckg = interpolate.splrep(dgs,np.arange(len(dgs)),k=1)
	tckz = interpolate.splrep(dzs,np.arange(len(dzs)),k=1)
	A = np.zeros(len(mw)) + interpolate.splev(t,tckt)
	B = np.zeros(len(mw)) + interpolate.splev(g,tckg)
	C = np.zeros(len(mw)) + interpolate.splev(z,tckz)
	D = np.arange(len(mw))

	coords = np.vstack((A,B,C,D))
	zi = ndimage.map_coordinates(data, coords, order=3, mode='nearest')
	if library == 'P':
		coords2 = np.array([[interpolate.splev(t,tckt),interpolate.splev(g,tckg),interpolate.splev(z,tckz)],]).T
		zi2 = ndimage.map_coordinates(data2, coords2, order=3, mode='nearest')
		mict = zi2[0]
	elif library == 'C' or library == 'R':
		if g>=2.75:
			mict = 1.
		elif g >=1.25 and g<=2.75:
			mict = 1.8
		elif g < 1.25:
			mict = 2.5
	return np.array(zi),mict

def get_model(t,g,z):
	if os.access(get_oriname(t,g,z),os.F_OK):
		if library == 'P':
			mf = pyfits.getdata(get_oriname(t,g,z))[Iwav]
			mhd = pyfits.getheader(get_oriname(t,g,z))
			if 'PHXXI_L' in mhd.keys():
				mict = float(mhd['PHXXI_L'])
			else:
				mict = 2.0
		elif library == 'C':
			mf = pyfits.getdata(get_oriname(t,g,z))[0][Iwav]
			if g>=2.75:
				mict = 1.
			elif g >=1.25 and g<=2.75:
				mict = 1.8
			elif g < 1.25:
				mict = 2.5
		elif library == 'R':
			mf = pyfits.getdata(get_oriname(t,g,z))[Iwav]
			if g>=2.75:
				mict = 1.
			elif g >=1.25 and g<=2.75:
				mict = 1.8
			elif g < 1.25:
				mict = 2.5
	
	else:
		if linear_interpolation:
			mf,mict = get_linear_interpol(t,g,z)
		else:
			mf,mict = get_interpolated(t,g,z)
	return mf,mict


def pixelization(x,y,rwav,nref=50):
	tck = interpolate.splrep(np.arange(len(rwav))+0.5,rwav,k=3)
	X = np.arange(nref*len(rwav))/float(nref)
	Y = interpolate.splev(X,tck)
	tck2 = interpolate.splrep(x,y,k=3)
	ny = interpolate.splev(Y,tck2)
	ny = ny.reshape((len(rwav),nref))
	ny = np.sum(ny,axis=1)/float(nref)
	ny[0]=ny[1]
	ny[-1]=ny[-2]

	return ny

def get_mict(t,g,z):
	mmhd = pyfits.getheader(get_oriname(t,g,z))
	if 'PHXXI_L' in mmhd.keys():
		return mmhd['PHXXI_L']
	else:
		return 2.0

def get_full_model(t,g,z,rot,R):
	mf,mict=get_model(t,g,z)
	vmac = get_vmac(t)
	sigma_mac = 0.297 * vmac
	sigma_mac = np.zeros(len(mw)) + lux / sigma_mac
	nmf = mf.copy()
	gn = integration.InstConvVarGau(mw.astype('double'),mf.astype('double'),nmf.astype('double'),sigma_mac.astype('double'),len(mw))
	mf = np.array(gn)
	ac,bc = rot_conv.get_ldcoef2(t, g, z, mict)
	#ac = np.array([ 0.5268,  0.4355,  0.2671,  0.1921,  0.1423])
 	#bc = np.array([ 0.2938,  0.3151,  0.3717,  0.364 ,  0.3541])
	rflx = rot_conv.conv2(mw,mf,ac,bc,rot)
	mf = np.array(rflx)
	nmf  = mf.copy()
	R = np.zeros(len(mw))+R
	rflx = integration.InstConvVarGau(mw.astype('double'),mf.astype('double'),nmf.astype('double'),R.astype('double'),len(mw))
	mf  = np.array(rflx)
	return mf

def get_cont(x,y,n=1,sl=1.,sh=5.):
	orilen = len(x)
	coef = np.polyfit(x,y,n)
	res = y - np.polyval(coef,x)
	IH = np.where(res>0)[0]
	IL = np.where(res<0)[0]
	dev = np.mean(res[IH])
	I = np.where((res>-sl*dev) & (res<sh*dev))[0]
	J1 = np.where(res<=-sl*dev)[0]
	J2 = np.where(res>=sh*dev)[0]
	J = np.unique(np.hstack((J1,J2)))
	cond = True
	if len(J)==0 or len(x)< .3*orilen:
		cond=False
	while cond:
		x = np.delete(x,J)
		y = np.delete(y,J)
		coef = np.polyfit(x,y,n)
		res = y - np.polyval(coef,x)
		IH = np.where(res>0)[0]
		IL = np.where(res<0)[0]
		dev = np.mean(res[IH])
		I = np.where((res>-sl*dev) & (res<sh*dev))[0]
		J1 = np.where(res<=-sl*dev)[0]
		J2 = np.where(res>=sh*dev)[0]
		J = np.unique(np.hstack((J1,J2)))
		cond = True
		if len(J)==0 or len(x)< .1*orilen:
			cond=False
	return coef

def get_ratio(sciw,rat,n=3):
	rat = scipy.signal.medfilt(rat,11)
	lori = len(sciw)
	coef = np.polyfit(sciw,rat,n)
	res = rat - np.polyval(coef,sciw)
	rms = np.sqrt(np.mean(res**2))
	I = np.where(res> 3*rms)[0]
	I2 = np.where(res< -3*rms)[0]
	I = np.sort(np.hstack((I,I2)))
	cond = True
	if len(I) == 0 or len(sciw) < .3 * lori:
		cond = False

	while cond:
		#imax = np.argmax(res**2)
		#sciw = np.delete(sciw,imax)
		#rat  = np.delete(rat,imax)
		sciw = np.delete(sciw,I)
		rat  = np.delete(rat,I)
		coef = np.polyfit(sciw,rat,n)
		res = rat - np.polyval(coef,sciw)
		rms = np.sqrt(np.mean(res**2))
		I = np.where(res> 3*rms)[0]
		I2 = np.where(res< -3*rms)[0]
		I = np.sort(np.hstack((I,I2)))
		if len(I) == 0 or len(sciw) < .3 * lori:
			cond = False

	return coef

def get_chis_comp(pars, save=False,plt=False):
	t = pars[0]
	g = pars[1]
	z = pars[2]
	rot = pars[3]
	mf = get_full_model(t,g,z,rot,RES_POW)
	mwav,mflx,wav,flx = mw.copy(),mf.copy(),sc[0].copy(),sc[3].copy()
	res = np.array([])
	Is = 0
	largs = 0
	tmpw,tmpf,tmpm,tmpo = np.array([]),np.array([]),np.array([]),np.array([])
	for i in ords:
		I = np.where((mwav>wav[i,0]) & (mwav<wav[i,-1]))[0]
		modw = mwav[I]
		modf = mflx[I]
		sciw = wav[i]
		scif = flx[i]/np.median(flx[i])
		modf = pixelization(modw,modf,sciw)

		if len(zmask)>0:
			II = np.where(zmask[i]!=0)[0]
			modf /= np.mean(modf[II])
		else:
			modf /= modf.mean()
		mscif = scipy.signal.medfilt(scif,11)
		INF = np.where(mscif!=0)[0]
                I0 = np.where(mscif==0)[0]
                mscif[I0] = 1.
                rat = modf/mscif
                rat[I0] = 0.
		#IB = np.where(bmask[i]!=0)[0]
		coef = get_ratio(sciw[INF],rat[INF])
		scif = scif * np.polyval(coef,sciw)
		mscif = mscif * np.polyval(coef,sciw)
		coef = get_cont(sciw,mscif)
		scif = scif / np.polyval(coef,sciw)
		coef = get_cont(sciw,modf)
		modf = modf / np.polyval(coef,sciw)
		
		if plt:
			#print len(tmpo)
			plot(sciw,scif,'r')
			plot(sciw,modf,'b')
			tmpw = np.hstack((tmpw,sciw))
			tmpf = np.hstack((tmpf,scif))
			tmpm = np.hstack((tmpm,modf))
			tmpo = np.hstack((tmpo,np.zeros(len(modf))+i))
			
		#	show()
		#	print fds
		#	#plot(sciw,modf)
		rest = scif - modf
		largI = len(scif)
		if len(zmask)>0:
			if sim_type == 'factors':
				#plot(sciw,modf,'b')
				modf = -(modf - 1.)*zmask[i]
				modf = -modf + 1.
				#plot(sciw,modf,'r')
				rest = scif - modf
			else:
				rest = scif - modf*zmask[i]
			II = np.where(zmask[i]==0)[0]
			rest[II] = 0.
			largI -= len(II) 
		largs += largI

		if len(res) == 0:
			res = rest.copy()
		else:
			res = np.vstack((res,rest))
	if plt:
		ftemp=open('temp_sp_ff.txt','w')
		outtt = np.vstack((tmpw,tmpf,tmpm,tmpo)).T
		np.savetxt('temp_sp_ff.txt',outtt)
		show()
		print vgbhnj
	#if plt:
	#	show()
	#	print gfds
	#res, largs = compare(mw,mf,sc[0],sc[3],save=save)
	chiss = np.sum(res**2) / float(largs)
	#print pars,np.sum(res**2)
	return chiss

def gauccf(params,vels):
	f = params[3] * np.exp(-0.5*(vels - params[0])**2/params[1]**2) + params[2]
	return f

def res_gauccf(params,ccf,vels):
	return ccf - gauccf(params,vels)

def get_fwhm(VELS,CCF,dplot=False):
	guess = np.array([0.,1.,CCF.min(),(CCF-CCF.min()).max()])
	gauss_fit = optimize.leastsq(res_gauccf,guess,args=(CCF,VELS))[0]
	im = np.argmax(CCF)
	base = gauss_fit[2]
	CCF -= base
	CCF /= CCF.max()
	#print VELS
	#print CCF
	if dplot:
		plot(VELS,CCF)
		show()
		print fds
	v1,c1,v2,c2 = VELS[:im],CCF[:im],VELS[im:],CCF[im:]
	try:
		I1 = np.where(c1 < 0.5)[0][-1]
	except:
		I1 = 1

	m = (v1[I1] - v1[I1+1]) / (c1[I1] - c1[I1+1])
	n = v1[I1] - m*c1[I1]
	vli = m*0.5 + n
	try:
		I1 = np.where(c2 < 0.5)[0][0]
	except:
		I1 = len(c2)- 2
	m = (v2[I1] - v2[I1-1]) / (c2[I1] - c2[I1-1])
	n = v2[I1] - m*c2[I1]
	vld = m*0.5 + n
	return vld-vli, gauss_fit[0]

def sim_rot(pars):
	rot = pars[0]
	t,g,z = pars[1],pars[2],pars[3]
	wav = pars[4]
	flx = pars[5]
	vels = pars[6]

	mf,mict = get_model(t,g,z)
	vmac = get_vmac(t)
	sigma_mac = 0.297 * vmac
	sigma_mac = np.zeros(len(mw)) + lux / sigma_mac
	nmf = mf.copy()
	gn = integration.InstConvVarGau(mw.astype('double'),mf.astype('double'),nmf.astype('double'),sigma_mac.astype('double'),len(mw))
	mf = np.array(gn)
	nmf  = mf.copy()
	R = np.zeros(len(mw)) + RES_POW
	rflx = integration.InstConvVarGau(mw.astype('double'),mf.astype('double'),nmf.astype('double'),R.astype('double'),len(mw))
	mf  = np.array(rflx)/np.median(np.array(rflx))

	ac,bc = rot_conv.get_ldcoef2(t, g, z, mict)
	rflx = rot_conv.conv2(mw,mf,ac,bc,rot)
	tmf = np.array(rflx)
	ccft = np.array([])
	nwav = wav.copy()
	nflx = flx.copy()
	for i in ords:
		sciw = wav[i]
		I = np.where((mw>sciw[0])&(mw<sciw[-1]))[0]
		I = np.hstack((I[0]-1,I,I[-1]+1))
		modf = pixelization(mw[I],mf[I],sciw)
		scif = pixelization(mw[I],tmf[I],sciw)			
		mscif = scipy.signal.medfilt(scif,11)
		rat = modf/mscif
		INF = np.where(mscif!=0)[0]
		coef = get_ratio(sciw[INF],rat[INF])
		nflx[i] = scif * np.polyval(coef,sciw)

		#III = np.where(mask_bin[i]!=0)[0]
		#print len(III)

	for v in vels:
		twav = mw*(1 + v/lux)
		ccf = np.zeros((2,len(ords),wav.shape[1]))
		for i in ords:
			sciw = nwav[i]
			scif = nflx[i]
			I = np.where((twav>sciw[0])&(twav<sciw[-1]))[0]
			I = np.hstack((I[0]-1,I,I[-1]+1))
			modf = pixelization(twav[I],mf[I],sciw)
			ccf[0,i-ords.min()] = scif / np.sum(scif*rot_mask[i])
			ccf[1,i-ords.min()] = modf / np.sum(modf*rot_mask[i])
		ccft = np.hstack((ccft,np.sum(ccf[0]*ccf[1]*rot_mask[ords[0]:ords[-1]+1])))
	fwhm,rv = get_fwhm(vels,ccft)
	return fwhm

def get_vsini(pars,dplot=False):
	t = pars[0]
	g = pars[1]
	z = pars[2]
	wav = pars[3]
	flx = pars[4]
	rots = pars[5]
	ur = rots.mean()

	mf,mict = get_model(t,g,z)
	vmac = get_vmac(t)
	sigma_mac = 0.297 * vmac
	sigma_mac = np.zeros(len(mw)) + lux / sigma_mac
	nmf = mf.copy()
	gn = integration.InstConvVarGau(mw.astype('double'),mf.astype('double'),nmf.astype('double'),sigma_mac.astype('double'),len(mw))
	mf = np.array(gn)
	nmf  = mf.copy()
	R = np.zeros(len(mw))+RES_POW
	rflx = integration.InstConvVarGau(mw.astype('double'),mf.astype('double'),nmf.astype('double'),R.astype('double'),len(mw))
	mf  = np.array(rflx)/np.median(np.array(rflx))

	vels = np.arange(-7*ur,7*ur,1.0)
	ccft = np.array([])

	coefs = []
	for i in ords:
		sciw = wav[i]
		scif = flx[i]
		I = np.where((mw>sciw[0])&(mw<sciw[-1]))[0]
		I = np.hstack((I[0]-1,I,I[-1]+1))
		#tck = interpolate.splrep(mw[I],mf[I],k=3)
		modf = pixelization(mw[I],mf[I],sciw)
		if sim_type == 'factors':
			modf = -(modf -1) * zmask2[i]
			modf = -modf + 1.
		else:
			modf = modf * zmask2[i]
		mscif = scipy.signal.medfilt(scif,11)
		I0 = np.where(mscif==0)[0]
                mscif[I0] = 1.
                rat = modf/mscif
                rat[I0] = 0.
		INF = np.where(mscif!=0)[0]
		coef = get_ratio(sciw[INF],rat[INF])
		coefs.append(coef)

	for v in vels:
		twav = mw*(1 + v/lux)
		ccf = np.zeros((2,len(ords),wav.shape[1]))
		for i in ords:
			sciw = wav[i]
			scif = flx[i]
			I = np.where((twav>sciw[0])&(twav<sciw[-1]))[0]
			if I[-1]+1 <= len(twav):
				I = np.hstack((I,I[-1]+1))
			if I[0]+-1 >= 0:
				I = np.hstack((I[0]-1,I))
			#I = np.hstack((I[0]-1,I,I[-1]+1))
			modf = pixelization(twav[I],mf[I],sciw)
			scif = scif * np.polyval(coefs[i-ords.min()],sciw)
			if sim_type == 'factors':
				modf = -(modf-1)*zmask2[i]
				modf = -modf + 1.
				scif = -(scif-1)*zmask2[i]
				scif = -scif + 1.
			else:
				modf *= zmask2[i]
				scif = scif * zmask2[i]
			ccf[0,i-ords.min()] = scif / np.sum(scif*rot_mask[i])
			ccf[1,i-ords.min()] = modf / np.sum(modf*rot_mask[i])
		ccft = np.hstack((ccft,np.sum(ccf[0]*ccf[1]*rot_mask[ords[0]:ords[-1]+1])))
	#plot(vels,ccft)
	#show()
	#print gvfcd

	obfwhm,NRV = get_fwhm(vels,ccft,dplot=dplot)
	fwhms = np.array([])
	pars = []
	for rot in rots:
		pars.append([rot,t,g,z,wav,flx,vels])
	#	print rot
	#	sim_rot(pars[-1])

	p = Pool(npools)
	t0 = time.time()
	fwhms = np.array((p.map(sim_rot, pars)))
	p.terminate()

	tck = interpolate.splrep(fwhms,rots,k=3)

	rout = interpolate.splev(obfwhm,tck)

	return rout,NRV

def get_rough_pars(spec,RV0=0,guess_vsini=5.,RESI=120000.,ncores=6,mask=[],trunc=0,errors=False,use_masks=False,very_rough=False,printing=False,fixG=-1,elim=0.1,zmin=500,nit=15):

	global namespec
	namespec = spec
	global sc,hd,RES_POW,npools,ords,zmask, zmask2, mask_bin, rot_mask

	RES_POW = RESI
	npools = ncores
	sc = pyfits.getdata(spec)
	sc = sc[:,:,trunc:sc.shape[2]-trunc]
	hd = pyfits.getheader(spec)
	
	"""
	zmask2 = np.ones(sc[0].shape)
	mask_bin = np.zeros(sc[0].shape)
	change = False
	if len(mask)>0:
		zmask = mask.copy()
		if not errors:
			change = True
		for i in range(zmask.shape[0]):
			II = np.where(zmask[i]==0)[0]
			III = np.where(zmask[i]!=0)[0]
			zmask2[i][II]=1.
			mask_bin[i][III] = 1.
	else:
		mask_bin = np.ones(sc[0].shape)
		zmask = np.ones(sc[0].shape)
	"""
	if not errors:
		mask_bin = np.ones(sc[0].shape)
		zmask    = np.ones(sc[0].shape)
		zmask2   = np.ones(sc[0].shape)
		rot_mask   = np.ones(sc[0].shape)
	else:
		zmask = mask.copy()
		zmask2 = zmask.copy()
		mask_bin = np.zeros(sc[0].shape)
		for i in range(zmask.shape[0]):
			II = np.where(zmask[i]==0)[0]
			III = np.where(zmask[i]!=0)[0]
			zmask2[i][II]=1.
			mask_bin[i][III] = 1.

	#clean data
	
	for i in range(sc.shape[1]):
		rms = np.sqrt(np.var(sc[3,i]))
		res = sc[3,i] - sc[3,i].mean()
		I = np.where(res > 5*rms)[0]
		J = np.where(res <= 5*rms)[0]
		cond = True
		if len(I) == 0:
			cond = False
		while cond:
			sc[3,i,I] = np.mean(sc[3,i,J])
			rms = np.sqrt(np.var(sc[3,i]))
			res = sc[3,i] - sc[3,i].mean()
			I = np.where(res > 5*rms)[0]
			J = np.where(res <= 5*rms)[0]
			if len(I) == 0:
				cond = False


	if 'RV' in hd.keys():
		sc[0] *= (1. - float(hd['RV'])/lux)
		RVTOT = float(hd['RV'])
	else:
		sc[0] *= (1. - RV0/lux)
		RVTOT = RV0

	# determine which echelle orders are going to be used
	ords = []
	for i in range(sc.shape[1]):
		J1 = np.where(mw > sc[0,i,-1]*(1. + 100./lux))[0]
		J2 = np.where(mw < sc[0,i,0]*(1. - 100./lux))[0]
		if len(J1)>0 and len(J2)>0:
			ords.append(i)

	ords = np.array(ords)
	#ords = np.array([6,7])
	#get_chis_comp([4000,4.00,0.3,2.0],plt=True)
	#print gfdx
	#get_chis_comp([5750,3.50,-0.5,10.],plt=True)
	#show()
	#print gfd
	if guess_vsini == -1:
		rots = np.arange(0.5,50.6,10.)
		curr,NRV = get_vsini([5700,4.0,0.0,sc[0],sc[3],rots])
		if curr < 0.5:
			curr = 0.5
	else:
		curr = guess_vsini
	print curr
	all_currs,all_curts,all_curgs,all_curzs = np.array([]),np.array([]),np.array([]),np.array([])
	#curr = guess_vsini
	ik = 0
	while ik < nit:
		if ik < 1  or very_rough:
			
			rts = np.arange(4000,7001,200)
			rgs = np.arange(1.5, 4.6, 0.5)
			rzs = np.arange(-1.0,0.6,0.5)
			rzs = np.array([-1.0,-0.5,0.0,0.5])
		else:
			if ik < 2 or very_rough:
				lowt,upt,ddt = np.around(curt-500,-2),np.around(curt+500,-2) + 1, 100
				lowg,upg,ddg = np.around(curg-.6,1),np.around(curg+.6,1)+.01, 0.2
				lowz,upz,ddz = np.around(curz-.4,1),np.around(curz+.4,1)+.01, 0.1
			elif ik < 3:
				lowt,upt,ddt = np.around(curt-300,-2),np.around(curt+300,-2) + 1, 75
				lowg,upg,ddg = np.around(curg-.6,1),np.around(curg+.6,1)+.01, 0.2
				lowz,upz,ddz = np.around(curz-.3,1),np.around(curz+.3,1)+.01, 0.075
			elif ik < 4:
				lowt,upt,ddt = np.around(curt-200,-2),np.around(curt+200,-2) + 1, 50
				lowg,upg,ddg = np.around(curg-.4,1),np.around(curg+.4,1)+.01, 0.1
				lowz,upz,ddz = np.around(curz-.2,1),np.around(curz+.2,1)+.01, 0.05
			else:
				lowt,upt,ddt = np.around(curt-50,-1),np.around(curt+50,-1) + 1, 10
				lowg,upg,ddg = np.around(curg-.2,1),np.around(curg+.2,1)+.01, 0.05
				lowz,upz,ddz = np.around(curz-.06,1),np.around(curz+.06,1)+.01, 0.02		

			if lowt < ts.min():
				lowt = ts.min()
			if  upt > ts.max():
				upt = ts.max()+1
			if lowg < gs.min():
				lowg = gs.min()
			if upg > gs.max():
				upg = gs.max()+0.01
			if lowz < zs.min():
				lowz = zs.min()
			if  upz > zs.max():
				upz = zs.max()+0.01

			rts = np.arange(lowt,upt,ddt)
			rgs = np.arange(lowg, upg, ddg)
			rzs = np.arange(lowz,upz,ddz)
		if fixG != -1:
			rgs = np.array([fixG])

		rrs = np.array([curr])
		nrts = np.repeat(rts, len(rgs)*len(rzs))
		nrgs = np.repeat(rgs, len(rzs))
		nrzs = np.tile(rzs,len(rgs))
		nrgs = np.tile(nrgs,len(rts))
		nrzs = np.tile(nrzs,len(rts))

		models = np.vstack((nrts,nrgs,nrzs,np.zeros(len(nrts))+curr)).T

		"""
		tckt = interpolate.splrep(rts,np.arange(len(rts)),k=1)
		tckg = interpolate.splrep(rgs,np.arange(len(rgs)),k=1)
		tckz = interpolate.splrep(rzs,np.arange(len(rzs)),k=1)
		tckit = interpolate.splrep(np.arange(len(rts)),rts,k=1)
		tckig = interpolate.splrep(np.arange(len(rgs)),rgs,k=1)
		tckiz = interpolate.splrep(np.arange(len(rzs)),rzs,k=1)

		if ik < 1:
			fts = interpolate.splev(np.arange(rts[0],rts[-1]+1,100.),tckt)
			fgs = interpolate.splev(np.arange(rgs[0],rgs[-1]+.001,.1),tckg)
			fzs = interpolate.splev(np.arange(rzs[0],rzs[-1]+.001,.1),tckz)
		else:

			fts = interpolate.splev(np.arange(rts[0],rts[-1]+1,1.),tckt)
			fgs = interpolate.splev(np.arange(rgs[0],rgs[-1]+.001,.01),tckg)
			fzs = interpolate.splev(np.arange(rzs[0],rzs[-1]+.001,.01),tckz)
		nfts = np.repeat(fts, len(fgs)*len(fzs))
		nfgs = np.repeat(fgs, len(fzs))
		nfzs = np.tile(fzs,len(fgs))
		nfgs = np.tile(nfgs,len(fts))
		nfzs = np.tile(nfzs,len(fts))
		coords = np.vstack((nfts,nfgs,nfzs))
		"""

		p = Pool(npools)
		t0 = time.time()

		final_chis = np.array((p.map(get_chis_comp, models)))
		p.terminate()
		III = np.argmin(final_chis)
		final_chis_o = final_chis.copy()
		"""
		final_chis = final_chis.reshape((len(rts),len(rgs),len(rzs)))
		zi = ndimage.map_coordinates(final_chis, coords,order=3, mode='nearest')
		I = np.where(zi==zi.min())[0][0]
		curt2,curg2,curz2 = interpolate.splev(nfts[I],tckit),interpolate.splev(nfgs[I],tckig),interpolate.splev(nfzs[I],tckiz)
		"""

		curt,curg,curz = models[III][0],models[III][1], models[III][2]
		if ik>0:
			IT = np.where((nrgs==curg)&(nrzs==curz))[0]
			tckt2 = interpolate.splrep(nrts[IT],final_chis_o[IT],k=3)
			temT = np.arange(lowt,upt,1.)
			curt3 = interpolate.splev(temT,tckt2)
			curt3 = temT[np.argmin(curt3)]
			if fixG == -1:
				IG = np.where((nrts==curt)&(nrzs==curz))[0]
				tckg2 = interpolate.splrep(nrgs[IG],final_chis_o[IG],k=3)
				temG = np.arange(lowg,upg,0.01)
				curg3 = interpolate.splev(temG,tckg2)
				curg3 = temG[np.argmin(curg3)]
			else:
				curg3 = fixG
			IZ = np.where((nrts==curt)&(nrgs==curg))[0]
			tckz2 = interpolate.splrep(nrzs[IZ],final_chis_o[IZ],k=3)
			temZ = np.arange(lowz,upz,0.01)
			curz3 = interpolate.splev(temZ,tckz2)
			#plot(nrzs[IZ],final_chis_o[IZ],'ro')
			#plot(temZ,curz3)
			#show()
			curz3 = temZ[np.argmin(curz3)]
		if ik >= 4:
			curt,curg,curz = curt3,curg3,curz3
		
		roundR = np.around(curr)
		if curr > roundR:
			roundR += .5
		else:
			roundR -= .5
		rots = np.arange(roundR-5,roundR+5,1.)
		I = np.where(rots >= 0.5)[0]
		rots = rots[I]
		past_curr = curr
		dplot=False
		#if ik == 3:
		#	dplot=True
		curr,NRV = get_vsini([curt,curg,curz,sc[0],sc[3],rots],dplot=dplot)
		if curr < 0.5:
			curr = 0.5
		sc[0] *= (1. - NRV/lux)
		RVTOT += NRV
		curr = np.around(curr,2)
		RVTOT = np.around(RVTOT,3)
		if printing:
			print 'Iteration', ik, ':', curt,curg,curz,curr,RVTOT,final_chis.min()
		#print fds
		#KK = np.where(np.absolute(all_currs-curr)<0.05)[0]
		dift,difg,difz,difr = np.absolute(all_curts - curt), np.absolute(all_curgs - curg), np.absolute(all_curzs - curz), np.absolute(all_currs - curr)
		KK1 = np.where(dift<10)[0]
		dift = np.zeros(len(dift))
		if len(KK1)>0:
			dift[KK1] = 1.
		KK1 = np.where(difg<0.03)[0]
		difg = np.zeros(len(difg))
		if len(KK1)>0:
			difg[KK1] = 1.
		KK1 = np.where(difz<0.02)[0]
		difz = np.zeros(len(difz))
		if len(KK1)>0:
			difz[KK1] = 1.
		#KK1 = np.where(difr<0.05)[0]
		#difr = np.zeros(len(difr))
		#if len(KK1)>0:
		#	difr[KK1] = 1.
		KK = dift+difg+difz
		#print KK
		KK = np.where(KK==3.)[0]
		if len(KK)>0 and ik > 3:
			break
		all_currs = np.hstack((all_currs,curr))
		all_curts = np.hstack((all_curts,curt))
		all_curgs = np.hstack((all_curgs,curg))
		all_curzs = np.hstack((all_curzs,curz))

		if not errors and use_masks:
			zi,zf = get_zones([curt,curg,curz,curr,RESI])
			ZO,ZI,ZF = get_zetas(zi,zf)
			ZO,ZI,ZF,BZO,BZI,BZF = the_good_zones(ZO,ZI,ZF,[curt,curg,curz,curr,RVTOT],th=5.0,zmin=zmin,limit=elim)
			zmask = np.zeros((sc.shape[1],sc.shape[2]))
			for i in range(len(ZI)):
				zmask[ZO[i],ZI[i]:ZF[i]] = 1.

			#zri,zrf = get_rot_zones([curt,curg,curz,curr,RESI])
			#ZRO,ZRI,ZRF = get_zetas(zri,zrf)
			#ZRO,ZRI,ZRF,BZRO,BZRI,BZRF = the_good_zones(ZRO,ZRI,ZRF,[curt,curg,curz,curr,RVTOT],th=3.0,zmin=400,limit=elim)
			#rot_mask = np.zeros((sc.shape[1],sc.shape[2]))
			#for i in range(len(ZRI)):
			#	rot_mask[ZRO[i],ZRI[i]:ZRF[i]] = 1.
			
			zmask2 = np.ones((sc.shape[1],sc.shape[2]))
			mask_bin = zmask.copy()
			rot_mask = mask_bin.copy()
			
		ik+=1
	#print curt,curg,curz
	if use_masks:
		if errors:
			print curt,curg,curz,curr,RVTOT
		else:
			print curt,curg,curz,curr,RVTOT,len(ZO),len(BZO)
	else:
		print curt,curg,curz,curr,RVTOT

	return np.array([curt,curg,curz,curr,RVTOT])

def get_zones(pars):
	lim = 0.05
	thrT,thrG,thrZ = 200,0.3,0.2
	rot2 = max(3.0,pars[3])
	R2 = pars[4]
	R = 127000.	
	rot = 3.0
	mod  = get_full_model(pars[0],pars[1],pars[2],rot,R)

	if pars[0]-thrT >= ts.min():
		mtl = get_full_model(pars[0]-thrT,pars[1],pars[2],rot,R)
	else:
		mtl = get_full_model(pars[0]+thrT,pars[1],pars[2],rot,R)
	if pars[0]+thrT <= ts.max():
		mtu = get_full_model(pars[0]+thrT,pars[1],pars[2],rot,R)
	else:
		mtu = get_full_model(pars[0]-thrT,pars[1],pars[2],rot,R)

	if pars[1]-thrG >= gs.min():
		mgl = get_full_model(pars[0],pars[1]-thrG,pars[2],rot,R)
	else:
		mgl = get_full_model(pars[0],pars[1]+thrG,pars[2],rot,R)
	if pars[1]+thrG <= gs.max():
		mgu = get_full_model(pars[0],pars[1]+thrG,pars[2],rot,R)
	else:
		mgu = get_full_model(pars[0],pars[1]-thrG,pars[2],rot,R)

	if pars[2]-thrZ >= zs.min():
		mzl = get_full_model(pars[0],pars[1],pars[2]-thrZ,rot,R)
	else:
		mzl = get_full_model(pars[0],pars[1],pars[2]+thrZ,rot,R)
	if pars[2]+thrZ <= zs.max():
		mzu = get_full_model(pars[0],pars[1],pars[2]+thrZ,rot,R)
	else:
		mzu = get_full_model(pars[0],pars[1],pars[2]-thrZ,rot,R)

	if library == 'P':
		w0 = mw[0]
		while w0 < mw[-1]:
			w1 = w0 + 100
			last = False
			if w1 + 100 > mw[-1]:
				w1 = mw[-1]
				last = True
			I = np.where((mw>=w0) & (mw<w1))[0]
			if last:
				I = np.hstack((I,I[-1]+1))
			mod[I] /= mod[I].max()
			mtl[I] /= mtl[I].max()
			mtu[I] /= mtu[I].max()
			mgl[I] /= mgl[I].max()
			mgu[I] /= mgu[I].max()
			mzl[I] /= mzl[I].max()
			mzu[I] /= mzu[I].max()
			w0 = w1
			if last:
				break

	difT = .5 * (np.absolute(mod-mtl) + np.absolute(mod-mtu))
	difG = .5 * (np.absolute(mod-mgl) + np.absolute(mod-mgu))
	difZ = .5 * (np.absolute(mod-mzl) + np.absolute(mod-mzu))

	#plot(mw,difT+difZ+difG)
	#plot(mw,np.zeros(len(mod))+0.05)
	#show()
	#print gfd
	diff = difT+difZ+difG
	I = np.where(diff>lim)[0]
	I2 = np.hstack((I[-1],I[:-1]))
	J = np.where(I-1 != I2)[0]
	JI = I[J]
	I2 = np.hstack((I[1:],I[0]))
	J = np.where(I+1!=I2)[0]
	JF = I[J]+1

	IF = np.where(JF >= len(mw))[0]
	if len(IF)>0:
		JF[IF] = len(mw)-1
	miw = mw[JI]
	mfw = mw[JF]
	devR = 300000./(2.355*R2)
	devr = max(3.0,rot2 * 0.66)
	dev = np.sqrt(devR**2+devr**2)
	extraw2 = 3*.5*(miw + mfw) * dev /300000.

	#extraw = .5*(miw + mfw) / (2.355*R2)
	#extraw2 = 3*.5*(miw + mfw) * rot2 * 0.66 / 300000.
	#print extraw
	#print extraw2
	#II = np.where(extraw2 < extraw)[0]
	#extraw2[II] = extraw[II]
	mmm = float(len(mw))/(mw[-1]-mw[0])
	extrap = (np.around(extraw2*mmm)).astype('int')
	#print extrap
	JI -= extrap
	IF = np.where(JI<0)[0]
	if len(IF)>0:
		JI[IF]=0
	JF += extrap
	IF = np.where(JF>=len(mw))[0]
	if len(IF)>0:
		JF[IF] = len(mw) - 1
	miw = mw[JI]
	mfw = mw[JF]

	"""
	i=0
	while i < len(miw)-1:
		if miw[i+1] < mfw[i]:
			mfw[i] = mfw[i+1]
			miw = np.delete(miw,i+1)
			mfw = np.delete(mfw,i+1)
		else:
			i+=1
	"""

	#for i in range(len(JF)):
	#	plot(mw[JI[i]:JF[i]],mod[JI[i]:JF[i]],'b')
	#show()
	#print J.shape
	#print len(miw)
	#print gfcd
	return miw,mfw	

def get_rot_zones(pars):
	lim = 0.01
	thrT,thrG,thrZ = 200,0.3,0.2
	rot2 = max(3.0,pars[3])
	R2 = pars[4]
	R = 127000.	
	rot = 3.0
	mod  = get_full_model(pars[0],pars[1],pars[2],rot,R)

	if pars[2]-thrZ >= zs.min():
		mzl = get_full_model(pars[0],pars[1],pars[2]-thrZ,rot,R)
	else:
		mzl = get_full_model(pars[0],pars[1],pars[2]+thrZ,rot,R)
	if pars[2]+thrZ <= zs.max():
		mzu = get_full_model(pars[0],pars[1],pars[2]+thrZ,rot,R)
	else:
		mzu = get_full_model(pars[0],pars[1],pars[2]-thrZ,rot,R)

	if library == 'P':
		w0 = mw[0]
		while w0 < mw[-1]:
			w1 = w0 + 100
			last = False
			if w1 + 100 > mw[-1]:
				w1 = mw[-1]
				last = True
			I = np.where((mw>=w0) & (mw<w1))[0]
			if last:
				I = np.hstack((I,I[-1]+1))
			mod[I] /= mod[I].max()
			mzl[I] /= mzl[I].max()
			mzu[I] /= mzu[I].max()
			w0 = w1
			if last:
				break

	difZ = .5 * (np.absolute(mod-mzl) + np.absolute(mod-mzu))

	diff = difZ
	I = np.where(diff>lim)[0]
	I2 = np.hstack((I[-1],I[:-1]))
	J = np.where(I-1 != I2)[0]
	JI = I[J]
	I2 = np.hstack((I[1:],I[0]))
	J = np.where(I+1!=I2)[0]
	JF = I[J]+1

	IF = np.where(JF >= len(mw))[0]
	if len(IF)>0:
		JF[IF] = len(mw)-1
	miw = mw[JI]
	mfw = mw[JF]
	devR = 300000./(2.355*R2)
	devr = rot2 * 0.66
	dev = np.sqrt(devR**2+devr**2)
	extraw2 = 3*.5*(miw + mfw) * dev /300000.

	#extraw = .5*(miw + mfw) / (2.355*R2)
	#extraw2 = 3*.5*(miw + mfw) * rot2 * 0.66 / 300000.
	#print extraw
	#print extraw2
	#II = np.where(extraw2 < extraw)[0]
	#extraw2[II] = extraw[II]
	mmm = float(len(mw))/(mw[-1]-mw[0])
	extrap = (np.around(extraw2*mmm)).astype('int')
	#print extrap
	JI -= extrap
	IF = np.where(JI<0)[0]
	if len(IF)>0:
		JI[IF]=0
	JF += extrap
	IF = np.where(JF>=len(mw))[0]
	if len(IF)>0:
		JF[IF] = len(mw) - 1
	miw = mw[JI]
	mfw = mw[JF]

	"""
	i=0
	while i < len(miw)-1:
		if miw[i+1] < mfw[i]:
			mfw[i] = mfw[i+1]
			miw = np.delete(miw,i+1)
			mfw = np.delete(mfw,i+1)
		else:
			i+=1
	"""
	nmiw,nmfw = [],[]
	for i in range(len(JF)):
		tf = mod[JI[i]:JF[i]]
		if tf.max() - tf.min() < 0.5 and tf.max() - tf.min() > 0.2:
			nmiw.append(miw[i])
			nmfw.append(mfw[i])
			#plot(mw[JI[i]:JF[i]],mod[JI[i]:JF[i]],'b')
	#show()
	#print J.shape
	#print len(miw)
	#print gfds
	nmiw,nmfw = np.array(nmiw),np.array(nmfw)
	return miw,mfw	

def make_mask(zi,zf):
	mask = np.zeros((sc.shape[1],sc.shape[2]))
	for i in range(sc.shape[1]):
		I = np.where((zi>sc[i,0]) & (zf < sc[i,-1]))[0]
		ZI,ZF = zi[I],zf[I]
		for o in range(len(ZI)):
			mask[ZI[o]:ZF[o]+1] = 1.
	return mask

def get_rats(ZO,ZI,ZF,pars):
	ords = []
	for i in range(sc.shape[1]):
		J1 = np.where(mw > sc[0,i,-1])[0]
		J2 = np.where(mw < sc[0,i,0])[0]
		if len(J1)>0 and len(J2)>0:
			ords.append(i)
	ords = np.array(ords)
	mf = get_full_model(pars[0],pars[1],pars[2],pars[3],RES_POW)
	tmodf = np.zeros((sc.shape[1],sc.shape[2]))
	tscif = np.zeros((sc.shape[1],sc.shape[2]))
	test_plot = np.zeros((4,sc.shape[1],sc.shape[2]))
	for i in ords:
		I = np.where((mw>sc[0,i,0]) & (mw<sc[0,i,-1]))[0]
		modw = mw[I]
		modf = mf[I]
		sciw = sc[0,i]
		scif = sc[3,i]/np.median(sc[3,i])
		modf = pixelization(modw,modf,sciw)
		#IMB = np.where(mask_bin[i]!=0)[0]
		#modf /= modf[IMB].mean()
		mscif = scipy.signal.medfilt(scif,11)
		rat = modf/mscif
		INF = np.where(mscif!=0)[0]
		coef = get_ratio(sciw[INF],rat[INF])
		scif = scif * np.polyval(coef,sciw)
		mscif = mscif * np.polyval(coef,sciw)
		coef = get_cont(sciw,mscif)
		scif = scif / np.polyval(coef,sciw)
		#plot(sciw,scif)
		coef = get_cont(sciw,modf)
		modf = modf / np.polyval(coef,sciw)
		#plot(sciw,modf)	
		tmodf[i] = modf
		tscif[i] = scif
		test_plot[0,i] = sc[0,i]
		test_plot[1,i] = scif
		test_plot[2,i] = modf
		test_plot[3,i] = mask_bin[i]
	#show()
	#print vcdx
	hdu = pyfits.PrimaryHDU(test_plot)
	os.system('rm example.fits')
	hdu.writeto('example.fits')
	rat = tscif/tmodf

	nejx = np.arange(100)/100.
	ratsout = []

	for i in range(len(ZI)):
		ejy = rat[ZO[i],ZI[i]:ZF[i]]
		ejx = np.arange(len(ejy))/float(len(ejy))
		tck = interpolate.splrep(ejx,ejy,k=3)
		if len(ratsout)==0:
			ratsout = interpolate.splev(nejx,tck)
		else:
			ratsout = np.vstack((ratsout,interpolate.splev(nejx,tck)))
			#plot(interpolate.splev(nejx,tck))
	#show()
	return ratsout

def res_f(p,s,m):
	return s - p[0]*m

def get_factors(ZO,ZI,ZF,pars):
	ords = np.sort(np.unique(ZO))
	mf = get_full_model(pars[0],pars[1],pars[2],pars[3],RES_POW)
	tmodf = np.zeros((sc.shape[1],sc.shape[2]))
	tscif = np.zeros((sc.shape[1],sc.shape[2]))
	test_plot = np.zeros((4,sc.shape[1],sc.shape[2]))
	for i in ords:
		I = np.where((mw>sc[0,i,0]) & (mw<sc[0,i,-1]))[0]
		modw = mw[I]
		modf = mf[I]
		sciw = sc[0,i]
		scif = sc[3,i]/np.median(sc[3,i])
		modf = pixelization(modw,modf,sciw)
		#IMB = np.where(mask_bin[i]!=0)[0]
		#modf /= modf[IMB].mean()
		mscif = scipy.signal.medfilt(scif,11)
		rat = modf/mscif
		INF = np.where(mscif!=0)[0]
		coef = get_ratio(sciw[INF],rat[INF])
		scif = scif * np.polyval(coef,sciw)
		mscif = mscif * np.polyval(coef,sciw)
		coef = get_cont(sciw,mscif)
		scif = scif / np.polyval(coef,sciw)
		#plot(sciw,scif)
		coef = get_cont(sciw,modf)
		modf = modf / np.polyval(coef,sciw)
		#plot(sciw,modf)	
		tmodf[i] = modf
		tscif[i] = scif
		test_plot[0,i] = sc[0,i]
		test_plot[1,i] = scif
		test_plot[2,i] = modf
		test_plot[3,i] = mask_bin[i]
	#show()
	#print vcdx
	#print 'HI'
	#hdu = pyfits.PrimaryHDU(test_plot)
	#os.system('rm example.fits')
	#hdu.writeto('example.fits')

	factors = np.array([])
	for i in range(len(ZI)):
		ejs = tscif[ZO[i],ZI[i]:ZF[i]]
		ejm = tmodf[ZO[i],ZI[i]:ZF[i]]
		ejs = -(ejs - 1)
		ejm = -(ejm - 1)
		f = optimize.leastsq(res_f,[1.],args=(ejs,ejm))[0]
		factors = np.append(factors,f)
	#hist(factors,30)
	#show()
	return factors

def the_good_zones(ZO,ZI,ZF,pars, th=5.0,zmin=400,limit=0.001):
	
	BZO,BZI,BZF = np.array([]),np.array([]),np.array([])
	ords = []
	for i in range(sc.shape[1]):
		J1 = np.where(mw > sc[0,i,-1])[0]
		J2 = np.where(mw < sc[0,i,0])[0]
		if len(J1)>0 and len(J2)>0:
			ords.append(i)
	ords = np.array(ords)
	mf = get_full_model(pars[0],pars[1],pars[2],pars[3],RES_POW)
	tmodf = np.zeros((sc.shape[1],sc.shape[2]))
	tscif = np.zeros((sc.shape[1],sc.shape[2]))

	mwav,mflx,wav,flx = mw.copy(),mf.copy(),sc[0].copy(),sc[3].copy()

	for i in ords:
		I = np.where((mwav>wav[i,0]) & (mwav<wav[i,-1]))[0]
		modw = mwav[I]
		modf = mflx[I]
		sciw = wav[i]
		scif = flx[i]/np.median(flx[i])
		modf = pixelization(modw,modf,sciw)
		modf /= modf.mean()
		mscif = scipy.signal.medfilt(scif,11)
		INF = np.where(mscif!=0)[0]
                I0 = np.where(mscif==0)[0]
                mscif[I0] = 1.
		rat = modf/mscif
		rat[I0] = 0.
		coef = get_ratio(sciw[INF],rat[INF])
		scif = scif * np.polyval(coef,sciw)
		mscif = mscif * np.polyval(coef,sciw)
		coef = get_cont(sciw,mscif)
		scif = scif / np.polyval(coef,sciw)
		coef = get_cont(sciw,modf)
		modf = modf / np.polyval(coef,sciw)
		#plot(sciw,scif)
		#plot(sciw,modf)
		tmodf[i] = modf
		tscif[i] = scif
	
	#show()
	#print gfds
	dif = tscif - tmodf
	
	res = dif**2
	resout = np.array([])
	difout = np.array([])
	deltzones = np.array([])
	lenzones = np.array([])
	factors = np.array([])
	#print len(ZO)
	for i in range(len(ZI)):
		resout = np.hstack((resout,np.sum(res[ZO[i],ZI[i]:ZF[i]])/float(ZF[i]-ZI[i])))
		difout = np.hstack((difout,np.mean(dif[ZO[i],ZI[i]:ZF[i]])))
		deltzones = np.hstack((deltzones,tscif[ZO[i],ZI[i]:ZF[i]].max()-tscif[ZO[i],ZI[i]:ZF[i]].min()))
		lenzones = np.hstack((lenzones,ZF[i]-ZI[i]+1))
		ejs = tscif[ZO[i],ZI[i]:ZF[i]]
		ejm = tmodf[ZO[i],ZI[i]:ZF[i]]
		ejs = -(ejs - 1)
		ejm = -(ejm - 1)
		f = optimize.leastsq(res_f,[1.],args=(ejs,ejm))[0]
		factors = np.append(factors,f)
	#print len(ZO)
	IB1 = np.where(deltzones<0.1)[0]
	IB2 = np.where(deltzones>1.5)[0]
	IB3 = np.where(factors>2)[0]
	IB4 = np.where(factors<0.5)[0]
	#IB3 = np.where(lenzones<5)[0]
	IB = np.unique(np.hstack((IB1,IB2,IB3,IB4)))
	IG = np.where((deltzones>=0.1)&(deltzones<=1.5)&(factors<=2)&(factors>=.5))[0]
	BZO = np.hstack((BZO,ZO[IB]))
	BZI = np.hstack((BZI,ZI[IB]))
	BZF = np.hstack((BZF,ZF[IB]))
	ZO  = ZO[IG]
	ZI = ZI[IG]
	ZF = ZF[IG]
	resout = resout[IG]
	difout = difout[IG]

	#hist(difout,30)
	#print len(difout)
	#show()

	
	meddif = np.median(difout)
	res  = difout - meddif
	dev  = np.sqrt(np.mean(res**2))
	I = np.where(np.absolute(res) > th*dev)[0]
	cond = True
	if len(I) == 0:
		cond = False
	while cond:
		im = np.argmax(res**2)
		difout = np.delete(difout,im)
		BZO = np.hstack((BZO,ZO[im]))
		BZI = np.hstack((BZI,ZI[im]))
		BZF = np.hstack((BZF,ZF[im]))
		ZO = np.delete(ZO,im)
		ZI = np.delete(ZI,im)
		ZF = np.delete(ZF,im)
		meddif = np.median(difout)
		res  = difout - meddif
		dev  = np.sqrt(np.mean(res**2))
		I = np.where(np.absolute(res) > th*dev)[0]
		if (len(I) == 0 and dev < limit/3.):
			cond = False
		if len(I) == 0 and len(difout)<zmin:
			cond = False
	
	#print len(difout)
	#hist(difout,30)
	#show()
	#print gfds
	"""
	devr = np.mean(resout)
	I = np.where(resout > th*devr)[0]
	cond = True
	if len(I) == 0:
		cond = False
	hist(difout,30)
	show()
	print gvfcd
	
	while cond:
		im = np.argmax(resout)
		resout = np.delete(resout,im)
		BZO = np.hstack((BZO,ZO[im]))
		BZI = np.hstack((BZI,ZI[im]))
		BZF = np.hstack((BZF,ZF[im]))
		ZO = np.delete(ZO,im)
		ZI = np.delete(ZI,im)
		ZF = np.delete(ZF,im)
		devr = np.mean(resout)
		I = np.where(resout > th*devr)[0]
		if (len(I) == 0 and th*devr < limit):
			cond = False
		if len(I) == 0 and len(resout)<zmin:
			cond = False
	#print th*devr
	#hist(resout,30)
	#show()
	#print gfds
	#print len(ZO)
	"""
	return ZO,ZI,ZF,BZO,BZI,BZF


def get_zetas(zi,zf):
	ZO,ZI,ZF = np.array([]),np.array([]),np.array([])
	for i in range(sc.shape[1]):
		J1 = np.where(mw > sc[0,i,-1])[0]
		J2 = np.where(mw < sc[0,i,0])[0]
		if len(J1)>0 and len(J2)>0:
			J3 = np.where((zi>sc[0,i,0]) & (zf<sc[0,i,-1]))
			zit,zft = zi[J3],zf[J3]
			tck = interpolate.splrep(sc[0,i],np.arange(sc.shape[2]))
			zpi = np.around(interpolate.splev(zit,tck)).astype('int')
			zpf = np.around(interpolate.splev(zft,tck)).astype('int')
			ZI = np.hstack((ZI,zpi))
			ZF = np.hstack((ZF,zpf))
			ZO = np.hstack((ZO,np.zeros(len(zpi))+i))
	return ZO,ZI,ZF
	
def get_final_zetas(ZO,ZI,ZF):
	NZO,NZI,NZF = np.array([]),np.array([]),np.array([])
	oos = np.unique(ZO)
	for o in oos:
		I = np.where(ZO==o)[0]
		tzo,tzi,tzf = ZO[I],ZI[I],ZF[I]
		
		i=0
		while i < len(tzi)-1:
			if tzi[i+1] < tzf[i]:
				tzf[i] = tzf[i+1]
				tzi = np.delete(tzi,i+1)
				tzf = np.delete(tzf,i+1)
			else:
				i+=1
		tzo = np.zeros(len(tzi))+o
		NZO = np.hstack((NZO,tzo))
		NZI = np.hstack((NZI,tzi))
		NZF = np.hstack((NZF,tzf))
	return NZO,NZI,NZF
	
		
def get_precise_parameters(spec,rough_parameters,RESI=120000.,ncores=14,trunc=0,fixG=-1,elim=0.1,zmin=500,nsim=50,efixG=-1):
	global sc, ords, RES_POW, mask_bin, rot_mask
	RES_POW = RESI
	rt,rg,rz,rr,rv = rough_parameters[0],rough_parameters[1],rough_parameters[2],rough_parameters[3],rough_parameters[4]

	sc = pyfits.getdata(spec)
	sc = sc[:,:,trunc:sc.shape[2]-trunc]
	sc[0] *= (1. - rv/lux)
	hd = pyfits.getheader(spec)

	zi,zf = get_zones([rt,rg,rz,rr,RESI])
	ZO,ZI,ZF = get_zetas(zi,zf)
	ZO,ZI,ZF,BZO,BZI,BZF = the_good_zones(ZO,ZI,ZF,rough_parameters,th=5.0,zmin=zmin,limit=elim)
	print len(ZO)
	mask = np.zeros((sc.shape[1],sc.shape[2]))
	mask_bin = np.ones((sc.shape[1],sc.shape[2]))
	for i in range(len(ZI)):
		mask[ZO[i],ZI[i]:ZF[i]] = 1.
	mask_bin = mask.copy()	
	rot_mask = mask_bin.copy()

	#for i in range(sc.shape[1]):
	#	plot(sc[0,i],sc[3,i])
	#	plot(sc[0,i],sc[3,i]*mask_bin[i])
	#show()
	#print gfds
	#zri,zrf = get_rot_zones([rt,rg,rz,rr,RESI])
	#ZRO,ZRI,ZRF = get_zetas(zri,zrf)
	#ZRO,ZRI,ZRF,BZRO,BZRI,BZRF = the_good_zones(ZRO,ZRI,ZRF,rough_parameters,th=3.0,zmin=400,limit=0.005)
	#rot_mask = np.zeros((sc.shape[1],sc.shape[2]))
	#for i in range(len(ZRI)):
	#	rot_mask[ZRO[i],ZRI[i]:ZRF[i]] = 1.

	#rpars = get_rough_pars(spec,RV0=0.,guess_vsini=rr,RESI=RESI,ncores=ncores,mask=mask,trunc=trunc,errors=False,use_masks=True)
	rpars = rough_parameters
	print len(ZO)
	ZO,ZI,ZF = get_final_zetas(ZO,ZI,ZF)
	print len(ZO)
	factors  = get_factors(ZO,ZI,ZF,rpars)
	#print hist(factors,20)
	#show()
	#print gfd
	#rats  = get_rats(ZO,ZI,ZF,rpars)
	
	#drats = {'rats':rats}
	#pickle.dump(drats,open('test_mask.pkl','w'))
	#print gfcdx
	print 'Starting Simulation ...'
	TPARS = []
	count = 0
	watch = []
	while count < nsim:
		nmask = mask.copy()
		for j in np.arange(len(ZI)):
			ran_val = np.random.randint(len(ZI))
			ran_val2 = np.random.randint(2)
			
			if sim_type == 'factors':
				factor = factors[ran_val]
				#print ZO[j],ZI[j],ZF[j],factor
				if ran_val2 == 0:
					nmask[ZO[j],ZI[j]:ZF[j]] *= factor
					watch.append(factor)
				else:
					nmask[ZO[j],ZI[j]:ZF[j]] /= factor
					watch.append(1./factor)
			else:
				rat = rats[ran_val]
				tck = interpolate.splrep(np.arange(len(rat))/float(len(rat)),rat,k=1)
				nejx = np.arange(ZF[j]-ZI[j])/float(ZF[j]-ZI[j])
				nejy = interpolate.splev(nejx,tck)
				if ran_val2 == 0:
					nmask[ZO[j],ZI[j]:ZF[j]] *= nejy
				else:
					nmask[ZO[j],ZI[j]:ZF[j]] /= nejy
		#hist(watch)
		#show()
		#print fdss
		#ords = np.unique(ZO)
		#global zmask
		#zmask = nmask
		#get_chis_comp([6000.,4.5,0.0,6.0])
		#print fds
		#for i in range(sc.shape[1]):
		#	plot(sc[0,i],nmask[i])
		#show()
		#print gfds
		try:
			if efixG !=-1:
				if fixG != -1:
					fixG2 = min(np.random.normal(loc=fixG,scale=efixG),5.0)
				else:
					print 'Problem with fixG'
			else:
				fixG2 = fixG
			tpars = get_rough_pars(spec,RV0=0.,guess_vsini=rpars[3],RESI=RESI,ncores=ncores,mask=nmask,trunc=trunc,errors=True,use_masks=True,fixG=fixG2)
			if len(TPARS)==0:
				TPARS = np.array(tpars)
			else:
				TPARS = np.vstack((TPARS,tpars))
			count +=1
		except:
			print 'Oops ...'
	"""
	et = np.around(np.sqrt(np.var(TPARS[:,0])))
	eg = np.around(np.sqrt(np.var(TPARS[:,1])),3)
	ez = np.around(np.sqrt(np.var(TPARS[:,2])),3)
	er = np.around(np.sqrt(np.var(TPARS[:,3])),3)
	ev = np.around(np.sqrt(np.var(TPARS[:,4])),3)

	hdu = pyfits.PrimaryHDU(TPARS)
	if os.access('results',os.F_OK) == False:
		os.system('mkdir results')
	date = str(time.localtime()).split('=')
	odat = '_'+date[1].split(',')[0]+'_'+date[2].split(',')[0]+'_'+date[3].split(',')[0]+\
		'_'+date[4].split(',')[0]+'_'+date[4].split(',')[0]+'_'+date[6].split(',')[0] + '_' 
	rout = 'results/' + spec.split('/')[-1][:-5] + odat + 'zaspe_out.fits'

	hdu.header.update('DATE',odat)
	hdu.header.update('SPEC',spec)
	hdu.header.update('RES',RESI)
	hdu.header.update('TRUNC',trunc)
	hdu.header.update('ELIM',elim)
	hdu.header.update('ZMIN',zmin)
	hdu.header.update('FIXG',False)
	if fixG != -1:
		hdu.header['FIXG'] = True
	hdu.header.update('TEFF',rt)
	hdu.header.update('LOGG',rg)
	hdu.header.update('FEH',rz)
	hdu.header.update('VSINI',rr)
	hdu.header.update('RV',rv)
	hdu.header.update('ETEFF',et)
	hdu.header.update('ELOGG',eg)
	hdu.header.update('EFEH',ez)
	hdu.header.update('EVSINI',er)
	hdu.header.update('ERV',ev)
	hdu.header.update('LIB',library)
	hdu.header.update('WI',wi)
	hdu.header.update('WF',wf)
	hdu.writeto(rout)
	"""
	return TPARS


