#!/usr/bin/python

import numpy as np
import scipy
import Cfunctions
from scipy import optimize
from scipy import interpolate
from scipy.interpolate import griddata

def conv2(lam,flu,ac,bc,rot):
	
	lux = 299792.458
	
	WAV = np.array([3530,4860,6260,7670,8350])
		
	"""Determining Limb Darkening Coefficients"""
	tck = interpolate.splrep(WAV,ac,k=3,s=0)
	A = interpolate.splev(lam,tck,der=0)
	tck = interpolate.splrep(WAV,bc,k=3,s=0)
	B = interpolate.splev(lam,tck,der=0)
		
	F = flu.copy()
	
	MIN = lam * ( 1. - float(rot) / lux )
	MAX = lam * ( 1. + float(rot) / lux )
	
	lar = len(F)
	
	if rot > 0.0:
		H = Cfunctions.Conv(lam.astype("double"),flu.astype("double"),F.astype("double"), \
		    A.astype("double"),B.astype("double"),MIN.astype("double"),MAX.astype("double"),rot,lar)
	else:
		H = flu
	
	Fnu = np.array(H)

	return Fnu

def get_ldcoef(TEFF, LOG_G, FEH, MT, path = 'qua_coe_sloan2.dat'):
	
	f = open(path,'r')
	lines = f.readlines()
	f.close()
	if MT == 1:
		FEH =0.
	else:
		if FEH < -2.5:	
			FEH = 2.5
		elif FEH > 0.5:
			FEH = 0.5
	if TEFF < 3500:
		TEFF = 3500
	if LOG_G > 5.0:
		LOG_G = 5.0

	th = np.loadtxt(path)
	n1 = th[:,0]
	qc = th[:,1]
	tur = th[:,2]
	lg = th[:,3]
	te = th[:,4]
	me = th[:,5]
	su = th[:,6]
	sg = th[:,7]
	sr = th[:,8]
	si = th[:,9]
	sz = th[:,10]
	
	I = np.where(tur == MT)[0]
	qc = qc[I]
	te = te[I]
	lg = lg[I]
	me = me[I]
	su = su[I]
	sg = sg[I]
	sr = sr[I]
	si = si[I]
	sz = sz[I]

	I = np.where( (np.absolute(lg - LOG_G) <= 1) & (np.absolute(me - FEH) <= 1) & (np.absolute(te - TEFF) <= 500) )[0]
	qc = qc[I]
	te = te[I]
	lg = lg[I]
	me = me[I]
	su = su[I]
	sg = sg[I]
	sr = sr[I]
	si = si[I]
	sz = sz[I]


	I = np.where(qc == 0.)[0]
	W = np.where(qc == 1.)[0]

	if MT == 2.:
		apoints = np.zeros((len(I),3))
		apoints[:,0] = te[I]
		apoints[:,1] = lg[I]
		apoints[:,2] = me[I]
		bpoints = np.zeros((len(W),3))
		bpoints[:,0] = te[W]
		bpoints[:,1] = lg[W]
		bpoints[:,2] = me[W]
	
	

		asu = griddata(apoints, su[I], (TEFF, LOG_G, FEH), method='linear')
		asg = griddata(apoints, sg[I], (TEFF, LOG_G, FEH), method='linear')
		asr = griddata(apoints, sr[I], (TEFF, LOG_G, FEH), method='linear')	
		asi = griddata(apoints, si[I], (TEFF, LOG_G, FEH), method='linear')
		asz = griddata(apoints, sz[I], (TEFF, LOG_G, FEH), method='linear')

		bsu = griddata(bpoints, su[W], (TEFF, LOG_G, FEH), method='linear')
		bsg = griddata(bpoints, sg[W], (TEFF, LOG_G, FEH), method='linear')
		bsr = griddata(bpoints, sr[W], (TEFF, LOG_G, FEH), method='linear')	
		bsi = griddata(bpoints, si[W], (TEFF, LOG_G, FEH), method='linear')
		bsz = griddata(bpoints, sz[W], (TEFF, LOG_G, FEH), method='linear')

	else:
		apoints = np.zeros((len(I),2))
		apoints[:,0] = te[I]
		apoints[:,1] = lg[I]
		
		bpoints = np.zeros((len(W),2))
		bpoints[:,0] = te[W]
		bpoints[:,1] = lg[W]
			

		asu = griddata(apoints, su[I], (TEFF, LOG_G), method='linear')
		asg = griddata(apoints, sg[I], (TEFF, LOG_G), method='linear')
		asr = griddata(apoints, sr[I], (TEFF, LOG_G), method='linear')	
		asi = griddata(apoints, si[I], (TEFF, LOG_G), method='linear')
		asz = griddata(apoints, sz[I], (TEFF, LOG_G), method='linear')

		bsu = griddata(bpoints, su[W], (TEFF, LOG_G), method='linear')
		bsg = griddata(bpoints, sg[W], (TEFF, LOG_G), method='linear')
		bsr = griddata(bpoints, sr[W], (TEFF, LOG_G), method='linear')	
		bsi = griddata(bpoints, si[W], (TEFF, LOG_G), method='linear')
		bsz = griddata(bpoints, sz[W], (TEFF, LOG_G), method='linear')

	ac = np.array([asu,asg,asr,asi,asz])
	bc = np.array([bsu,bsg,bsr,bsi,bsz])

	return ac,bc

def get_ldcoef2(TEFF, LOG_G, FEH, MT, path = 'qua_coe_sloan2.dat'):
	d = np.loadtxt(path)
	ldn = d[:,1].astype('int')

	ms = d[:,2]

	ts = d[:,4]
	ats = np.unique(ts)
	gs = d[:,3]
	ags = np.unique(gs)
	zs = d[:,5]
	azs = np.unique(zs)
	
	TT = ats[np.argmin((ats - TEFF)**2)]
	GG = ags[np.argmin((ags - LOG_G)**2)]
	ZZ = azs[np.argmin((azs - FEH)**2)]
	cond = True
	I = np.where((ts==TT) & (gs==GG) & (zs==ZZ))[0]
	if len(I)>0:
		cond = False
	deltz = 0.1
	while cond:
		ZZ1 = ZZ - deltz
		ZZ2 = ZZ + deltz
		I = np.where((ts==TT) & (gs==GG) & (zs==ZZ1))[0]
		if len(I)>0:
			cond = False
			ZZ = ZZ1
			break
		I = np.where((ts==TT) & (gs==GG) & (zs==ZZ2))[0]
		if len(I)>0:
			cond = False
			ZZ = ZZ2
			break		
		
	ams = np.unique(ms[I])
	MM = ams[np.argmin((ams - MT)**2)]

	Ia = np.where((ts==TT) & (gs==GG) & (zs==ZZ) & (ms==MM) & (ldn==0))[0][0]
	Ib = np.where((ts==TT) & (gs==GG) & (zs==ZZ) & (ms==MM) & (ldn==1))[0][0]

	ac = d[Ia][6:]
	bc = d[Ib][6:]

	return ac,bc
