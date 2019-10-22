import os
import numpy as np
from .TraceField import TraceField
from scipy.interpolate import InterpolatedUnivariateSpline

ocbfile = os.path.dirname(__file__)+'/__data/OCB.bin'

def GetOCB():
	'''
	Tries to calculate where the open-closed field line boundary is for 
	MLT 0.0-24.0 then saves it to a file.
	
	'''


	mlt = np.zeros(24,dtype='float32')-1
	lctn = np.zeros(24,dtype='float32')-1
	lcts = np.zeros(24,dtype='float32')-1
	mlat = np.zeros(24,dtype='float32')-1
	latn = np.zeros(24,dtype='float32')-1
	lats = np.zeros(24,dtype='float32')-1

	T = np.arange(24)*15.0
	L = np.arange(280)*0.25+10.0

	for i in range(0,24):
		print('{0} of {1}'.format(i+1,24))
		Z = np.sin(np.pi*L/180.0)
		XY = np.cos(np.pi*L/180.0)
		X = -XY*np.cos(np.pi*T[i]/180.0)
		Y = -XY*np.sin(np.pi*T[i]/180.0)
		
		Tr = TraceField(X,Y,Z)
		
		j = np.where((Tr.Lshell == np.nanmax(Tr.Lshell)) & np.isfinite(Tr.FlLen))[0]
		mlt[i] = T[i]/15
		if j.size > 0:
			#mlt[i]=trace.mltn[j]
			#mlat[i]=trace.mlatn[j]
			
			mlat[i] = L[j]
			latn[i] = Tr.GlatN[j]
			lats[i] = Tr.GlatS[j]
			lctn[i] = Tr.GltN[j]
			lcts[i] = Tr.GltS[j]
		else:
			lctn[i] = T[i]/15
			lcts[i] = T[i]/15

	f = open(ocbfile,'wb')
	mlt.tofile(f)
	lctn.tofile(f)
	lcts.tofile(f)
	mlat.tofile(f)
	latn.tofile(f)
	lats.tofile(f)
	f.close()	

def ReadOCBFile():
	'''
	Reads the file created in GetOCB()
	
	Returns
	=======
	numpy.recarray with the following fields:-
	
	Mlt:	Magnetic local time of the boundary.
	LctN:	Local time in the northern hemisphere.
	LctS:	Local time in the southern hemisphere.
	Mlat:	Invariant latitude of the boundary.
	LatN:	Latitude of the northern boundary.
	LatS:	Latitude of the southern boundary.
	'''
	if not os.path.isfile(ocbfile):
		GetOCB()
	
	dtype = [('Mlt','float32'),('LctN','float32'),('Lcts','float32'),
			('Mlat','float32'),('LatN','float32'),('Lats','float32')]
	data = np.recarray(24,dtype=dtype)
	f = open(ocbfile,'rb')
	data.Mlt = np.fromfile(f,dtype='float32',count=24)
	data.LctN = np.fromfile(f,dtype='float32',count=24)
	data.LctS = np.fromfile(f,dtype='float32',count=24)
	data.Mlat = np.fromfile(f,dtype='float32',count=24)
	data.LatN = np.fromfile(f,dtype='float32',count=24)
	data.LatS = np.fromfile(f,dtype='float32',count=24)
	f.close()
	
	return data
