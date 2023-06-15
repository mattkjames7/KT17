import os
import numpy as np
from .TraceField import TraceField
from scipy.interpolate import InterpolatedUnivariateSpline

ocbfile = os.path.dirname(__file__)+'/__data/OCB.bin'

def GetOCB(**kwargs):
	'''
	Tries to calculate where the open-closed field line boundary is for 
	MLT 0.0-24.0 then saves it to a file.
	
	'''


	dtype = [
		('MLT','float32'),
		('LTN','float32'),
		('LTS','float32'),
		('MLat','float32'),
		('LatN','float32'),
		('LatS','float32')
	]
	
	ocb = np.recarray(24,dtype=dtype)
	ocb.fill(np.nan)

	T = np.arange(24)*15.0
	L = np.arange(320)*0.25+10.0

	kwargs['alpha'] = kwargs.get('alpha',[])

	TraceOut = []

	for i in range(0,24):
		print('{0} of {1}'.format(i+1,24))
		Z = np.sin(np.pi*L/180.0)
		XY = np.cos(np.pi*L/180.0)
		X = -XY*np.cos(np.pi*T[i]/180.0)
		Y = -XY*np.sin(np.pi*T[i]/180.0)


		Tr = TraceField(X,Y,Z,**kwargs)
		TraceOut.append(Tr)
		
		#j = np.where((Tr.Lshell == np.nanmax(Tr.Lshell)) & np.isfinite(Tr.FlLen))[0]

		closed = np.isfinite(Tr.MLatN) & np.isfinite(Tr.MLatS) & np.isfinite(Tr.Lshell)
		inds = np.where(closed)[0]
		lmax = np.argmax(Tr.Lshell[inds])
		j = inds[lmax]

		ocb.MLT[i] = T[i]/15
		if j.size > 0:
			#mlt[i]=trace.mltn[j]
			#mlat[i]=trace.mlatn[j]
			
			ocb.MLat[i] = L[j]
			ocb.LatN[i] = Tr.LatN[j]
			ocb.LatS[i] = Tr.LatS[j]
			ocb.LTN[i] = Tr.LTN[j]
			ocb.LTS[i] = Tr.LTS[j]
		else:
			ocb.LTN[i] = T[i]/15
			ocb.LTS[i] = T[i]/15

	if (not os.path.isfile(ocbfile)) or kwargs.get('Overwrite',False):
		f = open(ocbfile,'wb')
		ocb.MLT.tofile(f)
		ocb.LTN.tofile(f)
		ocb.LTS.tofile(f)
		ocb.MLat.tofile(f)
		ocb.LatN.tofile(f)
		ocb.LatS.tofile(f)
		f.close()	

	if kwargs.get('ReturnTraces',False):
		return ocb,TraceOut
	else:
		return ocb

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
	
	dtype = [
		('MLT','float32'),
		('LTN','float32'),
		('LTS','float32'),
		('MLat','float32'),
		('LatN','float32'),
		('LatS','float32')
	]
	data = np.recarray(24,dtype=dtype)
	f = open(ocbfile,'rb')
	data.MLT = np.fromfile(f,dtype='float32',count=24)
	data.LTN = np.fromfile(f,dtype='float32',count=24)
	data.LTS = np.fromfile(f,dtype='float32',count=24)
	data.MLat = np.fromfile(f,dtype='float32',count=24)
	data.LatN = np.fromfile(f,dtype='float32',count=24)
	data.LatS = np.fromfile(f,dtype='float32',count=24)
	f.close()
	
	return data
