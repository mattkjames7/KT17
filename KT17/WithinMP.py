import numpy as np
from ._CFunctions import _CWithinMP

def WithinMP(x, y, z, Rsm):
	'''
	Tests whether a position in MSM coordinates is within the 
	magnetopause.
	
	Inputs
	======
	x : float
		Scalar or array of x MSM coordinate(s) (Rm)
	y : float
		Scalar or array of y MSM coordinate(s) (Rm)
	z : float
		Scalar or array of z MSM coordinate(s) (Rm)
	Rsm : float
		Subsolar magnetopause standoff distance (Rm).
			
	Returns
	=======
	res : bool
		True if within the magnetopause.

	'''
	#Convert input variables to appropriate numpy dtype:
	if np.size(x) == 1:
		_x = np.float64(x)
		_y = np.float64(y)
		_z = np.float64(z)
		_Rsm = np.float64(Rsm)
		res = _CWithinMP(_x, _y, _z, _Rsm)
	else:
		n = np.size(x)
		res = np.zeros(n,dtype='bool')
		rsm = np.zeros(n,dtype='float64') + Rsm
		for i in range(0,n):
			_x = np.float64(x[i])
			_y = np.float64(y[i])
			_z = np.float64(z[i])
			_Rsm = np.float64(rsm[i])
			res[i] = _CWithinMP(_x, _y, _z, _Rsm)			

	return res
