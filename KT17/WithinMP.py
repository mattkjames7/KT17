import numpy as np
from ._CFunctions import _CWithinMP

###### File created automatically using PopulateCtypes ######

def WithinMP(x, y, z, Rsm):
	'''
	Tests whether a position in MSM coordinates is within the 
	magnetopause.
	
	Inputs
	======
	x,y,z: Scalars containing a position in MSM coordinates.
	Rsm: Subsolar distance of the magnetopause in R_m.
	
	Returns
	=======
	Boolean - True if within the magnetopause.
	
	'''
	#Convert input variables to appropriate numpy dtype:
	_x = np.float64(x)
	_y = np.float64(y)
	_z = np.float64(z)
	_Rsm = np.float64(Rsm)
	res = _CWithinMP(_x, _y, _z, _Rsm)

	return res
