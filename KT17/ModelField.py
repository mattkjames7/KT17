import numpy as np
from ._CFunctions import _Ckt17Barray

###### File created automatically using PopulateCtypes ######

def ModelField(x, y, z, Params=[1.42,7.37,2.16]):

	#Convert input variables to appropriate numpy dtype:
	_n = np.int32(np.size(x))
	_x = np.array([x]).astype("float64").flatten()
	_y = np.array([y]).astype("float64").flatten()
	_z = np.array([z]).astype("float64").flatten()
	_Bx = np.zeros(_n,dtype="float64")
	_By = np.zeros(_n,dtype="float64")
	_Bz = np.zeros(_n,dtype="float64")
	
	Par = np.array(Params)	
	shP = Par.shape
	if np.size(shP) == 2:
		if shP[1] == _n:
			Par = Par.T
			shP = Par.shape
		_nP = np.int32(shP[1]) # number of parameters for each trace
	else:
		_nP = np.int32(shP[0])
	_Plen = np.int32(np.size(Par)) # total number of parameters
	
	_Params = Par.flatten().astype('float64')
	
	_Ckt17Barray(_n, _x, _y, _z, _Bx, _By, _Bz, _nP, _Plen, _Params)

	return (_Bx,_By,_Bz)
