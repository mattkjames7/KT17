import numpy as np
from ._CFunctions import _CModelField
from .ct import ctDoublePtrPtr

def ModelField(x, y, z, **kwargs):
	'''
	Returns the KT14/KT17 model field at one or more positions within
	Mercury's magnetosphere.
	
	Inputs
	=======
	x : float
		Scalar or array of x MSM coordinate(s) (Rm)
	y : float
		Scalar or array of y MSM coordinate(s) (Rm)
	z : float
		Scalar or array of z MSM coordinate(s) (Rm)
	
	**kwargs
	========
	Rsm : float
		Subsolar magnetopause standoff distance (Rm).
	t1 : float
		Disk current magntitude
	t2 : float
		Quasi-harris sheet magnitude
	Rsun : float
		Distance of Mercury from the Sun (AU).
	DistIndex : float
		Anderson et al 2013 disturbance index (0.0-97.0).
		
	NOTE : this model will default to using the KT14 paper's model
	parameters (Rsm = 1.42, t1 = 7.37, t2 = 2.16). The KT17 parameters
	(Rsun and DistIndex) are used to calculate the KT14 ones at 
	different activity levels. If both sets of parameters are provided
	then KT14 ones will be prioritized.
	
	Returns
	=======
	Bx : float
			Magnetic field vectors in MSM coordinate system (nT).
	By : float
			Magnetic field vectors in MSM coordinate system (nT).
	Bz : float
			Magnetic field vectors in MSM coordinate system (nT).
	'''

	#Convert input variables to appropriate numpy dtype:
	_n = np.int32(np.size(x))
	_x = np.array([x]).astype("float64").flatten()
	_y = np.array([y]).astype("float64").flatten()
	_z = np.array([z]).astype("float64").flatten()
	_Bx = np.zeros(_n,dtype="float64")
	_By = np.zeros(_n,dtype="float64")
	_Bz = np.zeros(_n,dtype="float64")
	
	#figure out what parameters we have and which to pass to the C code
	kt14 = np.array(['Rsm' in kwargs,'t1' in kwargs,'t2' in kwargs])
	kt17 = np.array(['Rsun' in kwargs,'DistIndex' in kwargs])
	
	if kt14.all():
		#use just kt14 parameters
		_Params = np.zeros((_n,3),dtype='float64')
		_Params[:,0] = kwargs['Rsm']
		_Params[:,1] = kwargs['t1']
		_Params[:,2] = kwargs['t2']
	elif kt17.all():
		#use only kt17 parameters
		_Params = np.zeros((_n,2),dtype='float64')
		_Params[:,0] = kwargs['Rsun']
		_Params[:,1] = kwargs['DistIndex']
	elif kt14.any():
		#use some kt14 parameters + defaults
		_Params = np.zeros((_n,3),dtype='float64')
		_Params[:,0] = kwargs.get('Rsm',1.42)
		_Params[:,1] = kwargs.get('t1',7.37)
		_Params[:,2] = kwargs.get('t2',2.16)		
	elif kt17.any():
		#use some kt17 parameters + defaults
		_Params = np.zeros((_n,2),dtype='float64')
		_Params[:,0] = kwargs.get('Rsun',0.427)
		_Params[:,1] = kwargs.get('DistIndex',50.0)		
	else:
		#use defaults
		_Params = np.zeros((_n,3),dtype='float64')
		_Params[:,0] = 1.42
		_Params[:,1] = 7.37
		_Params[:,2] = 2.16		

	_lP,_nP = _Params.shape

	#do some ctypes magic
	_Params2D = ctDoublePtrPtr(_Params)
	
	#call the model code
	_CModelField(_n, _x, _y, _z,_lP, _nP, _Params2D, _Bx, _By, _Bz)

	return (_Bx,_By,_Bz)
