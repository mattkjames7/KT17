import numpy as np
from .ct import c_char_p,c_bool,c_bool_ptr,c_int,c_int_ptr
from .ct import c_float,c_float_ptr,c_double,c_double_ptr,c_double_ptr_ptr
from ._CppLib import  getLibFilename,addWindowsSearchPaths
import os
import platform
import ctypes

modulePath = os.path.dirname(__file__)


#try loading the C++ library
try:
	if platform.system() == 'Darwin':
		cwd = os.getcwd()
		os.chdir(modulePath + '/__data/libkt/lib/')
		libkt17 = ctypes.CDLL(getLibFilename(True))
		os.chdir(cwd)
	elif platform.system() == 'Windows':
		addWindowsSearchPaths()
		libkt17 = ctypes.CDLL(getLibFilename())
	else:
		libkt17 = ctypes.CDLL(getLibFilename())
except:
	print('Importing '+getLibFilename(isShort=True)+' failed, please reinstall')
	raise SystemError



_CWithinMP = libkt17.WithinMP
_CWithinMP.argtypes = [	c_double,	#x MSM
						c_double,	#y MSM
						c_double,	#z MSM
						c_double]	#Rsm
_CWithinMP.restype = c_bool


_CWithinMPRT = libkt17.WithinMPRT
_CWithinMPRT.argtypes = [	c_double,	#r MSM
							c_double,	#theta MSM
							c_double]	#Rsm
_CWithinMPRT.restype = c_bool


_CModelField = libkt17.ModelField
_CModelField.argtypes = [	c_int,				#n
							c_double_ptr,		#x MSM
							c_double_ptr,		#y MSM
							c_double_ptr,		#z MSM
							c_int,				#lP
							c_int,				#nP
							c_double_ptr_ptr,	#Params	
							c_double_ptr,		#Bx
							c_double_ptr,		#By
							c_double_ptr]		#Bz
_CModelField.restype = None

_CTraceField = libkt17.TraceField
_CTraceField.argtypes = [	c_int,				#n
							c_double_ptr,		#x0
							c_double_ptr,		#y0
							c_double_ptr,		#z0
							c_int,				#nP	
							c_double_ptr,		#P0
							c_double_ptr,		#P1
							c_double_ptr,		#P2
							c_bool,				#BoundMP
							c_double,			#BoundTail
							c_int,				#BoundSurface
							c_int,				#MaxLen
							c_double,			#MaxStep
							c_double,			#InitStep
							c_double,			#MinStep
							c_double,			#ErrMax
							c_double,			#Delta
							c_bool,				#Verbose
							c_int,				#TraceDir
							c_int_ptr,			#nstep
							c_double_ptr_ptr,	#x
							c_double_ptr_ptr,	#y
							c_double_ptr_ptr,	#z
							c_double_ptr_ptr,	#Bx
							c_double_ptr_ptr,	#By
							c_double_ptr_ptr,	#Bz
							c_double_ptr_ptr,	#Rmsm
							c_double_ptr_ptr,	#Rmso
							c_double_ptr_ptr,	#S
							c_double_ptr_ptr,	#Rnorm
							c_double_ptr_ptr,	#FP
							c_int,				#nalpha
							c_double_ptr,		#alpha
							c_double_ptr]		#halpha
_CTraceField.restype = None

_Cdipole = libkt17.dipoleWrapper
_Cdipole.restype = None
_Cdipole.argtypes = [
    c_int,
    c_double_ptr,
    c_double_ptr,
    c_double_ptr,
    c_double_ptr,
    c_double_ptr,
    c_double_ptr,
]