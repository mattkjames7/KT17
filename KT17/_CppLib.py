import numpy as np
import os
import subprocess
import ctypes

def _LibPath():
	'''
	Return a path to the C++ library
	
	Returns
	=======
	path : str
		path to the library file.
	
	'''
	return os.path.dirname(__file__)+"/__data/libkt/"

def _LibName(WithPath=False):
	'''
	Return the name of the library.
	
	Inputs
	======
	WithPath : bool
		If True then the fuull path to the library will be included.
		
	Returns
	=======
	libpath : str
		Library name.
	
	'''
	if WithPath:
		path = _LibPath()
	else:
		path = ''
		
	if(os.name == 'posix'):
		ext = "so"
	elif(os.name == 'nt'):
		ext = "dll"
	else:
		raise Exception("The Operating System is not supported")
	
	return path + 'libkt.' + ext


def _LibExists():
	'''
	Check if the library file exists.
	
	Returns
	=======
	exists : bool
		True if the file exists
	'''
	return os.path.isfile(_LibName(True))
	
def _GetLib():
	'''	
	Return an instance of the C++ library
	
	Returns
	=======
	lib : ctypes.CDLL
		C++ library containing the KT14/17 code
	'''
	fname = _LibName(True)
	
	try:
		lib = ctypes.CDLL(fname)
	except:
		print("Importing C++ library failed. Attempting recompilation...")
		_CompileSource()
		lib = ctypes.CDLL(fname)
		
	return lib


def _CompileSource():
	'''
	Attempt to recompile the library if needed.
	
	'''
	
	if(os.name=='posix'):
		#check if we need root or not!
		path = os.path.dirname(__file__)
		if '/usr/local/' in path:
			sudo = 'sudo '
		else:
			sudo = ''

		CWD = os.getcwd()
		os.chdir(os.path.dirname(__file__)+"/__data/libkt/")
		os.system(sudo+'make')
		os.chdir(CWD)
	elif(os.name=='nt'):
		CWD = os.getcwd()
		os.chdir(os.path.dirname(__file__)+"/__data/libkt/")
		compile = subprocess.Popen("compile.bat")
		compile.communicate()
		comperr = compile.returncode
		if(comperr==6):
			raise Exception("There is no GCC compiler in PATH. Unable to compile C source files.")
		if(comperr==7):
			raise Exception("There is no GFORTRAN compiler in PATH. Unable to compile FORTRAN source files.")
		if(comperr==8):
			raise Exception("An error occurred during compilation.")
		os.chdir(CWD)
	else:
		raise Exception("The Operating System is not supported")	
