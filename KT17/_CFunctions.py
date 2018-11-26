import numpy as np
import ctypes as ct
import os
import platform

Arch = platform.architecture()[0]
if Arch == '64bit':
	libkt17 = ct.CDLL(os.path.dirname(__file__)+"/__data/libkt17/libkt17_64.so")
elif Arch == '32bit':
	libkt17 = ct.CDLL(os.path.dirname(__file__)+"/__data/libkt17/libkt17_32.so")
else:
	print('Architecture ({:s}) not found'.format(Arch))



_CWithinMP = libkt17.WithinMP
_CWithinMP.argtypes = [ct.c_double, ct.c_double, ct.c_double, ct.c_double]
_CWithinMP.restype = ct.c_bool

_Ckt17DipoleField = libkt17.kt17DipoleField
_Ckt17DipoleField.argtypes = [ct.c_double, ct.c_double, ct.c_double, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")]
_Ckt17DipoleField.restype = None

_Ckt17DipoleShield = libkt17.kt17DipoleShield
_Ckt17DipoleShield.argtypes = [ct.c_double, ct.c_double, ct.c_double, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")]
_Ckt17DipoleShield.restype = None

_Ckt17DipoleB = libkt17.kt17DipoleB
_Ckt17DipoleB.argtypes = [ct.c_double, ct.c_double, ct.c_double, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")]
_Ckt17DipoleB.restype = None

_Ckt17DiskThickness = libkt17.kt17DiskThickness
_Ckt17DiskThickness.argtypes = [ct.c_double, ct.c_double]
_Ckt17DiskThickness.restype = ct.c_double

_Ckt17DiskField = libkt17.kt17DiskField
_Ckt17DiskField.argtypes = [ct.c_double, ct.c_double, ct.c_double, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")]
_Ckt17DiskField.restype = None

_Ckt17DiskShield = libkt17.kt17DiskShield
_Ckt17DiskShield.argtypes = [ct.c_double, ct.c_double, ct.c_double, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")]
_Ckt17DiskShield.restype = None

_Ckt17DiskB = libkt17.kt17DiskB
_Ckt17DiskB.argtypes = [ct.c_double, ct.c_double, ct.c_double, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_double]
_Ckt17DiskB.restype = None

_Ckt17QuasiHarrisThickness = libkt17.kt17QuasiHarrisThickness
_Ckt17QuasiHarrisThickness.argtypes = [ct.c_double, ct.c_double]
_Ckt17QuasiHarrisThickness.restype = ct.c_double

_Ckt17QuasiHarrisField = libkt17.kt17QuasiHarrisField
_Ckt17QuasiHarrisField.argtypes = [ct.c_double, ct.c_double, ct.c_double, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")]
_Ckt17QuasiHarrisField.restype = None

_Ckt17QuasiHarrisShield = libkt17.kt17QuasiHarrisShield
_Ckt17QuasiHarrisShield.argtypes = [ct.c_double, ct.c_double, ct.c_double, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")]
_Ckt17QuasiHarrisShield.restype = None

_Ckt17QuasiHarrisB = libkt17.kt17QuasiHarrisB
_Ckt17QuasiHarrisB.argtypes = [ct.c_double, ct.c_double, ct.c_double, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_double]
_Ckt17QuasiHarrisB.restype = None

_Ckt17B = libkt17.kt17B
_Ckt17B.argtypes = [ct.c_double, ct.c_double, ct.c_double, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_double, ct.c_double, ct.c_double]
_Ckt17B.restype = None

_Ckt17Barray = libkt17.kt17Barray
_Ckt17Barray.argtypes = [ct.c_int, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_int, ct.c_int, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")]
_Ckt17Barray.restype = None

_Ckt17 = libkt17.kt17
_Ckt17.argtypes = [ct.c_double, ct.c_double, ct.c_double, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_double, ct.c_double]
_Ckt17.restype = None

# _Ckt17array = libkt17.kt17array
# _Ckt17array.argtypes = [ct.c_int, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")]
# _Ckt17array.restype = None

_CUnitVector = libkt17.UnitVector
_CUnitVector.argtypes = [ct.c_double, ct.c_double, ct.c_double, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")]
_CUnitVector.restype = None

_CReverseElements = libkt17.ReverseElements
_CReverseElements.argtypes = [np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_int]
_CReverseElements.restype = None

_CInsideTraceLims = libkt17.InsideTraceLims
_CInsideTraceLims.argtypes = [ct.c_double, ct.c_double, ct.c_double,np.ctypeslib.ndpointer(ct.c_bool,flags="C_CONTIGUOUS"), ct.c_int, ct.c_double]
_CInsideTraceLims.restype = ct.c_bool

_CNorthSouthFLs = libkt17.NorthSouthFLs
_CNorthSouthFLs.argtypes = [np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_int, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_int,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_int,flags="C_CONTIGUOUS")]
_CNorthSouthFLs.restype = None

_Clinterp = libkt17.linterp
_Clinterp.argtypes = [ct.c_double, ct.c_double, ct.c_double, ct.c_double, ct.c_double]
_Clinterp.restype = ct.c_double

_CEqFootprint = libkt17.EqFootprint
_CEqFootprint.argtypes = [np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_int, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_int, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")]
_CEqFootprint.restype = None

_Cmin = libkt17.min
_Cmin.argtypes = [ct.c_double, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_int]
_Cmin.restype = ct.c_double

_CPlanetFootprints = libkt17.PlanetFootprints
_CPlanetFootprints.argtypes = [np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_int, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")]
_CPlanetFootprints.restype = None

_CFieldLineLength = libkt17.FieldLineLength
_CFieldLineLength.argtypes = [np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_int, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")]
_CFieldLineLength.restype = None

_CRvecs = libkt17.Rvecs
_CRvecs.argtypes = [ct.c_double, ct.c_double, ct.c_double, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_double, ct.c_double, ct.c_double, ct.c_double]
_CRvecs.restype = None

_Ckt17StepRKM = libkt17.kt17StepRKM
_Ckt17StepRKM.argtypes = [ct.c_double, ct.c_double, ct.c_double, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_double, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_double, ct.c_double, ct.c_double]
_Ckt17StepRKM.restype = None

_Ckt17Trace = libkt17.kt17Trace
_Ckt17Trace.argtypes = [ct.c_double, ct.c_double, ct.c_double, ct.c_int, ct.c_double, ct.c_double, ct.c_int, ct.c_int, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_int,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")]
_Ckt17Trace.restype = None

_Ckt17MultiTrace = libkt17.kt17MultiTrace
_Ckt17MultiTrace.argtypes = [np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"),np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"),np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_int, ct.c_int, ct.c_double, ct.c_double, ct.c_int, ct.c_int, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_int,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")]
_Ckt17MultiTrace.restype = None


# _Ckt17Trace = libkt17.kt17Trace
# _Ckt17Trace.argtypes = [ct.c_double, ct.c_double, ct.c_double, ct.c_int, ct.c_double, ct.c_double, np.ctypeslib.ndpointer(ct.c_int,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_int, ct.c_double, ct.c_double, ct.c_double]
# _Ckt17Trace.restype = None

# _Ckt17TraceScaled = libkt17.kt17TraceScaled
# _Ckt17TraceScaled.argtypes = [ct.c_double, ct.c_double, ct.c_double, ct.c_int, ct.c_double, ct.c_double, np.ctypeslib.ndpointer(ct.c_int,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_int, ct.c_double, ct.c_double]
# _Ckt17TraceScaled.restype = None

# _Ckt17MultiTrace = libkt17.kt17MultiTrace
# _Ckt17MultiTrace.argtypes = [np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_int, ct.c_int, ct.c_double, ct.c_double, np.ctypeslib.ndpointer(ct.c_int,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_int, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")]
# _Ckt17MultiTrace.restype = None

# _Ckt17MultiTraceScaled = libkt17.kt17MultiTraceScaled
# _Ckt17MultiTraceScaled.argtypes = [np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_int, ct.c_int, ct.c_double, ct.c_double, np.ctypeslib.ndpointer(ct.c_int,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), ct.c_int, np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")]
# _Ckt17MultiTraceScaled.restype = None

