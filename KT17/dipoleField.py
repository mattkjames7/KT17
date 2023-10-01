import numpy as np
from ._CFunctions import _Cdipole

def dipoleField(x,y,z):

    n = np.int32(np.size(x))
    if n == 1:
        _x = np.zeros(1,dtype='float64') + x
        _y = np.zeros(1,dtype='float64') + y
        _z = np.zeros(1,dtype='float64') + z
    else:
        _x = np.array(x,dtype='float64')
        _y = np.array(y,dtype='float64')
        _z = np.array(z,dtype='float64')

    _bx = np.zeros(n,dtype='float64')
    _by = np.zeros(n,dtype='float64')
    _bz = np.zeros(n,dtype='float64')

    _Cdipole(n,_x,_y,_z,_bx,_by,_bz)

    return _bx,_by,_bz