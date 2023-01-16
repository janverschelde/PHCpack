"""
Exports the dimmension of the system of polynomials.
"""

import ctypes
from ctypes import c_int, c_double, pointer
from ctypes import create_string_buffer, sizeof
from version import getPHCmod

def setDoubleDimension(dim, vrblvl=0):
    """
    Sets the number of polynomials in double precision
    to the value of the first parameter dim.
    """
    if(vrblvl > 0):
        print("in setDoubleDimension, dim = ", dim)
    modPHCpack = getPHCmod()
    f = modPHCpack._ada_use_c2phc
    a = pointer(c_int(dim))
    b = pointer(c_int(0))
    c = pointer(c_double(0.0))
    v = c_int(vrblvl)
    r = f(23, a, b, c, v)
    if(vrblvl > 0):
        print('the return value of setDoubleDimension :', r)
    return r

def getDoubleDimension(vrblvl=0):
    """
    Returns the number of polynomials in double precision.
    """
    if(vrblvl > 0):
        print("in getDoubleDimension ...")
    modPHCpack = getPHCmod()
    f = modPHCpack._ada_use_c2phc
    a = pointer(c_int(0))
    b = pointer(c_int(0))
    c = pointer(c_double(0.0))
    v = c_int(vrblvl)
    r = f(22, a, b, c, v)
    if(vrblvl > 0):
        print('the return value of setDoubleDimension :', r)
        print('the retrieved dimension :', a[0])
    return a[0]

def showDimension():
    """
    Test on setting and getting the dimension.
    """
    setDoubleDimension(12345, 10)
    dim = getDoubleDimension(10)
    print('the dimension :', dim)

if __name__=="__main__":
    showDimension()
