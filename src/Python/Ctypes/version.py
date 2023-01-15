"""
Tests the retrieval of the version string of the PHCpack library.
The function getPHCmod() returns the library module for the three
supported platforms.
"""

import ctypes
from ctypes import c_int, c_double, pointer
from ctypes import create_string_buffer, sizeof

def getPHCmod():
    """
    Returns the proper module according to the platform.
    """
    import sys
    if 'linux' in sys.platform:
        LIBPHCPACK = "../../lib/libPHCpack.so"
        modPHCpack = ctypes.CDLL(LIBPHCPACK)
        return modPHCpack
    elif 'darwin' in sys.platform:
        LIBPHCPACK = "../../lib/libPHCpack.dylib"
        modPHCpack = ctypes.CDLL(LIBPHCPACK)
        return modPHCpack
    elif 'win' in sys.platform:
        LIBPHCPACK = "../../lib/libPHCpack"
        modPHCpack = ctypes.WinDLL(LIBPHCPACK, winmode=0)
        return modPHCpack
    else:
        print('The platform', sys.platform, 'is not supported.')
        return None

def version(verbose=True):
    """
    Returns the version string of PHCpack.
    """
    modPHCpack = getPHCmod()
    f = modPHCpack._ada_use_c2phc
    a = pointer(c_int(0))
    b = create_string_buffer(30*4)
    c = pointer(c_double(0.0))
    if verbose:
        r = f(999, a, b, c, 2)
    else:
        r = f(999, a, b, c, 0)
    if verbose:
        print('return value :', r)
        print('a :', a[0])
        print('sizeof(b) :', sizeof(b))
    result = ""
    for k in range(30):
        result = result + b[4*k].decode()
    if verbose:
         print(result)
    return result

if __name__=="__main__":
    print(version(False))
