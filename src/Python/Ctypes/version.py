"""
Tests the retrieval of the version string of the PHCpack library.
"""

import ctypes
from ctypes import c_int, c_double, pointer
from ctypes import create_string_buffer, sizeof
import sys

if 'linux' in sys.platform:
    LIBPHCPACK = "../../lib/libPHCpack.so"
    modPHCpack = ctypes.CDLL(LIBPHCPACK)
elif 'darwin' in sys.platform:
    LIBPHCPACK = "../../lib/libPHCpack.dylib"
    modPHCpack = ctypes.CDLL(LIBPHCPACK)
elif 'win' in sys.platform:
    LIBPHCPACK = "../../lib/libPHCpack"
    modPHCpack = ctypes.WinDLL(LIBPHCPACK, winmode=0)
else:
    print('The platform', sys.platform, 'is not supported.')

def version(verbose=True):
    """
    Returns the version string of PHCpack.
    """
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

print(version(False))
