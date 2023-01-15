import ctypes
from ctypes import c_int, c_double, pointer
from ctypes import create_string_buffer, sizeof

LIBPHCPACK = "../lib/libPHCpack"

def version(verbose=True):
    """
    Returns the version string of PHCpack.
    """
    d = ctypes.WinDLL(LIBPHCPACK, winmode=0)
    f = d._ada_use_c2phc
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
