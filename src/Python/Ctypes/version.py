"""
Tests the retrieval of the version string of the PHCpack library,
using ctypes, in a platform independent manner.
The function getPHCmod() returns the library module for the three
supported platforms.  The functions int4a2str() and str2int4a()
convert integer arrays to strings and string to integer arrays,
using the ctypes string buffer type.
The function int4a2nbr() converts a list of integers to the
string buffer representation of the 32-bit integer array version.
"""

import ctypes
import sys
from ctypes import create_string_buffer
from ctypes import c_int, c_double, pointer, sizeof

# relative location of the PHCpack library
LOCATION = "../../lib"

def getPHCmod():
    """
    Returns the proper module according to the platform.
    For the correct execution, the file libPHCpack,
    with the extension .so, .dylib, or .dll must be present.
    """
    if 'linux' in sys.platform:
        LIBPHCPACK = LOCATION + "/libPHCpack.so"
        modPHCpack = ctypes.CDLL(LIBPHCPACK)
        return modPHCpack
    if 'darwin' in sys.platform:
        LIBPHCPACK = LOCATION + "/libPHCpack.dylib"
        modPHCpack = ctypes.CDLL(LIBPHCPACK)
        return modPHCpack
    if 'win' in sys.platform:
        LIBPHCPACK = LOCATION + "/libPHCpack.dll"
        modPHCpack = ctypes.WinDLL(LIBPHCPACK, winmode=0)
        return modPHCpack
    print('The platform', sys.platform, 'is not supported.')
    return None

def int4a2str(data, verbose=False):
    """
    Given in data is an integer array stored as a ctypes string buffer,
    returns the string which represent the data.
    The '4a' refers to the 4-byte (or 32-bit) integer array.
    If verbose, then results of intermediate steps are printed,
    otherwise, the function remains silent.
    """
    if verbose:
        print('-> int4a2str, data is', data)
    szd = sizeof(data)
    if verbose:
        print('-> int4a2str, size of the data is', szd)
    dim = szd//4
    if verbose:
        print('-> int4a2str, number of characters :', dim)
    result = ""
    for k in range(dim):
        result = result + data[4*k].decode()
    if verbose:
        print('-> int4a2str returns', result)
    return result

def str2int4a(data, verbose=False):
    """
    Given in data is a string, stored as a Python string,
    returns the 32-bit integer array representation,
    stored as a ctypes string buffer.
    If verbose, then results of intermediate steps are printed,
    otherwise, the function remains silent.
    """
    if verbose:
        print('-> str2int4a, data is', data)
    dim = len(data)
    if verbose:
        print('-> str2int4a, size of the data is', dim)
    szd = 4*dim
    result = create_string_buffer(szd)
    for k in range(dim):
        result[4*k] = data[k].encode()
    if verbose:
        print('-> str2int4a returns', result)
    return result

def int4a2nbr(data, verbose=False):
    """
    Given in data is a Python list of integers,
    returns the encoding in a ctypes string buffer,
    for conversion into an 32-bit integer array.
    If verbose, then results of intermediate steps are shown,
    otherwise the function remains silent.
    """
    if verbose:
        print('-> int4a2nbr, data :', data)
    dim = len(data)
    szd = 4*dim
    if verbose:
        print('-> int4a2nbr size of result :', szd)
    result = create_string_buffer(szd)
    for k in range(dim):
        result[4*k] = data[k] # no encode() because plain integer
    if verbose:
        print('-> int4anbr returns', result)
    return result

def version(verbose=True):
    """
    Returns the version string of PHCpack.
    If verbose, then the conversions between strings and integer arrays
    are verified via the ctypes string buffer types.
    """
    modPHCpack = getPHCmod()
    phc = modPHCpack._ada_use_c2phc
    aaa = pointer(c_int(0))
    name = create_string_buffer(30*4)
    ccc = pointer(c_double(0.0))
    if verbose:
        retval = phc(999, aaa, name, ccc, 2)
    else:
        retval = phc(999, aaa, name, ccc, 0)
    if verbose:
        print('return value :', retval)
        print('a :', aaa[0])
        print('sizeof(name) :', sizeof(name))
    result = int4a2str(name, True)
    if verbose:
        print(result)
        print('version checks conversions ...')
        int4aresult = str2int4a(result, True)
        strres = int4a2str(int4aresult)
        print(strres)
    return result

if __name__=="__main__":
    print(version(True))
