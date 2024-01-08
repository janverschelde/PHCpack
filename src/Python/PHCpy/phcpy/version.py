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
from ctypes import c_int32, c_double, pointer, sizeof
from struct import unpack

# relative location of the PHCpack library
LOCATION = "."

def get_phcfun_fromlib(vrblvl=0):
    """
    Returns the proper function according to the platform.
    For the correct execution, the file libPHCpack,
    with the extension .so, .dylib, or .dll must be present.
    """
    if vrblvl > 0:
        print('in get_phcfun_from_lib ...')
    if 'linux' in sys.platform:
        libphcpack = LOCATION + "/libPHCpack.so"
        phcpack = ctypes.CDLL(libphcpack)
        return phcpack._ada_use_c2phc
    if 'darwin' in sys.platform:
        libphcpack = LOCATION + "/libPHCpack.dylib"
        phcpack = ctypes.CDLL(libphcpack)
        return phcpack._ada_use_c2phc
    if 'win' in sys.platform:
        libphcpack = LOCATION + "/libPHCpack.dll"
        phcpack = ctypes.WinDLL(libphcpack, winmode=0)
        return phcpack._ada_use_c2phc
    print('The platform', sys.platform, 'is not supported.')
    return None

def get_phcfun(vrblvl=0):
    """
    Returns phcpy.phc, or if that not works,
    returns get_phcfun_fromlib()
    """
    if vrblvl > 0:
        print('in get_phcfun ...')
    try:
        import phcpy
        return phcpy.phc
    except:
        return get_phcfun_fromlib()

def int4a2str(data, vrblvl=0):
    """
    Given in data is an integer array stored as a ctypes string buffer,
    returns the string which represent the data.
    The '4a' refers to the 4-byte (or 32-bit) integer array.
    If vrblvl > 0, then results of intermediate steps are printed,
    otherwise, the function remains silent.
    """
    if vrblvl > 0:
        print('-> int4a2str, data is', data)
    szd = sizeof(data)
    if vrblvl > 0:
        print('-> int4a2str, size of the data is', szd)
    dim = szd//4
    if vrblvl > 0:
        print('-> int4a2str, number of characters :', dim)
    result = ""
    for k in range(dim):
        result = result + data[4*k].decode()
    if vrblvl > 0:
        print('-> int4a2str returns', result)
    return result

def str2int4a(data, vrblvl=0):
    """
    Given in data is a string, stored as a Python string,
    returns the 32-bit integer array representation,
    stored as a ctypes string buffer.
    If vrblvl > 0, then results of intermediate steps are printed,
    otherwise, the function remains silent.
    """
    if vrblvl > 0:
        print('-> str2int4a, data is', data)
    dim = len(data)
    if vrblvl > 0:
        print('-> str2int4a, size of the data is', dim)
    szd = 4*dim
    result = create_string_buffer(b"", szd)
    for k in range(dim):
        result[4*k] = data[k].encode()
    if vrblvl > 0:
        print('-> str2int4a returns', result)
    return result

def int4a2nbr(data, vrblvl=0):
    """
    Given in data is a Python list of integers,
    returns the encoding in a ctypes string buffer,
    for conversion into an 32-bit integer array.
    If vrblvl > 0, then results of intermediate steps are shown,
    otherwise the function remains silent.
    """
    if vrblvl > 0:
        print('-> int4a2nbr, data :', data)
    dim = len(data)
    szd = 4*dim
    if vrblvl > 0:
        print('-> int4a2nbr size of result :', szd)
    result = create_string_buffer(b"", szd)
    for k in range(dim):
        if data[k] < 256:
            result[4*k] = data[k] # no encode() because plain integer
        elif data[k] < 65536:
            result[4*k+1], result[4*k] = divmod(data[k], 256)
        elif data[k] < 16777216:
            result[4*k+2], rest = divmod(data[k], 65536)
            result[4*k+1], result[4*k] = divmod(rest, 256)
        else:
            result[4*k+3], rest = divmod(data[k], 16777216)
            result[4*k+2], rest = divmod(rest, 65536)
            result[4*k+1], result[4*k] = divmod(rest, 256)
    if vrblvl > 0:
        print('-> int4a2nbr return type :', type(result))
        print('-> int4a2nbr return value :', nbr2int4a(result))
    return result

def nbr2int4a(data, vrblvl=0):
    """
    Given in data is a 32-bit integer array,
    in a ctypes string buffer.
    Returns the list of python integers.
    """
    if vrblvl > 0:
        print('-> nbr2int4a, data type :', type(data))
    fmt = str(len(data)//4) + 'I'
    if vrblvl > 0:
        print('format :', fmt)
    result = list(unpack(fmt, data))
    if vrblvl > 0:
        print('-> nbr2int4a returns', result)
    return result

def version_string(vrblvl=0):
    """
    Returns the version string of PHCpack.
    If vrblvl > 0, then the conversions between strings and integer arrays
    are verified via the ctypes string buffer types.
    """
    if vrblvl > 0:
        print('in version_string ...')
    phc = get_phcfun(vrblvl)
    aaa = pointer(c_int32(0))
    name = create_string_buffer(b"", 30*4)
    ccc = pointer(c_double(0.0))
    if vrblvl > 0:
        retval = phc(999, aaa, name, ccc, 2)
    else:
        retval = phc(999, aaa, name, ccc, 0)
    if vrblvl > 0:
        print('return value :', retval)
        print('a :', aaa[0])
        print('sizeof(name) :', sizeof(name))
    result = int4a2str(name, vrblvl=vrblvl-1)
    if vrblvl > 0:
        print(result)
        print('version checks conversions ...')
        int4aresult = str2int4a(result, vrblvl=vrblvl-1)
        strres = int4a2str(int4aresult, vrblvl=vrblvl-1)
        print(strres)
    return result

def test_byte_strings(vrblvl=0):
    """
    Tests the conversion of a string into a string buffer
    and then backwards.
    """
    if vrblvl > 0:
        print('in test_byte_strings ...')
    greeting = 'Hello World!'
    if vrblvl > 0:
        print(greeting)
    coded = str2int4a(greeting, vrblvl-1)
    decoded = int4a2str(coded, vrblvl-1)
    if vrblvl > 0:
        print(decoded)
    if greeting == decoded:
        if vrblvl > 0:
            print('Test passed.')
        return 0
    if vrblvl > 0:
        print('Test failed!')
    return 1

def test_integer_encodings(vrblvl=0):
    """
    Tests the encoding of a list of integers as a ctypes string buffer.
    """
    if vrblvl > 0:
        print('in test_integer_encodings ...')
    first = [12, 1033, 129129, 20123543]
    print('L =', first)
    bfirst = int4a2nbr(first, vrblvl)
    second = nbr2int4a(bfirst, vrblvl)
    print('K =', second)
    if first == second:
        if vrblvl > 0:
            print('Test passed.')
        return 0
    if vrblvl > 0:
        print('Test failed!')
    return 1

def main():
    """
    Prints the version string and runs two tests.
    """
    lvl = 1
    print(version_string(lvl))
    fail = test_byte_strings(lvl)
    fail = fail + test_integer_encodings(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=="__main__":
    main()
