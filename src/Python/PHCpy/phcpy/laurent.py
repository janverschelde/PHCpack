"""
Exports working with real powered Laurent homotopies.
"""
from ctypes import c_int32, c_double, pointer
from phcpy.version import get_phcfun

def initialize_series_coefficients(dims, vrblvl=0):
    """
    Sets the number of polynomials in double precision
    to the value of the first parameter dim.
    """
    if vrblvl > 0:
        print('in initialize_series_coefficients, dims :', dims)
    phc = get_phcfun(vrblvl-1)
    alen = pointer(c_int32(len(dims)))
    bpar = (c_int32 * len(dims))()
    for (idx, dim) in enumerate(dims):
        bpar[idx] = dim
    bdim = pointer(bpar)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> initialize_series_coefficients calls phc', end='')
    retval = phc(930, alen, bdim, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def test_initialization(vrblvl=0):
    """
    Tests initializing the series coefficients.
    """
    if vrblvl > 0:
        print("in test_initialization ...")
    dims = [3, 2, 4]
    fail = initialize_series_coefficients(dims, vrblvl)
    return fail

def test_laurent(vrblvl=0):
    """
    Tests operations in this module.
    The verbose level is defined by vrblvl.
    """
    if vrblvl > 0:
        print("in test_laurent ...")
    fail = test_initialization(vrblvl)
    if vrblvl > 0:
        if fail == 0:
            print('=> All tests on the laurent module passed.')
        else:
            print('Number of failed tests on laurent module :', fail)
    return fail

def main():
    """
    Runs tests on the laurent module.
    """
    lvl = 3
    fail = test_laurent(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__== '__main__':
    main()
