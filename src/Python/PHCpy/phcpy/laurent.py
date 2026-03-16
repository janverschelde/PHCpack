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
    bval = (c_int32 * len(dims))()
    for (idx, dim) in enumerate(dims):
        bval[idx] = dim
    bdim = pointer(bval)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> initialize_series_coefficients calls phc', end='')
    retval1 = phc(930, alen, bdim, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval1)
    if vrblvl > 0:
        print('-> initialize_series_coefficients calls phc', end='')
    retval2 = phc(931, alen, bdim, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval2)
    return retval1 + retval2

def set_series_term(adx, vdx, pwr, cff, vrblvl=0):
    """
    Sets the powers and coefficients of one term, at the position adx
    in the array and vector with index vdx in the array.
    The powers in pwr are doubles, whereas cff contains complex numbers.
    """
    if vrblvl > 0:
        print('in set_series_term, adx :', adx, ', vdx :', vdx)
        print('pwr : ', pwr)
        print('cff : ', cff)
    phc = get_phcfun(vrblvl-1)
    alen = pointer(c_int32(len(pwr)))
    bval = (c_int32 * 2)()
    bval[0] = adx
    bval[1] = vdx
    bdim = pointer(bval)
    cval = (c_double * len(pwr))()
    for (idx, val) in enumerate(pwr):
        cval[idx] = val
    ccc = pointer(cval)
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_series_term calls phc', end='')
    retval1 = phc(932, alen, bdim, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval1)
    cval = (c_double * len(pwr))()
    for (idx, val) in enumerate(cff):
        cval[2*idx] = val.real
        cval[2*idx+1] = val.imag
    ccc = pointer(cval)
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_series_term calls phc', end='')
    retval2 = phc(933, alen, bdim, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval2)
    return retval1 + retval2

def get_series_term(adx, vdx, vrblvl=0):
    """
    Gets the powers and coefficients of one term, at the position adx
    in the array and vector with index vdx in the array.
    """
    if vrblvl > 0:
        print('in set_series_term, adx :', adx, ', vdx :', vdx)
        print('pwr : ', pwr)
        print('cff : ', cff)

def clear_series_terms(vrblvl=0):
    """
    Deallocates the space occupied by the real powers
    and the complex coefficients.
    """
    if vrblvl > 0:
        print('in clear_series_terms ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_series_terms calls phc', end='')
    retval1 = phc(936, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval1)
    if vrblvl > 0:
        print('-> clear_series_terms calls phc', end='')
    retval2 = phc(937, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval2)
    return retval1 + retval2

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
