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

def power_dimension(vrblvl=0):
    """
    Returns the number of arrays of power vectors.
    """
    if vrblvl > 0:
        print('in power_dimension ...')
    phc = get_phcfun(vrblvl-1)
    dim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> power_dimension calls phc', end='')
    retval = phc(934, dim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return dim[0]

def coefficient_dimension(vrblvl=0):
    """
    Returns the number of arrays of coefficient vectors.
    """
    if vrblvl > 0:
        print('in coefficient_dimension ...')
    phc = get_phcfun(vrblvl-1)
    if vrblvl > 0:
        print('in power_dimension ...')
    phc = get_phcfun(vrblvl-1)
    dim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> coefficient_dimension calls phc', end='')
    retval = phc(935, dim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return dim[0]

def size_power_array(adx, vrblvl=0):
    """
    Returns the size of the power array with index adx.
    """
    if vrblvl > 0:
        print('in size_power_array, adx :', adx, '...')
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(adx))
    size = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> size_power_array calls phc', end='')
    retval = phc(936, aidx, size, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return size[0]

def size_coefficient_array(adx, vrblvl=0):
    """
    Returns the size of the coefficient array with index adx.
    """
    if vrblvl > 0:
        print('in size_coefficient_array, adx :', adx, '...')
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(adx))
    size = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> size_coefficient_array calls phc', end='')
    retval = phc(937, aidx, size, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return size[0]

def size_power_vector(adx, vdx, vrblvl=0):
    """
    Returns the size of the power vector at position vdx in array adx.
    """
    if vrblvl > 0:
        print('in size_power_vector, adx :', adx, ', vdx :', vdx, '...')
    phc = get_phcfun(vrblvl-1)

def size_coefficient_vector(adx, vdx, vrblvl=0):
    """
    Returns the size of the coefficient vector at position vdx in array adx.
    """
    if vrblvl > 0:
        print('in size_coefficient_vector, adx :', adx, ', vdx :', vdx, '...')
    phc = get_phcfun(vrblvl-1)

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
    retval1 = phc(942, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval1)
    if vrblvl > 0:
        print('-> clear_series_terms calls phc', end='')
    retval2 = phc(943, aaa, bbb, ccc, vrb)
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

def test_dimensions(vrblvl=0):
    """
    Tests retrieving the dimensions.
    """
    if vrblvl > 0:
        print("in test_dimensions ...")
    powdim = power_dimension(vrblvl)
    print('the power dimension :', powdim)
    cffdim = coefficient_dimension(vrblvl)
    print('the coefficient dimension :', cffdim)
    powsizes = [size_power_array(i+1,vrblvl) for i in range(powdim)]
    print('size power arrays :', powsizes)
    cffsizes = [size_power_array(i+1,vrblvl) for i in range(cffdim)]
    print('size coefficient arrays :', cffsizes)
    return 0

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
    fail = fail + test_dimensions(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__== '__main__':
    main()
