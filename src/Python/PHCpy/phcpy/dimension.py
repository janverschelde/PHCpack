"""
Exports the dimension of the system of polynomials.
Setting the dimension allocates memory to store polynomials.
"""
from ctypes import c_int32, c_double, pointer
from phcpy.version import get_phcfun

def set_double_dimension(dim, vrblvl=0):
    """
    Sets the number of polynomials in double precision
    to the value of the first parameter dim.
    """
    if vrblvl > 0:
        print('in set_double_dimension, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_dimension calls phc', end='')
    retval = phc(23, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_double_dimension(dim, vrblvl=0):
    """
    Sets the number of polynomials in double double precision
    to the value of the first parameter dim.
    """
    if vrblvl > 0:
        print('in set_double_double_dimension, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_double_dimension calls phc', end='')
    retval = phc(333, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_quad_double_dimension(dim, vrblvl=0):
    """
    Sets the number of polynomials in quad double precision
    to the value of the first parameter dim.
    """
    if vrblvl > 0:
        print('in set_quad_double_dimension, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_quad_double_dimension calls phc', end='')
    retval = phc(383, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_laurent_dimension(dim, vrblvl=0):
    """
    Sets the number of Laurent polynomials in double precision
    to the value of the first parameter dim.
    """
    if vrblvl > 0:
        print('in set_double_laurent_dimension, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_laurent_dimension calls phc', end='')
    retval = phc(123, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_double_laurent_dimension(dim, vrblvl=0):
    """
    Sets the number of Laurent polynomials in double double precision
    to the value of the first parameter dim.
    """
    if vrblvl > 0:
        print('in set_double_double_laurent_dimension, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_double_laurent_dimension calls phc', end='')
    retval = phc(553, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_quad_double_laurent_dimension(dim, vrblvl=0):
    """
    Sets the number of Laurent polynomials in quad double precision
    to the value of the first parameter dim.
    """
    if vrblvl > 0:
        print('in set_quad_double_laurent_dimension, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_quad_double_laurent_dimension calls phc', end='')
    retval = phc(563, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def get_double_dimension(vrblvl=0):
    """
    Returns the number of polynomials in double precision.
    """
    if vrblvl > 0:
        print("in get_double_dimension ...")
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> get_double_dimension calls phc', end='')
    retval = phc(22, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the retrieved dimension :', adim[0])
    return adim[0]

def get_double_double_dimension(vrblvl=0):
    """
    Returns the number of polynomials in double double precision.
    """
    if vrblvl > 0:
        print("in get_double_double_dimension ...")
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> get_double_double_dimension calls phc', end='')
    retval = phc(332, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the retrieved dimension :', adim[0])
    return adim[0]

def get_quad_double_dimension(vrblvl=0):
    """
    Returns the number of polynomials in quad double precision.
    """
    if vrblvl > 0:
        print("in get_quad_double_dimension ...")
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> get_quad_double_dimension calls phc', end='')
    retval = phc(382, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the retrieved dimension :', adim[0])
    return adim[0]

def get_double_laurent_dimension(vrblvl=0):
    """
    Returns the number of Laurent polynomials in double precision.
    """
    if vrblvl > 0:
        print("in get_double_laurent_dimension ...")
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> get_double_laurent_dimension calls phc', end='')
    retval = phc(122, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the retrieved dimension :', adim[0])
    return adim[0]

def get_double_double_laurent_dimension(vrblvl=0):
    """
    Returns the number of Laurent polynomials in double double precision.
    """
    if vrblvl > 0:
        print("in get_double_double_laurent_dimension ...")
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> get_double_double_laurent_dimension calls phc', end='')
    retval = phc(552, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the retrieved dimension :', adim[0])
    return adim[0]

def get_quad_double_laurent_dimension(vrblvl=0):
    """
    Returns the number of Laurent polynomials in quad double precision.
    """
    if vrblvl > 0:
        print("in get_quad_double_laurent_dimension ...")
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> get_quad_double_laurent_dimension calls phc', end='')
    retval = phc(562, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the retrieved dimension :', adim[0])
    return adim[0]

def set_seed(seed, vrblvl=0):
    """
    Sets the seed for the random number generators to seed.
    """
    if vrblvl > 0:
        print('in set_seed, seed :', seed)
    phc = get_phcfun(vrblvl-1)
    aseed = pointer(c_int32(seed))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_seed calls phc', end='')
    retval = phc(998, aseed, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def get_seed(vrblvl=0):
    """
    Returns the seed used to generate random numbers.
    """
    if vrblvl > 0:
        print('in get_seed ...')
    phc = get_phcfun(vrblvl-1)
    aseed = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> get_seed calls phc', end='')
    retval = phc(997, aseed, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the retrieved seed :', aseed[0])
    return aseed[0]

def get_core_count(vrblvl=0):
    """
    Returns the number of available cores.
    """
    if vrblvl > 0:
        print('in get_core_count ...')
    phc = get_phcfun(vrblvl-1)
    acores = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> get_core_count calls phc', end='')
    retval = phc(994, acores, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the core count :', acores[0])
    return acores[0]

def test_double_dimension(vrblvl=0):
    """
    Tests setting and getting the dimension
    for systems in double precision.
    The verbose level is defined by vrblvl.
    """
    if vrblvl > 0:
        print("in test_double_dimension ...")
    set_double_dimension(12345, vrblvl)
    dim = get_double_dimension(vrblvl)
    if vrblvl > 0:
        print('the dimension :', dim)
    fail = int(dim != 12345)
    return fail

def test_double_double_dimension(vrblvl=0):
    """
    Tests setting and getting the dimension
    for systems in double double precision.
    The verbose level is defined by vrblvl.
    """
    if vrblvl > 0:
        print("in test_double_double_dimension ...")
    set_double_double_dimension(12345, vrblvl)
    dim = get_double_double_dimension(vrblvl)
    if vrblvl > 0:
        print('the dimension :', dim)
    fail = int(dim != 12345)
    return fail

def test_quad_double_dimension(vrblvl=0):
    """
    Tests setting and getting the dimension
    for systems in quad double precision.
    The verbose level is defined by vrblvl.
    """
    if vrblvl > 0:
        print("in test_quad_double_dimension ...")
    set_quad_double_dimension(12345, vrblvl)
    dim = get_quad_double_dimension(vrblvl)
    if vrblvl > 0:
        print('the dimension :', dim)
    fail = int(dim != 12345)
    return fail

def test_double_laurent_dimension(vrblvl=0):
    """
    Tests setting and getting the dimension
    for Laurent systems in double precision.
    The verbose level is defined by vrblvl.
    """
    if vrblvl > 0:
        print("in test_double_laurent_dimension ...")
    set_double_laurent_dimension(12345, vrblvl)
    dim = get_double_laurent_dimension(vrblvl)
    if vrblvl > 0:
        print('the dimension :', dim)
    fail = int(dim != 12345)
    return fail

def test_double_double_laurent_dimension(vrblvl=0):
    """
    Tests setting and getting the dimension
    for Laurent systems in double double precision.
    The verbose level is defined by vrblvl.
    """
    if vrblvl > 0:
        print("in test_double_double_laurent_dimension ...")
    set_double_double_laurent_dimension(12345, vrblvl)
    dim = get_double_double_laurent_dimension(vrblvl)
    if vrblvl > 0:
        print('the dimension :', dim)
    fail = int(dim != 12345)
    return fail

def test_quad_double_laurent_dimension(vrblvl=0):
    """
    Tests setting and getting the dimension
    for Laurent systems in quad double precision.
    The verbose level is defined by vrblvl.
    """
    if vrblvl > 0:
        print("in test_quad_double_laurent_dimension ...")
    set_quad_double_laurent_dimension(12345, vrblvl)
    dim = get_quad_double_laurent_dimension(vrblvl)
    if vrblvl > 0:
        print('the dimension :', dim)
    fail = int(dim != 12345)
    return fail

def test_dimension(vrblvl=0):
    """
    Tests setting and getting the dimension.
    The verbose level is defined by vrblvl.
    """
    if vrblvl > 0:
        print("in test_dimension ...")
    fail = test_double_dimension(vrblvl)
    fail = fail + test_double_double_dimension(vrblvl)
    fail = fail + test_quad_double_dimension(vrblvl)
    fail = fail + test_double_laurent_dimension(vrblvl)
    fail = fail + test_double_double_laurent_dimension(vrblvl)
    fail = fail + test_quad_double_laurent_dimension(vrblvl)
    if vrblvl > 0:
        if fail == 0:
            print('=> All tests on set/get dimension passed.')
        else:
            print('Number of failed tests on set/get dimension :', fail)
    return fail

def test_seed(vrblvl=0):
    """
    Tests the setting and getting of the seed.
    The verbose level is defined by vrblvl.
    """
    if vrblvl > 0:
        print('in test_seed ...')
    seed = get_seed(vrblvl)
    if vrblvl > 0:
        print('The current seed :', seed)
        print('Setting the seed to 12345 ...')
    set_seed(12345, vrblvl)
    seed = get_seed(vrblvl)
    if vrblvl > 0:
        print('The current seed :', seed)
    fail = int(seed != 12345)
    if vrblvl > 0:
        if fail == 0:
            print('=> Tests on set/get seed passed.')
        else:
            print('Test on set/get seed failed!')
    return fail

def test_core_count(vrblvl=0):
    """
    Tests if the number of available cores is positive.
    The verbose level is defined by vrblvl.
    """
    if vrblvl > 0:
        print('in test_core_count ...')
    cores = get_core_count(vrblvl)
    if vrblvl > 0:
        print('The number of available cores :', cores)
    fail = int(cores <= 0)
    if vrblvl > 0:
        if fail == 0:
            print('=> Test on get core count passed.')
        else:
            print('Test on get core count failed!')
    return fail

def main():
    """
    Runs tests on set/get dimension and set/get seed.
    """
    lvl = 1
    fail = test_dimension(lvl)
    fail = fail + test_seed(lvl)
    fail = fail + test_core_count(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__== '__main__':
    main()
