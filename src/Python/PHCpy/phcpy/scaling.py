"""
Some application areas turn out systems with coefficients that
have a large variation in magnitudes.  Scaling a system reduces
the variability in the magnitude of the coefficients and leads
to more accurate solutions.
"""
from ctypes import c_int32, c_double, pointer
from phcpy.version import get_phcfun
from phcpy.polynomials import set_double_system, get_double_system
from phcpy.polynomials import set_double_double_system
from phcpy.polynomials import get_double_double_system
from phcpy.polynomials import set_quad_double_system
from phcpy.polynomials import get_quad_double_system

def double_scale_polynomials(mode, dim, vrblvl=0):
    """
    Applies equation and variable scaling to the system set in
    double precision.  The mode is either 0, 1, or 2, with
    0 : only scaling of the equations,
    1 : variable scaling without variability reduction,
    2 : variable scaling with variability reduction.
    The number dim is the dimension of the system.
    Returns the scaling coefficients (as complex numbers)
    and the estimated condition number of the scaling problem 
    if the mode > 0, as the last double in the coefficients.
    """
    if vrblvl > 0:
        print('in double_scale_polynomials, mode :', mode)
    phc = get_phcfun(vrblvl-1)
    amode = pointer(c_int32(mode))
    bbb = pointer(c_int32(0))
    lencffs = 4*dim + 2
    cffs = (c_double * lencffs)()
    pcff = pointer(cffs)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_scale_polynomials calls phc', end='')
    retval = phc(590, amode, bbb, pcff, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    allcffs = pcff[0][:lencffs] 
    result = []
    for idx in range(lencffs):
        result.append(allcffs[idx])
    if vrblvl > 0:
        print('the scaling coefficients :', result)
    return result

def double_double_scale_polynomials(mode, dim, vrblvl=0):
    """
    Applies equation and variable scaling to the system set in
    double double precision.  The mode is either 0, 1, or 2, with
    0 : only scaling of the equations,
    1 : variable scaling without variability reduction,
    2 : variable scaling with variability reduction.
    The number dim is the dimension of the system.
    Returns the scaling coefficients (as complex numbers)
    and the estimated condition number of the scaling problem 
    if the mode > 0, as the last double in the coefficients.
    """
    if vrblvl > 0:
        print('in double_double_scale_polynomials, mode :', mode)
    phc = get_phcfun(vrblvl-1)
    amode = pointer(c_int32(mode))
    bbb = pointer(c_int32(0))
    lencffs = 8*dim + 4
    cffs = (c_double * lencffs)()
    pcff = pointer(cffs)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_scale_polynomials calls phc', end='')
    retval = phc(591, amode, bbb, pcff, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    allcffs = pcff[0][:lencffs] 
    result = []
    for idx in range(lencffs):
        result.append(allcffs[idx])
    if vrblvl > 0:
        print('the scaling coefficients :', result)
    return result

def quad_double_scale_polynomials(mode, dim, vrblvl=0):
    """
    Applies equation and variable scaling to the system set in
    quad double precision.  The mode is either 0, 1, or 2, with
    0 : only scaling of the equations,
    1 : variable scaling without variability reduction,
    2 : variable scaling with variability reduction.
    The number dim is the dimension of the system.
    Returns the scaling coefficients (as complex numbers)
    and the estimated condition number of the scaling problem 
    if the mode > 0, as the last double in the coefficients.
    """
    if vrblvl > 0:
        print('in quad_double_scale_polynomials, mode :', mode)
    phc = get_phcfun(vrblvl-1)
    amode = pointer(c_int32(mode))
    bbb = pointer(c_int32(0))
    lencffs = 16*dim + 8
    cffs = (c_double * lencffs)()
    pcff = pointer(cffs)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_scale_polynomials calls phc', end='')
    retval = phc(592, amode, bbb, pcff, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    allcffs = pcff[0][:lencffs] 
    result = []
    for idx in range(lencffs):
        result.append(allcffs[idx])
    if vrblvl > 0:
        print('the scaling coefficients :', result)
    return result

def double_scale_system(pols, vrblvl=0):
    r"""
    Applies equation and variable scaling in double precision
    to the polynomials in the list *pols*.
    On return is the list of scaled polynomials and the scaling coefficients.
    """
    if vrblvl > 0:
        print('in double_scale_system, pols :')
        for pol in pols:
            print(pol)
    set_double_system(len(pols), pols, vrblvl-1)
    cffs = double_scale_polynomials(2, vrblvl)
    spol = get_double_system(vrblvl-1)
    return (spol, cffs)

def double_double_scale_system(pols, vrblvl=0):
    r"""
    Applies equation and variable scaling in double double precision
    to the polynomials in the list *pols*.
    On return is the list of scaled polynomials and the scaling coefficients.
    """
    if vrblvl > 0:
        print('in double_double_scale_system, pols :')
        for pol in pols:
            print(pol)
    set_double_double_system(len(pols), pols, vrblvl-1)
    cffs = double_double_scale_polynomials(2, vrblvl)
    spol = get_double_double_system(vrblvl-1)
    return (spol, cffs)

def quad_double_scale_system(pols, vrblvl=0):
    r"""
    Applies equation and variable scaling in quad double precision
    to the polynomials in the list *pols*.
    On return is the list of scaled polynomials and the scaling coefficients.
    """
    if vrblvl > 0:
        print('in quad_double_scale_system, pols :')
        for pol in pols:
            print(pol)
    set_quad_double_system(len(pols), pols, vrblvl-1)
    cffs = quad_double_scale_polynomials(2, vrblvl)
    spol = get_quad_double_system(vrblvl-1)
    return (spol, cffs)

def test_double_scale_polynomial(vrblvl=0):
    """
    Runs a test on scaling a polynomial in double precision.
    """
    if vrblvl > 0:
        print('in test_double_scale_polynomial ...')
    orgp = ['100*x^2 + 10*x + 1;']
    print('a badly scaled problem :', orgp)
    print('scaling a polynomial in double precision')
    (scp, cff) = double_scale_system(orgp, vrblvl)
    print('the scaled polynomial', scp)
    print('the scaling coefficients :', cff[0:-2])
    print('estimated inverse condition number :', cff[-2:-1])
    return int(cff[0] != -1.0) # set x = 10^(-1)*y

def test_double_double_scale_polynomial(vrblvl=0):
    """
    Runs a test on scaling a polynomial in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_scale_polynomial ...')
    orgp = ['100*x^2 + 10*x + 1;']
    print('a badly scaled problem :', orgp)
    print('scaling a polynomial in double double precision')
    (scp, cff) = double_double_scale_system(orgp, vrblvl)
    print('the scaled polynomial', scp)
    print('the scaling coefficients :', cff[0:-4])
    print('estimated inverse condition number :', cff[-4:-2])
    return int(cff[0] != -1.0) # set x = 10^(-1)*y

def test_quad_double_scale_polynomial(vrblvl=0):
    """
    Runs a test on scaling a polynomial in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_scale_polynomial ...')
    orgp = ['100*x^2 + 10*x + 1;']
    print('a badly scaled problem :', orgp)
    print('scaling a polynomial in double double precision')
    (scp, cff) = quad_double_scale_system(orgp, vrblvl)
    print('the scaled polynomial', scp)
    print('the scaling coefficients :', cff[0:-8])
    print('estimated inverse condition number :', cff[-8:-4])
    return int(cff[0] != -1.0) # set x = 10^(-1)*y

def main():
    """
    Runs some tests.
    """
    lvl = 1
    fail = test_double_scale_polynomial(lvl)
    fail = fail + test_double_double_scale_polynomial(lvl)
    fail = fail + test_quad_double_scale_polynomial(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()
