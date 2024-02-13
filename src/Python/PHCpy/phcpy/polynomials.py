"""
Exports the definition of polynomial systems, in double, double double,
and quad double precision.  Also polynomials with negative exponents,
the so-called Laurent polynomials, are supported.
"""
from ctypes import c_int32, c_double, pointer, create_string_buffer
from phcpy.version import get_phcfun, int4a2nbr, str2int4a, int4a2str
from phcpy.dimension import set_double_dimension, get_double_dimension
from phcpy.dimension import set_double_double_dimension
from phcpy.dimension import get_double_double_dimension
from phcpy.dimension import set_quad_double_dimension
from phcpy.dimension import get_quad_double_dimension
from phcpy.dimension import set_double_laurent_dimension
from phcpy.dimension import get_double_laurent_dimension
from phcpy.dimension import set_double_double_laurent_dimension
from phcpy.dimension import get_double_double_laurent_dimension
from phcpy.dimension import set_quad_double_laurent_dimension
from phcpy.dimension import get_quad_double_laurent_dimension

def set_double_polynomial(idx, nvr, pol, vrblvl=0):
    """
    Sets the polynomial by the string pol, in double precision,
    with a number of variables no more than nvr, at index idx.
    The index starts at one, not at zero.
    This function does not set the dimension,
    which must be set to a value at least idx.
    """
    if vrblvl > 0:
        print('in set_double_polynomial, pol :', pol, end='')
        print(', idx :', idx, ', nvr :', nvr)
    phc = get_phcfun(vrblvl-1)
    apars = int4a2nbr([len(pol), nvr, idx], vrblvl=vrblvl-1)
    bpol = str2int4a(pol, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vlvl = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_polynomial calls phc', end='')
    retval = phc(76, apars, bpol, ccc, vlvl)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_double_polynomial(idx, nvr, pol, vrblvl=0):
    """
    Sets the polynomial by the string pol, in double double precision,
    with a number of variables no more than nvr, at index idx.
    The index starts at one, not at zero.
    This function does not set the dimension,
    which must be set to a value at least idx.
    """
    if vrblvl > 0:
        print('in set_double_double_polynomial, pol :', pol, end='')
        print(', idx :', idx, ', nvr :', nvr)
    phc = get_phcfun(vrblvl-1)
    apars = int4a2nbr([len(pol), nvr, idx], vrblvl=vrblvl-1)
    bpol = str2int4a(pol, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vlvl = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_double_polynomial calls phc', end='')
    retval = phc(338, apars, bpol, ccc, vlvl)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_quad_double_polynomial(idx, nvr, pol, vrblvl=0):
    """
    Sets the polynomial by the string pol, in quad double precision,
    with a number of variables no more than nvr, at index idx.
    The index starts at one, not at zero.
    This function does not set the dimension,
    which must be set to a value at least idx.
    """
    if vrblvl > 0:
        print('in set_quad_double_polynomial, pol :', pol, end='')
        print(', idx :', idx, ', nvr :', nvr)
    phc = get_phcfun(vrblvl-1)
    apars = int4a2nbr([len(pol), nvr, idx], vrblvl=vrblvl-1)
    bpol = str2int4a(pol, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vlvl = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_quad_double_polynomial calls phc', end='')
    retval = phc(388, apars, bpol, ccc, vlvl)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_laurent_polynomial(idx, nvr, pol, vrblvl=0):
    """
    Sets the polynomial by the string pol, in double precision,
    with a number of variables no more than nvr, at index idx.
    The index starts at one, not at zero.
    This function does not set the dimension,
    which must be set to a value at least idx.
    """
    if vrblvl > 0:
        print('in set_double_laurent_polynomial, pol :', pol, end='')
        print(', idx :', idx, ', nvr :', nvr)
    phc = get_phcfun(vrblvl-1)
    apars = int4a2nbr([len(pol), nvr, idx], vrblvl=vrblvl-1)
    bpol = str2int4a(pol, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vlvl = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_laurent_polynomial calls phc', end='')
    retval = phc(74, apars, bpol, ccc, vlvl)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_double_laurent_polynomial(idx, nvr, pol, vrblvl=0):
    """
    Sets the polynomial by the string pol, in double double precision,
    with a number of variables no more than nvr, at index idx.
    The index starts at one, not at zero.
    This function does not set the dimension,
    which must be set to a value at least idx.
    """
    if vrblvl > 0:
        print('in set_double_double_laurent_polynomial, pol :', pol, end='')
        print(', idx :', idx, ', nvr :', nvr)
    phc = get_phcfun(vrblvl-1)
    apars = int4a2nbr([len(pol), nvr, idx], vrblvl=vrblvl-1)
    bpol = str2int4a(pol, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vlvl = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_double_laurent_polynomial calls phc', end='')
    retval = phc(558, apars, bpol, ccc, vlvl)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_quad_double_laurent_polynomial(idx, nvr, pol, vrblvl=0):
    """
    Sets the polynomial by the string pol, in quad double precision,
    with a number of variables no more than nvr, at index idx.
    The index starts at one, not at zero.
    This function does not set the dimension,
    which must be set to a value at least idx.
    """
    if vrblvl > 0:
        print('in set_quad_double_laurent_polynomial, pol :', pol, end='')
        print(', idx :', idx, ', nvr :', nvr)
    phc = get_phcfun(vrblvl-1)
    apars = int4a2nbr([len(pol), nvr, idx], vrblvl=vrblvl-1)
    bpol = str2int4a(pol, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vlvl = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_quad_double_laurent_polynomial calls phc', end='')
    retval = phc(568, apars, bpol, ccc, vlvl)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_system(nvr, pols, vrblvl=0):
    """
    Sets the system defines by the strings in pols, in double precision,
    with a number of variables no more than nvr.
    The dimension of the system is set to len(pols).
    Returns the sum of the return values of set_double_polynomial.
    """
    if vrblvl > 0:
        print('in set_double_system, nvr :', nvr)
    dim = len(pols)
    set_double_dimension(dim, vrblvl)
    result = 0
    for k in range(dim):
        if vrblvl > 0:
            print(f"-> set_double_system, pols[{k}] = {pols[k]}")
        retval = set_double_polynomial(k+1, nvr, pols[k], vrblvl)
        if vrblvl > 0:
            print('-> set_double_polynomial returns', retval)
        result = result + retval
    return result

def set_double_double_system(nvr, pols, vrblvl=0):
    """
    Sets the system defines by the strings in pols, in double double precision,
    with a number of variables no more than nvr.
    The dimension of the system is set to len(pols).
    Returns the sum of the return values of set_double_double_polynomial.
    """
    if vrblvl > 0:
        print('in set_double_double_system, nvr :', nvr)
    dim = len(pols)
    set_double_double_dimension(dim, vrblvl)
    result = 0
    for k in range(dim):
        if vrblvl > 0:
            print(f"-> set_double_double_system, pols[{k}] = {pols[k]}")
        retval = set_double_double_polynomial(k+1, nvr, pols[k], vrblvl)
        if vrblvl > 0:
            print('-> set_double_double_polynomial returns', retval)
        result = result + retval
    return result

def set_quad_double_system(nvr, pols, vrblvl=0):
    """
    Sets the system defines by the strings in pols, in quad double precision,
    with a number of variables no more than nvr.
    The dimension of the system is set to len(pols).
    Returns the sum of the return values of set_quad_double_polynomial.
    """
    if vrblvl > 0:
        print('in set_quad_double_system, nvr :', nvr)
    dim = len(pols)
    set_quad_double_dimension(dim, vrblvl)
    result = 0
    for k in range(dim):
        if vrblvl > 0:
            print(f"-> set_quad_double_system, pols[{k}] = {pols[k]}")
        retval = set_quad_double_polynomial(k+1, nvr, pols[k], vrblvl)
        if vrblvl > 0:
            print('-> set_quad_double_polynomial returns', retval)
        result = result + retval
    return result

def set_double_laurent_system(nvr, pols, vrblvl=0):
    """
    Sets the laurent system defines by the strings in pols,
    in double precision,
    with a number of variables no more than nvr.
    The dimension of the system is set to len(pols).
    Returns the sum of the return values of set_double_laurent_polynomial.
    """
    if vrblvl > 0:
        print('in set_double_laurent_system, nvr :', nvr)
    dim = len(pols)
    set_double_laurent_dimension(dim, vrblvl)
    result = 0
    for k in range(dim):
        if vrblvl > 0:
            print(f"-> set_double_laurent_system, pols[{k}] = {pols[k]}")
        retval = set_double_laurent_polynomial(k+1, nvr, pols[k], vrblvl)
        if vrblvl > 0:
            print('-> set_double_laurent_polynomial returns', retval)
        result = result + retval
    return result

def set_double_double_laurent_system(nvr, pols, vrblvl=0):
    """
    Sets the laurent system defines by the strings in pols,
    in double double precision,
    with a number of variables no more than nvr.
    The dimension of the system is set to len(pols).
    Returns the sum of the return values 
    of set_double_double_laurent_polynomial.
    """
    if vrblvl > 0:
        print('in set_double_double_laurent_system, nvr :', nvr)
    dim = len(pols)
    set_double_double_laurent_dimension(dim, vrblvl)
    result = 0
    for k in range(dim):
        if vrblvl > 0:
            print(f"-> set_double_double_laurent_system, pols[{k}] = {pols[k]}")
        retval = set_double_double_laurent_polynomial(k+1, nvr, pols[k], vrblvl)
        if vrblvl > 0:
            print('-> set_double_double_laurent_polynomial returns', retval)
        result = result + retval
    return result

def set_quad_double_laurent_system(nvr, pols, vrblvl=0):
    """
    Sets the laurent system defines by the strings in pols,
    in quad double precision,
    with a number of variables no more than nvr.
    The dimension of the system is set to len(pols).
    Returns the sum of the return values of set_quad_double_polynomial.
    """
    if vrblvl > 0:
        print('in set_quad_double_laurent_system, nvr :', nvr)
    dim = len(pols)
    set_quad_double_laurent_dimension(dim, vrblvl)
    result = 0
    for k in range(dim):
        if vrblvl > 0:
            print(f"-> set_quad_double_laurent_system, pols[{k}] = {pols[k]}")
        retval = set_quad_double_laurent_polynomial(k+1, nvr, pols[k], vrblvl)
        if vrblvl > 0:
            print('-> set_quad_double_laurent_polynomial returns', retval)
        result = result + retval
    return result

def get_double_polynomial(idx, vrblvl=0):
    """
    Returns the string representation of the polynomial
    in double precision, at index idx.
    """
    if vrblvl > 0:
        print('in get_double_polynomial idx :', idx)
    phc = get_phcfun(vrblvl-1)
    adx = pointer(c_int32(idx))
    bsz = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_double_polynomial calls phc', end='')
    retval = phc(600, adx, bsz, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('-> size of the polynomial :', bsz[0])
    szd = bsz[0]
    poldata = create_string_buffer(b"", 4*szd)
    if vrblvl > 0:
        print('-> get_double_polynomial calls phc', end='')
    retval = phc(67, adx, poldata, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    strpol = int4a2str(poldata, vrblvl=vrblvl-1)
    pols = strpol.split(';')
    result = pols[0] + ';'
    return result

def get_double_double_polynomial(idx, vrblvl=0):
    """
    Returns the string representation of the polynomial
    in double double precision, at index idx.
    """
    if vrblvl > 0:
        print('in get_double_double_polynomial idx :', idx)
    phc = get_phcfun(vrblvl-1)
    adx = pointer(c_int32(idx))
    bsz = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_double_double_polynomial calls phc', end='')
    retval = phc(601, adx, bsz, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('-> size of the polynomial :', bsz[0])
    szd = bsz[0]
    poldata = create_string_buffer(b"", 4*szd)
    if vrblvl > 0:
        print('-> get_double_double_polynomial calls phc', end='')
    retval = phc(106, adx, poldata, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    strpol = int4a2str(poldata, vrblvl=vrblvl-1)
    pols = strpol.split(';')
    result = pols[0] + ';'
    return result

def get_quad_double_polynomial(idx, vrblvl=0):
    """
    Returns the string representation of the polynomial
    in quad double precision, at index idx.
    """
    if vrblvl > 0:
        print('in get_quad_double_polynomial idx :', idx)
    phc = get_phcfun(vrblvl-1)
    adx = pointer(c_int32(idx))
    bsz = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_quad_double_polynomial calls phc', end='')
    retval = phc(602, adx, bsz, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('-> size of the polynomial :', bsz[0])
    szd = bsz[0]
    poldata = create_string_buffer(b"", 4*szd)
    if vrblvl > 0:
        print('-> get_quad_double_polynomial calls phc', end='')
    retval = phc(107, adx, poldata, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    strpol = int4a2str(poldata, vrblvl=vrblvl-1)
    pols = strpol.split(';')
    result = pols[0] + ';'
    return result

def get_double_laurent_polynomial(idx, vrblvl=0):
    """
    Returns the string representation of the Laurent polynomial
    in double precision, at index idx.
    """
    if vrblvl > 0:
        print('in get_double_laurent_polynomial idx :', idx)
    phc = get_phcfun(vrblvl-1)
    adx = pointer(c_int32(idx))
    bsz = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_double_laurent_polynomial calls phc', end='')
    retval = phc(604, adx, bsz, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('-> size of the polynomial :', bsz[0])
    szd = bsz[0]
    poldata = create_string_buffer(b"", 4*szd)
    if vrblvl > 0:
        print('-> get_double_laurent_polynomial calls phc', end='')
    retval = phc(128, adx, poldata, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    strpol = int4a2str(poldata, vrblvl=vrblvl-1)
    pols = strpol.split(';')
    result = pols[0] + ';'
    return result

def get_double_double_laurent_polynomial(idx, vrblvl=0):
    """
    Returns the string representation of the Laurent polynomial
    in double double precision, at index idx.
    """
    if vrblvl > 0:
        print('in get_double_double_laurent_polynomial idx :', idx)
    phc = get_phcfun(vrblvl-1)
    adx = pointer(c_int32(idx))
    bsz = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_double_double_laurent_polynomial calls phc', end='')
    retval = phc(605, adx, bsz, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('-> size of the polynomial :', bsz[0])
    szd = bsz[0]
    poldata = create_string_buffer(b"", 4*szd)
    if vrblvl > 0:
        print('-> get_double_double_laurent_polynomial calls phc', end='')
    retval = phc(559, adx, poldata, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    strpol = int4a2str(poldata, vrblvl=vrblvl-1)
    pols = strpol.split(';')
    result = pols[0] + ';'
    return result

def get_quad_double_laurent_polynomial(idx, vrblvl=0):
    """
    Returns the string representation of the Laurent polynomial
    in quad double precision, at index idx.
    """
    if vrblvl > 0:
        print('in get_quad_double_laurent_polynomial idx :', idx)
    phc = get_phcfun(vrblvl-1)
    adx = pointer(c_int32(idx))
    bsz = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_quad_double_laurent_polynomial calls phc', end='')
    retval = phc(606, adx, bsz, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('-> size of the polynomial :', bsz[0])
    szd = bsz[0]
    poldata = create_string_buffer(b"", 4*szd)
    if vrblvl > 0:
        print('-> get_quad_double_laurent_polynomial calls phc', end='')
    retval = phc(569, adx, poldata, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    strpol = int4a2str(poldata, vrblvl=vrblvl-1)
    pols = strpol.split(';')
    result = pols[0] + ';'
    return result

def get_double_system(vrblvl=0):
    """
    Returns the string representation of the polynomials
    in double precision.
    """
    dim = get_double_dimension(vrblvl)
    if vrblvl > 0:
        print('in get_double_system, dim :', dim)
    result = []
    for k in range(1, dim+1):
        pol = get_double_polynomial(k, vrblvl)
        if vrblvl > 0:
            print('get_double_system retrieved pol :', pol)
        result.append(pol)
    return result

def get_double_double_system(vrblvl=0):
    """
    Returns the string representation of the polynomials
    in double double precision.
    """
    dim = get_double_double_dimension(vrblvl)
    if vrblvl > 0:
        print('in get_double_double_system, dim :', dim)
    result = []
    for k in range(1, dim+1):
        pol = get_double_double_polynomial(k, vrblvl)
        if vrblvl > 0:
            print('-> get_double_double_system retrieved pol :', pol)
        result.append(pol)
    return result

def get_quad_double_system(vrblvl=0):
    """
    Returns the string representation of the polynomials
    in quad double precision.
    """
    dim = get_quad_double_dimension(vrblvl)
    if vrblvl > 0:
        print('in get_quad_double_system, dim :', dim)
    result = []
    for k in range(1, dim+1):
        pol = get_quad_double_polynomial(k, vrblvl)
        if vrblvl > 0:
            print('get_quad_double_system retrieved pol :', pol)
        result.append(pol)
    return result

def get_double_laurent_system(vrblvl=0):
    """
    Returns the string representation of the Laurent polynomials
    in double precision.
    """
    dim = get_double_laurent_dimension(vrblvl)
    if vrblvl > 0:
        print('in get_double_laurent_system, dim :', dim)
    result = []
    for k in range(1, dim+1):
        pol = get_double_laurent_polynomial(k, vrblvl)
        if vrblvl > 0:
            print('get_double_laurent_system retrieved pol :', pol)
        result.append(pol)
    return result

def get_double_double_laurent_system(vrblvl=0):
    """
    Returns the string representation of the Laurent polynomials
    in double double precision.
    """
    dim = get_double_double_laurent_dimension(vrblvl)
    if vrblvl > 0:
        print('in get_double_double_laurent_system, dim :', dim)
    result = []
    for k in range(1, dim+1):
        pol = get_double_double_laurent_polynomial(k, vrblvl)
        if vrblvl > 0:
            print('get_double_double_laurent_system retrieved pol :', pol)
        result.append(pol)
    return result

def get_quad_double_laurent_system(vrblvl=0):
    """
    Returns the string representation of the Laurent polynomials
    in quad double precision.
    """
    dim = get_quad_double_laurent_dimension(vrblvl)
    if vrblvl > 0:
        print('in get_quad_double_laurent_system, dim :', dim)
    result = []
    for k in range(1, dim+1):
        pol = get_quad_double_laurent_polynomial(k, vrblvl)
        if vrblvl > 0:
            print('get_quad_double_laurent_system retrieved pol :', pol)
        result.append(pol)
    return result

def get_double_number_terms(equ, vrblvl=0):
    """
    Returns the number of terms in the equation with index equ,
    set as a polynomial system in double precision.
    """
    if vrblvl > 0:
        print('in get_double_number_terms, equ :', equ)
    phc = get_phcfun(vrblvl-1)
    idx = (c_int32 * 2)()
    idx[0] = equ
    idx[1] = equ
    aidx = pointer(idx)
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_double_number_of_terms calls phc', end='')
    retval = phc(24, aidx, bbb, ccc, vrb)
    result = aidx[0][0]
    if vrblvl > 0:
        print(', return value :', retval)
        print('number of terms :', result)
    return result

def get_double_number_laurent_terms(equ, vrblvl=0):
    """
    Returns the number of terms in the equation with index equ,
    set as a Laurent system in double precision.
    """
    if vrblvl > 0:
        print('in get_double_number_laurent_terms, equ :', equ)
    phc = get_phcfun(vrblvl-1)
    idx = (c_int32 * 2)()
    idx[0] = equ
    idx[1] = equ
    aidx = pointer(idx)
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_double_number_of_laurent_terms calls phc', end='')
    retval = phc(124, aidx, bbb, ccc, vrb)
    result = aidx[0][0]
    if vrblvl > 0:
        print(', return value :', retval)
        print('number of terms :', result)
    return result

def get_double_double_number_terms(equ, vrblvl=0):
    """
    Returns the number of terms in the equation with index equ,
    set as a polynomial system in double double precision.
    """
    if vrblvl > 0:
        print('in get_double_double_number_terms, equ :', equ)
    phc = get_phcfun(vrblvl-1)
    idx = (c_int32 * 2)()
    idx[0] = equ
    idx[1] = equ
    aidx = pointer(idx)
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_double_double_number_of_terms calls phc', end='')
    retval = phc(334, aidx, bbb, ccc, vrb)
    result = aidx[0][0]
    if vrblvl > 0:
        print(', return value :', retval)
        print('number of terms :', result)
    return result

def get_double_double_number_laurent_terms(equ, vrblvl=0):
    """
    Returns the number of terms in the equation with index equ,
    set as a Laurent system in double double precision.
    """
    if vrblvl > 0:
        print('in get_double_double_number_laurent_terms, equ :', equ)
    phc = get_phcfun(vrblvl-1)
    idx = (c_int32 * 2)()
    idx[0] = equ
    idx[1] = equ
    aidx = pointer(idx)
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_double_double_number_of_laurent_terms calls phc', end='')
    retval = phc(554, aidx, bbb, ccc, vrb)
    result = aidx[0][0]
    if vrblvl > 0:
        print(', return value :', retval)
        print('number of terms :', result)
    return result

def get_quad_double_number_terms(equ, vrblvl=0):
    """
    Returns the number of terms in the equation with index equ,
    set as a polynomial system in quad double precision.
    """
    if vrblvl > 0:
        print('in get_quad_double_number_terms, equ :', equ)
    phc = get_phcfun(vrblvl-1)
    idx = (c_int32 * 2)()
    idx[0] = equ
    idx[1] = equ
    aidx = pointer(idx)
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_quad_double_number_of_terms calls phc', end='')
    retval = phc(384, aidx, bbb, ccc, vrb)
    result = aidx[0][0]
    if vrblvl > 0:
        print(', return value :', retval)
        print('number of terms :', result)
    return result

def get_quad_double_number_laurent_terms(equ, vrblvl=0):
    """
    Returns the number of terms in the equation with index equ,
    set as a Laurent system in quad double precision.
    """
    if vrblvl > 0:
        print('in get_quad_double_number_laurent_terms, equ :', equ)
    phc = get_phcfun(vrblvl-1)
    idx = (c_int32 * 2)()
    idx[0] = equ
    idx[1] = equ
    aidx = pointer(idx)
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_quad_double_number_of_laurent_terms calls phc', end='')
    retval = phc(564, aidx, bbb, ccc, vrb)
    result = aidx[0][0]
    if vrblvl > 0:
        print(', return value :', retval)
        print('number of terms :', result)
    return result

def string_of_symbols(maxlen=100, vrblvl=0):
    """
    Returns the list of all symbols (as strings),
    defined when storing polynomials.
    The maxlen on entry equals the maximum number of characters
    in the symbol string, that is: the sequence of all string
    representations of the symbols, separated by one space.
    """
    if vrblvl > 0:
        print('in string_of_symbols, maxlen :', maxlen)
    phc = get_phcfun(vrblvl-1)
    slen = pointer(c_int32(0))
    ssym = create_string_buffer(b"", maxlen*4)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> string_of_symbols calls phc', end='')
    retval = phc(295, slen, ssym, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('number of characters :', slen[0])
    variables = int4a2str(ssym, vrblvl=vrblvl-1)
    result = variables[0:slen[0]].split(' ')
    return result

def number_of_symbols(pols, vrblvl=0):
    """
    Returns the number of symbols that appear as variables
    in the polynomials, givein in the list of strings pols.
    Useful to determine whether a system is square or not.
    """
    if vrblvl > 0:
        print('in number_of_symbols, pols :')
        for pol in pols:
            print(pol)
    inpols = ''.join(pols)
    phc = get_phcfun(vrblvl-1)
    slen = pointer(c_int32(len(inpols)))
    bpol = str2int4a(inpols)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> number_of_symbols calls phc', end='')
    retval = phc(439, slen, bpol, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('-> number_of_symbols, result :', slen[0])
    return slen[0]

def check_semicolons(pols, vrblvl=0):
    """
    Counts the number of semicolons and writes a warning
    if the number of semicolons in pols does not match len(pols).
    This warning may prevent unwanted concatenation, e.g., as in
    pols = ['x^2 + 4*y^2 - 4;' '2*y^2 - x;']
    where the omitted comma results in one polynomial in two variables.
    """
    if vrblvl > 0:
        print('in check_semicolons, pols :')
        for pol in pols:
            print(pol)
    checksum = sum(pol.count(';') for pol in pols)
    if checksum != len(pols):
        print('#semicolons :', checksum, '!=', len(pols), 'polynomials.')
    return int(checksum != len(pols))

def is_square(pols, vrblvl=0):
    r"""
    Given in the list *pols* are string representations 
    of Laurent polynomials.
    A system is square if it has as many unknowns as equations.
    Returns True if the system is square, False otherwise.
    """
    if vrblvl > 0:
        print('in is_square, pols :')
        for pol in pols:
            print(pol)
    nbrvar = number_of_symbols(pols)
    nbreqs = len(pols)
    return nbrvar == nbreqs

def degree_of_double_polynomial(idx, vrblvl=0):
    """
    Returns the degree of the polynomial with double precision
    coefficients stored at position idx.
    """
    if vrblvl > 0:
        print('in degree_of_double_polynomial, idx :', idx)
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(idx))
    bdeg = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> degree_of_double_system calls phc', end='')
    retval = phc(119, aidx, bdeg, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('-> degree of polynomial :', bdeg[0])
    return bdeg[0]

def clear_double_system(vrblvl=0):
    """
    Clears the system set in double precision.
    """
    if vrblvl > 0:
        print('in clear_double_system ...')
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_double_system calls phc', end='')
    retval = phc(27, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_double_double_system(vrblvl=0):
    """
    Clears the system set in double double precision.
    """
    if vrblvl > 0:
        print('in clear_double_double_system ...')
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_double_double_system calls phc', end='')
    retval = phc(337, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_quad_double_system(vrblvl=0):
    """
    Clears the system set in quad double precision.
    """
    if vrblvl > 0:
        print('in clear_quad_double_system ...')
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_quad_double_system calls phc', end='')
    retval = phc(387, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_double_laurent_system(vrblvl=0):
    """
    Clears the Laurent system set in double precision.
    """
    if vrblvl > 0:
        print('in clear_double_laurent_system ...')
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_double_laurent_system calls phc', end='')
    retval = phc(127, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_double_double_laurent_system(vrblvl=0):
    """
    Clears the Laurent system set in double double precision.
    """
    if vrblvl > 0:
        print('in clear_double_double_laurent_system ...')
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_double_double_laurent_system calls phc', end='')
    retval = phc(557, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_quad_double_laurent_system(vrblvl=0):
    """
    Clears the Laurent system set in quad double precision.
    """
    if vrblvl > 0:
        print('in clear_quad_double_laurent_system ...')
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_quad_double_laurent_system calls phc', end='')
    retval = phc(567, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_symbol_table(vrblvl=0):
    """
    Clears the table of symbols used to represent polynomials.
    """
    if vrblvl > 0:
        print('in clear_symbol_table ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_symbol_table calls phc', end='')
    retval = phc(29, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_double_syspool(dim, vrblvl=0):
    """
    Initialize the systems pool in double precision with dim.
    """
    if vrblvl > 0:
        print('in initialize_double_syspool, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_double_syspool calls phc', end='')
    retval = phc(300, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_double_double_syspool(dim, vrblvl=0):
    """
    Initialize the systems pool in double double precision with dim.
    """
    if vrblvl > 0:
        print('in initialize_double_double_syspool, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_double_double_syspool calls phc', end='')
    retval = phc(318, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_quad_double_syspool(dim, vrblvl=0):
    """
    Initialize the systems pool in quad double precision with dim.
    """
    if vrblvl > 0:
        print('in initialize_quad_double_syspool, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_quad_double_syspool calls phc', end='')
    retval = phc(319, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def size_double_syspool(vrblvl=0):
    """
    Returns the size of the systems pool in double precision.
    """
    if vrblvl > 0:
        print('in size_double_syspool ...')
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> size_double_syspool calls phc', end='')
    retval = phc(301, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('number of systems in the pool :', adim[0])
    return adim[0]

def size_double_double_syspool(vrblvl=0):
    """
    Returns the size of the systems pool in double double precision.
    """
    if vrblvl > 0:
        print('in size_double_double_syspool ...')
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> size_double_double_syspool calls phc', end='')
    retval = phc(316, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('number of systems in the pool :', adim[0])
    return adim[0]

def size_quad_double_syspool(vrblvl=0):
    """
    Returns the size of the systems pool in quad double precision.
    """
    if vrblvl > 0:
        print('in size_quad_double_syspool ...')
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> size_quad_double_syspool calls phc', end='')
    retval = phc(317, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('number of systems in the pool :', adim[0])
    return adim[0]

def copy_to_double_syspool(idx, vrblvl=0):
    """
    Copies the system set in double precision to position idx
    in the systems pool.
    """
    if vrblvl > 0:
        print('in copy_to_double_syspool, idx :', idx)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(idx))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_to_double_syspool calls phc', end='')
    retval = phc(304, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_to_double_double_syspool(idx, vrblvl=0):
    """
    Copies the system set in double double precision to position idx
    in the systems pool.
    """
    if vrblvl > 0:
        print('in copy_to_double_double_syspool, idx :', idx)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(idx))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_to_double_double_syspool calls phc', end='')
    retval = phc(608, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_to_quad_double_syspool(idx, vrblvl=0):
    """
    Copies the system set in quad double precision to position idx
    in the systems pool.
    """
    if vrblvl > 0:
        print('in copy_to_quad_double_syspool, idx :', idx)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(idx))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_to_quad_double_syspool calls phc', end='')
    retval = phc(609, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_from_double_syspool(idx, vrblvl=0):
    """
    Copies the system in double precision at position idx
    in the systems pool to the defined system in double precision.
    """
    if vrblvl > 0:
        print('in copy_from_double_syspool, idx :', idx)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(idx))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_from_double_syspool calls phc', end='')
    retval = phc(313, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_from_double_double_syspool(idx, vrblvl=0):
    """
    Copies the system in double double precision at position idx
    in the systems pool to the defined system in double double precision.
    """
    if vrblvl > 0:
        print('in copy_from_double_double_syspool, idx :', idx)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(idx))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_from_double_double_syspool calls phc', end='')
    retval = phc(314, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_from_quad_double_syspool(idx, vrblvl=0):
    """
    Copies the system in quad double precision at position idx
    in the systems pool to the defined system in quad double precision.
    """
    if vrblvl > 0:
        print('in copy_from_quad_double_syspool, idx :', idx)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(idx))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_from_quad_double_syspool calls phc', end='')
    retval = phc(315, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_double_syspool(vrblvl=0):
    """
    Clears the systems pool in double precision.
    """
    if vrblvl > 0:
        print('in clear_double_syspool ...')
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_double_syspool calls phc', end='')
    retval = phc(697, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_double_double_syspool(vrblvl=0):
    """
    Clears the systems pool in double double precision.
    """
    if vrblvl > 0:
        print('in clear_double_double_syspool ...')
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_double_double_syspool calls phc', end='')
    retval = phc(698, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_quad_double_syspool(vrblvl=0):
    """
    Clears the systems pool in quad double precision.
    """
    if vrblvl > 0:
        print('in clear_quad_double_syspool ...')
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_quad_double_syspool calls phc', end='')
    retval = phc(699, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def test_double_polynomial(vrblvl=0):
    """
    Tests the setting and getting of a polynomial, in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_double_polynomial ...')
    set_double_dimension(2, vrblvl)
    dim = get_double_dimension(vrblvl)
    print('the dimension :', dim)
    org = "x*y - 1;"
    idx = 1
    set_double_polynomial(idx, dim, org, vrblvl)
    pol = get_double_polynomial(idx, vrblvl)
    print('the retrieved polynomial :', pol)
    smb = string_of_symbols(100, vrblvl)
    print('the list of symbols :', smb)
    return int(len(smb) != 2)

def test_double_double_polynomial(vrblvl=0):
    """
    Tests the setting and getting of a polynomial,
    in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_double_double_polynomial ...')
    set_double_double_dimension(2, vrblvl)
    dim = get_double_double_dimension(vrblvl)
    print('the dimension :', dim)
    org = "x*y - 1;"
    idx = 1
    set_double_double_polynomial(idx, dim, org, vrblvl)
    pol = get_double_double_polynomial(idx, vrblvl)
    print('the retrieved polynomial :', pol)
    smb = string_of_symbols(100, vrblvl)
    print('the list of symbols :', smb)
    return int(len(smb) != 2)

def test_quad_double_polynomial(vrblvl=0):
    """
    Tests the setting and getting of a polynomial,
    in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_quad_double_polynomial ...')
    set_quad_double_dimension(2, vrblvl)
    dim = get_quad_double_dimension(vrblvl)
    print('the dimension :', dim)
    org = "x*y - 1;"
    idx = 1
    set_quad_double_polynomial(idx, dim, org, vrblvl)
    pol = get_quad_double_polynomial(idx, vrblvl)
    print('the retrieved polynomial :', pol)
    smb = string_of_symbols(100, vrblvl)
    print('the list of symbols :', smb)
    return int(len(smb) != 2)

def test_double_laurent_polynomial(vrblvl=0):
    """
    Tests the setting and getting of a Laurent polynomial,
    in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_double_laurent_polynomial ...')
    set_double_laurent_dimension(2, vrblvl)
    dim = get_double_laurent_dimension(vrblvl)
    print('the dimension :', dim)
    org = "x*y^(-3) - 1;"
    idx = 1
    set_double_laurent_polynomial(idx, dim, org, vrblvl)
    pol = get_double_laurent_polynomial(idx, vrblvl)
    print('the retrieved polynomial :', pol)
    smb = string_of_symbols(100, vrblvl)
    print('the list of symbols :', smb)
    return int(len(smb) != 2)

def test_double_double_laurent_polynomial(vrblvl=0):
    """
    Tests the setting and getting of a Laurent polynomial,
    in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_double_double_laurent_polynomial ...')
    set_double_double_laurent_dimension(2, vrblvl)
    dim = get_double_double_laurent_dimension(vrblvl)
    print('the dimension :', dim)
    org = "x*y^(-3) - 1;"
    idx = 1
    set_double_double_laurent_polynomial(idx, dim, org, vrblvl)
    pol = get_double_double_laurent_polynomial(idx, vrblvl)
    print('the retrieved polynomial :', pol)
    smb = string_of_symbols(100, vrblvl)
    print('the list of symbols :', smb)
    return int(len(smb) != 2)

def test_quad_double_laurent_polynomial(vrblvl=0):
    """
    Tests the setting and getting of a Laurent polynomial,
    in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_quad_double_laurent_polynomial ...')
    set_quad_double_laurent_dimension(2, vrblvl)
    dim = get_quad_double_laurent_dimension(vrblvl)
    print('the dimension :', dim)
    org = "x*y^(-3) - 1;"
    idx = 1
    set_quad_double_laurent_polynomial(idx, dim, org, vrblvl)
    pol = get_quad_double_laurent_polynomial(idx, vrblvl)
    print('the retrieved polynomial :', pol)
    smb = string_of_symbols(100, vrblvl)
    print('the list of symbols :', smb)
    return int(len(smb) != 2)

def test_double_system(vrblvl=0):
    """
    Tests the setting and getting of a system,
    in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_double_system ...')
    polynomials = ["x^3 + 2*x*y - 1;", "x + y - 1;", "x - 1;"]
    dim = number_of_symbols(polynomials, vrblvl)
    print('number of symbols :', dim)
    if dim == len(polynomials):
        print('The system is square.')
    else:
        print('number of polynomials :', len(polynomials))
        print('  number of variables :', dim)
        print('The system is not square.')
    set_double_system(dim, polynomials, vrblvl)
    pols = get_double_system(vrblvl)
    print('the retrieved polynomials :')
    for pol in pols:
        print(pol)
    return int(len(pols) != 3)

def test_double_double_system(vrblvl=0):
    """
    Tests the setting and getting of a system,
    in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_double_double_system ...')
    polynomials = ["x^3 + 2*x*y - 1;", "x + y - 1/3;", "x - 1;"]
    dim = number_of_symbols(polynomials, vrblvl)
    print('number of symbols :', dim)
    if dim == len(polynomials):
        print('The system is square.')
    else:
        print('number of polynomials :', len(polynomials))
        print('  number of variables :', dim)
        print('The system is not square.')
    set_double_double_system(dim, polynomials, vrblvl)
    pols = get_double_double_system(vrblvl)
    print('the retrieved polynomials :')
    for pol in pols:
        print(pol)
    return int(len(pols) != 3)

def test_quad_double_system(vrblvl=0):
    """
    Tests the setting and getting of a system,
    in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_quad_double_system ...')
    polynomials = ["x^3 + 2*x*y - 1;", "x + y - 1/3;", "y - 1;"]
    dim = number_of_symbols(polynomials, vrblvl)
    print('number of symbols :', dim)
    if dim == len(polynomials):
        print('The system is square.')
    else:
        print('number of polynomials :', len(polynomials))
        print('  number of variables :', dim)
        print('The system is not square.')
    set_quad_double_system(dim, polynomials, vrblvl)
    pols = get_quad_double_system(vrblvl)
    print('the retrieved polynomials :')
    for pol in pols:
        print(pol)
    return int(len(pols) != 3)

def test_double_laurent_system(vrblvl=0):
    """
    Tests the setting and getting of a laurent system,
    in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_double_laurent_system ...')
    polynomials = ["x^(-3) + 2*x*y - 1;", "x + y^(-1) - 1;", "x - 1;"]
    dim = number_of_symbols(polynomials, vrblvl)
    print('number of symbols :', dim)
    if dim == len(polynomials):
        print('The system is square.')
    else:
        print('number of polynomials :', len(polynomials))
        print('  number of variables :', dim)
        print('The system is not square.')
    set_double_laurent_system(dim, polynomials, vrblvl)
    pols = get_double_laurent_system(vrblvl)
    print('the retrieved polynomials :')
    for pol in pols:
        print(pol)
    return int(len(pols) != 3)

def test_double_double_laurent_system(vrblvl=0):
    """
    Tests the setting and getting of a laurent system,
    in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('test_double_double_laurent_system ...')
    polynomials = ["x^(-3) + 2*x*y - 1;", "x + y**(-1) - 1/3;", "x - 1;"]
    dim = number_of_symbols(polynomials, vrblvl)
    print('number of symbols :', dim)
    if dim == len(polynomials):
        print('The system is square.')
    else:
        print('number of polynomials :', len(polynomials))
        print('  number of variables :', dim)
        print('The system is not square.')
    set_double_double_laurent_system(dim, polynomials, vrblvl)
    pols = get_double_double_laurent_system(vrblvl)
    print('the retrieved polynomials :')
    for pol in pols:
        print(pol)
    return len(pols) != 3

def test_quad_double_laurent_system(vrblvl=0):
    """
    Tests the setting and getting of a laurent system,
    in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_quad_double_laurent_system ...')
    polynomials = ["x^(-3) + 2*x*y - 1;", "x + y**(-1) - 1/3;", "y - 1;"]
    dim = number_of_symbols(polynomials, vrblvl)
    print('number of symbols :', dim)
    if dim == len(polynomials):
        print('The system is square.')
    else:
        print('number of polynomials :', len(polynomials))
        print('  number of variables :', dim)
        print('The system is not square.')
    set_quad_double_laurent_system(dim, polynomials, vrblvl)
    pols = get_quad_double_laurent_system(vrblvl)
    print('the retrieved polynomials :')
    for pol in pols:
        print(pol)
    return int(len(pols) != 3)

def test_double_syspool(vrblvl=0):
    """
    Tests the systems pool in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_double_syspool ...')
    initialize_double_syspool(3, vrblvl)
    dim = size_double_syspool(vrblvl)
    print('The size of the systems pool :', dim)
    pol1 = ['t - 1;']
    set_double_system(1, pol1, vrblvl)
    copy_to_double_syspool(1)
    pol2 = ['t - 2;']
    set_double_system(1, pol2, vrblvl)
    copy_to_double_syspool(2)
    pol3 = ['t - 3;']
    set_double_system(1, pol3, vrblvl)
    copy_to_double_syspool(3)
    for i in range(1, dim+1):
        clear_double_system(vrblvl)
        copy_from_double_syspool(i)
        pols = get_double_system(vrblvl)
        print('system at', i, 'in the pool :', pols)
    clear_double_syspool(vrblvl)
    return int(dim != 3)

def test_double_double_syspool(vrblvl=0):
    """
    Tests the systems pool in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_double_double_syspool ...')
    initialize_double_double_syspool(3, vrblvl)
    dim = size_double_double_syspool(vrblvl)
    print('The size of the systems pool :', dim)
    pol1 = ['t - 1/3;']
    set_double_double_system(1, pol1, vrblvl)
    copy_to_double_double_syspool(1)
    pol2 = ['t - 2/3;']
    set_double_double_system(1, pol2, vrblvl)
    copy_to_double_double_syspool(2)
    pol3 = ['t - 1;']
    set_double_double_system(1, pol3, vrblvl)
    copy_to_double_double_syspool(3)
    for i in range(1, dim+1):
        clear_double_double_system(vrblvl)
        copy_from_double_double_syspool(i)
        pols = get_double_double_system(vrblvl)
        print('system at', i, 'in the pool :', pols)
    clear_double_double_syspool(vrblvl)
    return int(dim != 3)

def test_quad_double_syspool(vrblvl=0):
    """
    Tests the systems pool in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_quad_double_syspool ...')
    initialize_quad_double_syspool(3, vrblvl)
    dim = size_quad_double_syspool(vrblvl)
    print('The size of the systems pool :', dim)
    pol1 = ['t - 1/3;']
    set_quad_double_system(1, pol1, vrblvl)
    copy_to_quad_double_syspool(1)
    pol2 = ['t - 2/3;']
    set_quad_double_system(1, pol2, vrblvl)
    copy_to_quad_double_syspool(2)
    pol3 = ['t - 1;']
    set_quad_double_system(1, pol3, vrblvl)
    copy_to_quad_double_syspool(3)
    for i in range(1, dim+1):
        clear_quad_double_system(vrblvl)
        copy_from_quad_double_syspool(i)
        pols = get_quad_double_system(vrblvl)
        print('system at', i, 'in the pool :', pols)
    clear_quad_double_syspool(vrblvl)
    return int(dim != 3)

def test_degree_of_double_polynomial(vrblvl=0):
    """
    Tests the degree of a polynomial in double precision.
    """
    if vrblvl > 0:
        print('in test_degree_of_double_polynomial ...')
    pols = ['x^2*y + y^2 + 1;', 'x^2 + y^2 - 1;']
    if vrblvl > 0:
        print('in test_degree_of_double_polynomial ...')
    set_double_system(len(pols), pols, vrblvl)
    degsum = 0
    for (idx, pol) in enumerate(pols):
        print('index :', idx)
        deg = degree_of_double_polynomial(idx+1, vrblvl)
        print('degree of', pol, ':', deg)
        degsum = degsum + deg
    return int(degsum != 5)

def test_double_number_terms(vrblvl=0):
    """
    Tests the number of terms in a polynomial system
    set in double precision.
    """
    if vrblvl > 0:
        print('in test_double_number_terms ...')
    pols = ['x^2*y;', 'x^2 + y^2;', 'x + z + y;']
    set_double_system(len(pols), pols, vrblvl)
    fail = 0
    for idx in range(len(pols)):
        nbterms = get_double_number_terms(idx+1, vrblvl)
        if vrblvl > 0:
            print('number of terms in equation', idx+1, ':', nbterms)
        fail = fail + int(nbterms != idx+1)
    return fail

def test_double_number_laurent_terms(vrblvl=0):
    """
    Tests the number of terms in a Laurent system
    set in double precision.
    """
    if vrblvl > 0:
        print('in test_double_number_laurent_terms ...')
    pols = ['x^2*y;', 'x^2 + y^2;', 'x + z + y;']
    set_double_laurent_system(len(pols), pols, vrblvl)
    fail = 0
    for idx in range(len(pols)):
        nbterms = get_double_number_laurent_terms(idx+1, vrblvl)
        if vrblvl > 0:
            print('number of terms in equation', idx+1, ':', nbterms)
        fail = fail + int(nbterms != idx+1)
    return fail

def test_double_double_number_terms(vrblvl=0):
    """
    Tests the number of terms in a polynomial system
    set in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_number_terms ...')
    pols = ['x^2*y;', 'x^2 + y^2;', 'x + z + y;']
    set_double_double_system(len(pols), pols, vrblvl)
    fail = 0
    for idx in range(len(pols)):
        nbterms = get_double_double_number_terms(idx+1, vrblvl)
        if vrblvl > 0:
            print('number of terms in equation', idx+1, ':', nbterms)
        fail = fail + int(nbterms != idx+1)
    return fail

def test_double_double_number_laurent_terms(vrblvl=0):
    """
    Tests the number of terms in a Laurent system
    set in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_number_laurent_terms ...')
    pols = ['x^2*y;', 'x^2 + y^2;', 'x + z + y;']
    set_double_double_laurent_system(len(pols), pols, vrblvl)
    fail = 0
    for idx in range(len(pols)):
        nbterms = get_double_double_number_laurent_terms(idx+1, vrblvl)
        if vrblvl > 0:
            print('number of terms in equation', idx+1, ':', nbterms)
        fail = fail + int(nbterms != idx+1)
    return fail

def test_quad_double_number_terms(vrblvl=0):
    """
    Tests the number of terms in a polynomial system
    set in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_number_terms ...')
    pols = ['x^2*y;', 'x^2 + y^2;', 'x + z + y;']
    set_quad_double_system(len(pols), pols, vrblvl)
    fail = 0
    for idx in range(len(pols)):
        nbterms = get_quad_double_number_terms(idx+1, vrblvl)
        if vrblvl > 0:
            print('number of terms in equation', idx+1, ':', nbterms)
        fail = fail + int(nbterms != idx+1)
    return fail

def test_quad_double_number_laurent_terms(vrblvl=0):
    """
    Tests the number of terms in a Laurent system
    set in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_number_laurent_terms ...')
    pols = ['x^2*y;', 'x^2 + y^2;', 'x + z + y;']
    set_quad_double_laurent_system(len(pols), pols, vrblvl)
    fail = 0
    for idx in range(len(pols)):
        nbterms = get_quad_double_number_laurent_terms(idx+1, vrblvl)
        if vrblvl > 0:
            print('number of terms in equation', idx+1, ':', nbterms)
        fail = fail + int(nbterms != idx+1)
    return fail

def test_is_square(vrblvl=0):
    """
    Tests on catching an omitted comma between the polynomials.
    """
    pols = ['x^2 + 4*y^2 - 4;' '2*y^2 - x;']
    retval = check_semicolons(pols, vrblvl)
    fail = int(retval != 1)
    retval = is_square(pols, vrblvl)
    if retval:
        print('The system is square.')
    else:
        print('The system is not square.')
    fail = fail + int(retval)
    return fail

def main():
    """
    Runs tests on the set/get polynomials.
    """
    lvl = 1
    fail = test_double_polynomial(lvl)
    fail = fail + test_double_double_polynomial(lvl)
    fail = fail + test_quad_double_polynomial(lvl)
    fail = fail + test_double_laurent_polynomial(lvl)
    fail = fail + test_double_double_laurent_polynomial(lvl)
    fail = fail + test_quad_double_laurent_polynomial(lvl)
    fail = fail + test_double_system(lvl)
    fail = fail + test_double_double_system(lvl)
    fail = fail + test_quad_double_system(lvl)
    fail = fail + test_double_laurent_system(lvl)
    fail = fail + test_double_double_laurent_system(lvl)
    fail = fail + test_quad_double_laurent_system(lvl)
    fail = fail + test_double_syspool(lvl)
    fail = fail + test_double_double_syspool(lvl)
    fail = fail + test_quad_double_syspool(lvl)
    fail = fail + test_degree_of_double_polynomial(lvl)
    fail = fail + test_double_number_terms(lvl)
    fail = fail + test_double_number_laurent_terms(lvl)
    fail = fail + test_double_double_number_terms(lvl)
    fail = fail + test_double_double_number_laurent_terms(lvl)
    fail = fail + test_quad_double_number_terms(lvl)
    fail = fail + test_quad_double_number_laurent_terms(lvl)
    fail = fail + test_is_square(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=="__main__":
    main()
