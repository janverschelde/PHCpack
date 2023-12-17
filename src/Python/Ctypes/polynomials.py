"""
Exports the definition of polynomial systems.
"""
from ctypes import c_int, c_double, pointer, create_string_buffer
from version import get_phcfun, int4a2nbr, str2int4a, int4a2str
from dimension import set_double_dimension, get_double_dimension
from dimension import set_double_double_dimension, get_double_double_dimension
from dimension import set_quad_double_dimension, get_quad_double_dimension

def set_double_polynomial(idx, nvr, pol, vrblvl=0):
    """
    Sets the polynomial by the string pol, in double precision,
    with a number of variables no more than nvr, at index idx.
    The index starts at one, not at zero.
    This function does not set the dimension,
    which must be set to a value at least idx.
    """
    if vrblvl > 0:
        print("-> set_double_polynomial, pol = ", pol)
    phc = get_phcfun()
    vrb = (vrblvl > 0)
    apars = int4a2nbr([len(pol), nvr, idx], vrb)
    bpol = str2int4a(pol)
    ccc = pointer(c_double(0.0))
    vlvl = c_int(vrblvl)
    retval = phc(76, apars, bpol, ccc, vlvl)
    if vrblvl > 0:
        print('-> set_double_polynomial return value :', retval)
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
        print("-> set_double_double_polynomial, pol = ", pol)
    phc = get_phcfun()
    vrb = (vrblvl > 0)
    apars = int4a2nbr([len(pol), nvr, idx], vrb)
    bpol = str2int4a(pol)
    ccc = pointer(c_double(0.0))
    vlvl = c_int(vrblvl)
    retval = phc(338, apars, bpol, ccc, vlvl)
    if vrblvl > 0:
        print('-> set_double_double_polynomial return value :', retval)
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
        print("-> set_quad_double_polynomial, pol = ", pol)
    phc = get_phcfun()
    vrb = (vrblvl > 0)
    apars = int4a2nbr([len(pol), nvr, idx], vrb)
    bpol = str2int4a(pol)
    ccc = pointer(c_double(0.0))
    vlvl = c_int(vrblvl)
    retval = phc(388, apars, bpol, ccc, vlvl)
    if vrblvl > 0:
        print('-> set_quad_double_polynomial return value :', retval)
    return retval

def set_double_system(nvr, pols, vrblvl=0):
    """
    Sets the system defines by the strings in pols, in double precision,
    with a number of variables no more than nvr.
    The dimension of the system is set to len(pols).
    Returns the sum of the return values of set_double_polynomial.
    """
    if vrblvl > 0:
        print('-> set_double_system, nvr =', nvr)
    dim = len(pols)
    set_double_dimension(dim, vrblvl)
    result = 0
    for k in range(dim):
        if vrblvl > 0:
            print('-> set_double_system, pols[%d] = %s' % (k, pols[k]))
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
        print('-> set_double_double_system, nvr =', nvr)
    dim = len(pols)
    set_double_double_dimension(dim, vrblvl)
    result = 0
    for k in range(dim):
        if vrblvl > 0:
            print('-> set_double_double_system, pols[%d] = %s' % (k, pols[k]))
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
        print('-> set_quad_double_system, nvr =', nvr)
    dim = len(pols)
    set_quad_double_dimension(dim, vrblvl)
    result = 0
    for k in range(dim):
        if vrblvl > 0:
            print('-> set_quad_double_system, pols[%d] = %s' % (k, pols[k]))
        retval = set_quad_double_polynomial(k+1, nvr, pols[k], vrblvl)
        if vrblvl > 0:
            print('-> set_quad_double_polynomial returns', retval)
        result = result + retval
    return result

def get_double_polynomial(idx, vrblvl=0):
    """
    Returns the string representation of the polynomial
    in double precision, at index idx.
    """
    if vrblvl > 0:
        print("-> get_double_polynomial idx = ", idx)
    phc = get_phcfun()
    adx = pointer(c_int(idx))
    bsz = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(600, adx, bsz, ccc, vrb)
    if vrblvl > 0:
        print('-> the return value of getting the size :', retval)
        print('-> size of the polynomial :', bsz[0])
    szd = bsz[0]
    poldata = create_string_buffer(szd)
    retval = phc(67, adx, poldata, ccc, vrb)
    if vrblvl > 0:
        print('-> the return value of get_double_polynomial :', retval)
    result = int4a2str(poldata, True)
    return result

def get_double_double_polynomial(idx, vrblvl=0):
    """
    Returns the string representation of the polynomial
    in double double precision, at index idx.
    """
    if vrblvl > 0:
        print("-> get_double_double_polynomial idx = ", idx)
    phc = get_phcfun()
    adx = pointer(c_int(idx))
    bsz = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(601, adx, bsz, ccc, vrb)
    if vrblvl > 0:
        print('-> the return value of getting the size :', retval)
        print('-> size of the polynomial :', bsz[0])
    szd = bsz[0]
    poldata = create_string_buffer(szd)
    retval = phc(106, adx, poldata, ccc, vrb)
    if vrblvl > 0:
        print('-> the return value of get_double_double_polynomial :', retval)
    result = int4a2str(poldata, True)
    return result

def get_quad_double_polynomial(idx, vrblvl=0):
    """
    Returns the string representation of the polynomial
    in quad double precision, at index idx.
    """
    if vrblvl > 0:
        print("-> get_quad_double_polynomial idx = ", idx)
    phc = get_phcfun()
    adx = pointer(c_int(idx))
    bsz = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(602, adx, bsz, ccc, vrb)
    if vrblvl > 0:
        print('-> the return value of getting the size :', retval)
        print('-> size of the polynomial :', bsz[0])
    szd = bsz[0]
    poldata = create_string_buffer(szd)
    retval = phc(107, adx, poldata, ccc, vrb)
    if vrblvl > 0:
        print('-> the return value of get_quad_double_polynomial :', retval)
    result = int4a2str(poldata, True)
    return result

def get_double_system(vrblvl=0):
    """
    Returns the string representation of the polynomials
    in double precision.
    """
    dim = get_double_dimension(vrblvl)
    if vrblvl > 0:
        print('-> get_double_system, dim =', dim)
    result = []
    for k in range(1, dim+1):
        pol = get_double_polynomial(k, vrblvl)
        print('-> get_double_system, pol =', pol)
        result.append(pol)
    return result

def get_double_double_system(vrblvl=0):
    """
    Returns the string representation of the polynomials
    in double double precision.
    """
    dim = get_double_double_dimension(vrblvl)
    if vrblvl > 0:
        print('-> get_double_double_system, dim =', dim)
    result = []
    for k in range(1, dim+1):
        pol = get_double_double_polynomial(k, vrblvl)
        print('-> get_double_double_system, pol =', pol)
        result.append(pol)
    return result

def get_quad_double_system(vrblvl=0):
    """
    Returns the string representation of the polynomials
    in quad double precision.
    """
    dim = get_quad_double_dimension(vrblvl)
    if vrblvl > 0:
        print('-> get_quad_double_system, dim =', dim)
    result = []
    for k in range(1, dim+1):
        pol = get_quad_double_polynomial(k, vrblvl)
        print('-> get_quad_double_system, pol =', pol)
        result.append(pol)
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
        print('-> string_of_symbols, maxlen :', maxlen)
    phc = get_phcfun()
    slen = pointer(c_int(0))
    ssym = create_string_buffer(maxlen*4)
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(295, slen, ssym, ccc, vrb)
    if vrblvl > 0:
        print('-> string_of_symbols, return value :', retval)
        print('number of characters :', slen[0])
    vars = int4a2str(ssym, (vrblvl > 0))
    result = vars[0:slen[0]].split(' ')
    return result

def number_of_symbols(pols, vrblvl=0):
    """
    Returns the number of symbols that appear as variables
    in the polynomials, givein in the list of strings pols.
    Useful to determine whether a system is square or not.
    """
    if vrblvl > 0:
        print('-> string_of_symbols, pols :')
        for pol in pols:
            print(pol)
    inpols = ''.join(pols)
    phc = get_phcfun()
    slen = pointer(c_int(len(inpols)))
    bpol = str2int4a(inpols)
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(439, slen, bpol, ccc, vrb)
    if vrblvl > 0:
        print('-> number_of_symbols, return value :', retval)
        print('-> number_of_symbols, result :', slen[0])
    return slen[0]

def test_double_polynomial():
    """
    Tests the setting and getting of a polynomial,
    in double precision.
    """
    lvl = 10
    set_double_dimension(2, lvl)
    dim = get_double_dimension(lvl)
    print('the dimension :', dim)
    org = "x*y - 1;"
    idx = 1
    set_double_polynomial(idx, dim, org, lvl)
    pol = get_double_polynomial(idx, lvl)
    print('the retrieved polynomial :', pol)
    smb = string_of_symbols(100,lvl)
    print('the list of symbols :', smb)

def test_double_double_polynomial():
    """
    Tests the setting and getting of a polynomial,
    in double double precision.
    """
    lvl = 10
    set_double_double_dimension(2, lvl)
    dim = get_double_double_dimension(lvl)
    print('the dimension :', dim)
    org = "x*y - 1;"
    idx = 1
    set_double_double_polynomial(idx, dim, org, lvl)
    pol = get_double_double_polynomial(idx, lvl)
    print('the retrieved polynomial :', pol)

def test_quad_double_polynomial():
    """
    Tests the setting and getting of a polynomial,
    in quad double precision.
    """
    lvl = 10
    set_quad_double_dimension(2, lvl)
    dim = get_quad_double_dimension(lvl)
    print('the dimension :', dim)
    org = "x*y - 1;"
    idx = 1
    set_quad_double_polynomial(idx, dim, org, lvl)
    pol = get_quad_double_polynomial(idx, lvl)
    print('the retrieved polynomial :', pol)

def test_double_system():
    """
    Tests the setting and getting of a system,
    in double precision.
    """
    lvl = 10
    polynomials = ["x^3 + 2*x*y - 1;", "x + y - 1;", "x - 1;"]
    dim = number_of_symbols(polynomials, lvl)
    print('number of symbols :', dim)
    if dim == len(polynomials):
        print('The system is square.')
    else:
        print('number of polynomials :', len(polynomials))
        print('  number of variables :', dim)
        print('The system is not square.')
    set_double_system(dim, polynomials, lvl)
    pols = get_double_system(lvl)
    print('the retrieved polynomials :')
    for pol in pols:
        print(pol)

def test_double_double_system():
    """
    Tests the setting and getting of a system,
    in double double precision.
    """
    lvl = 10
    polynomials = ["x^3 + 2*x*y - 1;", "x + y - 1/3;", "x - 1;"]
    dim = number_of_symbols(polynomials, lvl)
    print('number of symbols :', dim)
    if dim == len(polynomials):
        print('The system is square.')
    else:
        print('number of polynomials :', len(polynomials))
        print('  number of variables :', dim)
        print('The system is not square.')
    set_double_double_system(dim, polynomials, lvl)
    pols = get_double_double_system(lvl)
    print('the retrieved polynomials :')
    for pol in pols:
        print(pol)

def test_quad_double_system():
    """
    Tests the setting and getting of a system,
    in quad double precision.
    """
    lvl = 10
    polynomials = ["x^3 + 2*x*y - 1;", "x + y - 1/3;", "y - 1;"]
    dim = number_of_symbols(polynomials, lvl)
    print('number of symbols :', dim)
    if dim == len(polynomials):
        print('The system is square.')
    else:
        print('number of polynomials :', len(polynomials))
        print('  number of variables :', dim)
        print('The system is not square.')
    set_quad_double_system(dim, polynomials, lvl)
    pols = get_quad_double_system(lvl)
    print('the retrieved polynomials :')
    for pol in pols:
        print(pol)

if __name__=="__main__":
    # test_double_polynomial()
    # test_double_system()
    test_double_double_system()
    # test_quad_double_system()
