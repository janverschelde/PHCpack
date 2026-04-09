"""
Exports working with real powered Laurent homotopies.
The terms of the Laurent polynomials are managed by the polynomials module.
This module manages the arrays of real powered series with complex 
coefficients which correspond to the terms in the Laurent polynomials.
"""
from ctypes import c_int32, c_double, pointer
from phcpy.version import get_phcfun
from phcpy.solver import random_monomial
from phcpy.polynomials import clear_double_laurent_system
from phcpy.polynomials import set_double_laurent_system
from phcpy.polynomials import get_double_laurent_system

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
    szcff = 2*len(pwr)
    cval = (c_double * szcff)()
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
    aval = (c_int32 * 2)()
    aval[0] = adx
    aval[1] = vdx
    idxs = pointer(aval)
    size = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> size_power_vector calls phc', end='')
    retval = phc(938, idxs, size, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return size[0]

def size_coefficient_vector(adx, vdx, vrblvl=0):
    """
    Returns the size of the coefficient vector at position vdx in array adx.
    """
    if vrblvl > 0:
        print('in size_coefficient_vector, adx :', adx, ', vdx :', vdx, '...')
    phc = get_phcfun(vrblvl-1)
    aval = (c_int32 * 2)()
    aval[0] = adx
    aval[1] = vdx
    idxs = pointer(aval)
    size = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> size_power_vector calls phc', end='')
    retval = phc(939, idxs, size, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return size[0]

def to_rps_string(pwrs, cffs, tsb='t', vrblvl=0):
    """
    Given in pwrs the real powers of t and in cffs the corresponding
    list of complex coefficiens, returns the string representation
    of the real powered series.  
    The 'j' in the string of a complex number is replaced by '*i'.
    """
    if vrblvl > 0:
        print('in to_rps_string, pwrs =', pwrs, ', cffs =', cffs, '...')
    res = ''
    for (pwr, cff) in zip(pwrs, cffs):
        jtrm = str(cff)
        itrm = jtrm.replace('j', '*i')
        if pwr == 0.0:
            trm = itrm
        else:
            trm = itrm + '*' + tsb + '**' + str(pwr)
        if res == '':
            res = trm
        else:
            res = res + ' + ' + trm
    return res

def from_rps_string(strrep, tsb='t', vrblvl=0):
    """
    Given in strrep the string representation of a real powered series,
    returns the tuple of two lists, powers and coefficients.
    """
    if vrblvl > 0:
        print('in from_rps_string, strrep =', strrep, '...')
    splitsym = '*' + tsb + '**'
    jstrrep = strrep.replace('*i', 'j')
    data = jstrrep.split(splitsym)
    if vrblvl > 0:
        print('data :', data)
    if len(data) == 0:
        return ([0.0], [eval(data[0])])
    else:
        cff0 = data[0].split(' + ')
        cffs = [eval(cff0[0]), eval(cff0[1])]
        pwrs = [0.0]
        for item in data[1:-1]:
            cff = item.split(' + ')
            pwrs.append(eval(cff[0]))
            cffs.append(eval(cff[1]))
        pwrs.append(eval(data[-1]))
        return (pwrs, cffs)

def get_series_term(adx, vdx, vrblvl=0):
    """
    Gets the powers and coefficients of one term, at the position adx
    in the array and vector with index vdx in the array.
    """
    if vrblvl > 0:
        print('in get_series_term, adx :', adx, ', vdx :', vdx, '...')
    phc = get_phcfun(vrblvl-1)
    powsize = size_power_vector(adx, vdx, vrblvl-1)
    size = pointer(c_int32(powsize))
    bval = (c_int32 * 2)()
    bval[0] = adx
    bval[1] = vdx
    idxs = pointer(bval)
    cval = (c_double * powsize)()
    powers = pointer(cval)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_series_term calls phc', end='')
    retval = phc(940, size, idxs, powers, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    vals = powers[0:powsize]
    pwrs = []
    for idx in range(powsize):
        pwrs.append(float(vals[0][idx]))
    cffsize = size_coefficient_vector(adx, vdx, vrblvl-1)
    size = pointer(c_int32(cffsize))
    bval = (c_int32 * 2)()
    bval[0] = adx
    bval[1] = vdx
    idxs = pointer(bval)
    cffsize = 2*cffsize
    cval = (c_double * cffsize)()
    coefficients = pointer(cval)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_series_term calls phc', end='')
    retval = phc(941, size, idxs, coefficients, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    vals = coefficients[0:cffsize]
    cffs = []
    for idx in range(cffsize//2):
        cff = complex(float(vals[0][2*idx]), float(vals[0][2*idx+1]))
        cffs.append(cff)
    return (pwrs, cffs)

def solve_linear_system(nbr, vrblvl=0):
    """
    Computes the first nbr terms of the solution series of
    a linear system defined by real powered series set in
    the vectors of vectors containers and by a corresponding
    system of Laurent polynomials, in double precision.
    """
    if vrblvl > 0:
        print('in solve_linear_system, nbr :', nbr)
    phc = get_phcfun(vrblvl-1)
    anbr = pointer(c_int32(nbr))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> solve_linear_system calls phc', end='')
    retval = phc(944, anbr, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def run_newton_steps(nbr, vrblvl=0):
    """
    Runs nbr steps of Newton's method to compute the solution series
    real powered Laurent homotopy, defined by the vectors of vectors
    containers and by a corresponding system of Laurent polynomials,
    in double precision.
    """
    if vrblvl > 0:
        print('in run_newton_steps, nbr :', nbr)
    phc = get_phcfun(vrblvl-1)
    anbr = pointer(c_int32(nbr))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> run_newton_steps calls phc', end='')
    retval = phc(945, anbr, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

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

def random_real_powers(deg, vrblvl=0):
    """
    Returns a list of range 0..deg,
    of increasing real powers, starting at zero.
    The second exponent is larger than one,
    and all other powers form an increasing sequence.
    """
    if vrblvl > 0:
        print("in random_real_powers ...")
    from random import random
    result = [0.0, 1.0 + random()]
    for idx in range(deg-1):
        inc = random()
        lst = result[-1]
        result.append(lst + inc)
    return result

def random_complex_coefficients(deg, vrblvl=0):
    """
    Returns a list of range 0..deg,
    of random complex numbers on the unit circle.
    """
    if vrblvl > 0:
        print("in random_complex_coefficients ...")
    from math import pi, sin, cos
    from random import random
    result = []
    for idx in range(deg+1):
        angle = 2*pi*random()
        nbr = complex(sin(angle), cos(angle))
        result.append(nbr)
    return result

def random_real_powered_series(deg, vrblvl=0):
    """
    Returns a tuple of two random vectors
    to represent a real powered series with complex coefficients, 
    truncated at degree deg.
    The first tuple item is the vector of real exponents,
    the second tuple item are the complex coefficients.
    """
    if vrblvl > 0:
        print("in random_real_powered_series ...")
    return (random_real_powers(deg, vrblvl-1),
            random_complex_coefficients(deg, vrblvl-1))

def random_real_powered_polynomial(nvr, nbt, deg, lowexp=-9, uppexp=9, \
    vrblvl=0):
    """
    Returns the string representation of a Laurent polynomial in nvr
    variables with nbt monomials, with exponents in [lowexp, uppexp],
    with coefficients real powered series truncated at degree deg.
    """
    if vrblvl > 0:
        print('in random_real_powered_polynomial, nvr :', nvr, end=', ')
        print('nbt :', nbt, ', lowexp :', lowexp, end=', ')
        print('uppexp :', uppexp)
    pol = ''
    for tix in range(nbt):
        (pwr, cff) = random_real_powered_series(deg, vrblvl-1)
        strrep = to_rps_string(pwr, cff, vrblvl=vrblvl-1)
        if vrblvl > 0:
            print('strrep :', strrep)
        mon = random_monomial(nvr, lowexp, uppexp, vrblvl-1)
        idx = mon.find(')')
        trm = '(' + strrep + ')' + mon[idx+1:]
        if pol == '':
            pol = trm
        else:
            pol = pol + ' + ' + trm
    return pol   

def parse_real_powered_polynomial(pol, vsb='x', vrblvl=0):
    """
    Given the string representation of a real powered polynomial in pol,
    returns a tuple of two lists: the list of string representations of
    the coefficients and the list of monomials.
    The variables are expected to be x1, x2, etc, or whatever the
    symbol of the variable in vsb is.
    """
    if vrblvl > 0:
        print('in parse_real_powered_polynomial, pol :', pol)
    cffs, mons = [], []
    idx = pol.find(vsb)
    c0 = pol[1:idx-2]
    cffs = [c0]
    if vrblvl > 0:
        print('the first coefficient :\n', c0)
    rest = pol[idx:]
    idx = rest.find('+')
    m0 = rest[:idx].strip()
    mons = [m0]
    if vrblvl > 0:
        print('the first monomial :\n', m0)
    rest = rest[idx+1:]
    while(rest != ''):
        if vrblvl > 0:
            print('the rest :\n', rest)
        idx0 = rest.find('(')
        rest = rest[idx0+1:]
        idx1 = rest.find(vsb)
        c1 = rest[:idx1-2]
        cffs.append(c1)
        if vrblvl > 0:
            print('the next coefficient :\n', c1)
        rest = rest[idx1:]
        idx = rest.find('+')
        if idx != -1:
            m1 = rest[:idx].strip()
            mons.append(m1)
            if vrblvl > 0:
                print('the next monomial :\n', m1)
            rest = rest[idx+1:]
        else:
            m1 = rest.strip()
            mons.append(m1)
            if vrblvl > 0:
                print('the last monomial :\n', m1)
            rest = ''
    return (cffs, mons)

def random_real_powered_system(nbq, nvr, nbt, deg, lowexp=-9, uppexp=9, \
    vrblvl=0):
    """
    Returns a list of string representations of a Laurent system
    of nbq polynomials in nvr variables, where the k-th polynomial
    has nbt[k] monomials, with coefficients as power series truncated
    at degree deg, and exponents in [lowexp,uppexp].
    """
    if vrblvl > 0:
        print('in random_real_powered_system, nbq :', nbq, end=', ')
        print('nvr :', nvr, 'nbt :', nbt, ', lowexp :', lowexp, end=', ')
        print('uppexp :', uppexp)
    res = []
    for qix in range(nbq):
        pol = random_real_powered_polynomial\
                  (nvr, nbt[qix], deg, lowexp, uppexp, vrblvl-1)
        res.append(pol)
    return res 

def parse_real_powered_system(pols, vsb='x', vrblvl=0):
    """
    Given the string representations of a real powered system in pols,
    returns a tuple of two lists of lists: 
    (1) string representations of the coefficients; and
    (2) the list of monomial lists for each polynomial.
    The variables are expected to be x1, x2, etc, or whatever the
    symbol of the variable in vsb is.
    """
    if vrblvl > 0:
        print('in parse_real_powered_system, pols :', pols)
    cffs, mons = [], []
    for pol in pols:
        (cff, mon) = parse_real_powered_polynomial(pol, vsb, vrblvl-1)
        cffs.append(cff)
        mons.append(mon)
    return (cffs, mons)

def test_additions(deg, vrblvl=0):
    """
    Fills up initialized arrays with random real powered series
    truncated at degree deg.
    """
    if vrblvl > 0:
        print("in test_additions ...")
    powdim = power_dimension(vrblvl-1)
    powsizes = [size_power_array(i+1,vrblvl-1) for i in range(powdim)]
    print('size power arrays :', powsizes)
    fail = 0
    for (adx, dim) in enumerate(powsizes):
        for vdx in range(dim):
            if vrblvl > 0:
                print('adding series', vdx+1, 'to array', adx+1, '...')
            (pwr, cff) = random_real_powered_series(deg, vrblvl-1)
            fail = fail + set_series_term(adx+1, vdx+1, pwr, cff, vrblvl)
    return fail

def test_vector_dimensions(vrblvl=0):
    """
    Tests retrieving the dimensions of all vectors
    added to the initialized arrays of vectors.
    """
    if vrblvl > 0:
        print("in test_vector_dimensions ...")
    powdim = power_dimension(vrblvl-1)
    cffdim = coefficient_dimension(vrblvl-1)
    powsizes = [size_power_array(i+1,vrblvl-1) for i in range(powdim)]
    if vrblvl > 0:
        print('size power arrays :', powsizes)
    cffsizes = [size_power_array(i+1,vrblvl-1) for i in range(cffdim)]
    if vrblvl > 0:
        print('size coefficient arrays :', cffsizes)
    fail = 0
    for (adx, dim) in enumerate(powsizes):
        for vdx in range(dim):
            size = size_power_vector(adx+1, vdx+1, vrblvl)
            if vrblvl > 0:
                print('power vector', vdx+1, 'of array', adx+1, end=' ')
                print('has size', size)
            size = size_coefficient_vector(adx+1, vdx+1, vrblvl)
            if vrblvl > 0:
                print('coefficient vector', vdx+1, end=' ')
                print('of array', adx+1, 'has size', size)
    return fail

def test_retrievals(vrblvl=0):
    """
    Retrieves all vectors stored in the arrays.
    """
    if vrblvl > 0:
        print("in test_retrievals ...")
    powdim = power_dimension(vrblvl-1)
    powsizes = [size_power_array(i+1,vrblvl-1) for i in range(powdim)]
    print('size power arrays :', powsizes)
    fail = 0
    for (adx, dim) in enumerate(powsizes):
        for vdx in range(dim):
            if vrblvl > 0:
                print('getting series', vdx+1, 'to array', adx+1, '...')
            (pwr, cff) = get_series_term(adx+1, vdx+1, vrblvl)
            if vrblvl > 0:
                print('powers :', pwr)
                print('coefficients :', cff)
    return fail

def test_string_representations(vrblvl=0):
    """
    Retrieves all vectors stored in the arrays
    and writes the string representations of the series.
    """
    if vrblvl > 0:
        print("in test_retrievals ...")
    powdim = power_dimension(vrblvl-1)
    powsizes = [size_power_array(i+1,vrblvl-1) for i in range(powdim)]
    print('size power arrays :', powsizes)
    fail = 0
    for (adx, dim) in enumerate(powsizes):
        for vdx in range(dim):
            if vrblvl > 0:
                print('getting series', vdx+1, 'to array', adx+1, '...')
            (pwr, cff) = get_series_term(adx+1, vdx+1, vrblvl)
            strrep = to_rps_string(pwr, cff, vrblvl=vrblvl)
            print('strrep :', strrep)
            (pwrs, cffs) = from_rps_string(strrep, vrblvl=vrblvl)
            print('powers :', pwrs)
            print('coefficients :', cffs)
    return fail

def test_random_system(deg, vrblvl=0):
    """
    Generates a random real powered system,
    where the coefficients are series truncated at degree deg.
    """
    if vrblvl > 0:
        print('in test_random_system ...')
    pols = random_real_powered_system(3, 2, [1, 2, 3], 2) 
    print('a random real powered system :')
    for (idx, pol) in enumerate(pols):
        print('polynomial', idx+1, ':', pol)
    return 0

def test_parse_polynomial(deg, vrblvl=0):
    """
    Generates a random real powered polynomial where the coefficient
    series are truncated at degree deg and then extracts the series
    coefficients and the monomials.
    """
    pol = random_real_powered_polynomial(3, 3, deg, vrblvl=vrblvl-1)
    print('a random polynomial :\n', pol)
    (cffs, mons) = parse_real_powered_polynomial(pol, vrblvl=vrblvl)
    print('the coefficients :\n', cffs)
    print('the monomials :\n', mons)
    return 0

def test_parse_system(deg, vrblvl=0):
    """
    Generates a random real powered system where the coefficient
    series are truncated at degree deg and then extracts the series
    coefficients and the monomials.
    """
    pols = random_real_powered_system(3, 3, [2, 2, 2], deg, vrblvl=vrblvl-1)
    print('a random system :\n', pols)
    (cffs, mons) = parse_real_powered_system(pols, vrblvl=vrblvl)
    for (idx, pol) in enumerate(pols):
        print('polynomial', idx+1, ':\n', pol)
        print('has coefficients :\n', cffs[idx])
        print('and monomials :\n', mons[idx])
    return 0

def linear_laurent_system(dim, vrblvl=0):
    """
    Sets the linear Laurent system of dimension dim.
    """
    if vrblvl > 0:
        print('in linear_laurent_system, dim :', dim, '...')
    clear_double_laurent_system(vrblvl-1)
    variables = ['x' + str(k) for k in range(1, dim+1)]
    poly = ' + '.join(variables) + ' + 1;'
    pols = [poly for _ in range(dim)]
    if vrblvl > 0:
        print('the polynomials :')
        for pol in pols:
            print(pol)
    set_double_laurent_system(dim, pols, vrblvl-1)
    if vrblvl > 0:
        gpol = get_double_laurent_system(vrblvl-1)
        print('the stored polynomials :')
        for pol in gpol:
            print(pol)

def random_linear_matrix(dim, deg, vrblvl=0):
    """
    Generates a dim-by-dim matrix of random real powered series,
    returned as a list of rows, where each element on the row
    is a tuple of two lists: the powers and coefficients,
    truncated at the term of index deg.
    """
    if vrblvl > 0:
        print('in random_linear_matrix, dim :', dim, ', deg :', deg, '...')
    res = []
    for idx in range(dim):
        row = [random_real_powered_series(deg, vrblvl-1) for _ in range(dim)]
        res.append(row)
    return res

def random_series_vector(dim, deg, vrblvl=0):
    """
    Generates a vector of size dim with powers and coefficients
    of a series truncated at term at index deg.
    """
    if vrblvl > 0:
        print('in random_linear_matrix, dim :', dim, ', deg :', deg, '...')
    res = [random_real_powered_series(deg, vrblvl-1) for _ in range(dim)]
    return res

def sort_series_powers(pwr, cff, vrblvl=0):
    """
    Given in pwr the unsorted exponents of a series
    and in cff the corresponding coefficients.
    Sorts the exponents in pwr, swapping the coefficients in cff.
    The sorting happens in place, the function modifies pwr and cff,
    returning nothing.
    """
    if vrblvl > 0:
        print('in sort_series_powers ...')
    for idx1 in range(len(pwr)):
        minidx = idx1
        for idx2 in range(idx1+1,len(pwr)):
            if pwr[idx2] < pwr[minidx]:
                minidx = idx2
        if not(minidx == idx1):
            pwr[minidx], pwr[idx1] = pwr[idx1], pwr[minidx]
            cff[minidx], cff[idx1] = cff[idx1], cff[minidx]

def normalize_series_powers(pwr, cff, vrblvl=0):
    """
    Given in pwr are sorted exponents of power series,
    with corresponding coefficients in cff.
    For equal powers, the corresponding coefficients need
    to be added and the lists need be shortened.
    """
    if vrblvl > 0:
        print('in normalize_series_powers ...')
    for idx in range(len(pwr)):
        if idx == len(pwr)-1:
            break
        if pwr[idx] == pwr[idx+1]:
            cff[idx] = cff[idx] + cff[idx+1]
            del cff[idx+1]
            del pwr[idx+1]

def series_product(mat, vec, vrblvl=0):
    """
    Given in mat a square matrix of real powered series
    and in vec a real powered series vector of compatible dimension,
    all series given as tuples of powers and coefficients,
    returns the vector of real powered series which is the product
    of the matrix with the vector.
    """
    if vrblvl > 0:
        print('in series_product ...')
    res = []
    dim = len(mat)
    for rowidx in range(dim):
        respwr = []
        rescff = []
        for colidx in range(dim):
            (apwr, acff) = mat[rowidx][colidx]
            (vpwr, vcff) = vec[colidx]
            for adx in range(len(apwr)):
                for vdx in range(len(vpwr)):
                    respwr.append(apwr[adx] + vpwr[vdx])
                    rescff.append(acff[adx] * vcff[vdx])
        res.append((respwr, rescff))
    for (pwr, cff) in res:
        if vrblvl > 0:
            print('-> before the sort ...')
            print('powers :', pwr)
            print('coeffs :', cff)
        sort_series_powers(pwr, cff, vrblvl-1)
        if vrblvl > 0:
            print('-> after sorting ...')
            print('powers :', pwr)
            print('coeffs :', cff)
        normalize_series_powers(pwr, cff, vrblvl-1)
        if vrblvl > 0:
            print('-> after normalizing ...')
            print('powers :', pwr)
            print('coeffs :', cff)
            print('len(pwr) =', len(pwr), '==', len(cff), end=' = ')
            print('len(cff) ?', len(pwr)==len(cff))
    return res

def random_linear_system(dim, deg, vrblvl=0):
    """
    Makes a random linear system of dimension dim,
    of real powered series truncated at index deg,
    and sets the coefficients and the Laurent system.
    """
    if vrblvl > 0:
        print('in random_linear_system, dim :', dim, '...')
    linear_laurent_system(dim, vrblvl)
    mat = random_linear_matrix(dim, deg, vrblvl)
    if vrblvl > 0:
        print('the matrix :')
        for (idx, row) in enumerate(mat):
            print('series at row', idx+1, ':')
            for (pwr, cff) in row:
                print('powers :', pwr)
                print('coefficients :', cff)
    clear_series_terms(vrblvl-1)
    dims = [dim+1 for _ in range(dim)]
    initialize_series_coefficients(dims, vrblvl-1)
    fail = 0
    for (adx, row) in enumerate(mat):
        vdx = 0
        for (pwr, cff) in row:
            fail = fail + set_series_term(adx+1, vdx+1, pwr, cff, vrblvl-1)
            vdx = vdx + 1
    sol = random_series_vector(dim, deg, vrblvl)
    rhs = series_product(mat, sol, vrblvl)
    for adx in range(len(mat)): 
        (pwr, cff) = rhs[adx]
        fail = fail + set_series_term(adx+1, dim+1, pwr, cff, vrblvl-1)
    return fail

def test_linear_solver(dim, deg, vrblvl=0):
    """
    Tests the linear solver on a generated system of dimension dim,
    with real powered series truncated at index deg.
    """
    if vrblvl > 0:
        print('in test_linear_system, dim :', dim, '...')
    fail = random_linear_system(dim, deg, vrblvl)
    if fail != 0:
        print(fail, 'failures occurred in defining the linear system!')
    res = solve_linear_system(3, vrblvl)
    return res

def test_newton_steps(vrblvl=0):
    """
    Tests newton's method.
    """
    if vrblvl > 0:
        print("in test_newton_steps ...")
    res = run_newton_steps(3, vrblvl)
    return res

def test_laurent(deg, vrblvl=0):
    """
    Tests operations in this module.
    The verbose level is defined by vrblvl.
    """
    if vrblvl > 0:
        print("in test_laurent ...")
    fail = test_initialization(vrblvl)
    fail = fail + test_dimensions(vrblvl)
    fail = fail + test_additions(deg, vrblvl)
    fail = fail + test_vector_dimensions(vrblvl)
    fail = fail + test_retrievals(vrblvl)
    fail = fail + test_string_representations(vrblvl)
    fail = fail + test_random_system(deg, vrblvl)
    fail = fail + test_parse_polynomial(deg, vrblvl)
    fail = fail + test_parse_system(deg, vrblvl)
    fail = fail + test_linear_solver(deg, deg, vrblvl)
    fail = fail + test_newton_steps(vrblvl)
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
    deg = 2
    fail = test_laurent(deg, lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__== '__main__':
    main()
