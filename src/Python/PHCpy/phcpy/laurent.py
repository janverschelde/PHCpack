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
        print('in laurent.initialize_series_coefficients, dims :', dims)
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
        print('in laurent.set_series_term, adx :', adx, ', vdx :', vdx)
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
        print('in laurent.power_dimension ...')
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
        print('in laurent.coefficient_dimension ...')
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
        print('in laurent.size_power_array, adx :', adx, '...')
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
        print('in laurent.size_coefficient_array, adx :', adx, '...')
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
        print('in laurent.size_power_vector, adx :', adx, end='')
        print(', vdx :', vdx, '...')
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
        print('in laurent.size_coefficient_vector, adx :', adx, end='')
        print(', vdx :', vdx, '...')
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

def get_series_term(adx, vdx, vrblvl=0):
    """
    Gets the powers and coefficients of one term, at the position adx
    in the array and vector with index vdx in the array.
    """
    if vrblvl > 0:
        print('in laurent.get_series_term, adx :', adx, ', vdx :', vdx, '...')
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
        print('in laurent.solve_linear_system, nbr :', nbr)
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

def split_powers_coefficients(dim, numbers, vrblvl=0):
    """
    Given a sequence of doubles in the list numbers, as
    [cre0, cim0, cre1, cim1, pwr1, cre2, cim2, pwr2, ...],
    with real and imaginary parts in cre, cim, powers in pwr,
    for the dim many components of the series,
    returns two lists: real powers and complex coefficients.
    """
    if vrblvl > 0:
        print('in laurent.split_powers_coefficients, dim :', dim)
    pwrs, cffs = [0.0 for _ in range(dim)], []
    recff, imcff = 0.0, 0.0
    for (idx, number) in enumerate(numbers):
        if idx < 2*dim:
            if idx % 2 == 0: 
                recff = number
            else:
                imcff = number
                cff = complex(recff, imcff)
                cffs.append(cff)
        else:
            if (idx - 2*dim) % 3 == 0:
                recff = number
            elif (idx - 2*dim) % 3 == 1:
                imcff = number
                cff = complex(recff, imcff)
                cffs.append(cff)
            else:
                pwrs.append(number)
    return pwrs, cffs

def distribute_powers(dim, powers, vrblvl=0):
    """
    Given one list of powers for every component
    of a dim-dimensional vector of series,
    returns dim many lists with the powers for each component,
    using the values in powers.
    """
    if vrblvl > 0:
        print('in laurent.distribute_powers, dim :', dim)
        if not (len(powers) % dim == 0):
            print('size of powers :', len(powers), end=' ')
            print('is not a multiple of', dim, '!')
    result = [[] for _ in range(dim)]
    pwridx = 0
    size = len(powers)//dim
    for column in range(size):
        for row in range(dim):
            result[row].append(powers[pwridx])
            pwridx = pwridx + 1
            if pwridx > len(powers):
                if vrblvl > 0:
                     print('exceeding', len(powers), ', returning ...')
                return result
    return result

def distribute_coefficients(dim, coefficients, vrblvl=0):
    """
    Given one list of coefficients for every component
    of a dim-dimensional vector of series,
    returns dim many lists with coefficients for each component,
    using the values in coefficients.
    """
    if vrblvl > 0:
        print('in laurent.distribute_coefficients, dim :', dim)
        if not (len(coefficients) % (dim) == 0):
            print('size of coefficients :', len(coefficients), end=' ')
            print('is not a multiple of', dim, '!')
    result = [[] for _ in range(dim)]
    cffidx = 0
    size = len(coefficients)//dim
    for column in range(size):
        for row in range(dim):
            result[row].append(coefficients[cffidx])
            cffidx = cffidx + 1
            if cffidx > len(coefficients):
                if vrblvl > 0:
                     print('exceeding', len(coefficients), ', returning ...')
                return result
    return result

def run_newton_steps(dim, nbr, vrblvl=0):
    """
    Runs nbr steps of Newton's method to compute the solution series
    real powered Laurent homotopy of dimension dim,
    defined by the vectors of vectors containers and by a corresponding 
    system of Laurent polynomials, in double precision.
    """
    if vrblvl > 0:
        print('in laurent.run_newton_steps, nbr :', nbr)
    phc = get_phcfun(vrblvl-1)
    anbr = pointer(c_int32(nbr))
    bbb = pointer(c_int32(0))
    # for every term we have three doubles, except for constant
    outsize = 2*dim + 3*dim*nbr
    cval = (c_double * outsize)()
    for idx in range(outsize):
        cval[idx] = 0.0
    result = pointer(cval)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> run_newton_steps calls phc', end='')
    retval = phc(945, anbr, bbb, result, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    vals = result[0:outsize]
    numbers = []
    for idx in range(outsize):
        numbers.append(float(vals[0][idx]))
    if vrblvl > 0:
        print('returned numbers :', numbers)
    pwrs, cffs = split_powers_coefficients(dim, numbers, vrblvl-1)
    if vrblvl > 0:
        print('coefficients :\n', cffs)
        print('powers :\n', pwrs)
    powers = distribute_powers(dim, pwrs, vrblvl-1)
    coefficients = distribute_coefficients(dim, cffs, vrblvl-1)
    if vrblvl > 0:
        print('distributed powers :\n', powers)
        print('distributed coefficients :\n', coefficients)
    return (powers, coefficients)

def clear_series_terms(vrblvl=0):
    """
    Deallocates the space occupied by the real powers
    and the complex coefficients.
    """
    if vrblvl > 0:
        print('in laurent.clear_series_terms ...')
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
        print("in laurent.test_initialization ...")
    dims = [3, 2, 4]
    fail = initialize_series_coefficients(dims, vrblvl)
    return fail

def test_dimensions(vrblvl=0):
    """
    Tests retrieving the dimensions.
    """
    if vrblvl > 0:
        print("in laurent.test_dimensions ...")
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
        print("in laurent.random_real_powers ...")
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
        print("in laurent.random_complex_coefficients ...")
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
        print("in laurent.random_real_powered_series ...")
    return (random_real_powers(deg, vrblvl-1),
            random_complex_coefficients(deg, vrblvl-1))

def to_rps_string(pwrs, cffs, tsb='t', vrblvl=0):
    """
    Given in pwrs the real powers of t and in cffs the corresponding
    list of complex coefficiens, returns the string representation
    of the real powered series.  
    The 'j' in the string of a complex number is replaced by '*i'.
    """
    if vrblvl > 0:
        print('in laurent.to_rps_string, pwrs =', pwrs, end='')
        print(', cffs =', cffs, '...')
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
        print('in laurent.from_rps_string, strrep =', strrep, '...')
    splitsym = '*' + tsb + '**'
    jstrrep = strrep.replace('*i', 'j')
    data = jstrrep.split(splitsym)
    if vrblvl > 0:
        print('data :', data)
    if len(data) == 0:
        return ([0.0], [eval(data[0])])
    else:
        cff0 = data[0].split(' + ')
        if len(cff0) == 1:
            cffs = [eval(cff0[0])]
        else:
            cffs = [eval(cff0[0]), eval(cff0[1])]
        pwrs = [0.0]
        for item in data[1:-1]:
            cff = item.split(' + ')
            pwrs.append(eval(cff[0]))
            cffs.append(eval(cff[1]))
        pwrs.append(eval(data[-1]))
        return (pwrs, cffs)

def random_real_powered_polynomial\
    (nvr, nbt, deg, zerocst=False, xsb='x', lowexp=-9, uppexp=9, vrblvl=0):
    """
    Returns the string representation of a Laurent polynomial in nvr
    variables with nbt monomials, with exponents in [lowexp, uppexp],
    with coefficients real powered series truncated at degree deg.
    Every variable name starts with xsb.  If zerocst is True,
    then the constant coefficient of all series coefficients is zero.
    """
    if vrblvl > 0:
        print('in laurent.random_real_powered_polynomial, nvr : ', end='')
        print(nvr, ', nbt :', nbt, ', lowexp :', lowexp, end=', ')
        print('uppexp :', uppexp)
    pol = ''
    for tix in range(nbt):
        (pwr, cff) = random_real_powered_series(deg, vrblvl-1)
        if zerocst:
            cff[0] = complex(0.0, 0.0)
        strrep = to_rps_string(pwr, cff, vrblvl=vrblvl-1)
        if vrblvl > 0:
            print('strrep :', strrep)
        mon = random_monomial(nvr, xsb, lowexp, uppexp, vrblvl-1)
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
        print('in laurent.parse_real_powered_polynomial, pol :', pol)
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

def random_real_powered_system\
    (nbq, nvr, nbt, deg, \
    zerocst=False, xsb='x', lowexp=-9, uppexp=9, vrblvl=0):
    """
    Returns a list of string representations of a Laurent system
    of nbq polynomials in nvr variables, where the k-th polynomial
    has nbt[k] monomials, with coefficients as power series truncated
    at degree deg, and exponents in [lowexp,uppexp].
    """
    if vrblvl > 0:
        print('in laurent.random_real_powered_system, nbq :', nbq, end=', ')
        print('nvr :', nvr, 'nbt :', nbt, ', deg :', deg, end=', ')
        print('zerocst :', zerocst, 'xsb :', xsb, end=', ')
        print('lowexp :', lowexp, ', uppexp :', uppexp)
    res = []
    for qix in range(nbq):
        pol = random_real_powered_polynomial\
                  (nvr, nbt[qix], deg, zerocst, xsb, lowexp, uppexp, vrblvl-1)
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
        print('in laurent.parse_real_powered_system, pols :', pols)
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
        print('in laurent.test_additions ...')
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
        print('in laurent.test_vector_dimensions ...')
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
        print('in laurent.test_retrievals ...')
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
        print('in laurent.test_retrievals ...')
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
        print('in laurent.test_random_system ...')
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
    if vrblvl > 0:
        print('in laurent.test_parse polynomial ...')
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
    if vrblvl > 0:
        print('in laurent.test_parse_system ...')
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
        print('in laurent.linear_laurent_system, dim :', dim, '...')
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
        print('in laurent.random_linear_matrix, dim :', dim, end='')
        print(', deg :', deg, '...')
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
        print('in laurent.random_linear_matrix, dim :', dim, end='')
        print(', deg :', deg, '...')
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
        print('in laurent.sort_series_powers ...')
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
        print('in laurent.normalize_series_powers ...')
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
        print('in laurent.series_product ...')
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
        print('in laurent.random_linear_system, dim :', dim, '...')
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

def degree_monomial(mon, xsb='x', vrblvl=0):
    """
    Given the string representation of a monomial using xsb plus index
    as variable name, returns the degree of the monomials.
    """
    if vrblvl > 0:
        print('in laurent.degree_monomial, mon :', mon)
    if not xsb in mon:
        return 0                   # the monomial is a constant
    elif not '*' in mon:
        if not '^' in mon:
            return 1               # only one variable without exponent
        else:
            L = mon.split('^')
            return int(L[1])       # only one exponent
    else:
        factors = mon.split('*')   # we have at least two factors
        result = 0
        for factor in factors:     # add exponents of the factors
            L = factor.split('^')
            result = result + int(L[1])
        return result

def support_monomial(mon, dim, xsb='x', vrblvl=0):
    """
    Given the string representation of a monomial using xsb plus index
    as variable name, in dim many variables,
    returns the list of degrees of the variables.
    """
    if vrblvl > 0:
        print('in laurent.support_monomial, mon :', mon)
    if not xsb in mon:
        return [0 for _ in range(dim)]     # monomial is a constant
    elif not '*' in mon:
        result = [0 for _ in range(dim)]   # only one variable
        idxvar = int(mon[1:])
        if not '^' in mon:
            result[idxvar-1] = 1
            return result          # only one variable without exponent
        else:
            L = mon.split('^')     
            result[idxvar-1] = int(L[1])  # only one exponent
            return result
    else:
        factors = mon.split('*')   # we have at least two factors
        result = [0 for _ in range(dim)]   # exponent space
        for factor in factors:     # add exponents of the factors
            L = factor.split('^')
            idxvar = int(L[0][1:])
            result[idxvar-1] = int(L[1])
        return result

def multiply_monomial(supmon, idxvar, xsb='x', vrblvl=0):
    """
    Given in supmon is the support of a monomial and in idxvar
    the index of a variable, where each varable starts with xsb,
    followed by the index plus one, return the string representation
    of the monomial multiplied by the variable.
    """
    if vrblvl > 0:
        print('in laurent.multiply_monomial, supmon :', supmon)
    result = ''
    for (idx, degree) in enumerate(supmon):
        if idx > 0:
            result = result + '*'
        result = result + xsb + str(idx+1) + '^'
        if idx == idxvar:
            result = result + str(degree+1)
        else:
            result = result + str(degree)
    return result

def sort_monomial_series(mons, cffs, xsb='x', vrblvl=0):
    """
    Sorts the monomials in mons according to decreasing degrees,
    swapping also the corresponding coefficients in cffs.
    """
    if vrblvl > 0:
        print('in laurent.sort_monomial_series ...')
        print('before the sort :')
        print('monomials :', mons)
        print('coefficients :', cffs)
    for idx1 in range(len(mons)):
        maxidx = idx1
        maxdeg = degree_monomial(mons[idx1], xsb)
        for idx2 in range(idx1+1, len(mons)):
            mondeg = degree_monomial(mons[idx2], xsb)
            if mondeg > maxdeg:
                maxidx, maxdeg = idx2, mondeg
        if not(maxidx == idx1):
            mons[maxidx], mons[idx1] = mons[idx1], mons[maxidx]
            cffs[maxidx], cffs[idx1] = cffs[idx1], cffs[maxidx]
    if vrblvl > 0:
        print('after the sort :')
        print('monomials :', mons)
        print('coefficients :', cffs)


def random_product_homotopy(dim, nbt, deg, sol, xsb='x', vrblvl=0):
    """
    Returns the coefficients and monomials, in lists of strings,
    for a Laurent homotopy where the k-th polynomial is the
    product of x[k] - sol[k] times a random polynomial of nbt[k] terms,
    with coefficients that are power series truncated at the first index.
    The index at which the series in sol are truncated equals deg.
    """
    if vrblvl > 0:
        print('in laurent.random_product_homotopy ...')
    pols = random_real_powered_system\
              (dim, dim, nbt, 1, True, xsb, vrblvl=vrblvl-1)
    if vrblvl > 0:
        print('a random system :\n', pols)
    (cffs, mons) = parse_real_powered_system(pols, vrblvl=vrblvl)
    if vrblvl > 0:
        for (idx, pol) in enumerate(pols):
            print('polynomial', idx+1, ':\n', pol)
            print('has coefficients :\n', cffs[idx])
            print('and monomials :\n', mons[idx])
            sup = [support_monomial(mon, dim, xsb, vrblvl) for mon in mons[idx]]
            print('support', idx+1, ':', sup)
    expmons = [ [] for _ in range(len(pols)) ] 
    expcffs = [ [] for _ in range(len(pols)) ] 
    for (idx, pol) in enumerate(pols):
        (solpwr, solcff) = sol[idx]
        for (mon, polcff) in zip(mons[idx], cffs[idx]):
            supmon = support_monomial(mon, dim, xsb, vrblvl)
            mulmon = multiply_monomial(supmon, idx, xsb, vrblvl)
            expmons[idx].append(mulmon)
            expcffs[idx].append(polcff)   # copy coefficient
            if vrblvl > 0:
                print('multiplied monomial :', mulmon)
            (pwrpolcff, cffpolcff) = from_rps_string(polcff)
            if vrblvl > 0:
                print('multiplying solution powers', solpwr)
                print('  with coefficients power :', pwrpolcff[1])
                print('multiplying solution coefficients', solcff)
                print('  with coefficient :', cffpolcff[1])
            prodpwrs = [(pwrpolcff[1] + power) for power in solpwr]
            prodcffs = [(-cffpolcff[1]*coeff) for coeff in solcff]
            if vrblvl > 0:
                print('the product of powers :\n', prodpwrs)
                print('the corresponding coefficients :\n', prodcffs)
                print('len check :', (len(prodpwrs) == len(prodcffs)))
                if not (len(prodpwrs) == len(prodcffs)):
                    print('lengths do not match!')
            strcff = to_rps_string(prodpwrs, prodcffs, vrblvl=vrblvl-1)
            expcffs[idx].append(strcff)
            expmons[idx].append(mon)
        if vrblvl > 0:
            print('expanded monomials :', expmons[idx])
            print('expanded coefficients :', expcffs[idx])
    return (expcffs, expmons)

def random_binomial_homotopy(dim, nbt, deg, sol, xsb='x', vrblvl=0):
    """
    Makes a random Laurent homotopy of dimension dim,
    with number of terms in the polynomials given in the list nbt,
    which is a list of length dim, with series truncated at the
    terms with index deg.  The solution is given in sol.
    Returns the coefficients and the monomials, in lists of strings.
    """
    if vrblvl > 0:
        print('in laurent.random_binomial_homotopy ...')
    ssol = []
    for ksol in sol:
        (pwr, cff) = ksol
        cff[1] = complex(0.0, 0.0) # only a constant
        pwr[1] = 1.0
        flippedsign = [-nbr for nbr in cff]
        ssol.append(to_rps_string(pwr, flippedsign))
    if vrblvl > 0:
        print('the solution series strings :', ssol)
    pols = random_real_powered_system\
              (dim, dim, nbt, deg, True, xsb, vrblvl=vrblvl-1)
    if vrblvl > 0:
        print('a random system :\n', pols)
    (cffs, mons) = parse_real_powered_system(pols, vrblvl=vrblvl)
    for idx in range(len(pols)):
        mons[idx].append(xsb + str(idx+1))
        cffs[idx].append('(1.0+0.0*i) + (0.0+0.0*i)*t**1.0')
        mons[idx].append('1')
        if vrblvl > 0:
            print('appending', ssol[idx])
        cffs[idx].append(ssol[idx])
    for idx in range(len(pols)):
        sort_monomial_series(mons[idx], cffs[idx], xsb, vrblvl)
    if vrblvl > 0:
        for (idx, pol) in enumerate(pols):
            print('polynomial', idx+1, ':\n', pol)
            print('has coefficients :\n', cffs[idx])
            print('and monomials :\n', mons[idx])
            degmon = [degree_monomial(mon, xsb, vrblvl) for mon in mons[idx]]
            print('monomial degrees :', degmon)
    return (cffs, mons)

def store_laurent_homotopy(cffs, mons, vrblvl=0):
    """
    Stores the coefficients and the monomials given as string representations
    of a Laurent homotopy.
    """
    if vrblvl > 0:
        print('in laurent.store_binomial_homotopy ...')
    clear_series_terms(vrblvl-1)
    dims = [len(cff) for cff in cffs]
    if vrblvl > 0:
        print('the dimensions :', dims)
    initialize_series_coefficients(dims, vrblvl-1)
    fail = 0
    for (adx, polcff) in enumerate(cffs):
        if vrblvl > 0:
            print('parsing coefficients of polynomial ', adx+1, '...')
        for (vdx, coefficient) in enumerate(polcff):
            (serpwrs, sercffs) = from_rps_string(coefficient, vrblvl=vrblvl)
            if vrblvl > 0:
                print('powers :', serpwrs)
                print('coeffs :', sercffs)
            fail = fail \
                 + set_series_term(adx+1, vdx+1, serpwrs, sercffs, vrblvl)
    clear_double_laurent_system(vrblvl-1)
    pols = []
    for monomials in mons:
        poly = ' + '.join(monomials) + ';'
        pols.append(poly)
    if vrblvl > 0:
        print('the polynomials :')
        for pol in pols:
            print(pol)
    set_double_laurent_system(len(pols), pols, vrblvl-1)
    if vrblvl > 0:
        gpol = get_double_laurent_system(vrblvl-1)
        print('the stored polynomials :')
        for pol in gpol:
            print(pol)
    return fail

def laurent_homotopy_strings(cffs, mons, vrblvl=0):
    """
    Returns the list of string representations of the polynomials
    in the Laurent homotopy, defined by coefficients and monomials,
    given by their string representations.
    """
    if vrblvl > 0:
        print("in laurent.laurent_homotopy_strings ...")
    result = []
    for (idx, (cff, mon)) in enumerate(zip(cffs, mons)):
        if vrblvl > 0:
            print('polynomial', idx+1, 'has monomials', mon)
            print('and coefficients :\n', cff)
        pol = ''
        for (coefficient, monomial) in zip(cff, mon):
            newcff = coefficient.replace('*i', 'j')
            if pol == '':
                pol = '(' + newcff + ')'
            else:
                pol = pol + ' + (' + newcff + ')'
            if not(monomial == '1'):
                newmon = monomial.replace('^', '**')
                pol = pol + '*' + newmon
        if vrblvl > 0:
            print('polynomial', idx+1, ' : ', pol)
        result.append(pol)
    return result

def evaluate_series(pwr, cff, tval, vrblvl=0):
    """
    Given a series with powers in pwr and corresponding
    coefficients in cff, evaluates the series at t = tval.
    """
    if vrblvl > 0:
        print('in laurent.evaluate_series, tval :', tval)
    value = 0.0
    for (power, coefficient) in zip(pwr, cff):
        value = value + coefficient*tval**power
    return value

def evaluate_series_vector(pwrs, cffs, tval, vrblvl=0):
    """
    Given a series vector with powers in pwrs and corresponding
    coefficients in cffs, evaluates the series at t = tval.
    """
    if vrblvl > 0:
        print('in laurent.evaluate_series_vector, tval :', tval)
    values = []
    for (powers, coefficients) in zip(pwrs, cffs):
        value = evaluate_series(powers, coefficients, tval, vrblvl-1)
        values.append(value)
    return values

def evaluate_laurent_homotopy(lhom, xsol, tval, xsb='x', vrblvl=0):
    """
    Given in hom are the string representations of a Laurent homotopy,
    in xsol are the complex values for solution series evaluated at tval,
    for variables with names starting with xsb, returns the values
    obtained by replacing symbols in hom by values in xsol and tval.
    """
    if vrblvl > 0:
        print('in laurent.evaluate_laurent_homotopy, tval :', tval)
        print('xsol :\n', xsol)
        print('Laurent homotopy :\n', lhom)
    thom = [pol.replace('t', str(tval)) for pol in lhom]
    pval = []
    for pol in thom:
        for (idx, xval) in enumerate(xsol):
            var = xsb + str(idx+1)
            pol = pol.replace(var, str(xval))
        pval.append(pol)
    if vrblvl > 0:
        print('evaluating the expressions :\n', pval)
    result = [eval(value) for value in pval]
    return result

def test_linear_solver(dim, deg, vrblvl=0):
    """
    Tests the linear solver on a generated system of dimension dim,
    with real powered series truncated at index deg.
    """
    if vrblvl > 0:
        print('in laurent.test_linear_system, dim :', dim, '...')
    fail = random_linear_system(dim, deg, vrblvl)
    if fail != 0:
        print(fail, 'failures occurred in defining the linear system!')
    res = solve_linear_system(3, vrblvl)
    return res

def test_newton_product(vrblvl=0):
    """
    Tests newton's method on a random product homotopy.
    """
    if vrblvl > 0:
        print("in laurent.test_newton_product ...")
    dim, deg = 2, 2
    sol = random_series_vector(dim, deg, vrblvl-1)
    if vrblvl > 0:
        print('the solution series :', sol)
    nbt = [2 for _ in range(dim)]
    (cffs, mons) = random_product_homotopy(dim, nbt, deg, sol, vrblvl=vrblvl)
    fail = store_laurent_homotopy(cffs, mons, vrblvl)
    if fail != 0:
        print(fail, 'failures occurred in storing Laurent homotopy!')
    lauhom = laurent_homotopy_strings(cffs, mons, vrblvl)
    print('the polynomials in the Laurent homotopy :')
    for (idx, pol) in enumerate(lauhom):
        print('polynomial', idx+1, ':\n', pol)
    return 0

def test_newton_binomial(vrblvl=0):
    """
    Tests newton's method on a random binomial homotopy.
    """
    if vrblvl > 0:
        print("in laurent.test_newton_binomial ...")
    dim, deg = 2, 2
    sol = random_series_vector(dim, 1, vrblvl-1)
    if vrblvl > 0:
        print('the solution series :', sol)
    nbt = [2 for _ in range(dim)]
    (cffs, mons) = random_binomial_homotopy(dim, nbt, deg, sol, vrblvl=vrblvl)
    fail = store_laurent_homotopy(cffs, mons, vrblvl)
    if fail != 0:
        print(fail, 'failures occurred in storing Laurent homotopy!')
    lauhom = laurent_homotopy_strings(cffs, mons, vrblvl)
    print('the polynomials in the Laurent homotopy :')
    for (idx, pol) in enumerate(lauhom):
        print('polynomial', idx+1, ':\n', pol)
    (pwrs, cffs) = run_newton_steps(dim, 4, vrblvl)
    for (idx, (pwr, cff)) in enumerate(zip(pwrs, cffs)):
        print('series component', idx+1, ':', to_rps_string(pwr, cff))
    for tval in [1.0e-4, 1.0e-6, 1.0e-8]:
        xval = evaluate_series_vector(pwrs, cffs, tval, vrblvl)
        print('the values at t =', tval, ':\n', xval)
        residuals = evaluate_laurent_homotopy(lauhom, xval, tval, 'x', vrblvl)
        print('the residuals :', residuals)
        backward = sum([abs(nbr) for nbr in residuals])
        print('backward error :', backward, 'for t =', tval)
    return fail

def test_newton_steps(product=False, vrblvl=0):
    """
    Tests newton's method, either on the product (if True),
    or on the binomial,
    """
    if vrblvl > 0:
        print("in laurent.test_newton_steps ...")
    #if not product:
    #    return test_newton_binomial(vrblvl)
    #else:
    return test_newton_product(1)

def test_laurent(deg, vrblvl=0):
    """
    Tests operations in this module.
    The verbose level is defined by vrblvl.
    """
    if vrblvl > 0:
        print("in laurent.test_laurent ...")
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

def main(vrblvl=0):
    """
    Runs tests on the laurent module.
    """
    if vrblvl > 0:
        print("in laurent.main ...")
    deg = 2
    fail = test_laurent(deg, vrblvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__== '__main__':
    main()
