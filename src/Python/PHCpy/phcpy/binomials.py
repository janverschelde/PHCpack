"""
A binomial system is a polynomial system with exactly two monomials 
with nonzero coefficient in every equation.  Binomials systems can be
solved much more efficiently than other polynomial systems.
The positive dimensional solution sets of binomial systems
are monomial maps.
"""
from ctypes import c_int32, c_double, pointer
from phcpy.version import get_phcfun
from phcpy.dimension import get_double_laurent_dimension
from phcpy.polynomials import string_of_symbols
from phcpy.polynomials import set_double_laurent_system
from phcpy.polynomials import get_double_number_laurent_terms

def is_binomial_system(vrblvl=0):
    r"""
    Returns True if the system stored as a Laurent system
    in double precision is a binomial system, return False otherwise.
    The verbose level is given by the value of vrblvl.
    """
    if vrblvl > 0:
        print('in is_binomial_system ...')
    nbequ = get_double_laurent_dimension(vrblvl-1)
    if vrblvl > 0:
        print('number of Laurent polynomials :', nbequ)
    for i in range(1, nbequ+1):
        nbterms = get_double_number_laurent_terms(i, vrblvl-1)
        if vrblvl > 0:
            print('  -> number of terms in polynomial', i, ':', nbterms)
        if nbterms != 2:
            if vrblvl > 0:
                print('  the system is not a binomial system')
            return False
    if vrblvl > 0:
        print('  the system is a binomial system')
    return True

def double_map_top_dimension(vrblvl=0):
    """
    Returns the top dimension of the solution maps,
    in double precision.
    """
    if vrblvl > 0:
        print('in double_map_top_dimension ...')
    phc = get_phcfun(vrblvl-1)
    dim = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_map_top_dimension calls phc', end='')
    retval = phc(433, dim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('-> top dimension :', dim[0])
    return dim[0]

def number_double_maps(dim, vrblvl=0):
    """
    Returns the number of maps of dimension dim,
    in double precision.
    """
    if vrblvl > 0:
        print('in number_double_maps, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    nbr = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> number_double_maps calls phc', end='')
    retval = phc(434, adim, nbr, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('-> number of maps :', nbr[0])
    return nbr[0]

def degree_double_map(dim, idx, vrblvl=0):
    """
    Returns the degree of map with index idx (starting at one),
    of dimension dim, in double precision.
    """
    if vrblvl > 0:
        print('in degree_double_map, dim :', dim, end='')
        print(', idx :', idx)
    phc = get_phcfun(vrblvl-1)
    dims = (c_int32 * 2)()
    dims[0] = dim
    dims[1] = idx
    adim = pointer(dims)
    deg = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> degree_double_map calls phc', end='')
    retval = phc(435, adim, deg, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('-> degree of map :', deg[0])
    return deg[0]

def coefficients_double_map(dim, idx, nvr, vrblvl=0):
    """
    Returns the complex coefficients of a map of dimension dim,
    of index idx, in nvr variables, in double precision.
    """
    if vrblvl > 0:
        print('in coefficients_double_map, dim :', dim, end='')
        print(', idx :', idx, ', nvr :', nvr)
    phc = get_phcfun(vrblvl-1)
    dims = (c_int32 * 3)()
    dims[0] = dim
    dims[1] = idx
    dims[2] = nvr
    adim = pointer(dims)
    bbb = pointer(c_int32(0))
    # complex coefficients
    ccnvr = 2*nvr
    cffs = (c_double * ccnvr)()
    ccff = pointer(cffs)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> coefficients_double_map calls phc', end='')
    retval = phc(436, adim, bbb, ccff, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    vals = ccff[:2*nvr]
    result = []
    for cidx in range(0, 2*nvr, 2):
        coefficient = complex(vals[0][cidx], vals[0][cidx+1])
        result.append(coefficient)
    if vrblvl > 0:
        print('the coefficients :', result)
    return result

def exponents_double_map(dim, idx, nvr, vrblvl=0):
    """
    Returns the exponents of a map of dimension dim,
    of index idx, in nvr variables, in double precision.
    """
    if vrblvl > 0:
        print('in exponents_double_map, dim :', dim, end='')
        print(', idx :', idx, ', nvr :', nvr)
    phc = get_phcfun(vrblvl-1)
    dims = (c_int32 * 3)()
    dims[0] = dim
    dims[1] = idx
    dims[2] = nvr
    adim = pointer(dims)
    nbexp = dim*nvr
    exp = (c_int32 * nbexp)()
    bexp = pointer(exp)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> exponents_double_map calls phc', end='')
    retval = phc(437, adim, bexp, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    vals = bexp[:nvr]
    result = []
    for cidx in range(0, nbexp, dim):
        expdim = []
        for cnt in range(dim):
            expdim.append(vals[0][cidx+cnt])
        result.append(expdim)
    if vrblvl > 0:
        print('the exponents :', result)
    return result

def double_string_map(dim, nvr, cff, exp, vrblvl=0):
    """
    Returns the string representation of a map of dimension dim,
    in nvr variables, with double complex coefficients in cff,
    and exponents in exp.
    """
    if vrblvl > 0:
        print('in double_string_map, dim :', dim, end='')
        print(', nvr :', nvr)
        print('the coefficients :', cff)
        print('the exponents :', exp)
    smb = string_of_symbols(100, vrblvl)
    if vrblvl > 0:
        print('the list of symbols :', smb)
    result = []
    for idx in range(nvr):
        strvar = smb[idx]
        if cff[idx] != complex(0.0, 0.0):
            strvar = strvar + ' - ' + str(cff[idx])
            for (expidx, pwr) in enumerate(exp[idx]):
                if pwr != 0:
                    strvar = strvar + '*' + 't' + str(expidx+1)
                    if pwr != 1:
                        strvar = strvar + '**' + str(pwr)
        result.append(strvar)
    return result

def double_solve(nvr, pols, puretopdim=False, vrblvl=0):
    """
    If the polynomial system in pols in number of variables nvr
    is a binomial system, then the solution maps are returned.
    When the flag puretopdim is True, then the solver will only
    look for the solution sets of the expected top dimension,
    which will lead to a faster execution.
    """
    if vrblvl > 0:
        print('in double_solve, nvr :', nvr)
        print('the polynomials :')
        for pol in pols:
            print(pol)
    set_double_laurent_system(nvr, pols, vrblvl-1)
    result = []
    isbin = is_binomial_system(vrblvl)
    if not isbin:
        if vrblvl > 0:
            print('The system is not binomial.')
    else:
        phc = get_phcfun(vrblvl-1)
        adim = pointer(c_int32(int(puretopdim)))
        bbb = pointer(c_int32(0))
        ccc = pointer(c_double(0.0))
        vrb = c_int32(vrblvl-1)
        if vrblvl > 0:
            print('-> double_solve calls phc', end='')
        retval = phc(430, adim, bbb, ccc, vrb)
        if vrblvl > 0:
            print(', return value :', retval)
        topdim = double_map_top_dimension(vrblvl)
        if vrblvl > 0:
            print('the top dimension :', topdim)
        for dim in range(topdim, -1, -1):
            nbr = number_double_maps(dim, vrblvl)
            if vrblvl > 0:
                print('at dimension', dim, 'number of maps :', nbr)
            for idx in range(nbr):
                deg = degree_double_map(dim, idx+1, vrblvl)
                if vrblvl > 0:
                    print('map', idx+1, 'has degree', deg)
                cff = coefficients_double_map(dim, idx+1, nvr, vrblvl)
                if vrblvl > 0:
                    print('the coefficients :')
                    for coeff in cff:
                        print(coeff)
                exp = exponents_double_map(dim, idx+1, nvr, vrblvl)
                if vrblvl > 0:
                    print('the exponents :', exp)
                strmap = double_string_map(dim, nvr, cff, exp, vrblvl)
                if vrblvl > 0:
                    print('the string representation :')
                    for component in strmap:
                        print(component)
                result.append(strmap)
    return result

def test_is_binomial_system(vrblvl=0):
    """
    Tests on the binomial system test.
    """
    if vrblvl > 0:
        print('in test_is_binomial_system ...')
    twisted_cubic = ['x^2 - y;', 'x^3 - z;']
    if vrblvl > 0:
        print('the twisted cubic :')
        for pol in twisted_cubic:
            print(pol)
    set_double_laurent_system(3, twisted_cubic, vrblvl-1)
    isbin = is_binomial_system(vrblvl)
    fail = int(not isbin)
    return fail

def test_double_solve(vrblvl=0):
    """
    Tests on solving a binomial system.
    """
    if vrblvl > 0:
        print('in test_double_solve ...')
    twisted_cubic = ['x**2*y - z*x;', 'x**2*z - y**2*x;']
    if vrblvl > 0:
        print('the twisted cubic :')
        for pol in twisted_cubic:
            print(pol)
    sols = double_solve(3, twisted_cubic, vrblvl=vrblvl)
    if vrblvl > 0:
        for (idxmap, solmap) in enumerate(sols):
            print('Solution map', idxmap+1, ':')
            print(solmap)
    topdim = double_map_top_dimension(vrblvl)
    fail = int(topdim != 2)
    dims = {}
    for dim in range(topdim, -1, -1):
        dims[dim] = number_double_maps(dim, vrblvl)
    if vrblvl > 0:
        print('the dimensions :')
        print(dims)
    fail = fail + int(dims[2] != 1) # one two dimensional set
    fail = fail + int(dims[1] != 2) # two one dimensional sets
    fail = fail + int(dims[0] != 0) # no isolated solutions
    if vrblvl > 0:
        print('the degrees :')
    degs = []
    for dim in range(topdim, -1, -1):
        for nbr in range(dims[dim]):
            deg = degree_double_map(dim, nbr+1, vrblvl)
            if vrblvl > 0:
                print('dimension :', dim, end='')
                print(', degree :', deg)
            degs.append(deg)
    if vrblvl > 0:
        print('degrees :', degs)
    fail = fail + int(degs[0] != 1) # two dimensional set is linear
    fail = fail + int(degs[1] != 3) # twisted cubic
    fail = fail + int(degs[2] != 1) # linear one dimensional set
    return fail

def main():
    """
    Runs some tests.
    """
    lvl = 1
    fail = test_is_binomial_system(lvl)
    fail = fail + test_double_solve(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=="__main__":
    main()
