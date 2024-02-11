"""
Numerical Schubert calculus defines homotopies for enumerative geometry.
"""
from ctypes import c_int32, c_double, pointer, create_string_buffer
from random import uniform
from math import pi, sin, cos
from phcpy.version import get_phcfun
from phcpy.version import int4a2nbr, nbr2int4a, int4a2str
from phcpy.polynomials import get_double_polynomial
from phcpy.polynomials import get_double_system
from phcpy.solutions import get_double_solutions, clear_double_solutions
from phcpy.solutions import verify

def pieri_root_count(mdim, pdim, qdeg=0, vrblvl=0):
    r"""
    Computes the number of *pdim*-plane producing maps of degree *qdeg*
    that meet *mdim*-planes at mdim*pdim + qdeg*(mdim+pdim) points.
    """
    if vrblvl > 0:
        print('in pieri_root_count, m :', mdim, end='')
        print(', p :', pdim, end='')
        print(', q :', qdeg)
    phc = get_phcfun(vrblvl-1)
    apars = int4a2nbr([mdim, pdim, qdeg], vrblvl=vrblvl-1)
    roco = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> pieri_root_count calls phc ...')
        print('apars =', nbr2int4a(apars))
    retval = phc(223, apars, roco, ccc, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
        print('the Pieri root count :', roco[0])
    return roco[0]

def pieri_localization_poset(mdim, pdim, qdeg=0, size=10240, vrblvl=0):
    r"""
    Returns the string representation of the localization poset
    used to compute the Pieri root count, for the number of *pdim*-plane
    producing maps of degree *qdeg* that meet *mdim*-planes at
    mdim*pdim + qdeg*(mdim+pdim) points.
    The size of the string representation is given on input in size.
    """
    if vrblvl > 0:
        print('in pieri_root_count, m :', mdim, end='')
        print(', p :', pdim, end='')
        print(', q :', qdeg)
    phc = get_phcfun(vrblvl-1)
    apars = int4a2nbr([mdim, pdim, qdeg], vrblvl=vrblvl-1)
    poset = create_string_buffer(b"", 4*size)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> pieri_localization_poset calls phc ...')
        print('apars =', nbr2int4a(apars))
    retval = phc(224, apars, poset, ccc, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
    allresult = int4a2str(poset, (vrblvl > 0))
    endidx = allresult.find('\0')
    result = allresult[:endidx]
    if vrblvl > 0:
        print('the localization poset :')
        print(result)
    return result

def resolve_schubert_conditions(ndim, kdim, brackets, vrblvl=0):
    r"""
    In n-dimensional space we consider k-dimensional planes,
    subject to intersection conditions represented by brackets.
    The *brackets* is a list of brackets.  A bracket is a list
    of as many natural numbers (in the range 1..*ndim*) as *kdim*.
    On return is the formal root count, which is sharp for general flags.
    and the coordinates of the flags, stored row wise in a list
    of real and imaginary parts.
    """
    if vrblvl > 0:
        print('in resolve_schubert_conditions, n :', ndim, end='')
        print(', k :', kdim, end='')
        print(', #brackets :', len(brackets))
    size = kdim*len(brackets)
    cds = (c_int32 * size)()
    if vrblvl > 0:
        print('conditions :', cds[:])
    idx = 0
    for bracket in brackets:
        for nbr in bracket:
            cds[idx] = c_int32(nbr)
            idx = idx + 1
    if vrblvl > 0:
        print('conditions :', cds[:])
    phc = get_phcfun(vrblvl-1)
    apars = int4a2nbr([ndim, kdim, len(brackets)], vrblvl=vrblvl-1)
    brk = pointer(cds)
    roco = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> resolve_schubert_conditions calls phc ...')
        print('apars =', nbr2int4a(apars))
    retval = phc(228, apars, brk, roco, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
        print('the root count :', int(roco[0]))
    return int(roco[0])

def real_osculating_planes(mdim, pdim, qdeg=0, vrblvl=0):
    """
    Returns m*p + qdeg*(m+p) real m-planes osculating
    a rational normal curve.
    The value of the verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in real_osculating planes, m :', mdim, end='')
        print(', p :', pdim, end='')
        print(', q :', qdeg)
    dim = mdim*pdim + qdeg*(mdim+pdim)
    pts = [uniform(-1, +1) for _ in range(dim)]
    if vrblvl > 0:
        print('the points :', pts)
    size = dim*mdim*(mdim+pdim)
    cpts = (c_double*size)()
    for (idx, point) in enumerate(pts):
        cpts[idx] = c_double(point)
    phc = get_phcfun(vrblvl-1)
    apars = int4a2nbr([mdim, pdim, qdeg], vrblvl=vrblvl-1)
    bbb = pointer(c_int32(0))
    cff = pointer(cpts)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> real_osculating_planes calls phc ...')
        print('apars =', nbr2int4a(apars))
    retval = phc(226, apars, bbb, cff, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
    idx = 0
    planes = []
    for k in range(0, dim):
        plane = []
        for _ in range(0, mdim+pdim):
            row = []
            for _ in range(0, mdim):
                row.append(cff[0][idx])
                idx = idx + 1
            plane.append(row)
        planes.append(plane)
        if vrblvl > 0:
            print('the coefficients of the plane', k, ':')
            for row in plane:
                print(row)
    return planes

def random_complex_matrix(nbrows, nbcols):
    r"""
    Returns a random *nbrows*-by-*nbcols* matrix
    with randomly generated complex coefficients
    on the unit circle, as a list of *rows*.
    """
    result = []
    for _ in range(0, nbrows):
        angles = [uniform(0, 2*pi) for _ in range(nbcols)]
        cols = [complex(cos(a), sin(a)) for a in angles]
        result.append(cols)
    return result

def random_complex_matrices(nbr, nbrows, nbcols):
    r"""
    Returns a list of matrix of length *nbr*,
    all of dimension *nbrows* by *nbcols*.
    """
    result = []
    for _ in range(nbr):
        result.append(random_complex_matrix(nbrows, nbcols))
    return result

def make_pieri_system(mdim, pdim, qdeg, planes, is_real=False, vrblvl=0):
    """
    Makes the polynomial system defined by the mdim-planes
    in the list planes.
    """
    if vrblvl > 0:
        print('in make_pieri_system, m :', mdim, end='')
        print(', p :', pdim, end='')
        print(', q :', qdeg)
    phc = get_phcfun(vrblvl-1)
    apars = int4a2nbr([mdim, pdim, qdeg], vrblvl=vrblvl-1)
    bbb = pointer(c_int32(0))
    dim = mdim*pdim + qdeg*(mdim+pdim)
    size = 2*dim*mdim*(mdim+pdim)
    cpts = (c_double*size)()
    idx = 0
    for plane in planes:
        for row in plane:
            for nbr in row:
                if is_real:
                    cpts[idx] = c_double(nbr)
                    idx = idx + 2
                else:
                    cpts[idx] = c_double(nbr.real)
                    idx = idx + 1
                    cpts[idx] = c_double(nbr.imag)
                    idx = idx + 1
    cff = pointer(cpts)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> make_pieri_system calls phc ...')
        print('apars =', nbr2int4a(apars))
    retval = phc(227, apars, bbb, cff, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
    result = []
    if qdeg == 0:
        for i in range(1, mdim*pdim+1):
            result.append(get_double_polynomial(i, vrblvl))
    return result

def run_pieri_homotopies(mdim, pdim, qdeg, planes, vrblvl=0, *pts):
    r"""
    Computes the number of *pdim*-plane producing maps of degree *qdeg*
    that meet *mdim*-planes at mdim*pdim + qdeq*(mdim+pdim) points.
    For *qdeg* = 0, there are no interpolation points in *pts*.
    The planes and pts are assumed to be complex.
    """
    if vrblvl > 0:
        print('in run_pieri_homotopies, m :', mdim, end='')
        print(', p :', pdim, end='')
        print(', q :', qdeg)
    roco = pieri_root_count(mdim, pdim, qdeg, vrblvl)
    if vrblvl > 0:
        print('Pieri root count for', (mdim, pdim, qdeg), 'is', roco)
    clear_double_solutions(vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    apars = int4a2nbr([mdim, pdim, qdeg], vrblvl=vrblvl-1)
    broco = pointer(c_int32(0))
    dim = mdim*pdim + qdeg*(mdim+pdim)
    if qdeg == 0:
        size = 2*dim*mdim*(mdim+pdim)
    else:
        size = 2*len(pts[0]) + 2*dim*mdim*(mdim+pdim)
    cpts = (c_double*size)()
    idx = 0
    if qdeg > 0:            # points are in rows of a matrix
        for row in pts[0]:
            for cffpts in row:
                cpts[idx] = c_double(cffpts.real)
                idx = idx + 1
                cpts[idx] = c_double(cffpts.imag)
                idx = idx + 1
    for plane in planes:
        for row in plane:
            for nbr in row:
                cpts[idx] = c_double(nbr.real)
                idx = idx + 1
                cpts[idx] = c_double(nbr.imag)
                idx = idx + 1
    cff = pointer(cpts)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> make_pieri_system calls phc ...')
        print('apars =', nbr2int4a(apars))
    retval = phc(225, apars, broco, cff, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
        print('the root count :', broco[0])
    if vrblvl > 0:
        print('making the system ...')
    pols = []
    if qdeg == 0:
        for i in range(1, mdim*pdim+1):
            pols.append(get_double_polynomial(i, vrblvl-1))
    else:
        for i in range(1, mdim*pdim+qdeg*(mdim+pdim)+1):
            pols.append(get_double_polynomial(i, vrblvl-1))
    if vrblvl > 0:
        print('the system :')
        for pol in pols:
            print(pol)
        print('root count :', roco)
    sols = get_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('the solutions :')
        for sol in sols:
            print(sol)
    return (pols, sols)

def double_littlewood_richardson_homotopies(ndim, kdim, brackets, \
    verbose=True, vrfcnd=False, minrep=True, tosqr=False, outputfilename='',
    vrblvl=0):
    r"""
    In n-dimensional space we consider k-dimensional planes,
    subject to intersection conditions represented by brackets.
    The parameters *ndim* and *kdim* give values for n and k respectively.
    The parameter brackets is a list of brackets.  A bracket is a list
    of as many natural numbers (in the range 1..*ndim*) as *kdim*.
    The Littlewood-Richardson homotopies compute k-planes that
    meet the flags at spaces of dimensions prescribed by the brackets,
    in standard double precision.  Four options are passed as Booleans:

    *verbose*: for adding extra output during computations,

    *vrfcnd*: for extra diagnostic verification of Schubert conditions,

    *minrep*: for a minimial representation of the problem formulation,

    *tosqr*: to square the overdetermined systems.

    On return is a 4-tuple.  The first item of the tuple is the
    formal root count, sharp for general flags, then as second
    item the coordinates of the flags.  The coordinates of the
    flags are stored row wise in a list of real and imaginary parts.
    The third and fourth item of the tuple on return are respectively
    the polynomial system that has been solved and its solutions.
    The length of the list of solution should match the root count.
    """
    if vrblvl > 0:
        print('in double_littlewood_richardson_homotopies ...')
        print('ndim :', ndim, '  kdim :', kdim)
        print('the brackets :')
        for bracket in brackets:
            print(bracket)
    clear_double_solutions(vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    adims = (c_int32 * 8)()
    adims[0] = c_int32(ndim)
    adims[1] = c_int32(kdim)
    adims[2] = c_int32(len(brackets))
    adims[3] = c_int32(int(verbose))
    adims[4] = c_int32(int(vrfcnd))
    adims[5] = c_int32(len(outputfilename))
    adims[6] = c_int32(int(minrep))
    adims[7] = c_int32(int(tosqr))
    pdims = pointer(adims)
    nbc = len(brackets)*kdim
    bcnd = (c_int32 * nbc)()
    idx = 0
    for bracket in brackets:
        for num in bracket:
            bcnd[idx] = c_int32(num)
            idx = idx + 1
    pcnd = pointer(bcnd)
    size = 4*(len(brackets) - 2)*ndim*ndim+1
    name = (c_double * size)()
    for (idx, letter) in enumerate(outputfilename):
        name[idx] = letter
    pname = pointer(name)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_littlewood_richardson_homotopies calls phc ...')
    retval = phc(229, pdims, pcnd, pname, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
    roco = int(pname[0][0])
    if vrblvl > 0:
        print('the root count :', roco)
    allcffs = pname[0][1:size]
    flgs = []
    for idx in range(size, 2):
        realnbr = allcffs[0][idx]
        imagnbr = allcffs[0][idx+1]
        flgs.append(complex(realnbr, imagnbr))
    fsys = get_double_system(vrblvl-1)
    sols = get_double_solutions(vrblvl-1)
    return (roco, flgs, fsys, sols)

def test_pieri_count(vrblvl=0):
    """
    Tests the Pieri root count.
    """
    if vrblvl > 0:
        print('in test_pieri ...')
    roco = pieri_root_count(2, 2, 1, vrblvl)
    poset = pieri_localization_poset(2, 2, 1, vrblvl=vrblvl)
    if vrblvl > 0:
        print('The pieri localization poset :')
        print(poset)
    if roco == 8:
        if vrblvl > 0:
            print('Pieri root count', roco, 'is correct.')
        return 0
    if vrblvl > 0:
        print('Pieri root count', roco, 'is wrong.')
    return 1

def test_littlewood_richardson_rule(vrblvl=0):
    """
    Tests the Littlewood-Richardson rule.
    """
    if vrblvl > 0:
        print('in test_littlewood_richardson_rule ...')
    brk = [[2, 4, 6], [2, 4, 6], [2, 4, 6]]
    roco = resolve_schubert_conditions(6, 3, brk, vrblvl)
    if roco == 2:
        if vrblvl > 0:
            print('Littlewood-Richardson root count :', roco, 'is correct.')
        return 0
    if vrblvl > 0:
        print('Littlewood-Richardson root count :', roco, 'is wrong.')
    return 1

def test_pieri_problem(vrblvl=0):
    """
    Tests the making of a Pieri problem.
    """
    if vrblvl > 0:
        print('in test_pieri_problem ...')
    planes = real_osculating_planes(2, 2, 0, vrblvl)
    pols = make_pieri_system(2, 2, 0, planes, is_real=True, vrblvl=vrblvl)
    if vrblvl > 0:
        for pol in pols:
            print(pol)
    if vrblvl > 0:
        print('making a complex problem ...')
    planes = random_complex_matrices(2*2, 4, 2)
    if vrblvl > 0:
        for (idx, plane) in enumerate(planes):
            print('plane', idx, ':')
            for row in plane:
                print(row)
    pols = make_pieri_system(2, 2, 0, planes, is_real=False, vrblvl=vrblvl)
    if vrblvl > 0:
        for pol in pols:
            print(pol)
    return 0

def test_pieri_homotopies(vrblvl=0):
    """
    Tests running the Pieri homotopies on the four lines problem.
    """
    if vrblvl > 0:
        print('in test_pieri_homotopies ...')
    planes = random_complex_matrices(2*2, 4, 2)
    if vrblvl > 0:
        for (idx, plane) in enumerate(planes):
            print('plane', idx, ':')
            for row in plane:
                print(row)
    pols, sols = run_pieri_homotopies(2, 2, 0, planes, vrblvl=vrblvl)
    if vrblvl > 0:
        for pol in pols:
            print(pol)
        for sol in sols:
            print(sol)
    err = verify(pols, sols, vrblvl-1)
    if vrblvl > 0:
        print('the error sum :', err)
    if len(sols) == 2 and err < 1.0e-10:
        if vrblvl > 0:
            print('Found 2 solutions and error is okay.')
        return 0
    if len(sols) != 2:
        if vrblvl > 0:
            print('Number of solutions is not 2 :', len(sols))
        return 1
    if err >= 1.0e-10:
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def test_littlewood_richardson_homotopies(vrblvl=0):
    """
    Performs a test on the Littlewood-Richardson homotopies.
    """
    if vrblvl > 0:
        print('in test_littlewood_richardson_homotopies ...')
    brk = [[2, 4, 6], [2, 4, 6], [2, 4, 6]]
    (roco, flags, fsys, sols) \
        = double_littlewood_richardson_homotopies\
             (6, 3, brk, vrblvl=vrblvl)
    if vrblvl > 0:
        print('the root count :', roco)
    if vrblvl > 0:
        print('the flags :', flags)
    if vrblvl > 0:
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(fsys, sols, vrblvl-1)
    if vrblvl > 0:
        print('the error sum :', err)
    if len(sols) == 2 and err < 1.0e-10:
        if vrblvl > 0:
            print('Found 2 solutions and error is okay.')
        return 0
    if len(sols) != 2:
        if vrblvl > 0:
            print('Number of solutions is not 2 :', len(sols))
        return 1
    if err >= 1.0e-10:
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def test_pieri_curves(vrblvl=0):
    """
    Computes all line producing curves that meet given lines in 3-space,
    at given interpolation points, computed with Pieri homotopies.
    """
    if vrblvl > 0:
        print('in test_pieri_curves ...')
    (mdim, pdim, qdeg) = (2, 2, 1)
    dim = mdim*pdim + qdeg*(mdim+pdim)
    roco = pieri_root_count(mdim, pdim, qdeg)
    if vrblvl > 0:
        print('the root count :', roco)
    pls = [random_complex_matrix(mdim+pdim, mdim) for _ in range(dim)]
    points = random_complex_matrix(dim, 1)
    (fsys, fsols) = run_pieri_homotopies(mdim, pdim, qdeg, pls, vrblvl, points)
    if vrblvl > 0:
        print('number of solutions :', len(fsols))
    fail = int(roco != len(fsols))
    err = verify(fsys, fsols, vrblvl-1)
    if vrblvl > 0:
        print('the error :', err)
    fail = fail + int(err > 1.0e-8)
    return fail

def main():
    """
    Runs some tests.
    """
    lvl = 1
    fail = test_pieri_count(lvl)
    fail = fail + test_littlewood_richardson_rule(lvl)
    fail = fail + test_pieri_problem(lvl)
    fail = fail + test_pieri_homotopies(lvl)
    fail = fail + test_littlewood_richardson_homotopies(lvl)
    fail = fail + test_pieri_curves(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=="__main__":
    main()
