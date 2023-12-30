"""
Numerical Schubert calculus defines homotopies for enumerative geometry.
"""
from ctypes import c_int32, c_double, pointer, create_string_buffer
from random import uniform
from math import pi, sin, cos
from phcpy.version import get_phcfun
from phcpy.version import int4a2nbr, nbr2int4a, int4a2str, str2int4a
from phcpy.polynomials import get_double_polynomial

def pieri_root_count(mdim, pdim, qdeg, vrblvl=0):
    r"""
    Computes the number of *pdim*-plane producing maps of degree *qdeg*
    that meet *mdim*-planes at mdim*pdim + qdeg*(mdim+pdim) points.
    """
    if vrblvl > 0:
        print('in pieri_root_count, m :', mdim, end='')
        print(', p :', pdim, end='')
        print(', q :', qdeg)
    phc = get_phcfun()
    apars = int4a2nbr([mdim, pdim, qdeg], (vrblvl > 0))
    roco = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> pieri_root_count calls phc ...')
        print('apars =', nbr2int4a(apars))
    retval = phc(223, apars, roco, ccc, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
        print('the Pieri root count :', roco[0])
    return roco[0]

def pieri_localization_poset(mdim, pdim, qdeg, size=10240, vrblvl=0):
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
    phc = get_phcfun()
    apars = int4a2nbr([mdim, pdim, qdeg], (vrblvl > 0))
    poset = create_string_buffer(b"", 4*size)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> pieri_localization_poset calls phc ...')
        print('apars =', nbr2int4a(apars))
    retval = phc(224, apars, poset, ccc, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
    result = int4a2str(poset, (vrblvl > 0))
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
    phc = get_phcfun()
    apars = int4a2nbr([ndim, kdim, len(brackets)], (vrblvl > 0))
    brk = pointer(cds)
    roco = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> resolve_schubert_conditions calls phc ...')
        print('apars =', nbr2int4a(apars))
    retval = phc(228, apars, brk, roco, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
        print('the root count :', int(roco[0]))
    return int(roco[0])

def real_osculating_planes(mdim, pdim, qdeg, vrblvl=0):
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
    for i in range(len(pts)):
        cpts[i] = c_double(pts[i])
    phc = get_phcfun()
    apars = int4a2nbr([mdim, pdim, qdeg], (vrblvl > 0))
    bbb = pointer(c_int32(0))
    cff = pointer(cpts)
    vrb = c_int32(vrblvl)
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
        for i in range(0, mdim+pdim):
            row = []
            for j in range(0, mdim):
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
    from math import pi, sin, cos
    from random import uniform as u
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
    phc = get_phcfun()
    apars = int4a2nbr([mdim, pdim, qdeg], (vrblvl > 0))
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
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> make_pieri_system calls phc ...')
        print('apars =', nbr2int4a(apars))
    retval = phc(227, apars, bbb, cff, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
    result = []
    if(qdeg == 0):
        for i in range(1, mdim*pdim+1):
            result.append(get_double_polynomial(i))
    return result

def test_pieri(vrblvl=0):
    """
    Tests the Pieri root count.
    """
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
    lvl = 10
    print('making a real problem ...')
    planes = real_osculating_planes(2, 2, 0, lvl)
    pols = make_pieri_system(2, 2, 0, planes, is_real=True, vrblvl=lvl)
    for pol in pols:
        print(pol)
    print('making a complex problem ...')
    planes = random_complex_matrices(2*2, 4, 2)
    for (idx, plane) in enumerate(planes):
        print('plane', idx, ':')
        for row in plane:
            print(row)
    pols = make_pieri_system(2, 2, 0, planes, is_real=False, vrblvl=lvl)
    for pol in pols:
        print(pol)
    return 0

def main():
    """
    Runs some tests.
    """
    lvl = 10
    fail = test_pieri(lvl)
    fail = fail + test_littlewood_richardson_rule(lvl)
    fail = fail + test_pieri_problem(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=="__main__":
    main()
