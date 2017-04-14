"""
PHCpack offers numerical Schubert calculus, exported here.
"""

def prompt_for_dimensions():
    """
    Returns the triplet (m,p,q),
    where m is the dimension of the input planes,
    p is the dimension of the output planes, and
    q is the degree of the maps.
    """
    mdim = input('give the dimension of the input planes : ')
    pdim = input('give the dimension of the output planes : ')
    qdeg = input('give the degree of the solution maps : ')
    return (int(mdim), int(pdim), int(qdeg))

def pieri_root_count(mdim, pdim, qdeg, verbose=True):
    r"""
    Computes the number of *pdim*-plane producing maps of
    degree qdeg that meet *mdim*-planes at mdim*pdim + qdeg*(mdim+pdim) points.
    """
    from phcpy.phcpy2c3 import py2c_schubert_pieri_count
    from phcpy.phcpy2c3 import py2c_schubert_localization_poset
    root_count = py2c_schubert_pieri_count(mdim, pdim, qdeg)
    if verbose:
        print('Pieri root count for', (mdim, pdim, qdeg), 'is', root_count)
        poset = py2c_schubert_localization_poset(mdim, pdim, qdeg)
        print('the localization poset :')
        print(poset)
    return root_count

def resolve_schubert_conditions(ndim, kdim, brackets, verbose=True):
    r"""
    In n-dimensional space we consider k-dimensional planes,
    subject to intersection conditions represented by brackets.
    The *brackets* is a list of brackets.  A bracket is a list
    of as many natural numbers (in the range 1..*ndim*) as *kdim*.
    On return is the formal root count, which is sharp for general flags.
    and the coordinates of the flags, stored row wise in a list
    of real and imaginary parts.
    """
    from phcpy.phcpy2c3 import py2c_schubert_resolve_conditions as resolve
    nbc = len(brackets)
    cds = ''
    for bracket in brackets:
        for num in bracket:
            cds = cds + ' ' + str(num)
    # print 'the condition string :', cds
    roco = resolve(ndim, kdim, nbc, len(cds), cds, int(verbose))
    return roco

def standard_littlewood_richardson_homotopies(ndim, kdim, brackets, \
    verbose=True, vrfcnd=False, minrep=True, tosqr=False, outputfilename=''):
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
    from phcpy.phcpy2c3 import py2c_solcon_clear_standard_solutions
    from phcpy.phcpy2c3 \
       import py2c_schubert_standard_littlewood_richardson_homotopies \
       as stlrhom
    from phcpy.interface import load_standard_solutions, load_standard_system
    py2c_solcon_clear_standard_solutions()
    nbc = len(brackets)
    cds = ''
    for bracket in brackets:
        for num in bracket:
            cds = cds + ' ' + str(num)
    # print 'the condition string :', cds
    (roco, sflags) = stlrhom(ndim, kdim, nbc, len(cds), cds, \
        int(verbose), int(vrfcnd), int(minrep), int(tosqr), \
        len(outputfilename), outputfilename)
    rflags = eval(sflags)
    flgs = []
    for k in range(len(rflags)//2):
        flgs.append(complex(rflags[2*k], rflags[2*k+1]))
    fsys = load_standard_system()
    sols = load_standard_solutions()
    return (roco, flgs, fsys, sols)

def dobldobl_littlewood_richardson_homotopies(ndim, kdim, brackets, \
    verbose=True, vrfcnd=False, minrep=True, tosqr=False, outputfilename=''):
    r"""
    In n-dimensional space we consider k-dimensional planes,
    subject to intersection conditions represented by brackets.
    The parameters *ndim* and *kdim* give values for n and k respectively.
    The parameter *brackets* is a list of brackets.  A bracket is a list
    of as many natural numbers (in the range 1..*ndim*) as *kdim*.
    The Littlewood-Richardson homotopies compute k-planes that
    meet the flags at spaces of dimensions prescribed by the brackets,
    in double double precision.  Four options are passed as Booleans:

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
    from phcpy.phcpy2c3 import py2c_solcon_clear_dobldobl_solutions
    from phcpy.phcpy2c3 \
       import py2c_schubert_dobldobl_littlewood_richardson_homotopies \
       as ddlrhom
    from phcpy.interface import load_dobldobl_solutions, load_dobldobl_system
    py2c_solcon_clear_dobldobl_solutions()
    nbc = len(brackets)
    cds = ''
    for bracket in brackets:
        for num in bracket:
            cds = cds + ' ' + str(num)
    # print 'the condition string :', cds
    (roco, sflags) = ddlrhom(ndim, kdim, nbc, len(cds), cds, \
        int(verbose), int(vrfcnd), int(minrep), int(tosqr), \
        len(outputfilename), outputfilename)
    rflags = eval(sflags)
    flgs = []
    for k in range(len(rflags)//4):
        flgs.append(complex(rflags[2*k], rflags[2*k+2]))
    fsys = load_dobldobl_system()
    sols = load_dobldobl_solutions()
    return (roco, flgs, fsys, sols)

def quaddobl_littlewood_richardson_homotopies(ndim, kdim, brackets, \
    verbose=True, vrfcnd=False, minrep=True, tosqr=False, outputfilename=''):
    r"""
    In n-dimensional space we consider k-dimensional planes,
    subject to intersection conditions represented by brackets.
    The parameters *ndim* and *kdim* give values for n and k respectively.
    The parameter *brackets* is a list of brackets.  A bracket is a list
    of as many natural numbers (in the range 1..ndim) as kdim.
    The Littlewood-Richardson homotopies compute k-planes that
    meet the flags at spaces of dimensions prescribed by the brackets,
    in quad double precision.  Four options are passed as Booleans:

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
    from phcpy.phcpy2c3 import py2c_solcon_clear_quaddobl_solutions
    from phcpy.phcpy2c3 \
       import py2c_schubert_quaddobl_littlewood_richardson_homotopies \
       as qdlrhom
    from phcpy.interface import load_quaddobl_solutions, load_quaddobl_system
    py2c_solcon_clear_quaddobl_solutions()
    nbc = len(brackets)
    cds = ''
    for bracket in brackets:
        for num in bracket:
            cds = cds + ' ' + str(num)
    # print 'the condition string :', cds
    (roco, sflags) = qdlrhom(ndim, kdim, nbc, len(cds), cds, \
        int(verbose), int(vrfcnd), int(minrep), int(tosqr), \
        len(outputfilename), outputfilename)
    rflags = eval(sflags)
    flgs = []
    for k in range(len(rflags)//8):
        flgs.append(complex(rflags[2*k], rflags[2*k+4]))
    fsys = load_quaddobl_system()
    sols = load_quaddobl_solutions()
    return (roco, flgs, fsys, sols)

def littlewood_richardson_homotopies(ndim, kdim, brackets, \
    verbose=True, vrfcnd=False, minrep=True, tosqr=False, \
    precision='d', outputfilename=''):
    """
    In n-dimensional space we consider k-dimensional planes,
    subject to intersection conditions represented by brackets.
    The parameters *ndim* and *kdim* give values for n and k respectively.
    The parameter *brackets* is a list of brackets.  A bracket is a list
    of as many natural numbers (in the range 1..*ndim*) as *kdim*.
    The Littlewood-Richardson homotopies compute k-planes that
    meet the flags at spaces of dimensions prescribed by the brackets.
    Four options are passed as Booleans:

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
    if(precision == 'd'):
        return standard_littlewood_richardson_homotopies(ndim, kdim, \
                  brackets, verbose, vrfcnd, minrep, tosqr, outputfilename)
    elif(precision == 'dd'):
        return dobldobl_littlewood_richardson_homotopies(ndim, kdim, \
                  brackets, verbose, vrfcnd, minrep, tosqr, outputfilename)
    elif(precision == 'qd'):
        return quaddobl_littlewood_richardson_homotopies(ndim, kdim, \
                  brackets, verbose, vrfcnd, minrep, tosqr, outputfilename)
    else:
        print('wrong level of precision, use d, dd, or qd')

def random_complex_matrix(nbrows, nbcols):
    r"""
    Returns a random *nbrows*-by-*nbcols* matrix
    with randomly generated complex coefficients
    on the unit circle, as a list of *rows*.
    """
    from math import pi, sin, cos
    from random import uniform as u
    result = []
    for i in range(0, nbrows):
        angles = [u(0, 2*pi) for _ in range(nbcols)]
        cols = [complex(cos(a), sin(a)) for a in angles]
        result.append(cols)
    return result

def random_complex_matrices(nbr, nbrows, nbcols):
    r"""
    Returns a list of matrix of length *nbr*,
    all of dimension *nbrows* by *nbcols*.
    """
    result = []
    for i in range(nbr):
        result.append(random_complex_matrix(nbrows, nbcols))
    return result

def planes_to_string(planes):
    r"""
    Returns one long string with all numbers
    in *planes*, a list of lists of rows.
    The numbers are the real and imaginary parts,
    separated by space.
    """
    result = ""
    for plane in planes:
        for row in plane:
            for cff in row:
                (cffre, cffim) = (cff.real, cff.imag)
                result = result + ' ' + ('%.17e' % cffre)
                result = result + ' ' + ('%.17e' % cffim)
    return result

def points_to_string(pts):
    r"""
    Returns one long string with all numbers in *pts*,
    as sequences of real and imaginary parts,
    every number is separated by one space.
    """
    result = ""
    for row in pts:
        for cff in row:
            (xre, yre) = (cff.real, cff.imag)
            result = result + ' ' + ('%.17e' % xre)
            result = result + ' ' + ('%.17e' % yre)
    return result

def run_pieri_homotopies(mdim, pdim, qdeg, planes, *pts, **opt):
    r"""
    Computes the number of *pdim*-plane producing maps of degree *qdeg*
    that meet *mdim*-planes at mdim*pdim + qdeq*(mdim+pdim) points.
    For *qdeg* = 0, there are no interpolation points in *pts*.
    """
    from phcpy.phcpy2c3 import py2c_schubert_pieri_count
    from phcpy.phcpy2c3 import py2c_schubert_pieri_homotopies
    from phcpy.phcpy2c3 import py2c_syscon_load_standard_polynomial
    from phcpy.phcpy2c3 import py2c_solcon_clear_standard_solutions
    from phcpy.phcpy2c3 import py2c_solcon_number_of_standard_solutions
    from phcpy.phcpy2c3 import py2c_solcon_length_standard_solution_string
    from phcpy.phcpy2c3 import py2c_solcon_write_standard_solution_string
    if 'verbose' not in opt:
        verbose = True        # by default, the function is verbose
    else:
        verbose = opt['verbose']
    root_count = py2c_schubert_pieri_count(mdim, pdim, qdeg)
    if verbose:
        print('Pieri root count for', (mdim, pdim, qdeg), 'is', root_count)
    strplanes = planes_to_string(planes)
    # print 'length of input data :', len(strplanes)
    if(qdeg > 0):
        strpts = points_to_string(pts[0])
        # print 'the interpolation points :', strpts
    else:
        strpts = ''
    # print 'calling py2c_pieri_homotopies ...'
    py2c_solcon_clear_standard_solutions()
    if verbose:
        print('passing %d characters for (m, p, q) = (%d, %d, %d)' \
            % (len(strplanes), mdim, pdim, qdeg))
    py2c_schubert_pieri_homotopies(mdim, pdim, qdeg, \
        len(strplanes), strplanes, strpts)
    # print 'making the system ...'
    pols = []
    if(qdeg == 0):
        for i in range(1, mdim*pdim+1):
            pols.append(py2c_syscon_load_standard_polynomial(i))
    else:
        for i in range(1, mdim*pdim+qdeg*(mdim+pdim)+1):
            pols.append(py2c_syscon_load_standard_polynomial(i))
    if verbose:
        print('the system :')
        for poly in pols:
            print(poly)
        print('root count :', root_count)
    nbsols = py2c_solcon_number_of_standard_solutions()
    sols = []
    for k in range(1, nbsols+1):
        lns = py2c_solcon_length_standard_solution_string(k)
        sol = py2c_solcon_write_standard_solution_string(k, lns)
        sols.append(sol)
    if verbose:
        print('the solutions :')
        for solution in sols:
            print(solution)
    return (pols, sols)

def verify(pols, sols):
    r"""
    Verifies whether the solutions in *sols*
    satisfy the polynomials of the system in *pols*.
    """
    from phcpy.solutions import strsol2dict, evaluate
    dictsols = [strsol2dict(sol) for sol in sols]
    checksum = 0
    for sol in dictsols:
        sumeval = sum(evaluate(pols, sol))
        print(sumeval)
        checksum = checksum + sumeval
    print('the total check sum :', checksum)

def real_osculating_planes(mdim, pdim, qdeg):
    """
    Returns m*p + qdeg*(m+p) real m-planes osculating
    a rational normal curves.
    """
    from phcpy.phcpy2c3 import py2c_schubert_osculating_planes
    dim = mdim*pdim + qdeg*(mdim+pdim)
    from random import uniform as u
    pts = ""
    for k in range(dim):
        cff = '%.17lf' % u(-1, +1)
        pts = pts + ' ' + cff
    # print 'the points :', pts
    osc = py2c_schubert_osculating_planes(mdim, pdim, qdeg, len(pts), pts)
    # print 'the coefficients of the planes :'
    # print osc
    items = osc.split(' ')
    ind = 0
    planes = []
    for k in range(0, dim):
        plane = []
        for i in range(0, mdim+pdim):
            row = []
            for j in range(0, mdim):
                row.append(eval(items[ind]))
                ind = ind + 1
            plane.append(row)
        planes.append(plane)
    return planes

def make_pieri_system(mdim, pdim, qdeg, planes, is_real=False):
    """
    Makes the polynomial system defined by the mdim-planes
    in the list planes.
    """
    from phcpy.phcpy2c3 import py2c_schubert_pieri_system
    from phcpy.phcpy2c3 import py2c_syscon_load_standard_polynomial
    strplanes = planes_to_string(planes)
    # print 'the string of planes :', strplanes
    if is_real:
        py2c_schubert_pieri_system(mdim, pdim, qdeg, \
            len(strplanes), strplanes, 1)
    else:
        py2c_schubert_pieri_system(mdim, pdim, qdeg, \
            len(strplanes), strplanes, 0)
    result = []
    if(qdeg == 0):
        for i in range(1, mdim*pdim+1):
            result.append(py2c_syscon_load_standard_polynomial(i))
    return result

def cheater(mdim, pdim, qdeg, start, startsols):
    r"""
    Generates a random Pieri problem of dimensions (*mdim, pdim, qdeg*)
    and solves it with a Cheater's homotopy, starting from
    the Pieri system in *start*, at the solutions in *startsols*.
    """
    dim = mdim*pdim + qdeg*(mdim+pdim)
    planes = [random_complex_matrix(mdim+pdim, mdim) for _ in range(0, dim)]
    pols = make_pieri_system(mdim, pdim, qdeg, planes)
    from phcpy.trackers import track
    print('cheater homotopy with %d paths' % len(startsols))
    sols = track(pols, start, startsols)
    for sol in sols:
        print(sol)
    verify(pols, sols)

def osculating_input(mdim, pdim, qdeg, start, startsols):
    r"""
    Generates real *mdim*-planes osculating a rational normal curve
    and solves this Pieri problem using the system in *start*,
    with corresponding solutions in *startsols*.
    """
    target_planes = real_osculating_planes(mdim, pdim, qdeg)
    # print 'real osculating planes :', target_planes
    target_system = make_pieri_system(mdim, pdim, qdeg, target_planes, False)
    print('the start system of length %d :' % len(start))
    for pol in start:
        print(pol)
    print('the target system of length %d :' % len(target_system))
    for pol in target_system:
        print(pol)
    from phcpy.trackers import track
    target_solutions = track(target_system, start, startsols)
    for sol in target_solutions:
        print(sol)
    verify(target_system, target_solutions)

def test_lrhom(prc='d'):
    """
    Performs a test on the Littlewood-Richardson homotopies.
    """
    brk = [[2, 4, 6], [2, 4, 6], [2, 4, 6]]
    (roco, flags, fsys, sols) \
        = littlewood_richardson_homotopies(6, 3, brk, precision=prc)
    print('the root count :', roco)
    print('the flags :', flags)
    print('the solutions :')
    for sol in sols:
        print(sol)
    verify(fsys, sols)

def test_pieri():
    """
    Does a test on the Pieri homotopies.
    """
    (mdim, pdim, qdeg) = prompt_for_dimensions()
    pieri_root_count(mdim, pdim, qdeg)
    dim = mdim*pdim + qdeg*(mdim+pdim)
    planes = [random_complex_matrix(mdim+pdim, mdim) for k in range(0, dim)]
    # print '%d random %d-planes :' % (dim, m)
    # for A in planes:
    #    for row in A: print row
    #    print ''
    if(qdeg > 0):
        points = random_complex_matrix(dim, 1)
        print('interpolation points :')
        for point in points:
            print(point)
        (system, sols) = run_pieri_homotopies(mdim, pdim, qdeg, planes, points)
    else:
        (system, sols) = run_pieri_homotopies(mdim, pdim, qdeg, planes)
    print('evaluation of the solutions :')
    verify(system, sols)
    from phcpy.solver import newton_step
    print('verification with one Newton step :')
    newton_step(system, sols)
    # cheater(m, p, qdeg, system, sols)
    if(qdeg == 0):
        osculating_input(mdim, pdim, qdeg, system, sols)

def main():
    """
    Tests the Pieri homotopies and the Littlewood-Richardson homotopies.
    """
    print("\nTesting the Pieri homotopies ...\n")
    test_pieri()
    print("\nTesting the Littlewood-Richardson homotopies ...")
    test_lrhom()

if __name__ == "__main__":
    main()
