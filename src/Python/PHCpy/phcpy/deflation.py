"""
Deflation restores the quadratic convergence of Newton's method
at a singular, isolated solution.
"""
from ctypes import c_int32, c_double, pointer
from phcpy.version import get_phcfun
from phcpy.polynomials import number_of_symbols
from phcpy.polynomials import set_double_system
from phcpy.polynomials import set_double_double_system
from phcpy.polynomials import set_quad_double_system
from phcpy.solutions import set_double_solutions
from phcpy.solutions import get_double_solutions
from phcpy.solutions import set_double_double_solutions
from phcpy.solutions import get_double_double_solutions
from phcpy.solutions import set_quad_double_solutions
from phcpy.solutions import get_quad_double_solutions
from phcpy.solutions import make_solution, diagnostics

def double_newton_step(pols, sols, vrblvl=0):
    r"""
    Applies one Newton step to the *sols* of the *pols*,
    in double precision.
    """
    if vrblvl > 0:
        print('in double_newton_step')
        print('the polynomial system :')
        for pol in pols:
            print(pol)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    nvr = number_of_symbols(pols, vrblvl-1)
    set_double_system(nvr, pols, vrblvl-1)
    set_double_solutions(nvr, sols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_newton_step calls phc', end='')
    retval = phc(199, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    result = get_double_solutions(vrblvl-1)
    return result

def double_double_newton_step(pols, sols, vrblvl=0):
    r"""
    Applies one Newton step to the *sols* of the *pols*,
    in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_newton_step')
        print('the polynomial system :')
        for pol in pols:
            print(pol)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    nvr = number_of_symbols(pols, vrblvl-1)
    set_double_double_system(nvr, pols, vrblvl-1)
    set_double_double_solutions(nvr, sols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_newton_step calls phc', end='')
    retval = phc(198, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    result = get_double_double_solutions(vrblvl-1)
    return result

def quad_double_newton_step(pols, sols, vrblvl=0):
    r"""
    Applies one Newton step to the *sols* of the *pols*,
    in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_newton_step')
        print('the polynomial system :')
        for pol in pols:
            print(pol)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    nvr = number_of_symbols(pols, vrblvl-1)
    set_quad_double_system(nvr, pols, vrblvl-1)
    set_quad_double_solutions(nvr, sols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_newton_step calls phc', end='')
    retval = phc(197, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    result = get_quad_double_solutions(vrblvl-1)
    return result

def double_deflate(pols, sols, maxitr=3, maxdef=3, \
    tolerr=1.0e-8, tolres=1.0e-8, tolrnk=1.0e-6, vrblvl=0):
    r"""
    Applies deflation on the solutions sols of the system in pols,
    in double precision.  The parameters are as follows:

    *maxitr*: the maximum number of iterations per root,

    *maxdef*: the maximum number of deflations per root,

    *tolerr*: tolerance on the forward error on each root,

    *tolres*: tolerance on the backward error on each root,

    *tolrnk*: tolerance on the numerical rank of the Jacobian matrices.
    """
    if vrblvl > 0:
        print('in double_deflate, maxitr :', maxitr, end='')
        print(', maxdef :', maxdef)
        print(f'tolerr : {tolerr:.3e}', end='')
        print(f', tolres : {tolres:.3e}, tolrnk : {tolrnk:.3e}')
        print('the polynomial system :')
        for pol in pols:
            print(pol)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    nvr = number_of_symbols(pols, vrblvl-1)
    set_double_system(nvr, pols, vrblvl-1)
    set_double_solutions(nvr, sols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    amaxitr = pointer(c_int32(maxitr))
    bmaxdef = pointer(c_int32(maxdef))
    tols = (c_double * 3)()
    tols[0] = tolerr
    tols[1] = tolres
    tols[2] = tolrnk
    ctols = pointer(tols)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_deflate calls phc', end='')
    retval = phc(196, amaxitr, bmaxdef, ctols, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    result = get_double_solutions(vrblvl-1)
    return result

def double_double_deflate(pols, sols, maxitr=3, maxdef=3, \
    tolerr=1.0e-8, tolres=1.0e-8, tolrnk=1.0e-6, vrblvl=0):
    r"""
    Applies deflation on the solutions sols of the system in pols,
    in double double precision.  The parameters are as follows:

    *maxitr*: the maximum number of iterations per root,

    *maxdef*: the maximum number of deflations per root,

    *tolerr*: tolerance on the forward error on each root,

    *tolres*: tolerance on the backward error on each root,

    *tolrnk*: tolerance on the numerical rank of the Jacobian matrices.
    """
    if vrblvl > 0:
        print('in double_double_deflate, maxitr :', maxitr, end='')
        print(', maxdef :', maxdef)
        print(f'tolerr : {tolerr:.3e}', end='')
        print(f', tolres : {tolres:.3e}, tolrnk : {tolrnk:.3e}')
        print('the polynomial system :')
        for pol in pols:
            print(pol)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    nvr = number_of_symbols(pols, vrblvl-1)
    set_double_double_system(nvr, pols, vrblvl-1)
    set_double_double_solutions(nvr, sols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    amaxitr = pointer(c_int32(maxitr))
    bmaxdef = pointer(c_int32(maxdef))
    tols = (c_double * 3)()
    tols[0] = tolerr
    tols[1] = tolres
    tols[2] = tolrnk
    ctols = pointer(tols)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_deflate calls phc', end='')
    retval = phc(249, amaxitr, bmaxdef, ctols, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    result = get_double_double_solutions(vrblvl-1)
    return result

def quad_double_deflate(pols, sols, maxitr=3, maxdef=3, \
    tolerr=1.0e-8, tolres=1.0e-8, tolrnk=1.0e-6, vrblvl=0):
    r"""
    Applies deflation on the solutions sols of the system in pols,
    in quad double precision.  The parameters are as follows:

    *maxitr*: the maximum number of iterations per root,

    *maxdef*: the maximum number of deflations per root,

    *tolerr*: tolerance on the forward error on each root,

    *tolres*: tolerance on the backward error on each root,

    *tolrnk*: tolerance on the numerical rank of the Jacobian matrices.
    """
    if vrblvl > 0:
        print('in quad_double_deflate, maxitr :', maxitr, end='')
        print(', maxdef :', maxdef)
        print(f'tolerr : {tolerr:.3e}', end='')
        print(f', tolres : {tolres:.3e}, tolrnk : {tolrnk:.3e}')
        print('the polynomial system :')
        for pol in pols:
            print(pol)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    nvr = number_of_symbols(pols, vrblvl-1)
    set_quad_double_system(nvr, pols, vrblvl-1)
    set_quad_double_solutions(nvr, sols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    amaxitr = pointer(c_int32(maxitr))
    bmaxdef = pointer(c_int32(maxdef))
    tols = (c_double * 3)()
    tols[0] = tolerr
    tols[1] = tolres
    tols[2] = tolrnk
    ctols = pointer(tols)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_deflate calls phc', end='')
    retval = phc(250, amaxitr, bmaxdef, ctols, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    result = get_quad_double_solutions(vrblvl-1)
    return result

def double_multiplicity(system, solution, \
    order=5, tol=1.0e-8, vrblvl=0):
    r"""
    Computes the multiplicity structure in double precision
    of an isolated solution (in the string *solution*)
    of a polynomial system (in the list *system*).
    The other parameters are

    *order*: the maximum order of differentiation,

    *tol*: tolerance on the numerical rank,

    *vrblvl*: is the verbose level.

    On return is the computed multiplicity.
    """
    if vrblvl > 0:
        print('in double_multiplicity ...')
        print('the polynomial system :')
        for pol in system:
            print(pol)
        print('the solution :')
        print(solution)
    dim = number_of_symbols(system, vrblvl-1)
    set_double_system(dim, system, vrblvl-1)
    set_double_solutions(dim, [solution], vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    amlt = pointer(c_int32(order))
    size = order + 1
    bhlb = (c_int32 * size)()
    bhlb[0] = vrblvl
    phlb = pointer(bhlb)
    ctol = pointer(c_double(tol))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_multiplicity calls phc', end='')
    retval = phc(732, amlt, phlb, ctol, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    vals = phlb[0:size]
    hlbfun = []
    for idx in range(size):
        hlbfun.append(int(vals[0][idx]))
    if vrblvl > 0:
        print('The values of the Hilbert function :', hlbfun)
        print('The multiplicity :', amlt[0])
    return (amlt[0], hlbfun)

def double_double_multiplicity(system, solution, \
    order=5, tol=1.0e-8, vrblvl=0):
    r"""
    Computes the multiplicity structure in double double precision
    of an isolated solution (in the string *solution*)
    of a polynomial system (in the list *system*).
    The other parameters are

    *order*: the maximum order of differentiation,

    *tol*: tolerance on the numerical rank,

    *vrblvl*: is the verbose level.

    On return is the computed multiplicity.
    """
    if vrblvl > 0:
        print('in double_double_multiplicity ...')
        print('the polynomial system :')
        for pol in system:
            print(pol)
        print('the solution :')
        print(solution)
    dim = number_of_symbols(system, vrblvl-1)
    set_double_double_system(dim, system, vrblvl-1)
    set_double_double_solutions(dim, [solution], vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    amlt = pointer(c_int32(order))
    size = order + 1
    bhlb = (c_int32 * size)()
    bhlb[0] = vrblvl
    phlb = pointer(bhlb)
    ctol = pointer(c_double(tol))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_multiplicity calls phc', end='')
    retval = phc(733, amlt, phlb, ctol, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    vals = phlb[0:size]
    hlbfun = []
    for idx in range(size):
        hlbfun.append(int(vals[0][idx]))
    if vrblvl > 0:
        print('The values of the Hilbert function :', hlbfun)
        print('The multiplicity :', amlt[0])
    return (amlt[0], hlbfun)

def quad_double_multiplicity(system, solution, \
    order=5, tol=1.0e-8, vrblvl=0):
    r"""
    Computes the multiplicity structure in quad double precision
    of an isolated solution (in the string *solution*)
    of a polynomial system (in the list *system*).
    The other parameters are

    *order*: the maximum order of differentiation,

    *tol*: tolerance on the numerical rank,

    *vrblvl*: is the verbose level.

    On return is the computed multiplicity.
    """
    if vrblvl > 0:
        print('in quad_double_multiplicity ...')
        print('the polynomial system :')
        for pol in system:
            print(pol)
        print('the solution :')
        print(solution)
    dim = number_of_symbols(system, vrblvl-1)
    set_quad_double_system(dim, system, vrblvl-1)
    set_quad_double_solutions(dim, [solution], vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    amlt = pointer(c_int32(order))
    size = order + 1
    bhlb = (c_int32 * size)()
    bhlb[0] = vrblvl
    phlb = pointer(bhlb)
    ctol = pointer(c_double(tol))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_multiplicity calls phc', end='')
    retval = phc(734, amlt, phlb, ctol, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    vals = phlb[0:size]
    hlbfun = []
    for idx in range(size):
        hlbfun.append(int(vals[0][idx]))
    if vrblvl > 0:
        print('The values of the Hilbert function :', hlbfun)
        print('The multiplicity :', amlt[0])
    return (amlt[0], hlbfun)

def test_double_newton_step(vrblvl=0):
    """
    Tests Newton's method in double precision.
    """
    if vrblvl > 0:
        print('in test_double_newton_step ...')
    pols = ['x^2 - 1;', 'x*y - 1;']
    sol = make_solution(['x', 'y'], [0.9999, 0.9999])
    sols = [sol]
    newsols = double_newton_step(pols, sols, vrblvl)
    if vrblvl > 0:
        print('solution after a first Newton step :')
        print(newsols[0])
    newsols = double_newton_step(pols, newsols, vrblvl)
    if vrblvl > 0:
        print('solution after a second Newton step :')
        print(newsols[0])
    err, rco, res = diagnostics(newsols[0], vrblvl-1)
    if vrblvl > 0:
        print(f'forward error : {err:.3e}')
        print(f'estimate for inverse condition number : {rco:.3e}')
        print(f'backward error : {res:.3e}')
    return int(res > 1.0e-12)

def test_double_double_newton_step(vrblvl=0):
    """
    Tests Newton's method in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_newton_step ...')
    pols = ['x^2 - 1;', 'x*y - 1;']
    sol = make_solution(['x', 'y'], [0.9999, 0.9999])
    sols = [sol]
    newsols = double_double_newton_step(pols, sols, vrblvl)
    if vrblvl > 0:
        print('solution after a first Newton step :')
        print(newsols[0])
    newsols = double_double_newton_step(pols, newsols, vrblvl)
    if vrblvl > 0:
        print('solution after a second Newton step :')
        print(newsols[0])
    newsols = double_double_newton_step(pols, newsols, vrblvl)
    if vrblvl > 0:
        print('solution after a third Newton step :')
        print(newsols[0])
    err, rco, res = diagnostics(newsols[0], vrblvl-1)
    if vrblvl > 0:
        print(f'forward error : {err:.3e}')
        print(f'estimate for inverse condition number : {rco:.3e}')
        print(f'backward error : {res:.3e}')
    return int(res > 1.0e-24)

def test_quad_double_newton_step(vrblvl=0):
    """
    Tests Newton's method in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_newton_step ...')
    pols = ['x^2 - 1;', 'x*y - 1;']
    sol = make_solution(['x', 'y'], [0.9999, 0.9999])
    sols = [sol]
    newsols = quad_double_newton_step(pols, sols, vrblvl)
    if vrblvl > 0:
        print('solution after a first Newton step :')
        print(newsols[0])
    newsols = quad_double_newton_step(pols, newsols, vrblvl)
    if vrblvl > 0:
        print('solution after a second Newton step :')
        print(newsols[0])
    newsols = quad_double_newton_step(pols, newsols, vrblvl)
    if vrblvl > 0:
        print('solution after a third Newton step :')
        print(newsols[0])
    newsols = quad_double_newton_step(pols, newsols, vrblvl)
    if vrblvl > 0:
        print('solution after a fourth Newton step :')
        print(newsols[0])
    err, rco, res = diagnostics(newsols[0], vrblvl-1)
    if vrblvl > 0:
        print(f'forward error : {err:.3e}')
        print(f'estimate for inverse condition number : {rco:.3e}')
        print(f'backward error : {res:.3e}')
    return int(res > 1.0e-48)

def test_double_deflate(vrblvl=0):
    """
    Tests deflation in double precision,
    on the 'ojika1' demonstration test system, in the paper by
    T. Ojika: "Modified deflation algorithm for the solution of 
    singular problems. I. A system of nonlinear algebraic equations."
    J. Math. Anal. Appl. 123, 199-221, 1987.
    """
    pols = ['x^2 + y - 3;', 'x + 0.125*y^2 - 1.5;']
    regsol = make_solution(['x', 'y'], [-2.999, -5.999]) # regular
    sinsol = make_solution(['x', 'y'], [0.9999, 1.999])  # singular
    sols = [regsol, sinsol]
    defsols = double_deflate(pols, sols, vrblvl=vrblvl)
    print('the first solution after deflation :')
    print(defsols[0])
    print('the second solution after deflation :')
    print(defsols[1])
    err, rco, res = diagnostics(defsols[1], vrblvl-1)
    if vrblvl > 0:
        print(f'forward error : {err:.3e}')
        print(f'estimate for inverse condition number : {rco:.3e}')
        print(f'backward error : {res:.3e}')
    return int(res > 1.0e-12)

def test_double_double_deflate(vrblvl=0):
    """
    Tests deflation in double double precision,
    on the 'ojika1' demonstration test system, in the paper by
    T. Ojika: "Modified deflation algorithm for the solution of 
    singular problems. I. A system of nonlinear algebraic equations."
    J. Math. Anal. Appl. 123, 199-221, 1987.
    """
    pols = ['x^2 + y - 3;', 'x + 0.125*y^2 - 1.5;']
    regsol = make_solution(['x', 'y'], [-2.999, -5.999]) # regular
    sinsol = make_solution(['x', 'y'], [0.9999, 1.999])  # singular
    sols = [regsol, sinsol]
    defsols = double_double_deflate(pols, sols, vrblvl=vrblvl)
    print('the first solution after deflation :')
    print(defsols[0])
    print('the second solution after deflation :')
    print(defsols[1])
    err, rco, res = diagnostics(defsols[1], vrblvl-1)
    if vrblvl > 0:
        print(f'forward error : {err:.3e}')
        print(f'estimate for inverse condition number : {rco:.3e}')
        print(f'backward error : {res:.3e}')
    return int(res > 1.0e-12)

def test_quad_double_deflate(vrblvl=0):
    """
    Tests deflation in quad double precision,
    on the 'ojika1' demonstration test system, in the paper by
    T. Ojika: "Modified deflation algorithm for the solution of 
    singular problems. I. A system of nonlinear algebraic equations."
    J. Math. Anal. Appl. 123, 199-221, 1987.
    """
    pols = ['x^2 + y - 3;', 'x + 0.125*y^2 - 1.5;']
    regsol = make_solution(['x', 'y'], [-2.999, -5.999]) # regular
    sinsol = make_solution(['x', 'y'], [0.9999, 1.999])  # singular
    sols = [regsol, sinsol]
    defsols = quad_double_deflate(pols, sols, vrblvl=vrblvl)
    print('the first solution after deflation :')
    print(defsols[0])
    print('the second solution after deflation :')
    print(defsols[1])
    err, rco, res = diagnostics(defsols[1], vrblvl-1)
    if vrblvl > 0:
        print(f'forward error : {err:.3e}')
        print(f'estimate for inverse condition number : {rco:.3e}')
        print(f'backward error : {res:.3e}')
    return int(res > 1.0e-12)

def test_double_multiplicity(vrblvl=0):
    """
    Tests the multiplicity structure in double precision.
    """
    if vrblvl > 0:
        print('in test_double_multiplicity ...')
    pols = [ 'x**2+y-3;', 'x+0.125*y**2-1.5;']
    sol = make_solution(['x', 'y'], [1, 2])
    multiplicity, hilfun = double_multiplicity(pols, sol, vrblvl=vrblvl)
    if vrblvl > 0:
        print('the multiplicity :', multiplicity)
        print('the Hibert function :', hilfun)
    return int(multiplicity != 3)

def test_double_double_multiplicity(vrblvl=0):
    """
    Tests the multiplicity structure in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_multiplicity ...')
    pols = [ 'x**2+y-3;', 'x+0.125*y**2-1.5;']
    sol = make_solution(['x', 'y'], [1, 2])
    multiplicity, hilfun = double_double_multiplicity(pols, sol, \
        vrblvl=vrblvl)
    if vrblvl > 0:
        print('the multiplicity :', multiplicity)
        print('the Hibert function :', hilfun)
    return int(multiplicity != 3)

def test_quad_double_multiplicity(vrblvl=0):
    """
    Tests the multiplicity structure in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_multiplicity ...')
    pols = [ 'x**2+y-3;', 'x+0.125*y**2-1.5;']
    sol = make_solution(['x', 'y'], [1, 2])
    multiplicity, hilfun = quad_double_multiplicity(pols, sol, \
        vrblvl=vrblvl)
    if vrblvl > 0:
        print('the multiplicity :', multiplicity)
        print('the Hibert function :', hilfun)
    return int(multiplicity != 3)

def main():
    """
    Runs some tests.
    """
    lvl = 1
    fail = test_double_newton_step(lvl)
    fail = fail + test_double_double_newton_step(lvl)
    fail = fail + test_quad_double_newton_step(lvl)
    fail = fail + test_double_deflate(lvl)
    fail = fail + test_double_double_deflate(lvl)
    fail = fail + test_quad_double_deflate(lvl)
    fail = fail + test_double_multiplicity(lvl)
    fail = fail + test_double_double_multiplicity(lvl)
    fail = fail + test_quad_double_multiplicity(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()
