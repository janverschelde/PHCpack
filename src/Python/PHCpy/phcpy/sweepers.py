"""
The module sweepers exports the definition of sweep homotopies and
the tracking of solution paths defined by sweep homotopies.
A sweep homotopy is a polynomial system where some of the variables
are considered as parameters.  Given solutions for some parameters
and new values for the parameters, we can track the solution paths
starting at the given solutions and ending at the new solutions for
the new values of the parameters.
The sweep is controlled by a convex linear combination between the
list of start and target values for the parameters.
We distinguish between a complex and a real sweep.
In a complex sweep, with a randomly generated gamma we avoid singularities
along the solution paths, in a complex convex combination between the
start and target values for the parameters.  This complex sweep is
applicable only when the parameter space is convex.
In a real sweep, arc-length parameter continuation is applied.
The algorithms applied in this module are described in the paper by
Kathy Piret and Jan Verschelde: Sweeping Algebraic Curves for Singular 
Solutions.  Journal of Computational and Applied Mathematics,
volume 234, number 4, pages 1228-1237, 2010. 
"""
from ctypes import c_int32, c_double, pointer
from math import sqrt
from phcpy.version import get_phcfun
from phcpy.polynomials import set_double_system
from phcpy.polynomials import set_double_double_system
from phcpy.polynomials import set_quad_double_system
from phcpy.solutions import set_double_solutions, get_double_solutions
from phcpy.solutions import set_double_double_solutions
from phcpy.solutions import get_double_double_solutions
from phcpy.solutions import set_quad_double_solutions
from phcpy.solutions import get_quad_double_solutions
from phcpy.solutions import make_solution, verify

def set_double_start(start, vrblvl=0):
    """
    Sets the values of all parameters to start,
    which contains the consecutive values of all real and imaginary parts
    of the start values of all parameters, in double precision.
    """
    if vrblvl > 0:
        print('in set_double_start, start :')
        for cff in start:
            print(cff)
    size = len(start)
    phc = get_phcfun(vrblvl-1)
    pars = (c_int32 * 2)()
    pars[0] = 0 # double precision
    pars[1] = 0 # start values
    apar = pointer(pars)
    bnum = pointer(c_int32(size))
    vals = (c_double * size)()
    for (idx, coefficient) in enumerate(start):
        vals[idx] = c_double(coefficient)
    pval = pointer(vals)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_start calls phc', end='')
    retval = phc(618, apar, bnum, pval, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_double_start(start, vrblvl=0):
    """
    Sets the values of all parameters to start,
    which contains the consecutive values of all real and imaginary parts
    of the start values of all parameters, in double double precision.
    The length of start must be twice the size of the start
    in double precision.
    """
    if vrblvl > 0:
        print('in set_double_double_start, start :')
        for cff in start:
            print(cff)
    size = len(start)
    phc = get_phcfun(vrblvl-1)
    pars = (c_int32 * 2)()
    pars[0] = 1 # double double precision
    pars[1] = 0 # start values
    apar = pointer(pars)
    bnum = pointer(c_int32(size))
    vals = (c_double * size)()
    for (idx, coefficient) in enumerate(start):
        vals[idx] = c_double(coefficient)
    pval = pointer(vals)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_double_start calls phc', end='')
    retval = phc(618, apar, bnum, pval, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_quad_double_start(start, vrblvl=0):
    """
    Sets the values of all parameters to start,
    which contains the consecutive values of all real and imaginary parts
    of the start values of all parameters, in quad double precision.
    The length of start must be four times the size of the start
    in double precision.
    """
    if vrblvl > 0:
        print('in set_quad_double_start, start :')
        for cff in start:
            print(cff)
    size = len(start)
    phc = get_phcfun(vrblvl-1)
    pars = (c_int32 * 2)()
    pars[0] = 2 # quad double precision
    pars[1] = 0 # start values
    apar = pointer(pars)
    bnum = pointer(c_int32(size))
    vals = (c_double * size)()
    for (idx, coefficient) in enumerate(start):
        vals[idx] = c_double(coefficient)
    pval = pointer(vals)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_quad_double_start calls phc', end='')
    retval = phc(618, apar, bnum, pval, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_target(target, vrblvl=0):
    """
    Sets the values of all parameters to target,
    which contains the consecutive values of all real and imaginary parts
    of the target values of all parameters, in double precision.
    """
    if vrblvl > 0:
        print('in set_double_target, target :')
        for cff in target:
            print(cff)
    size = len(target)
    phc = get_phcfun(vrblvl-1)
    pars = (c_int32 * 2)()
    pars[0] = 0 # double precision
    pars[1] = 1 # target values
    apar = pointer(pars)
    bnum = pointer(c_int32(size))
    vals = (c_double * size)()
    for (idx, coefficient) in enumerate(target):
        vals[idx] = c_double(coefficient)
    pval = pointer(vals)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_target calls phc', end='')
    retval = phc(618, apar, bnum, pval, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_double_target(target, vrblvl=0):
    """
    Sets the values of all parameters to target,
    which contains the consecutive values of all real and imaginary parts
    of the target values of all parameters, in double double precision.
    The length of start must be twice the size of the start
    in double precision.
    """
    if vrblvl > 0:
        print('in set_double_double_target, target :')
        for cff in target:
            print(cff)
    size = len(target)
    phc = get_phcfun(vrblvl-1)
    pars = (c_int32 * 2)()
    pars[0] = 1 # double double precision
    pars[1] = 1 # target values
    apar = pointer(pars)
    bnum = pointer(c_int32(size))
    vals = (c_double * size)()
    for (idx, coefficient) in enumerate(target):
        vals[idx] = c_double(coefficient)
    pval = pointer(vals)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_double_target calls phc', end='')
    retval = phc(618, apar, bnum, pval, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_quad_double_target(target, vrblvl=0):
    """
    Sets the values of all parameters to target,
    which contains the consecutive values of all real and imaginary parts
    of the target values of all parameters, in double double precision.
    The length of start must be four times the size of the start
    in double precision.
    """
    if vrblvl > 0:
        print('in set_quad_double_target, target :')
        for cff in target:
            print(cff)
    size = len(target)
    phc = get_phcfun(vrblvl-1)
    pars = (c_int32 * 2)()
    pars[0] = 2 # quad double precision
    pars[1] = 1 # target values
    apar = pointer(pars)
    bnum = pointer(c_int32(size))
    vals = (c_double * size)()
    for (idx, coefficient) in enumerate(target):
        vals[idx] = c_double(coefficient)
    pval = pointer(vals)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_quad_double_target calls phc', end='')
    retval = phc(618, apar, bnum, pval, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_parameter_names(neq, nvr, pars, vrblvl=0):
    """
    Defines which variables serve as parameters,
    by providing the names of the parameters in the list pars.
    The number of equations is given in neq and
    the number of variables is in nvr.
    """
    if vrblvl > 0:
        print('in set_parameter_names, neq :', neq, end='')
        print(', nvr :', nvr)
        print('parameters :', pars)
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 4)()
    apars[0] = neq
    apars[1] = nvr
    apars[2] = len(pars)
    parnames = ' '.join(pars)
    nbc = len(parnames)
    apars[3] = nbc
    apar = pointer(apars)
    size = nbc
    names = (c_int32 *size)()
    for (idx, letter) in enumerate(parnames):
        names[idx] = c_int32(ord(letter))
    bnames = pointer(names)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_parameter_names calls phc', end='')
    retval = phc(611, apar, bnames, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_complex_sweep_run(gchoice, regamma, imgamma, vrblvl=0):
    """
    Starts the trackers in a complex convex parameter homotopy,
    in double precision, where the indices to the parameters,
    start and target values are already defined,
    and the systems and solution are set in double precision.
    The input parameter gchoice is 0, 1, or 2, for respectively
    a randomly generated gamma (0), or no gamma (1), or a user given
    gamma with real and imaginary parts in regamma and imgamma.
    With a random gamma, this is known as cheater's homotopy.
    """
    if vrblvl > 0:
        print('in double_complex_sweep_run, gchoice :', gchoice)
        print('regamma :', regamma, '  imgamma :', imgamma)
    phc = get_phcfun(vrblvl-1)
    pars = (c_int32 * 2)()
    pars[0] = 0  # double precision
    pars[1] = gchoice
    apar = pointer(pars)
    bbb = pointer(c_int32(0))
    gamma = (c_double * 2)()
    gamma[0] = regamma
    gamma[1] = imgamma
    cgamma = pointer(gamma)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_complex_sweep_run calls phc', end='')
    retval = phc(620, apar, bbb, cgamma, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_double_complex_sweep_run(gchoice, regamma, imgamma, vrblvl=0):
    """
    Starts the trackers in a complex convex parameter homotopy,
    in double double precision, where the indices to the parameters,
    start and target values are already defined,
    and the systems and solution are set in double double precision.
    The input parameter gchoice is 0, 1, or 2, for respectively
    a randomly generated gamma (0), or no gamma (1), or a user given
    gamma with real and imaginary parts in regamma and imgamma.
    With a random gamma, this is known as cheater's homotopy.
    """
    if vrblvl > 0:
        print('in double_double_complex_sweep_run, gchoice :', gchoice)
        print('regamma :', regamma, '  imgamma :', imgamma)
    phc = get_phcfun(vrblvl-1)
    pars = (c_int32 * 2)()
    pars[0] = 1  # double double precision
    pars[1] = gchoice
    apar = pointer(pars)
    bbb = pointer(c_int32(0))
    gamma = (c_double * 4)()
    gamma[0] = regamma
    gamma[1] = 0.0
    gamma[2] = imgamma
    gamma[3] = 0.0
    cgamma = pointer(gamma)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_complex_sweep_run calls phc', end='')
    retval = phc(620, apar, bbb, cgamma, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def quad_double_complex_sweep_run(gchoice, regamma, imgamma, vrblvl=0):
    """
    Starts the trackers in a complex convex parameter homotopy,
    in quad double precision, where the indices to the parameters,
    start and target values are already defined,
    and the systems and solution are set in double double precision.
    The input parameter gchoice is 0, 1, or 2, for respectively
    a randomly generated gamma (0), or no gamma (1), or a user given
    gamma with real and imaginary parts in regamma and imgamma.
    With a random gamma, this is known as cheater's homotopy.
    """
    if vrblvl > 0:
        print('in quad_double_complex_sweep_run, gchoice :', gchoice)
        print('regamma :', regamma, '  imgamma :', imgamma)
    phc = get_phcfun(vrblvl-1)
    pars = (c_int32 * 2)()
    pars[0] = 2  # quad double precision
    pars[1] = gchoice
    apar = pointer(pars)
    bbb = pointer(c_int32(0))
    gamma = (c_double * 8)()
    gamma[0] = regamma
    gamma[1] = 0.0
    gamma[2] = 0.0
    gamma[3] = 0.0
    gamma[4] = imgamma
    gamma[5] = 0.0
    gamma[6] = 0.0
    gamma[7] = 0.0
    cgamma = pointer(gamma)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_complex_sweep_run calls phc', end='')
    retval = phc(620, apar, bbb, cgamma, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_real_sweep_run(vrblvl=0):
    """
    Starts a sweep with a natural parameter in a family of n equations
    in n+1 variables, where the last variable is the artificial parameter s
    that moves the one natural parameter from a start to target value.
    The last equation is of the form (1-s)*(A - v[0]) + s*(A - v[1]),
    where A is the natural parameter, going from the start value v[0]
    to the target value v[1], already set in double precision.
    Also set before calling the function are the homotopy and the start
    solutions, where every solution has the value v[0] for the A variable.
    The sweep stops when s reaches the value v[1],
    or when a singularity is encountered on the path.
    """
    if vrblvl > 0:
        print('in double_real_sweep ...')
    phc = get_phcfun(vrblvl-1)
    aprc = pointer(c_int32(0)) # double precision
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_real_sweep_run calls phc', end='')
    retval = phc(621, aprc, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_double_real_sweep_run(vrblvl=0):
    """
    Starts a sweep with a natural parameter in a family of n equations
    in n+1 variables, where the last variable is the artificial parameter s
    that moves the one natural parameter from a start to target value.
    The last equation is of the form (1-s)*(A - v[0]) + s*(A - v[1]),
    where A is the natural parameter, going from the start value v[0]
    to the target value v[1], already set in double double precision.
    Also set before calling the function are the homotopy and the start
    solutions, where every solution has the value v[0] for the A variable.
    The sweep stops when s reaches the value v[1],
    or when a singularity is encountered on the path.
    """
    if vrblvl > 0:
        print('in double_double_real_sweep ...')
    phc = get_phcfun(vrblvl-1)
    aprc = pointer(c_int32(1)) # double double precision
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_real_sweep_run calls phc', end='')
    retval = phc(621, aprc, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def quad_double_real_sweep_run(vrblvl=0):
    """
    Starts a sweep with a natural parameter in a family of n equations
    in n+1 variables, where the last variable is the artificial parameter s
    that moves the one natural parameter from a start to target value.
    The last equation is of the form (1-s)*(A - v[0]) + s*(A - v[1]),
    where A is the natural parameter, going from the start value v[0]
    to the target value v[1], already set in quad double precision.
    Also set before calling the function are the homotopy and the start
    solutions, where every solution has the value v[0] for the A variable.
    The sweep stops when s reaches the value v[1],
    or when a singularity is encountered on the path.
    """
    if vrblvl > 0:
        print('in quad_double_real_sweep ...')
    phc = get_phcfun(vrblvl-1)
    aprc = pointer(c_int32(2)) # quad double precision
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_real_sweep_run calls phc', end='')
    retval = phc(621, aprc, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_complex_sweep(pols, sols, nvar, pars, start, target, vrblvl=0):
    r"""
    For the polynomials in the list of strings *pols*
    and the solutions in *sols* for the values in the list *start*,
    a sweep through the parameter space will be performed
    in double precision to the target values of
    the parameters in the list *target*.
    The number of variables in the polynomials and the solutions
    must be the same and be equal to the value of *nvar*.
    The list of symbols in *pars* contains the names of the variables
    in the polynomials *pols* that serve as parameters.
    The size of the lists *pars*, *start*, and *target* must be same.
    """
    if vrblvl > 0:
        print('in double_complex_sweep, nvar :', nvar)
        print('the polynomials :')
        for pol in pols:
            print(pol)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
        print('parameters :', pars)
        print('start :', start, '  target :', target)
    set_double_system(nvar, pols, vrblvl-1)
    set_double_solutions(nvar, sols, vrblvl-1)
    set_parameter_names(len(pols), nvar, pars, vrblvl)
    set_double_start(start, vrblvl)
    set_double_target(target, vrblvl)
    double_complex_sweep_run(0, 0.0, 0.0)
    result = get_double_solutions(vrblvl-1)
    return result

def double_double_complex_sweep(pols, sols, nvar, pars, start, target, \
    vrblvl=0):
    r"""
    For the polynomials in the list of strings *pols*
    and the solutions in *sols* for the values in the list *start*,
    a sweep through the parameter space will be performed
    in double double precision to the target values of
    the parameters in the list *target*.
    The number of variables in the polynomials and the solutions
    must be the same and be equal to the value of *nvar*.
    The list of symbols in *pars* contains the names of the variables
    in the polynomials *pols* that serve as parameters.
    The size of the lists *pars*, *start*, and *target* must be same.
    """
    if vrblvl > 0:
        print('in double_double_complex_sweep, nvar :', nvar)
        print('the polynomials :')
        for pol in pols:
            print(pol)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
        print('parameters :', pars)
        print('start :', start, '  target :', target)
    set_double_double_system(nvar, pols, vrblvl-1)
    set_double_double_solutions(nvar, sols, vrblvl-1)
    set_parameter_names(len(pols), nvar, pars, vrblvl)
    set_double_double_start(start, vrblvl)
    set_double_double_target(target, vrblvl)
    double_double_complex_sweep_run(0, 0.0, 0.0)
    result = get_double_double_solutions(vrblvl-1)
    return result

def quad_double_complex_sweep(pols, sols, nvar, pars, start, target, \
    vrblvl=0):
    r"""
    For the polynomials in the list of strings *pols*
    and the solutions in *sols* for the values in the list *start*,
    a sweep through the parameter space will be performed
    in double double precision to the target values of
    the parameters in the list *target*.
    The number of variables in the polynomials and the solutions
    must be the same and be equal to the value of *nvar*.
    The list of symbols in *pars* contains the names of the variables
    in the polynomials *pols* that serve as parameters.
    The size of the lists *pars*, *start*, and *target* must be same.
    """
    if vrblvl > 0:
        print('in quad_double_complex_sweep, nvar :', nvar)
        print('the polynomials :')
        for pol in pols:
            print(pol)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
        print('parameters :', pars)
        print('start :', start, '  target :', target)
    set_quad_double_system(nvar, pols, vrblvl-1)
    set_quad_double_solutions(nvar, sols, vrblvl-1)
    set_parameter_names(len(pols), nvar, pars, vrblvl)
    set_quad_double_start(start, vrblvl)
    set_quad_double_target(target, vrblvl)
    quad_double_complex_sweep_run(0, 0.0, 0.0)
    result = get_quad_double_solutions(vrblvl-1)
    return result

def double_real_sweep(pols, sols, par='s', start=0.0, target=1.0, vrblvl=0):
    r"""
    A real sweep homotopy is a family of n equations in n+1 variables,
    where one of the variables is the artificial parameter s which moves
    from 0.0 to 1.0.  The last equation can then be of the form

    (1 - s)*(lambda - L[0]) + s*(lambda - L[1]) = 0 so that,

    at s = 0, the natural parameter lambda has the value L[0], and

    at s = 1, the natural parameter lambda has the value L[1].

    Thus: as s moves from 0 to 1, lambda goes from L[0] to L[1].

    All solutions in the list *sols* must have then the value L[0]
    for the variable lambda.
    The sweep stops when the target value for s is reached
    or when a singular solution is encountered.
    Computations happen in double precision.
    """
    if vrblvl > 0:
        print('in double_real_sweep ...')
        print('the polynomials :')
        for pol in pols:
            print(pol)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
        print('parameters :', par)
        print('start :', start, '  target :', target)
    nvar = len(pols) + 1
    set_double_system(nvar, pols, vrblvl-1)
    set_double_solutions(nvar, sols, vrblvl-1)
    set_parameter_names(len(pols), nvar, [par], vrblvl)
    set_double_start([start], vrblvl)
    set_double_target([target], vrblvl)
    double_real_sweep_run(vrblvl)
    result = get_double_solutions(vrblvl-1)
    return result

def double_double_real_sweep(pols, sols, par='s', start=0.0, target=1.0, \
    vrblvl=0):
    r"""
    A real sweep homotopy is a family of n equations in n+1 variables,
    where one of the variables is the artificial parameter s which moves
    from 0.0 to 1.0.  The last equation can then be of the form

    (1 - s)*(lambda - L[0]) + s*(lambda - L[1]) = 0 so that,

    at s = 0, the natural parameter lambda has the value L[0], and

    at s = 1, the natural parameter lambda has the value L[1].

    Thus: as s moves from 0 to 1, lambda goes from L[0] to L[1].

    All solutions in the list *sols* must have then the value L[0]
    for the variable lambda.
    The sweep stops when the target value for s is reached
    or when a singular solution is encountered.
    Computations happen in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_real_sweep ...')
        print('the polynomials :')
        for pol in pols:
            print(pol)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
        print('parameters :', par)
        print('start :', start, '  target :', target)
    nvar = len(pols) + 1
    set_double_double_system(nvar, pols, vrblvl-1)
    set_double_double_solutions(nvar, sols, vrblvl-1)
    set_parameter_names(len(pols), nvar, [par], vrblvl)
    set_double_double_start([start], vrblvl)
    set_double_double_target([target], vrblvl)
    double_double_real_sweep_run(vrblvl)
    result = get_double_double_solutions(vrblvl-1)
    return result

def quad_double_real_sweep(pols, sols, par='s', start=0.0, target=1.0, \
    vrblvl=0):
    r"""
    A real sweep homotopy is a family of n equations in n+1 variables,
    where one of the variables is the artificial parameter s which moves
    from 0.0 to 1.0.  The last equation can then be of the form

    (1 - s)*(lambda - L[0]) + s*(lambda - L[1]) = 0 so that,

    at s = 0, the natural parameter lambda has the value L[0], and

    at s = 1, the natural parameter lambda has the value L[1].

    Thus: as s moves from 0 to 1, lambda goes from L[0] to L[1].

    All solutions in the list *sols* must have then the value L[0]
    for the variable lambda.
    The sweep stops when the target value for s is reached
    or when a singular solution is encountered.
    Computations happen in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_real_sweep ...')
        print('the polynomials :')
        for pol in pols:
            print(pol)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
        print('parameters :', par)
        print('start :', start, '  target :', target)
    nvar = len(pols) + 1
    set_quad_double_system(nvar, pols, vrblvl-1)
    set_quad_double_solutions(nvar, sols, vrblvl-1)
    set_parameter_names(len(pols), nvar, [par], vrblvl)
    set_quad_double_start([start], vrblvl)
    set_quad_double_target([target], vrblvl)
    quad_double_real_sweep_run(vrblvl)
    result = get_quad_double_solutions(vrblvl-1)
    return result

def test_double_complex_sweep(vrblvl=0):
    """
    Runs a complex sweep on two points on the unit circle.
    Although we start at two points with real coordinates
    and we end at two points that have nonzero imaginary parts,
    the sweep does not encounter a singularity because of
    the random complex gamma constant.
    """
    if vrblvl > 0:
        print('in test_double_complex_sweep ...')
    circle = ['x^2 + y^2 - 1;']
    first = make_solution(['x', 'y'], [0, 1])
    second = make_solution(['x', 'y'], [0, -1])
    startsols = [first, second]
    xpar = ['x']
    start = [0, 0]  # real and imaginary parts of the start value
    target = [2, 0]
    newsols = double_complex_sweep\
                 (circle, startsols, 2, xpar, start, target, vrblvl)
    if vrblvl > 0:
        for (idx, sol) in enumerate(newsols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(circle, newsols, vrblvl-1)
    if vrblvl > 0:
        print('the error :', err)
    fail = int(err > 1.0e-8)
    return fail

def test_double_double_complex_sweep(vrblvl=0):
    """
    Runs a complex sweep on two points on the unit circle,
    in double double precision.
    Although we start at two points with real coordinates
    and we end at two points that have nonzero imaginary parts,
    the sweep does not encounter a singularity because of
    the random complex gamma constant.
    """
    if vrblvl > 0:
        print('in test_double_double_complex_sweep ...')
    circle = ['x^2 + y^2 - 1;']
    first = make_solution(['x', 'y'], [0, 1])
    second = make_solution(['x', 'y'], [0, -1])
    startsols = [first, second]
    xpar = ['x']
    # real and imaginary parts of the start value
    start = [0.0, 0.0, 0.0, 0.0] # as double doubles
    target = [2.0, 0.0, 0.0, 0.0]
    newsols = double_double_complex_sweep\
                 (circle, startsols, 2, xpar, start, target, vrblvl)
    if vrblvl > 0:
        for (idx, sol) in enumerate(newsols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(circle, newsols, vrblvl-1)
    if vrblvl > 0:
        print('the error :', err)
    fail = int(err > 1.0e-8)
    return fail

def test_quad_double_complex_sweep(vrblvl=0):
    """
    Runs a complex sweep on two points on the unit circle,
    in quad double precision.
    Although we start at two points with real coordinates
    and we end at two points that have nonzero imaginary parts,
    the sweep does not encounter a singularity because of
    the random complex gamma constant.
    """
    if vrblvl > 0:
        print('in test_quad_double_complex_sweep ...')
    circle = ['x^2 + y^2 - 1;']
    first = make_solution(['x', 'y'], [0, 1])
    second = make_solution(['x', 'y'], [0, -1])
    startsols = [first, second]
    xpar = ['x']
    # real and imaginary parts of the start value
    start = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # as quad doubles
    target = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    newsols = quad_double_complex_sweep\
                 (circle, startsols, 2, xpar, start, target, vrblvl)
    if vrblvl > 0:
        for (idx, sol) in enumerate(newsols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(circle, newsols, vrblvl-1)
    if vrblvl > 0:
        print('the error :', err)
    fail = int(err > 1.0e-8)
    return fail

def test_double_real_sweep(vrblvl=0):
    """
    Runs a real sweep on two points on the unit circle: (1,0), (-1,0),
    moving the second coordinate from 0 to 2.
    The sweep will stop at the quadratic turning point: (0,1).
    We can also run the sweep starting at two complex points:
    (2*j, sqrt(5)) and (-2*j, sqrt(5)), moving the second coordinate
    from sqrt(5) to 0.  This sweep will also stop at (0,1).
    """
    if vrblvl > 0:
        print('in test_double_real_sweep ...')
    circle = ['x^2 + y^2 - 1;', 'y*(1-s) + (y-2)*s;']
    first = make_solution(['x', 'y', 's'], [1, 0, 0])
    second = make_solution(['x', 'y', 's'], [-1, 0, 0])
    startsols = [first, second]
    newsols = double_real_sweep(circle, startsols, vrblvl=vrblvl)
    if vrblvl > 0:
        for (idx, sol) in enumerate(newsols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(circle, newsols, vrblvl-1)
    if vrblvl > 0:
        print('the error :', err)
    fail = int(err > 1.0e-8)
    sqrt5 = sqrt(5)
    sweepline = f'(y - {sqrt5:.15e})*(1-s) + y*s;'
    circle = ['x^2 + y^2 - 1;', sweepline]
    first = make_solution(['x', 'y', 's'], [complex(0,2), sqrt5, 0])
    second = make_solution(['x', 'y', 's'], [complex(0,-2), sqrt5, 0])
    startsols = [first, second]
    newsols = double_real_sweep(circle, startsols, vrblvl=vrblvl)
    if vrblvl > 0:
        for (idx, sol) in enumerate(newsols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(circle, newsols, vrblvl-1)
    if vrblvl > 0:
        print('the error :', err)
    fail = fail + int(err > 1.0e-8)
    return fail

def test_double_double_real_sweep(vrblvl=0):
    """
    Runs a real sweep on two points on the unit circle: (1,0), (-1,0),
    moving the second coordinate from 0 to 2.
    The sweep will stop at the quadratic turning point: (0,1).
    We can also run the sweep starting at two complex points:
    (2*j, sqrt(5)) and (-2*j, sqrt(5)), moving the second coordinate
    from sqrt(5) to 0.  This sweep will also stop at (0,1).
    """
    if vrblvl > 0:
        print('in test_double_double_real_sweep ...')
    circle = ['x^2 + y^2 - 1;', 'y*(1-s) + (y-2)*s;']
    first = make_solution(['x', 'y', 's'], [1, 0, 0])
    second = make_solution(['x', 'y', 's'], [-1, 0, 0])
    startsols = [first, second]
    newsols = double_double_real_sweep(circle, startsols, vrblvl=vrblvl)
    if vrblvl > 0:
        for (idx, sol) in enumerate(newsols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(circle, newsols, vrblvl-1)
    if vrblvl > 0:
        print('the error :', err)
    fail = int(err > 1.0e-8)
    sqrt5 = sqrt(5)
    sweepline = f'(y - {sqrt5:.15e})*(1-s) + y*s;'
    circle = ['x^2 + y^2 - 1;', sweepline]
    first = make_solution(['x', 'y', 's'], [complex(0,2), sqrt5, 0])
    second = make_solution(['x', 'y', 's'], [complex(0,-2), sqrt5, 0])
    startsols = [first, second]
    newsols = double_double_real_sweep(circle, startsols, vrblvl=vrblvl)
    if vrblvl > 0:
        for (idx, sol) in enumerate(newsols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(circle, newsols, vrblvl-1)
    if vrblvl > 0:
        print('the error :', err)
    fail = fail + int(err > 1.0e-8)
    return fail

def test_quad_double_real_sweep(vrblvl=0):
    """
    Runs a real sweep on two points on the unit circle: (1,0), (-1,0),
    moving the second coordinate from 0 to 2.
    The sweep will stop at the quadratic turning point: (0,1).
    We can also run the sweep starting at two complex points:
    (2*j, sqrt(5)) and (-2*j, sqrt(5)), moving the second coordinate
    from sqrt(5) to 0.  This sweep will also stop at (0,1).
    """
    if vrblvl > 0:
        print('in test_quad_double_real_sweep ...')
    circle = ['x^2 + y^2 - 1;', 'y*(1-s) + (y-2)*s;']
    first = make_solution(['x', 'y', 's'], [1, 0, 0])
    second = make_solution(['x', 'y', 's'], [-1, 0, 0])
    startsols = [first, second]
    newsols = quad_double_real_sweep(circle, startsols, vrblvl=vrblvl)
    if vrblvl > 0:
        for (idx, sol) in enumerate(newsols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(circle, newsols, vrblvl-1)
    if vrblvl > 0:
        print('the error :', err)
    fail = int(err > 1.0e-8)
    sqrt5 = sqrt(5)
    sweepline = f'(y - {sqrt5:.15e})*(1-s) + y*s;'
    circle = ['x^2 + y^2 - 1;', sweepline]
    first = make_solution(['x', 'y', 's'], [complex(0,2), sqrt5, 0])
    second = make_solution(['x', 'y', 's'], [complex(0,-2), sqrt5, 0])
    startsols = [first, second]
    newsols = quad_double_real_sweep(circle, startsols, vrblvl=vrblvl)
    if vrblvl > 0:
        for (idx, sol) in enumerate(newsols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(circle, newsols, vrblvl-1)
    if vrblvl > 0:
        print('the error :', err)
    fail = fail + int(err > 1.0e-8)
    return fail

def main():
    """
    Runs some tests.
    """
    lvl = 1
    fail = test_double_complex_sweep(lvl)
    fail = test_double_double_complex_sweep(lvl)
    fail = test_quad_double_complex_sweep(lvl)
    fail = test_double_real_sweep(lvl)
    fail = test_double_double_real_sweep(lvl)
    fail = test_quad_double_real_sweep(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__== '__main__':
    main()
