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
The algorithms applied in this module are described in the paper by
Kathy Piret and Jan Verschelde: Sweeping Algebraic Curves for Singular 
Solutions.  Journal of Computational and Applied Mathematics,
volume 234, number 4, pages 1228-1237, 2010. 
"""
from ctypes import c_int32, c_double, pointer
from phcpy.version import get_phcfun
from phcpy.polynomials import set_double_system
from phcpy.solutions import set_double_solutions, get_double_solutions
from phcpy.solutions import make_solution

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
    aprc = pointer(c_int32(0))
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

def double_complex_sweep(pols, sols, nvar, pars, start, target, vrblvl=0):
    r"""
    For the polynomials in the list of strings *pols*
    and the solutions in *sols* for the values in the list *start*,
    a sweep through the parameter space will be performed
    in standard double precision to the target values of
    the parameters in the list *target*.
    The number of variables in the polynomials and the solutions
    must be the same and be equal to the value of *nvar*.
    The list of symbols in *pars* contains the names of the variables
    in the polynomials *pols* that serve as parameters.
    The size of the lists *pars*, *start*, and *target* must be same.
    """
    if vrblvl > 0:
        print('in double_complex_sweep, nvar :', nvar)
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

def test_double_complex_sweep(vrblvl=0):
    """
    Runs a complex sweep on two points on the unit circle.
    Although we start at two points with real coordinates
    and we end at two points that have nonzero imaginary parts,
    the sweep does not encounter a singularity because of
    the random complex gamma constant.
    """
    circle = ['x^2 + y^2 - 1;']
    first = make_solution(['x', 'y'], [0, 1])
    second = make_solution(['x', 'y'], [0, -1])
    startsols = [first, second]
    xpar = ['x']
    ststart = [0, 0]  # real and imaginary parts of the start value
    sttarget = [2, 0]
    newsols = double_complex_sweep\
                 (circle, startsols, 2, xpar, ststart, sttarget, vrblvl)
    if vrblvl > 0:
        for (idx, sol) in enumerate(newsols):
            print('Solution', idx+1, ':')
            print(sol)
    return 0

def main():
    """
    Runs some tests.
    """
    lvl = 1
    fail = test_double_complex_sweep(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__== '__main__':
    main()
