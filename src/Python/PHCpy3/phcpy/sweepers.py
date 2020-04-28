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

def standard_complex_sweep(pols, sols, nvar, pars, start, target):
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
    from phcpy.interface import store_standard_solutions as storesols
    from phcpy.interface import store_standard_system as storesys
    storesys(pols, nbvar=nvar)
    storesols(nvar, sols)
    from phcpy.interface import load_standard_solutions as loadsols
    from phcpy.phcpy2c3 \
    import py2c_sweep_define_parameters_symbolically as define
    from phcpy.phcpy2c3 \
    import py2c_sweep_set_standard_start as set_start
    from phcpy.phcpy2c3 \
    import py2c_sweep_set_standard_target as set_target
    from phcpy.phcpy2c3 import py2c_sweep_standard_complex_run as run
    (nbq, nbp) = (len(pols), len(pars))
    parnames = ' '.join(pars)
    nbc = len(parnames)
    define(nbq, nvar, nbp, nbc, parnames)
    print('setting the start and the target ...')
    set_start(nbp, str(start))
    set_target(nbp, str(target))
    print('calling run in standard double precision ...')
    run(0, 0.0, 0.0)
    result = loadsols()
    return result

def dobldobl_complex_sweep(pols, sols, nvar, pars, start, target):
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
    The size of the lists *pars, *start*, and *target* must be same.
    """
    from phcpy.interface import store_dobldobl_solutions as storesols
    from phcpy.interface import store_dobldobl_system as storesys
    storesys(pols, nbvar=nvar)
    storesols(nvar, sols)
    from phcpy.interface import load_dobldobl_solutions as loadsols
    from phcpy.phcpy2c3 \
    import py2c_sweep_define_parameters_symbolically as define
    from phcpy.phcpy2c3 \
    import py2c_sweep_set_dobldobl_start as set_start
    from phcpy.phcpy2c3 \
    import py2c_sweep_set_dobldobl_target as set_target
    from phcpy.phcpy2c3 import py2c_sweep_dobldobl_complex_run as run
    (nbq, nbp) = (len(pols), len(pars))
    parnames = ' '.join(pars)
    nbc = len(parnames)
    define(nbq, nvar, nbp, nbc, parnames)
    print('setting the start and the target ...')
    set_start(nbp, str(start))
    set_target(nbp, str(target))
    print('calling run in double double precision ...')
    run(0, 0.0, 0.0)
    result = loadsols()
    return result

def quaddobl_complex_sweep(pols, sols, nvar, pars, start, target):
    r"""
    For the polynomials in the list of strings p*ols*
    and the solutions in *sols* for the values in the list *start*,
    a sweep through the parameter space will be performed
    in quad double precision to the target values of
    the parameters in the list *target*.
    The number of variables in the polynomials and the solutions
    must be the same and be equal to the value of *nvar*.
    The list of symbols in *pars* contains the names of the variables
    in the polynomials pols that serve as parameters.
    The size of the lists *pars*, *start*, and *target* must be same.
    """
    from phcpy.interface import store_quaddobl_solutions as storesols
    from phcpy.interface import store_quaddobl_system as storesys
    storesys(pols, nbvar=nvar)
    storesols(nvar, sols)
    from phcpy.interface import load_quaddobl_solutions as loadsols
    from phcpy.phcpy2c3 \
    import py2c_sweep_define_parameters_symbolically as define
    from phcpy.phcpy2c3 \
    import py2c_sweep_set_quaddobl_start as set_start
    from phcpy.phcpy2c3 \
    import py2c_sweep_set_quaddobl_target as set_target
    from phcpy.phcpy2c3 import py2c_sweep_quaddobl_complex_run as run
    (nbq, nbp) = (len(pols), len(pars))
    parnames = ' '.join(pars)
    nbc = len(parnames)
    define(nbq, nvar, nbp, nbc, parnames)
    print('setting the start and the target ...')
    set_start(nbp, str(start))
    set_target(nbp, str(target))
    print('calling run in quad double precision ...')
    run(0, 0.0, 0.0)
    result = loadsols()
    return result

def standard_real_sweep(pols, sols, par='s', start=0.0, target=1.0):
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
    Computations happend in standard double precision.
    """
    from phcpy.interface import store_standard_solutions as storesols
    from phcpy.interface import store_standard_system as storesys
    nvar = len(pols) + 1
    storesys(pols, nbvar=nvar)
    storesols(nvar, sols)
    from phcpy.interface import load_standard_solutions as loadsols
    from phcpy.phcpy2c3 \
        import py2c_sweep_define_parameters_symbolically as define
    from phcpy.phcpy2c \
    import py2c_sweep_set_standard_start as set_start
    from phcpy.phcpy2c3 \
        import py2c_sweep_set_standard_target as set_target
    (nbq, nbp) = (len(pols), 1)
    pars = [par]
    parnames = ' '.join(pars)
    nbc = len(parnames)
    define(nbq, nvar, nbp, nbc, parnames)
    set_start(nbp, str([start, 0.0]))
    set_target(nbp, str([target, 0.0]))
    from phcpy.phcpy2c3 import py2c_sweep_standard_real_run as run
    run()
    result = loadsols()
    return result

def dobldobl_real_sweep(pols, sols, par='s', start=0.0, target=1.0):
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
    from phcpy.interface import store_dobldobl_solutions as storesols
    from phcpy.interface import store_dobldobl_system as storesys
    nvar = len(pols) + 1
    storesys(pols, nbvar=nvar)
    storesols(nvar, sols)
    # print('done storing system and solutions ...')
    from phcpy.interface import load_dobldobl_solutions as loadsols
    from phcpy.phcpy2c3 \
        import py2c_sweep_define_parameters_symbolically as define
    from phcpy.phcpy2c3 \
        import py2c_sweep_set_dobldobl_start as set_start
    from phcpy.phcpy2c3 \
        import py2c_sweep_set_dobldobl_target as set_target
    (nbq, nbp) = (len(pols), 1)
    pars = [par]
    parnames = ' '.join(pars)
    nbc = len(parnames)
    # print('defining the parameters ...')
    define(nbq, nvar, nbp, nbc, parnames)
    set_start(nbp, str([start, 0.0, 0.0, 0.0]))  # double doubles !
    set_target(nbp, str([target, 0.0, 0.0, 0.0]))
    from phcpy.phcpy2c3 import py2c_sweep_dobldobl_real_run as run
    run()
    result = loadsols()
    return result

def quaddobl_real_sweep(pols, sols, par='s', start=0.0, target=1.0):
    r"""
    A real sweep homotopy is a family of n equations in n+1 variables,
    where one of the variables is the artificial parameter s which moves
    from 0.0 to 1.0.  The last equation can then be of the form

    (1 - s)*(lambda - L[0]) + s*(lambda - L[1]) = 0 so that,

    at s = 0, the natural parameter lambda has the value L[0], and

    at s = 1, the natural parameter lambda has the value L[1].

    Thus: as s moves from 0 to 1, lambda goes from L[0] to L[1].
    All solutions in the list sols must have then the value L[0]
    for the variable lambda.
    The sweep stops when the target value for s is reached
    or when a singular solution is encountered.
    Computations happen in quad double precision.
    """
    from phcpy.interface import store_quaddobl_solutions as storesols
    from phcpy.interface import store_quaddobl_system as storesys
    nvar = len(pols) + 1
    storesys(pols, nbvar=nvar)
    storesols(nvar, sols)
    # print('done storing system and solutions ...')
    from phcpy.interface import load_quaddobl_solutions as loadsols
    from phcpy.phcpy2c3 \
        import py2c_sweep_define_parameters_symbolically as define
    from phcpy.phcpy2c3 \
        import py2c_sweep_set_quaddobl_start as set_start
    from phcpy.phcpy2c3 \
    import py2c_sweep_set_quaddobl_target as set_target
    (nbq, nbp) = (len(pols), 1)
    pars = [par]
    parnames = ' '.join(pars)
    nbc = len(parnames)
    # print('defining the parameters ...')
    define(nbq, nvar, nbp, nbc, parnames)
    set_start(nbp, str([start, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
    set_target(nbp, str([target, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
    from phcpy.phcpy2c3 import py2c_sweep_quaddobl_real_run as run
    run()
    result = loadsols()
    return result

def complex_sweep_test(precision='d'):
    """
    Runs a complex sweep on two points on the unit circle.
    Although we start at two points with real coordinates
    and we end at two points that have nonzero imaginary parts,
    the sweep does not encounter a singularity because of
    the random complex gamma constant.
    """
    from phcpy.solutions import make_solution as makesol
    circle = ['x^2 + y^2 - 1;']
    first = makesol(['x', 'y'], [0, 1])
    second = makesol(['x', 'y'], [0, -1])
    startsols = [first, second]
    xpar = ['x']
    if(precision == 'd'):
        ststart = [0, 0]  # real and imaginary parts of the start value
        sttarget = [2, 0]
        newsols = standard_complex_sweep(circle, startsols, 2, xpar, \
                                         ststart, sttarget)
    elif(precision == 'dd'):
        ddstart = [0, 0, 0, 0]  # double doubles
        ddtarget = [2, 0, 0, 0]
        newsols = dobldobl_complex_sweep(circle, startsols, 2, xpar, \
                                         ddstart, ddtarget)
    elif(precision == 'qd'):
        qdstart = [0, 0, 0, 0, 0, 0, 0, 0]  # quad doubles
        qdtarget = [2, 0, 0, 0, 0, 0, 0, 0]
        newsols = quaddobl_complex_sweep(circle, startsols, 2, xpar, \
                                         qdstart, qdtarget)
    else:
        print('wrong precision given as input parameter to test')
    for sol in newsols:
        print(sol)

def real_sweep_test(precision='d'):
    """
    Runs a real sweep on two points on the unit circle: (1,0), (-1,0),
    moving the second coordinate from 0 to 2.
    The sweep will stop at the quadratic turning point: (0,1).
    We can also run the sweep starting at two complex points:
    (2*j, sqrt(5)) and (-2*j, sqrt(5)), moving the second coordinate
    from sqrt(5) to 0.  This sweep will also stop at (0,1).
    """
    from phcpy.solutions import make_solution as makesol
    rcircle = ['x^2 + y^2 - 1;', 'y*(1-s) + (y-2)*s;']
    rfirst = makesol(['x', 'y', 's'], [1, 0, 0])
    rsecond = makesol(['x', 'y', 's'], [-1, 0, 0])
    rstartsols = [rfirst, rsecond]
    if(precision == 'd'):
       rnewsols = standard_real_sweep(rcircle, rstartsols)
    elif(precision == 'dd'):
       rnewsols = dobldobl_real_sweep(rcircle, rstartsols)
    elif(precision == 'qd'):
       rnewsols = quaddobl_real_sweep(rcircle, rstartsols)
    else:
       print('wrong precision given as input parameter to test')
    print('after the sweep that started at real solutions :')
    for sol in rnewsols:
        print(sol)
    from math import sqrt
    sqrt5 = sqrt(5)
    sweepline = '(y - %.15e)*(1-s) + y*s;' % sqrt5
    ccircle = ['x^2 + y^2 - 1;', sweepline]
    cfirst = makesol(['x', 'y', 's'], [complex(0,2), sqrt5, 0])
    csecond = makesol(['x', 'y', 's'], [complex(0,-2), sqrt5, 0])
    cstartsols = [cfirst, csecond]
    if(precision == 'd'):
        cnewsols = standard_real_sweep(ccircle, cstartsols)
    elif(precision == 'dd'):
        cnewsols = dobldobl_real_sweep(ccircle, cstartsols)
    elif(precision == 'qd'):
        cnewsols = quaddobl_real_sweep(ccircle, cstartsols)
    else:
        print('wrong precision given as input parameter to test')
    print('after the sweep that started at complex solutions :')
    for sol in cnewsols:
        print(sol)

if __name__ == "__main__":
    real_sweep_test('qd')
    #complex_sweep_test('d')
    #complex_sweep_test('dd')
    #complex_sweep_test('qd')
