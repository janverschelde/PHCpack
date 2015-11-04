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
"""

def standard_complex_sweep(pols, sols, nvar, pars, start, target):
    """
    For the polynomials in the list of strings pols
    and the solutions in sols for the values in the list start,
    a sweep through the parameter space will be performed
    in standard double precision to the target values of
    the parameters in the list target.
    The number of variables in the polynomials and the solutions
    must be the same and be equal to the value of nvar.
    The list of symbols in pars contains the names of the variables
    in the polynomials pols that serve as parameters.
    The size of the lists pars, start, and target must be same.
    """
    from phcpy.interface import store_standard_solutions as storesols
    from phcpy.interface import store_standard_system as storesys
    storesys(pols, nbvar=nvar)
    storesols(nvar, sols)
    from phcpy.interface import load_standard_solutions as loadsols
    from phcpy.phcpy2c \
    import py2c_sweep_define_parameters_symbolically as define
    from phcpy.phcpy2c \
    import py2c_sweep_set_standard_start as set_start
    from phcpy.phcpy2c \
    import py2c_sweep_set_standard_target as set_target
    from phcpy.phcpy2c import py2c_sweep_standard_complex_run as run
    (nbq, nbp) = (len(pols), len(pars))
    parnames = ' '.join(pars)
    nbc = len(parnames)
    print 'defining the parameters in the sweep homotopy ...'
    define(nbq, nvar, nbp, nbc, parnames)  
    print 'setting the start and the target ...'
    set_start(nbp, str(start));
    set_target(nbp, str(target));
    print 'calling run in standard double precision ...'
    run(0, 0.0, 0.0);
    result = loadsols()
    return result

def dobldobl_complex_sweep(pols, sols, nvar, pars, start, target):
    """
    For the polynomials in the list of strings pols
    and the solutions in sols for the values in the list start,
    a sweep through the parameter space will be performed
    in double double precision to the target values of
    the parameters in the list target.
    The number of variables in the polynomials and the solutions
    must be the same and be equal to the value of nvar.
    The list of symbols in pars contains the names of the variables
    in the polynomials pols that serve as parameters.
    The size of the lists pars, start, and target must be same.
    """
    from phcpy.interface import store_dobldobl_solutions as storesols
    from phcpy.interface import store_dobldobl_system as storesys
    storesys(pols, nbvar=nvar)
    storesols(nvar, sols)
    from phcpy.interface import load_dobldobl_solutions as loadsols
    from phcpy.phcpy2c \
    import py2c_sweep_define_parameters_symbolically as define
    from phcpy.phcpy2c \
    import py2c_sweep_set_dobldobl_start as set_start
    from phcpy.phcpy2c \
    import py2c_sweep_set_dobldobl_target as set_target
    from phcpy.phcpy2c import py2c_sweep_dobldobl_complex_run as run
    (nbq, nbp) = (len(pols), len(pars))
    parnames = ' '.join(pars)
    nbc = len(parnames)
    print 'defining the parameters in the sweep homotopy ...'
    define(nbq, nvar, nbp, nbc, parnames)  
    print 'setting the start and the target ...'
    set_start(nbp, str(start));
    set_target(nbp, str(target));
    print 'calling run in double double precision ...'
    run(0, 0.0, 0.0);
    result = loadsols()
    return result

def quaddobl_complex_sweep(pols, sols, nvar, pars, start, target):
    """
    For the polynomials in the list of strings pols
    and the solutions in sols for the values in the list start,
    a sweep through the parameter space will be performed
    in quad double precision to the target values of
    the parameters in the list target.
    The number of variables in the polynomials and the solutions
    must be the same and be equal to the value of nvar.
    The list of symbols in pars contains the names of the variables
    in the polynomials pols that serve as parameters.
    The size of the lists pars, start, and target must be same.
    """
    from phcpy.interface import store_quaddobl_solutions as storesols
    from phcpy.interface import store_quaddobl_system as storesys
    storesys(pols, nbvar=nvar)
    storesols(nvar, sols)
    from phcpy.interface import load_quaddobl_solutions as loadsols
    from phcpy.phcpy2c \
    import py2c_sweep_define_parameters_symbolically as define
    from phcpy.phcpy2c \
    import py2c_sweep_set_quaddobl_start as set_start
    from phcpy.phcpy2c \
    import py2c_sweep_set_quaddobl_target as set_target
    from phcpy.phcpy2c import py2c_sweep_quaddobl_complex_run as run
    (nbq, nbp) = (len(pols), len(pars))
    parnames = ' '.join(pars)
    nbc = len(parnames)
    print 'defining the parameters in the sweep homotopy ...'
    define(nbq, nvar, nbp, nbc, parnames)  
    print 'setting the start and the target ...'
    set_start(nbp, str(start));
    set_target(nbp, str(target));
    print 'calling run in quad double precision ...'
    run(0, 0.0, 0.0);
    result = loadsols()
    return result

def test(precision='d'):
    """
    Runs a sweep on two points on the unit circle.
    """
    from solutions import make_solution as makesol
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
        print 'wrong precision given as input parameter to test'
    for sol in newsols:
        print sol

if __name__ == "__main__":
    test('d')
    test('dd')
    test('qd')
