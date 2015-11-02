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
"""

def standard_sweep(pols, sols, pars, start, target):
    """
    For the polynomials in the list of strings pols
    and the solutions in sols for the values in the list start,
    a sweep through the parameter space will be performed
    in standard double precision to the target values of
    the parameters in the list target.
    The list of symbols in pars contains the names of the variables
    in the polynomials pols that serve as parameters.
    The size of the lists pars, start, and target must be same.
    """
    from phcpy.interface import store_standard_solutions as storesols
    from phcpy.interface import store_standard_system as storesys
    storesys(pols)
    storesols(len(pols), sols)
    from phcpy.interface import load_standard_solutions as loadsols
    from phcpy.phcpy2c \
    import py2c_sweep_define_parameters_symbolically as define
    from phcpy.phcpy2c \
    import py2c_sweep_set_standard_start as set_start
    from phcpy.phcpy2c \
    import py2c_sweep_set_standard_target as set_target
    from phcpy.phcpy2c import py2c_sweep_standard_run as run
    (nbq, nbv, nbp) = (len(pols), len(pols), len(pars))
    parnames = ' '.join(pars)
    nbc = len(parnames)
    print 'defining the parameters in the sweep homotopy ...'
    define(nbq, nbv, nbp, nbc, parnames)  
    print 'setting the start and the target ...'
    set_start(nbp, str(start));
    set_target(nbp, str(target));
    print 'calling run in standard double precision ...'
    run(0, 0.0, 0.0);
    result = loadsols()
    return result

def dobldobl_sweep(pols, sols, pars, start, target):
    """
    For the polynomials in the list of strings pols
    and the solutions in sols for the values in the list start,
    a sweep through the parameter space will be performed
    in double double precision to the target values of
    the parameters in the list target.
    The list of symbols in pars contains the names of the variables
    in the polynomials pols that serve as parameters.
    The size of the lists pars, start, and target must be same.
    """
    from phcpy.interface import store_dobldobl_solutions as storesols
    from phcpy.interface import store_dobldobl_system as storesys
    storesys(pols)
    storesols(len(pols), sols)
    from phcpy.interface import load_dobldobl_solutions as loadsols
    from phcpy.phcpy2c \
    import py2c_sweep_define_parameters_symbolically as define
    from phcpy.phcpy2c \
    import py2c_sweep_set_dobldobl_start as set_start
    from phcpy.phcpy2c \
    import py2c_sweep_set_dobldobl_target as set_target
    from phcpy.phcpy2c import py2c_sweep_dobldobl_run as run
    (nbq, nbv, nbp) = (len(pols), len(pols), len(pars))
    parnames = ' '.join(pars)
    nbc = len(parnames)
    print 'defining the parameters in the sweep homotopy ...'
    define(nbq, nbv, nbp, nbc, parnames)  
    print 'setting the start and the target ...'
    set_start(nbp, str(start));
    set_target(nbp, str(target));
    print 'calling run in double double precision ...'
    run(0, 0.0, 0.0);
    result = loadsols()
    return result

def quaddobl_sweep(pols, sols, pars, start, target):
    """
    For the polynomials in the list of strings pols
    and the solutions in sols for the values in the list start,
    a sweep through the parameter space will be performed
    in quad double precision to the target values of
    the parameters in the list target.
    The list of symbols in pars contains the names of the variables
    in the polynomials pols that serve as parameters.
    The size of the lists pars, start, and target must be same.
    """
    from phcpy.interface import store_quaddobl_solutions as storesols
    from phcpy.interface import store_quaddobl_system as storesys
    storesys(pols)
    storesols(len(pols), sols)
    from phcpy.interface import load_quaddobl_solutions as loadsols
    from phcpy.phcpy2c \
    import py2c_sweep_define_parameters_symbolically as define
    from phcpy.phcpy2c \
    import py2c_sweep_set_quaddobl_start as set_start
    from phcpy.phcpy2c \
    import py2c_sweep_set_quaddobl_target as set_target
    from phcpy.phcpy2c import py2c_sweep_quaddobl_run as run
    (nbq, nbv, nbp) = (len(pols), len(pols), len(pars))
    parnames = ' '.join(pars)
    nbc = len(parnames)
    print 'defining the parameters in the sweep homotopy ...'
    define(nbq, nbv, nbp, nbc, parnames)  
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
    circle = ['x^2 + y^2 - 1;', 'x;']  # this is a bug ...
    first = makesol(['x', 'y'], [0, 1])
    second = makesol(['x', 'y'], [0, -1])
    startsols = [first, second]
    xpar = ['x']
    if(precision == 'd'):
        stxstart = [0, 0]  # real and imaginary parts of the start value
        stxtarget = [2, 0]
        newsols = standard_sweep(circle, startsols, xpar, stxstart, stxtarget)
    elif(precision == 'dd'):
        ddxstart = [0, 0, 0, 0]  # double doubles
        ddxtarget = [2, 0, 0, 0] 
        newsols = dobldobl_sweep(circle, startsols, xpar, ddxstart, ddxtarget)
    elif(precision == 'qd'):
        qdxstart = [0, 0, 0, 0, 0, 0, 0, 0]  # quad doubles
        qdxtarget = [2, 0, 0, 0, 0, 0, 0, 0]
        newsols = quaddobl_sweep(circle, startsols, xpar, qdxstart, qdxtarget)
    else:
        print 'wrong precision given as input parameter to test'
    for sol in newsols:
        print sol

if __name__ == "__main__":
    test('qd')
