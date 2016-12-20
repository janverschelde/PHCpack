"""
The module series exports functions to compute power series solutions with
Newton's method in double, double double, or quad double precision.
"""

def standard_newton_series(pols, sols, idx=1, nbr=4, verbose=True):
    r"""
    Computes series in standard double precision for the polynomials
    in *pols*, where the leading coefficients are the solutions in *sols*.
    On entry are the following five parameters:

    *pols*: a list of string representations of polynomials,

    *sols*: a list of solutions of the polynomials in *pols*,

    *idx*: index of the series parameter, by default equals 1,

    *nbr*: number of steps with Newton's method,

    *verbose*: whether to write intermediate output to screen or not.

    On return is a list of lists of strings.  Each lists of strings
    represents the series solution for the variables in the list *pols*.
    """
    from phcpy.solver import number_of_symbols
    from phcpy.interface import store_standard_solutions
    from phcpy.interface import store_standard_system, load_standard_system
    from phcpy.phcpy2c2 import py2c_standard_Newton_series as newton
    from phcpy.phcpy2c2 import py2c_syspool_standard_size as poolsize
    from phcpy.phcpy2c2 import py2c_syspool_copy_to_standard_container
    from phcpy.phcpy2c2 import py2c_syspool_standard_clear
    nbsym = number_of_symbols(pols)
    if verbose:
        print "the polynomials :"
        for pol in pols:
            print pol
        print "Number of variables :", nbsym
    store_standard_system(pols, nbvar=nbsym)
    store_standard_solutions(nbsym, sols)
    fail = newton(idx, nbr, int(verbose))
    size = (-1 if fail else poolsize())
    if verbose:
        if size == -1:
            print "An error occurred in the execution of Newton's method."
        else:
            print "Computed %d series solutions." % size
    result = []
    for k in range(1, size+1):
        py2c_syspool_copy_to_standard_container(k)
        sersol = load_standard_system()
        result.append(sersol)
    py2c_syspool_standard_clear()
    return result

def dobldobl_newton_series(pols, sols, idx=1, nbr=4, verbose=True):
    r"""
    Computes series in double double precision for the polynomials
    in *pols*, where the leading coefficients are the solutions in *sols*.
    On entry are the following five parameters:

    *pols*: a list of string representations of polynomials,

    *sols*: a list of solutions of the polynomials in *pols*,

    *idx*: index of the series parameter, by default equals 1,

    *nbr*: number of steps with Newton's method,

    *verbose*: whether to write intermediate output to screen or not.

    On return is a list of lists of strings.  Each lists of strings
    represents the series solution for the variables in the list *pols*.
    """
    from phcpy.solver import number_of_symbols
    from phcpy.interface import store_dobldobl_solutions
    from phcpy.interface import store_dobldobl_system, load_dobldobl_system
    from phcpy.phcpy2c2 import py2c_dobldobl_Newton_series as newton
    from phcpy.phcpy2c2 import py2c_syspool_dobldobl_size as poolsize
    from phcpy.phcpy2c2 import py2c_syspool_copy_to_dobldobl_container
    from phcpy.phcpy2c2 import py2c_syspool_dobldobl_clear
    nbsym = number_of_symbols(pols)
    if verbose:
        print "the polynomials :"
        for pol in pols:
            print pol
        print "Number of variables :", nbsym
    store_dobldobl_system(pols, nbvar=nbsym)
    store_dobldobl_solutions(nbsym, sols)
    fail = newton(idx, nbr, int(verbose))
    size = (-1 if fail else poolsize())
    if verbose:
        if size == -1:
            print "An error occurred in the execution of Newton's method."
        else:
            print "Computed %d series solutions." % size
    result = []
    for k in range(1, size+1):
        py2c_syspool_copy_to_dobldobl_container(k)
        sersol = load_dobldobl_system()
        result.append(sersol)
    py2c_syspool_dobldobl_clear()
    return result

def quaddobl_newton_series(pols, sols, idx=1, nbr=4, verbose=True):
    r"""
    Computes series in quad double precision for the polynomials
    in *pols*, where the leading coefficients are the solutions in *sols*.
    On entry are the following five parameters:

    *pols*: a list of string representations of polynomials,

    *sols*: a list of solutions of the polynomials in *pols*,

    *idx*: index of the series parameter, by default equals 1,

    *nbr*: number of steps with Newton's method,

    *verbose*: whether to write intermediate output to screen or not.

    On return is a list of lists of strings.  Each lists of strings
    represents the series solution for the variables in the list *pols*.
    """
    from phcpy.solver import number_of_symbols
    from phcpy.interface import store_quaddobl_solutions
    from phcpy.interface import store_quaddobl_system, load_quaddobl_system
    from phcpy.phcpy2c2 import py2c_quaddobl_Newton_series as newton
    from phcpy.phcpy2c2 import py2c_syspool_quaddobl_size as poolsize
    from phcpy.phcpy2c2 import py2c_syspool_copy_to_quaddobl_container
    from phcpy.phcpy2c2 import py2c_syspool_quaddobl_clear
    nbsym = number_of_symbols(pols)
    if verbose:
        print "the polynomials :"
        for pol in pols:
            print pol
        print "Number of variables :", nbsym
    store_quaddobl_system(pols, nbvar=nbsym)
    store_quaddobl_solutions(nbsym, sols)
    fail = newton(idx, nbr, int(verbose))
    size = (-1 if fail else poolsize())
    if verbose:
        if size == -1:
            print "An error occurred in the execution of Newton's method."
        else:
            print "Computed %d series solutions." % size
    result = []
    for k in range(1, size+1):
        py2c_syspool_copy_to_quaddobl_container(k)
        sersol = load_quaddobl_system()
        result.append(sersol)
    py2c_syspool_quaddobl_clear()
    return result

def standard_newton_power_series(pols, lser, idx=1, nbr=4, verbose=True):
    r"""
    Computes series in standard double precision for the polynomials
    in *pols*, where the leading terms are given in the list *lser*.
    On entry are the following five parameters:

    *pols*: a list of string representations of polynomials,

    *lser*: a list of polynomials in the series parameter (e.g.: t),
    for use as start terms in Newton's method,

    *idx*: index of the series parameter, by default equals 1,

    *nbr*: number of steps with Newton's method,

    *verbose*: whether to write intermediate output to screen or not.

    On return is a list of lists of strings.  Each lists of strings
    represents the series solution for the variables in the list *pols*.
    """
    from phcpy.solver import number_of_symbols
    from phcpy.interface import store_standard_system, load_standard_system
    from phcpy.phcpy2c2 import py2c_standard_Newton_power_series as newton
    from phcpy.phcpy2c2 import py2c_syspool_standard_init
    from phcpy.phcpy2c2 import py2c_syspool_standard_create
    from phcpy.phcpy2c2 import py2c_syspool_standard_size as poolsize
    from phcpy.phcpy2c2 import py2c_syspool_copy_to_standard_container
    from phcpy.phcpy2c2 import py2c_syspool_standard_clear
    nbsym = number_of_symbols(pols)
    if verbose:
        print "the polynomials :"
        for pol in pols:
            print pol
        print "Number of variables :", nbsym
    store_standard_system(lser, nbvar=1)
    py2c_syspool_standard_init(1);
    py2c_syspool_standard_create(1);
    store_standard_system(pols, nbvar=nbsym)
    fail = newton(idx, nbr, int(verbose))
    size = (-1 if fail else poolsize())
    if verbose:
        if size == -1:
            print "An error occurred in the execution of Newton's method."
        else:
            print "Computed one series solution."
    py2c_syspool_copy_to_standard_container(1)
    result = load_standard_system()
    py2c_syspool_standard_clear()
    return result

def dobldobl_newton_power_series(pols, lser, idx=1, nbr=4, verbose=True):
    r"""
    Computes series in double double precision for the polynomials
    in *pols*, where the leading terms are given in the list *lser*.
    On entry are the following five parameters:

    *pols*: a list of string representations of polynomials,

    *lser*: a list of polynomials in the series parameter (e.g.: t),
    for use as start terms in Newton's method,

    *idx*: index of the series parameter, by default equals 1,

    *nbr*: number of steps with Newton's method,

    *verbose*: whether to write intermediate output to screen or not.

    On return is a list of lists of strings.  Each lists of strings
    represents the series solution for the variables in the list *pols*.
    """
    from phcpy.solver import number_of_symbols
    from phcpy.interface import store_dobldobl_system, load_dobldobl_system
    from phcpy.phcpy2c2 import py2c_dobldobl_Newton_power_series as newton
    from phcpy.phcpy2c2 import py2c_syspool_dobldobl_init
    from phcpy.phcpy2c2 import py2c_syspool_dobldobl_create
    from phcpy.phcpy2c2 import py2c_syspool_dobldobl_size as poolsize
    from phcpy.phcpy2c2 import py2c_syspool_copy_to_dobldobl_container
    from phcpy.phcpy2c2 import py2c_syspool_dobldobl_clear
    nbsym = number_of_symbols(pols)
    if verbose:
        print "the polynomials :"
        for pol in pols:
            print pol
        print "Number of variables :", nbsym
    store_dobldobl_system(lser, nbvar=1)
    py2c_syspool_dobldobl_init(1);
    py2c_syspool_dobldobl_create(1);
    store_dobldobl_system(pols, nbvar=nbsym)
    fail = newton(idx, nbr, int(verbose))
    size = (-1 if fail else poolsize())
    if verbose:
        if size == -1:
            print "An error occurred in the execution of Newton's method."
        else:
            print "Computed one series solution."
    py2c_syspool_copy_to_dobldobl_container(1)
    result = load_dobldobl_system()
    py2c_syspool_dobldobl_clear()
    return result

def quaddobl_newton_power_series(pols, lser, idx=1, nbr=4, verbose=True):
    r"""
    Computes series in quad double precision for the polynomials
    in *pols*, where the leading terms are given in the list *lser*.
    On entry are the following five parameters:

    *pols*: a list of string representations of polynomials,

    *lser*: a list of polynomials in the series parameter (e.g.: t),
    for use as start terms in Newton's method,

    *idx*: index of the series parameter, by default equals 1,

    *nbr*: number of steps with Newton's method,

    *verbose*: whether to write intermediate output to screen or not.

    On return is a list of lists of strings.  Each lists of strings
    represents the series solution for the variables in the list *pols*.
    """
    from phcpy.solver import number_of_symbols
    from phcpy.interface import store_quaddobl_system, load_quaddobl_system
    from phcpy.phcpy2c2 import py2c_quaddobl_Newton_power_series as newton
    from phcpy.phcpy2c2 import py2c_syspool_quaddobl_init
    from phcpy.phcpy2c2 import py2c_syspool_quaddobl_create
    from phcpy.phcpy2c2 import py2c_syspool_quaddobl_size as poolsize
    from phcpy.phcpy2c2 import py2c_syspool_copy_to_quaddobl_container
    from phcpy.phcpy2c2 import py2c_syspool_quaddobl_clear
    nbsym = number_of_symbols(pols)
    if verbose:
        print "the polynomials :"
        for pol in pols:
            print pol
        print "Number of variables :", nbsym
    store_quaddobl_system(lser, nbvar=1)
    py2c_syspool_quaddobl_init(1);
    py2c_syspool_quaddobl_create(1);
    store_quaddobl_system(pols, nbvar=nbsym)
    fail = newton(idx, nbr, int(verbose))
    size = (-1 if fail else poolsize())
    if verbose:
        if size == -1:
            print "An error occurred in the execution of Newton's method."
        else:
            print "Computed one series solution."
    py2c_syspool_copy_to_quaddobl_container(1)
    result = load_quaddobl_system()
    py2c_syspool_quaddobl_clear()
    return result

def viviani(prc='d'):
    """
    Returns the system which stores the Viviani curve,
    with some solutions intersected with a plane,
    in double ('d'), double double ('dd'), or quad double('qd') precision.
    """
    from phcpy.solver import solve
    pols = ['(1-s)*y + s*(y-1);', \
            'x^2 + y^2 + z^2 - 4;' , \
            '(x-1)^2 + y^2 - 1;', \
            's;']
    sols = solve(pols, silent=True, precision=prc)
    print "The solutions on the Viviani curve :"
    for sol in sols:
        print sol
    return (pols[:3], sols)

def viviani2(precision='d'):
    """
    Computes the power series expansion for the Viviani curve,
    from a natural paramter perspective.
    The default precision is double ('d').  Other precisions
    are double double ('dd') and quad double ('qd').
    """
    pols = [ '2*t^2 - x;', \
             'x^2 + y^2 + z^2 - 4;' , \
             '(x-1)^2 + y^2 - 1;']
    lser = [ '2*t^2;', '2*t;', '2;']
    if precision == 'd':
        nser = standard_newton_power_series(pols, lser, nbr=8)
    elif precision == 'dd':
        nser = dobldobl_newton_power_series(pols, lser, nbr=8)
    elif precision == 'qd':
        nser = quaddobl_newton_power_series(pols, lser, nbr=8)
    else:
        print 'invalid argument for the precision'
    print nser

def apollonius(precision='d'):
    """
    Test on computing the power series at a double solution
    for the problem of Apolonius.
    The parameter t is the fourth variable, whence we call
    Newton's method with idx equal to four.
    """
    pols = [ 'x1^2 + 3*x2^2 - r^2 - 2*r - 1;', \
             'x1^2 + 3*x2^2 - r^2 - 4*x1 - 2*r + 3;', \
       '3*t^2 + x1^2 - 6*t*x2 + 3*x2^2 - r^2 + 6*t - 2*x1 - 6*x2 + 2*r + 3;']
 
    lser1 = ['1;', '1 + 0.536*t;', '1 + 0.904*t;']
    lser2 = ['1;', '1 + 7.464*t;', '1 + 11.196*t;']
    if precision == 'd':
        nser1 = standard_newton_power_series(pols, lser1, idx=4, nbr=8)
        nser2 = standard_newton_power_series(pols, lser2, idx=4, nbr=8)
    elif precision == 'dd':
        nser1 = dobldobl_newton_power_series(pols, lser1, idx=4, nbr=8)
        nser2 = dobldobl_newton_power_series(pols, lser2, idx=4, nbr=8)
    elif precision == 'qd':
        nser1 = quaddobl_newton_power_series(pols, lser1, idx=4, nbr=8)
        nser2 = quaddobl_newton_power_series(pols, lser2, idx=4, nbr=8)
    else:
        print 'invalid argument for the precision'
    print nser1
    print nser2

def test(precision='d'):
    """
    Tests the application of Newton's method to compute power
    series solutions of a polynomial system.
    """
    (pols, sols) = viviani(precision)
    if precision == 'd':
        sersols = standard_newton_series(pols, sols)
    elif precision == 'dd':
        sersols = dobldobl_newton_series(pols, sols)
    elif precision == 'qd':
        sersols = quaddobl_newton_series(pols, sols)
    for series in sersols:
        print series

if __name__ == "__main__":
    #test('d')
    #test('dd')
    #test('qd')
    #viviani2('qd')
    apollonius()
