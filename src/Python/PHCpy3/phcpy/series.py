"""
The module series exports functions to compute power series solutions with
Newton's method in double, double double, or quad double precision.
"""

def replace_symbol(pol, idx):
    """
    In the polynomial pol, 
    replaces the first symbol by the symbol at place idx.
    """
    from phcpy.phcpy2c3 import py2c_syscon_string_of_symbols
    sbl = py2c_syscon_string_of_symbols()
    var = sbl.split(' ')
    result = pol.replace(var[0], var[idx-1])
    return result

def substitute_symbol(pols, idx):
    """
    Given in pols is a list of polynomials,
    replaces the first symbol by the symbol at place idx.
    """
    if idx == 1:
        return pols
    else:
        result = []
        for pol in pols:
            result.append(replace_symbol(pol, idx))
        return result

def standard_newton_series(pols, sols, idx=1, maxdeg=4, nbr=4, verbose=True):
    r"""
    Computes series in standard double precision for the polynomials
    in *pols*, where the leading coefficients are the solutions in *sols*.
    On entry are the following five parameters:

    *pols*: a list of string representations of polynomials,

    *sols*: a list of solutions of the polynomials in pols,

    *idx*: index of the series parameter, by default equals 1,

    *maxdeg*: maximal degree of the series,

    *nbr*: number of steps with Newton's method,

    *verbose*: whether to write intermediate output to screen or not.

    On return is a list of lists of strings.  Each lists of strings
    represents the series solution for the variables in the list *pols*.
    """
    from phcpy.solver import number_of_symbols
    from phcpy.interface import store_standard_solutions
    from phcpy.interface import store_standard_system, load_standard_system
    from phcpy.phcpy2c3 import py2c_standard_Newton_series as newton
    from phcpy.phcpy2c3 import py2c_syspool_standard_size as poolsize
    from phcpy.phcpy2c3 import py2c_syspool_copy_to_standard_container
    nbsym = number_of_symbols(pols)
    if verbose:
        print("the polynomials :")
        for pol in pols:
            print(pol)
        print("Number of variables :", nbsym)
    store_standard_system(pols, nbvar=nbsym)
    store_standard_solutions(nbsym, sols)
    fail = newton(idx, maxdeg, nbr, int(verbose))
    size = (-1 if fail else poolsize())
    if verbose:
        if size == -1:
            print("An error occurred in the execution of Newton's method.")
        else:
            print("Computed %d series solutions." % size)
    result = []
    for k in range(1, size+1):
        py2c_syspool_copy_to_standard_container(k)
        sersol = load_standard_system()
        result.append(substitute_symbol(sersol, idx))
    return result

def dobldobl_newton_series(pols, sols, idx=1, maxdeg=4, nbr=4, verbose=True):
    r"""
    Computes series in double double precision for the polynomials
    in *pols*, where the leading coefficients are the solutions in *sols*.
    On entry are the following five parameters:

    *pols*: a list of string representations of polynomials,

    *sols*: a list of solutions of the polynomials in pols,

    *idx*: index of the series parameter, by default equals 1,

    *maxdeg*: maximal degree of the series,

    *nbr*: number of steps with Newton's method,

    *verbose*: whether to write intermediate output to screen or not.

    On return is a list of lists of strings.  Each lists of strings
    represents the series solution for the variables in the list *pols*.
    """
    from phcpy.solver import number_of_symbols
    from phcpy.interface import store_dobldobl_solutions
    from phcpy.interface import store_dobldobl_system, load_dobldobl_system
    from phcpy.phcpy2c3 import py2c_dobldobl_Newton_series as newton
    from phcpy.phcpy2c3 import py2c_syspool_dobldobl_size as poolsize
    from phcpy.phcpy2c3 import py2c_syspool_copy_to_dobldobl_container
    nbsym = number_of_symbols(pols)
    if verbose:
        print("the polynomials :")
        for pol in pols:
            print(pol)
        print("Number of variables :", nbsym)
    store_dobldobl_system(pols, nbvar=nbsym)
    store_dobldobl_solutions(nbsym, sols)
    fail = newton(idx, maxdeg, nbr, int(verbose))
    size = (-1 if fail else poolsize())
    if verbose:
        if size == -1:
            print("An error occurred in the execution of Newton's method.")
        else:
            print("Computed %d series solutions." % size)
    result = []
    for k in range(1, size+1):
        py2c_syspool_copy_to_dobldobl_container(k)
        sersol = load_dobldobl_system()
        result.append(substitute_symbol(sersol, idx))
    return result

def quaddobl_newton_series(pols, sols, idx=1, maxdeg=4, nbr=4, verbose=True):
    """
    Computes series in quad double precision for the polynomials
    in *pols*, where the leading coefficients are the solutions in *sols*.
    On entry are the following five parameters:

    *pols*: a list of string representations of polynomials,

    *sols*: a list of solutions of the polynomials in pols,

    *idx*: index of the series parameter, by default equals 1,

    *maxdeg*: maximal degree of the series,

    *nbr*: number of steps with Newton's method,

    *verbose*: whether to write intermediate output to screen or not.

    On return is a list of lists of strings.  Each lists of strings
    represents the series solution for the variables in the list *pols*.
    """
    from phcpy.solver import number_of_symbols
    from phcpy.interface import store_quaddobl_solutions
    from phcpy.interface import store_quaddobl_system, load_quaddobl_system
    from phcpy.phcpy2c3 import py2c_quaddobl_Newton_series as newton
    from phcpy.phcpy2c3 import py2c_syspool_quaddobl_size as poolsize
    from phcpy.phcpy2c3 import py2c_syspool_copy_to_quaddobl_container
    nbsym = number_of_symbols(pols)
    if verbose:
        print("the polynomials :")
        for pol in pols:
            print(pol)
        print("Number of variables :", nbsym)
    store_quaddobl_system(pols, nbvar=nbsym)
    store_quaddobl_solutions(nbsym, sols)
    fail = newton(idx, maxdeg, nbr, int(verbose))
    size = (-1 if fail else poolsize())
    if verbose:
        if size == -1:
            print("An error occurred in the execution of Newton's method.")
        else:
            print("Computed %d series solutions." % size)
    result = []
    for k in range(1, size+1):
        py2c_syspool_copy_to_quaddobl_container(k)
        sersol = load_quaddobl_system()
        result.append(substitute_symbol(sersol, idx))
    return result

def checkin_newton_power_series(nbsym, lser, idx):
    """
    Given in nbsym the number of symbols in the polynomial system,
    in lser the list of leading terms in the series and 
    in idx the index of the parameter, returns True
    if nbsym = len(lser) if idx == 0, or otherwise
    if nbsym = len(lser) + 1 if idx != 0.
    An error message is written and False is returned
    if the above conditions are not satisfied.
    """
    if idx == 0:
        okay = (nbsym == len(lser))
    else:
        okay = (nbsym == len(lser) + 1)
    if not okay:
        if idx == 0:
            dim = nbsym
        else:
            dim = nbsym - 1
        print('Wrong length of list of leading terms, should be', \
            str(dim) + '.')
       
    return okay   

def standard_newton_power_series(pols, lser, idx=1, maxdeg=4, nbr=4, \
    checkin=True, verbose=True):
    r"""
    Computes series in standard double precision for the polynomials
    in *pols*, where the leading terms are given in the list *lser*.
    On entry are the following five parameters:

    *pols*: a list of string representations of polynomials,

    *lser*: a list of polynomials in the series parameter (e.g.: t),
    for use as start terms in Newton's method,

    *idx*: index of the series parameter, by default equals 1,

    *maxdeg*: maximal degree of the series,

    *nbr*: number of steps with Newton's method,

    *checkin*: checks whether the number of symbols in pols matches
    the length of the list lser if idx == 0, or is one less than 
    the length of the list lser if idx != 0.  If the conditions are
    not satisfied, then an error message is printed and lser is returned.

    *verbose*: whether to write intermediate output to screen or not.

    On return is a list of lists of strings.  Each lists of strings
    represents the series solution for the variables in the list *pols*.
    """
    from phcpy.solver import number_of_symbols
    from phcpy.interface import store_standard_system, load_standard_system
    from phcpy.phcpy2c3 import py2c_standard_Newton_power_series as newton
    from phcpy.phcpy2c3 import py2c_syspool_standard_init
    from phcpy.phcpy2c3 import py2c_syspool_standard_create
    from phcpy.phcpy2c3 import py2c_syspool_standard_size as poolsize
    from phcpy.phcpy2c3 import py2c_syspool_copy_to_standard_container
    from phcpy.phcpy2c3 import py2c_syspool_standard_clear
    nbsym = number_of_symbols(pols)
    if verbose:
        print("the polynomials :")
        for pol in pols:
            print(pol)
        print("Number of variables :", nbsym)
    if checkin:
        if not checkin_newton_power_series(nbsym, lser, idx):
            return lser
    store_standard_system(lser, nbvar=1)
    py2c_syspool_standard_init(1);
    py2c_syspool_standard_create(1);
    store_standard_system(pols, nbvar=nbsym)
    fail = newton(idx, maxdeg, nbr, int(verbose))
    size = (-1 if fail else poolsize())
    if verbose:
        if size == -1:
            print("An error occurred in the execution of Newton's method.")
        else:
            print("Computed one series solution.")
    py2c_syspool_copy_to_standard_container(1)
    result = load_standard_system()
    result = substitute_symbol(result, idx)
    py2c_syspool_standard_clear()
    return result

def dobldobl_newton_power_series(pols, lser, idx=1, maxdeg=4, nbr=4, \
    checkin=True, verbose=True):
    r"""
    Computes series in double double precision for the polynomials
    in *pols*, where the leading terms are given in the list *lser*.
    On entry are the following five parameters:

    *pols*: a list of string representations of polynomials,

    *lser*: a list of polynomials in the series parameter (e.g.: t),
    for use as start terms in Newton's method,

    *idx*: index of the series parameter, by default equals 1,

    *maxdeg*: maximal degree of the series,

    *nbr*: number of steps with Newton's method,

    *checkin*: checks whether the number of symbols in pols matches
    the length of the list lser if idx == 0, or is one less than 
    the length of the list lser if idx != 0.  If the conditions are
    not satisfied, then an error message is printed and lser is returned.

    *verbose*: whether to write intermediate output to screen or not.

    On return is a list of lists of strings.  Each lists of strings
    represents the series solution for the variables in the list *pols*.
    """
    from phcpy.solver import number_of_symbols
    from phcpy.interface import store_dobldobl_system, load_dobldobl_system
    from phcpy.phcpy2c3 import py2c_dobldobl_Newton_power_series as newton
    from phcpy.phcpy2c3 import py2c_syspool_dobldobl_init
    from phcpy.phcpy2c3 import py2c_syspool_dobldobl_create
    from phcpy.phcpy2c3 import py2c_syspool_dobldobl_size as poolsize
    from phcpy.phcpy2c3 import py2c_syspool_copy_to_dobldobl_container
    from phcpy.phcpy2c3 import py2c_syspool_dobldobl_clear
    nbsym = number_of_symbols(pols)
    if verbose:
        print("the polynomials :")
        for pol in pols:
            print(pol)
        print("Number of variables :", nbsym)
    if checkin:
        if not checkin_newton_power_series(nbsym, lser, idx):
            return lser
    store_dobldobl_system(lser, nbvar=1)
    py2c_syspool_dobldobl_init(1);
    py2c_syspool_dobldobl_create(1);
    store_dobldobl_system(pols, nbvar=nbsym)
    fail = newton(idx, maxdeg, nbr, int(verbose))
    size = (-1 if fail else poolsize())
    if verbose:
        if size == -1:
            print("An error occurred in the execution of Newton's method.")
        else:
            print("Computed one series solution.")
    py2c_syspool_copy_to_dobldobl_container(1)
    result = load_dobldobl_system()
    result = substitute_symbol(result, idx)
    py2c_syspool_dobldobl_clear()
    return result

def quaddobl_newton_power_series(pols, lser, idx=1, maxdeg=4, nbr=4, \
    checkin=True, verbose=True):
    r"""
    Computes series in quad double precision for the polynomials
    in *pols*, where the leading terms are given in the list *lser*.
    On entry are the following five parameters:

    *pols*: a list of string representations of polynomials,

    *lser*: a list of polynomials in the series parameter (e.g.: t),
    for use as start terms in Newton's method,

    *idx*: index of the series parameter, by default equals 1,

    *maxdeg*: maximal degree of the series,

    *nbr*: number of steps with Newton's method,

    *checkin*: checks whether the number of symbols in pols matches
    the length of the list lser if idx == 0, or is one less than 
    the length of the list lser if idx != 0.  If the conditions are
    not satisfied, then an error message is printed and lser is returned.

    *verbose*: whether to write intermediate output to screen or not.

    On return is a list of lists of strings.  Each lists of strings
    represents the series solution for the variables in the list *pols*.
    """
    from phcpy.solver import number_of_symbols
    from phcpy.interface import store_quaddobl_system, load_quaddobl_system
    from phcpy.phcpy2c3 import py2c_quaddobl_Newton_power_series as newton
    from phcpy.phcpy2c3 import py2c_syspool_quaddobl_init
    from phcpy.phcpy2c3 import py2c_syspool_quaddobl_create
    from phcpy.phcpy2c3 import py2c_syspool_quaddobl_size as poolsize
    from phcpy.phcpy2c3 import py2c_syspool_copy_to_quaddobl_container
    from phcpy.phcpy2c3 import py2c_syspool_quaddobl_clear
    nbsym = number_of_symbols(pols)
    if verbose:
        print("the polynomials :")
        for pol in pols:
            print(pol)
        print("Number of variables :", nbsym)
    if checkin:
        if not checkin_newton_power_series(nbsym, lser, idx):
            return lser
    store_quaddobl_system(lser, nbvar=1)
    py2c_syspool_quaddobl_init(1);
    py2c_syspool_quaddobl_create(1);
    store_quaddobl_system(pols, nbvar=nbsym)
    fail = newton(idx, maxdeg, nbr, int(verbose))
    size = (-1 if fail else poolsize())
    if verbose:
        if size == -1:
            print("An error occurred in the execution of Newton's method.")
        else:
            print("Computed one series solution.")
    py2c_syspool_copy_to_quaddobl_container(1)
    result = load_quaddobl_system()
    result = substitute_symbol(result, idx)
    py2c_syspool_quaddobl_clear()
    return result

def make_fractions(pols):
    """
    Given a list of string representations for the numerator and
    denominator polynomials in its even and odd numbered indices,
    returns a list of string representations for the fractions.
    """
    result = []
    (nbr, idx) = (len(pols)//2, 0)
    for k in range(nbr):
        (num, den) = (pols[idx], pols[idx+1])
        idx = idx + 2
        frac = '(' + num[:-1] + ')/(' + den[:-1] + ')'
        result.append(frac)
    return result

def rational_forms(pols):
    """
    Given a list of lists of string representations for the numerators
    and denominators, returns the proper rational representations for
    the Pade approximants.
    """
    result = []
    for pol in pols:
        result.append(make_fractions(pol))
    return result

def standard_pade_approximants(pols, sols, idx=1, numdeg=2, dendeg=2, \
    nbr=4, verbose=True):
    r"""
    Computes Pade approximants based on the series in standard double 
    precision for the polynomials in *pols*, where the leading 
    coefficients of the series are the solutions in *sols*.
    On entry are the following seven parameters:

    *pols*: a list of string representations of polynomials,

    *sols*: a list of solutions of the polynomials in *pols*,

    *idx*: index of the series parameter, by default equals 1,

    *numdeg*: the degree of the numerator,

    *dendeg*: the degree of the denominator,

    *nbr*: number of steps with Newton's method,

    *verbose*: whether to write intermediate output to screen or not.

    On return is a list of lists of strings.  Each lists of strings
    represents the series solution for the variables in the list *pols*.
    """
    from phcpy.solver import number_of_symbols
    from phcpy.interface import store_standard_solutions
    from phcpy.interface import store_standard_system, load_standard_system
    from phcpy.phcpy2c3 \
        import py2c_standard_Pade_approximant as Pade_approximants
    from phcpy.phcpy2c3 import py2c_syspool_standard_size as poolsize
    from phcpy.phcpy2c3 import py2c_syspool_copy_to_standard_container
    from phcpy.phcpy2c3 import py2c_syspool_standard_clear
    nbsym = number_of_symbols(pols)
    if verbose:
        print("the polynomials :")
        for pol in pols:
            print(pol)
        print("Number of variables :", nbsym)
    store_standard_system(pols, nbvar=nbsym)
    store_standard_solutions(nbsym, sols)
    fail = Pade_approximants(idx, numdeg, dendeg, nbr, int(verbose))
    size = (-1 if fail else poolsize())
    if verbose:
        if size == -1:
            print("An error occurred in the Pade constructor.")
        else:
            print("Computed %d Pade approximants." % size)
    result = []
    for k in range(1, size+1):
        py2c_syspool_copy_to_standard_container(k)
        sersol = load_standard_system()
        substsersol = substitute_symbol(sersol, idx)
        result.append(make_fractions(substsersol))
    py2c_syspool_standard_clear()
    return result

def dobldobl_pade_approximants(pols, sols, idx=1, numdeg=2, dendeg=2, \
    nbr=4, verbose=True):
    r"""
    Computes Pade approximants based on the series in double double 
    precision for the polynomials in *pols*, where the leading 
    coefficients of the series are the solutions in *sols*.
    On entry are the following seven parameters:

    *pols*: a list of string representations of polynomials,

    *sols*: a list of solutions of the polynomials in *pols*,

    *idx*: index of the series parameter, by default equals 1,

    *numdeg*: the degree of the numerator,

    *dendeg*: the degree of the denominator,

    *nbr*: number of steps with Newton's method,

    *verbose*: whether to write intermediate output to screen or not.

    On return is a list of lists of strings.  Each lists of strings
    represents the series solution for the variables in the list *pols*.
    """
    from phcpy.solver import number_of_symbols
    from phcpy.interface import store_dobldobl_solutions
    from phcpy.interface import store_dobldobl_system, load_dobldobl_system
    from phcpy.phcpy2c3 \
        import py2c_dobldobl_Pade_approximant as Pade_approximants
    from phcpy.phcpy2c3 import py2c_syspool_dobldobl_size as poolsize
    from phcpy.phcpy2c3 import py2c_syspool_copy_to_dobldobl_container
    from phcpy.phcpy2c3 import py2c_syspool_dobldobl_clear
    nbsym = number_of_symbols(pols)
    if verbose:
        print("the polynomials :")
        for pol in pols:
            print(pol)
        print("Number of variables :", nbsym)
    store_dobldobl_system(pols, nbvar=nbsym)
    store_dobldobl_solutions(nbsym, sols)
    fail = Pade_approximants(idx, numdeg, dendeg, nbr, int(verbose))
    size = (-1 if fail else poolsize())
    if verbose:
        if size == -1:
            print("An error occurred in the Pade constructor.")
        else:
            print("Computed %d Pade approximants." % size)
    result = []
    for k in range(1, size+1):
        py2c_syspool_copy_to_dobldobl_container(k)
        sersol = load_dobldobl_system()
        substsersol = substitute_symbol(sersol, idx)
        result.append(make_fractions(substsersol))
    py2c_syspool_dobldobl_clear()
    return result

def quaddobl_pade_approximants(pols, sols, idx=1, numdeg=2, dendeg=2, \
    nbr=4, verbose=True):
    r"""
    Computes Pade approximants based on the series in quad double 
    precision for the polynomials in *pols*, where the leading 
    coefficients of the series are the solutions in *sols*.
    On entry are the following seven parameters:

    *pols*: a list of string representations of polynomials,

    *sols*: a list of solutions of the polynomials in *pols*,

    *idx*: index of the series parameter, by default equals 1,

    *numdeg*: the degree of the numerator,

    *dendeg*: the degree of the denominator,

    *nbr*: number of steps with Newton's method,

    *verbose*: whether to write intermediate output to screen or not.

    On return is a list of lists of strings.  Each lists of strings
    represents the series solution for the variables in the list *pols*.
    """
    from phcpy.solver import number_of_symbols
    from phcpy.interface import store_quaddobl_solutions
    from phcpy.interface import store_quaddobl_system, load_quaddobl_system
    from phcpy.phcpy2c3 \
        import py2c_quaddobl_Pade_approximant as Pade_approximants
    from phcpy.phcpy2c3 import py2c_syspool_quaddobl_size as poolsize
    from phcpy.phcpy2c3 import py2c_syspool_copy_to_quaddobl_container
    from phcpy.phcpy2c3 import py2c_syspool_quaddobl_clear
    nbsym = number_of_symbols(pols)
    if verbose:
        print("the polynomials :")
        for pol in pols:
            print(pol)
        print("Number of variables :", nbsym)
    store_quaddobl_system(pols, nbvar=nbsym)
    store_quaddobl_solutions(nbsym, sols)
    fail = Pade_approximants(idx, numdeg, dendeg, nbr, int(verbose))
    size = (-1 if fail else poolsize())
    if verbose:
        if size == -1:
            print("An error occurred in the Pade constructor.")
        else:
            print("Computed %d Pade approximants." % size)
    result = []
    for k in range(1, size+1):
        py2c_syspool_copy_to_quaddobl_container(k)
        sersol = load_quaddobl_system()
        substsersol = substitute_symbol(sersol, idx)
        result.append(make_fractions(substsersol))
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
    sols = solve(pols, verbose=False, precision=prc)
    print("The solutions on the Viviani curve :")
    for sol in sols:
        print(sol)
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
        nser = standard_newton_power_series(pols, lser, maxdeg=12, nbr=8)
    elif precision == 'dd':
        nser = dobldobl_newton_power_series(pols, lser, maxdeg=12, nbr=8)
    elif precision == 'qd':
        nser = quaddobl_newton_power_series(pols, lser, maxdeg=12, nbr=8)
    else:
        print('invalid argument for the precision')
    print(nser)

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
        nser1 = standard_newton_power_series(pols, lser1, idx=4, nbr=7)
        nser2 = standard_newton_power_series(pols, lser2, idx=4, nbr=7)
    elif precision == 'dd':
        nser1 = dobldobl_newton_power_series(pols, lser1, idx=4, nbr=7)
        nser2 = dobldobl_newton_power_series(pols, lser2, idx=4, nbr=7)
    elif precision == 'qd':
        nser1 = quaddobl_newton_power_series(pols, lser1, idx=4, nbr=7)
        nser2 = quaddobl_newton_power_series(pols, lser2, idx=4, nbr=7)
    else:
        print('invalid argument for the precision')
    print(nser1)
    print(nser2)

def example4pade(prc='d'):
    """
    The function f(z) = ((1 + 1/2*z)/(1 + 2*z))^(1/2) is
    a solution x(s) of (1-s)*(x^2 - 1) + s*(3*x^2 - 3/2) = 0
    """
    pols = ['(x^2 - 1)*(1-s) + (3*x^2 - 3/2)*s;', 's;']
    from phcpy.solver import solve
    sols = solve(pols, verbose=False, precision=prc)
    for sol in sols:
        print(sol)
    if prc == 'd':
        sers = standard_newton_series(pols[:1], sols, idx=2)
    else:
        sers = dobldobl_newton_series(pols[:1], sols, idx=2)
    print('the series solutions :')
    for ser in sers:
        print(ser)
    if prc == 'd':
        pade = standard_pade_approximants(pols[:1], sols, idx=2)
    elif prc == 'dd':
        pade = dobldobl_pade_approximants(pols[:1], sols, idx=2)
    elif prc == 'qd':
        pade = quaddobl_pade_approximants(pols[:1], sols, idx=2)
    else:
        print('wrong value for the precision')
    print('the Pade approximants :')
    for pad in pade:
        print(pad)

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
        print(series)

if __name__ == "__main__":
    #test('d')
    #test('dd')
    #test('qd')
    #viviani2('d')
    #viviani2('dd')
    #viviani2('qd')
    #apollonius()
    example4pade('qd')
