"""
Exports functions to compute series of solution curves defined by
polynomial homotopies using Newton's method.
"""

from ctypes import c_int, c_double, pointer
from version import get_phcfun, int4a2nbr, int4a2str, nbr2int4a
from polynomials import string_of_symbols, number_of_symbols
from polynomials import set_double_system, get_double_system
from polynomials import set_double_double_system, get_double_double_system
from polynomials import set_quad_double_system, get_quad_double_system
from polynomials import copy_from_double_syspool, size_double_syspool
from polynomials import copy_from_double_double_syspool
from polynomials import size_double_double_syspool
from polynomials import copy_from_quad_double_syspool
from polynomials import size_quad_double_syspool
from solutions import set_double_solutions, get_double_solutions
from solutions import set_double_double_solutions, get_double_double_solutions
from solutions import set_quad_double_solutions, get_quad_double_solutions
from solver import solve_double_system, write_double_solutions
from solver import solve_double_double_system, write_double_double_solutions
from solver import solve_quad_double_system, write_quad_double_solutions

def replace_symbol(pol, idx):
    """
    In the polynomial pol, 
    replaces the first symbol by the symbol at place idx.
    """
    var = string_of_symbols()
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

def double_newton_series(pols, sols, idx=1, maxdeg=4, nbr=4, vrblvl=0):
    r"""
    Computes series in double precision for the polynomials
    in *pols*, where the leading coefficients are the solutions in *sols*.
    On entry are the following five parameters:

    *pols*: a list of string representations of polynomials,

    *sols*: a list of solutions of the polynomials in pols,

    *idx*: index of the series parameter, by default equals 1,

    *maxdeg*: maximal degree of the series,

    *nbr*: number of steps with Newton's method,

    *vrblvl*: the verbose level.

    On return is a list of lists of strings.  Each lists of strings
    represents the series solution for the variables in the list *pols*.
    """
    nbsym = number_of_symbols(pols, vrblvl)
    if vrblvl > 0:
        print('-> double_newton_series, idx :', idx)
        print('-> double_newton_series, the polynomials :')
        for pol in pols:
            print(pol)
        print('number of variables :', nbsym)
    set_double_system(nbsym, pols, vrblvl)
    set_double_solutions(nbsym, sols)
    syspol = get_double_system(vrblvl)
    for pol in syspol:
        print(pol)
    write_double_solutions(vrblvl)
    phc = get_phcfun()
    apars = int4a2nbr([idx, maxdeg, nbr], (vrblvl > 0))
    bbb = pointer(c_int(vrblvl))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    if vrblvl > 0:
        print('-> double_newton_series calls phc ...')
        print('apars =', nbr2int4a(apars));
    retval = phc(691, apars, bbb, ccc, vrb)
    fail = (retval > 0)
    size = (-1 if fail else size_double_syspool(vrblvl))
    if vrblvl > 0:
        print('the return value of double_newton_series :', retval)
        if size == -1:
            print("An error occurred in the execution of Newton's method.")
        else:
            print('Computed %d series solutions.' % size)
    result = []
    for k in range(1, size+1):
        if vrblvl > 0:
            print('')
            print('retrieving series', k, '...')
            print('')
        copy_from_double_syspool(k)
        sersol = get_double_system(vrblvl)
        result.append(substitute_symbol(sersol, idx))
    return result

def double_double_newton_series(pols, sols, idx=1, maxdeg=4, nbr=4, vrblvl=0):
    r"""
    Computes series in double double precision for the polynomials
    in *pols*, where the leading coefficients are the solutions in *sols*.
    On entry are the following five parameters:

    *pols*: a list of string representations of polynomials,

    *sols*: a list of solutions of the polynomials in pols,

    *idx*: index of the series parameter, by default equals 1,

    *maxdeg*: maximal degree of the series,

    *nbr*: number of steps with Newton's method,

    *vrblvl*: the verbose level.

    On return is a list of lists of strings.  Each lists of strings
    represents the series solution for the variables in the list *pols*.
    """
    nbsym = number_of_symbols(pols, vrblvl)
    if vrblvl > 0:
        print('-> double_double_newton_series, idx :', idx)
        print('-> double_double_newton_series, the polynomials :')
        for pol in pols:
            print(pol)
        print('number of variables :', nbsym)
    set_double_double_system(nbsym, pols, vrblvl)
    set_double_double_solutions(nbsym, sols)
    syspol = get_double_double_system(vrblvl)
    for pol in syspol:
        print(pol)
    write_double_double_solutions(vrblvl)
    phc = get_phcfun()
    apars = int4a2nbr([idx, maxdeg, nbr], (vrblvl > 0))
    bbb = pointer(c_int(vrblvl))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    if vrblvl > 0:
        print('-> double_double_newton_series calls phc ...')
        print('apars =', nbr2int4a(apars));
    retval = phc(692, apars, bbb, ccc, vrb)
    fail = (retval > 0)
    size = (-1 if fail else size_double_double_syspool(vrblvl))
    if vrblvl > 0:
        print('the return value of double_double_newton_series :', retval)
        if size == -1:
            print("An error occurred in the execution of Newton's method.")
        else:
            print('Computed %d series solutions.' % size)
    result = []
    for k in range(1, size+1):
        if vrblvl > 0:
            print('')
            print('retrieving series', k, '...')
            print('')
        copy_from_double_double_syspool(k)
        sersol = get_double_double_system(vrblvl)
        result.append(substitute_symbol(sersol, idx))
    return result

def quad_double_newton_series(pols, sols, idx=1, maxdeg=4, nbr=4, vrblvl=0):
    r"""
    Computes series in quad double precision for the polynomials
    in *pols*, where the leading coefficients are the solutions in *sols*.
    On entry are the following five parameters:

    *pols*: a list of string representations of polynomials,

    *sols*: a list of solutions of the polynomials in pols,

    *idx*: index of the series parameter, by default equals 1,

    *maxdeg*: maximal degree of the series,

    *nbr*: number of steps with Newton's method,

    *vrblvl*: the verbose level.

    On return is a list of lists of strings.  Each lists of strings
    represents the series solution for the variables in the list *pols*.
    """
    nbsym = number_of_symbols(pols, vrblvl)
    if vrblvl > 0:
        print('-> quad_double_newton_series, idx :', idx)
        print('-> quad_double_newton_series, the polynomials :')
        for pol in pols:
            print(pol)
        print('number of variables :', nbsym)
    set_quad_double_system(nbsym, pols, vrblvl)
    set_quad_double_solutions(nbsym, sols)
    syspol = get_quad_double_system(vrblvl)
    for pol in syspol:
        print(pol)
    write_quad_double_solutions(vrblvl)
    phc = get_phcfun()
    apars = int4a2nbr([idx, maxdeg, nbr], (vrblvl > 0))
    bbb = pointer(c_int(vrblvl))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    if vrblvl > 0:
        print('-> quad_double_newton_series calls phc ...')
        print('apars =', nbr2int4a(apars));
    retval = phc(693, apars, bbb, ccc, vrb)
    fail = (retval > 0)
    size = (-1 if fail else size_quad_double_syspool(vrblvl))
    if vrblvl > 0:
        print('the return value of quad_double_newton_series :', retval)
        if size == -1:
            print("An error occurred in the execution of Newton's method.")
        else:
            print('Computed %d series solutions.' % size)
    result = []
    for k in range(1, size+1):
        if vrblvl > 0:
            print('')
            print('retrieving series', k, '...')
            print('')
        copy_from_quad_double_syspool(k)
        sersol = get_quad_double_system(vrblvl)
        result.append(substitute_symbol(sersol, idx))
    return result

def test_double_viviani(vrblvl=0):
    """
    Returns the system which stores the Viviani curve,
    with some solutions intersected with a plane,
    in double precision.
    """
    pols = ['(1-s)*y + s*(y-1);', \
            'x^2 + y^2 + z^2 - 4;' , \
            '(x-1)^2 + y^2 - 1;', \
            's;']
    set_double_system(len(pols), pols, vrblvl)
    nbr, roco = solve_double_system(vrblvl=vrblvl)
    sols = get_double_solutions(vrblvl)
    print("The solutions on the Viviani curve :")
    for sol in sols:
        print(sol)
    sersols = double_newton_series(pols[:3], sols, vrblvl=vrblvl)
    for srs in sersols:
        print(srs)

def test_double_double_viviani(vrblvl=0):
    """
    Returns the system which stores the Viviani curve,
    with some solutions intersected with a plane,
    in double double precision.
    """
    pols = ['(1-s)*y + s*(y-1);', \
            'x^2 + y^2 + z^2 - 4;' , \
            '(x-1)^2 + y^2 - 1;', \
            's;']
    set_double_double_system(len(pols), pols, vrblvl)
    nbr, roco = solve_double_double_system(vrblvl=vrblvl)
    sols = get_double_double_solutions(vrblvl)
    print("The solutions on the Viviani curve :")
    for sol in sols:
        print(sol)
    sersols = double_double_newton_series(pols[:3], sols, vrblvl=vrblvl)
    for srs in sersols:
        print(srs)

def test_quad_double_viviani(vrblvl=0):
    """
    Returns the system which stores the Viviani curve,
    with some solutions intersected with a plane,
    in quad double precision.
    """
    pols = ['(1-s)*y + s*(y-1);', \
            'x^2 + y^2 + z^2 - 4;' , \
            '(x-1)^2 + y^2 - 1;', \
            's;']
    set_quad_double_system(len(pols), pols, vrblvl)
    nbr, roco = solve_quad_double_system(vrblvl=vrblvl)
    sols = get_quad_double_solutions(vrblvl)
    print("The solutions on the Viviani curve :")
    for sol in sols:
        print(sol)
    sersols = quad_double_newton_series(pols[:3], sols, vrblvl=vrblvl)
    for srs in sersols:
        print(srs)

if __name__=="__main__":
    # test_double_viviani(10)
    # test_double_double_viviani(10)
    test_quad_double_viviani(10)
