"""
A homotopy is a family of polynomial systems which connects a given target
system to a start system.  This module exports several start systems.
"""
from ctypes import c_int, c_double, pointer
from phcpy.version import get_phcfun
from phcpy.polynomials import set_double_system, string_of_symbols
from phcpy.polynomials import degree_of_double_polynomial
from phcpy.examples import noon3
from phcpy.solver import solve_checkin, solve

def total_degree(pols, vrblvl=0):
    """
    Given in pols a list of string representations of polynomials,
    returns the product of the degrees of the polynomials,
    the so-called total degree which bounds the number of
    isolated solutions of the polynomial system.
    The system is assumed to be square.
    The value of the verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in total degree, pols :')
        for pol in pols:
            print(pol)
    set_double_system(len(pols), pols, vrblvl)
    phc = get_phcfun()
    deg = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    if vrblvl > 0:
        print('-> total_degree calls phc', end='')
    retval = phc(28, deg, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the total degree :', deg[0])
    return deg[0]

def total_degree_start_system(pols, checkin=True, vrblvl=0):
    r"""
    Returns the system and solutions of the total degree start system
    for the polynomials represented by the strings in the list *pols*.
    If *checkin*, then the list *pols* is tested to see if *pols* defines
    a square polynomial system.  If the input system is not square,
    then an error message is printed and None is returned.
    """
    if vrblvl > 0:
        print('in total degree, pols :')
        for pol in pols:
            print(pol)
    if checkin:
        errmsg = 'Start systems are defined only for square systems,'
        if not solve_checkin(pols, errmsg):
            return None
    dim = len(pols)
    set_double_system(dim, pols, vrblvl)
    svars = string_of_symbols(200, vrblvl)
    degrees = [degree_of_double_polynomial(k+1) for k in range(dim)]
    result = []
    for ind in range(dim):
        result.append(svars[ind]+'^'+str(degrees[ind])+' - 1;')
    return (result, solve(result))

def main():
    """
    Runs some tests.
    """
    lvl = 10
    pols = noon3()
    totdeg = total_degree(pols, lvl)
    print('the total degree of noon3 :', totdeg)
    (start, startsols) = total_degree_start_system(pols, vrblvl=lvl)
    print('the start system :')
    for pol in start:
        print(pol)
    print('the start solutions :')
    for (idx, sol) in enumerate(startsols):
        print('Solution', idx+1, ':')
        print(sol)

if __name__=='__main__':
    main()
