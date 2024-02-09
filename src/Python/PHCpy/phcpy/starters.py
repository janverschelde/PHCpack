"""
A an artificial parameter homotopy is a family of polynomial systems 
which connects a given target system to a start system.
The module starters exports several functions to construct start systems.
Analoguous to the culinary starters, start systems may be viewed as
systems that are easier to solve (or digest) than general systems.
"""
from ctypes import c_int32, c_double, pointer, create_string_buffer
from phcpy.version import int4a2str, str2int4a
from phcpy.version import get_phcfun
from phcpy.polynomials import set_double_system, string_of_symbols
from phcpy.polynomials import get_double_system
from phcpy.polynomials import degree_of_double_polynomial
from phcpy.solutions import get_double_solutions
from phcpy.examples import noon3, game4two
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
    set_double_system(len(pols), pols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    deg = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
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
    set_double_system(dim, pols, vrblvl-1)
    svars = string_of_symbols(200, vrblvl-1)
    degrees = [degree_of_double_polynomial(k+1) for k in range(dim)]
    result = []
    for ind in range(dim):
        result.append(svars[ind]+'^'+str(degrees[ind])+' - 1;')
    return (result, solve(result))

def m_homogeneous_bezout_number(pols, vrblvl=0):
    r"""
    Given in *pols* a list of string representations of polynomials,
    in as many variables as the elements in the list,
    this function applies a heuristic to generate a partition of the
    set of unknowns to exploit the product structure of the system.
    On return are the m-homogeneous Bezout number and the partition
    of the set of unknowns.  If the partition equals the entire
    set of unknowns, then the 1-homogeneous Bezout number equals
    the total degree of the system.
    """
    if vrblvl > 0:
        print('in m_homogeneous_bezout_number, pols :')
        for pol in pols:
            print(pol)
    dim = len(pols)
    set_double_system(dim, pols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    deg = pointer(c_int32(0))
    pbuffer = create_string_buffer(b"", 4*256)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> m_homogeneous_bezout_number calls phc', end='')
    retval = phc(530, deg, pbuffer, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('an m-homogeneous Bezout number :', deg[0])
    wholepartition = int4a2str(pbuffer, vrblvl=vrblvl-1)
    endidx = wholepartition.find('\0')
    partition = wholepartition[:endidx]
    if vrblvl > 0:
        print('the partition :', partition)
    return (deg[0], partition)

def m_partition_bezout_number(pols, partition, vrblvl=0):
    r"""
    There are as many m-homogeneous Bezout numbers as there are
    partitions of the set of unknowns of a polynomial system.
    Given in *pols* the string representations of a polynomial system
    in as many variables as equations, and a string representation of
    a *partition* of the set of unknowns, this function returns the
    m-homogeneous Bezout number corresponding to the given partition.
    """
    if vrblvl > 0:
        print('in m_partition_bezout_number, pols :')
        for pol in pols:
            print(pol)
        print('the partition :', partition)
    dim = len(pols)
    set_double_system(dim, pols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    deg = pointer(c_int32(len(partition)))
    pbuffer = str2int4a(partition, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> m_partition_bezout_number calls phc', end='')
    retval = phc(531, deg, pbuffer, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the m-homogeneous Bezout number :', deg[0])
    return deg[0]

def m_homogeneous_start_system(pols, partition, checkin=True, vrblvl=0):
    r"""
    For an m-homogeneous Bezout number of a polynomial system defined by
    a *partition* of the set of unknowns, one can define a linear-product
    system that has exactly as many regular solutions as the Bezount number.
    This linear-product system can then be used as start system in a
    homotopy to compute all isolated solutions of any polynomial system
    with the same m-homogeneous structure.
    This function returns a linear-product start system with random
    coefficients and its solutions for the given polynomials in pols
    and the partition.
    If *checkin*, then the list *pols* is tested to see if *pols* defines
    a square polynomial system.  If the input system is not square,
    then an error message is printed and None is returned.
    """
    if vrblvl > 0:
        print('in m_homogeneous_start_system, pols :')
        for pol in pols:
            print(pol)
        print('the partition :', partition)
    if checkin:
        errmsg = 'Start systems are defined only for square systems,'
        if not solve_checkin(pols, errmsg):
            return None
    dim = len(pols)
    set_double_system(dim, pols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    deg = pointer(c_int32(len(partition)))
    pbuffer = str2int4a(partition, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> m_partition_start_system calls phc', end='')
    retval = phc(532, deg, pbuffer, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the m-homogeneous Bezout number :', deg[0])
    startsys = get_double_system(vrblvl-1)
    if vrblvl > 0:
        print('the start system :')
        for pol in startsys:
            print(pol)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> m_partition_start_system calls phc', end='')
        print(' to solve system')
    retval = phc(114, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    startsols = get_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('solution', idx+1, ':')
            print(sol)
    return (startsys, startsols)

def linear_product_root_count(pols, checkin=True, vrblvl=0):
    r"""
    Given in *pols* a list of string representations of polynomials,
    returns a linear-product root count based on a supporting
    set structure of the polynomials in *pols*.  This root count is
    an upper bound for the number of isolated solutions.
    """
    if vrblvl > 0:
        print('in linear_product_root_count ...')
        print('the polynomials :')
        for pol in pols:
            print(pol)
    if checkin:
        errmsg = 'Root counts are defined only for square systems,'
        if not solve_checkin(pols, errmsg):
            return None
    dim = len(pols)
    set_double_system(dim, pols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    roco = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> linear_product_root_count calls phc', end='')
    retval = phc(110, roco, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    if vrblvl > 0:
        print('-> linear_product_root_count calls phc', end='')
    retval = phc(112, roco, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('linear product root count :', roco[0])
    lprc = roco[0]
    strsets = create_string_buffer(b"", 4*1024)
    if vrblvl > 0:
        print('-> linear_product_root_count calls phc', end='')
    retval = phc(116, roco, strsets, ccc, vrb)
    wholesets = int4a2str(strsets, vrblvl=vrblvl-1)
    endidx = wholesets.find('\0')
    sets = wholesets[:endidx]
    if vrblvl > 0:
        print(', return value :', retval)
        print('supporting set structure :')
        print(sets)
    return (lprc, sets)

def random_linear_product_system(pols, checkin=True, tosolve=True, vrblvl=0):
    r"""
    Given in *pols* a list of string representations of polynomials,
    returns a random linear-product system based on a supporting
    set structure and its solutions as well (if *tosolve*).
    If *checkin*, then the list *pols* is tested to see if *pols* defines
    a square polynomial system.  If the input system is not square,
    then an error message is printed and None is returned.
    """
    if vrblvl > 0:
        print('in random_linear_product_system ...')
        print('the polynomials :')
        for pol in pols:
            print(pol)
    if checkin:
        errmsg = 'Root counts are defined only for square systems,'
        if not solve_checkin(pols, errmsg):
            return None
    dim = len(pols)
    set_double_system(dim, pols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    roco = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> random_linear_product_system calls phc', end='')
    retval = phc(110, roco, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    if vrblvl > 0:
        print('-> random_linear_product_system calls phc', end='')
    retval = phc(113, roco, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    result = get_double_system()
    if not tosolve:
        return result
    if vrblvl > 0:
        print('-> random_linear_product_system calls phc', end='')
    retval = phc(114, roco, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    sols = get_double_solutions()
    return (result, sols)

def test_total_degree(vrblvl=0):
    """
    Tests the total degree and the start system.
    """
    if vrblvl > 0:
        print('in test_total_degree ...')
    pols = noon3()
    totdeg = total_degree(pols, vrblvl)
    if vrblvl > 0:
        print('the total degree of noon3 :', totdeg)
    start, startsols = total_degree_start_system(pols, vrblvl)
    if vrblvl > 0:
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    return int(len(startsols) != 27)

def test_m_homogeneous_degree(vrblvl=0):
    """
    Tests m-homogeneous Bezout number.
    """
    if vrblvl > 0:
        print('in test_m_homogeneous_degree ...')
    pols = game4two()
    deg, partition = m_homogeneous_bezout_number(pols, vrblvl)
    fail = int(deg != 9)
    deg = m_partition_bezout_number(pols, partition, vrblvl=vrblvl)
    fail = fail + int(deg != 9)
    start, sols = m_homogeneous_start_system(pols, partition, vrblvl=vrblvl)
    if vrblvl > 0:
        print('the start system :')
        for pol in start:
            print(pol)
    fail = fail + int(len(sols) != 9)
    return fail

def test_linear_product_root_count(vrblvl=0):
    """
    Tests the linear product root count.
    """
    if vrblvl > 0:
        print('in test_linear_product_root_count ...')
    pols = noon3()
    lprc, sets = linear_product_root_count(pols, vrblvl=vrblvl)
    if vrblvl > 0:
        print('linear product root count of noon3 :', lprc)
        print('the supporting set structure :')
        print(sets)
    fail = int(lprc != 21)
    prodsys, prodsols = random_linear_product_system(pols, vrblvl=vrblvl)
    if vrblvl > 0:
        print('a random linear-product system :')
        for pol in prodsys:
            print(pol)
        print('the solutions :')
        for (idx, sol) in enumerate(prodsols):
            print('Solution', idx+1, ':')
            print(sol)
    fail = fail + int(len(prodsols) != 21)
    return fail

def main():
    """
    Runs some tests.
    """
    lvl = 1
    fail = test_total_degree(lvl)
    fail = fail + test_m_homogeneous_degree(lvl)
    fail = fail + test_linear_product_root_count(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()
