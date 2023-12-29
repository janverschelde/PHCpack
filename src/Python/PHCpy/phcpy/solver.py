"""
Exports the blackbox solver for isolated solutions of square systems.
"""
from ctypes import c_int, c_double, pointer, create_string_buffer
from phcpy.version import get_phcfun, int4a2nbr, int4a2str
from phcpy.dimension import set_seed, get_core_count
from phcpy.polynomials import set_double_system, set_double_Laurent_system
from phcpy.polynomials import set_double_double_system
from phcpy.polynomials import set_double_double_Laurent_system
from phcpy.polynomials import set_quad_double_system
from phcpy.polynomials import set_quad_double_Laurent_system
from phcpy.polynomials import is_square
from phcpy.solutions import clear_double_solutions
from phcpy.solutions import clear_double_double_solutions
from phcpy.solutions import clear_quad_double_solutions
from phcpy.solutions import get_double_solutions
from phcpy.solutions import get_double_double_solutions
from phcpy.solutions import get_quad_double_solutions
from phcpy.solutions import write_double_solutions
from phcpy.solutions import write_double_double_solutions
from phcpy.solutions import write_quad_double_solutions
from phcpy.solutions import formdictlist
from phcpy.volumes import stable_mixed_volume

def random_trinomials():
    """
    Returns a system of two trinomials equations for testing.
    A trinomial consists of three monomials in two variables.
    Exponents are uniform between 0 and 5 and coefficients are
    on the complex unit circle.
    """
    from random import randint as r
    exponents = [(r(0, 5), r(0, 5)) for _ in range(0, 6)]
    makemonf = lambda e: 'x^%d*y^%d' % e
    monomials = [makemonf(e) for e in exponents]
    from random import uniform as u
    from math import cos, sin, pi
    angles = [u(0, 2*pi) for _ in range(0, 6)]
    makecff = lambda a: '(' + str(cos(a)) + '%+.14f' % sin(a) + '*i)'
    cff = [makecff(a) for a in angles]
    one = '+'.join(cff[i] + '*' + monomials[i] for i in range(0, 3)) + ';'
    two = '+'.join(cff[i] + '*' + monomials[i] for i in range(3, 6)) + ';'
    return [one, two]

def real_random_trinomials(sys):
    r"""
    On input in sys are two random trinonials with complex coefficients,
    in the format what **random_trinomials()** returns.
    On return is a list of two real random trinomials with the same
    monomial structure but with random real coefficients in [-1,+1].
    """
    from random import uniform as u
    result = []
    for pol in sys:
        terms = pol.split(')')
        rpol = ''
        for i in range(1, len(terms)-1):
            rterm = terms[i].split('+')
            cff = '%+.17f' % u(-1, +1)
            rpol = rpol + cff + rterm[0]
        cff = '%+.17f' % u(-1, +1)
        rpol = rpol + cff + terms[len(terms)-1]
        result.append(rpol)
    return result

def solve_double_system(nbtasks=0, mvfocus=0, vrblvl=0):
    """
    Solves the system stored in double precision, where
    nbtasks equals the number of tasks, no multitasking if zero,
    mvfocus equals zero by default and all root counts are computed,
    otherwise, the focus is on mixed volumes and polyhedral homotopies,
    and vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print('in solve_double_system, nbtasks :', nbtasks)
    phc = get_phcfun()
    vrb = (vrblvl > 0)
    apars = create_string_buffer(b"", 8)
    broco = create_string_buffer(b"", 8096)
    bpars = int4a2nbr([0, nbtasks, mvfocus], vrb)
    for i in range(12):
        broco[i] = bpars[i]
    ccc = pointer(c_double(0.0))
    vlvl = c_int(vrblvl)
    if vrblvl > 0:
        print('-> solve_double_system calls phc', end='')
    retval = phc(77, apars, broco, ccc, vlvl)
    if vrblvl > 0:
        print(', return value :', retval)
    roco = int.from_bytes(apars[0], "big")
    result = int4a2str(broco, vrb)
    return (roco, result)

def solve_double_double_system(nbtasks=0, mvfocus=0, vrblvl=0):
    """
    Solves the system stored in double double precision, where
    nbtasks equals the number of tasks, no multitasking if zero,
    mvfocus equals zero by default and all root counts are computed,
    otherwise, the focus is on mixed volumes and polyhedral homotopies,
    and vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print("in solve_double_double_system, nbtasks :", nbtasks)
    phc = get_phcfun()
    vrb = (vrblvl > 0)
    apars = create_string_buffer(b"", 8)
    broco = create_string_buffer(b"", 8096)
    bpars = int4a2nbr([0, nbtasks, mvfocus], vrb)
    for i in range(12):
        broco[i] = bpars[i]
    ccc = pointer(c_double(0.0))
    vlvl = c_int(vrblvl)
    if vrblvl > 0:
        print('-> solve_double_double_system calls phc', end='')
    retval = phc(700, apars, broco, ccc, vlvl)
    if vrblvl > 0:
        print(', return value :', retval)
    roco = int.from_bytes(apars[0], "big")
    result = int4a2str(broco, vrb)
    return (roco, result)

def solve_quad_double_system(nbtasks=0, mvfocus=0, vrblvl=0):
    """
    Solves the system stored in quad double precision, where
    nbtasks equals the number of tasks, no multitasking if zero,
    mvfocus equals zero by default and all root counts are computed,
    otherwise, the focus is on mixed volumes and polyhedral homotopies,
    and vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print("in solve_quad_double_system, nbtasks :", nbtasks)
    phc = get_phcfun()
    vrb = (vrblvl > 0)
    apars = create_string_buffer(b"", 8)
    broco = create_string_buffer(b"", 8096)
    bpars = int4a2nbr([0, nbtasks, mvfocus], vrb)
    for i in range(12):
        broco[i] = bpars[i]
    ccc = pointer(c_double(0.0))
    vlvl = c_int(vrblvl)
    if vrblvl > 0:
        print('-> solve_quad_double_system calls phc', end='')
    retval = phc(702, apars, broco, ccc, vlvl)
    if vrblvl > 0:
        print(', return value :', retval)
    roco = int.from_bytes(apars[0], "big")
    result = int4a2str(broco, vrb)
    return (roco, result)

def solve_double_Laurent_system(nbtasks=0, mvfocus=0, vrblvl=0):
    """
    Solves the Laurent system stored in double precision, where
    nbtasks equals the number of tasks, no multitasking if zero,
    mvfocus equals zero by default and all root counts are computed,
    otherwise, the focus is on mixed volumes and polyhedral homotopies,
    and vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print('in solve_double_Laurent_system, nbtasks :', nbtasks)
    phc = get_phcfun()
    vrb = (vrblvl > 0)
    apars = create_string_buffer(b"", 8)
    broco = create_string_buffer(b"", 8096)
    bpars = int4a2nbr([0, nbtasks, mvfocus], vrb)
    for i in range(12):
        broco[i] = bpars[i]
    ccc = pointer(c_double(0.0))
    vlvl = c_int(vrblvl)
    if vrblvl > 0:
        print('-> solve_double_Laurent_system calls phc', end='')
    retval = phc(75, apars, broco, ccc, vlvl)
    if vrblvl > 0:
        print(', return value :', retval)
    roco = int.from_bytes(apars[0], "big")
    result = int4a2str(broco, vrb)
    return (roco, result)

def solve_double_double_Laurent_system(nbtasks=0, mvfocus=0, vrblvl=0):
    """
    Solves the Laurent system stored in double double precision, where
    nbtasks equals the number of tasks, no multitasking if zero,
    mvfocus equals zero by default and all root counts are computed,
    otherwise, the focus is on mixed volumes and polyhedral homotopies,
    and vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print('in solve_double_double_Laurent_system, nbtasks :', nbtasks)
    phc = get_phcfun()
    vrb = (vrblvl > 0)
    apars = create_string_buffer(b"", 8)
    broco = create_string_buffer(b"", 8096)
    bpars = int4a2nbr([0, nbtasks, mvfocus], vrb)
    for i in range(12):
        broco[i] = bpars[i]
    ccc = pointer(c_double(0.0))
    vlvl = c_int(vrblvl)
    if vrblvl > 0:
        print('-> solve_double_double_Laurent_system calls phc', end='')
    retval = phc(701, apars, broco, ccc, vlvl)
    if vrblvl > 0:
        print(', return value :', retval)
    roco = int.from_bytes(apars[0], "big")
    result = int4a2str(broco, vrb)
    return (roco, result)

def solve_quad_double_Laurent_system(nbtasks=0, mvfocus=0, vrblvl=0):
    """
    Solves the Laurent system stored in quad double precision, where
    nbtasks equals the number of tasks, no multitasking if zero,
    mvfocus equals zero by default and all root counts are computed,
    otherwise, the focus is on mixed volumes and polyhedral homotopies,
    and vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print('in solve_quad_double_Laurent_system, nbtasks :', nbtasks)
    phc = get_phcfun()
    vrb = (vrblvl > 0)
    apars = create_string_buffer(b"", 8)
    broco = create_string_buffer(b"", 8096)
    bpars = int4a2nbr([0, nbtasks, mvfocus], vrb)
    for i in range(12):
        broco[i] = bpars[i]
    ccc = pointer(c_double(0.0))
    vlvl = c_int(vrblvl)
    if vrblvl > 0:
        print('-> solve_quad_double_Laurent_system calls phc', end='')
    retval = phc(703, apars, broco, ccc, vlvl)
    if vrblvl > 0:
        print(', return value :', retval)
    roco = int.from_bytes(apars[0], "big")
    result = int4a2str(broco, vrb)
    return (roco, result)

def solve_checkin(pols, msg):
    r"""
    Checks whether the system defined by the list of strings in *pols*
    is square.  If so, True is returned.  Otherwise, the error message
    in the string *msg* is printed to help the user.
    """
    if is_square(pols):
        return True
    else:
        print(msg)
        dim = number_of_symbols(pols)
        neq = len(pols)
        print('got %d polynomials in %d variables.' % (neq, dim))
        print('Either correct the input, or use phcpy.factor.solve')
        print('to solve polynomial systems that are not square.')

def solve(pols, tasks=0, mvfocus=0, precision='d', 
    checkin=True, dictionary_output=False, verbose_level=0):
    r"""
    Calls the blackbox solver to compute all isolated solutions.
    To compute all solutions, also all positive dimensional solution sets,
    with a numerical irreducible decomposition, use solve in phcpy.factor.
    On input in *pols* is a list of strings.
    The number of tasks for multithreading is given by *tasks*.
    The default zero value for *tasks* indicates no multithreading.
    If the *mvfocus* option is set to one, then only mixed volumes
    and polyhedral homotopies will be applied in the solver, and no
    degree bounds will be computed, as is already the case when the
    input system is genuinely Laurent and has negative exponents.
    Three levels of precision are supported:

    *d*: standard double precision (1.1e-15 or 2^(-53)),

    *dd*: double double precision (4.9e-32 or 2^(-104)),

    *qd*: quad double precision (1.2e-63 or 2^(-209)).

    If *checkin* (by default), the input *pols* is checked for being square.
    If *dictionary_output*, then on return is a list of dictionaries,
    else the returned list is a list of strings.
    If *verbose_level* is larger than 0, then the names of the procedures
    called in the running of the blackbox solver will be listed.
    """
    if checkin:
        errmsg = 'The blackbox solver accepts only square systems,'
        if not solve_checkin(pols, errmsg):
            return None
        if tasks < 0:
            print('The number of tasks must be a nonnegative integer.')
            return None
    if(precision == 'd'):
        set_double_Laurent_system(len(pols), pols, vrblvl=verbose_level)
        nbr, roco = solve_double_Laurent_system(nbtasks=tasks, \
            mvfocus=mvfocus, vrblvl=verbose_level)
        sols = get_double_solutions(vrblvl=verbose_level)
        clear_double_solutions(vrblvl=verbose_level)
        if dictionary_output:
            return formdictlist(sols)
        else:
            return sols
    elif(precision == 'dd'):
        set_double_double_Laurent_system(len(pols), pols, \
            vrblvl=verbose_level)
        nbr, roco = solve_double_double_Laurent_system(nbtasks=tasks, \
            mvfocus=mvfocus, vrblvl=verbose_level)
        sols = get_double_double_solutions(vrblvl=verbose_level)
        clear_double_double_solutions(vrblvl=verbose_level)
        if dictionary_output:
            return formdictlist(sols)
        else:
            return sols
    elif(precision == 'qd'):
        set_quad_double_Laurent_system(len(pols), pols, \
            vrblvl=verbose_level)
        nbr, roco = solve_quad_double_Laurent_system(nbtasks=tasks, \
            mvfocus=mvfocus, vrblvl=verbose_level)
        sols = get_quad_double_solutions(vrblvl=verbose_level)
        clear_quad_double_solutions(vrblvl=verbose_level)
        if dictionary_output:
            return formdictlist(sols)
        else:
            return sols
    else:
        print('wrong level of precision, use d, dd, or qd')

def test_trinomial_solve(vrblvl=0):
    """
    Generates a random trinomial system and solves it.
    The verbose level is given by vrblvl.
    """
    from random import seed
    seed(12345678)
    pols = random_trinomials()
    print('two random trinomials :')
    print(pols)
    set_seed(12345678)
    set_double_system(2, pols)
    (mixvol, stable_mixvol) = stable_mixed_volume()
    print('its mixed volume :', mixvol)
    print('its stable mixed volume :', stable_mixvol)
    nbr, roco = solve_double_system()
    print('number of computed solutions :', nbr)
    clear_double_solutions(vrblvl)
    return nbr != 28

def test_double_solve(vrblvl=0):
    """
    Solves a simple system in double precision.
    The verbose level is given by vrblvl.
    """
    polynomials = ["x^3 + 2*x*y - x^2;", "x + y - x^3;"]
    set_double_system(2, polynomials, vrblvl)
    nbr, roco = solve_double_system(vrblvl)
    if vrblvl > 0:
        print('number of solutions :', nbr)
        print('root counts :\n', roco)
        write_double_solutions(vrblvl)
    sols = get_double_solutions(vrblvl)
    if vrblvl > 0:
        print('number of retrieved solutions :', len(sols))
    clear_double_solutions(vrblvl)
    return len(sols) != 3

def test_double_double_solve(vrblvl=0):
    """
    Solves a simple system in double double precision.
    The verbose level is given by vrblvl.
    """
    polynomials = ["x^3 + 2*x*y - x^2;", "x + y - x^3;"]
    set_double_double_system(2, polynomials, vrblvl)
    nbr, roco = solve_double_double_system(vrblvl)
    if vrblvl > 0:
        print('number of solutions :', nbr)
        print('root counts :\n', roco)
        write_double_double_solutions(vrblvl)
    sols = get_double_double_solutions(vrblvl)
    if vrblvl > 0:
        print('number of retrieved solutions :', len(sols))
    clear_double_double_solutions(vrblvl)
    return len(sols) != 3 # is known bug!

def test_quad_double_solve(vrblvl=0):
    """
    Solves a simple system in quad double precision.
    The verbose level is given by vrblvl.
    """
    polynomials = ["x^3 + 2*x*y - x^2;", "x + y - x^3;"]
    set_quad_double_system(2, polynomials, vrblvl)
    nbr, roco = solve_quad_double_system(vrblvl)
    if vrblvl > 0:
        print("number of solutions :", nbr)
        print("root counts :\n", roco)
        write_quad_double_solutions(vrblvl)
    sols = get_quad_double_solutions(vrblvl)
    if vrblvl > 0:
        print('number of retrieved solutions :', len(sols))
    clear_quad_double_solutions(vrblvl)
    return len(sols) != 3 # is known bug!

def test_double_Laurent_solve(vrblvl=0):
    """
    Solves a simple Laurent system in double precision.
    The verbose level is given by vrblvl.
    """
    polynomials = ["x^(-3) + 2*x*y - x;", "x + y^(-2) - x^3;"]
    set_double_Laurent_system(2, polynomials, vrblvl)
    nbr, roco = solve_double_Laurent_system(vrblvl)
    if vrblvl > 0:
        print('number of solutions :', nbr)
        print('root counts :\n', roco)
        write_double_solutions(vrblvl)
    sols = get_double_solutions(vrblvl)
    if vrblvl > 0:
        print('number of retrieved solutions :', len(sols))
    clear_double_solutions(vrblvl)
    return len(sols) != 10

def test_double_double_Laurent_solve(vrblvl=0):
    """
    Solves a simple Laurent system in double double precision.
    The verbose level is given by vrblvl.
    """
    polynomials = ["x^(-3) + 2*x*y - x;", "x + y^(-2) - x^3;"]
    set_double_double_Laurent_system(2, polynomials, vrblvl)
    nbr, roco = solve_double_double_Laurent_system(vrblvl)
    if vrblvl > 0:
        print('number of solutions :', nbr)
        print('root counts :\n', roco)
        write_double_double_solutions(vrblvl)
    sols = get_double_double_solutions(vrblvl)
    if vrblvl > 0:
        print('number of retrieved solutions :', len(sols))
    clear_double_double_solutions(vrblvl)
    return len(sols) != 10

def test_quad_double_Laurent_solve(vrblvl=0):
    """
    Solves a simple Laurent system in quad double precision.
    The verbose level is given by vrblvl.
    """
    polynomials = ["x^(-3) + 2*x*y - x;", "x + y^(-2) - x^3;"]
    set_quad_double_Laurent_system(2, polynomials, vrblvl)
    nbr, roco = solve_quad_double_Laurent_system(vrblvl)
    if vrblvl > 0:
        print('number of solutions :', nbr)
        print('root counts :\n', roco)
        write_quad_double_solutions(vrblvl)
    sols = get_quad_double_solutions(vrblvl)
    if vrblvl > 0:
        print('number of retrieved solutions :', len(sols))
    clear_quad_double_solutions(vrblvl)
    return len(sols) != 10

def test_solve(vrblvl=0):
    """
    Tests the solve function on a simple system.
    """
    polynomials = ["x^3 + 2*x*y - x^2;", "x + y - x^3;"]
    sols = solve(polynomials, verbose_level=vrblvl)
    if vrblvl > 0:
        print('the solutions :')
        for sol in sols:
            print(sol)
    fail = int(len(sols) != 3)
    """
    polynomials = ["x^(-3) + 2*x*y - x;", "x + y^(-2) - x^3;"]
    sols = solve(polynomials, dictionary_output=True, \
        verbose_level=vrblvl)
    if vrblvl > 0:
        print('the solutions :')
        for (k, sol) in zip(range(len(sols)), sols):
            print('solution', k+1, ':')
            print(sol)
    fail = fail + int(len(sols) != 10)
    polynomials = ["x^(-3) + 2*x*y - x;", "x + y^(-2) - x^3;"]
    nbr_cores = get_core_count()
    if vrblvl > 0:
        print('number of tasks :', nbr_cores)
    sols = solve(polynomials, tasks=nbr_cores, precision='dd', \
        verbose_level=vrblvl)
    if vrblvl > 0:
        print('the solutions :')
        for (k, sol) in zip(range(len(sols)), sols):
            print('solution', k+1, ':')
            print(sol)
    fail = fail + int(len(sols) != 10)
    polynomials = ["x^(-3) + 2*x*y - x;", "x + y^(-2) - x^3;"]
    nbr_cores = get_core_count()
    if vrblvl > 0:
        print('number of tasks :', nbr_cores)
    sols = solve(polynomials, tasks=nbr_cores, precision='qd', \
        verbose_level=vrblvl)
    if vrblvl > 0:
        print('the solutions :')
        for (k, sol) in zip(range(len(sols)), sols):
            print('solution', k+1, ':')
            print(sol)
    fail = fail + int(len(sols) != 10)
    """
    return fail

def main():
    """
    Runs tests on the blackbox solver.
    """
    lvl = 10
    fail = test_trinomial_solve(lvl)
    fail = fail + test_double_solve(lvl)
    fail = fail + test_double_double_solve(lvl)
    fail = fail + test_quad_double_solve(lvl)
    #fail = fail + test_double_Laurent_solve(lvl)
    #fail = fail + test_double_double_Laurent_solve(lvl)
    #fail = fail + test_quad_double_Laurent_solve(lvl)
    # fail = fail + test_solve(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()
