"""
Exports the blackbox solver for isolated solutions of square systems,
in double, double double, and quad double precision.
Also Laurent systems are accepted on input.
"""
from random import randint, uniform, seed
from math import cos, sin, pi
from ctypes import c_int32, c_double, pointer, create_string_buffer
from phcpy.version import get_phcfun, int4a2nbr, int4a2str
from phcpy.dimension import set_seed
from phcpy.polynomials import set_double_system, set_double_laurent_system
from phcpy.polynomials import set_double_double_system
from phcpy.polynomials import set_double_double_laurent_system
from phcpy.polynomials import set_quad_double_system
from phcpy.polynomials import set_quad_double_laurent_system
from phcpy.polynomials import check_semicolons
from phcpy.polynomials import is_square, number_of_symbols
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

def random_trinomials(vrblvl=0):
    """
    Returns a system of two trinomials equations for testing.
    A trinomial consists of three monomials in two variables.
    Exponents are uniform between 0 and 5 and coefficients are
    on the complex unit circle.
    """
    if vrblvl > 0:
        print('in random trinomials ...')
    exponents = [(randint(0, 5), randint(0, 5)) for _ in range(0, 6)]
    monomials = [f'x^{xe}*y^{ye}' for (xe, ye) in exponents]
    angles = [uniform(0, 2*pi) for _ in range(0, 6)]
    cff = [f'({cos(a):.14f} {sin(a):+.14f}*i)' for a in angles]
    one = '+'.join(cff[i] + '*' + monomials[i] for i in range(0, 3)) + ';'
    two = '+'.join(cff[i] + '*' + monomials[i] for i in range(3, 6)) + ';'
    return [one, two]

def real_random_trinomials(sys, vrblvl=0):
    r"""
    On input in sys are two random trinonials with complex coefficients,
    in the format what **random_trinomials()** returns.
    On return is a list of two real random trinomials with the same
    monomial structure but with random real coefficients in [-1,+1].
    """
    if vrblvl > 0:
        print('in real_random_trinomials, the system :')
        for pol in sys:
            print(pol)
    result = []
    for pol in sys:
        terms = pol.split(')')
        rpol = ''
        for i in range(1, len(terms)-1):
            rterm = terms[i].split('+')
            cff = f'{uniform(-1, +1):+.17f}'
            rpol = rpol + cff + rterm[0]
        cff = f'{uniform(-1, +1):+.17f}'
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
    clear_double_solutions(vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    apars = create_string_buffer(b"", 8)
    broco = create_string_buffer(b"", 8096)
    bpars = int4a2nbr([0, nbtasks, mvfocus], vrblvl=vrblvl-1)
    for i in range(12):
        broco[i] = bpars[i]
    ccc = pointer(c_double(0.0))
    vlvl = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> solve_double_system calls phc', end='')
    retval = phc(77, apars, broco, ccc, vlvl)
    if vrblvl > 0:
        print(', return value :', retval)
    roco = int.from_bytes(apars[0], "big")
    result = int4a2str(broco, vrblvl=vrblvl-1)
    resultend = result.find('\0')
    return (roco, result[:resultend])

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
    clear_double_double_solutions(vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    apars = create_string_buffer(b"", 8)
    broco = create_string_buffer(b"", 8096)
    bpars = int4a2nbr([0, nbtasks, mvfocus], vrblvl=vrblvl-1)
    for i in range(12):
        broco[i] = bpars[i]
    ccc = pointer(c_double(0.0))
    vlvl = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> solve_double_double_system calls phc', end='')
    retval = phc(700, apars, broco, ccc, vlvl)
    if vrblvl > 0:
        print(', return value :', retval)
    roco = int.from_bytes(apars[0], "big")
    result = int4a2str(broco, vrblvl=vrblvl-1)
    resultend = result.find('\0')
    return (roco, result[:resultend])

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
    clear_quad_double_solutions(vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    apars = create_string_buffer(b"", 8)
    broco = create_string_buffer(b"", 8096)
    bpars = int4a2nbr([0, nbtasks, mvfocus], vrblvl=vrblvl-1)
    for i in range(12):
        broco[i] = bpars[i]
    ccc = pointer(c_double(0.0))
    vlvl = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> solve_quad_double_system calls phc', end='')
    retval = phc(702, apars, broco, ccc, vlvl)
    if vrblvl > 0:
        print(', return value :', retval)
    roco = int.from_bytes(apars[0], "big")
    result = int4a2str(broco, vrblvl=vrblvl-1)
    resultend = result.find('\0')
    return (roco, result[:resultend])

def solve_double_laurent_system(nbtasks=0, mvfocus=0, vrblvl=0):
    """
    Solves the laurent system stored in double precision, where
    nbtasks equals the number of tasks, no multitasking if zero,
    mvfocus equals zero by default and all root counts are computed,
    otherwise, the focus is on mixed volumes and polyhedral homotopies,
    and vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print('in solve_double_laurent_system, nbtasks :', nbtasks)
    clear_double_solutions(vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    apars = create_string_buffer(b"", 8)
    broco = create_string_buffer(b"", 8096)
    bpars = int4a2nbr([0, nbtasks, mvfocus], vrblvl=vrblvl-1)
    for i in range(12):
        broco[i] = bpars[i]
    ccc = pointer(c_double(0.0))
    vlvl = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> solve_double_laurent_system calls phc', end='')
    retval = phc(75, apars, broco, ccc, vlvl)
    if vrblvl > 0:
        print(', return value :', retval)
    roco = int.from_bytes(apars[0], "big")
    result = int4a2str(broco, vrblvl=vrblvl-1)
    resultend = result.find('\0')
    return (roco, result[:resultend])

def solve_double_double_laurent_system(nbtasks=0, mvfocus=0, vrblvl=0):
    """
    Solves the laurent system stored in double double precision, where
    nbtasks equals the number of tasks, no multitasking if zero,
    mvfocus equals zero by default and all root counts are computed,
    otherwise, the focus is on mixed volumes and polyhedral homotopies,
    and vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print('in solve_double_double_laurent_system, nbtasks :', nbtasks)
    clear_double_double_solutions(vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    apars = create_string_buffer(b"", 8)
    broco = create_string_buffer(b"", 8096)
    bpars = int4a2nbr([0, nbtasks, mvfocus], vrblvl=vrblvl-1)
    for i in range(12):
        broco[i] = bpars[i]
    ccc = pointer(c_double(0.0))
    vlvl = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> solve_double_double_laurent_system calls phc', end='')
    retval = phc(701, apars, broco, ccc, vlvl)
    if vrblvl > 0:
        print(', return value :', retval)
    roco = int.from_bytes(apars[0], "big")
    result = int4a2str(broco, vrblvl=vrblvl-1)
    resultend = result.find('\0')
    return (roco, result[:resultend])

def solve_quad_double_laurent_system(nbtasks=0, mvfocus=0, vrblvl=0):
    """
    Solves the laurent system stored in quad double precision, where
    nbtasks equals the number of tasks, no multitasking if zero,
    mvfocus equals zero by default and all root counts are computed,
    otherwise, the focus is on mixed volumes and polyhedral homotopies,
    and vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print('in solve_quad_double_laurent_system, nbtasks :', nbtasks)
    clear_quad_double_solutions(vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    apars = create_string_buffer(b"", 8)
    broco = create_string_buffer(b"", 8096)
    bpars = int4a2nbr([0, nbtasks, mvfocus], vrblvl=vrblvl-1)
    for i in range(12):
        broco[i] = bpars[i]
    ccc = pointer(c_double(0.0))
    vlvl = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> solve_quad_double_laurent_system calls phc', end='')
    retval = phc(703, apars, broco, ccc, vlvl)
    if vrblvl > 0:
        print(', return value :', retval)
    roco = int.from_bytes(apars[0], "big")
    result = int4a2str(broco, vrblvl=vrblvl-1)
    resultend = result.find('\0')
    return (roco, result[:resultend])

def solve_checkin(pols, msg, vrblvl=0):
    r"""
    Checks whether the system defined by the list of strings in *pols*
    is square.  If so, True is returned.  Otherwise, the error message
    in the string *msg* is printed to help the user.
    """
    if vrblvl > 0:
        print('in solve_checkin ...')
    mismatch = check_semicolons(pols, vrblvl-1)
    if mismatch:
        print('Warning: mismatched number of semicolons.')
        return False
    if is_square(pols, vrblvl-1):
        return True
    print(msg)
    dim = number_of_symbols(pols, vrblvl-1)
    neq = len(pols)
    print(f'got {neq} polynomials in {dim} variables.')
    print('Either correct the input, or consider phcpy.decomposition')
    print('to solve polynomial systems that are not square.')
    return False

def solve(pols, tasks=0, mvfocus=0, precision='d',
    checkin=True, dictionary_output=False, vrblvl=0):
    r"""
    Calls the blackbox solver to compute all isolated solutions.
    To compute all solutions, also all positive dimensional solution sets,
    with a numerical irreducible decomposition, see phcpy.decomposition.
    On input in *pols* is a list of strings.
    The number of tasks for multithreading is given by *tasks*.
    The default zero value for *tasks* indicates no multithreading.
    If the *mvfocus* option is set to one, then only mixed volumes
    and polyhedral homotopies will be applied in the solver, and no
    degree bounds will be computed, as is already the case when the
    input system is genuinely laurent and has negative exponents.
    Three levels of precision are supported:

    *d*: standard double precision (1.1e-15 or 2^(-53)),

    *dd*: double double precision (4.9e-32 or 2^(-104)),

    *qd*: quad double precision (1.2e-63 or 2^(-209)).

    If *checkin* (by default), the input *pols* is checked for being square.
    If *dictionary_output*, then on return is a list of dictionaries,
    else the returned list is a list of strings.
    If *vrblvl* is larger than 0, then the names of the procedures
    called in the running of the blackbox solver will be listed.
    """
    if vrblvl > 0:
        print('in solve, tasks :', tasks, end='')
        print(', mvfocus :', mvfocus, end='')
        print(', precision :', precision)
        print('the polynomials on input :')
        for pol in pols:
            print(pol)
    if checkin:
        errmsg = 'The blackbox solver accepts only square systems,'
        if not solve_checkin(pols, errmsg):
            return None
        if tasks < 0:
            print('The number of tasks must be a nonnegative integer.')
            return None
    if precision == 'd':
        set_double_laurent_system(len(pols), pols, vrblvl=vrblvl-1)
        nbr, roco = solve_double_laurent_system(nbtasks=tasks, \
            mvfocus=mvfocus, vrblvl=vrblvl-1)
        if vrblvl > 0:
            print('nbr :', nbr, end='')
            print(', roco :')
            print(roco)
        sols = get_double_solutions(vrblvl=vrblvl-1)
        clear_double_solutions(vrblvl=vrblvl-1)
        if dictionary_output:
            return formdictlist(sols)
        return sols
    if precision == 'dd':
        set_double_double_laurent_system(len(pols), pols, \
            vrblvl=vrblvl-1)
        nbr, roco = solve_double_double_laurent_system(nbtasks=tasks, \
            mvfocus=mvfocus, vrblvl=vrblvl-1)
        if vrblvl > 0:
            print('nbr :', nbr, end='')
            print(', roco :')
            print(roco)
        sols = get_double_double_solutions(vrblvl=vrblvl-1)
        clear_double_double_solutions(vrblvl=vrblvl-1)
        if dictionary_output:
            return formdictlist(sols)
        return sols
    if precision == 'qd':
        set_quad_double_laurent_system(len(pols), pols, \
            vrblvl=vrblvl-1)
        nbr, roco = solve_quad_double_laurent_system(nbtasks=tasks, \
            mvfocus=mvfocus, vrblvl=vrblvl-1)
        if vrblvl > 0:
            print('nbr :', nbr, end='')
            print(', roco :')
            print(roco)
        sols = get_quad_double_solutions(vrblvl=vrblvl-1)
        clear_quad_double_solutions(vrblvl=vrblvl-1)
        if dictionary_output:
            return formdictlist(sols)
        return sols
    print('wrong level of precision, use d, dd, or qd')
    return None

def test_trinomial_solve(vrblvl=0):
    """
    Generates a random trinomial system and solves it.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_trinomial_solve ...')
    seed(12345678)
    pols = random_trinomials(vrblvl)
    if vrblvl > 0:
        print('two random trinomials :')
        print(pols)
    set_seed(12345678)
    (mixvol, stable_mixvol) = stable_mixed_volume(pols, vrblvl-1)
    if vrblvl > 0:
        print('its mixed volume :', mixvol)
        print('its stable mixed volume :', stable_mixvol)
    nbr, roco = solve_double_system(vrblvl=vrblvl)
    if vrblvl > 0:
        print('number of computed solutions :', nbr)
        print('the root counts :')
        print(roco)
    clear_double_solutions(vrblvl-1)
    return int(nbr != stable_mixvol) + int(nbr != 28)

def test_double_solve(vrblvl=0):
    """
    Solves a simple system in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_double_solve ...')
    polynomials = ["x^3 + 2*x*y - x^2;", "x + y - x^3;"]
    set_double_system(2, polynomials, vrblvl-1)
    nbr, roco = solve_double_system(vrblvl=vrblvl)
    if vrblvl > 0:
        print('number of solutions :', nbr)
        print('root counts :\n', roco)
        write_double_solutions(vrblvl-1)
    sols = get_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('number of retrieved solutions :', len(sols))
    clear_double_solutions(vrblvl-1)
    return int(len(sols) != 3)

def test_double_double_solve(vrblvl=0):
    """
    Solves a simple system in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_double_double_solve ...')
    polynomials = ["x^3 + 2*x*y - x^2;", "x + y - x^3;"]
    set_double_double_system(2, polynomials, vrblvl-1)
    nbr, roco = solve_double_double_system(vrblvl=vrblvl)
    if vrblvl > 0:
        print('number of solutions :', nbr)
        print('root counts :\n', roco)
        write_double_double_solutions(vrblvl-1)
    sols = get_double_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('number of retrieved solutions :', len(sols))
    clear_double_double_solutions(vrblvl-1)
    return int(len(sols) != 3) # is known bug!

def test_quad_double_solve(vrblvl=0):
    """
    Solves a simple system in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_quad_double_solve ...')
    polynomials = ["x^3 + 2*x*y - x^2;", "x + y - x^3;"]
    set_quad_double_system(2, polynomials, vrblvl-1)
    nbr, roco = solve_quad_double_system(vrblvl=vrblvl)
    if vrblvl > 0:
        print("number of solutions :", nbr)
        print("root counts :\n", roco)
        write_quad_double_solutions(vrblvl-1)
    sols = get_quad_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('number of retrieved solutions :', len(sols))
    clear_quad_double_solutions(vrblvl-1)
    return len(sols) != 3 # is known bug!

def test_double_laurent_solve(vrblvl=0):
    """
    Solves a simple laurent system in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_double_laurent_solve ...')
    polynomials = ["x^(-3) + 2*x*y - x;", "x + y^(-2) - x^3;"]
    set_double_laurent_system(2, polynomials, vrblvl-1)
    nbr, roco = solve_double_laurent_system(vrblvl=vrblvl)
    if vrblvl > 0:
        print('number of solutions :', nbr)
        print('root counts :\n', roco)
        write_double_solutions(vrblvl-1)
    sols = get_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('number of retrieved solutions :', len(sols))
    clear_double_solutions(vrblvl-1)
    return int(len(sols) != 10)

def test_double_double_laurent_solve(vrblvl=0):
    """
    Solves a simple laurent system in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_double_double_laurent_solve ...')
    polynomials = ["x^(-3) + 2*x*y - x;", "x + y^(-2) - x^3;"]
    set_double_double_laurent_system(2, polynomials, vrblvl-1)
    nbr, roco = solve_double_double_laurent_system(vrblvl=vrblvl)
    if vrblvl > 0:
        print('number of solutions :', nbr)
        print('root counts :\n', roco)
        write_double_double_solutions(vrblvl-1)
    sols = get_double_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('number of retrieved solutions :', len(sols))
    clear_double_double_solutions(vrblvl-1)
    return int(len(sols) != 10)

def test_quad_double_laurent_solve(vrblvl=0):
    """
    Solves a simple laurent system in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_quad_double_laurent_solve ...')
    polynomials = ["x^(-3) + 2*x*y - x;", "x + y^(-2) - x^3;"]
    set_quad_double_laurent_system(2, polynomials, vrblvl-1)
    nbr, roco = solve_quad_double_laurent_system(vrblvl=vrblvl)
    if vrblvl > 0:
        print('number of solutions :', nbr)
        print('root counts :\n', roco)
        write_quad_double_solutions(vrblvl-1)
    sols = get_quad_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('number of retrieved solutions :', len(sols))
    clear_quad_double_solutions(vrblvl-1)
    return int(len(sols) != 10)

def test_solve(vrblvl=0):
    """
    Tests the solve function on a simple system.
    """
    if vrblvl > 0:
        print('in test_solve ...')
    polynomials = ["x^3 + 2*x*y - x^2;", "x + y - x^3;"]
    sols = solve(polynomials, vrblvl=vrblvl)
    if vrblvl > 0:
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    fail = int(len(sols) != 3)
    return fail

def main():
    """
    Runs tests on the blackbox solver.
    """
    lvl = 1
    fail = test_trinomial_solve(lvl)
    fail = fail + test_double_solve(lvl)
    fail = fail + test_double_double_solve(lvl)
    fail = fail + test_quad_double_solve(lvl)
    fail = fail + test_double_laurent_solve(lvl)
    fail = fail + test_double_double_laurent_solve(lvl)
    fail = fail + test_quad_double_laurent_solve(lvl)
    fail = fail + test_solve(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()
