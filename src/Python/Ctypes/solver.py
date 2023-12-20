"""
Exports the blackbox solver.
"""
from ctypes import c_int, c_double, pointer, create_string_buffer
from version import get_phcfun, int4a2nbr, int4a2str
from polynomials import set_double_system
from polynomials import set_double_double_system
from polynomials import set_quad_double_system
from volumes import stable_mixed_volume

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

def write_double_solutions(vrblvl=0):
    """
    Writes the solutions stored in double precision.
    """
    if vrblvl > 0:
        print("-> write_double_solutions, vrblvl =", vrblvl)
    phc = get_phcfun()
    aaa = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(31, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print('-> write_double_solutions return value :', retval)
    return retval

def write_double_double_solutions(vrblvl=0):
    """
    Writes the solutions stored in double double precision.
    """
    if vrblvl > 0:
        print("-> write_double_double_solutions, vrblvl =", vrblvl)
    phc = get_phcfun()
    aaa = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(341, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print('-> write_double_double_solutions return value :', retval)
    return retval

def write_quad_double_solutions(vrblvl=0):
    """
    Writes the solutions stored in quad double precision.
    """
    if vrblvl > 0:
        print("-> write_quad_double_solutions, vrblvl =", vrblvl)
    phc = get_phcfun()
    aaa = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(391, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print('-> write_quad_double_solutions return value :', retval)
    return retval

def solve_double_system(nbtasks=0, mvfocus=0, vrblvl=0):
    """
    Solves the system stored in double precision, where
    nbtasks equals the number of tasks, no multitasking if zero,
    mvfocus equals zero by default and all root counts are computed,
    otherwise, the focus is on mixed volumes and polyhedral homotopies,
    and vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print("-> solve_double_system, nbtasks =", nbtasks)
    phc = get_phcfun()
    vrb = (vrblvl > 0)
    apars = int4a2nbr([1, nbtasks, mvfocus], vrb)
    broco = create_string_buffer(b"", 2048)
    ccc = pointer(c_double(0.0))
    vlvl = c_int(vrblvl)
    retval = phc(77, apars, broco, ccc, vlvl)
    if vrblvl > 0:
        print('-> solve_double_system return value :', retval)
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
        print("-> solve_double_double_system, nbtasks =", nbtasks)
    phc = get_phcfun()
    vrb = (vrblvl > 0)
    apars = int4a2nbr([1, nbtasks, mvfocus], vrb)
    broco = create_string_buffer(b"", 2048)
    ccc = pointer(c_double(0.0))
    vlvl = c_int(vrblvl)
    retval = phc(700, apars, broco, ccc, vlvl)
    if vrblvl > 0:
        print('-> solve_double_double_system return value :', retval)
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
        print("-> solve_quad_double_system, nbtasks =", nbtasks)
    phc = get_phcfun()
    vrb = (vrblvl > 0)
    apars = int4a2nbr([1, nbtasks, mvfocus], vrb)
    broco = create_string_buffer(b"", 2048)
    ccc = pointer(c_double(0.0))
    vlvl = c_int(vrblvl)
    retval = phc(702, apars, broco, ccc, vlvl)
    if vrblvl > 0:
        print('-> solve_quad_double_system return value :', retval)
    roco = int.from_bytes(apars[0], "big")
    result = int4a2str(broco, vrb)
    return (roco, result)

def test_trinomial_solve():
    """
    Generates a random trinomial system and solves it.
    """
    pols = random_trinomials()
    print('two random trinomials :')
    print(pols)
    set_double_system(2, pols)
    (mixvol, stable_mixvol) = stable_mixed_volume()
    print('its mixed volume :', mixvol)
    print('its stable mixed volume :', stable_mixvol)
    nbr, roco = solve_double_system()
    print('number of computed solutions :', nbr)

def test_double_solve():
    """
    Solves a simple system in double precision.
    """
    lvl = 10
    polynomials = ["x^3 + 2*x*y - x^2;", "x + y - x^3;"]
    set_double_system(2, polynomials, lvl)
    nbr, roco = solve_double_system(lvl)
    print("number of solutions :", nbr)
    print("root counts :\n", roco)
    write_double_solutions(lvl)

def test_double_double_solve():
    """
    Solves a simple system in double double precision.
    """
    lvl = 10
    polynomials = ["x^3 + 2*x*y - x^2;", "x + y - x^3;"]
    set_double_double_system(2, polynomials, lvl)
    nbr, roco = solve_double_double_system(lvl)
    print("number of solutions :", nbr)
    print("root counts :\n", roco)
    write_double_double_solutions(lvl)

def test_quad_double_solve():
    """
    Solves a simple system in quad double precision.
    """
    lvl = 10
    polynomials = ["x^3 + 2*x*y - x^2;", "x + y - x^3;"]
    set_quad_double_system(2, polynomials, lvl)
    nbr, roco = solve_quad_double_system(lvl)
    print("number of solutions :", nbr)
    print("root counts :\n", roco)
    write_quad_double_solutions(lvl)

if __name__=="__main__":
    test_trinomial_solve()
    # test_double_solve()
    # test_double_double_solve()
    # test_quad_double_solve()
