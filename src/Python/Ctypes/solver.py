"""
Exports the blackbox solver.
"""
from ctypes import c_int, c_double, pointer, create_string_buffer
from version import get_phcfun, int4a2nbr, int4a2str
from polynomials import set_double_system
from polynomials import set_double_double_system
from polynomials import set_quad_double_system

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
    broco = create_string_buffer(1024)
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
    broco = create_string_buffer(1024)
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
    broco = create_string_buffer(1024)
    ccc = pointer(c_double(0.0))
    vlvl = c_int(vrblvl)
    retval = phc(702, apars, broco, ccc, vlvl)
    if vrblvl > 0:
        print('-> solve_quad_double_system return value :', retval)
    roco = int.from_bytes(apars[0], "big")
    result = int4a2str(broco, vrb)
    return (roco, result)

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
    # test_double_solve()
    # test_double_double_solve()
    test_quad_double_solve()
