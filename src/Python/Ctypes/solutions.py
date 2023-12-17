"""
Exports operations on solutions.
"""
from ctypes import c_int, c_double, pointer, create_string_buffer
from version import get_phcfun, int4a2str
from polynomials import set_double_system
from polynomials import set_double_double_system
from polynomials import set_quad_double_system
from solver import solve_double_system, write_double_solutions
from solver import solve_double_double_system, write_double_double_solutions
from solver import solve_quad_double_system, write_quad_double_solutions

def number_double_solutions(vrblvl=0):
    """
    Returns the number of solutions in double precision.
    The vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print("-> number_double_solutions, vrblvl =", vrblvl)
    phc = get_phcfun()
    aaa = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(32, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print('-> number_double_solutions, return value :', retval)
    return bbb[0]

def number_double_double_solutions(vrblvl=0):
    """
    Returns the number of solutions in double double precision.
    The vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print("-> number_double_double_solutions, vrblvl =", vrblvl)
    phc = get_phcfun()
    aaa = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(342, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print('-> number_double_double_solutions, return value :', retval)
    return bbb[0]

def number_quad_double_solutions(vrblvl=0):
    """
    Returns the number of solutions in quad double precision.
    The vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print("-> number_quad_double_solutions, vrblvl =", vrblvl)
    phc = get_phcfun()
    aaa = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(392, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print('-> number_quad_double_solutions, return value :', retval)
    return bbb[0]

def get_next_double_solution(idx, vrblvl=0):
    """
    Returns the string representation of the next solution
    in double precision, at the index idx.
    The vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print("-> get_next_double_solution, vrblvl =", vrblvl)
    phc = get_phcfun()
    aaa = pointer(c_int(idx)) # at the given index
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(525, aaa, bbb, ccc, vrb)
    size = bbb[0]
    if vrblvl > 0:
        print('-> get_next_double_solution, size :', size)
    soldata = create_string_buffer(4*size)
    retval = phc(533, bbb, soldata, ccc, vrb)
    if vrblvl > 0:
        print("-> get_next_double_solution, return value :", retval)
    result = int4a2str(soldata, False)
    return result

def get_next_double_double_solution(idx, vrblvl=0):
    """
    Returns the string representation of the next solution
    in double double precision, at the index idx.
    The vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print("-> get_next_double_double_solution, vrblvl =", vrblvl)
    phc = get_phcfun()
    aaa = pointer(c_int(idx)) # at the given index
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(526, aaa, bbb, ccc, vrb)
    size = bbb[0]
    if vrblvl > 0:
        print('-> get_next_double_double_solution, size :', size)
    soldata = create_string_buffer(4*size)
    retval = phc(534, bbb, soldata, ccc, vrb)
    if vrblvl > 0:
        print("-> get_next_double_double_solution, return value :", retval)
    result = int4a2str(soldata, False)
    return result

def get_next_quad_double_solution(idx, vrblvl=0):
    """
    Returns the string representation of the next solution
    in quad double precision, at the index idx.
    The vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print("-> get_next_quad_double_solution, vrblvl =", vrblvl)
    phc = get_phcfun()
    aaa = pointer(c_int(idx)) # at the given index
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(527, aaa, bbb, ccc, vrb)
    size = bbb[0]
    if vrblvl > 0:
        print('-> get_next_quad_double_solution, size :', size)
    soldata = create_string_buffer(4*size)
    retval = phc(535, bbb, soldata, ccc, vrb)
    if vrblvl > 0:
        print("-> get_next_quad_double_solution, return value :", retval)
    result = int4a2str(soldata, False)
    return result

def move_double_solution_cursor(idx, vrblvl=0):
    """
    Moves the cursor to the next solution, following the index,
    in double precision.
    """
    if vrblvl > 0:
        print("-> move_double_solution_cursor, vrblvl =", vrblvl)
    phc = get_phcfun()
    aaa = pointer(c_int(idx)) # at the given index
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(454, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print("-> move_double_solution_cursor, return value :", retval)
    return aaa[0]

def move_double_double_solution_cursor(idx, vrblvl=0):
    """
    Moves the cursor to the next solution, following the index,
    in double double precision.
    """
    if vrblvl > 0:
        print("-> move_double_double_solution_cursor, vrblvl =", vrblvl)
    phc = get_phcfun()
    aaa = pointer(c_int(idx)) # at the given index
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(455, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print("-> move_double_double_solution_cursor, return value :", retval)
    return aaa[0]

def move_quad_double_solution_cursor(idx, vrblvl=0):
    """
    Moves the cursor to the next solution, following the index,
    in quad double precision.
    """
    if vrblvl > 0:
        print("-> move_quad_double_solution_cursor, vrblvl =", vrblvl)
    phc = get_phcfun()
    aaa = pointer(c_int(idx)) # at the given index
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(456, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print("-> move_quad_double_solution_cursor, return value :", retval)
    return aaa[0]

def test_double_solutions():
    """
    Solves a simple system and tests the operations
    on the solutions in double precision.
    """
    lvl = 10
    polynomials = ["x^3 + 2*x*y - x^2;", "x + y - x^3;"]
    set_double_system(2, polynomials, lvl)
    nbr, roco = solve_double_system(lvl)
    print("number of solutions :", nbr)
    print("root counts :\n", roco)
    write_double_solutions(lvl)
    nbrsols = number_double_solutions(lvl)
    print("number of solutions retrieved :", nbrsols)
    sol = get_next_double_solution(1, lvl)
    print("the first solution :")
    print(sol)
    idx = 1
    for _ in range(1, nbrsols):
        idx = move_double_solution_cursor(idx, lvl)
        print("the next index :", idx)
        sol = get_next_double_solution(idx, lvl)
        print("the solution at index", idx, ":")
        print(sol)

def test_double_double_solutions():
    """
    Solves a simple system and tests the operations
    on the solutions in double double precision.
    """
    lvl = 10
    polynomials = ["x^3 + 2*x*y - x^2;", "x + y - x^3;"]
    set_double_double_system(2, polynomials, lvl)
    nbr, roco = solve_double_double_system(lvl)
    print("number of solutions :", nbr)
    print("root counts :\n", roco)
    write_double_double_solutions(lvl)
    nbrsols = number_double_double_solutions(lvl)
    print("number of solutions retrieved :", nbrsols)
    sol = get_next_double_double_solution(1, lvl)
    print("the first solution :")
    print(sol)
    idx = 1
    for _ in range(1, nbrsols):
        idx = move_double_double_solution_cursor(idx, lvl)
        print("the next index :", idx)
        sol = get_next_double_double_solution(idx, lvl)
        print("the solution at index", idx, ":")
        print(sol)

def test_quad_double_solutions():
    """
    Solves a simple system and tests the operations
    on the solutions in quad double precision.
    """
    lvl = 10
    polynomials = ["x^3 + 2*x*y - x^2;", "x + y - x^3;"]
    set_quad_double_system(2, polynomials, lvl)
    nbr, roco = solve_quad_double_system(lvl)
    print("number of solutions :", nbr)
    print("root counts :\n", roco)
    write_quad_double_solutions(lvl)
    nbrsols = number_quad_double_solutions(lvl)
    print("number of solutions retrieved :", nbrsols)
    sol = get_next_quad_double_solution(1, lvl)
    print("the first solution :")
    print(sol)
    idx = 1
    for _ in range(1, nbrsols):
        idx = move_quad_double_solution_cursor(idx, lvl)
        print("the next index :", idx)
        sol = get_next_quad_double_solution(idx, lvl)
        print("the solution at index", idx, ":")
        print(sol)

if __name__=="__main__":
    # test_double_solutions()
    # test_double_double_solutions()
    test_quad_double_solutions()
