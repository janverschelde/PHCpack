"""
Exports functions on solutions.
"""
from ctypes import c_int, c_double, pointer, create_string_buffer
from version import get_phcfun, int4a2nbr, int4a2str, str2int4a
from polynomials import set_double_system
from polynomials import set_double_double_system
from polynomials import set_quad_double_system
from solver import solve_double_system, write_double_solutions
from solver import solve_double_double_system, write_double_double_solutions
from solver import solve_quad_double_system, write_quad_double_solutions

def append_double_solution_string(nvr, sol, vrblvl=0):
    """
    Appends the string in sol to the solutions in double precision,
    where the number of variables equals nvr.
    """
    if vrblvl > 0:
        print('-> append_double_solution_string, nvr =', nvr)
        print('The solution string :')
        print(sol)
    phc = get_phcfun()
    apars = int4a2nbr([nvr, len(sol)], (vrblvl > 0))
    bsol = str2int4a(sol)
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(208, apars, bsol, ccc, vrblvl)
    if vrblvl > 0:
        print('-> append_double_solution_string, return value :', retval)
    return retval

def append_double_double_solution_string(nvr, sol, vrblvl=0):
    """
    Appends the string in sol to the solutions in double double precision,
    where the number of variables equals nvr.
    """
    if vrblvl > 0:
        print('-> append_double_double_solution_string, nvr =', nvr)
        print('The solution string :')
        print(sol)
    phc = get_phcfun()
    apars = int4a2nbr([nvr, len(sol)], vrblvl)
    bsol = str2int4a(sol)
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(378, apars, bsol, ccc, vrblvl)
    if vrblvl > 0:
        print('-> append_double_double_solution_string, return value :', retval)
    return retval

def append_quad_double_solution_string(nvr, sol, vrblvl=0):
    """
    Appends the string in sol to the solutions in quad double precision,
    where the number of variables equals nvr.
    """
    if vrblvl > 0:
        print('-> append_quad_double_solution_string, nvr =', nvr)
        print('The solution string :')
        print(sol)
    phc = get_phcfun()
    apars = int4a2nbr([nvr, len(sol)], vrblvl)
    bsol = str2int4a(sol)
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(428, apars, bsol, ccc, vrblvl)
    if vrblvl > 0:
        print('-> append_quad_double_solution_string, return value :', retval)
    return retval

def set_double_solutions(nvr, sols, vrblvl=0):
    """
    Sets the solutions in double precision, with the strings in sols,
    where the number of variables equals nvr.
    """
    if vrblvl > 0:
        print('-> set_double_solutions, nvr =', nvr)
    clear_double_solutions(vrblvl)
    fail = 0
    for ind in range(0, len(sols)):
        print('calling append_double_solution_string ...')
        print('the solution :')
        print(sols[ind])
        print('nvr =', nvr) 
        print('len(sol) =', len(sols[ind])) 
        fail = append_double_solution_string(nvr, sols[ind], vrblvl)
        if(fail != 0):
            print('Solution at position', ind, 'is not appended.')
    return fail

def set_double_double_solutions(nvr, sols, vrblvl=0):
    """
    Sets the solutions in double double precision, with the strings in sols,
    where the number of variables equals nvr.
    """
    clear_double_double_solutions(vrblvl)
    fail = 0
    for ind in range(0, len(sols)):
        fail = append_double_double_solution_string(nvr, sols[ind], vrblvl)
        if(fail != 0):
            print('Solution at position', ind, 'is not appended.')
    return fail

def set_quad_double_solutions(nvr, sols, vrblvl=0):
    """
    Sets the solutions in quad double precision, with the strings in sols,
    where the number of variables equals nvr.
    """
    clear_quad_double_solutions(vrblvl)
    fail = 0
    for ind in range(0, len(sols)):
        fail = append_quad_double_solution_string(nvr, sols[ind], vrblvl)
        if(fail != 0):
            print('Solution at position', ind, 'is not appended.')
    return fail

def number_double_solutions(vrblvl=0):
    """
    Returns the number of solutions in double precision.
    The vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print('-> number_double_solutions', end='')
    phc = get_phcfun()
    aaa = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(32, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return bbb[0]

def number_double_double_solutions(vrblvl=0):
    """
    Returns the number of solutions in double double precision.
    The vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print('-> number_double_double_solutions', end='')
    phc = get_phcfun()
    aaa = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(342, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return bbb[0]

def number_quad_double_solutions(vrblvl=0):
    """
    Returns the number of solutions in quad double precision.
    The vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print('-> number_quad_double_solutions', end='')
    phc = get_phcfun()
    aaa = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(392, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return bbb[0]

def get_next_double_solution(idx, vrblvl=0):
    """
    Returns the string representation of the next solution
    in double precision, at the index idx.
    The vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print('-> get_next_double_solution, idx =', idx)
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
        print('-> get_next_double_double_solution, idx =', idx)
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
        print('-> get_next_quad_double_solution, idx =', idx)
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
        print('-> move_double_solution_cursor, idx =', idx)
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
        print('-> move_double_double_solution_cursor, idx =', idx)
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
        print('-> move_quad_double_solution_cursor, idx =', idx)
    phc = get_phcfun()
    aaa = pointer(c_int(idx)) # at the given index
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(456, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print("-> move_quad_double_solution_cursor, return value :", retval)
    return aaa[0]

def clear_double_solutions(vrblvl=0):
    """
    Clears the solutions defined in double precision.
    """
    if vrblvl > 0:
        print('-> clear_double_solutions', end='')
    phc = get_phcfun()
    aaa = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(37, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_double_double_solutions(vrblvl=0):
    """
    Clears the solutions defined in double double precision.
    """
    if vrblvl > 0:
        print('-> clear_double_double_solutions', end='')
    phc = get_phcfun()
    aaa = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(347, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_quad_double_solutions(vrblvl=0):
    """
    Clears the solutions defined in quad double precision.
    """
    if vrblvl > 0:
        print('-> clear_quad_double_solutions', end='')
    phc = get_phcfun()
    aaa = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(397, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def get_double_solutions(vrblvl=0):
    """
    Returns the solution strings in double precision.
    """
    nbrsols = number_double_solutions(vrblvl)
    if vrblvl > 0:
        print("number of solutions retrieved :", nbrsols)
    result = []
    if nbrsols > 0:
        sol = get_next_double_solution(1, vrblvl)
        if vrblvl > 0:
            print("the first solution :")
            print(sol)
        result.append(sol)
        idx = 1
        for _ in range(1, nbrsols):
            idx = move_double_solution_cursor(idx, vrblvl)
            if vrblvl > 0:
                print("the next index :", idx)
            sol = get_next_double_solution(idx, vrblvl)
            if vrblvl > 0:
                print("the solution at index", idx, ":")
                print(sol)
        result.append(sol)
    return result

def get_double_double_solutions(vrblvl=0):
    """
    Returns the solution strings in double double precision.
    """
    nbrsols = number_double_double_solutions(vrblvl)
    if vrblvl > 0:
        print("number of solutions retrieved :", nbrsols)
    result = []
    if nbrsols > 0:
        sol = get_next_double_double_solution(1, vrblvl)
        if vrblvl > 0:
            print("the first solution :")
            print(sol)
        result.append(sol)
        idx = 1
        for _ in range(1, nbrsols):
            idx = move_double_double_solution_cursor(idx, vrblvl)
            if vrblvl > 0:
                print("the next index :", idx)
            sol = get_next_double_double_solution(idx, vrblvl)
            if vrblvl > 0:
                print("the solution at index", idx, ":")
                print(sol)
        result.append(sol)
    return result

def get_quad_double_solutions(vrblvl=0):
    """
    Returns the solution strings in quad double precision.
    """
    nbrsols = number_quad_double_solutions(vrblvl)
    if vrblvl > 0:
        print("number of solutions retrieved :", nbrsols)
    result = []
    if nbrsols > 0:
        sol = get_next_quad_double_solution(1, vrblvl)
        if vrblvl > 0:
            print("the first solution :")
            print(sol)
        result.append(sol)
        idx = 1
        for _ in range(1, nbrsols):
            idx = move_quad_double_solution_cursor(idx, vrblvl)
            if vrblvl > 0:
                print("the next index :", idx)
            sol = get_next_quad_double_solution(idx, vrblvl)
            if vrblvl > 0:
                print("the solution at index", idx, ":")
                print(sol)
        result.append(sol)
    return result

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
    sols = get_double_solutions(lvl)

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
    sols = get_double_double_solutions(lvl)

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
    sols = get_quad_double_solutions(lvl)

if __name__=="__main__":
    test_double_solutions()
    # test_double_double_solutions()
    # test_quad_double_solutions()
