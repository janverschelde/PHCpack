"""
Exports operations on solutions.
"""
from ctypes import c_int, c_double, pointer
from version import getPHCmod
from polynomials import set_double_system
from solver import solve_double_system, write_double_solutions

def number_double_solutions(vrblvl=0):
    """
    Returns the number of solutions in double precision.
    The vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print("-> number_double_solutions, vrblvl =", vrblvl)
    phcpack = getPHCmod()
    phc = phcpack._ada_use_c2phc
    aaa = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(32, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print('-> number_double_solutions, return value :', retval)
    return bbb[0]

def show_solutions():
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

if __name__=="__main__":
    show_solutions()
