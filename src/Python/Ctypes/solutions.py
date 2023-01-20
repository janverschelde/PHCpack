"""
Exports operations on solutions.
"""
from ctypes import c_int, c_double, pointer, create_string_buffer
from version import getPHCmod, int4a2nbr, str2int4a, int4a2str
from polynomials import setDoubleSystem
from solver import solveDoubleSystem, writeDoubleSolutions

def numberDoubleSolutions(vrblvl=0):
    """
    Returns the number of solutions in double precision.
    The vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print("-> numberDoubleSolutions, vrblvl =", vrblvl)
    modPHCpack = getPHCmod()
    phc = modPHCpack._ada_use_c2phc
    aaa = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(32, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print('-> numberDoubleSolutions, return value :', retval)
    return bbb[0]

def showSolutions():
    """
    Solves a simple system and tests the operations
    on the solutions in double precision.
    """
    lvl = 10
    polynomials = ["x^3 + 2*x*y - x^2;", "x + y - x^3;"]
    setDoubleSystem(2, polynomials, lvl)
    nbr, roco = solveDoubleSystem(lvl)
    print("number of solutions :", nbr)
    print("root counts :\n", roco)
    writeDoubleSolutions(lvl)
    nbrsols = numberDoubleSolutions(lvl)
    print("number of solutions retrieved :", nbrsols)

if __name__=="__main__":
    showSolutions()
