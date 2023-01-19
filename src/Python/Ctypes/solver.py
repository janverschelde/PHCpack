"""
Exports the blackbox solver.
"""
from ctypes import c_int, c_double, pointer, create_string_buffer
from version import getPHCmod, int4a2nbr, int4a2str
from polynomials import setDoubleSystem

def writeDoubleSolutions(vrblvl=0):
    """
    Writes the solutions stored in double precision.
    """
    if vrblvl > 0:
        print("-> solveDoubleSystem, vrblvl =", vrblvl)
    modPHCpack = getPHCmod()
    phc = modPHCpack._ada_use_c2phc
    aaa = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    retval = phc(31, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print('-> writeDoubleSolutions return value :', retval)
    return retval

def solveDoubleSystem(nbtasks=0, mvfocus=0, vrblvl=0):
    """
    Solves the system stored in double precision, where
    nbtasks equals the number of tasks, no multitasking if zero,
    mvfocus equals zero by default and all root counts are computed,
    otherwise, the focus is on mixed volumes and polyhedral homotopies,
    and vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print("-> solveDoubleSystem, nbtasks =", nbtasks)
    modPHCpack = getPHCmod()
    phc = modPHCpack._ada_use_c2phc
    vrb = (vrblvl > 0)
    apars = int4a2nbr([1, nbtasks, mvfocus], vrb)
    broco = create_string_buffer(1024)
    ccc = pointer(c_double(0.0))
    vlvl = c_int(vrblvl)
    retval = phc(77, apars, broco, ccc, vlvl)
    if vrblvl > 0:
        print('-> solveDoubleSystem return value :', retval)
    roco = int.from_bytes(apars[0], "big")
    result = int4a2str(broco, vrb)
    return (roco, result)

def showSolve():
    """
    Solves a simple system.
    """
    lvl = 10
    polynomials = ["x^3 + 2*x*y - x^2;", "x + y - x^3;"]
    setDoubleSystem(2, polynomials, lvl)
    nbr, roco = solveDoubleSystem(lvl)
    print("number of solutions :", nbr)
    print("root counts :\n", roco)
    writeDoubleSolutions(lvl)

if __name__=="__main__":
    showSolve()
