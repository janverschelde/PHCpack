"""
Exports the blackbox solver.
"""

def writeDoubleSolutions(vrblvl=0):
    """
    Writes the solutions stored in double precision.
    """
    from ctypes import c_int, c_double, pointer, create_string_buffer
    from version import getPHCmod, int4a2nbr, int4a2str
    if(vrblvl > 0):
        print("-> solveDoubleSystem, vrblvl =", vrblvl)
    modPHCpack = getPHCmod()
    f = modPHCpack._ada_use_c2phc
    vrb = (vrblvl > 0)
    a = pointer(c_int(0))
    b = pointer(c_int(0))
    c = pointer(c_double(0.0))
    v = c_int(vrblvl)
    r = f(31, a, b, c, v)
    if(vrblvl > 0):
        print('-> writeDoubleSolutions return value :', r)
    return r

def solveDoubleSystem(nbtasks=0, mvfocus=0, vrblvl=0):
    """
    Solves the system stored in double precision, where
    nbtasks equals the number of tasks, no multitasking if zero,
    mvfocus equals zero by default and all root counts are computed,
    otherwise, the focus is on mixed volumes and polyhedral homotopies,
    and vrblvl is the verbose level.
    """
    from ctypes import c_int, c_double, pointer, create_string_buffer
    from version import getPHCmod, int4a2nbr, int4a2str
    if(vrblvl > 0):
        print("-> solveDoubleSystem, nbtasks =", nbtasks)
    modPHCpack = getPHCmod()
    f = modPHCpack._ada_use_c2phc
    vrb = (vrblvl > 0)
    a = int4a2nbr([1, nbtasks, mvfocus], vrb)
    print('a[0] :', a[0])
    b = create_string_buffer(1024)
    c = pointer(c_double(0.0))
    v = c_int(vrblvl)
    r = f(77, a, b, c, v)
    if(vrblvl > 0):
        print('-> solveDoubleSystem return value :', r)
    print('a[0] :', a[0])
    roco = a[0].decode()
    result = int4a2str(b, vrb)
    return (roco, result)

def showSolve():
    """
    Solves a simple system.
    """
    from polynomials import setDoubleSystem
    lvl = 10
    polynomials = ["x^3 + 2*x*y - x^2;", "x + y - x^3;"]
    r = setDoubleSystem(2, polynomials, lvl)
    nbr, roco = solveDoubleSystem(lvl)
    print("number of solutions : ", nbr)
    print("root counts :\n", roco)
    writeDoubleSolutions(lvl)

if __name__=="__main__":
    showSolve()
