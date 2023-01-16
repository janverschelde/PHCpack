"""
Exports functions to compute mixed volumes and stable mixed volumes.
"""

def mixedVolume(demics=True, vrblvl=0):
    """
    Returns the mixed volume of the polynomial system in double precision.
    The demics flag indicates if dynamic enumeration as implemented by
    the software DEMiCs will be used, otherwise MixedVol is called.
    The verbose level is given by vrblvl.
    """
    from ctypes import c_int, c_double, pointer
    from version import getPHCmod
    if(vrblvl > 0):
        print('-> mixedVolume, demics flag :', demics)
    modPHCpack = getPHCmod()
    f = modPHCpack._ada_use_c2phc
    a = pointer(c_int(0))
    b = pointer(c_int(0))
    c = pointer(c_double(0.0))
    v = c_int(vrblvl)
    if demics:
        r = f(843, a, b, c, v)
    else:
        r = f(78, a, b, c, v)
    if(vrblvl > 0):
        print('-> mixedVolume, return value :', r)
    return a[0]

def stableMixedVolume(demics=True, vrblvl=0):
    """
    Returns the stable mixed volume of the polynomial system
    in double precision.
    The demics flag indicates if dynamic enumeration as implemented by
    the software DEMiCs will be used, otherwise MixedVol is called.
    The verbose level is given by vrblvl.
    """
    from ctypes import c_int, c_double, pointer
    from version import getPHCmod
    if(vrblvl > 0):
        print('-> stableMixedVolume, demics flag :', demics)
    modPHCpack = getPHCmod()
    f = modPHCpack._ada_use_c2phc
    a = pointer(c_int(0))
    b = pointer(c_int(0))
    c = pointer(c_double(0.0))
    v = c_int(vrblvl)
    if demics:
        r = f(844, a, b, c, v)
    else:
        r = f(79, a, b, c, v)
    if(vrblvl > 0):
        print('-> stableMixedVolume, return value :', r)
    return (a[0], b[0])

def showMixedVolume():
    """
    Computes the mixed volume of a simple example.
    """
    from polynomials import setDoubleSystem
    lvl = 10
    polynomials = ["x^3 + 2*x*y - 1;", "x + y - 1;"]
    setDoubleSystem(2, polynomials, lvl)
    mv = mixedVolume(True, lvl)
    print('the mixed volume by DEMiCs :', mv)
    mv = mixedVolume(False, lvl)
    print('the mixed volume by MixedVol :', mv)

def showStableMixedVolume():
    """
    Computes the stable mixed volume of a simple example.
    """
    from polynomials import setDoubleSystem
    lvl = 10
    polynomials = ["x^3 + 2*x*y - x^2*y;", "x + y - x^3;"]
    setDoubleSystem(2, polynomials, lvl)
    mv, smv = stableMixedVolume(True, lvl)
    print('the mixed volume by DEMiCs :', mv)
    print('the stable mixed volume by DEMiCs :', smv)
    mv, smv = stableMixedVolume(False, lvl)
    print('the mixed volume by MixedVol :', mv)
    print('the stable mixed volume by MixedVol :', smv)

if __name__=="__main__":
    showMixedVolume()
    showStableMixedVolume()
