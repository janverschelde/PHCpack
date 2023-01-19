"""
Exports functions to compute mixed volumes and stable mixed volumes.
"""
from ctypes import c_int, c_double, pointer
from version import getPHCmod
from polynomials import setDoubleSystem

def mixedVolume(demics=True, vrblvl=0):
    """
    Returns the mixed volume of the polynomial system in double precision.
    The demics flag indicates if dynamic enumeration as implemented by
    the software DEMiCs will be used, otherwise MixedVol is called.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('-> mixedVolume, demics flag :', demics)
    modPHCpack = getPHCmod()
    phc = modPHCpack._ada_use_c2phc
    mixvol = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    if demics:
        retval = phc(843, mixvol, bbb, ccc, vrb)
    else:
        retval = phc(78, mixvol, bbb, ccc, vrb)
    if vrblvl > 0:
        print('-> mixedVolume, return value :', retval)
    return mixvol[0]

def stableMixedVolume(demics=True, vrblvl=0):
    """
    Returns the stable mixed volume of the polynomial system
    in double precision.
    The demics flag indicates if dynamic enumeration as implemented by
    the software DEMiCs will be used, otherwise MixedVol is called.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('-> stableMixedVolume, demics flag :', demics)
    modPHCpack = getPHCmod()
    phc = modPHCpack._ada_use_c2phc
    mixvol = pointer(c_int(0))
    stablemv = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    if demics:
        retval = phc(844, mixvol, stablemv, ccc, vrb)
    else:
        retval = phc(79, mixvol, stablemv, ccc, vrb)
    if vrblvl > 0:
        print('-> stableMixedVolume, return value :', retval)
    return (mixvol[0], stablemv[0])

def showMixedVolume():
    """
    Computes the mixed volume of a simple example.
    """
    lvl = 10
    polynomials = ["x^3 + 2*x*y - 1;", "x + y - 1;"]
    setDoubleSystem(2, polynomials, lvl)
    mvl = mixedVolume(True, lvl)
    print('the mixed volume by DEMiCs :', mvl)
    mvl = mixedVolume(False, lvl)
    print('the mixed volume by MixedVol :', mvl)

def showStableMixedVolume():
    """
    Computes the stable mixed volume of a simple example.
    """
    lvl = 10
    polynomials = ["x^3 + 2*x*y - x^2*y;", "x + y - x^3;"]
    setDoubleSystem(2, polynomials, lvl)
    mvl, smv = stableMixedVolume(True, lvl)
    print('the mixed volume by DEMiCs :', mvl)
    print('the stable mixed volume by DEMiCs :', smv)
    mvl, smv = stableMixedVolume(False, lvl)
    print('the mixed volume by MixedVol :', mvl)
    print('the stable mixed volume by MixedVol :', smv)

if __name__=="__main__":
    showMixedVolume()
    showStableMixedVolume()
