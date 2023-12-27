"""
Exports functions to compute mixed volumes and stable mixed volumes.
"""
from ctypes import c_int, c_double, pointer
from version import get_phcfun
from polynomials import set_double_system

def mixed_volume(demics=True, vrblvl=0):
    """
    Returns the mixed volume of the polynomial system in double precision.
    The demics flag indicates if dynamic enumeration as implemented by
    the software DEMiCs will be used, otherwise MixedVol is called.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in mixed_volume, demics flag :', demics)
    phc = get_phcfun()
    mixvol = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    if vrblvl > 0:
        print('-> mixed_volume calls phc', end='')
    if demics:
        retval = phc(843, mixvol, bbb, ccc, vrb)
    else:
        retval = phc(78, mixvol, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return mixvol[0]

def stable_mixed_volume(demics=True, vrblvl=0):
    """
    Returns the mixed and the stable mixed volume of the polynomial system
    in double precision.
    The demics flag indicates if dynamic enumeration as implemented by
    the software DEMiCs will be used, otherwise MixedVol is called.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in stable_mixed_volume, demics flag :', demics)
    phc = get_phcfun()
    mixvol = pointer(c_int(0))
    stablemv = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    if vrblvl > 0:
        print('-> stable_mixed_volume calls phc', end='')
    if demics:
        retval = phc(844, mixvol, stablemv, ccc, vrb)
    else:
        retval = phc(79, mixvol, stablemv, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return (mixvol[0], stablemv[0])

def test_mixed_volume(vrblvl=0):
    """
    Computes the mixed volume of a simple example.
    The verbose level is given by vrblvl.
    """
    polynomials = ["x^3 + 2*x*y - 1;", "x + y - 1;"]
    set_double_system(2, polynomials, vrblvl)
    mvl = mixed_volume(True, vrblvl)
    if vrblvl > 0:
        print('the mixed volume by DEMiCs :', mvl)
    mvl = mixed_volume(False, vrblvl)
    if vrblvl > 0:
        print('the mixed volume by MixedVol :', mvl)
    return int(mvl != 3)

def test_stable_mixed_volume(vrblvl=0):
    """
    Computes the stable mixed volume of a simple example.
    The verbose level is given by vrblvl.
    """
    polynomials = ["x^3 + 2*x*y - x^2*y;", "x + y - x^3;"]
    set_double_system(2, polynomials, vrblvl)
    mvl, smv = stable_mixed_volume(True, vrblvl)
    if vrblvl > 0:
        print('the mixed volume by DEMiCs :', mvl)
        print('the stable mixed volume by DEMiCs :', smv)
    mvl, smv = stable_mixed_volume(False, vrblvl)
    if vrblvl > 0:
        print('the mixed volume by MixedVol :', mvl)
        print('the stable mixed volume by MixedVol :', smv)
    return int(mvl != 3) + int(smv != 5)

def main():
    """
    Runs tests on mixed volumes and stable mixed volumes.
    """
    lvl = 1
    fail = test_mixed_volume(lvl)
    fail = fail + test_stable_mixed_volume(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()
