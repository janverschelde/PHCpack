"""
Exports functions to compute mixed volumes and stable mixed volumes.
"""
from ctypes import c_int32, c_double, pointer
from phcpy.version import get_phcfun
from phcpy.polynomials import set_double_system

def mixed_volume(demics=True, vrblvl=0):
    """
    Returns the mixed volume of the polynomial system in double precision.
    The demics flag indicates if dynamic enumeration as implemented by
    the software DEMiCs will be used, otherwise MixedVol is called.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in mixed_volume, demics flag :', demics)
    phc = get_phcfun(vrblvl-1)
    mixvol = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
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
    phc = get_phcfun(vrblvl-1)
    mixvol = pointer(c_int32(0))
    stablemv = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> stable_mixed_volume calls phc', end='')
    if demics:
        retval = phc(844, mixvol, stablemv, ccc, vrb)
    else:
        retval = phc(79, mixvol, stablemv, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return (mixvol[0], stablemv[0])

def clear_cells(vrblvl=0):
    """
    Clears the computed mixed cells.
    """
    if vrblvl > 0:
        print('in clear_cells ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_cells calls phc', end='')
    retval = phc(94, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def number_of_cells(vrblvl=0):
    """
    Returns the number of cells computed by the mixed_volume function.
    """
    if vrblvl > 0:
        print('in number of cells ...')
    phc = get_phcfun(vrblvl-1)
    alen = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> number_of_cells calls phc', end='')
    retval = phc(82, alen, bbb, ccc, vrb)
    result = alen[0]
    if vrblvl > 0:
        print(', return value :', retval)
        print('the number of cells :', result)
    return result

def cell_mixed_volume(idx, vrblvl=0):
    """
    Returns the mixed volume of cell with index idx.
    """
    if vrblvl > 0:
        print('in cell_mixed_volume, idx :', idx)
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(idx))
    bmv = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> cell_mixed_volume calls phc', end='')
    retval = phc(90, aidx, bmv, ccc, vrb)
    result = bmv[0]
    if vrblvl > 0:
        print(', return value :', retval)
        print('the mixed volume of cell', idx, 'is', result)
    return result

def test_mixed_volume(vrblvl=0):
    """
    Computes the mixed volume of a simple example.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_mixed_volume ...')
    polynomials = ["x^3 + 2*x*y - 1;", "x + y - 1;"]
    set_double_system(2, polynomials, vrblvl)
    mvl = mixed_volume(True, vrblvl)
    nbr = number_of_cells(vrblvl)
    if vrblvl > 0:
        print('the mixed volume by DEMiCs :', mvl)
        print('the number of cells :', nbr)
    fail = int(mvl != 3)
    mvcells = 0
    for idx in range(1, nbr+1):
        mvcells = mvcells + cell_mixed_volume(idx, vrblvl)
    if vrblvl > 0:
        print('sum of mixed volumes of cells :', mvcells)
    fail = fail + int(mvcells != 3)
    clear_cells(vrblvl)
    mvl = mixed_volume(False, vrblvl)
    nbr = number_of_cells(vrblvl)
    if vrblvl > 0:
        print('the mixed volume by MixedVol :', mvl)
        print('the number of cells :', nbr)
    fail = fail + int(mvl != 3)
    mvcells = 0
    for idx in range(1, nbr+1):
        mvcells = mvcells + cell_mixed_volume(idx, vrblvl)
    if vrblvl > 0:
        print('sum of mixed volumes of cells :', mvcells)
    fail = fail + int(mvcells != 3)
    return fail

def test_stable_mixed_volume(vrblvl=0):
    """
    Computes the stable mixed volume of a simple example.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_stable_mixed_volume ...')
    polynomials = ["x^3 + 2*x*y - x^2*y;", "x + y - x^3;"]
    set_double_system(2, polynomials, vrblvl)
    mvl, smv = stable_mixed_volume(True, vrblvl)
    nbr = number_of_cells(vrblvl)
    if vrblvl > 0:
        print('the mixed volume by DEMiCs :', mvl)
        print('the stable mixed volume by DEMiCs :', smv)
        print('the number of cells :', nbr)
    mvcells = 0
    for idx in range(1, nbr+1):
        mvcells = mvcells + cell_mixed_volume(idx, vrblvl)
    if vrblvl > 0:
        print('sum of mixed volumes of cells :', mvcells)
    fail = int(mvl != 3) + int(smv != 5) + int(mvcells != 5)
    clear_cells(vrblvl)
    mvl, smv = stable_mixed_volume(False, vrblvl)
    nbr = number_of_cells(vrblvl)
    if vrblvl > 0:
        print('the mixed volume by MixedVol :', mvl)
        print('the stable mixed volume by MixedVol :', smv)
        print('the number of cells :', nbr)
    mvcells = 0
    for idx in range(1, nbr+1):
        mvcells = mvcells + cell_mixed_volume(idx, vrblvl)
    if vrblvl > 0:
        print('sum of mixed volumes of cells :', mvcells)
    fail = fail + int(mvl != 3) + int(smv != 5) + int(mvcells != 5)
    return fail

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
