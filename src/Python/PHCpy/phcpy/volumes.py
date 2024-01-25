"""
Exports functions to compute mixed volumes and stable mixed volumes.
"""
from ctypes import c_int32, c_double, pointer
from phcpy.version import get_phcfun
from phcpy.polynomials import set_double_system, get_double_system
from phcpy.solutions import get_double_solutions, verify
from phcpy.solutions import clear_double_solutions

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

def are_cells_stable(vrblvl=0):
    """
    Returns True if stable mixed cells were computed,
    returns False otherwise.
    """
    if vrblvl > 0:
        print('in are_cells_stable ...')
    phc = get_phcfun(vrblvl-1)
    aflag = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> are_cells_stable calls phc', end='')
    retval = phc(879, aflag, bbb, ccc, vrb)
    result = aflag[0] == 1
    if vrblvl > 0:
        print(', return value :', retval)
        print('the flag :', result)
    return result

def number_of_cells(vrblvl=0):
    """
    Returns the number of cells computed by the mixed_volume function.
    """
    if vrblvl > 0:
        print('in number_of_cells ...')
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

def number_of_original_cells(vrblvl=0):
    """
    Returns the number of original cells,
    the cells without artificial origin,
    computed by the mixed_volume function.
    """
    if vrblvl > 0:
        print('in number_of_mixed_original_cells ...')
    phc = get_phcfun(vrblvl-1)
    alen = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> number_of_original_cells calls phc', end='')
    retval = phc(880, alen, bbb, ccc, vrb)
    result = alen[0]
    if vrblvl > 0:
        print(', return value :', retval)
        print('the number of original cells :', result)
    return result

def number_of_stable_cells(vrblvl=0):
    """
    Returns the number of stable cells computed 
    by the mixed_volume function.
    """
    if vrblvl > 0:
        print('in number_of_mixed_cells ...')
    phc = get_phcfun(vrblvl-1)
    alen = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> number_of_cells calls phc', end='')
    retval = phc(881, alen, bbb, ccc, vrb)
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

def double_random_coefficient_system(vrblvl=0):
    """
    Makes a random coefficient system using the type of mixture
    and the supports used to compute the mixed volume.
    Sets the system in double precision.
    """
    if vrblvl > 0:
        print('in double_random_coefficient_system ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_random_coefficient_system calls phc', end='')
    retval = phc(96, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_coefficient_system(vrblvl=0):
    """
    Sets the constructed random coefficient system to the system
    in double precision.
    """
    if vrblvl > 0:
        print('in set_double_coefficient_system ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_coefficient_system calls phc', end='')
    retval = phc(99, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_polyhedral_homotopy(vrblvl=0):
    """
    Sets the polyhedral homotopy in double precision based
    on the lifted supports used to compute the mixed volume.
    Initializes the data to store the solutions.
    """
    if vrblvl > 0:
        print('in set_double_polyhedral_homotopy ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_polyhedral_homotopy calls phc', end='')
    retval = phc(101, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_solve_start_system(idx, vrblvl=0):
    """
    Solves the start system corresponding to the cell with index idx,
    using double precision arithmetic.  Return the number of solutions
    found, which must equal the mixed volume of the cell.
    """
    if vrblvl > 0:
        print('in double_solve_start_system, idx :', idx)
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(idx))
    bnbr = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_solve_start_system calls phc', end='')
    retval = phc(102, aidx, bnbr, ccc, vrb)
    result = bnbr[0]
    if vrblvl > 0:
        print(', return value :', retval)
        print('number of solutions found :', result)
    return result

def double_track_path(cellidx, pathidx, vrblvl=0):
    """
    Tracks path with index pathidx defined by cell with index cellidx.
    Prior to calling this function, the start system corresponding
    to the cell must have been solved.
    """
    if vrblvl > 0:
        print('in double_track_path, cell :', cellidx, end='')
        print(', path :', pathidx)
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(cellidx))
    nbrs = (c_int32 * 2)()
    nbrs[0] = pathidx
    nbrs[1] = vrblvl
    bnbrs = pointer(nbrs)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_track_path calls phc', end='')
    retval = phc(103, aidx, bnbrs, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_start_solution(cellidx, solidx, vrblvl=0):
    """
    Sets the start solution of index solidx, with cell with index cellidx,
    to the solutions in double precision.
    """ 
    if vrblvl > 0:
        print('in set_double_start_solution, cell :', cellidx, end='')
        print(', solidx :', solidx)
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(cellidx))
    bidx = pointer(c_int32(solidx))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_start_solution calls phc', end='')
    retval = phc(104, aidx, bidx, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_polyhedral_homotopies(vrblvl=0):
    r"""
    Runs the polyhedral homotopies and returns a random coefficient
    system based on the contents of the cell container,
    in standard double precision arithmetic.
    For this to work, the function **mixed_volume()** must be called first.
    """
    if vrblvl > 0:
        print('in double_random_coefficient_system ...')
    double_random_coefficient_system(vrblvl)
    set_double_coefficient_system(vrblvl)
    rndcffsys = get_double_system(vrblvl-1)
    set_double_polyhedral_homotopy(vrblvl)
    is_stable = are_cells_stable(vrblvl)
    if is_stable:
        nbcells = number_of_original_cells(vrblvl)
    else:
        nbcells = number_of_cells(vrblvl)
    if vrblvl > 0:
         print('the number of cells :', nbcells)
    clear_double_solutions(vrblvl-1)
    for cell in range(1, nbcells+1):
        nbsols = double_solve_start_system(cell, vrblvl)
        if vrblvl > 0:
            print(f'system {cell} has {nbsols} solutions')
        for idxsol in range(1, nbsols+1):
            if vrblvl > 0:
                print(f'-> tracking path {idxsol} out of {nbsols}')
            double_track_path(cell, idxsol, vrblvl)
            set_double_start_solution(cell, idxsol, vrblvl)
    sols = get_double_solutions(vrblvl-1)
    return (rndcffsys, sols)

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
    stbnbr = number_of_stable_cells(vrblvl)
    if vrblvl > 0:
        print('the mixed volume by DEMiCs :', mvl)
        print('the stable mixed volume by DEMiCs :', smv)
        print('the number of cells :', nbr)
        print('the number of stable cells :', stbnbr)
    mvcells = 0
    for idx in range(1, nbr+1):
        mvcells = mvcells + cell_mixed_volume(idx, vrblvl)
    if vrblvl > 0:
        print('sum of mixed volumes of cells :', mvcells)
    fail = int(mvl != 3) + int(smv != 5) + int(mvcells != 5)
    clear_cells(vrblvl)
    mvl, smv = stable_mixed_volume(False, vrblvl)
    nbr = number_of_cells(vrblvl)
    stbnbr = number_of_stable_cells(vrblvl)
    if vrblvl > 0:
        print('the mixed volume by MixedVol :', mvl)
        print('the stable mixed volume by MixedVol :', smv)
        print('the number of cells :', nbr)
        print('the number of stable cells :', stbnbr)
    mvcells = 0
    for idx in range(1, nbr+1):
        mvcells = mvcells + cell_mixed_volume(idx, vrblvl)
    if vrblvl > 0:
        print('sum of mixed volumes of cells :', mvcells)
    fail = fail + int(mvl != 3) + int(smv != 5) + int(mvcells != 5)
    return fail

def test_double_polyhedral_homotopies(vrblvl=0):
    """
    Test the polyhedral homotopies in double precision
    to compute and solve a random coefficient system.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_double_polyhedral_homotopies ...')
    polynomials = ["x^3 + 2*x*y - 1;", "x + y - 1;"]
    mv = mixed_volume(True, vrblvl)
    rndcffsys, sols = double_polyhedral_homotopies(vrblvl)
    fail = int(mv != len(sols))
    if vrblvl > 0:
        print('a random coefficient system :')
        for pol in rndcffsys:
            print(pol)
        print('its solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(rndcffsys, sols, vrblvl-1)
    if vrblvl > 0:
        print('the sum of errors :', err)
    fail = fail + int(abs(err.real) > 1.0e-8) + int(abs(err.imag) > 1.0e-8)
    return fail

def main():
    """
    Runs tests on mixed volumes and stable mixed volumes.
    """
    lvl = 1
    fail = test_mixed_volume(lvl)
    fail = fail + test_stable_mixed_volume(lvl)
    fail = fail + test_double_polyhedral_homotopies(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()
