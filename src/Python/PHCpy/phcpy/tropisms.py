"""
The module tropisms exports functions to manage numerically computed
tropisms in double, double double, or quad double precision,
in a polyhedral end game with aposteriori step size control.
"""
from ctypes import c_int32, c_double, pointer
from random import randint, uniform
from phcpy.version import get_phcfun
from phcpy.trackers import set_parameter_value, show_parameters
from phcpy.trackers import double_track
from phcpy.trackers import double_double_track, quad_double_track
from phcpy.solver import solve

def double_initialize_tropisms(nbt, dim, wnd, dirs, errs, vrblvl=0):
    r"""
    Initializes the tropisms, given in double precision,
    along with estimates for their winding numbers and errors.
    On entry are the following five parameters:

    *nbt*: the number of direction vectors;

    *dim*: the number of coordinates in each vector;

    *wnd*: a list of integer values for the winding numbers, as many as *nbt*;

    *dirs*: a list of lists of doubles with the coordinates of the directions,
    each inner list has *dim* doubles and *nbt* vectors are given;

    *errs*: a list of *nbt* doubles.
    """
    if vrblvl > 0:
        print('in double_initialize_tropisms, nbt :', nbt, end='')
        print(', dim :', dim)
        print('the winding numbers :', wnd)
        print('the coordinates of tropisms :')
        for direction in dirs:
            print(direction)
        print('the errors :', errs)
    phc = get_phcfun(vrblvl-1)
    pars = (c_int32 * 2)()
    pars[0] = nbt
    pars[1] = dim
    apars = pointer(pars)
    bwnd = (c_int32 * nbt)()
    for (idx, nbr) in enumerate(wnd):
        bwnd[idx] = c_int32(nbr)
    pwnd = pointer(bwnd)
    size = nbt*(dim+1)
    cffs = (c_double * size)()
    idx = 0
    for direction in dirs:
        for nbr in direction:
            cffs[idx] = c_double(nbr)
            idx = idx + 1
    for (idx, nbr) in enumerate(errs):
        cffs[nbt*dim+idx] = c_double(nbr)
    pcffs = pointer(cffs)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_initialize_tropisms calls phc ...')
    retval = phc(711, apars, pwnd, pcffs, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
    return retval

def double_double_initialize_tropisms(nbt, dim, wnd, dirs, errs, vrblvl=0):
    r"""
    Initializes the tropisms, given in double double precision,
    along with estimates for their winding numbers and errors.
    On entry are the following five parameters:

    *nbt*: the number of direction vectors;

    *dim*: the number of coordinates in each vector;

    *wnd*: a list of integer values for the winding numbers, as many as *nbt*;

    *dirs*: a list of lists of doubles with the coordinates of the directions,
    each inner list has *dim* doubles and *nbt* vectors are given;

    *errs*: a list of *nbt* doubles.
    """
    if vrblvl > 0:
        print('in double_double_initialize_tropisms, nbt :', nbt, end='')
        print(', dim :', dim)
        print('the winding numbers :', wnd)
        print('the coordinates of tropisms :')
        for direction in dirs:
            print(direction)
        print('the errors :', errs)
    phc = get_phcfun(vrblvl-1)
    pars = (c_int32 * 2)()
    pars[0] = nbt
    pars[1] = dim
    apars = pointer(pars)
    bwnd = (c_int32 * nbt)()
    for (idx, nbr) in enumerate(wnd):
        bwnd[idx] = c_int32(nbr)
    pwnd = pointer(bwnd)
    size = 2*nbt*(dim+1)
    cffs = (c_double * size)()
    idx = 0
    for direction in dirs:
        for nbr in direction:
            cffs[idx] = c_double(nbr)
            idx = idx + 1
    for (idx, nbr) in enumerate(errs):
        cffs[2*nbt*dim+idx] = c_double(nbr)
    pcffs = pointer(cffs)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_initialize_tropisms calls phc ...')
    retval = phc(712, apars, pwnd, pcffs, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
    return retval

def quad_double_initialize_tropisms(nbt, dim, wnd, dirs, errs, vrblvl=0):
    r"""
    Initializes the tropisms, given in quad double precision,
    along with estimates for their winding numbers and errors.
    On entry are the following five parameters:

    *nbt*: the number of direction vectors;

    *dim*: the number of coordinates in each vector;

    *wnd*: a list of integer values for the winding numbers, as many as *nbt*;

    *dirs*: a list of lists of doubles with the coordinates of the directions,
    each inner list has *dim* doubles and *nbt* vectors are given;

    *errs*: a list of *nbt* doubles.
    """
    if vrblvl > 0:
        print('in quad_double_initialize_tropisms, nbt :', nbt, end='')
        print(', dim :', dim)
        print('the winding numbers :', wnd)
        print('the coordinates of tropisms :')
        for direction in dirs:
            print(direction)
        print('the errors :', errs)
    phc = get_phcfun(vrblvl-1)
    pars = (c_int32 * 2)()
    pars[0] = nbt
    pars[1] = dim
    apars = pointer(pars)
    bwnd = (c_int32 * nbt)()
    for (idx, nbr) in enumerate(wnd):
        bwnd[idx] = c_int32(nbr)
    pwnd = pointer(bwnd)
    size = 4*nbt*(dim+1)
    cffs = (c_double * size)()
    idx = 0
    for direction in dirs:
        for nbr in direction:
            cffs[idx] = c_double(nbr)
            idx = idx + 1
    for (idx, nbr) in enumerate(errs):
        cffs[4*nbt*dim+idx] = c_double(nbr)
    pcffs = pointer(cffs)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_initialize_tropisms calls phc ...')
    retval = phc(713, apars, pwnd, pcffs, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
    return retval

def double_tropisms_number(vrblvl=0):
    """
    Returns the number of tropisms in double precision.
    """
    if vrblvl > 0:
        print('in double_tropisms_number')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_tropisms_number calls phc ...')
    retval = phc(720, apar, bbb, ccc, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
        print('the number of tropisms :', apar[0])
    return apar[0]

def double_double_tropisms_number(vrblvl=0):
    """
    Returns the number of tropisms in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_tropisms_number')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_tropisms_number calls phc ...')
    retval = phc(721, apar, bbb, ccc, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
        print('the number of tropisms :', apar[0])
    return apar[0]

def quad_double_tropisms_number(vrblvl=0):
    """
    Returns the number of tropisms in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_tropisms_number')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_tropisms_number calls phc ...')
    retval = phc(722, apar, bbb, ccc, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
        print('the number of tropisms :', apar[0])
    return apar[0]

def double_tropisms_dimension(vrblvl=0):
    """
    Returns the dimension of tropisms in double precision.
    """
    if vrblvl > 0:
        print('in double_tropisms_dimension')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_tropisms_dimension calls phc ...')
    retval = phc(729, apar, bbb, ccc, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
        print('the dimension of tropisms :', apar[0])
    return apar[0]

def double_double_tropisms_dimension(vrblvl=0):
    """
    Returns the dimension of tropisms in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_tropisms_dimension')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_tropisms_dimension calls phc ...')
    retval = phc(730, apar, bbb, ccc, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
        print('the dimension of tropisms :', apar[0])
    return apar[0]

def quad_double_tropisms_dimension(vrblvl=0):
    """
    Returns the dimension of tropisms in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_tropisms_dimension')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_tropisms_dimension calls phc ...')
    retval = phc(731, apar, bbb, ccc, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
        print('the dimension of tropisms :', apar[0])
    return apar[0]

def get_double_tropisms(nbt, dim, vrblvl=0):
    r"""
    Given on input the number of tropisms in *nbt* and the dimension in *dim*,
    returns a tuple of three lists: the winding numbers, coordinates of
    the direction vectors, and the errors; in double precision.
    """
    if vrblvl > 0:
        print('in get_double_tropisms, nbt :', nbt, end='')
        print(', dim :', dim)
    phc = get_phcfun(vrblvl-1)
    pars = (c_int32 * 2)()
    pars[0] = nbt
    pars[1] = dim
    apars = pointer(pars)
    bwnd = (c_int32 * nbt)()
    for idx in range(nbt):
        bwnd[idx] = c_int32(0)
    pwnd = pointer(bwnd)
    size = nbt*(dim+1)
    cffs = (c_double * size)()
    for idx in range(size):
        cffs[idx] = c_double(0.0)
    pcffs = pointer(cffs)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_double_tropisms calls phc ...')
    retval = phc(717, apars, pwnd, pcffs, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
    vals = pwnd[0:nbt]
    wndnbrs = []
    for idx in range(nbt):
        wndnbrs.append(int(vals[0][idx]))
    vals = pcffs[0:size]
    directions = []
    idx = 0
    for _ in range(nbt):
        direction = []
        for _ in range(dim):
            direction.append(float(vals[0][idx]))
            idx = idx + 1
        directions.append(direction)
    errors = []
    for idx in range(size-dim+1, size):
        errors.append(float(vals[0][idx]))
    if vrblvl > 0:
        print('winding numbers :', wndnbrs)
        print('directions :')
        for direction in directions:
            print(direction)
        print('errors :', errors)
    return (wndnbrs, directions, errors)

def get_double_double_tropisms(nbt, dim, vrblvl=0):
    r"""
    Given on input the number of tropisms in *nbt* and the dimension in *dim*,
    returns a tuple of three lists: the winding numbers, coordinates of
    the direction vectors, and the errors; in double double precision.
    """
    if vrblvl > 0:
        print('in get_double_double_tropisms, nbt :', nbt, end='')
        print(', dim :', dim)
    phc = get_phcfun(vrblvl-1)
    pars = (c_int32 * 2)()
    pars[0] = nbt
    pars[1] = dim
    apars = pointer(pars)
    bwnd = (c_int32 * nbt)()
    for idx in range(nbt):
        bwnd[idx] = c_int32(0)
    pwnd = pointer(bwnd)
    size = 2*nbt*(dim+1)
    cffs = (c_double * size)()
    for idx in range(size):
        cffs[idx] = c_double(0.0)
    pcffs = pointer(cffs)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_double_double_tropisms calls phc ...')
    retval = phc(718, apars, pwnd, pcffs, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
    vals = pwnd[0:nbt]
    wndnbrs = []
    for idx in range(nbt):
        wndnbrs.append(int(vals[0][idx]))
    vals = pcffs[0:size]
    directions = []
    idx = 0
    for _ in range(nbt):
        direction = []
        for _ in range(2*dim):
            direction.append(float(vals[0][idx]))
            idx = idx + 1
        directions.append(direction)
    errors = []
    for idx in range(size-2*dim+2, size):
        errors.append(float(vals[0][idx]))
    if vrblvl > 0:
        print('winding numbers :', wndnbrs)
        print('directions :')
        for direction in directions:
            print(direction)
        print('errors :', errors)
    return (wndnbrs, directions, errors)

def get_quad_double_tropisms(nbt, dim, vrblvl=0):
    r"""
    Given on input the number of tropisms in *nbt* and the dimension in *dim*,
    returns a tuple of three lists: the winding numbers, coordinates of
    the direction vectors, and the errors; in quad double precision.
    """
    if vrblvl > 0:
        print('in get_quad_double_tropisms, nbt :', nbt, end='')
        print(', dim :', dim)
    phc = get_phcfun(vrblvl-1)
    pars = (c_int32 * 2)()
    pars[0] = nbt
    pars[1] = dim
    apars = pointer(pars)
    bwnd = (c_int32 * nbt)()
    for idx in range(nbt):
        bwnd[idx] = c_int32(0)
    pwnd = pointer(bwnd)
    size = 4*nbt*(dim+1)
    cffs = (c_double * size)()
    for idx in range(size):
        cffs[idx] = c_double(0.0)
    pcffs = pointer(cffs)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_quad_double_tropisms calls phc ...')
    retval = phc(719, apars, pwnd, pcffs, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
    vals = pwnd[0:nbt]
    wndnbrs = []
    for idx in range(nbt):
        wndnbrs.append(int(vals[0][idx]))
    vals = pcffs[0:size]
    directions = []
    idx = 0
    for _ in range(nbt):
        direction = []
        for _ in range(4*dim):
            direction.append(float(vals[0][idx]))
            idx = idx + 1
        directions.append(direction)
    errors = []
    for idx in range(size-4*dim+4, size):
        errors.append(float(vals[0][idx]))
    if vrblvl > 0:
        print('winding numbers :', wndnbrs)
        print('directions :')
        for direction in directions:
            print(direction)
        print('errors :', errors)
    return (wndnbrs, directions, errors)

def clear_double_tropisms(vrblvl=0):
    """
    Clears the tropisms in double precision.
    """
    if vrblvl > 0:
        print('in clear_double_tropisms')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_double_tropisms calls phc ...')
    retval = phc(726, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
    return retval

def clear_double_double_tropisms(vrblvl=0):
    """
    Clears the tropisms in double double precision.
    """
    if vrblvl > 0:
        print('in clear_double_double_tropisms')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_double_double_tropisms calls phc ...')
    retval = phc(727, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
    return retval

def clear_quad_double_tropisms(vrblvl=0):
    """
    Clears the tropisms in quad double precision.
    """
    if vrblvl > 0:
        print('in clear_quad_double_tropisms')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_quad_double_tropisms calls phc ...')
    retval = phc(728, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print('the return value of phc :', retval)
    return retval

def test_double_tropisms_data(vrblvl=0):
    """
    Tests tropisms data in double precision.
    """
    if vrblvl > 0:
        print('in test_double_tropisms_data ...')
    nbr = 2
    dim = 3
    wnd = [randint(10, 99) for _ in range(nbr)]
    dirs = []
    for _ in range(nbr):
        dirs.append([uniform(-1, +1) for _ in range(dim)])
    errs = [uniform(0,1)/1000.0 for _ in range(nbr)]
    if vrblvl > 0:
        print('random winding numbers :', wnd)
        for k in range(len(dirs)):
            print('direction', k+1, ':', dirs[k])
        print('random errors :', errs)
    double_initialize_tropisms(nbr, dim, wnd, dirs, errs, vrblvl)
    size = double_tropisms_number(vrblvl)
    if vrblvl > 0:
        print('number of tropisms :', size)
    fail = int(size != nbr)
    tdim = double_tropisms_dimension(vrblvl)
    if vrblvl > 0:
        print('dimension of tropisms :', tdim)
    fail = fail + int(tdim != dim)
    (rwnd, rdirs, rerrs) = get_double_tropisms(nbr, dim, vrblvl)
    if vrblvl > 0:
        print('retrieved winding numbers :', rwnd)
    for (awnd, bwnd) in zip(wnd, rwnd):
        fail = fail + int(awnd != bwnd)
    if vrblvl > 0:
        print('retrieved directions :')
        for (idx, direction) in enumerate(rdirs):
            print('direction', idx+1, ':', direction)
    for (adir, bdir) in zip(dirs, rdirs):
        for (anbr, bnbr) in zip(adir, bdir):
            fail = fail + int(anbr != bnbr)
    if vrblvl > 0:
        print('retrieved errors :', rerrs)
    for (aerr, berr) in zip(errs, rerrs):
        fail = fail + int(aerr != berr)
    clear_double_tropisms(vrblvl)
    return fail

def test_double_double_tropisms_data(vrblvl=0):
    """
    Tests tropisms data in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_tropisms_data ...')
    nbr = 2
    dim = 3
    wnd = [randint(10, 99) for _ in range(nbr)]
    dirs = []
    for _ in range(nbr):
        dirs.append([uniform(-1, +1) for _ in range(2*dim)])
    errs = [uniform(0,1)/1000.0 for _ in range(2*nbr)]
    if vrblvl > 0:
        print('random winding numbers :', wnd)
        for k in range(len(dirs)):
            print('direction', k+1, ':', dirs[k])
        print('random errors :', errs)
    double_double_initialize_tropisms(nbr, dim, wnd, dirs, errs, vrblvl)
    size = double_double_tropisms_number(vrblvl)
    if vrblvl > 0:
        print('number of tropisms :', size)
    fail = int(size != nbr)
    tdim = double_double_tropisms_dimension(vrblvl)
    if vrblvl > 0:
        print('dimension of tropisms :', tdim)
    fail = fail + int(tdim != dim)
    (rwnd, rdirs, rerrs) = get_double_double_tropisms(nbr, dim, vrblvl)
    if vrblvl > 0:
        print('retrieved winding numbers :', rwnd)
    for (awnd, bwnd) in zip(wnd, rwnd):
        fail = fail + int(awnd != bwnd)
    if vrblvl > 0:
        print('retrieved directions :')
        for (idx, direction) in enumerate(rdirs):
            print('direction', idx+1, ':', direction)
    for (adir, bdir) in zip(dirs, rdirs):
        for (anbr, bnbr) in zip(adir, bdir):
            fail = fail + int(anbr != bnbr)
    if vrblvl > 0:
        print('retrieved errors :', rerrs)
    for (aerr, berr) in zip(errs, rerrs):
        fail = fail + int(aerr != berr)
    clear_double_double_tropisms(vrblvl)
    return fail

def test_quad_double_tropisms_data(vrblvl=0):
    """
    Tests tropisms data in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_tropisms_data ...')
    nbr = 2
    dim = 3
    wnd = [randint(10, 99) for _ in range(nbr)]
    dirs = []
    for _ in range(nbr):
        dirs.append([uniform(-1, +1) for _ in range(4*dim)])
    errs = [uniform(0,1)/1000.0 for _ in range(4*nbr)]
    if vrblvl > 0:
        print('random winding numbers :', wnd)
        for k in range(len(dirs)):
            print('direction', k+1, ':', dirs[k])
        print('random errors :', errs)
    quad_double_initialize_tropisms(nbr, dim, wnd, dirs, errs, vrblvl)
    size = quad_double_tropisms_number(vrblvl)
    if vrblvl > 0:
        print('number of tropisms :', size)
    fail = int(size != nbr)
    tdim = quad_double_tropisms_dimension(vrblvl)
    if vrblvl > 0:
        print('dimension of tropisms :', tdim)
    fail = fail + int(tdim != dim)
    (rwnd, rdirs, rerrs) = get_quad_double_tropisms(nbr, dim, vrblvl)
    if vrblvl > 0:
        print('retrieved winding numbers :', rwnd)
    for (awnd, bwnd) in zip(wnd, rwnd):
        fail = fail + int(awnd != bwnd)
    if vrblvl > 0:
        print('retrieved directions :')
        for (idx, direction) in enumerate(rdirs):
            print('direction', idx+1, ':', direction)
    for (adir, bdir) in zip(dirs, rdirs):
        for (anbr, bnbr) in zip(adir, bdir):
            fail = fail + int(anbr != bnbr)
    if vrblvl > 0:
        print('retrieved errors :', rerrs)
    for (aerr, berr) in zip(errs, rerrs):
        fail = fail + int(aerr != berr)
    clear_quad_double_tropisms(vrblvl)
    return fail

def test_double_endgame(vrblvl=0):
    """
    Tests the numerical computation of a tropism,
    in double precision.
    """
    if vrblvl > 0:
        print('in test_double_endgame ...')
    pols = ['x*z-y*z-z^2+x;', 'z^3-x*y-y*z-z^2-x;', 'x;']
    start = [ 'x*z-y*z-z^2+x;', 'z^3-x*y-y*z-z^2-x;', \
              '(-8.16462369561530E-01 - 5.77398648327108E-01*i)*x' \
           +  '+( 6.63704141654579E-01 - 7.47995195406065E-01*i);']
    startsols = solve(start)
    if vrblvl > 0:
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    if vrblvl > 0:
        print('setting the order of the extrapolator ...')
    set_parameter_value(5, 4, vrblvl-1)
    if vrblvl > 0:
        print('settings of the parameter :')
        show_parameters(vrblvl-1)
    gmm = complex(-0.9669413930172692, 0.25499871072188346)
    gamma, sols = double_track(pols, start, startsols, \
        gamma=gmm, pwt=1, vrblvl=vrblvl-1)
    if vrblvl > 0:
        print('gamma :', gamma)
        print('the solutions at the end :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    size = double_tropisms_number(vrblvl-1)
    if vrblvl > 0:
        print('the number of tropisms :', size)
    (wnd, dirs, errs) = get_double_tropisms(len(sols), len(pols), vrblvl)
    if vrblvl > 0:
        print('the winding numbers :', wnd)
        print('the directions :')
        for direction in dirs:
            print(direction)
        print('the errors :', errs)
    return int(sum(wnd) != 10)

def test_double_double_endgame(vrblvl=0):
    """
    Tests the numerical computation of a tropism,
    in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_endgame ...')
    pols = ['x*z-y*z-z^2+x;', 'z^3-x*y-y*z-z^2-x;', 'x;']
    start = [ 'x*z-y*z-z^2+x;', 'z^3-x*y-y*z-z^2-x;', \
              '(-8.16462369561530E-01 - 5.77398648327108E-01*i)*x' \
           +  '+( 6.63704141654579E-01 - 7.47995195406065E-01*i);']
    startsols = solve(start)
    if vrblvl > 0:
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    if vrblvl > 0:
        print('setting the order of the extrapolator ...')
    set_parameter_value(5, 4, vrblvl-1)
    if vrblvl > 0:
        print('settings of the parameter :')
        show_parameters(vrblvl-1)
    gmm = complex(-0.9669413930172692, 0.25499871072188346)
    gamma, sols = double_double_track(pols, start, startsols, \
        gamma=gmm, pwt=1, vrblvl=vrblvl-1)
    if vrblvl > 0:
        print('gamma :', gamma)
        print('the solutions at the end :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    size = double_double_tropisms_number(vrblvl-1)
    if vrblvl > 0:
        print('the number of tropisms :', size)
    (wnd, dirs, errs) = get_double_double_tropisms(len(sols), len(pols), \
        vrblvl)
    if vrblvl > 0:
        print('the winding numbers :', wnd)
        print('the directions :')
        for direction in dirs:
            print(direction)
        print('the errors :', errs)
    return int(sum(wnd) != 10)

def test_quad_double_endgame(vrblvl=0):
    """
    Tests the numerical computation of a tropism,
    in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_endgame ...')
    pols = ['x*z-y*z-z^2+x;', 'z^3-x*y-y*z-z^2-x;', 'x;']
    start = [ 'x*z-y*z-z^2+x;', 'z^3-x*y-y*z-z^2-x;', \
              '(-8.16462369561530E-01 - 5.77398648327108E-01*i)*x' \
           +  '+( 6.63704141654579E-01 - 7.47995195406065E-01*i);']
    startsols = solve(start)
    if vrblvl > 0:
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    if vrblvl > 0:
        print('setting the order of the extrapolator ...')
    set_parameter_value(5, 4, vrblvl-1)
    if vrblvl > 0:
        print('settings of the parameter :')
        show_parameters(vrblvl-1)
    gmm = complex(0.1401011540077633, -0.990137195870195)
    gamma, sols = quad_double_track(pols, start, startsols, \
        gamma=gmm, pwt=1, vrblvl=vrblvl-1)
    print(gamma)
    if vrblvl > 0:
        print('the solutions at the end :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    size = quad_double_tropisms_number(vrblvl-1)
    if vrblvl > 0:
        print('the number of tropisms :', size)
    (wnd, dirs, errs) = get_quad_double_tropisms(len(sols), len(pols), \
        vrblvl)
    if vrblvl > 0:
        print('the winding numbers :', wnd)
        print('the directions :')
        for direction in dirs:
            print(direction)
        print('the errors :', errs)
    return int(sum(wnd) != 10)

def main():
    """
    Runs some tests.
    """
    lvl = 1
    fail = test_double_tropisms_data(lvl)
    fail = fail + test_double_double_tropisms_data(lvl)
    fail = fail + test_quad_double_tropisms_data(lvl)
    fail = fail + test_double_endgame(lvl)
    fail = fail + test_double_double_endgame(lvl)
    fail = fail + test_quad_double_endgame(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__== '__main__':
    main()
