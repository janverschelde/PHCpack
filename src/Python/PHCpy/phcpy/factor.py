"""
Given a witness set representation of a pure dimensional solution set,
the functions in this module separate the generic points in the witness set
according to the irreducible components of the solution set.
"""
from ctypes import c_int32, c_double, pointer
from phcpy.version import get_phcfun, int4a2nbr
from phcpy.solutions import make_solution, endmultiplicity
from phcpy.solutions import get_double_solutions
from phcpy.solutions import set_double_solutions
from phcpy.solutions import clear_double_solutions
from phcpy.solutions import get_double_double_solutions
from phcpy.solutions import set_double_double_solutions
from phcpy.solutions import clear_double_double_solutions
from phcpy.solutions import get_quad_double_solutions
from phcpy.solutions import set_quad_double_solutions
from phcpy.solutions import clear_quad_double_solutions
from phcpy.sets import set_double_witness_set
from phcpy.sets import set_double_double_witness_set
from phcpy.sets import set_quad_double_witness_set
from phcpy.sets import set_double_Laurent_witness_set
from phcpy.sets import set_double_double_Laurent_witness_set
from phcpy.sets import set_quad_double_Laurent_witness_set

def set_double_verbose(verbose=True, vrblvl=0):
    """
    Sets the state of the monodromy algorithm in double precision
    to verbose if vrblvl > 0, otherwise it is set to remain mute.
    """
    if vrblvl > 0:
        print('in set_double_verbose, verbose :', verbose)
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_monodromy_verbose calls phc', end='')
    if verbose:
        retval = phc(630, aaa, bbb, ccc, vrb)
    else:
        retval = phc(39, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_double_verbose(verbose=True, vrblvl=0):
    """
    Sets the state of the monodromy algorithm in double double precision
    to verbose if vrblvl > 0, otherwise it is set to remain mute.
    """
    if vrblvl > 0:
        print('in set_double_double_verbose, verbose :', verbose)
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_double_verbose calls phc', end='')
    if verbose:
        retval = phc(660, aaa, bbb, ccc, vrb)
    else:
        retval = phc(658, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_quad_double_monodromy_verbose(verbose=True, vrblvl=0):
    """
    Sets the state of the monodromy algorithm in quad double precision
    to verbose if vrblvl > 0, otherwise it is set to remain mute.
    """
    if vrblvl > 0:
        print('in set_quad_double_verbose, verbose :', verbose)
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_quad_double_verbose calls phc', end='')
    if verbose:
        retval = phc(690, aaa, bbb, ccc, vrb)
    else:
        retval = phc(688, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_assign_labels(nvr, vrblvl=0):
    """
    Assigns a unique label to the multiplicity field
    for each solution set in double precision.
    """
    if vrblvl > 0:
        print('in double_assign_labels, nvr :', nvr)
    sols = get_double_solutions(vrblvl-1)
    newsols = []
    for (idx, sol) in enumerate(sols):
        data = sol.split('m : 1')
        idxm = f"m : {idx+1}"
        newsol = idxm.join(data)
        newsols.append(newsol)
    clear_double_solutions(vrblvl-1)
    set_double_solutions(nvr, newsols, vrblvl-1)
    return 0

def double_double_assign_labels(nvr, vrblvl=0):
    """
    Assigns a unique label to the multiplicity field
    for each solution set in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_assign_labels, nvr :', nvr)
    sols = get_double_double_solutions(vrblvl-1)
    newsols = []
    for (idx, sol) in enumerate(sols):
        data = sol.split('m : 1')
        idxm = f"m : {idx+1}"
        newsol = idxm.join(data)
        newsols.append(newsol)
    clear_double_double_solutions(vrblvl-1)
    set_double_double_solutions(nvr, newsols, vrblvl-1)
    return 0

def quad_double_assign_labels(nvr, vrblvl=0):
    """
    Assigns a unique label to the multiplicity field
    for each solution set in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_assign_labels, nvr :', nvr)
    sols = get_quad_double_solutions(vrblvl-1)
    newsols = []
    for (idx, sol) in enumerate(sols):
        data = sol.split('m : 1')
        idxm = f"m : {idx+1}"
        newsol = idxm.join(data)
        newsols.append(newsol)
    clear_quad_double_solutions(vrblvl-1)
    set_quad_double_solutions(nvr, newsols, vrblvl-1)
    return 0

def initialize_double_sampler(dim, vrblvl=0):
    """
    Initializes the sampling machine with a witness set,
    set in double precision.
    """
    if vrblvl > 0:
        print('in initialize_double_sampler, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_double_sampler calls phc', end='')
    retval = phc(42, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_double_double_sampler(dim, vrblvl=0):
    """
    Initializes the sampling machine with a witness set,
    set in double double precision.
    """
    if vrblvl > 0:
        print('in initialize_double_double_sampler, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_double_double_sampler calls phc', end='')
    retval = phc(632, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_quad_double_sampler(dim, vrblvl=0):
    """
    Initializes the sampling machine with a witness set,
    set in quad double precision.
    """
    if vrblvl > 0:
        print('in initialize_quad_double_sampler, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_quad_double_sampler calls phc', end='')
    retval = phc(662, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_double_Laurent_sampler(dim, vrblvl=0):
    """
    Initializes the sampling machine with a witness set,
    for a Laurent system set in double precision.
    """
    if vrblvl > 0:
        print('in initialize_double_Laurent_sampler, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_double_Laurent_sampler calls phc', end='')
    retval = phc(804, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_double_double_Laurent_sampler(dim, vrblvl=0):
    """
    Initializes the sampling machine with a witness set,
    for a Laurent system set in double double precision.
    """
    if vrblvl > 0:
        print('in initialize_double_double_Laurent_sampler, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_double_double_Laurent_sampler calls phc', end='')
    retval = phc(805, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_quad_double_Laurent_sampler(dim, vrblvl=0):
    """
    Initializes the sampling machine with a witness set,
    for a Laurent system set in quad double precision.
    """
    if vrblvl > 0:
        print('in initialize_quad_double_Laurent_sampler, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_quad_double_Laurent_sampler calls phc', end='')
    retval = phc(806, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_double_monodromy(nbloops, deg, dim, vrblvl=0):
    """
    Initializes to run nbloops monodromy loops to factor
    a positive dimension solution set of dimension dim and of degree deg
    in double precision.
    """
    if vrblvl > 0:
        print('in initialize_double_monodromy', end='')
        print(', nbloops :', nbloops, ', deg :', deg, ', dim :', dim)
    phc = get_phcfun(vrblvl-1)
    anbloops = pointer(c_int32(nbloops))
    bpars = int4a2nbr([deg, dim], vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_double_monodromy calls phc', end='')
    retval = phc(50, anbloops, bpars, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_double_double_monodromy(nbloops, deg, dim, vrblvl=0):
    """
    Initializes to run nbloops monodromy loops to factor
    a positive dimension solution set of dimension dim and of degree deg
    in double double precision.
    """
    if vrblvl > 0:
        print('in initialize_double_double_monodromy', end='')
        print(', nbloops :', nbloops, ', deg :', deg, ', dim :', dim)
    phc = get_phcfun(vrblvl-1)
    anbloops = pointer(c_int32(nbloops))
    bpars = int4a2nbr([deg, dim], vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_double_double_monodromy calls phc', end='')
    retval = phc(640, anbloops, bpars, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_quad_double_monodromy(nbloops, deg, dim, vrblvl=0):
    """
    Initializes to run nbloops monodromy loops to factor
    a positive dimension solution set of dimension dim and of degree deg
    in double double precision.
    """
    if vrblvl > 0:
        print('in initialize_quad_double_monodromy', end='')
        print(', nbloops :', nbloops, ', deg :', deg, ', dim :', dim)
    phc = get_phcfun(vrblvl-1)
    anbloops = pointer(c_int32(nbloops))
    bpars = int4a2nbr([deg, dim], vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_double_double_monodromy calls phc', end='')
    retval = phc(670, anbloops, bpars, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_monodromy_breakup(embsys, esols, dim, \
    islaurent=False, verbose=True, nbloops=20, vrblvl=0):
    r"""
    Applies the monodromy breakup algorithm in standard double precision
    to factor the *dim*-dimensional algebraic set represented by the
    embedded system *embsys* and its solutions *esols*.
    If the embedded polynomial system is a Laurent system,
    then islaurent must equal True.
    If *verbose* is False, then no output is written.
    The value of *nbloops* equals the maximum number of loops.
    """
    if vrblvl > 0:
        print('in double_monodromy_breakup', end='')
        print(', dim :', dim, ', nbloops :', nbloops)
        print(', islaurent :', islaurent, ', verbose :', verbose)
        print('the embedded polynomials :')
        for pol in embsys:
            print(pol)
        print('the generic points :')
        for (idx, sol) in enumerate(esols):
            print('Solution', idx+1, ':')
            print(sol)
    deg = len(esols)
    nvr = len(embsys)
    set_double_verbose(verbose, vrblvl)
    if islaurent:
        set_double_Laurent_witness_set(nvr, dim, embsys, esols, vrblvl-1)
        double_assign_labels(nvr, vrblvl)
        initialize_double_Laurent_sampler(dim, vrblvl)
    else:
        set_double_witness_set(nvr, dim, embsys, esols, vrblvl-1)
        double_assign_labels(nvr, vrblvl)
        initialize_double_sampler(dim, vrblvl)

def test_double_assign_labels(vrblvl=0):
    """
    Tests the assigning of labels for solutions
    set in double precision.
    """
    if vrblvl > 0:
        print('in test_double_assign_labels ...')
    names = ['x', 'y']
    sols = []
    for idx in range(3):
        sol = make_solution(names, [idx, 0.0])
        sols.append(sol)
    if vrblvl > 0:
        print('the test solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    set_double_solutions(2, sols, vrblvl-1)
    double_assign_labels(2, vrblvl)
    newsols = get_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('after assigning labels :')
    fail = 0
    for (idx, sol) in enumerate(newsols):
        if vrblvl > 0:
            print('Solution', idx+1, ':')
            print(sol)
        (_, mval) = endmultiplicity(sol, vrblvl-1)
        if vrblvl > 0:
            print('mval :', mval, mval == idx+1)
        fail = fail + int(mval != idx+1)
    return fail

def test_double_double_assign_labels(vrblvl=0):
    """
    Tests the assigning of labels for solutions
    set in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_assign_labels ...')
    names = ['x', 'y']
    sols = []
    for idx in range(3):
        sol = make_solution(names, [idx, 0.0])
        sols.append(sol)
    if vrblvl > 0:
        print('the test solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    set_double_double_solutions(2, sols, vrblvl-1)
    double_double_assign_labels(2, vrblvl)
    newsols = get_double_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('after assigning labels :')
    fail = 0
    for (idx, sol) in enumerate(newsols):
        if vrblvl > 0:
            print('Solution', idx+1, ':')
            print(sol)
        (_, mval) = endmultiplicity(sol, vrblvl-1)
        if vrblvl > 0:
            print('mval :', mval, mval == idx+1)
        fail = fail + int(mval != idx+1)
    return fail

def test_quad_double_assign_labels(vrblvl=0):
    """
    Tests the assigning of labels for solutions
    set in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_assign_labels ...')
    names = ['x', 'y']
    sols = []
    for idx in range(3):
        sol = make_solution(names, [idx, 0.0])
        sols.append(sol)
    if vrblvl > 0:
        print('the test solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    set_quad_double_solutions(2, sols, vrblvl-1)
    quad_double_assign_labels(2, vrblvl)
    newsols = get_quad_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('after assigning labels :')
    fail = 0
    for (idx, sol) in enumerate(newsols):
        if vrblvl > 0:
            print('Solution', idx+1, ':')
            print(sol)
        (_, mval) = endmultiplicity(sol, vrblvl-1)
        if vrblvl > 0:
            print('mval :', mval, mval == idx+1)
        fail = fail + int(mval != idx+1)
    return fail

def main():
    """
    Runs some tests.
    """
    lvl = 1
    fail = test_double_assign_labels(lvl)
    fail = fail + test_double_double_assign_labels(lvl)
    fail = fail + test_quad_double_assign_labels(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()
