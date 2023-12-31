"""
The module trackers offers functions to track paths starting
at the known solutions of a start system and leading to
the desired solutions of a target system,
with aposteriori step control algorithms.
The tuning of the parameters for the predictor, corrector,
and the settings of the tolerances is handled by the tuning module.
"""
from ctypes import c_int32, c_double, pointer
from random import uniform
from cmath import exp, pi
from phcpy.version import get_phcfun
from phcpy.polynomials import set_double_system
from phcpy.solutions import set_double_solutions, get_double_solutions
from phcpy.solutions import clear_double_solutions, verify
from phcpy.homotopies import total_degree_start_system

def set_double_target_system(pols, vrblvl=0):
    """
    Sets the target system in an artificial parameter homotopy
    in double precision to the list of polynomials in pols,
    which is assumed to be square.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_double_target_system, with pols :')
        for pol in pols:
            print(pol)
    set_double_system(len(pols), pols, vrblvl)
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_target_system calls phc', end='')
    retval = phc(2, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_start_system(pols, vrblvl=0):
    """
    Sets the start system in an artificial parameter homotopy
    in double precision to the list of polynomials in pols,
    which is assumed to be square.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_double_start_system, with pols :')
        for pol in pols:
            print(pol)
    set_double_system(len(pols), pols, vrblvl)
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_start_system calls phc', end='')
    retval = phc(4, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_start_solutions(nvr, sols, vrblvl=0):
    """
    Sets the start solutions in an artificial parameter homotopy
    in double precision to the list of solutions in sols,
    where the number of variables in nvr must match the
    dimension of the start system.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_double_start_solutions, with nvr :', nvr)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx, ':')
            print(sol)
    set_double_solutions(nvr, sols, vrblvl)
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_start_solutions calls phc', end='')
    retval = phc(8, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def get_double_target_solutions(vrblvl=0):
    """
    Returns the list of target solutions computed in double precision.
    """
    if vrblvl > 0:
        print('in get_double_target_solutions ...')
    clear_double_solutions(vrblvl)
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_target_solutions calls phc', end='')
    retval = phc(5, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    sols = get_double_solutions(vrblvl)
    if vrblvl > 0:
        print('the target solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx, ':')
            print(sol)
    return sols

def set_double_homotopy(gamma=0, pwt=2, vrblvl=0):
    """
    After the target and start system are set in double precision,
    the homotopy is constructed with either a random gamma constant,
    or with the given complex value of gamma.
    The power of the continuation parameter is given by pwt.
    The gamma used to make the homotopy is returned.
    """
    if vrblvl > 0:
        print('in set_double_homotopy, with power :', pwt)
        print('gamma :', gamma)
    phc = get_phcfun()
    apwt = pointer(c_int32(pwt))
    bbb = pointer(c_int32(0))
    if gamma != 0:
        usegamma = gamma
    else:
        angle = uniform(0, 2*pi)
        usegamma = exp(angle*complex(0, 1))
        if(vrblvl > 0):
            print('random gamma :', usegamma)
    c_gamma = (c_double*2)()
    c_gamma[0] = c_double(usegamma.real)
    c_gamma[1] = c_double(usegamma.imag)
    ptr_gamma = pointer(c_gamma)
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_homotopy calls phc', end='')
    retval = phc(153, apwt, bbb, ptr_gamma, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return usegamma 

def clear_double_homotopy(vrblvl=0):
    """
    Clears the homotopy set in double precision.
    """
    if vrblvl > 0:
        print('in clear_double_homotopy ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> clear_double_homotopy calls phc', end='')
    retval = phc(154, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_double_track_data(vrblvl=0):
    """
    Clears the data allocated for path tracking double precision
    with an artificial parameter homotopy.
    """
    if vrblvl > 0:
         print('in clear_double_track_data ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> clear_double_track_data calls phc', end='')
    retval = phc(18, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_track(target, start, startsols, gamma=0, pwt=2, tasks=0, vrblvl=0):
    r"""
    Does path tracking in double precision.
    On input are a target system, a start system with solutions,
    optionally: a (random) gamma constant and the number of tasks.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in sols.
    The *startsols* is a list of strings representing start solutions.
    By default, a random *gamma* constant is generated,
    otherwise *gamma* must be a nonzero complex constant.
    The *pwt* is the power of t in the homotopy.
    Changing the default of *pwt* can only be done if a nonzero complex
    value for *gamma* is provided as well.
    The number of tasks in the multithreading is defined by *tasks*.
    The default zero value for *tasks* indicates no multithreading.
    On return is a tuple, with first the gamma used in the homotopy
    and then second, the string representations of the solutions
    computed at the end of the paths.
    """
    if vrblvl > 0:
        print('in double_track, with target :')
        for pol in target:
            print(pol)
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    set_double_target_system(target, vrblvl)
    set_double_start_system(start, vrblvl)
    set_double_start_solutions(len(target), startsols, vrblvl)
    usedgamma = set_double_homotopy(gamma, pwt, vrblvl)
    phc = get_phcfun()
    nbtasks = pointer(c_int32(tasks))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_track calls phc', end='')
    retval = phc(16, nbtasks, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    sols = get_double_target_solutions(vrblvl)
    clear_double_homotopy(vrblvl)
    clear_double_track_data(vrblvl)
    return (usedgamma, sols)

def test_double_track(vrblvl=0):
    """
    Runs on the mickey mouse example of two quadrics.
    """
    mickey = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    start, startsols = total_degree_start_system(mickey, vrblvl=vrblvl)
    print('the start system :')
    for pol in start:
        print(pol)
    print('the start solutions :')
    for (idx, sol) in enumerate(startsols):
        print('Solution', idx+1, ':')
        print(sol)
    gamma, sols = double_track(mickey, start, startsols, vrblvl=vrblvl)
    print('the solutions :')
    for (idx, sol) in enumerate(sols):
        print('Solution', idx+1, ':')
        print(sol)
    err = verify(mickey, sols, vrblvl)
    if vrblvl > 0:
        print('the error sum :', err)
    if len(sols) == 4 and abs(err.real + err.imag) < 1.0e-10:
        if vrblvl > 0:
            print('Found 4 solutions and error is okay.')
        return 0
    if len(sols) != 4:
        if vrblvl > 0:
            print('Number of solutions is not 4 :', len(sols))
        return 1
    if abs(err.real + err.imag) >= 1.0e-10:
        if verblvl > 0:
            print('The error is too large.')
    return 1

def main():
    """
    Runs some tests on solutions.
    """
    lvl = 10
    fail = test_double_track(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()
