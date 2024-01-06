"""
A polynomial homotopy is a family of polynomial systems with one parameter.
In an artificial parameter homotopy, there is a start and a target system.
There is only one system in a natural parameter homotopy,
where one variable plays the role of the parameter in the homotopy.
The module homotopies exports several functions to set start and
target functions in a homotopy.
"""
from ctypes import c_int32, c_double, pointer
from random import uniform
from cmath import exp, pi
from phcpy.version import get_phcfun
from phcpy.polynomials import number_of_symbols
from phcpy.polynomials import set_double_system
from phcpy.polynomials import set_double_double_system
from phcpy.polynomials import get_double_double_system
from phcpy.polynomials import set_quad_double_system
from phcpy.solutions import set_double_solutions, get_double_solutions
from phcpy.solutions import set_double_double_solutions
from phcpy.solutions import get_double_double_solutions
from phcpy.solutions import write_double_double_solutions
from phcpy.solutions import set_quad_double_solutions
from phcpy.solutions import get_quad_double_solutions
from phcpy.solutions import clear_double_solutions
from phcpy.solutions import clear_double_double_solutions
from phcpy.solutions import clear_quad_double_solutions
from phcpy.solutions import get_next_double_solution
from phcpy.solutions import get_next_double_double_solution
from phcpy.solutions import get_next_quad_double_solution


def copy_double_target_system(vrblvl=0):
    """
    Copies the system set in double precision to the target
    in an artificial-parameter homotopy in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_target_system ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> copyt_double_target_system calls phc', end='')
    retval = phc(2, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_target_system(vrblvl=0):
    """
    Copies the system set in double double precision to the target
    in an artificial-parameter homotopy in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_double_target_system ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_double_target_system calls phc', end='')
    retval = phc(252, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_target_system(vrblvl=0):
    """
    Copies the system set in quad double precision to the target
    in an artificial-parameter homotopy in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_quad_double_target_system ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> copy_quad_double_target_system calls phc', end='')
    retval = phc(262, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_start_system(vrblvl=0):
    """
    Copies the system set in double precision to the start
    in an artificial-parameter homotopy in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_start_system ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> copy_double_start_system calls phc', end='')
    retval = phc(4, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_start_system(vrblvl=0):
    """
    Copies the system set in double double precision to the start
    in an artificial-parameter homotopy in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_double_start_system ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> copy_double_double_start_system calls phc', end='')
    retval = phc(254, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_start_system(vrblvl=0):
    """
    Copies the system set in quad double precision to the start
    in an artificial-parameter homotopy in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_quad_double_start_system ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> copy_quad_double_start_system calls phc', end='')
    retval = phc(264, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_start_solutions(vrblvl=0):
    """
    Copies the solutions set in double precision to the start
    solutions in an artificial-parameter homotopy in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_start_solutions ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> copy_double_start_solutions calls phc', end='')
    retval = phc(8, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_start_solutions(vrblvl=0):
    """
    Copies the solutions set in double double precision to the start
    solutions in an artificial-parameter homotopy in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_double_start_solutions ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> copy_double_double_start_solutions calls phc', end='')
    retval = phc(258, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_start_solutions(vrblvl=0):
    """
    Copies the solutions set in quad double precision to the start
    solutions in an artificial-parameter homotopy in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_quad_double_start_solutions ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> copy_quad_double_start_solutions calls phc', end='')
    retval = phc(268, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_target_solutions(vrblvl=0):
    """
    Copies the solutions set in double precision to the target
    solutions in an artificial-parameter homotopy in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_target_solutions ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> copy_double_target_solutions calls phc', end='')
    retval = phc(6, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_target_solutions(vrblvl=0):
    """
    Copies the solutions set in double double precision to the target
    solutions in an artificial-parameter homotopy in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_double_target_solutions ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> copy_double_double_target_solutions calls phc', end='')
    retval = phc(256, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_target_solutions(vrblvl=0):
    """
    Copies the solutions set in quad double precision to the target
    solutions in an artificial-parameter homotopy in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_quad_double_target_solutions ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> copy_quad_double_target_solutions calls phc', end='')
    retval = phc(266, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

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
    nvr = number_of_symbols(pols, vrblvl)
    set_double_system(nvr, pols, vrblvl)
    return copy_double_target_system(vrblvl)

def set_double_double_target_system(pols, vrblvl=0):
    """
    Sets the target system in an artificial parameter homotopy
    in double double precision to the list of polynomials in pols,
    which is assumed to be square.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_double_double_target_system, with pols :')
        for pol in pols:
            print(pol)
    nvr = number_of_symbols(pols, vrblvl)
    set_double_double_system(nvr, pols, vrblvl)
    return copy_double_double_target_system(vrblvl)

def set_quad_double_target_system(pols, vrblvl=0):
    """
    Sets the target system in an artificial parameter homotopy
    in quad double precision to the list of polynomials in pols,
    which is assumed to be square.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_quad_double_target_system, with pols :')
        for pol in pols:
            print(pol)
    nvr = number_of_symbols(pols, vrblvl)
    set_quad_double_system(nvr, pols, vrblvl)
    return copy_quad_double_target_system(vrblvl)

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
    nvr = number_of_symbols(pols, vrblvl)
    set_double_system(nvr, pols, vrblvl)
    return copy_double_start_system(vrblvl)

def set_double_double_start_system(pols, vrblvl=0):
    """
    Sets the start system in an artificial parameter homotopy
    in double double precision to the list of polynomials in pols,
    which is assumed to be square.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_double_double_start_system, with pols :')
        for pol in pols:
            print(pol)
    nvr = number_of_symbols(pols, vrblvl)
    set_double_double_system(nvr, pols, vrblvl)
    return copy_double_double_start_system(vrblvl)

def set_quad_double_start_system(pols, vrblvl=0):
    """
    Sets the start system in an artificial parameter homotopy
    in quad double precision to the list of polynomials in pols,
    which is assumed to be square.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_quad_double_start_system, with pols :')
        for pol in pols:
            print(pol)
    nvr = number_of_symbols(pols, vrblvl)
    set_quad_double_system(nvr, pols, vrblvl)
    return copy_quad_double_start_system(vrblvl)

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
            print('Solution', idx+1, ':')
            print(sol)
    clear_double_solutions(vrblvl)
    set_double_solutions(nvr, sols, vrblvl)
    return copy_double_start_solutions(vrblvl)

def set_double_double_start_solutions(nvr, sols, vrblvl=0):
    """
    Sets the start solutions in an artificial parameter homotopy
    in double double precision to the list of solutions in sols,
    where the number of variables in nvr must match the
    dimension of the start system.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_double_double_start_solutions, with nvr :', nvr)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    clear_double_double_solutions(vrblvl)
    set_double_double_solutions(nvr, sols, vrblvl)
    return copy_double_double_start_solutions(vrblvl)

def set_quad_double_start_solutions(nvr, sols, vrblvl=0):
    """
    Sets the start solutions in an artificial parameter homotopy
    in quad double precision to the list of solutions in sols,
    where the number of variables in nvr must match the
    dimension of the start system.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_quad_double_start_solutions, with nvr :', nvr)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    clear_quad_double_solutions(vrblvl)
    set_quad_double_solutions(nvr, sols, vrblvl)
    return copy_quad_double_start_solutions(vrblvl)

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
            print('Solution', idx+1, ':')
            print(sol)
    return sols

def get_double_double_target_solutions(vrblvl=0):
    """
    Returns the list of target solutions computed
    in double double precision.
    """
    if vrblvl > 0:
        print('in get_double_double_target_solutions ...')
    clear_double_double_solutions(vrblvl)
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_double_target_solutions calls phc', end='')
    retval = phc(255, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    sols = get_double_double_solutions(vrblvl)
    if vrblvl > 0:
        print('the target solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    return sols

def get_quad_double_target_solutions(vrblvl=0):
    """
    Returns the list of target solutions computed
    in quad double precision.
    """
    if vrblvl > 0:
        print('in get_quad_double_target_solutions ...')
    clear_quad_double_solutions(vrblvl)
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_quad_double_target_solutions calls phc', end='')
    retval = phc(265, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    sols = get_quad_double_solutions(vrblvl)
    if vrblvl > 0:
        print('the target solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
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
    if vrblvl > 0:
        print('-> set_double_homotopy calls phc', end='')
        print(' to set gamma', end='')
    aprc = pointer(c_int32(1))
    retval = phc(996, aprc, bbb, ptr_gamma, vrb) # set gamma
    if vrblvl > 0:
        print(', return value :', retval)
    return usegamma 

def set_double_double_homotopy(gamma=0, pwt=2, vrblvl=0):
    """
    After the target and start system are set in double double precision,
    the homotopy is constructed with either a random gamma constant,
    or with the given complex value of gamma.
    The power of the continuation parameter is given by pwt.
    The gamma used to make the homotopy is returned.
    """
    if vrblvl > 0:
        print('in set_double_double_homotopy, with power :', pwt)
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
        print('-> set_double_double_homotopy calls phc', end='')
    retval = phc(173, apwt, bbb, ptr_gamma, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    if vrblvl > 0:
        print('-> set_double_double_homotopy calls phc', end='')
        print(' to set gamma', end='')
    aprc = pointer(c_int32(2))
    retval = phc(996, aprc, bbb, ptr_gamma, vrb) # set gamma
    if vrblvl > 0:
        print(', return value :', retval)
    return usegamma 

def set_quad_double_homotopy(gamma=0, pwt=2, vrblvl=0):
    """
    After the target and start system are set in quad double precision,
    the homotopy is constructed with either a random gamma constant,
    or with the given complex value of gamma.
    The power of the continuation parameter is given by pwt.
    The gamma used to make the homotopy is returned.
    """
    if vrblvl > 0:
        print('in set_quad_double_homotopy, with power :', pwt)
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
        print('-> set_quad_double_homotopy calls phc', end='')
    retval = phc(183, apwt, bbb, ptr_gamma, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    if vrblvl > 0:
        print('-> set_quad_double_homotopy calls phc', end='')
        print(' to set gamma', end='')
    aprc = pointer(c_int32(3))
    retval = phc(996, aprc, bbb, ptr_gamma, vrb) # set gamma
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

def clear_double_double_homotopy(vrblvl=0):
    """
    Clears the homotopy set in double double precision.
    """
    if vrblvl > 0:
        print('in clear_double_double_homotopy ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> clear_double_double_homotopy calls phc', end='')
    retval = phc(174, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_quad_double_homotopy(vrblvl=0):
    """
    Clears the homotopy set in quad double precision.
    """
    if vrblvl > 0:
        print('in clear_quad_double_homotopy ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> clear_quad_double_homotopy calls phc', end='')
    retval = phc(184, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval
