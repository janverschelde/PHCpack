"""
The module trackers offers functions to track paths starting
at the known solutions of a start system and leading to
the desired solutions of a target system,
with aposteriori step control algorithms.
For small problems, the default values of the parameters and tolerances 
for predictor and corrector suffice, otherwise they must be tuned.
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
from phcpy.solutions import clear_double_solutions, verify
from phcpy.solutions import clear_double_double_solutions
from phcpy.solutions import clear_quad_double_solutions
from phcpy.solutions import get_next_double_solution
from phcpy.solutions import get_next_double_double_solution
from phcpy.solutions import get_next_quad_double_solution
from phcpy.homotopies import total_degree_start_system

def show_parameters(vrblvl=0):
    """
    Displays the current values of the continuation parameters.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in show_parameters ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> show_parameters calls phc', end='')
    retval = phc(194, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def autotune_parameters(difficulty_level, digits_of_precision, vrblvl=0):
    """
    Tunes the parameters given the difficulty level of the homotopy
    and the digits of working precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in autotune_parameters, difficulty :', difficulty_level, end='')
        print(', digits :', digits_of_precision)
    phc = get_phcfun()
    aaa = pointer(c_int32(difficulty_level))
    bbb = pointer(c_int32(digits_of_precision))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_parameter_value calls phc', end='')
    retval = phc(193, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def interactive_tune(vrblvl=0):
    """
    Interactive tuning of the parameters.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in interactive_tune ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_parameter_value calls phc', end='')
    retval = phc(70, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_parameter_value(idx, value, vrblvl=0):
    """
    Sets the parameter with index idx to the given value.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_parameter_value, idx :', idx, end='')
        print(', value :', value)
    phc = get_phcfun()
    aaa = pointer(c_int32(idx))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(value))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_parameter_value calls phc', end='')
    retval = phc(190, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def get_parameter_value(idx, vrblvl=0):
    """
    Returns the value of the parameter with index idx.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in get_parameter_value, idx :', idx)
    phc = get_phcfun()
    aaa = pointer(c_int32(idx))
    bbb = pointer(c_int32(0))
    value = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_parameter_value calls phc', end='')
    retval = phc(189, aaa, bbb, value, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('Value of parameter', idx, ':', value[0])
    return value[0]

def set_condition_level(level, vrblvl=0):
    r"""
    Sets the parameter that represents the difficulty level of the
    homotopy to the value of *level*.  The default level equals zero,
    higher values lead to smaller tolerances.
    On return is the failure code, which is zero if all went well.
    """
    return set_parameter_value(1, level, vrblvl)

def get_condition_level(vrblvl=0):
    """
    Returns the level of difficulty.
    The verbose level is given by vrblvl.
    """
    return int(get_parameter_value(1, vrblvl))

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

def clear_double_track_data(vrblvl=0):
    """
    Clears the data allocated for path tracking in double precision
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

def clear_double_double_track_data(vrblvl=0):
    """
    Clears the data allocated for path tracking in double double precision
    with an artificial parameter homotopy.
    """
    if vrblvl > 0:
         print('in clear_double_double_track_data ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> clear_double_double_track_data calls phc', end='')
    retval = phc(238, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_quad_double_track_data(vrblvl=0):
    """
    Clears the data allocated for path tracking in quad double precision
    with an artificial parameter homotopy.
    """
    if vrblvl > 0:
         print('in clear_quad_double_track_data ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> clear_quad_double_track_data calls phc', end='')
    retval = phc(248, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def do_double_track(tasks=0, vrblvl=0):
    """
    Calls the path trackers in double precision with a number
    of tasks equal to tasks (no multithreading if zero).
    The verbose level is given by vrblvl.
    """
    phc = get_phcfun()
    nbtasks = pointer(c_int32(tasks))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> do_double_track calls phc', end='')
    retval = phc(16, nbtasks, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def do_double_double_track(tasks=0, vrblvl=0):
    """
    Calls the path trackers in double double precision with a number
    of tasks equal to tasks (no multithreading if zero).
    The verbose level is given by vrblvl.
    """
    phc = get_phcfun()
    nbtasks = pointer(c_int32(tasks))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> do_double_double_track calls phc', end='')
    retval = phc(236, nbtasks, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def do_quad_double_track(tasks=0, vrblvl=0):
    """
    Calls the path trackers in quad double precision with a number
    of tasks equal to tasks (no multithreading if zero).
    The verbose level is given by vrblvl.
    """
    phc = get_phcfun()
    nbtasks = pointer(c_int32(tasks))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> do_quad_double_track calls phc', end='')
    retval = phc(246, nbtasks, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_track(target, start, startsols, gamma=0, pwt=2, tasks=0, vrblvl=0):
    r"""
    Track paths in double precision.
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
        print('in double_track, with gamma :', gamma, end='')
        print(', pwt :', pwt, end=''),
        print(', tasks :', tasks)
        print('the target system :')
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
    nvr = number_of_symbols(start)
    set_double_start_solutions(nvr, startsols, vrblvl)
    usedgamma = set_double_homotopy(gamma, pwt, vrblvl)
    do_double_track(tasks, vrblvl)
    sols = get_double_target_solutions(vrblvl)
    clear_double_solutions(vrblvl)
    clear_double_homotopy(vrblvl)
    clear_double_track_data(vrblvl)
    return (usedgamma, sols)

def double_double_track(target, start, startsols, \
    gamma=0, pwt=2, tasks=1, vrblvl=0):
    r"""
    Tracks paths in double double precision.
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
    Note: tasks=0 does not work ...
    """
    if vrblvl > 0:
        print('in double_double_track, with gamma :', gamma, end='')
        print(', pwt :', pwt, end=''),
        print(', tasks :', tasks)
        print('the target system :')
        for pol in target:
            print(pol)
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    set_double_double_target_system(target, vrblvl)
    set_double_double_start_system(start, vrblvl)
    nvr = number_of_symbols(start)
    set_double_double_start_solutions(nvr, startsols, vrblvl)
    usedgamma = set_double_double_homotopy(gamma, pwt, vrblvl)
    do_double_double_track(tasks, vrblvl)
    sols = get_double_double_target_solutions(vrblvl)
    clear_double_double_solutions(vrblvl)
    clear_double_double_homotopy(vrblvl)
    clear_double_double_track_data(vrblvl)
    return (usedgamma, sols)

def quad_double_track(target, start, startsols, \
    gamma=0, pwt=2, tasks=1, vrblvl=0):
    r"""
    Tracks paths in quad double precision.
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
    Note: tasks=0 does not work ...
    """
    if vrblvl > 0:
        print('in quad_double_track, with gamma :', gamma, end='')
        print(', pwt :', pwt, end=''),
        print(', tasks :', tasks)
        print('the target system :')
        for pol in target:
            print(pol)
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    set_quad_double_target_system(target, vrblvl)
    set_quad_double_start_system(start, vrblvl)
    nvr = number_of_symbols(start)
    set_quad_double_start_solutions(nvr, startsols, vrblvl)
    usedgamma = set_quad_double_homotopy(gamma, pwt, vrblvl)
    do_quad_double_track(tasks, vrblvl)
    sols = get_quad_double_target_solutions(vrblvl)
    clear_quad_double_homotopy(vrblvl)
    clear_quad_double_track_data(vrblvl)
    return (usedgamma, sols)

def initialize_double_tracker(target, start, fixedgamma=True, \
    regamma=0.0, imgamma=0.0, vrblvl=0):
    r"""
    Initializes a path tracker with a generator for a *target*
    and *start* system given in standard double precision.
    If *fixedgamma*, then gamma will be a fixed default value,
    otherwise, a random complex constant for gamma is generated,
    but only if *regamma* and *imgamma* are both equal to 0.0.
    If not *fixedgamma* and moreover: *regamma* and *imgamma* are not
    both zero, then the complex number with real part in *regamma*
    and imaginary part in *imgamma* will be the gamma constant.
    """
    if vrblvl > 0:
        print('in initialize_double_tracker', end='')
        print(', fixedgamma :', fixedgamma, end='')
        print(', regamma :', regamma, end='')
        print(', imgamma :', imgamma)
        print('the target system :')
        for pol in target:
            print(pol)
        print('the start system :')
        for pol in start:
            print(pol)
    set_double_target_system(target, vrblvl)
    set_double_start_system(start, vrblvl)
    phc = get_phcfun()
    afix = pointer(c_int32(int(fixedgamma)))
    bbb = pointer(c_int32(0))
    c_gamma = (c_double*2)()
    c_gamma[0] = c_double(regamma)
    c_gamma[1] = c_double(imgamma)
    ptr_gamma = pointer(c_gamma)
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> initialize_double_tracker calls phc', end='')
    retval = phc(500, afix, bbb, ptr_gamma, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_double_double_tracker(target, start, fixedgamma=True, \
    regamma=0.0, imgamma=0.0, vrblvl=0):
    r"""
    Initializes a path tracker with a generator for a *target*
    and *start* system given in double double precision.
    If *fixedgamma*, then gamma will be a fixed default value,
    otherwise, a random complex constant for gamma is generated,
    but only if *regamma* and *imgamma* are both equal to 0.0.
    If not *fixedgamma* and moreover: *regamma* and *imgamma* are not
    both zero, then the complex number with real part in *regamma*
    and imaginary part in *imgamma* will be the gamma constant.
    """
    if vrblvl > 0:
        print('in initialize_double_double_tracker', end='')
        print(', fixedgamma :', fixedgamma, end='')
        print(', regamma :', regamma, end='')
        print(', imgamma :', imgamma)
        print('the target system :')
        for pol in target:
            print(pol)
        print('the start system :')
        for pol in start:
            print(pol)
    set_double_double_target_system(target, vrblvl)
    set_double_double_start_system(start, vrblvl)
    phc = get_phcfun()
    afix = pointer(c_int32(int(fixedgamma)))
    bbb = pointer(c_int32(0))
    c_gamma = (c_double*2)()
    c_gamma[0] = c_double(regamma)
    c_gamma[1] = c_double(imgamma)
    ptr_gamma = pointer(c_gamma)
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> initialize_double_double_tracker calls phc', end='')
    retval = phc(501, afix, bbb, ptr_gamma, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_quad_double_tracker(target, start, fixedgamma=True, \
    regamma=0.0, imgamma=0.0, vrblvl=0):
    r"""
    Initializes a path tracker with a generator for a *target*
    and *start* system given in quad double precision.
    If *fixedgamma*, then gamma will be a fixed default value,
    otherwise, a random complex constant for gamma is generated,
    but only if *regamma* and *imgamma* are both equal to 0.0.
    If not *fixedgamma* and moreover: *regamma* and *imgamma* are not
    both zero, then the complex number with real part in *regamma*
    and imaginary part in *imgamma* will be the gamma constant.
    """
    if vrblvl > 0:
        print('in initialize_quad_double_tracker', end='')
        print(', fixedgamma :', fixedgamma, end='')
        print(', regamma :', regamma, end='')
        print(', imgamma :', imgamma)
        print('the target system :')
        for pol in target:
            print(pol)
        print('the start system :')
        for pol in start:
            print(pol)
    set_quad_double_target_system(target, vrblvl)
    set_quad_double_start_system(start, vrblvl)
    phc = get_phcfun()
    afix = pointer(c_int32(int(fixedgamma)))
    bbb = pointer(c_int32(0))
    c_gamma = (c_double*2)()
    c_gamma[0] = c_double(regamma)
    c_gamma[1] = c_double(imgamma)
    ptr_gamma = pointer(c_gamma)
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> initialize_quad_double_tracker calls phc', end='')
    retval = phc(502, afix, bbb, ptr_gamma, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_double_solution(nvr, sol, vrblvl=0):
    r"""
    A double precision path tracker with a generator is
    initialized with a start solution *sol* in a number of
    variables equal to the value of *nvr*.
    """
    if vrblvl > 0:
        print('in initialize_double_solution, nvr :', nvr)
        print('the solution :')
        print(sol)
    set_double_solutions(nvr, [sol], vrblvl)
    phc = get_phcfun()
    aidx = pointer(c_int32(1))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> initialize_double_solution calls phc', end='')
    retval = phc(503, aidx, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_double_double_solution(nvr, sol, vrblvl=0):
    r"""
    A double double precision path tracker with a generator is
    initialized with a start solution *sol* in a number of
    variables equal to the value of *nvr*.
    """
    if vrblvl > 0:
        print('in initialize_double_double_solution, nvr :', nvr)
        print('the solution :')
        print(sol)
    set_double_double_solutions(nvr, [sol], vrblvl)
    phc = get_phcfun()
    aidx = pointer(c_int32(1))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> initialize_double_double_solution calls phc', end='')
    retval = phc(504, aidx, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_quad_double_solution(nvr, sol, vrblvl=0):
    r"""
    A quad double precision path tracker with a generator is
    initialized with a start solution *sol* in a number of
    variables equal to the value of *nvr*.
    """
    if vrblvl > 0:
        print('in initialize_quad_double_solution, nvr :', nvr)
        print('the solution :')
        print(sol)
    set_quad_double_solutions(nvr, [sol], vrblvl)
    phc = get_phcfun()
    aidx = pointer(c_int32(1))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> initialize_quad_double_solution calls phc', end='')
    retval = phc(505, aidx, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def next_double_solution(vrblvl=0):
    r"""
    Returns the next solution on a path tracked with
    double precision arithmetic, provided the functions
    **initialize_double_tracker()** and
    **initialize_double_solution()** have been executed properly.
    """
    if vrblvl > 0:
        print('in next_double_solution ...')
    phc = get_phcfun()
    aidx = pointer(c_int32(1))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> next_double_solution calls phc', end='')
    retval = phc(506, aidx, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    sol = get_next_double_solution(1, vrblvl)
    return sol

def next_double_double_solution(vrblvl=0):
    r"""
    Returns the next solution on a path tracked with
    double double precision arithmetic, provided the functions
    **initialize_double_double_tracker()** and
    **initialize_double_double_solution()** have been executed properly.
    """
    if vrblvl > 0:
        print('in next_double_double_solution ...')
    phc = get_phcfun()
    aidx = pointer(c_int32(1))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> next_double_double_solution calls phc', end='')
    retval = phc(507, aidx, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    sol = get_next_double_double_solution(1, vrblvl)
    return sol

def next_quad_double_solution(vrblvl=0):
    r"""
    Returns the next solution on a path tracked with
    quad double precision arithmetic, provided the functions
    **initialize_quad_double_tracker()** and
    **initialize_quad_double_solution()** have been executed properly.
    """
    if vrblvl > 0:
        print('in next_quad_double_solution ...')
    phc = get_phcfun()
    aidx = pointer(c_int32(1))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> next_quad_double_solution calls phc', end='')
    retval = phc(508, aidx, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    sol = get_next_quad_double_solution(1, vrblvl)
    return sol

def test_double_track(vrblvl=0):
    """
    Tests tracking the mickey mouse example of two quadrics,
    in double precision.
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
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def test_double_double_track(vrblvl=0):
    """
    Tests tracking the mickey mouse example of two quadrics,
    in double double precision.
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
    gamma, sols = double_double_track(mickey, start, startsols, \
        tasks=4, vrblvl=vrblvl)
    ddpols = get_double_double_system(vrblvl)
    for pol in ddpols:
        print(pol)
    write_double_double_solutions(vrblvl)
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
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def test_quad_double_track(vrblvl=0):
    """
    Tests tracking the mickey mouse example of two quadrics,
    in quad double precision.
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
    gamma, sols = quad_double_track(mickey, start, startsols, 
        tasks=2, vrblvl=vrblvl)
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
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def test_next_double_track(vrblvl=0):
    """
    Tests the step-by-step tracking on the mickey mouse example
    of two quadrics, in double precision.
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
    initialize_double_tracker(mickey, start, vrblvl)
    initialize_double_solution(2, startsols[0], vrblvl)
    while True:
        sol = next_double_solution(vrblvl)
        print('the next solution :')
        print(sol)
        answer = input('continue ? (y/n) ')
        if(answer != 'y'):
            break
    clear_double_solutions(vrblvl)
    return 0

def test_next_double_double_track(vrblvl=0):
    """
    Tests the step-by-step tracking on the mickey mouse example
    of two quadrics, in double double precision.
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
    initialize_double_double_tracker(mickey, start, vrblvl)
    initialize_double_double_solution(2, startsols[0], vrblvl)
    while True:
        sol = next_double_double_solution(vrblvl)
        print('the next solution :')
        print(sol)
        answer = input('continue ? (y/n) ')
        if(answer != 'y'):
            break
    clear_double_solutions(vrblvl)
    clear_double_double_solutions(vrblvl)
    return 0

def test_next_quad_double_track(vrblvl=0):
    """
    Tests the step-by-step tracking on the mickey mouse example
    of two quadrics, in quad double precision.
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
    initialize_quad_double_tracker(mickey, start, vrblvl)
    initialize_quad_double_solution(2, startsols[0], vrblvl)
    while True:
        sol = next_quad_double_solution(vrblvl)
        print('the next solution :')
        print(sol)
        answer = input('continue ? (y/n) ')
        if(answer != 'y'):
            break
    clear_double_solutions(vrblvl)
    clear_quad_double_solutions(vrblvl)
    return 0

def test_tuning(vrblvl=0):
    """
    Runs some tests on tuning the parameters.
    """
    show_parameters(vrblvl)
    print('setting the condition level to 2 ...')
    set_condition_level(2, vrblvl)
    level = get_condition_level(vrblvl)
    print('the condition level :', level)
    autotune_parameters(level, 14, vrblvl)
    show_parameters(vrblvl)
    # interactive_tune(vrblvl)
    autotune_parameters(0, 14, vrblvl)
    show_parameters(vrblvl)
    return 0

def main():
    """
    Runs some tests on tuning and tracking.
    """
    lvl = 10
    fail = test_tuning(lvl)
    fail = fail + test_double_track(lvl)
    fail = fail + test_double_double_track(lvl)
    fail = fail + test_quad_double_track(lvl)
    fail = fail + test_next_double_track(lvl)
    fail = fail + test_next_double_double_track(lvl)
    fail = fail + test_next_quad_double_track(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()
