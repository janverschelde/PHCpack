"""
A polynomial homotopy is a family of polynomial systems with one parameter.
In an artificial parameter homotopy, there is a start and a target system.
There is only one system in a natural parameter homotopy,
where one variable plays the role of the parameter in the homotopy.
The module homotopies exports several functions to set start and
target functions in a homotopy.  There are 36 copy functions:

    copy_{double, double_double, quad_double}
        _{laurent_system, system, solutions}
        _{from, into}_{start, target}

which take no input arguments (other than the verbose level),
and have no return arguments (other than the return value of the call).
Those copy functions are auxiliary to the set functions to define
the homotopies.
"""
from ctypes import c_int32, c_double, pointer
from random import uniform
from cmath import exp, pi
from phcpy.version import get_phcfun
from phcpy.polynomials import number_of_symbols
from phcpy.polynomials import set_double_system
from phcpy.polynomials import get_double_system
from phcpy.polynomials import set_double_laurent_system
from phcpy.polynomials import get_double_laurent_system
from phcpy.polynomials import clear_double_system
from phcpy.polynomials import clear_double_laurent_system
from phcpy.polynomials import set_double_double_system
from phcpy.polynomials import get_double_double_system
from phcpy.polynomials import set_double_double_laurent_system
from phcpy.polynomials import get_double_double_laurent_system
from phcpy.polynomials import clear_double_double_system
from phcpy.polynomials import clear_double_double_laurent_system
from phcpy.polynomials import set_quad_double_system
from phcpy.polynomials import get_quad_double_system
from phcpy.polynomials import set_quad_double_laurent_system
from phcpy.polynomials import get_quad_double_laurent_system
from phcpy.polynomials import clear_quad_double_system
from phcpy.polynomials import clear_quad_double_laurent_system
from phcpy.solutions import set_double_solutions, get_double_solutions
from phcpy.solutions import set_double_double_solutions
from phcpy.solutions import get_double_double_solutions
from phcpy.solutions import set_quad_double_solutions
from phcpy.solutions import get_quad_double_solutions
from phcpy.solutions import clear_double_solutions
from phcpy.solutions import clear_double_double_solutions
from phcpy.solutions import clear_quad_double_solutions
from phcpy.solutions import make_solution, verify

def copy_double_system_from_target(vrblvl=0):
    """
    Copies the target system set in an artificial-parameter homotopy
    in double precision into the system in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_system_from_target ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_system_from_target calls phc', end='')
    retval = phc(1, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_system_into_target(vrblvl=0):
    """
    Copies the system set in double precision to the target
    in an artificial-parameter homotopy in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_system_into_target ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_system_into_target calls phc', end='')
    retval = phc(2, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_laurent_system_from_target(vrblvl=0):
    """
    Copies the target Laurent system set in an artificial-parameter homotopy
    in double precision into the Laurent system in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_laurent_system_from_target ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_laurent_system_from_target calls phc', end='')
    retval = phc(786, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_laurent_system_into_target(vrblvl=0):
    """
    Copies the Laurent system set in double precision to the target
    in an artificial-parameter homotopy in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_laurent_system_into_target ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_laurent_system_into_target calls phc', end='')
    retval = phc(780, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_system_from_target(vrblvl=0):
    """
    Copies the target system set in an artificial-parameter homotopy
    in double double precision into the system in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_double_system_from_target ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_double_system_from_target calls phc', end='')
    retval = phc(251, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_system_into_target(vrblvl=0):
    """
    Copies the system set in double double precision to the target
    in an artificial-parameter homotopy in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_double_system_into_target ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_double_system_into_target calls phc', end='')
    retval = phc(252, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_laurent_system_from_target(vrblvl=0):
    """
    Copies the target Laurent system set in an artificial-parameter homotopy
    in double double precision into the Laurent system
    in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_double_laurent_system_from_target ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_double_laurent_system_from_target calls phc', \
            end='')
    retval = phc(787, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_laurent_system_into_target(vrblvl=0):
    """
    Copies the system set in double double precision to the target
    in an artificial-parameter homotopy in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_double_laurent_system_into_target ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_double_laurent_system_into_target calls phc', \
            end='')
    retval = phc(781, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_system_from_target(vrblvl=0):
    """
    Copies the target system set in an artificial-parameter homotopy
    in quad double precision into the system in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_quad_double_system_from_target ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_quad_double_system_from_target calls phc', end='')
    retval = phc(261, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_system_into_target(vrblvl=0):
    """
    Copies the system set in quad double precision to the target
    system in an artificial-parameter homotopy in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_quad_double_system_into_target ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_quad_double_system_into_target calls phc', end='')
    retval = phc(262, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_laurent_system_from_target(vrblvl=0):
    """
    Copies the target Laurent system set in an artificial-parameter homotopy
    in quad double precision into the Laurent system in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_quad_double_laurent_system_from_target ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_quad_double_laurent_system_from_target calls phc', \
            end='')
    retval = phc(788, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_laurent_system_into_target(vrblvl=0):
    """
    Copies the Laurent system set in quad double precision to the target
    system in an artificial-parameter homotopy in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_quad_double_laurent_system_into_target ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_quad_double_laurent_system_into_target calls phc', \
            end='')
    retval = phc(782, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_system_from_start(vrblvl=0):
    """
    Copies the start system set in an artificial-parameter homotopy
    in double precision into the system in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_system_from_start ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_system_from_start calls phc', end='')
    retval = phc(3, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_system_into_start(vrblvl=0):
    """
    Copies the system set in double precision to the start
    system in an artificial-parameter homotopy in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_system_into_start ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_system_into_start calls phc', end='')
    retval = phc(4, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_laurent_system_from_start(vrblvl=0):
    """
    Copies the start Laurent system set in an artificial-parameter homotopy
    in double precision into the Laurent system in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_laurent_system_from_start ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_laurent_system_from_start calls phc', end='')
    retval = phc(783, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_laurent_system_into_start(vrblvl=0):
    """
    Copies the Laurent system set in double precision to the start
    system in an artificial-parameter homotopy in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_laurent_system_into_start ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_laurent_system_into_start calls phc', end='')
    retval = phc(777, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_system_from_start(vrblvl=0):
    """
    Copies the start system set in an artificial-parameter homotopy
    in double double precision into the system in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_double_system_from_start ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_double_system_from_start calls phc', end='')
    retval = phc(253, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_system_into_start(vrblvl=0):
    """
    Copies the system set in double double precision to the start
    in an artificial-parameter homotopy in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_double_system_into_start ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_double_system_into_start calls phc', end='')
    retval = phc(254, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_laurent_system_from_start(vrblvl=0):
    """
    Copies the start Laurent system set in an artificial-parameter homotopy
    in double double precision into the Laurent system
    in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_double_laurent_system_from_start ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_double_laurent_system_from_start calls phc', \
            end='')
    retval = phc(784, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_laurent_system_into_start(vrblvl=0):
    """
    Copies the Laurent system set in double double precision to the start
    in an artificial-parameter homotopy in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_double_laurent_system_into_start ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_double_laurent_system_into_start calls phc', \
            end='')
    retval = phc(778, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_system_from_start(vrblvl=0):
    """
    Copies the start system set in an artificial-parameter homotopy
    in quad double precision into the system in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_quad_double_system_from_start ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_quad_double_system_from_start calls phc', end='')
    retval = phc(263, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_system_into_start(vrblvl=0):
    """
    Copies the system set in quad double precision to the start
    in an artificial-parameter homotopy in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_quad_double_system_into_start ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_quad_double_system_into_start calls phc', end='')
    retval = phc(264, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_laurent_system_from_start(vrblvl=0):
    """
    Copies the start Laurent system set in an artificial-parameter homotopy
    in quad double precision into the Laurent system
    in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_quad_double_laurent_system_from_start ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_quad_double_laurent_system_from_start calls phc', \
            end='')
    retval = phc(785, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_laurent_system_into_start(vrblvl=0):
    """
    Copies the Laurent system set in quad double precision to the start
    in an artificial-parameter homotopy in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_quad_double_laurent_system_into_start ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_quad_double_laurent_system_into_start calls phc', \
            end='')
    retval = phc(779, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_solutions_from_start(vrblvl=0):
    """
    Copies the start solutions in an artificial-parameter homotopy
    in double precision to the solutions in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_solutions_from_start ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_solutions_from_start calls phc', end='')
    retval = phc(7, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_solutions_into_start(vrblvl=0):
    """
    Copies the solutions set in double precision to the start
    solutions in an artificial-parameter homotopy in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_solutions_into_start ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_solutions_into_start calls phc', end='')
    retval = phc(8, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_solutions_from_start(vrblvl=0):
    """
    Copies the start solutions in an artificial-parameter homotopy
    in double double precision to the solutions in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_double_solutions_from_start ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_double_solutions_from_start calls phc', end='')
    retval = phc(257, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_solutions_into_start(vrblvl=0):
    """
    Copies the solutions set in double double precision to the start
    solutions in an artificial-parameter homotopy in double double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_double_solutions_into_start ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_double_solutions_into_start calls phc', end='')
    retval = phc(258, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_solutions_from_start(vrblvl=0):
    """
    Copies the start solutions in an artificial-parameter homotopy
    in quad double precision to the solutions in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_quad_double_solutions_from_start ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_quad_double_solutions_from_start calls phc', end='')
    retval = phc(267, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_solutions_into_start(vrblvl=0):
    """
    Copies the solutions set in quad double precision to the start
    solutions in an artificial-parameter homotopy in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_quad_double_solutions_into_start ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_quad_double_solutions_into_start calls phc', end='')
    retval = phc(268, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_solutions_from_target(vrblvl=0):
    """
    Copies the target solutions in an artificial-parameter homotopy
    in double precision to the solutions in double precision.
    """
    if vrblvl > 0:
        print('in copy_double_solutions_from_target ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_solutions_from_target calls phc', end='')
    retval = phc(5, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_solutions_into_target(vrblvl=0):
    """
    Copies the solutions set in double precision to the target
    solutions in an artificial-parameter homotopy in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_solutions_into_target ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_solutions_into_target calls phc', end='')
    retval = phc(6, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_solutions_from_target(vrblvl=0):
    """
    Copies the target solutions in an artificial-parameter homotopy
    in double precision to the solutions in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_double_solutions_from_target ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_double_solutions_from_target calls phc', end='')
    retval = phc(255, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_solutions_into_target(vrblvl=0):
    """
    Copies the solutions set in double double precision to the target
    solutions in an artificial-parameter homotopy in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_double_double_target_solutions ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_double_target_solutions calls phc', end='')
    retval = phc(256, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_solutions_from_target(vrblvl=0):
    """
    Copies the target solutions in an artificial-parameter homotopy
    in quad double precision to the solutions in quad double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_quad_double_solutions_from_target ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_quad_double_solutions_from_target calls phc', end='')
    retval = phc(265, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_solutions_into_target(vrblvl=0):
    """
    Copies the solutions set in quad double precision to the target
    solutions in an artificial-parameter homotopy in double precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in copy_quad_double_solutions_into_target ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_quad_double_solutions_into_target calls phc', end='')
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
    return copy_double_system_into_target(vrblvl)

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
    return copy_double_double_system_into_target(vrblvl)

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
    return copy_quad_double_system_into_target(vrblvl)

def set_double_laurent_target_system(pols, vrblvl=0):
    """
    Sets the target Laurent system in an artificial parameter homotopy
    in double precision to the list of polynomials in pols,
    which is assumed to be square.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_double_laurent_target_system, with pols :')
        for pol in pols:
            print(pol)
    nvr = number_of_symbols(pols, vrblvl)
    set_double_laurent_system(nvr, pols, vrblvl)
    return copy_double_laurent_system_into_target(vrblvl)

def set_double_double_laurent_target_system(pols, vrblvl=0):
    """
    Sets the target Laurent system in an artificial parameter homotopy
    in double double precision to the list of polynomials in pols,
    which is assumed to be square.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_double_double_laurent_target_system, with pols :')
        for pol in pols:
            print(pol)
    nvr = number_of_symbols(pols, vrblvl)
    set_double_double_laurent_system(nvr, pols, vrblvl)
    return copy_double_double_laurent_system_into_target(vrblvl)

def set_quad_double_laurent_target_system(pols, vrblvl=0):
    """
    Sets the target Laurent system in an artificial parameter homotopy
    in quad double precision to the list of polynomials in pols,
    which is assumed to be square.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_quad_double_target_system, with pols :')
        for pol in pols:
            print(pol)
    nvr = number_of_symbols(pols, vrblvl)
    set_quad_double_laurent_system(nvr, pols, vrblvl)
    return copy_quad_double_laurent_system_into_target(vrblvl)

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
    return copy_double_system_into_start(vrblvl)

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
    return copy_double_double_system_into_start(vrblvl)

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
    return copy_quad_double_system_into_start(vrblvl)

def set_double_laurent_start_system(pols, vrblvl=0):
    """
    Sets the start Laurent system in an artificial parameter homotopy
    in double precision to the list of polynomials in pols,
    which is assumed to be square.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_double_start_system, with pols :')
        for pol in pols:
            print(pol)
    nvr = number_of_symbols(pols, vrblvl)
    set_double_laurent_system(nvr, pols, vrblvl)
    return copy_double_laurent_system_into_start(vrblvl)

def set_double_double_laurent_start_system(pols, vrblvl=0):
    """
    Sets the start Laurent system in an artificial parameter homotopy
    in double double precision to the list of polynomials in pols,
    which is assumed to be square.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_double_double_laurent_start_system, with pols :')
        for pol in pols:
            print(pol)
    nvr = number_of_symbols(pols, vrblvl)
    set_double_double_laurent_system(nvr, pols, vrblvl)
    return copy_double_double_laurent_system_into_start(vrblvl)

def set_quad_double_laurent_start_system(pols, vrblvl=0):
    """
    Sets the start Laurent system in an artificial parameter homotopy
    in quad double precision to the list of polynomials in pols,
    which is assumed to be square.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_quad_double_laurent_start_system, with pols :')
        for pol in pols:
            print(pol)
    nvr = number_of_symbols(pols, vrblvl)
    set_quad_double_laurent_system(nvr, pols, vrblvl)
    return copy_quad_double_laurent_system_into_start(vrblvl)

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
    return copy_double_solutions_into_start(vrblvl)

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
    return copy_double_double_solutions_into_start(vrblvl)

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
    return copy_quad_double_solutions_into_start(vrblvl)

def get_double_target_solutions(vrblvl=0):
    """
    Returns the list of target solutions computed in double precision.
    """
    if vrblvl > 0:
        print('in get_double_target_solutions ...')
    clear_double_solutions(vrblvl)
    copy_double_solutions_from_target(vrblvl)
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
    copy_double_double_solutions_from_target(vrblvl)
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
    copy_quad_double_solutions_from_target(vrblvl)
    sols = get_quad_double_solutions(vrblvl)
    if vrblvl > 0:
        print('the target solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    return sols

def set_double_target_solutions(nvr, sols, vrblvl=0):
    """
    Sets the target solutions in an artificial parameter homotopy
    in double precision to the list of solutions in sols,
    where the number of variables in nvr must match the
    dimension of the start system.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_double_target_solutions, with nvr :', nvr)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    clear_double_solutions(vrblvl)
    set_double_solutions(nvr, sols, vrblvl)
    return copy_double_solutions_into_target(vrblvl)

def set_double_double_target_solutions(nvr, sols, vrblvl=0):
    """
    Sets the target solutions in an artificial parameter homotopy
    in double double precision to the list of solutions in sols,
    where the number of variables in nvr must match the
    dimension of the start system.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_double_double_target_solutions, with nvr :', nvr)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    clear_double_double_solutions(vrblvl)
    set_double_double_solutions(nvr, sols, vrblvl)
    return copy_double_double_solutions_into_target(vrblvl)

def set_quad_double_target_solutions(nvr, sols, vrblvl=0):
    """
    Sets the target solutions in an artificial parameter homotopy
    in quad double precision to the list of solutions in sols,
    where the number of variables in nvr must match the
    dimension of the start system.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_quad_double_target_solutions, with nvr :', nvr)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    clear_quad_double_solutions(vrblvl)
    set_quad_double_solutions(nvr, sols, vrblvl)
    return copy_quad_double_solutions_into_target(vrblvl)

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
    phc = get_phcfun(vrblvl-1)
    apwt = pointer(c_int32(pwt))
    bbb = pointer(c_int32(0))
    if gamma != 0:
        usegamma = gamma
    else:
        angle = uniform(0, 2*pi)
        usegamma = exp(angle*complex(0, 1))
        if vrblvl > 0:
            print('random gamma :', usegamma)
    c_gamma = (c_double*2)()
    c_gamma[0] = c_double(usegamma.real)
    c_gamma[1] = c_double(usegamma.imag)
    ptr_gamma = pointer(c_gamma)
    vrb = c_int32(vrblvl-1)
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
    phc = get_phcfun(vrblvl-1)
    apwt = pointer(c_int32(pwt))
    bbb = pointer(c_int32(0))
    if gamma != 0:
        usegamma = gamma
    else:
        angle = uniform(0, 2*pi)
        usegamma = exp(angle*complex(0, 1))
        if vrblvl > 0:
            print('random gamma :', usegamma)
    c_gamma = (c_double*2)()
    c_gamma[0] = c_double(usegamma.real)
    c_gamma[1] = c_double(usegamma.imag)
    ptr_gamma = pointer(c_gamma)
    vrb = c_int32(vrblvl-1)
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
    phc = get_phcfun(vrblvl-1)
    apwt = pointer(c_int32(pwt))
    bbb = pointer(c_int32(0))
    if gamma != 0:
        usegamma = gamma
    else:
        angle = uniform(0, 2*pi)
        usegamma = exp(angle*complex(0, 1))
        if vrblvl > 0:
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

def set_double_laurent_homotopy(gamma=0, pwt=2, vrblvl=0):
    """
    After the target and start Laurent system are set in double precision,
    the homotopy is constructed with either a random gamma constant,
    or with the given complex value of gamma.
    The power of the continuation parameter is given by pwt.
    The gamma used to make the Laurent homotopy is returned.
    """
    if vrblvl > 0:
        print('in set_double_laurent_homotopy, with power :', pwt)
        print('gamma :', gamma)
    phc = get_phcfun(vrblvl-1)
    apwt = pointer(c_int32(pwt))
    bbb = pointer(c_int32(0))
    if gamma != 0:
        usegamma = gamma
    else:
        angle = uniform(0, 2*pi)
        usegamma = exp(angle*complex(0, 1))
        if vrblvl > 0:
            print('random gamma :', usegamma)
    c_gamma = (c_double*2)()
    c_gamma[0] = c_double(usegamma.real)
    c_gamma[1] = c_double(usegamma.imag)
    ptr_gamma = pointer(c_gamma)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_laurent_homotopy calls phc', end='')
    retval = phc(921, apwt, bbb, ptr_gamma, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    if vrblvl > 0:
        print('-> set_double_laurent_homotopy calls phc', end='')
        print(' to set gamma', end='')
    aprc = pointer(c_int32(1))
    retval = phc(996, aprc, bbb, ptr_gamma, vrb) # set gamma
    if vrblvl > 0:
        print(', return value :', retval)
    return usegamma

def set_double_double_laurent_homotopy(gamma=0, pwt=2, vrblvl=0):
    """
    After the target and start Laurent system are set 
    in double double precision,
    the homotopy is constructed with either a random gamma constant,
    or with the given complex value of gamma.
    The power of the continuation parameter is given by pwt.
    The gamma used to make the Laurent homotopy is returned.
    """
    if vrblvl > 0:
        print('in set_double_double_laurent_homotopy, with power :', pwt)
        print('gamma :', gamma)
    phc = get_phcfun(vrblvl-1)
    apwt = pointer(c_int32(pwt))
    bbb = pointer(c_int32(0))
    if gamma != 0:
        usegamma = gamma
    else:
        angle = uniform(0, 2*pi)
        usegamma = exp(angle*complex(0, 1))
        if vrblvl > 0:
            print('random gamma :', usegamma)
    c_gamma = (c_double*2)()
    c_gamma[0] = c_double(usegamma.real)
    c_gamma[1] = c_double(usegamma.imag)
    ptr_gamma = pointer(c_gamma)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_double_laurent_homotopy calls phc', end='')
    retval = phc(922, apwt, bbb, ptr_gamma, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    if vrblvl > 0:
        print('-> set_double_double_laurent_homotopy calls phc', end='')
        print(' to set gamma', end='')
    aprc = pointer(c_int32(2))
    retval = phc(996, aprc, bbb, ptr_gamma, vrb) # set gamma
    if vrblvl > 0:
        print(', return value :', retval)
    return usegamma

def set_quad_double_laurent_homotopy(gamma=0, pwt=2, vrblvl=0):
    """
    After the target and start Laurent system are set 
    in quad double precision,
    the homotopy is constructed with either a random gamma constant,
    or with the given complex value of gamma.
    The power of the continuation parameter is given by pwt.
    The gamma used to make the Laurent homotopy is returned.
    """
    if vrblvl > 0:
        print('in set_quad_double_laurent_homotopy, with power :', pwt)
        print('gamma :', gamma)
    phc = get_phcfun(vrblvl-1)
    apwt = pointer(c_int32(pwt))
    bbb = pointer(c_int32(0))
    if gamma != 0:
        usegamma = gamma
    else:
        angle = uniform(0, 2*pi)
        usegamma = exp(angle*complex(0, 1))
        if vrblvl > 0:
            print('random gamma :', usegamma)
    c_gamma = (c_double*2)()
    c_gamma[0] = c_double(usegamma.real)
    c_gamma[1] = c_double(usegamma.imag)
    ptr_gamma = pointer(c_gamma)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_quad_double_laurent_homotopy calls phc', end='')
    retval = phc(923, apwt, bbb, ptr_gamma, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    if vrblvl > 0:
        print('-> set_quad_double_laurent_homotopy calls phc', end='')
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
    phc = get_phcfun(vrblvl-1)
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
    phc = get_phcfun(vrblvl-1)
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
    phc = get_phcfun(vrblvl-1)
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

def clear_double_laurent_homotopy(vrblvl=0):
    """
    Clears the homotopy set in double precision.
    """
    if vrblvl > 0:
        print('in clear_double_laurent_homotopy ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_double_laurent_homotopy calls phc', end='')
    retval = phc(924, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_double_double_laurent_homotopy(vrblvl=0):
    """
    Clears the Laurent homotopy set in double double precision.
    """
    if vrblvl > 0:
        print('in clear_double_double_laurent_homotopy ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_double_double_laurent_homotopy calls phc', end='')
    retval = phc(925, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_quad_double_laurent_homotopy(vrblvl=0):
    """
    Clears the Laurent homotopy set in quad double precision.
    """
    if vrblvl > 0:
        print('in clear_quad_double_laurent_homotopy ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_quad_double_laurent_homotopy calls phc', end='')
    retval = phc(926, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def example_start_system(vrblvl=0):
    """
    Returns an example of a start system.
    """
    if vrblvl > 0:
        print('in example_start_system ...')
    pols = ['x^2 - 1;', 'y - 1;']
    names = ['x', 'y']
    sol1 = make_solution(names, [1, 1], vrblvl=vrblvl-1)
    sol2 = make_solution(names, [-1, 1], vrblvl=vrblvl-1)
    sols = [sol1, sol2]
    return (pols, sols)

def example_target_system(vrblvl=0):
    """
    Returns an example of a target system.
    """
    if vrblvl > 0:
        print('in example_start_system ...')
    pols = ['(x - 2)*(y - 2);', 'x + y - 2;']
    names = ['x', 'y']
    sol1 = make_solution(names, [2, 0], vrblvl=vrblvl-1)
    sol2 = make_solution(names, [0, 2], vrblvl=vrblvl-1)
    sols = [sol1, sol2]
    return (pols, sols)

def test_double_start_system(vrblvl=0):
    """
    Tests the definition of a start system in double precision.
    """
    if vrblvl > 0:
        print('in test_double_start_system ...')
    (pols, sols) = example_start_system(vrblvl)
    if vrblvl > 0:
        print('the start system :')
        for pol in pols:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(pols, sols, vrblvl=vrblvl-1)
    if vrblvl > 0:
        print('the error of the start solutions :', err)
    fail = int(abs(err.real) + abs(err.imag) > 1.0e-8)
    set_double_start_system(pols, vrblvl-1)
    clear_double_system(vrblvl-1)
    copy_double_system_from_start(vrblvl-1)
    pols = get_double_system(vrblvl-1)
    if vrblvl > 0:
        print('retrieved polynomials :')
        for pol in pols:
            print(pol)
    fail = fail + int(len(pols) != 2)
    set_double_start_solutions(2, sols, vrblvl-1)
    clear_double_solutions(vrblvl-1)
    copy_double_solutions_from_start(vrblvl-1)
    sols = get_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('retrieved solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    fail = fail + int(len(sols) != 2)
    return fail

def test_double_double_start_system(vrblvl=0):
    """
    Tests the definition of a start system in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_start_system ...')
    (pols, sols) = example_start_system(vrblvl)
    err = verify(pols, sols, vrblvl-1)
    if vrblvl > 0:
        print('the error of the start solutions :', err)
    fail = int(abs(err.real) + abs(err.imag) > 1.0e-8)
    set_double_double_start_system(pols, vrblvl)
    clear_double_double_system(vrblvl-1)
    copy_double_double_system_from_start(vrblvl)
    pols = get_double_double_system(vrblvl-1)
    if vrblvl > 0:
        print('retrieved polynomials :')
        for pol in pols:
            print(pol)
    fail = fail + int(len(pols) != 2)
    set_double_double_start_solutions(2, sols, vrblvl)
    clear_double_double_solutions(vrblvl-1)
    copy_double_double_solutions_from_start(vrblvl)
    sols = get_double_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('retrieved solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    fail = fail + int(len(sols) != 2)
    return fail

def test_quad_double_start_system(vrblvl=0):
    """
    Tests the definition of a start system in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_start_system ...')
    (pols, sols) = example_start_system(vrblvl)
    err = verify(pols, sols, vrblvl-1)
    if vrblvl > 0:
        print('the error of the start solutions :', err)
    fail = int(abs(err.real) + abs(err.imag) > 1.0e-8)
    set_quad_double_start_system(pols, vrblvl)
    clear_quad_double_system(vrblvl-1)
    copy_quad_double_system_from_start(vrblvl)
    pols = get_quad_double_system(vrblvl-1)
    if vrblvl > 0:
        print('retrieved polynomials :')
        for pol in pols:
            print(pol)
    fail = fail + int(len(pols) != 2)
    set_quad_double_start_solutions(2, sols, vrblvl)
    clear_quad_double_solutions(vrblvl-1)
    copy_quad_double_solutions_from_start(vrblvl)
    sols = get_quad_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('retrieved solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    fail = fail + int(len(sols) != 2)
    return fail

def test_double_target_system(vrblvl=0):
    """
    Tests the definition of a target system in double precision.
    """
    if vrblvl > 0:
        print('in test_double_target_system ...')
    (pols, sols) = example_target_system(vrblvl)
    err = verify(pols, sols, vrblvl-1)
    if vrblvl > 0:
        print('the error of the target solutions :', err)
    fail = int(abs(err.real) + abs(err.imag) > 1.0e-8)
    set_double_target_system(pols, vrblvl)
    clear_double_system(vrblvl-1)
    copy_double_system_from_target(vrblvl)
    pols = get_double_system(vrblvl-1)
    if vrblvl > 0:
        print('retrieved polynomials :')
        for pol in pols:
            print(pol)
    fail = fail + int(len(pols) != 2)
    set_double_target_solutions(2, sols, vrblvl)
    clear_double_solutions(vrblvl)
    sols = get_double_target_solutions(vrblvl)
    if vrblvl > 0:
        print('retrieved solutions :')
        for (idx, sol) in enumerate(sols):
            print('solution', idx+1, ':')
            print(sol)
    fail = fail + int(len(sols) != 2)
    return fail

def test_double_double_target_system(vrblvl=0):
    """
    Tests the definition of a target system in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_target_system ...')
    (pols, sols) = example_target_system(vrblvl)
    err = verify(pols, sols, vrblvl-1)
    if vrblvl > 0:
        print('the error of the target solutions :', err)
    fail = int(abs(err.real) + abs(err.imag) > 1.0e-8)
    set_double_double_target_system(pols, vrblvl)
    clear_double_double_system(vrblvl-1)
    copy_double_double_system_from_target(vrblvl)
    pols = get_double_double_system(vrblvl-1)
    if vrblvl > 0:
        print('retrieved polynomials :')
        for pol in pols:
            print(pol)
    fail = fail + int(len(pols) != 2)
    set_double_double_target_solutions(2, sols, vrblvl)
    clear_double_double_solutions(vrblvl)
    sols = get_double_double_target_solutions(vrblvl)
    if vrblvl > 0:
        print('retrieved solutions :')
        for (idx, sol) in enumerate(sols):
            print('solution', idx+1, ':')
            print(sol)
    fail = fail + int(len(sols) != 2)
    return fail

def test_quad_double_target_system(vrblvl=0):
    """
    Tests the definition of a target system in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_target_system ...')
    (pols, sols) = example_target_system(vrblvl)
    err = verify(pols, sols, vrblvl-1)
    if vrblvl > 0:
        print('the error of the target solutions :', err)
    fail = int(abs(err.real) + abs(err.imag) > 1.0e-8)
    set_quad_double_target_system(pols, vrblvl)
    clear_quad_double_system(vrblvl-1)
    copy_quad_double_system_from_target(vrblvl)
    pols = get_quad_double_system(vrblvl-1)
    if vrblvl > 0:
        print('retrieved polynomials :')
        for pol in pols:
            print(pol)
    fail = fail + int(len(pols) != 2)
    set_quad_double_target_solutions(2, sols, vrblvl)
    clear_quad_double_solutions(vrblvl)
    sols = get_quad_double_target_solutions(vrblvl)
    if vrblvl > 0:
        print('retrieved solutions :')
        for (idx, sol) in enumerate(sols):
            print('solution', idx+1, ':')
            print(sol)
    fail = fail + int(len(sols) != 2)
    return fail

def test_double_laurent_start_system(vrblvl=0):
    """
    Tests the definition of a start system in double precision,
    set as a Laurent system.
    """
    if vrblvl > 0:
        print('in test_double_laurent_start_system ...')
    (pols, sols) = example_start_system(vrblvl)
    if vrblvl > 0:
        print('the start system :')
        for pol in pols:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(pols, sols, vrblvl=vrblvl-1)
    if vrblvl > 0:
        print('the error of the start solutions :', err)
    fail = int(abs(err.real) + abs(err.imag) > 1.0e-8)
    set_double_laurent_start_system(pols, vrblvl-1)
    clear_double_laurent_system(vrblvl-1)
    copy_double_laurent_system_from_start(vrblvl-1)
    pols = get_double_laurent_system(vrblvl-1)
    if vrblvl > 0:
        print('retrieved polynomials :')
        for pol in pols:
            print(pol)
    fail = fail + int(len(pols) != 2)
    set_double_start_solutions(2, sols, vrblvl-1)
    clear_double_solutions(vrblvl-1)
    copy_double_solutions_from_start(vrblvl-1)
    sols = get_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('retrieved solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    fail = fail + int(len(sols) != 2)
    return fail

def test_double_double_laurent_start_system(vrblvl=0):
    """
    Tests the definition of a start system in double double precision,
    set as Laurent system.
    """
    if vrblvl > 0:
        print('in test_double_double_laurent_start_system ...')
    (pols, sols) = example_start_system(vrblvl)
    err = verify(pols, sols, vrblvl-1)
    if vrblvl > 0:
        print('the error of the start solutions :', err)
    fail = int(abs(err.real) + abs(err.imag) > 1.0e-8)
    set_double_double_laurent_start_system(pols, vrblvl)
    clear_double_double_laurent_system(vrblvl-1)
    copy_double_double_laurent_system_from_start(vrblvl)
    pols = get_double_double_laurent_system(vrblvl-1)
    if vrblvl > 0:
        print('retrieved polynomials :')
        for pol in pols:
            print(pol)
    fail = fail + int(len(pols) != 2)
    set_double_double_start_solutions(2, sols, vrblvl)
    clear_double_double_solutions(vrblvl-1)
    copy_double_double_solutions_from_start(vrblvl)
    sols = get_double_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('retrieved solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    fail = fail + int(len(sols) != 2)
    return fail

def test_quad_double_laurent_start_system(vrblvl=0):
    """
    Tests the definition of a start system in quad double precision,
    set as Laurent system.
    """
    if vrblvl > 0:
        print('in test_quad_double_laurent_start_system ...')
    (pols, sols) = example_start_system(vrblvl)
    err = verify(pols, sols, vrblvl-1)
    if vrblvl > 0:
        print('the error of the start solutions :', err)
    fail = int(abs(err.real) + abs(err.imag) > 1.0e-8)
    set_quad_double_laurent_start_system(pols, vrblvl)
    clear_quad_double_laurent_system(vrblvl-1)
    copy_quad_double_laurent_system_from_start(vrblvl)
    pols = get_quad_double_laurent_system(vrblvl-1)
    if vrblvl > 0:
        print('retrieved polynomials :')
        for pol in pols:
            print(pol)
    fail = fail + int(len(pols) != 2)
    set_quad_double_start_solutions(2, sols, vrblvl)
    clear_quad_double_solutions(vrblvl-1)
    copy_quad_double_solutions_from_start(vrblvl)
    sols = get_quad_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('retrieved solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    fail = fail + int(len(sols) != 2)
    return fail

def test_double_laurent_target_system(vrblvl=0):
    """
    Tests the definition of a target system in double precision,
    set as Laurent system.
    """
    if vrblvl > 0:
        print('in test_double_laurent_target_system ...')
    (pols, sols) = example_target_system(vrblvl)
    err = verify(pols, sols, vrblvl-1)
    if vrblvl > 0:
        print('the error of the target solutions :', err)
    fail = int(abs(err.real) + abs(err.imag) > 1.0e-8)
    set_double_laurent_target_system(pols, vrblvl)
    clear_double_laurent_system(vrblvl-1)
    copy_double_laurent_system_from_target(vrblvl)
    pols = get_double_laurent_system(vrblvl-1)
    if vrblvl > 0:
        print('retrieved polynomials :')
        for pol in pols:
            print(pol)
    fail = fail + int(len(pols) != 2)
    set_double_target_solutions(2, sols, vrblvl)
    clear_double_solutions(vrblvl)
    sols = get_double_target_solutions(vrblvl)
    if vrblvl > 0:
        print('retrieved solutions :')
        for (idx, sol) in enumerate(sols):
            print('solution', idx+1, ':')
            print(sol)
    fail = fail + int(len(sols) != 2)
    return fail

def test_double_double_laurent_target_system(vrblvl=0):
    """
    Tests the definition of a target system in double double precision,
    set as Laurent system.
    """
    if vrblvl > 0:
        print('in test_double_double_laurent_target_system ...')
    (pols, sols) = example_target_system(vrblvl)
    err = verify(pols, sols, vrblvl-1)
    if vrblvl > 0:
        print('the error of the target solutions :', err)
    fail = int(abs(err.real) + abs(err.imag) > 1.0e-8)
    set_double_double_laurent_target_system(pols, vrblvl)
    clear_double_double_laurent_system(vrblvl-1)
    copy_double_double_laurent_system_from_target(vrblvl)
    pols = get_double_double_laurent_system(vrblvl-1)
    if vrblvl > 0:
        print('retrieved polynomials :')
        for pol in pols:
            print(pol)
    fail = fail + int(len(pols) != 2)
    set_double_double_target_solutions(2, sols, vrblvl)
    clear_double_double_solutions(vrblvl)
    sols = get_double_double_target_solutions(vrblvl)
    if vrblvl > 0:
        print('retrieved solutions :')
        for (idx, sol) in enumerate(sols):
            print('solution', idx+1, ':')
            print(sol)
    fail = fail + int(len(sols) != 2)
    return fail

def test_quad_double_laurent_target_system(vrblvl=0):
    """
    Tests the definition of a target system in quad double precision,
    set as Laurent system.
    """
    if vrblvl > 0:
        print('in test_quad_double_laurent_target_system ...')
    (pols, sols) = example_target_system(vrblvl)
    err = verify(pols, sols, vrblvl-1)
    if vrblvl > 0:
        print('the error of the target solutions :', err)
    fail = int(abs(err.real) + abs(err.imag) > 1.0e-8)
    set_quad_double_laurent_target_system(pols, vrblvl)
    clear_quad_double_laurent_system(vrblvl-1)
    copy_quad_double_laurent_system_from_target(vrblvl)
    pols = get_quad_double_laurent_system(vrblvl-1)
    if vrblvl > 0:
        print('retrieved polynomials :')
        for pol in pols:
            print(pol)
    fail = fail + int(len(pols) != 2)
    set_quad_double_target_solutions(2, sols, vrblvl)
    clear_quad_double_solutions(vrblvl)
    sols = get_quad_double_target_solutions(vrblvl)
    if vrblvl > 0:
        print('retrieved solutions :')
        for (idx, sol) in enumerate(sols):
            print('solution', idx+1, ':')
            print(sol)
    fail = fail + int(len(sols) != 2)
    return fail

def main():
    """
    Runs some tests.
    """
    lvl = 1
    fail = test_double_start_system(lvl)
    fail = fail + test_double_double_start_system(lvl)
    fail = fail + test_quad_double_start_system(lvl)
    fail = fail + test_double_target_system(lvl)
    fail = fail + test_double_double_target_system(lvl)
    fail = fail + test_quad_double_target_system(lvl)
    fail = fail + test_double_laurent_start_system(lvl)
    fail = fail + test_double_double_laurent_start_system(lvl)
    fail = fail + test_quad_double_laurent_start_system(lvl)
    fail = fail + test_double_laurent_target_system(lvl)
    fail = fail + test_double_double_laurent_target_system(lvl)
    fail = fail + test_quad_double_laurent_target_system(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()
