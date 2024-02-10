"""
The module trackers offers functions to track paths starting
at the known solutions of a start system and leading to
the desired solutions of a target system,
with aposteriori step size control algorithms.
An aposteriori step size control algorithm determines the step size
based on the performance of the corrector.
For small problems, the default values of the parameters and tolerances 
for predictor and corrector suffice, otherwise they must be tuned.
Reruns of paths must happen with the same value of the gamma constant.
"""
from ctypes import c_int32, c_double, pointer
from phcpy.version import get_phcfun
from phcpy.polynomials import number_of_symbols
from phcpy.solutions import set_double_solutions
from phcpy.solutions import set_double_double_solutions
from phcpy.solutions import set_quad_double_solutions
from phcpy.solutions import clear_double_solutions
from phcpy.solutions import clear_double_double_solutions
from phcpy.solutions import clear_quad_double_solutions
from phcpy.solutions import get_next_double_solution
from phcpy.solutions import get_next_double_double_solution
from phcpy.solutions import get_next_quad_double_solution
from phcpy.solutions import strsol2dict, verify
from phcpy.homotopies import set_double_target_system
from phcpy.homotopies import set_double_start_system
from phcpy.homotopies import set_double_laurent_target_system
from phcpy.homotopies import set_double_laurent_start_system
from phcpy.homotopies import set_double_start_solutions
from phcpy.homotopies import set_double_homotopy
from phcpy.homotopies import set_double_laurent_homotopy
from phcpy.homotopies import get_double_target_solutions
from phcpy.homotopies import clear_double_homotopy
from phcpy.homotopies import clear_double_laurent_homotopy
from phcpy.homotopies import set_double_double_target_system
from phcpy.homotopies import set_double_double_start_system
from phcpy.homotopies import set_double_double_laurent_target_system
from phcpy.homotopies import set_double_double_laurent_start_system
from phcpy.homotopies import set_double_double_start_solutions
from phcpy.homotopies import set_double_double_homotopy
from phcpy.homotopies import set_double_double_laurent_homotopy
from phcpy.homotopies import get_double_double_target_solutions
from phcpy.homotopies import clear_double_double_homotopy
from phcpy.homotopies import clear_double_double_laurent_homotopy
from phcpy.homotopies import set_quad_double_target_system
from phcpy.homotopies import set_quad_double_start_system
from phcpy.homotopies import set_quad_double_laurent_target_system
from phcpy.homotopies import set_quad_double_laurent_start_system
from phcpy.homotopies import set_quad_double_start_solutions
from phcpy.homotopies import set_quad_double_homotopy
from phcpy.homotopies import set_quad_double_laurent_homotopy
from phcpy.homotopies import get_quad_double_target_solutions
from phcpy.homotopies import clear_quad_double_homotopy
from phcpy.homotopies import clear_quad_double_laurent_homotopy
from phcpy.starters import total_degree_start_system

def show_parameters(vrblvl=0):
    """
    Displays the current values of the continuation parameters.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in show_parameters ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
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
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(difficulty_level))
    bbb = pointer(c_int32(digits_of_precision))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
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
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
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
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(idx))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(value))
    vrb = c_int32(vrblvl-1)
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
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(idx))
    bbb = pointer(c_int32(0))
    value = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_parameter_value calls phc', end='')
    retval = phc(189, aaa, bbb, value, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('Value of parameter', idx, ':', value[0])
    return value[0]

def write_parameters(vrblvl=0):
    """
    Writes the parameters with repeated calls to get_parameter_value(),
    as the show_parameters() does not work in a Jupyter notebook.
    """
    pars = [get_parameter_value(idx, vrblvl) for idx in range(35)]
    print('GLOBAL MONITOR :')
    print('  1. the condition of the homotopy           :', int(pars[1]))
    print('  2. number of paths tracked simultaneously  :', int(pars[2]))
    print('  3. maximum number of steps along a path    :', int(pars[3]))
    print('  4. distance from target to start end game  :', end='')
    print(f' {pars[4]:1.3e}')
    print('  5. order of extrapolator in end game       :', int(pars[5]))
    print('  6. maximum number of re-runs               :', int(pars[6]))
    print('STEP CONTROL (PREDICTOR) :                    along path : end game')
    s1s = '( x:Sec,t:Rea )' # secant for x, real for t
    s1c = '( x:Sec,t:Com )' # secant for x, complex for t
    s1g = '( x:Sec,t:Geo )' # second for x, geometric for t
    t1s = '( x:Tan,t:Rea )' # tangent for x, real for t
    t1c = '( x:Tan,t:Com )' # tangent for x, complex for t
    t1g = '( x:Tan,t:Geo )' # tangent for x, geometric for t
    h3s = '( x:Her,t:Rea )' # Hermite for x, real for t
    q2s = '( x:Qu2,t:Rea )' # quadratic for for x, real for t
    c3s = '( x:Cub,t:Rea )' # cubic for x, real for t
    ptp = [s1s, s1c, s1g, t1s, t1c, t1g, h3s, q2s, c3s]
    uuu = 'no predictor'
    if int(pars[7]) in range(9):
        if int(pars[8]) in range(9):
            print(f'  7: 8. type {ptp[int(pars[7])]}:{ptp[int(pars[8])]} :',
                  end=' ')
        else:
            print(f'  7: 8. type {ptp[int(pars[7])]}:{uuu} :', end=' ')
    else:
        if int(pars[8]) in range(9):
            print(f'  7: 8. type {uuu}:{ptp[int(pars[8])]} :',
                  end=' ')
        else:
            print(f'  7: 8. type {uuu}:{uuu} :', end=' ')
    print(int(pars[7]), '        :', int(pars[8]))
    print('  9:10. minimum step size                    :', end='')
    print(f' {pars[9]:1.3e} : {pars[10]:1.3e}')
    print(' 11:12. maximum step size                    :', end='')
    print(f' {pars[11]:1.3e} : {pars[12]:1.3e}')
    print(' 13:14. reduction factor for step size       :', end='')
    print(f' {pars[13]:1.3e} : {pars[14]:1.3e}')
    print(' 15:16. expansion factor for step size       :', end='')
    print(f' {pars[15]:1.3e} : {pars[16]:1.3e}')
    print(' 17:18. expansion threshold                  :', end=' ')
    print(int(pars[17]), '        :', int(pars[18]))
    print('PATH CLOSENESS (CORRECTOR) :                  along path : end game')
    print(' 19:20. maximum number of iterations         :', end=' ')
    print(int(pars[19]), '        :', int(pars[20]))
    print(' 21:22. relative precision for residuals     :', end='')
    print(f' {pars[21]:1.3e} : {pars[22]:1.3e}')
    print(' 23:24. absolute precision for residuals     :', end='')
    print(f' {pars[23]:1.3e} : {pars[24]:1.3e}')
    print(' 25:26. relative precision for corrections   :', end='')
    print(f' {pars[25]:1.3e} : {pars[26]:1.3e}')
    print(' 27:28. absolute precision for corrections   :', end='')
    print(f' {pars[27]:1.3e} : {pars[28]:1.3e}')
    print('SOLUTION TOLERANCES :                         along path : end game')
    print(' 29:30. inverse condition of Jacobian        :', end='')
    print(f' {pars[29]:1.3e} : {pars[30]:1.3e}')
    print(' 31:32. clustering of solutions              :', end='')
    print(f' {pars[31]:1.3e} : {pars[32]:1.3e}')
    print(' 33:34. solution at infinity                 :', end='')
    print(f' {pars[33]:1.3e} : {pars[34]:1.3e}')
    return 0

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

def clear_double_track_data(vrblvl=0):
    """
    Clears the data allocated for path tracking in double precision
    with an artificial parameter homotopy.
    """
    if vrblvl > 0:
        print('in clear_double_track_data ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
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
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
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
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_quad_double_track_data calls phc', end='')
    retval = phc(248, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_double_laurent_track_data(vrblvl=0):
    """
    Clears the data allocated for path tracking in double precision
    with an artificial parameter Laurent homotopy.
    """
    if vrblvl > 0:
        print('in clear_double_laurent_track_data ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_double_laurent_track_data calls phc', end='')
    retval = phc(771, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_double_double_laurent_track_data(vrblvl=0):
    """
    Clears the data allocated for path tracking in double double precision
    with an artificial parameter Laurent homotopy.
    """
    if vrblvl > 0:
        print('in clear_double_double_laurent_track_data ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_double_double_laurent_track_data calls phc', end='')
    retval = phc(772, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_quad_double_laurent_track_data(vrblvl=0):
    """
    Clears the data allocated for path tracking in quad double precision
    with an artificial parameter Laurent homotopy.
    """
    if vrblvl > 0:
        print('in clear_quad_double_laurent_track_data ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_quad_double_laurent_track_data calls phc', end='')
    retval = phc(773, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def do_double_track(tasks=0, vrblvl=0):
    """
    Calls the path trackers in double precision with a number
    of tasks equal to tasks (no multithreading if zero).
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in do_double_track, tasks :', tasks)
    phc = get_phcfun(vrblvl-1)
    nbtasks = pointer(c_int32(tasks))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
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
    if vrblvl > 0:
        print('in do_double_double_track, tasks :', tasks)
    phc = get_phcfun(vrblvl-1)
    nbtasks = pointer(c_int32(tasks))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
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
    if vrblvl > 0:
        print('in do_quad_double_track, tasks :', tasks)
    phc = get_phcfun(vrblvl-1)
    nbtasks = pointer(c_int32(tasks))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> do_quad_double_track calls phc', end='')
    retval = phc(246, nbtasks, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def do_double_laurent_track(tasks=0, vrblvl=0):
    """
    Calls the path trackers in double precision with a number
    of tasks equal to tasks (no multithreading if zero),
    for a homotopy of Laurent polynomial systems.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in do_double_laurent_track, tasks :', tasks)
    phc = get_phcfun(vrblvl-1)
    nbtasks = pointer(c_int32(tasks))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> do_double_laurent_track calls phc', end='')
    retval = phc(774, nbtasks, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def do_double_double_laurent_track(tasks=0, vrblvl=0):
    """
    Calls the path trackers in double double precision with a number
    of tasks equal to tasks (no multithreading if zero),
    for a homotopy of Laurent polynomial sysetms.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in do_double_double_laurent_track, tasks :', tasks)
    phc = get_phcfun(vrblvl-1)
    nbtasks = pointer(c_int32(tasks))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> do_double_double_laurent_track calls phc', end='')
    retval = phc(775, nbtasks, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def do_quad_double_laurent_track(tasks=0, vrblvl=0):
    """
    Calls the path trackers in quad double precision with a number
    of tasks equal to tasks (no multithreading if zero),
    for a homotopy of Laurent polynomial systems.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in do_quad_double_laurent_track, tasks :', tasks)
    phc = get_phcfun(vrblvl-1)
    nbtasks = pointer(c_int32(tasks))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> do_quad_double_laurent_track calls phc', end='')
    retval = phc(776, nbtasks, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_track(target, start, startsols, \
    gamma=0, pwt=2, tasks=0, vrblvl=0):
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
        print(', pwt :', pwt, end='')
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
    set_double_target_system(target, vrblvl-1)
    set_double_start_system(start, vrblvl-1)
    nvr = number_of_symbols(start, vrblvl-1)
    set_double_start_solutions(nvr, startsols, vrblvl-1)
    usedgamma = set_double_homotopy(gamma, pwt, vrblvl-1)
    do_double_track(tasks, vrblvl)
    sols = get_double_target_solutions(vrblvl-1)
    clear_double_solutions(vrblvl-1)
    clear_double_homotopy(vrblvl-1)
    clear_double_track_data(vrblvl)
    return (usedgamma, sols)

def double_double_track(target, start, startsols, \
    gamma=0, pwt=2, tasks=0, vrblvl=0):
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
        print(', pwt :', pwt, end='')
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
    set_double_double_target_system(target, vrblvl-1)
    set_double_double_start_system(start, vrblvl-1)
    nvr = number_of_symbols(start, vrblvl-1)
    set_double_double_start_solutions(nvr, startsols, vrblvl-1)
    usedgamma = set_double_double_homotopy(gamma, pwt, vrblvl-1)
    do_double_double_track(tasks, vrblvl)
    sols = get_double_double_target_solutions(vrblvl-1)
    clear_double_double_solutions(vrblvl-1)
    clear_double_double_homotopy(vrblvl-1)
    clear_double_double_track_data(vrblvl)
    return (usedgamma, sols)

def quad_double_track(target, start, startsols, \
    gamma=0, pwt=2, tasks=0, vrblvl=0):
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
        print(', pwt :', pwt, end='')
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
    set_quad_double_target_system(target, vrblvl-1)
    set_quad_double_start_system(start, vrblvl-1)
    nvr = number_of_symbols(start, vrblvl-1)
    set_quad_double_start_solutions(nvr, startsols, vrblvl-1)
    usedgamma = set_quad_double_homotopy(gamma, pwt, vrblvl-1)
    do_quad_double_track(tasks, vrblvl)
    sols = get_quad_double_target_solutions(vrblvl-1)
    clear_quad_double_homotopy(vrblvl-1)
    clear_quad_double_track_data(vrblvl)
    return (usedgamma, sols)

def double_laurent_track(target, start, startsols, \
    gamma=0, pwt=2, tasks=0, vrblvl=0):
    r"""
    Track paths in double precision, for Laurent systems.
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
        print('in double_laurent_track, with gamma :', gamma, end='')
        print(', pwt :', pwt, end='')
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
    set_double_laurent_target_system(target, vrblvl)
    set_double_laurent_start_system(start, vrblvl)
    nvr = number_of_symbols(start)
    set_double_start_solutions(nvr, startsols, vrblvl)
    usedgamma = set_double_laurent_homotopy(gamma, pwt, vrblvl)
    do_double_laurent_track(tasks, vrblvl)
    sols = get_double_target_solutions(vrblvl)
    clear_double_solutions(vrblvl)
    clear_double_laurent_homotopy(vrblvl)
    clear_double_laurent_track_data(vrblvl)
    return (usedgamma, sols)

def double_double_laurent_track(target, start, startsols, \
    gamma=0, pwt=2, tasks=0, vrblvl=0):
    r"""
    Tracks paths in double double precision, for Laurent systems.
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
        print('in double_double_laurent_track, with gamma :', gamma, end='')
        print(', pwt :', pwt, end='')
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
    set_double_double_laurent_target_system(target, vrblvl-1)
    set_double_double_laurent_start_system(start, vrblvl-1)
    nvr = number_of_symbols(start, vrblvl-1)
    set_double_double_start_solutions(nvr, startsols, vrblvl-1)
    usedgamma = set_double_double_laurent_homotopy(gamma, pwt, vrblvl-1)
    do_double_double_laurent_track(tasks, vrblvl)
    sols = get_double_double_target_solutions(vrblvl-1)
    clear_double_double_solutions(vrblvl-1)
    clear_double_double_laurent_homotopy(vrblvl-1)
    clear_double_double_laurent_track_data(vrblvl)
    return (usedgamma, sols)

def quad_double_laurent_track(target, start, startsols, \
    gamma=0, pwt=2, tasks=0, vrblvl=0):
    r"""
    Tracks paths in quad double precision, for Laurent systems.
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
        print('in quad_double_laurent_track, with gamma :', gamma, end='')
        print(', pwt :', pwt, end='')
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
    set_quad_double_laurent_target_system(target, vrblvl-1)
    set_quad_double_laurent_start_system(start, vrblvl-1)
    nvr = number_of_symbols(start, vrblvl-1)
    set_quad_double_start_solutions(nvr, startsols, vrblvl-1)
    usedgamma = set_quad_double_laurent_homotopy(gamma, pwt, vrblvl-1)
    do_quad_double_laurent_track(tasks, vrblvl)
    sols = get_quad_double_target_solutions(vrblvl-1)
    clear_quad_double_laurent_homotopy(vrblvl-1)
    clear_quad_double_laurent_track_data(vrblvl)
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
    phc = get_phcfun(vrblvl-1)
    afix = pointer(c_int32(int(fixedgamma)))
    bbb = pointer(c_int32(0))
    c_gamma = (c_double*2)()
    c_gamma[0] = c_double(regamma)
    c_gamma[1] = c_double(imgamma)
    ptr_gamma = pointer(c_gamma)
    vrb = c_int32(vrblvl-1)
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
    phc = get_phcfun(vrblvl-1)
    afix = pointer(c_int32(int(fixedgamma)))
    bbb = pointer(c_int32(0))
    c_gamma = (c_double*2)()
    c_gamma[0] = c_double(regamma)
    c_gamma[1] = c_double(imgamma)
    ptr_gamma = pointer(c_gamma)
    vrb = c_int32(vrblvl-1)
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
    phc = get_phcfun(vrblvl-1)
    afix = pointer(c_int32(int(fixedgamma)))
    bbb = pointer(c_int32(0))
    c_gamma = (c_double*2)()
    c_gamma[0] = c_double(regamma)
    c_gamma[1] = c_double(imgamma)
    ptr_gamma = pointer(c_gamma)
    vrb = c_int32(vrblvl-1)
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
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(1))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
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
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(1))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
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
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(1))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
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
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(1))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
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
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(1))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
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
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(1))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
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
    if vrblvl > 0:
        print('in test_double_track ...')
    mickey = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    start, startsols = total_degree_start_system(mickey, vrblvl=vrblvl-1)
    if vrblvl > 0:
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    gamma, sols = double_track(mickey, start, startsols, vrblvl=vrblvl-1)
    if vrblvl > 0:
        print('gamma :', gamma)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(mickey, sols, vrblvl-1)
    if vrblvl > 0:
        print('the error sum :', err)
    if len(sols) == 4 and err < 1.0e-10:
        if vrblvl > 0:
            print('Found 4 solutions and error is okay.')
        return 0
    if len(sols) != 4:
        if vrblvl > 0:
            print('Number of solutions is not 4 :', len(sols))
        return 1
    if err >= 1.0e-10:
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def test_double_double_track(vrblvl=0):
    """
    Tests tracking the mickey mouse example of two quadrics,
    in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_track ...')
    mickey = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    start, startsols = total_degree_start_system(mickey, vrblvl=vrblvl-1)
    if vrblvl > 0:
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    gamma, sols = double_double_track(mickey, start, startsols, \
        vrblvl=vrblvl)
    if vrblvl > 0:
        print('gamma :', gamma)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(mickey, sols, vrblvl-1)
    if vrblvl > 0:
        print('the error sum :', err)
    if len(sols) == 4 and err < 1.0e-10:
        if vrblvl > 0:
            print('Found 4 solutions and error is okay.')
        return 0
    if len(sols) != 4:
        if vrblvl > 0:
            print('Number of solutions is not 4 :', len(sols))
        return 1
    if err >= 1.0e-10:
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def test_quad_double_track(vrblvl=0):
    """
    Tests tracking the mickey mouse example of two quadrics,
    in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_track ...')
    mickey = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    start, startsols = total_degree_start_system(mickey, vrblvl=vrblvl)
    if vrblvl > 0:
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    gamma, sols = quad_double_track(mickey, start, startsols, \
        vrblvl=vrblvl)
    if vrblvl > 0:
        print('gamma :', gamma)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(mickey, sols, vrblvl-1)
    if vrblvl > 0:
        print('the error sum :', err)
    if len(sols) == 4 and err < 1.0e-10:
        if vrblvl > 0:
            print('Found 4 solutions and error is okay.')
        return 0
    if len(sols) != 4:
        if vrblvl > 0:
            print('Number of solutions is not 4 :', len(sols))
        return 1
    if err >= 1.0e-10:
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def test_next_double_track(vrblvl=0):
    """
    Tests the step-by-step tracking on the mickey mouse example
    of two quadrics, in double precision.
    """
    if vrblvl > 0:
        print('in test_next_double_track ...')
    mickey = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    start, startsols = total_degree_start_system(mickey, vrblvl=vrblvl-1)
    if vrblvl > 0:
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    initialize_double_tracker(mickey, start, vrblvl-1)
    initialize_double_solution(2, startsols[0], vrblvl-1)
    autotune_parameters(0, 15, vrblvl-1)
    fail = 1
    for step in range(1, 101):
        sol = next_double_solution(vrblvl-1)
        if vrblvl > 0:
            print('the next solution :')
            print(sol)
        dsol = strsol2dict(sol)
        if vrblvl > 0:
            print('step', step, end='')
            print(', t :', dsol['t'].real)
        if dsol['t'].real >= 1.0:
            if vrblvl > 0:
                print('the end solution :')
                print(sol)
            fail = 0
            break
    clear_double_solutions(vrblvl-1)
    return fail

def test_next_double_double_track(vrblvl=0):
    """
    Tests the step-by-step tracking on the mickey mouse example
    of two quadrics, in double double precision.
    """
    if vrblvl > 0:
        print('in test_next_double_double_track ...')
    mickey = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    start, startsols = total_degree_start_system(mickey, vrblvl=vrblvl-1)
    if vrblvl > 0:
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    initialize_double_double_tracker(mickey, start, vrblvl-1)
    initialize_double_double_solution(2, startsols[0], vrblvl-1)
    autotune_parameters(0, 15, vrblvl-1)
    fail = 1
    for step in range(1, 101):
        sol = next_double_double_solution(vrblvl-1)
        if vrblvl > 0:
            print('the next solution :')
            print(sol)
        dsol = strsol2dict(sol)
        if vrblvl > 0:
            print('step', step, end='')
            print(', t :', dsol['t'].real)
        if dsol['t'].real >= 1.0:
            if vrblvl > 0:
                print('the end solution :')
                print(sol)
            fail = 0
            break
    clear_double_solutions(vrblvl-1)
    clear_double_double_solutions(vrblvl-1)
    return fail

def test_next_quad_double_track(vrblvl=0):
    """
    Tests the step-by-step tracking on the mickey mouse example
    of two quadrics, in quad double precision.
    """
    if vrblvl > 0:
        print('in test_next_quad_double_track ...')
    mickey = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    start, startsols = total_degree_start_system(mickey, vrblvl=vrblvl-1)
    if vrblvl > 0:
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    initialize_quad_double_tracker(mickey, start, vrblvl-1)
    initialize_quad_double_solution(2, startsols[0], vrblvl-1)
    autotune_parameters(0, 15, vrblvl-1)
    fail = 1
    for step in range(1, 101):
        sol = next_quad_double_solution(vrblvl-1)
        if vrblvl > 0:
            print('the next solution :')
            print(sol)
        dsol = strsol2dict(sol)
        if vrblvl > 0:
            print('step', step, end='')
            print(', t :', dsol['t'].real)
        if dsol['t'].real >= 1.0:
            if vrblvl > 0:
                print('the end solution :')
                print(sol)
            fail = 0
            break
    clear_double_solutions(vrblvl-1)
    clear_quad_double_solutions(vrblvl-1)
    return fail

def test_double_laurent_track(vrblvl=0):
    """
    Tests tracking the mickey mouse example of two quadrics,
    set as Laurent system in double precision.
    """
    if vrblvl > 0:
        print('in test_double_laurent_track ...')
    mickey = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    start, startsols = total_degree_start_system(mickey, vrblvl=vrblvl-1)
    if vrblvl > 0:
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    gamma, sols = double_laurent_track(mickey, start, startsols, \
        vrblvl=vrblvl-1)
    if vrblvl > 0:
        print('gamma :', gamma)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(mickey, sols, vrblvl-1)
    if vrblvl > 0:
        print('the error sum :', err)
    if len(sols) == 4 and err < 1.0e-10:
        if vrblvl > 0:
            print('Found 4 solutions and error is okay.')
        return 0
    if len(sols) != 4:
        if vrblvl > 0:
            print('Number of solutions is not 4 :', len(sols))
        return 1
    if err >= 1.0e-10:
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def test_double_double_laurent_track(vrblvl=0):
    """
    Tests tracking the mickey mouse example of two quadrics,
    set as Laurent system in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_laurent_track ...')
    mickey = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    start, startsols = total_degree_start_system(mickey, vrblvl=vrblvl-1)
    if vrblvl > 0:
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    gamma, sols = double_double_laurent_track(mickey, start, startsols, \
        vrblvl=vrblvl)
    if vrblvl > 0:
        print('gamma :', gamma)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(mickey, sols, vrblvl-1)
    if vrblvl > 0:
        print('the error sum :', err)
    if len(sols) == 4 and err < 1.0e-10:
        if vrblvl > 0:
            print('Found 4 solutions and error is okay.')
        return 0
    if len(sols) != 4:
        if vrblvl > 0:
            print('Number of solutions is not 4 :', len(sols))
        return 1
    if err >= 1.0e-10:
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def test_quad_double_laurent_track(vrblvl=0):
    """
    Tests tracking the mickey mouse example of two quadrics,
    set as Laurent system in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_laurent_track ...')
    mickey = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    start, startsols = total_degree_start_system(mickey, vrblvl=vrblvl)
    if vrblvl > 0:
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    gamma, sols = quad_double_laurent_track(mickey, start, startsols, \
        vrblvl=vrblvl)
    if vrblvl > 0:
        print('gamma :', gamma)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(mickey, sols, vrblvl-1)
    if vrblvl > 0:
        print('the error sum :', err)
    if len(sols) == 4 and err < 1.0e-10:
        if vrblvl > 0:
            print('Found 4 solutions and error is okay.')
        return 0
    if len(sols) != 4:
        if vrblvl > 0:
            print('Number of solutions is not 4 :', len(sols))
        return 1
    if err >= 1.0e-10:
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def test_tuning(vrblvl=0):
    """
    Runs some tests on tuning the parameters.
    """
    if vrblvl > 0:
        print('in test_tuning ...')
        show_parameters(vrblvl-1)
        print('setting the condition level to 2 ...')
    fail = set_condition_level(2, vrblvl-1)
    level = get_condition_level(vrblvl-1)
    if vrblvl > 0:
        print('the condition level :', level)
    fail = fail + autotune_parameters(level, 14, vrblvl-1)
    if vrblvl > 0:
        show_parameters(vrblvl-1)
    # interactive_tune(vrblvl)
    fail = fail + autotune_parameters(0, 14, vrblvl-1)
    if vrblvl > 0:
        show_parameters(vrblvl)
    return fail

def test_write_parameters(vrblvl=0):
    """
    Tests the writing of the parameters.
    """
    if vrblvl > 0:
        show_parameters(vrblvl)
        write_parameters(vrblvl)
    return 0

def main():
    """
    Runs some tests on tuning and tracking.
    """
    lvl = 1
    fail = test_tuning(lvl)
    fail = fail + test_double_track(lvl)
    fail = fail + test_double_double_track(lvl)
    fail = fail + test_quad_double_track(lvl)
    fail = fail + test_next_double_track(lvl)
    fail = fail + test_next_double_double_track(lvl)
    fail = fail + test_next_quad_double_track(lvl)
    fail = fail + test_double_laurent_track(lvl)
    fail = fail + test_double_double_laurent_track(lvl)
    fail = fail + test_quad_double_laurent_track(lvl)
    fail = fail + test_write_parameters(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()
