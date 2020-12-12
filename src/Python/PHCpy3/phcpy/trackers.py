"""
The module trackers offers functions to track paths starting
at the known solutions of a start system and leading to
the desired solutions of a target system.
The path tracking functions in this module can track all solution paths
in several levels of precision: standard double, double double,
quad double, or arbitrary multiprecision arithmetic.
For standard double, double double, and quad double arithmetic,
multitasking is supported which could give a good speedup
if sufficiently many cores are available on the processor.
The tuning of the parameters for the predictor, corrector,
and the settings of the tolerances is handled by the tuning module.
"""

def standard_double_track(target, start, sols, gamma=0, pwt=2, tasks=0):
    r"""
    Does path tracking with standard double precision.
    On input are a target system, a start system with solutions,
    optionally: a (random) gamma constant and the number of tasks.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in sols.
    The *sols* is a list of strings representing start solutions.
    By default, a random *gamma* constant is generated,
    otherwise *gamma* must be a nonzero complex constant.
    The *pwt* is the power of t in the homotopy.
    Changing the default of *pwt* can only be done if a nonzero complex
    value for *gamma* is provided as well.
    The number of tasks in the multithreading is defined by *tasks*.
    The default zero value for *tasks* indicates no multithreading.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_start_system
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_start_solutions
    from phcpy.phcpy2c3 import py2c_create_standard_homotopy
    from phcpy.phcpy2c3 import py2c_create_standard_homotopy_with_gamma
    from phcpy.phcpy2c3 import py2c_clear_standard_homotopy
    from phcpy.phcpy2c3 import py2c_solve_by_standard_homotopy_continuation
    from phcpy.phcpy2c3 import py2c_clear_standard_operations_data
    from phcpy.phcpy2c3 import py2c_solcon_clear_standard_solutions
    from phcpy.phcpy2c3 import py2c_copy_standard_target_solutions_to_container
    from phcpy.interface import store_standard_system
    from phcpy.interface import store_standard_solutions
    from phcpy.interface import load_standard_solutions
    from phcpy.solver import number_of_symbols
    dim = number_of_symbols(start)
    store_standard_system(target, nbvar=dim)
    py2c_copy_standard_container_to_target_system()
    store_standard_system(start, nbvar=dim)
    py2c_copy_standard_container_to_start_system()
    if(gamma == 0):
        py2c_create_standard_homotopy()
    else:
        py2c_create_standard_homotopy_with_gamma(gamma.real, gamma.imag, pwt)
    store_standard_solutions(dim, sols)
    py2c_copy_standard_container_to_start_solutions()
    py2c_solve_by_standard_homotopy_continuation(tasks)
    py2c_solcon_clear_standard_solutions()
    py2c_copy_standard_target_solutions_to_container()
    py2c_clear_standard_homotopy()
    py2c_clear_standard_operations_data()
    return load_standard_solutions()

def ade_double_track(target, start, sols, gamma=0, verbose=1):
    r"""
    Does path tracking with algorithmic differentiation,
    in standard double precision.
    On input are a target system, a start system with solutions.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in sols.
    The *sols* is a list of strings representing start solutions.
    On return are the string representations of the solutions
    computed at the end of the paths.
    If *gamma* on input equals zero, then a random complex number is generated,
    otherwise the real and imaginary parts of *gamma* are used.
    """
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_start_system
    from phcpy.phcpy2c3 import py2c_ade_manypaths_d
    from phcpy.interface import store_standard_system
    from phcpy.interface import store_standard_solutions
    from phcpy.interface import load_standard_solutions
    store_standard_system(target)
    py2c_copy_standard_container_to_target_system()
    store_standard_system(start)
    py2c_copy_standard_container_to_start_system()
    dim = len(start)
    store_standard_solutions(dim, sols)
    if(gamma == 0):
        from random import uniform
        from cmath import exp, pi
        angle = uniform(0, 2*pi)
        gamma = exp(angle*complex(0, 1))
        if(verbose > 0):
            print('random gamma constant :', gamma)
    fail = py2c_ade_manypaths_d(verbose, gamma.real, gamma.imag)
    if(fail == 0):
        if(verbose > 0):
            print('Path tracking with AD was a success!')
    else:
        print('Path tracking with AD failed!')
    return load_standard_solutions()

def ade_tuned_double_track(target, start, sols, pars, gamma=0, verbose=1):
    r"""
    Does path tracking with algorithm differentiation,
    in standard double precision, with tuned parameters in pars.
    On input are a target system, a start system with solutions.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in *sols*.
    The *sols* is a list of strings representing start solutions.
    The *pars* is a tuple of 14 values for the path parameters.
    On return are the string representations of the solutions
    computed at the end of the paths.
    If *gamma* on input equals zero, then a random complex number is generated,
    otherwise the real and imaginary parts of *gamma* are used.
    """
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_start_system
    from phcpy.phcpy2c3 import py2c_ade_manypaths_d_pars
    from phcpy.interface import store_standard_system
    from phcpy.interface import store_standard_solutions
    from phcpy.interface import load_standard_solutions
    store_standard_system(target)
    py2c_copy_standard_container_to_target_system()
    store_standard_system(start)
    py2c_copy_standard_container_to_start_system()
    dim = len(start)
    store_standard_solutions(dim, sols)
    if(gamma == 0):
        from random import uniform
        from cmath import exp, pi
        angle = uniform(0, 2*pi)
        gamma = exp(angle*complex(0, 1))
        if(verbose > 0):
            print('random gamma constant :', gamma)
    if(verbose > 0):
        print('the path parameters :\n', pars)
    fail = py2c_ade_manypaths_d_pars(verbose, gamma.real, gamma.imag, \
        pars[0], pars[1], pars[2], pars[3], pars[4], pars[5], pars[6], \
        pars[7], pars[8], pars[9], pars[10], pars[11], pars[12], pars[13])
    if(fail == 0):
        if(verbose > 0):
            print('Path tracking with AD was a success!')
    else:
        print('Path tracking with AD failed!')
    return load_standard_solutions()

def gpu_double_track(target, start, sols, gamma=0, verbose=1):
    r"""
    GPU accelerated path tracking with algorithm differentiation,
    in standard double precision.
    On input are a target system, a start system with solutions.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in sols.
    The *sols* is a list of strings representing start solutions.
    On return are the string representations of the solutions
    computed at the end of the paths.
    If *gamma* on input equals zero, then a random complex number is generated,
    otherwise the real and imaginary parts of *gamma* are used.
    """
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_start_system
    from phcpy.phcpy2c3 import py2c_gpu_manypaths_d
    from phcpy.interface import store_standard_system
    from phcpy.interface import store_standard_solutions
    from phcpy.interface import load_standard_solutions
    store_standard_system(target)
    py2c_copy_standard_container_to_target_system()
    store_standard_system(start)
    py2c_copy_standard_container_to_start_system()
    dim = len(start)
    store_standard_solutions(dim, sols)
    if(gamma == 0):
        from random import uniform
        from cmath import exp, pi
        angle = uniform(0, 2*pi)
        gamma = exp(angle*complex(0, 1))
        if(verbose > 0):
            print('random gamma constant :', gamma)
    fail = py2c_gpu_manypaths_d(2, verbose, gamma.real, gamma.imag)
    if(fail == 0):
        if(verbose > 0):
            print('Path tracking on the GPU was a success!')
    else:
        print('Path tracking on the GPU failed!')
    return load_standard_solutions()

def double_double_track(target, start, sols, gamma=0, pwt=2, tasks=0):
    r"""
    Does path tracking in double double precision.
    On input are a target system, a start system with solutions,
    optionally a (random) gamma constant and the number of tasks.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system with known solutions in sols.
    The *sols* is a list of strings representing start solutions.
    By default, a random *gamma* constant is generated,
    otherwise *gamma* must be a nonzero complex constant.
    The *pwt* is the power of t in the homotopy.
    Changing the default of *pwt* can only be done if a nonzero complex
    value for *gamma* is provided as well.
    The number of tasks in the multithreading is defined by *tasks*.
    The default zero value for *tasks* indicates no multithreading.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_start_system
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_start_solutions
    from phcpy.phcpy2c3 import py2c_create_dobldobl_homotopy
    from phcpy.phcpy2c3 import py2c_create_dobldobl_homotopy_with_gamma
    from phcpy.phcpy2c3 import py2c_clear_dobldobl_homotopy
    from phcpy.phcpy2c3 import py2c_solve_by_dobldobl_homotopy_continuation
    from phcpy.phcpy2c3 import py2c_clear_dobldobl_operations_data
    from phcpy.phcpy2c3 import py2c_solcon_clear_dobldobl_solutions
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_target_solutions_to_container
    from phcpy.interface import store_dobldobl_system
    from phcpy.interface import store_dobldobl_solutions
    from phcpy.interface import load_dobldobl_solutions
    from phcpy.solver import number_of_symbols
    dim = number_of_symbols(start)
    store_dobldobl_system(target, nbvar=dim)
    py2c_copy_dobldobl_container_to_target_system()
    store_dobldobl_system(start, nbvar=dim)
    py2c_copy_dobldobl_container_to_start_system()
    if(gamma == 0):
        py2c_create_dobldobl_homotopy()
    else:
        py2c_create_dobldobl_homotopy_with_gamma(gamma.real, gamma.imag, pwt)
    # dim = len(start)
    store_dobldobl_solutions(dim, sols)
    py2c_copy_dobldobl_container_to_start_solutions()
    py2c_solve_by_dobldobl_homotopy_continuation(tasks)
    py2c_solcon_clear_dobldobl_solutions()
    py2c_copy_dobldobl_target_solutions_to_container()
    py2c_clear_dobldobl_homotopy()
    py2c_clear_dobldobl_operations_data()
    return load_dobldobl_solutions()

def ade_double_double_track(target, start, sols, gamma=0, verbose=1):
    r"""
    Does path tracking with algorithmic differentiation,
    in double double precision.
    On input are a target system, a start system with solutions.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in sols.
    The *sols* is a list of strings representing start solutions.
    On return are the string representations of the solutions
    computed at the end of the paths.
    If *gamma* on input equals zero, then a random complex number is generated,
    otherwise the real and imaginary parts of *gamma* are used.
    """
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_start_system
    from phcpy.phcpy2c3 import py2c_ade_manypaths_dd
    from phcpy.interface import store_dobldobl_system
    from phcpy.interface import store_dobldobl_solutions
    from phcpy.interface import load_dobldobl_solutions
    store_dobldobl_system(target)
    py2c_copy_dobldobl_container_to_target_system()
    store_dobldobl_system(start)
    py2c_copy_dobldobl_container_to_start_system()
    dim = len(start)
    store_dobldobl_solutions(dim, sols)
    if(gamma == 0):
        from random import uniform
        from cmath import exp, pi
        angle = uniform(0, 2*pi)
        gamma = exp(angle*complex(0, 1))
        if(verbose > 0):
            print('random gamma constant :', gamma)
    fail = py2c_ade_manypaths_dd(verbose, gamma.real, gamma.imag)
    if(fail == 0):
        if(verbose > 0):
            print('Path tracking with AD was a success!')
    else:
        print('Path tracking with AD failed!')
    return load_dobldobl_solutions()

def ade_tuned_double_double_track\
    (target, start, sols, pars, gamma=0, verbose=1):
    r"""
    Does path tracking with algorithm differentiation,
    in double double precision, with tuned path parameters.
    On input are a target system, a start system with solutions.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in *sols*.
    The *sols* is a list of strings representing start solutions.
    The *pars* is a tuple of 14 values for the path parameters.
    On return are the string representations of the solutions
    computed at the end of the paths.
    If *gamma* on input equals zero, then a random complex number is generated,
    otherwise the real and imaginary parts of *gamma* are used.
    """
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_start_system
    from phcpy.phcpy2c3 import py2c_ade_manypaths_dd_pars
    from phcpy.interface import store_dobldobl_system
    from phcpy.interface import store_dobldobl_solutions
    from phcpy.interface import load_dobldobl_solutions
    store_dobldobl_system(target)
    py2c_copy_dobldobl_container_to_target_system()
    store_dobldobl_system(start)
    py2c_copy_dobldobl_container_to_start_system()
    dim = len(start)
    store_dobldobl_solutions(dim, sols)
    if(gamma == 0):
        from random import uniform
        from cmath import exp, pi
        angle = uniform(0, 2*pi)
        gamma = exp(angle*complex(0, 1))
        if(verbose > 0):
            print('random gamma constant :', gamma)
    if(verbose > 0):
        print('the path parameters :\n', pars)
    fail = py2c_ade_manypaths_dd_pars(verbose, gamma.real, gamma.imag, \
        pars[0], pars[1], pars[2], pars[3], pars[4], pars[5], pars[6], \
        pars[7], pars[8], pars[9], pars[10], pars[11], pars[12], pars[13])
    if(fail == 0):
        if(verbose > 0):
            print('Path tracking with AD was a success!')
    else:
        print('Path tracking with AD failed!')
    return load_dobldobl_solutions()

def gpu_double_double_track(target, start, sols, gamma=0, verbose=1):
    r"""
    GPU accelerated path tracking with algorithmic differentiation,
    in double double precision.
    On input are a target system, a start system with solutions.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in sols.
    The *sols* is a list of strings representing start solutions.
    On return are the string representations of the solutions
    computed at the end of the paths.
    If *gamma* on input equals zero, then a random complex number is generated,
    otherwise the real and imaginary parts of *gamma* are used.
    """
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_start_system
    from phcpy.phcpy2c3 import py2c_gpu_manypaths_dd
    from phcpy.interface import store_dobldobl_system
    from phcpy.interface import store_dobldobl_solutions
    from phcpy.interface import load_dobldobl_solutions
    store_dobldobl_system(target)
    py2c_copy_dobldobl_container_to_target_system()
    store_dobldobl_system(start)
    py2c_copy_dobldobl_container_to_start_system()
    dim = len(start)
    store_dobldobl_solutions(dim, sols)
    if(gamma == 0):
        from random import uniform
        from cmath import exp, pi
        angle = uniform(0, 2*pi)
        gamma = exp(angle*complex(0, 1))
        if(verbose > 0):
            print('random gamma constant :', gamma)
    fail = py2c_gpu_manypaths_dd(2, verbose, gamma.real, gamma.imag)
    if(fail == 0):
        if(verbose > 0):
            print('Path tracking on the GPU was a success!')
    else:
        print('Path tracking on the GPU failed!')
    return load_dobldobl_solutions()

def quad_double_track(target, start, sols, gamma=0, pwt=2, tasks=0):
    r"""
    Does path tracking with quad double precision.
    On input are a target system, a start system with solutions,
    optionally a (random) gamma constant and the number of tasks.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system with known solutions in sols.
    The *sols* is a list of strings representing start solutions.
    By default, a random *gamma* constant is generated,
    otherwise *gamma* must be a nonzero complex constant.
    The *pwt* is the power of t in the homotopy.
    Changing the default of *pwt* can only be done if a nonzero complex
    value for *gamma* is provided as well.
    The number of tasks in the multithreading is defined by *tasks*.
    The default zero value for *tasks* indicates no multithreading.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_start_system
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_start_solutions
    from phcpy.phcpy2c3 import py2c_create_quaddobl_homotopy
    from phcpy.phcpy2c3 import py2c_create_quaddobl_homotopy_with_gamma
    from phcpy.phcpy2c3 import py2c_clear_quaddobl_homotopy
    from phcpy.phcpy2c3 import py2c_solve_by_quaddobl_homotopy_continuation
    from phcpy.phcpy2c3 import py2c_clear_quaddobl_operations_data
    from phcpy.phcpy2c3 import py2c_solcon_clear_quaddobl_solutions
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_target_solutions_to_container
    from phcpy.interface import store_quaddobl_system
    from phcpy.interface import store_quaddobl_solutions
    from phcpy.interface import load_quaddobl_solutions
    from phcpy.solver import number_of_symbols
    dim = number_of_symbols(start)
    store_quaddobl_system(target, nbvar=dim)
    py2c_copy_quaddobl_container_to_target_system()
    store_quaddobl_system(start, nbvar=dim)
    py2c_copy_quaddobl_container_to_start_system()
    if(gamma == 0):
        py2c_create_quaddobl_homotopy()
    else:
        py2c_create_quaddobl_homotopy_with_gamma(gamma.real, gamma.imag, pwt)
    store_quaddobl_solutions(dim, sols)
    py2c_copy_quaddobl_container_to_start_solutions()
    py2c_solve_by_quaddobl_homotopy_continuation(tasks)
    py2c_solcon_clear_quaddobl_solutions()
    py2c_copy_quaddobl_target_solutions_to_container()
    py2c_clear_quaddobl_homotopy()
    py2c_clear_quaddobl_operations_data()
    return load_quaddobl_solutions()

def ade_quad_double_track(target, start, sols, gamma=0, verbose=1):
    r"""
    Does path tracking with algorithmic differentiation,
    in quad double precision.
    On input are a target system, a start system with solutions.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in sols.
    The *sols* is a list of strings representing start solutions.
    On return are the string representations of the solutions
    computed at the end of the paths.
    If *gamma* on input equals zero, then a random complex number is generated,
    otherwise the real and imaginary parts of *gamma* are used.
    """
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_start_system
    from phcpy.phcpy2c3 import py2c_ade_manypaths_qd
    from phcpy.interface import store_quaddobl_system
    from phcpy.interface import store_quaddobl_solutions
    from phcpy.interface import load_quaddobl_solutions
    store_quaddobl_system(target)
    py2c_copy_quaddobl_container_to_target_system()
    store_quaddobl_system(start)
    py2c_copy_quaddobl_container_to_start_system()
    dim = len(start)
    store_quaddobl_solutions(dim, sols)
    if(gamma == 0):
        from random import uniform
        from cmath import exp, pi
        angle = uniform(0, 2*pi)
        gamma = exp(angle*complex(0, 1))
        if(verbose > 0):
            print('random gamma constant :', gamma)
    fail = py2c_ade_manypaths_qd(verbose, gamma.real, gamma.imag)
    if(fail == 0):
        if(verbose > 0):
            print('Path tracking with AD was a success!')
    else:
        print('Path tracking with AD failed!')
    return load_quaddobl_solutions()

def ade_tuned_quad_double_track(target, start, sols, pars, gamma=0, verbose=1):
    r"""
    Does path tracking with algorithm differentiation,
    in quad double precision, with tuned path parameters.
    On input are a target system, a start system with solutions.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in *sols*.
    The *sols* is a list of strings representing start solutions.
    The *pars* is a tuple of tuned values for the path parameters.
    On return are the string representations of the solutions
    computed at the end of the paths.
    If *gamma* on input equals zero, then a random complex number is generated,
    otherwise the real and imaginary parts of *gamma* are used.
    """
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_start_system
    from phcpy.phcpy2c3 import py2c_ade_manypaths_qd_pars
    from phcpy.interface import store_quaddobl_system
    from phcpy.interface import store_quaddobl_solutions
    from phcpy.interface import load_quaddobl_solutions
    store_quaddobl_system(target)
    py2c_copy_quaddobl_container_to_target_system()
    store_quaddobl_system(start)
    py2c_copy_quaddobl_container_to_start_system()
    dim = len(start)
    store_quaddobl_solutions(dim, sols)
    if(gamma == 0):
        from random import uniform
        from cmath import exp, pi
        angle = uniform(0, 2*pi)
        gamma = exp(angle*complex(0, 1))
        if(verbose > 0):
            print('random gamma constant :', gamma)
    if(verbose > 0):
        print('the path parameters :\n', pars)
    fail = py2c_ade_manypaths_qd_pars(verbose, gamma.real, gamma.imag, \
        pars[0], pars[1], pars[2], pars[3], pars[4], pars[5], pars[6], \
        pars[7], pars[8], pars[9], pars[10], pars[11], pars[12], pars[13])
    if(fail == 0):
        if(verbose > 0):
            print('Path tracking with AD was a success!')
    else:
        print('Path tracking with AD failed!')
    return load_quaddobl_solutions()

def gpu_quad_double_track(target, start, sols, gamma=0, verbose=1):
    r"""
    GPU accelerated path tracking with algorithmic differentiation,
    in quad double precision.
    On input are a target system, a start system with solutions.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in sols.
    The *sols* is a list of strings representing start solutions.
    On return are the string representations of the solutions
    computed at the end of the paths.
    If *gamma* on input equals zero, then a random complex number is generated,
    otherwise the real and imaginary parts of *gamma* are used.
    """
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_start_system
    from phcpy.phcpy2c3 import py2c_gpu_manypaths_qd
    from phcpy.interface import store_quaddobl_system
    from phcpy.interface import store_quaddobl_solutions
    from phcpy.interface import load_quaddobl_solutions
    store_quaddobl_system(target)
    py2c_copy_quaddobl_container_to_target_system()
    store_quaddobl_system(start)
    py2c_copy_quaddobl_container_to_start_system()
    dim = len(start)
    store_quaddobl_solutions(dim, sols)
    if(gamma == 0):
        from random import uniform
        from cmath import exp, pi
        angle = uniform(0, 2*pi)
        gamma = exp(angle*complex(0, 1))
        if(verbose > 0):
            print('random gamma constant :', gamma)
    fail = py2c_gpu_manypaths_qd(2, verbose, gamma.real, gamma.imag)
    if(fail == 0):
        if(verbose > 0):
            print('Path tracking on the GPU was a success!')
    else:
        print('Path tracking on the GPU failed!')
    return load_quaddobl_solutions()

def multiprecision_track(target, start, sols, gamma=0, pwt=2, decimals=80):
    r"""
    Does path tracking with multiprecision.
    On input are a target system, a start system with solutions,
    and optionally a (random) *gamma* constant.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system with known solutions in sols.
    The *sols* is a list of strings representing start solutions.
    By default, a random *gamma* constant is generated,
    otherwise *gamma* must be a nonzero complex constant.
    The *pwt* is the power of t in the homotopy.
    Changing the default of *pwt* can only be done if a nonzero complex
    value for *gamma* is provided as well.
    The number of decimal places in the working precision is
    given by the value of *decimals*.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy.phcpy2c3 import py2c_copy_multprec_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_multprec_container_to_start_system
    from phcpy.phcpy2c3 import py2c_copy_multprec_container_to_start_solutions
    from phcpy.phcpy2c3 import py2c_create_multprec_homotopy
    from phcpy.phcpy2c3 import py2c_create_multprec_homotopy_with_gamma
    from phcpy.phcpy2c3 import py2c_solve_by_multprec_homotopy_continuation
    from phcpy.phcpy2c3 import py2c_solcon_clear_multprec_solutions
    from phcpy.phcpy2c3 import py2c_copy_multprec_target_solutions_to_container
    from phcpy.interface import store_multprec_system
    from phcpy.interface import store_multprec_solutions
    from phcpy.interface import load_multprec_solutions
    store_multprec_system(target, decimals)
    py2c_copy_multprec_container_to_target_system()
    store_multprec_system(start, decimals)
    py2c_copy_multprec_container_to_start_system()
    # py2c_clear_multprec_homotopy()
    if(gamma == 0):
        py2c_create_multprec_homotopy()
    else:
        py2c_create_multprec_homotopy_with_gamma(gamma.real, gamma.imag, pwt=2)
    dim = len(start)
    store_multprec_solutions(dim, sols)
    py2c_copy_multprec_container_to_start_solutions()
    py2c_solve_by_multprec_homotopy_continuation(decimals)
    py2c_solcon_clear_multprec_solutions()
    py2c_copy_multprec_target_solutions_to_container()
    return load_multprec_solutions()

def track(target, start, sols, \
    precision='d', decimals=80, gamma=0, pwt=2, tasks=0):
    r"""
    Runs the path trackers to track solutions in sols
    at the start system in start to the target system
    in the target list using the current settings of
    the numerical continuation parameters as tuned by
    the function tune_track_parameters() of the tuning module.
    Four levels of precision are supported:

    1. *d*: standard double precision (1.1e-15 or 2^(-53)),

    2. *dd*: double double precision (4.9e-32 or 2^(-104)),

    3. *qd*: quad double precision (1.2e-63 or 2^(-209)).

    4. *mp*: arbitrary multiprecision, with as many decimal places
       in the working precision as the value of *decimals*.

    The next to last parameter is optional.  By default,
    a random complex number will be used for *gamma*,
    otherwise, *gamma* can be any nonzero complex number.
    The *pwt* is the power of t in the homotopy.
    Changing the default of *pwt* can only be done if a nonzero complex
    value for *gamma* is provided as well.
    The last parameter equals the number of *tasks*.
    By default, for *tasks* equal to 0 there is no multitasking.
    For positive values of *tasks*, the multitasking could give
    a speedup of up to the number of tasks, depending how may
    cores are available on the processor.
    """
    if(precision == 'd'):
        return standard_double_track(target, start, sols, gamma, pwt, tasks)
    elif(precision == 'dd'):
        return double_double_track(target, start, sols, gamma, pwt, tasks)
    elif(precision == 'qd'):
        return quad_double_track(target, start, sols, gamma, pwt, tasks)
    elif(precision == 'mp'):
        return multiprecision_track(target, start, sols, gamma, pwt, decimals)
    else:
        print('wrong argument for precision')
        return None

def initialize_standard_tracker(target, start, fixedgamma=True, \
    regamma=0.0, imgamma=0.0):
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
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_start_system
    from phcpy.phcpy2c3 import py2c_initialize_standard_homotopy
    from phcpy.interface import store_standard_system
    store_standard_system(target)
    py2c_copy_standard_container_to_target_system()
    store_standard_system(start)
    py2c_copy_standard_container_to_start_system()
    if fixedgamma:
        return py2c_initialize_standard_homotopy(1, regamma, imgamma)
    else:
        return py2c_initialize_standard_homotopy(0, regamma, imgamma)

def initialize_dobldobl_tracker(target, start, fixedgamma=True, \
    regamma=0.0, imgamma=0.0):
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
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_start_system
    from phcpy.phcpy2c3 import py2c_initialize_dobldobl_homotopy
    from phcpy.interface import store_dobldobl_system
    store_dobldobl_system(target)
    py2c_copy_dobldobl_container_to_target_system()
    store_dobldobl_system(start)
    py2c_copy_dobldobl_container_to_start_system()
    if fixedgamma:
        return py2c_initialize_dobldobl_homotopy(1, regamma, imgamma)
    else:
        return py2c_initialize_dobldobl_homotopy(0, regamma, imgamma)

def initialize_quaddobl_tracker(target, start, fixedgamma=True, \
    regamma=0.0, imgamma=0.0):
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
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_start_system
    from phcpy.phcpy2c3 import py2c_initialize_quaddobl_homotopy
    from phcpy.interface import store_quaddobl_system
    store_quaddobl_system(target)
    py2c_copy_quaddobl_container_to_target_system()
    store_quaddobl_system(start)
    py2c_copy_quaddobl_container_to_start_system()
    if fixedgamma:
        return py2c_initialize_quaddobl_homotopy(1, regamma, imgamma)
    else:
        return py2c_initialize_quaddobl_homotopy(0, regamma, imgamma)

def initialize_multprec_tracker(target, start, fixedgamma=True, decimals=100):
    r"""
    Initializes a path tracker with a generator for a *target*
    and *start* system given in arbitrary multiprecision,
    with the number of decimal places in the working precision
    given by the value of *decimals*.
    If *fixedgamma*, then gamma will be a fixed default value,
    otherwise, a random complex constant for gamma is generated.
    """
    from phcpy.phcpy2c3 import py2c_copy_multprec_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_multprec_container_to_start_system
    from phcpy.phcpy2c3 import py2c_initialize_multprec_homotopy
    from phcpy.interface import store_multprec_system
    store_multprec_system(target, decimals)
    py2c_copy_multprec_container_to_target_system()
    store_multprec_system(start, decimals)
    py2c_copy_multprec_container_to_start_system()
    if fixedgamma:
        return py2c_initialize_multprec_homotopy(1, decimals)
    else:
        return py2c_initialize_multprec_homotopy(0, decimals)

def initialize_varbprec_tracker(target, start, fixedgamma=True):
    r"""
    Initializes a path tracker in variable precision with a *target*
    and *start* system, given as lists of string representations of
    multivariate polynomials.
    If *fixedgamma*, then gamma will be a fixed default value,
    otherwise, a random complex constant for gamma is generated.
    """
    from phcpy.phcpy2c3 import py2c_initialize_varbprec_homotopy
    tgtsys = ''
    for pol in target:
        tgtsys = tgtsys + pol
    nct = len(tgtsys)
    strsys = ''
    for pol in start:
        strsys = strsys + pol
    ncs = len(strsys)
    if fixedgamma:
        return py2c_initialize_varbprec_homotopy(1, nct, tgtsys, ncs, strsys)
    else:
        return py2c_initialize_varbprec_homotopy(0, nct, tgtsys, ncs, strsys)

def initialize_standard_solution(nvar, sol):
    r"""
    A standard double precision path tracker with a generator is
    initialized with a start solution *sol* in a number of
    variables equal to the value of *nvar*.
    """
    from phcpy.phcpy2c3 import py2c_initialize_standard_solution
    from phcpy.interface import store_standard_solutions
    store_standard_solutions(nvar, [sol])
    return py2c_initialize_standard_solution(1)

def initialize_dobldobl_solution(nvar, sol):
    r"""
    A double double precision path tracker with a generator is
    initialized with a start solution *sol* in a number of
    variables equal to the value of *nvar*.
    """
    from phcpy.phcpy2c3 import py2c_initialize_dobldobl_solution
    from phcpy.interface import store_dobldobl_solutions
    store_dobldobl_solutions(nvar, [sol])
    return py2c_initialize_dobldobl_solution(1)

def initialize_quaddobl_solution(nvar, sol):
    r"""
    A quad double precision path tracker with a generator is
    initialized with a start solution *sol* in a number of
    variables equal to the value of *nvar*.
    """
    from phcpy.phcpy2c3 import py2c_initialize_quaddobl_solution
    from phcpy.interface import store_quaddobl_solutions
    store_quaddobl_solutions(nvar, [sol])
    return py2c_initialize_quaddobl_solution(1)

def initialize_multprec_solution(nvar, sol):
    r"""
    A multiprecision path tracker with a generator is
    initialized with a start solution *sol* in a number of
    variables equal to the value of *nvar*.
    """
    from phcpy.phcpy2c3 import py2c_initialize_multprec_solution
    from phcpy.interface import store_multprec_solutions
    store_multprec_solutions(nvar, [sol])
    return py2c_initialize_multprec_solution(1)

def initialize_varbprec_solution(nvar, sol):
    r"""
    A variable precision path tracker with a generator is
    initialized with a start solution *sol* in a number of
    variables equal to the value of *nvar*.
    """
    from phcpy.phcpy2c3 import py2c_initialize_varbprec_solution
    return py2c_initialize_varbprec_solution(nvar, len(sol), sol)

def next_standard_solution():
    r"""
    Returns the next solution on a path tracked with standard
    double precision arithmetic, provided the functions
    **initialize_standard_tracker()** and **initialize_standard_solution()**
    have been executed properly.
    """
    from phcpy.phcpy2c3 import py2c_next_standard_solution
    from phcpy.phcpy2c3 import py2c_solcon_length_standard_solution_string
    from phcpy.phcpy2c3 import py2c_solcon_write_standard_solution_string
    py2c_next_standard_solution(1)
    lns = py2c_solcon_length_standard_solution_string(1)
    sol = py2c_solcon_write_standard_solution_string(1, lns)
    return sol

def next_dobldobl_solution():
    r"""
    Returns the next solution on a path tracked with double
    double precision arithmetic, provided the functions
    **initialize_dobldobl_tracker()** and **initialize_dobldobl_solution()**
    have been executed properly.
    """
    from phcpy.phcpy2c3 import py2c_next_dobldobl_solution
    from phcpy.phcpy2c3 import py2c_solcon_length_dobldobl_solution_string
    from phcpy.phcpy2c3 import py2c_solcon_write_dobldobl_solution_string
    py2c_next_dobldobl_solution(1)
    lns = py2c_solcon_length_dobldobl_solution_string(1)
    sol = py2c_solcon_write_dobldobl_solution_string(1, lns)
    return sol

def next_quaddobl_solution():
    r"""
    Returns the next solution on a path tracked with quad
    double precision arithmetic, provided the functions
    **initialize_quaddobl_tracker()** and **initialize_quaddobl_solution()**
    have been executed properly.
    """
    from phcpy.phcpy2c3 import py2c_next_quaddobl_solution
    from phcpy.phcpy2c3 import py2c_solcon_length_quaddobl_solution_string
    from phcpy.phcpy2c3 import py2c_solcon_write_quaddobl_solution_string
    py2c_next_quaddobl_solution(1)
    lns = py2c_solcon_length_quaddobl_solution_string(1)
    sol = py2c_solcon_write_quaddobl_solution_string(1, lns)
    return sol

def next_multprec_solution():
    r"""
    Returns the next solution on a path tracked with arbitrary
    multiprecision arithmetic, provided the functions
    **initialize_multprec_tracker()** and **initialize_multprec_solution()**
    have been executed properly.
    """
    from phcpy.phcpy2c3 import py2c_next_multprec_solution
    from phcpy.phcpy2c3 import py2c_solcon_length_multprec_solution_string
    from phcpy.phcpy2c3 import py2c_solcon_write_multprec_solution_string
    py2c_next_multprec_solution(1)
    lns = py2c_solcon_length_multprec_solution_string(1)
    sol = py2c_solcon_write_multprec_solution_string(1, lns)
    return sol

def next_varbprec_solution(wanted, maxprec, maxit, verbose):
    r"""
    Returns the next solution on a path tracked with variable
    precision arithmetic, provided the functions
    **initialize_varbprec_tracker()** and **initialize_varbprec_solution()**
    have been executed properly.  The four input parameters are

    1. *wanted*: the number of correct decimal places in the solution,

    2. *maxprec*: upper bound on the number of decimal places in the precision,

    3. *maxit*: maximum number of iterations, and

    4. *verbose*: flag to indicate if intermediate output is wanted.
    """
    from phcpy.phcpy2c3 import py2c_next_varbprec_solution
    sol = py2c_next_varbprec_solution(wanted, maxprec, maxit, verbose)
    return sol

def standard_double_crude_track(target, start, sols, gamma=0, verbose=True):
    r"""
    A crude path tracker does not refine or postprocess the solutions
    at the end of the paths, computed in standard double precision.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in *sols*.
    The *sols* is a list of strings representing start solutions.
    By default, a random *gamma* constant is generated,
    otherwise *gamma* must be a nonzero complex constant.
    If *verbose*, then the solution vectors are written to screen.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_start_system
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_start_solutions
    from phcpy.phcpy2c3 import py2c_create_standard_homotopy
    from phcpy.phcpy2c3 import py2c_create_standard_homotopy_with_gamma
    from phcpy.phcpy2c3 import py2c_standard_crude_tracker
    from phcpy.phcpy2c3 import py2c_solcon_clear_standard_solutions
    from phcpy.interface import store_standard_system
    from phcpy.interface import store_standard_solutions
    from phcpy.interface import load_standard_solutions
    from phcpy.solver import number_of_symbols
    dim = number_of_symbols(start)
    store_standard_system(target, nbvar=dim)
    py2c_copy_standard_container_to_target_system()
    store_standard_system(start, nbvar=dim)
    py2c_copy_standard_container_to_start_system()
    if(gamma == 0):
        py2c_create_standard_homotopy()
    else:
        py2c_create_standard_homotopy_with_gamma(gamma.real, gamma.imag)
    store_standard_solutions(dim, sols)
    py2c_copy_standard_container_to_start_solutions()
    py2c_standard_crude_tracker(int(verbose))
    return load_standard_solutions()

def double_double_crude_track(target, start, sols, gamma=0, verbose=True):
    r"""
    A crude path tracker does not refine or postprocess the solutions
    at the end of the paths, computed in double double precision.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in *sols*.
    The *sols* is a list of strings representing start solutions.
    By default, a random *gamma* constant is generated,
    otherwise *gamma* must be a nonzero complex constant.
    If *verbose*, then the solution vectors are written to screen.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_start_system
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_start_solutions
    from phcpy.phcpy2c3 import py2c_create_dobldobl_homotopy
    from phcpy.phcpy2c3 import py2c_create_dobldobl_homotopy_with_gamma
    from phcpy.phcpy2c3 import py2c_dobldobl_crude_tracker
    from phcpy.phcpy2c3 import py2c_solcon_clear_dobldobl_solutions
    from phcpy.interface import store_dobldobl_system
    from phcpy.interface import store_dobldobl_solutions
    from phcpy.interface import load_dobldobl_solutions
    from phcpy.solver import number_of_symbols
    dim = number_of_symbols(start)
    store_dobldobl_system(target, nbvar=dim)
    py2c_copy_dobldobl_container_to_target_system()
    store_dobldobl_system(start, nbvar=dim)
    py2c_copy_dobldobl_container_to_start_system()
    if(gamma == 0):
        py2c_create_dobldobl_homotopy()
    else:
        py2c_create_dobldobl_homotopy_with_gamma(gamma.real, gamma.imag)
    store_dobldobl_solutions(dim, sols)
    py2c_copy_dobldobl_container_to_start_solutions()
    py2c_dobldobl_crude_tracker(int(verbose))
    return load_dobldobl_solutions()

def quad_double_crude_track(target, start, sols, gamma=0, verbose=True):
    r"""
    A crude path tracker does not refine or postprocess the solutions
    at the end of the paths, computed in quad double precision.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in *sols*.
    The *sols* is a list of strings representing start solutions.
    By default, a random *gamma* constant is generated,
    otherwise *gamma* must be a nonzero complex constant.
    If *verbose*, then the solution vectors are written to screen.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_start_system
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_start_solutions
    from phcpy.phcpy2c3 import py2c_create_quaddobl_homotopy
    from phcpy.phcpy2c3 import py2c_create_quaddobl_homotopy_with_gamma
    from phcpy.phcpy2c3 import py2c_quaddobl_crude_tracker
    from phcpy.phcpy2c3 import py2c_solcon_clear_quaddobl_solutions
    from phcpy.interface import store_quaddobl_system
    from phcpy.interface import store_quaddobl_solutions
    from phcpy.interface import load_quaddobl_solutions
    from phcpy.solver import number_of_symbols
    dim = number_of_symbols(start)
    store_quaddobl_system(target, nbvar=dim)
    py2c_copy_quaddobl_container_to_target_system()
    store_quaddobl_system(start, nbvar=dim)
    py2c_copy_quaddobl_container_to_start_system()
    if(gamma == 0):
        py2c_create_quaddobl_homotopy()
    else:
        py2c_create_quaddobl_homotopy_with_gamma(gamma.real, gamma.imag)
    store_quaddobl_solutions(dim, sols)
    py2c_copy_quaddobl_container_to_start_solutions()
    py2c_quaddobl_crude_tracker(int(verbose))
    return load_quaddobl_solutions()

def test_track(verbose=True, precision='d', decimals=80):
    """
    Tests the path tracking on a small random system.
    Two random trinomials are generated and random constants
    are added to ensure there are no singular solutions
    so we can use this generated system as a start system.
    The target system has the same monomial structure as
    the start system, but with random real coefficients.
    Because all coefficients are random, the number of
    paths tracked equals the mixed volume of the system.
    """
    from phcpy.solver import random_trinomials, real_random_trinomials
    from phcpy.solver import solve, mixed_volume, newton_step
    pols = random_trinomials()
    real_pols = real_random_trinomials(pols)
    from random import uniform as u
    qone = pols[0][:-1] + ('%+.17f' % u(-1, +1)) + ';'
    qtwo = pols[1][:-1] + ('%+.17f' % u(-1, +1)) + ';'
    rone = real_pols[0][:-1] + ('%+.17f' % u(-1, +1)) + ';'
    rtwo = real_pols[1][:-1] + ('%+.17f' % u(-1, +1)) + ';'
    start = [qone, qtwo]
    target = [rone, rtwo]
    start_sols = solve(start, verbose)
    sols = track(target, start, start_sols, precision, decimals)
    mixvol = mixed_volume(target)
    print('mixed volume of the target is', mixvol)
    print('number of solutions found :', len(sols))
    newton_step(target, sols, precision, decimals)
    # for s in sols: print s

def test_next_track(precision='d', decimals=80):
    """
    Tests the step-by-step tracking of a solution path.
    Three levels of precision are supported:
    d  : standard double precision (1.1e-15 or 2^(-53)),
    dd : double double precision (4.9e-32 or 2^(-104)),
    qd : quad double precision (1.2e-63 or 2^(-209)).
    mp : arbitrary multiprecision with as many decimal places
    in the working precision as the value set by decimals.
    """
    from phcpy.solver import total_degree_start_system
    quadrics = ['x**2 + 4*y**2 - 4;', '2*y**2 - x;']
    (startsys, startsols) = total_degree_start_system(quadrics)
    print('the first start solution :\n', startsols[0])
    if(precision == 'd'):
        initialize_standard_tracker(quadrics, startsys)
        initialize_standard_solution(2, startsols[0])
        while True:
            sol = next_standard_solution()
            print('the next solution :\n', sol)
            answer = input('continue ? (y/n) ')
            if(answer != 'y'):
                break
    elif(precision == 'dd'):
        initialize_dobldobl_tracker(quadrics, startsys)
        initialize_dobldobl_solution(2, startsols[0])
        while True:
            sol = next_dobldobl_solution()
            print('the next solution :\n', sol)
            answer = input('continue ? (y/n) ')
            if(answer != 'y'):
                break
    elif(precision == 'qd'):
        initialize_quaddobl_tracker(quadrics, startsys)
        initialize_quaddobl_solution(2, startsols[0])
        while True:
            sol = next_quaddobl_solution()
            print('the next solution :\n', sol)
            answer = input('continue ? (y/n) ')
            if(answer != 'y'):
                break
    elif(precision == 'mp'):
        initialize_multprec_tracker(quadrics, startsys, decimals)
        initialize_multprec_solution(2, startsols[0])
        while True:
            sol = next_multprec_solution()
            print('the next solution :\n', sol)
            answer = input('continue ? (y/n) ')
            if(answer != 'y'):
                break
    else:
        print('wrong argument for precision')

def test_monitored_track():
    """
    Often the number of paths to track can be huge
    and waiting on the outcome of track() without knowing
    how many paths that have been tracked so far can be annoying.
    This script illustrates how one can monitor the progress
    of the path tracking.  We must use the same gamma constant
    with each call of track.
    """
    from random import uniform
    from cmath import pi, exp
    angle = uniform(0, 2*pi)
    ourgamma = exp(complex(0, 1)*angle)
    from phcpy.solver import total_degree_start_system, newton_step
    quadrics = ['x**2 + 4*y**2 - 4;', '2*y**2 - x;']
    (startsys, startsols) = total_degree_start_system(quadrics)
    targetsols = []
    for ind in range(0, len(startsols)):
        print('tracking path', ind+1, '...', end=' ')
        endsol = track(quadrics, startsys, [startsols[ind]], gamma=ourgamma)
        print('found solution\n', endsol[0])
        targetsols.append(endsol[0])
    print('tracked', len(targetsols), 'paths, running newton_step...')
    newton_step(quadrics, targetsols)

def cyclic3homotopy():
    """
    Returns a triplet with the cyclic 3-roots system as first item,
    then in second place a random coefficient system, and the start
    solutions in the third item of the triplet.
    """
    c3 = ['x0 + x1 + x2;', 'x0*x1 + x1*x2 + x2*x0;', 'x0*x1*x2 - 1;']
    c3q = ['+( 9.32743666318680E-01 - 3.60540223750954E-01*i)*x0' \
        + ' +(-6.60494034825272E-01 + 7.50831292608554E-01*i)*x1' \
        + ' +(-6.37456882105859E-01 - 7.70486030668874E-01*i)*x2;', \
          ' +(-6.65114940988527E-01 + 7.46740996111656E-01*i)*x0*x1' \
        + ' +( 6.01027652324430E-01 + 7.99228228443780E-01*i)*x0*x2' \
        + ' +( 9.99999810085609E-01 + 6.16302478459703E-04*i)*x1*x2;', \
          ' +(-5.95090196510060E-01 - 8.03658918956057E-01*i)*x0*x1*x2' \
        + ' +(-4.08655719716833E-01 + 9.12688612146945E-01*i);']
    c3qs0 = 't :  1.00000000000000E+00   0.00000000000000E+00\n' \
          + 'm : 1\n' \
          + 'the solution for t :\n' \
          + ' x0 :  1.15429240714256E-01   6.09909582313413E-01\n' \
          + ' x1 :  7.25702143828498E-01  -7.09583364893093E-01\n' \
          + ' x2 :  1.43006908284017E+00   6.88642063770792E-01\n' \
          + '== err :  2.203E-16 = rco :  1.557E-01 = res :  2.220E-16 =='
    c3qs1 = ' t :  1.00000000000000E+00   0.00000000000000E+00\n' \
          + 'm : 1\n' \
          + 'the solution for t :\n' \
          + ' x0 : -5.85911812652100E-01  -2.04990136358612E-01\n' \
          + ' x1 :  2.51666148186012E-01   9.83268174582854E-01\n' \
          + ' x2 : -1.31141606276013E+00   8.94155123020900E-01\n' \
          + '== err :  2.232E-16 = rco :  1.540E-01 = res :  6.661E-16 =='
    c3qs2 = 't :  1.00000000000000E+00   0.00000000000000E+00\n' \
          + 'm : 1\n' \
          + 'the solution for t :\n' \
          + ' x0 :  4.70482571937844E-01  -4.04919445954801E-01\n' \
          + ' x1 : -9.77368292014510E-01  -2.73684809689761E-01\n' \
          + ' x2 : -1.18653020080034E-01  -1.58279718679169E+00\n' \
          + '== err :  1.721E-16 = rco :  1.506E-01 = res :  3.331E-16 =='
    c3qs3 = 't :  1.00000000000000E+00   0.00000000000000E+00\n' \
          + 'm : 1\n' \
          + 'the solution for t :\n' \
          + ' x0 :  1.33160482410366E+00   9.06706762952865E-01\n' \
          + ' x1 :  5.15653962776839E-01   8.39542362953949E-01\n' \
          + ' x2 :  5.33980349490553E-01  -3.34360029828785E-01\n' \
          + '== err :  1.040E-16 = rco :  2.696E-01 = res :  2.220E-16 =='
    c3qs4 = 't :  1.00000000000000E+00   0.00000000000000E+00\n' \
          + 'm : 1\n' \
          + 'the solution for t :\n' \
          + ' x0 : -1.45103350255217E+00   6.99850223999248E-01\n' \
          + ' x1 : -9.84891995259755E-01   2.67982498498837E-02\n' \
          + ' x2 :  2.25741050965737E-02   6.29620562694904E-01\n' \
          + '== err :  2.796E-16 = rco :  3.004E-01 = res :  4.441E-16 =='
    c3qs5 = 't :  1.00000000000000E+00   0.00000000000000E+00\n' \
          + 'm : 1\n' \
          + 'the solution for t :\n' \
          + ' x0 :  1.19428678448505E-01  -1.60655698695211E+00\n' \
          + ' x1 :  4.69238032482916E-01  -8.66340612803833E-01\n' \
          + ' x2 : -5.56554454587126E-01  -2.95260532866119E-01\n' \
          + '== err :  2.131E-16 = rco :  2.195E-01 = res :  5.551E-17 =='
    c3qsols = [c3qs0, c3qs1, c3qs2, c3qs3, c3qs4, c3qs5]
    return (c3, c3q, c3qsols)

def test_ade_double_track():
    """
    Tests the path tracker on the cyclic 3-roots problem,
    in standard double precision.
    """
    (c3, c3q, c3qsols) = cyclic3homotopy()
    ans = input('Tune the path parameters ? (y/n) ')
    if(ans != 'y'):
        sols = ade_double_track(c3, c3q, c3qsols)
    else:
        from phcpy.tuning import tune_path_parameters as tune
        pars = tune(16)
        sols = ade_tuned_double_track(c3, c3q, c3qsols, pars)
    for sol in sols:
        print(sol)

def test_ade_double_double_track():
    """
    Tests the path tracker on the cyclic 3-roots problem,
    in double double precision.
    """
    (c3, c3q, c3qsols) = cyclic3homotopy()
    ans = input('Tune the path parameters ? (y/n) ')
    if(ans != 'y'):
        sols = ade_double_double_track(c3, c3q, c3qsols)
    else:
        from phcpy.tuning import tune_path_parameters as tune
        pars = tune(32)
        sols = ade_tuned_double_double_track(c3, c3q, c3qsols, pars)
    for sol in sols:
        print(sol)

def test_ade_quad_double_track():
    """
    Tests the path tracker on the cyclic 3-roots problem,
    in quad double precision.
    """
    (c3, c3q, c3qsols) = cyclic3homotopy()
    ans = input('Tune the path parameters ? (y/n) ')
    if(ans != 'y'):
        sols = ade_quad_double_track(c3, c3q, c3qsols)
    else:
        from phcpy.tuning import tune_path_parameters as tune
        pars = tune(64)
        sols = ade_tuned_quad_double_track(c3, c3q, c3qsols, pars)
    for sol in sols:
        print(sol)

def test_crude_tracker(precision='d'):
    """
    Runs the crude path trackers.  Values for precision are 'd', 'dd', or
    'qd', respectively for double, double double, or quad double precision.
    """
    from phcpy.families import noon
    from phcpy.solver import random_linear_product_system as rlps
    target = noon(5)
    (start, startsols) = rlps(target)
    if(precision == 'd'):
        sols = standard_double_crude_track(target, start, startsols)
    elif(precision == 'dd'):
        sols = double_double_crude_track(target, start, startsols)
    elif(precision == 'qd'):
        sols = quad_double_crude_track(target, start, startsols)
    else:
        print('Wrong value for the precision.')
    print('Number of solutions returned :', len(sols))

def test():
    """
    Runs test_track(), test_next_track(), and test_monitored_track().
    """
    print('\ntesting path tracker...\n')
    test_track() # precision='mp', decimals=48)
    print('\ntesting step-by-step tracking...\n')
    test_next_track() # (precision='mp', decimals=48)
    print('\ntesting monitored tracking...\n')
    test_monitored_track()

if __name__ == "__main__":
    test()
