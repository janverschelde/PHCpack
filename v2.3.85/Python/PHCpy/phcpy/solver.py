"""
The main functionality of PHCpack is its blackbox solver
and the path tracking routines, respectively exported by
the functions solve and track.  An important task of the
blackbox solver is the mixed-volume computation, available
in the function mixed_volume.
"""

def random_trinomials():
    """
    Returns a system of two trinomials equations for testing.
    A trinomial consists of three monomials in two variables.
    Exponents are uniform between 0 and 5 and coefficients are
    on the complex unit circle.
    """
    from random import randint as r
    exponents = [(r(0, 5), r(0, 5)) for i in range(0, 6)]
    monomials = map(lambda e: 'x^%d*y^%d' % e, exponents)
    from random import uniform as u
    from math import cos, sin, pi
    angles = [u(0, 2*pi) for i in range(0, 6)]
    cff = map(lambda a: '(' + str(cos(a)) + '%+.14f' % sin(a) + '*i)', angles)
    one = '+'.join(cff[i] + '*' + monomials[i] for i in range(0, 3)) + ';'
    two = '+'.join(cff[i] + '*' + monomials[i] for i in range(3, 6)) + ';'
    return [one, two]

def real_random_trinomials(sys):
    """
    On input in sys are two random trinonials with complex coefficients,
    in the format what random_trinomials() returns.
    On return is a list of two real random trinomials with the same
    monomial structure but with random real coefficients in [-1,+1].
    """
    from random import uniform as u
    result = []
    for pol in sys:
        terms = pol.split(')')
        rpol = ''
        for i in range(1, len(terms)-1):
            rterm = terms[i].split('+')
            cff = '%+.17f' % u(-1, +1)
            rpol = rpol + cff + rterm[0]
        cff = '%+.17f' % u(-1, +1)
        rpol = rpol + cff + terms[len(terms)-1]
        result.append(rpol)
    return result

def store_system(polsys):
    """
    Stores the polynomials represented by the list of
    strings in polsys into the system container.
    """
    from phcpy2c import py2c_syscon_clear_system
    from phcpy2c import py2c_syscon_initialize_number
    from phcpy2c import py2c_syscon_store_polynomial
    py2c_syscon_clear_system()
    dim = len(polsys)
    py2c_syscon_initialize_number(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        py2c_syscon_store_polynomial(nchar, dim, cnt+1, pol)

def store_dobldobl_system(polsys):
    """
    Stores the polynomials represented by the list of strings in polsys
    into the systems container for double double arithmetic.
    """
    from phcpy2c import py2c_syscon_clear_dobldobl_system
    from phcpy2c\
    import py2c_syscon_initialize_number_of_dobldobl_polynomials
    from phcpy2c import py2c_syscon_store_dobldobl_polynomial
    py2c_syscon_clear_dobldobl_system()
    dim = len(polsys)
    py2c_syscon_initialize_number_of_dobldobl_polynomials(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        py2c_syscon_store_dobldobl_polynomial(nchar, dim, cnt+1, pol)

def store_quaddobl_system(polsys):
    """
    Stores the polynomials represented by the list of strings in polsys
    into the systems container for quad double arithmetic.
    """
    from phcpy2c import py2c_syscon_clear_quaddobl_system
    from phcpy2c\
    import py2c_syscon_initialize_number_of_quaddobl_polynomials
    from phcpy2c import py2c_syscon_store_quaddobl_polynomial
    py2c_syscon_clear_quaddobl_system()
    dim = len(polsys)
    py2c_syscon_initialize_number_of_quaddobl_polynomials(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        py2c_syscon_store_quaddobl_polynomial(nchar, dim, cnt+1, pol)

def store_multprec_system(polsys, decimals):
    """
    Stores the polynomials represented by the list of strings in polsys
    into the systems container for multiprecision arithmetic.
    The parameter decimals equals the number of decimal places
    in the working precision for the parsing of the strings in polsys.
    """
    from phcpy2c import py2c_syscon_clear_multprec_system
    from phcpy2c\
    import py2c_syscon_initialize_number_of_multprec_polynomials
    from phcpy2c import py2c_syscon_store_multprec_polynomial
    py2c_syscon_clear_multprec_system()
    dim = len(polsys)
    py2c_syscon_initialize_number_of_multprec_polynomials(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        py2c_syscon_store_multprec_polynomial(nchar, dim, cnt+1, decimals, pol)

def load_system():
    """
    Returns the polynomials stored in the system container.
    """
    from phcpy2c import py2c_syscon_number_of_polynomials
    from phcpy2c import py2c_syscon_load_polynomial
    dim = py2c_syscon_number_of_polynomials()
    result = []
    for ind in range(1, dim+1):
        result.append(py2c_syscon_load_polynomial(ind))
    return result

def load_dobldobl_system():
    """
    Returns the polynomials stored in the system container
    with double double complex coefficients.
    """
    from phcpy2c import py2c_syscon_number_of_dobldobl_polynomials
    from phcpy2c import py2c_syscon_load_dobldobl_polynomial
    dim = py2c_syscon_number_of_dobldobl_polynomials()
    result = []
    for ind in range(1, dim+1):
        result.append(py2c_syscon_load_dobldobl_polynomial(ind))
    return result

def load_quaddobl_system():
    """
    Returns the polynomials stored in the system container
    with quad double complex coefficients.
    """
    from phcpy2c import py2c_syscon_number_of_quaddobl_polynomials
    from phcpy2c import py2c_syscon_load_quaddobl_polynomial
    dim = py2c_syscon_number_of_quaddobl_polynomials()
    result = []
    for ind in range(1, dim+1):
        result.append(py2c_syscon_load_quaddobl_polynomial(ind))
    return result

def load_multprec_system():
    """
    Returns the polynomials stored in the system container
    with arbitrary multiprecision complex coefficients.
    """
    from phcpy2c import py2c_syscon_number_of_multprec_polynomials
    from phcpy2c import py2c_syscon_load_multprec_polynomial
    dim = py2c_syscon_number_of_multprec_polynomials()
    result = []
    for ind in range(1, dim+1):
        result.append(py2c_syscon_load_multprec_polynomial(ind))
    return result

def store_solutions(nvar, sols):
    """
    Stores the solutions in the list sols, represented as strings
    in PHCpack format into the solution container.
    The number nvar equals the number of variables.
    """
    from phcpy2c import py2c_solcon_clear_solutions
    from phcpy2c import py2c_solcon_append_solution_string
    py2c_solcon_clear_solutions()
    for ind in range(0, len(sols)):
        py2c_solcon_append_solution_string(nvar, len(sols[ind]), sols[ind])

def store_dobldobl_solutions(nvar, sols):
    """
    Stores the solutions in the list sols, represented as strings
    in PHCpack format into the solution container for processing
    with complex double double arithmetic.
    The number nvar equals the number of variables.
    """
    from phcpy2c import py2c_solcon_clear_dobldobl_solutions
    from phcpy2c import py2c_solcon_append_dobldobl_solution_string
    py2c_solcon_clear_dobldobl_solutions()
    for ind in range(0, len(sols)):
        py2c_solcon_append_dobldobl_solution_string\
        (nvar, len(sols[ind]), sols[ind])

def store_quaddobl_solutions(nvar, sols):
    """
    Stores the solutions in the list sols, represented as strings
    in PHCpack format into the solution container for processing
    with complex quad double arithmetic.
    The number n equals the number of variables.
    """
    from phcpy2c import py2c_solcon_clear_quaddobl_solutions
    from phcpy2c import py2c_solcon_append_quaddobl_solution_string
    py2c_solcon_clear_quaddobl_solutions()
    for ind in range(0, len(sols)):
        py2c_solcon_append_quaddobl_solution_string\
        (nvar, len(sols[ind]), sols[ind])

def store_multprec_solutions(nvar, sols):
    """
    Stores the solutions in the list sols, represented as strings
    in PHCpack format into the solution container for processing
    with complex multiprecision arithmetic.
    The number n equals the number of variables.
    """
    from phcpy2c import py2c_solcon_clear_multprec_solutions
    from phcpy2c import py2c_solcon_append_multprec_solution_string
    py2c_solcon_clear_multprec_solutions()
    for ind in range(0, len(sols)):
        py2c_solcon_append_multprec_solution_string\
        (nvar, len(sols[ind]), sols[ind])

def load_solutions():
    """
    Returns the list of solutions stored in the solutions container.
    """
    from phcpy2c import py2c_solcon_number_of_solutions
    from phcpy2c import py2c_solcon_length_solution_string
    from phcpy2c import py2c_solcon_write_solution_string
    nbsols = py2c_solcon_number_of_solutions()
    result = []
    for ind in range(1, nbsols+1):
        lns = py2c_solcon_length_solution_string(ind)
        sol = py2c_solcon_write_solution_string(ind, lns)
        result.append(sol)
    return result

def load_dobldobl_solutions():
    """
    Returns the list of solutions stored in the container
    for complex double double solutions.
    """
    from phcpy2c import py2c_solcon_number_of_dobldobl_solutions
    from phcpy2c import py2c_solcon_length_dobldobl_solution_string
    from phcpy2c import py2c_solcon_write_dobldobl_solution_string
    nbsols = py2c_solcon_number_of_dobldobl_solutions()
    result = []
    for ind in range(1, nbsols+1):
        lns = py2c_solcon_length_dobldobl_solution_string(ind)
        sol = py2c_solcon_write_dobldobl_solution_string(ind, lns)
        result.append(sol)
    return result

def load_quaddobl_solutions():
    """
    Returns the list of solutions stored in the container
    for complex quad double solutions.
    """
    from phcpy2c import py2c_solcon_number_of_quaddobl_solutions
    from phcpy2c import py2c_solcon_length_quaddobl_solution_string
    from phcpy2c import py2c_solcon_write_quaddobl_solution_string
    nbsols = py2c_solcon_number_of_quaddobl_solutions()
    result = []
    for ind in range(1, nbsols+1):
        lns = py2c_solcon_length_quaddobl_solution_string(ind)
        sol = py2c_solcon_write_quaddobl_solution_string(ind, lns)
        result.append(sol)
    return result

def load_multprec_solutions():
    """
    Returns the list of solutions stored in the container
    for complex multiprecision solutions.
    """
    from phcpy2c import py2c_solcon_number_of_multprec_solutions
    from phcpy2c import py2c_solcon_length_multprec_solution_string
    from phcpy2c import py2c_solcon_write_multprec_solution_string
    nbsols = py2c_solcon_number_of_multprec_solutions()
    result = []
    for ind in range(1, nbsols+1):
        lns = py2c_solcon_length_multprec_solution_string(ind)
        sol = py2c_solcon_write_multprec_solution_string(ind, lns)
        result.append(sol)
    return result

def solve(pols, silent=False):
    """
    Calls the blackbox solver of PHCpack.
    On input in pols is a list of strings.
    By default, the solver will print to screen the
    computed root counts.  To make the solver silent,
    set the flag silent to True.
    """
    from phcpy2c import py2c_syscon_clear_Laurent_system
    from phcpy2c import py2c_syscon_initialize_number_of_Laurentials
    from phcpy2c import py2c_syscon_store_Laurential
    from phcpy2c import py2c_solcon_clear_solutions
    from phcpy2c import py2c_solve_Laurent_system
    py2c_syscon_clear_Laurent_system()
    py2c_solcon_clear_solutions()
    dim = len(pols)
    py2c_syscon_initialize_number_of_Laurentials(dim)
    for ind in range(0, dim):
        pol = pols[ind]
        nchar = len(pol)
        py2c_syscon_store_Laurential(nchar, dim, ind+1, pol)
    py2c_solve_Laurent_system(silent)
    return load_solutions()

def tune_track_parameters(difficulty=0, digits=16, \
                          interactive=False, silent=True):
    """
    Tunes the numerical continuation parameters for use
    during path tracking based on two parameters:
    the difficulty of the solution path (difficulty)
    and the number of decimal places (digits) for the
    accuracy of the approximations along a path.
    Increasing difficulty will decrease the step size
    and increase the number of steps along a path.
    Increasing digits will decrease the tolerances on
    the numerical approximations.
    If interactive is True, then the user can interactively
    adjust specific numerical continuation parameters.
    If silent is False, then the current values of the
    numerical continuation parameters are shown.
    """
    from phcpy2c import py2c_autotune_continuation_parameters
    from phcpy2c import py2c_tune_continuation_parameters
    from phcpy2c import py2c_show_continuation_parameters
    py2c_autotune_continuation_parameters(difficulty, digits)
    if interactive:
        py2c_tune_continuation_parameters()
    elif not silent:
        py2c_show_continuation_parameters()

def standard_double_track(target, start, sols, gamma=0):
    """
    Does path tracking with standard double precision.
    On input are a target system, a start system with solutions,
    and optionally a (random) gamma constant.
    The target is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The start is a list of strings representing the polynomials
    of the start system, with known solutions in sols.
    The sols is a list of strings representing start solutions.
    By default, a random gamma constant is generated,
    otherwise gamma must be a nonzero complex constant.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy2c import py2c_copy_container_to_target_system
    from phcpy2c import py2c_copy_container_to_start_system
    from phcpy2c import py2c_copy_container_to_start_solutions
    from phcpy2c import py2c_create_homotopy
    from phcpy2c import py2c_create_homotopy_with_gamma
    from phcpy2c import py2c_solve_by_homotopy_continuation
    from phcpy2c import py2c_solcon_clear_solutions
    from phcpy2c import py2c_copy_target_solutions_to_container
    store_system(target)
    py2c_copy_container_to_target_system()
    store_system(start)
    py2c_copy_container_to_start_system()
    # py2c_clear_homotopy()
    if(gamma == 0):
        py2c_create_homotopy()
    else:
        py2c_create_homotopy_with_gamma(gamma.real, gamma.imag)
    dim = len(start)
    store_solutions(dim, sols)
    py2c_copy_container_to_start_solutions()
    py2c_solve_by_homotopy_continuation()
    py2c_solcon_clear_solutions()
    py2c_copy_target_solutions_to_container()
    return load_solutions()

def double_double_track(target, start, sols, gamma=0):
    """
    Does path tracking with double double precision.
    On input are a target system, a start system with solutions,
    and optionally a (random) gamma constant.
    The target is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The start is a list of strings representing the polynomials
    of the start system with known solutions in sols.
    The sols is a list of strings representing start solutions.
    By default, a random gamma constant is generated,
    otherwise gamma must be a nonzero complex constant.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy2c import py2c_copy_dobldobl_container_to_target_system
    from phcpy2c import py2c_copy_dobldobl_container_to_start_system
    from phcpy2c import py2c_copy_dobldobl_container_to_start_solutions
    from phcpy2c import py2c_create_dobldobl_homotopy
    from phcpy2c import py2c_create_dobldobl_homotopy_with_gamma
    from phcpy2c import py2c_solve_by_dobldobl_homotopy_continuation
    from phcpy2c import py2c_solcon_clear_dobldobl_solutions
    from phcpy2c import py2c_copy_dobldobl_target_solutions_to_container
    store_dobldobl_system(target)
    py2c_copy_dobldobl_container_to_target_system()
    store_dobldobl_system(start)
    py2c_copy_dobldobl_container_to_start_system()
    # py2c_clear_homotopy()
    if(gamma == 0):
        py2c_create_dobldobl_homotopy()
    else:
        py2c_create_dobldobl_homotopy_with_gamma(gamma.real, gamma.imag)
    dim = len(start)
    store_dobldobl_solutions(dim, sols)
    py2c_copy_dobldobl_container_to_start_solutions()
    py2c_solve_by_dobldobl_homotopy_continuation()
    py2c_solcon_clear_dobldobl_solutions()
    py2c_copy_dobldobl_target_solutions_to_container()
    return load_dobldobl_solutions()

def quad_double_track(target, start, sols, gamma=0):
    """
    Does path tracking with quad double precision.
    On input are a target system, a start system with solutions,
    and optionally a (random) gamma constant.
    The target is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The start is a list of strings representing the polynomials
    of the start system with known solutions in sols.
    The sols is a list of strings representing start solutions.
    By default, a random gamma constant is generated,
    otherwise gamma must be a nonzero complex constant.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy2c import py2c_copy_quaddobl_container_to_target_system
    from phcpy2c import py2c_copy_quaddobl_container_to_start_system
    from phcpy2c import py2c_copy_quaddobl_container_to_start_solutions
    from phcpy2c import py2c_create_quaddobl_homotopy
    from phcpy2c import py2c_create_quaddobl_homotopy_with_gamma
    from phcpy2c import py2c_solve_by_quaddobl_homotopy_continuation
    from phcpy2c import py2c_solcon_clear_quaddobl_solutions
    from phcpy2c import py2c_copy_quaddobl_target_solutions_to_container
    store_quaddobl_system(target)
    py2c_copy_quaddobl_container_to_target_system()
    store_quaddobl_system(start)
    py2c_copy_quaddobl_container_to_start_system()
    # py2c_clear_homotopy()
    if(gamma == 0):
        py2c_create_quaddobl_homotopy()
    else:
        py2c_create_quaddobl_homotopy_with_gamma(gamma.real, gamma.imag)
    dim = len(start)
    store_quaddobl_solutions(dim, sols)
    py2c_copy_quaddobl_container_to_start_solutions()
    py2c_solve_by_quaddobl_homotopy_continuation()
    py2c_solcon_clear_quaddobl_solutions()
    py2c_copy_quaddobl_target_solutions_to_container()
    return load_quaddobl_solutions()

def multiprecision_track(target, start, sols, gamma=0, decimals=80):
    """
    Does path tracking with multiprecision.
    On input are a target system, a start system with solutions,
    and optionally a (random) gamma constant.
    The target is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The start is a list of strings representing the polynomials
    of the start system with known solutions in sols.
    The sols is a list of strings representing start solutions.
    By default, a random gamma constant is generated,
    otherwise gamma must be a nonzero complex constant.
    The number of decimal places in the working precision is
    given by the value of decimals.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy2c import py2c_copy_multprec_container_to_target_system
    from phcpy2c import py2c_copy_multprec_container_to_start_system
    from phcpy2c import py2c_copy_multprec_container_to_start_solutions
    from phcpy2c import py2c_create_multprec_homotopy
    from phcpy2c import py2c_create_multprec_homotopy_with_gamma
    from phcpy2c import py2c_solve_by_multprec_homotopy_continuation
    from phcpy2c import py2c_solcon_clear_multprec_solutions
    from phcpy2c import py2c_copy_multprec_target_solutions_to_container
    store_multprec_system(target, decimals)
    py2c_copy_multprec_container_to_target_system()
    store_multprec_system(start, decimals)
    py2c_copy_multprec_container_to_start_system()
    # py2c_clear_multprec_homotopy()
    if(gamma == 0):
        py2c_create_multprec_homotopy()
    else:
        py2c_create_multprec_homotopy_with_gamma(gamma.real, gamma.imag)
    dim = len(start)
    store_multprec_solutions(dim, sols)
    py2c_copy_multprec_container_to_start_solutions()
    py2c_solve_by_multprec_homotopy_continuation(decimals)
    py2c_solcon_clear_multprec_solutions()
    py2c_copy_multprec_target_solutions_to_container()
    return load_multprec_solutions()

def track(target, start, sols, precision='d', gamma=0):
    """
    Runs the path trackers to track solutions in sols
    at the start system in start to the target system
    in the target list using the current settings of
    the numerical continuation parameters as tuned by
    the function tune_track_parameters().
    Three levels of precision are supported:
    d  : standard double precision (1.1e-15 or 2^(-53)),
    dd : double double precision (4.9e-32 or 2^(-104)),
    qd : quad double precision (1.2e-63 or 2^(-209)).
    The last parameter is optional.  By default,
    a random complex number will be used for gamma,
    otherwise, gamma can be any nonzero complex number.
    """
    if(precision == 'd'):
        return standard_double_track(target, start, sols, gamma)
    elif(precision == 'dd'):
        return double_double_track(target, start, sols, gamma)
    elif(precision == 'qd'):
        return quad_double_track(target, start, sols, gamma)
    else:
        print 'wrong argument for precision'
        return None

def initialize_standard_tracker(target, start):
    """
    Initializes a path tracker with a generator for a target
    and start system given in standard double precision.
    """
    from phcpy2c import py2c_copy_container_to_target_system
    from phcpy2c import py2c_copy_container_to_start_system
    from phcpy2c import py2c_initialize_standard_homotopy
    store_system(target)
    py2c_copy_container_to_target_system()
    store_system(start)
    py2c_copy_container_to_start_system()
    py2c_initialize_standard_homotopy()

def initialize_dobldobl_tracker(target, start):
    """
    Initializes a path tracker with a generator for a target
    and start system given in double double precision.
    """
    from phcpy2c import py2c_copy_dobldobl_container_to_target_system
    from phcpy2c import py2c_copy_dobldobl_container_to_start_system
    from phcpy2c import py2c_initialize_dobldobl_homotopy
    store_dobldobl_system(target)
    py2c_copy_dobldobl_container_to_target_system()
    store_dobldobl_system(start)
    py2c_copy_dobldobl_container_to_start_system()
    py2c_initialize_dobldobl_homotopy()

def initialize_quaddobl_tracker(target, start):
    """
    Initializes a path tracker with a generator for a target
    and start system given in quad double precision.
    """
    from phcpy2c import py2c_copy_quaddobl_container_to_target_system
    from phcpy2c import py2c_copy_quaddobl_container_to_start_system
    from phcpy2c import py2c_initialize_quaddobl_homotopy
    store_quaddobl_system(target)
    py2c_copy_quaddobl_container_to_target_system()
    store_quaddobl_system(start)
    py2c_copy_quaddobl_container_to_start_system()
    py2c_initialize_quaddobl_homotopy()

def initialize_multprec_tracker(target, start, decimals=100):
    """
    Initializes a path tracker with a generator for a target
    and start system given in arbitrary multiprecision, with
    the number of decimal places in the working precision
    given by the value of decimals.
    """
    from phcpy2c import py2c_copy_multprec_container_to_target_system
    from phcpy2c import py2c_copy_multprec_container_to_start_system
    from phcpy2c import py2c_initialize_multprec_homotopy
    store_multprec_system(target, decimals)
    py2c_copy_multprec_container_to_target_system()
    store_multprec_system(start, decimals)
    py2c_copy_multprec_container_to_start_system()
    py2c_initialize_multprec_homotopy(decimals)

def initialize_standard_solution(nvar, sol):
    """
    Initializes a path tracker with a generator
    for a start solution sol given in standard double precision.
    The value of nvar must equal the number of variables in the
    solution sol and sol is a PHCpack solution string.
    """
    from phcpy2c import py2c_initialize_standard_solution
    store_solutions(nvar, [sol])
    py2c_initialize_standard_solution(1)

def initialize_dobldobl_solution(nvar, sol):
    """
    Initializes a path tracker with a generator
    for a start solution sol given in double double precision.
    The value of nvar must equal the number of variables in the
    solution sol and sol is a PHCpack solution string.
    """
    from phcpy2c import py2c_initialize_dobldobl_solution
    store_dobldobl_solutions(nvar, [sol])
    py2c_initialize_dobldobl_solution(1)

def initialize_quaddobl_solution(nvar, sol):
    """
    Initializes a path tracker with a generator
    for a start solution sol given in quad double precision.
    The value of nvar must equal the number of variables in the
    solution sol and sol is a PHCpack solution string.
    """
    from phcpy2c import py2c_initialize_quaddobl_solution
    store_quaddobl_solutions(nvar, [sol])
    py2c_initialize_quaddobl_solution(1)

def initialize_multprec_solution(nvar, sol):
    """
    Initializes a path tracker with a generator
    for a start solution sol given in arbitrary multiprecision.
    The value of nvar must equal the number of variables in the
    solution sol and sol is a PHCpack solution string.
    """
    from phcpy2c import py2c_initialize_multprec_solution
    store_multprec_solutions(nvar, [sol])
    py2c_initialize_multprec_solution(1)

def next_standard_solution():
    """
    Returns the next solution on a path tracked with standard
    double precision arithmetic, provided the functions
    initialize_standard_tracker() and initialize_standard_solution()
    have been executed properly.
    """
    from phcpy2c import py2c_next_standard_solution
    from phcpy2c import py2c_solcon_length_solution_string
    from phcpy2c import py2c_solcon_write_solution_string
    py2c_next_standard_solution(1)
    lns = py2c_solcon_length_solution_string(1)
    sol = py2c_solcon_write_solution_string(1, lns)
    return sol

def next_dobldobl_solution():
    """
    Returns the next solution on a path tracked with double
    double precision arithmetic, provided the functions
    initialize_dobldobl_tracker() and initialize_dobldobl_solution()
    have been executed properly.
    """
    from phcpy2c import py2c_next_dobldobl_solution
    from phcpy2c import py2c_solcon_length_dobldobl_solution_string
    from phcpy2c import py2c_solcon_write_dobldobl_solution_string
    py2c_next_dobldobl_solution(1)
    lns = py2c_solcon_length_dobldobl_solution_string(1)
    sol = py2c_solcon_write_dobldobl_solution_string(1, lns)
    return sol

def next_quaddobl_solution():
    """
    Returns the next solution on a path tracked with quad
    double precision arithmetic, provided the functions
    initialize_quaddobl_tracker() and initialize_quaddobl_solution()
    have been executed properly.
    """
    from phcpy2c import py2c_next_quaddobl_solution
    from phcpy2c import py2c_solcon_length_quaddobl_solution_string
    from phcpy2c import py2c_solcon_write_quaddobl_solution_string
    py2c_next_quaddobl_solution(1)
    lns = py2c_solcon_length_quaddobl_solution_string(1)
    sol = py2c_solcon_write_quaddobl_solution_string(1, lns)
    return sol

def next_multprec_solution():
    """
    Returns the next solution on a path tracked with arbitrary
    multiprecision arithmetic, provided the functions
    initialize_multprec_tracker() and initialize_multprec_solution()
    have been executed properly.
    """
    from phcpy2c import py2c_next_multprec_solution
    from phcpy2c import py2c_solcon_length_multprec_solution_string
    from phcpy2c import py2c_solcon_write_multprec_solution_string
    py2c_next_multprec_solution(1)
    lns = py2c_solcon_length_multprec_solution_string(1)
    sol = py2c_solcon_write_multprec_solution_string(1, lns)
    return sol

def newton_step(system, solutions, precision='d', decimals=100):
    """
    Applies one Newton step to the solutions of the system.
    For each solution, prints its last line of diagnostics.
    Three levels of precision are supported:
    d  : standard double precision (1.1e-15 or 2^(-53)),
    dd : double double precision (4.9e-32 or 2^(-104)),
    qd : quad double precision (1.2e-63 or 2^(-209)).
    mp : arbitrary precision, where the number of decimal places
    in the working precision is determined by decimals.
    """
    if(precision == 'd'):
        store_system(system)
        store_solutions(len(system), solutions)
        from phcpy2c import py2c_Newton_step
        py2c_Newton_step()
        result = load_solutions()
    elif(precision == 'dd'):
        store_dobldobl_system(system)
        store_dobldobl_solutions(len(system), solutions)
        from phcpy2c import py2c_dobldobl_Newton_step
        py2c_dobldobl_Newton_step()
        result = load_dobldobl_solutions()
    elif(precision == 'qd'):
        store_quaddobl_system(system)
        store_quaddobl_solutions(len(system), solutions)
        from phcpy2c import py2c_quaddobl_Newton_step
        py2c_quaddobl_Newton_step()
        result = load_quaddobl_solutions()
    elif(precision == 'mp'):
        store_multprec_system(system, decimals)
        store_multprec_solutions(len(system), solutions)
        from phcpy2c import py2c_multprec_Newton_step
        py2c_multprec_Newton_step(decimals)
        result = load_multprec_solutions()
    else:
        print 'wrong argument for precision'
        return None
    for sol in result:
        strsol = sol.split('\n')
        print strsol[-1]
    return result

def deflate(system, solutions):
    """
    The deflation method augments the given system with
    derivatives to restore the quadratic convergence of
    Newton's method at isolated singular solutions.
    After application of deflation with default settings,
    the new approximate solutions are returned.
    """
    from phcpy2c import py2c_deflate
    store_system(system)
    store_solutions(len(system), solutions)
    py2c_deflate()
    result = load_solutions()
    return result

def total_degree(pols):
    """
    Given in pols a list of string representations of polynomials,
    returns the product of the degrees of the polynomials,
    the so-called total degree which bounds the number of
    isolated solutions of the polynomial system.
    """
    from phcpy2c import py2c_syscon_total_degree
    store_system(pols)
    return py2c_syscon_total_degree()

def total_degree_start_system(pols):
    """
    Returns the system and solutions of the total degree start system
    for the polynomials represented by the strings in the list pols.
    """
    from phcpy2c import py2c_syscon_number_of_polynomials
    from phcpy2c import py2c_syscon_string_of_symbols
    from phcpy2c import py2c_syscon_degree_of_polynomial
    store_system(pols)
    dim = py2c_syscon_number_of_polynomials()
    svars = py2c_syscon_string_of_symbols()
    nvars = svars.split(' ')
    degrees = [py2c_syscon_degree_of_polynomial(k+1) for k in range(dim)]
    result = []
    for ind in range(dim):
        result.append(nvars[ind]+'^'+str(degrees[ind])+' - 1;')
    return (result, solve(result))

def linear_product_root_count(pols, silent=False):
    """
    Given in pols a list of string representations of polynomials,
    return a linear-product root count based on a supporting
    set structure of the polynomials in pols.  This root count is
    an upper bound for the number of isolated solutions.
    """
    from phcpy2c import py2c_product_supporting_set_structure
    from phcpy2c import py2c_product_write_set_structure
    from phcpy2c import py2c_product_linear_product_root_count
    store_system(pols)
    py2c_product_supporting_set_structure()
    if not silent:
        print 'a supporting set structure :'
        py2c_product_write_set_structure()
    root_count = py2c_product_linear_product_root_count()
    if not silent:
        print 'the root count :', root_count
    return root_count

def random_linear_product_system(pols, tosolve=True):
    """
    Given in pols a list of string representations of polynomials,
    returns a random linear-product system based on a supporting
    set structure and its solutions as well (if tosolve).
    """
    from phcpy2c import py2c_product_supporting_set_structure
    from phcpy2c import py2c_product_random_linear_product_system
    from phcpy2c import py2c_product_solve_linear_product_system
    store_system(pols)
    py2c_product_supporting_set_structure()
    py2c_product_random_linear_product_system()
    result = load_system()
    if not tosolve:
        return result
    py2c_product_solve_linear_product_system()
    sols = load_solutions()
    return (result, sols)

def mixed_volume(pols, stable=False):
    """
    Given in pols a list of string representations of polynomials,
    this function returns the mixed volume of the system.
    This is an interface to Algorithm 846: MixedVol of ACM TOMS,
    developed by Tangan Gao, T.Y. Li, Mengnien Wu, and Li Xing.
    If the option stable is set to True, then on return is a tuple
    containing the mixed volume and the stable mixed volume.
    The mixed volume counts the solutions with all their coordinates
    nonzero, the stable mixed volume counts all affine roots.
    Note that the stable mixed volume does not apply to systems
    with negative exponents.
    """
    from phcpy2c import py2c_celcon_clear_container
    from phcpy2c import py2c_syscon_clear_Laurent_system
    from phcpy2c import py2c_syscon_initialize_number_of_Laurentials
    from phcpy2c import py2c_syscon_store_Laurential
    from phcpy2c import py2c_mixed_volume
    py2c_celcon_clear_container()
    if stable:
        store_system(pols)
    else:
        py2c_syscon_clear_Laurent_system()
        dim = len(pols)
        py2c_syscon_initialize_number_of_Laurentials(dim)
        for ind in range(0, dim):
            lpol = pols[ind]
            nchar = len(lpol)
            py2c_syscon_store_Laurential(nchar, dim, ind+1, lpol)
    return py2c_mixed_volume(stable)

def random_coefficient_system(silent=False):
    """
    Runs the polyhedral homotopies and returns a random coefficient
    system based on the contents of the cell container.
    For this to work, the mixed_volume function must be called first.
    """
    from phcpy2c import py2c_celcon_create_random_coefficient_system
    from phcpy2c import py2c_celcon_copy_into_systems_container
    from phcpy2c import py2c_celcon_create_polyhedral_homotopy
    from phcpy2c import py2c_celcon_number_of_cells
    from phcpy2c import py2c_solcon_clear_solutions
    from phcpy2c import py2c_celcon_solve_start_system
    from phcpy2c import py2c_celcon_track_solution_path
    from phcpy2c import py2c_celcon_copy_target_solution_to_container
    py2c_celcon_create_random_coefficient_system()
    py2c_celcon_copy_into_systems_container()
    # py2c_syscon_write_system()
    result = load_system()
    # print result
    py2c_celcon_create_polyhedral_homotopy()
    nbcells = py2c_celcon_number_of_cells()
    py2c_solcon_clear_solutions()
    for cell in range(1, nbcells+1):
        mixvol = py2c_celcon_solve_start_system(cell)
        if not silent:
            print 'system %d has %d solutions' % (cell, mixvol)
        for j in range(1, mixvol+1):
            if not silent:
                print '-> tracking path %d out of %d' % (j, mixvol)
            py2c_celcon_track_solution_path(cell, j, 0)
            py2c_celcon_copy_target_solution_to_container(cell, j)
    sols = load_solutions()
    # print sols
    # newton_step(result, sols)
    return (result, sols)

def permute_system(pols):
    """
    Permutes the equations in the list of polynomials in pols
    along the permutation used in the mixed volume computation.
    """
    from phcpy2c import py2c_celcon_permute_system
    store_system(pols)
    py2c_celcon_permute_system()
    return load_system()

def test_polyhedral_homotopy():
    """
    Test on jumpstarting a polyhedral homotopy.
    """
    from phcpy2c import py2c_syscon_clear_system
    from phcpy2c import py2c_syscon_clear_Laurent_system
    py2c_syscon_clear_system()
    py2c_syscon_clear_Laurent_system()
    qrt = random_trinomials()
    mixvol = mixed_volume(qrt)
    print 'the mixed volume is', mixvol
    (rqs, rsols) = random_coefficient_system()
    print 'found %d solutions' % len(rsols)
    newton_step(rqs, rsols)
    print 'tracking to target...'
    pqr = permute_system(qrt)
    qsols = track(pqr, rqs, rsols)
    newton_step(qrt, qsols)

def test_track(silent=True, precision='d'):
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
    pols = random_trinomials()
    real_pols = real_random_trinomials(pols)
    from random import uniform as u
    qone = pols[0][:-1] + ('%+.17f' % u(-1, +1)) + ';'
    qtwo = pols[1][:-1] + ('%+.17f' % u(-1, +1)) + ';'
    rone = real_pols[0][:-1] + ('%+.17f' % u(-1, +1)) + ';'
    rtwo = real_pols[1][:-1] + ('%+.17f' % u(-1, +1)) + ';'
    start = [qone, qtwo]
    target = [rone, rtwo]
    start_sols = solve(start, silent)
    sols = track(target, start, start_sols, precision)
    mixvol = mixed_volume(target)
    print 'mixed volume of the target is', mixvol
    print 'number of solutions found :', len(sols)
    newton_step(target, sols)
    # for s in sols: print s

def test_next_track(precision='d'):
    """
    Tests the step-by-step tracking of a solution path.
    Three levels of precision are supported:
    d  : standard double precision (1.1e-15 or 2^(-53)),
    dd : double double precision (4.9e-32 or 2^(-104)),
    qd : quad double precision (1.2e-63 or 2^(-209)).
    """
    quadrics = ['x**2 + 4*y**2 - 4;', '2*y**2 - x;']
    (startsys, startsols) = total_degree_start_system(quadrics)
    print 'the first start solution :\n', startsols[0]
    if(precision == 'd'):
        initialize_standard_tracker(quadrics, startsys)
        initialize_standard_solution(2, startsols[0])
        while True:
            sol = next_standard_solution()
            print 'the next solution :\n', sol
            answer = raw_input('continue ? (y/n) ')
            if(answer != 'y'):
                break
    elif(precision == 'dd'):
        initialize_dobldobl_tracker(quadrics, startsys)
        initialize_dobldobl_solution(2, startsols[0])
        while True:
            sol = next_dobldobl_solution()
            print 'the next solution :\n', sol
            answer = raw_input('continue ? (y/n) ')
            if(answer != 'y'):
                break
    elif(precision == 'qd'):
        initialize_quaddobl_tracker(quadrics, startsys)
        initialize_quaddobl_solution(2, startsols[0])
        while True:
            sol = next_quaddobl_solution()
            print 'the next solution :\n', sol
            answer = raw_input('continue ? (y/n) ')
            if(answer != 'y'):
                break
    else:
        print 'wrong argument for precision'

def test_solver():
    """
    Generates a random trinomial system and solves it.
    """
    pols = random_trinomials()
    print 'two random trinomials :'
    print pols
    (mixvol, stable_mixvol) = mixed_volume(pols, stable=True)
    print 'its mixed volume :', mixvol
    print 'its stable mixed volume :', stable_mixvol
    sols = solve(pols)
    print 'number of computed solutions :', len(sols)
    # newton_step(pols, sols)
    # for sol in sols: print sol

def test_deflate():
    """
    Applies the deflation method to a system used as example in
    the paper by T. Ojika on Modified deflation algorithm for
    the solution of singular problems. I. A system of nonlinear
    algebraic equations, which appeared in
    J. Math. Anal. Appl. 123, 199-221, 1987.
    The approximate solutions were computed via homotopy continuation.
    The "solve" automatically deflates.
    """
    pols = [ 'x**2+y-3;', 'x+0.125*y**2-1.5;']
    sols = [ \
    "t :  1.00000000000000E+00   0.00000000000000E+00\n" + \
    "m : 1\n" + \
    "the solution for t :\n" + \
    " x : -3.00000000000000E+00   0.00000000000000E+00\n" + \
    " y : -6.00000000000000E+00   0.00000000000000E+00\n" + \
    "== err :  0.000E+00 = rco :  1.965E-01 = res :  0.000E+00 ==", \
    "t :  1.00000000000000E+00   0.00000000000000E+00\n" + \
    "m : 1\n" + \
    "the solution for t :\n" + \
    " x :  9.99995996645892E-01   2.90042038160562E-08\n" + \
    " y :  2.00000800669703E+00  -5.80082882217080E-08\n" + \
    "== err :  6.675E-06 = rco :  2.922E-12 = res :  7.423E-12 ==", \
    "t :  1.00000000000000E+00   0.00000000000000E+00\n" + \
    "m : 1\n" + \
    "the solution for t :\n" + \
    " x :  9.99997474629442E-01  -5.20158433080490E-06\n" + \
    " y :  2.00000505076133E+00   1.04031431625426E-05\n" + \
    "== err :  3.885E-06 = rco :  9.307E-12 = res :  1.863E-12 ==", \
    "t :  1.00000000000000E+00   0.00000000000000E+00\n" + \
    "m : 1\n" + \
    "the solution for t :\n" + \
    " x :  9.99998182236421E-01  -3.14843605387998E-06\n" + \
    " y :  2.00000363553377E+00   6.29686066191760E-06\n" + \
    "== err :  8.602E-08 = rco :  3.611E-12 = res :  7.957E-16 ==" ]
    print 'the system :'
    print pols
    print 'the solutions before deflation :'
    for sol in sols:
        print sol
    result = deflate(pols, sols)
    print 'the solutions after deflation :'
    for sol in result:
        print sol

def test():
    """
    Runs test_solver and test_track.
    """
    print '\ntesting polyhedral homotopy...\n'
    test_polyhedral_homotopy()
    print '\ntesting solver...\n'
    test_solver()
    print '\ntesting path tracker...\n'
    test_track()
    print '\ntesting deflation...\n'
    test_deflate()

if __name__ == "__main__":
    test()
