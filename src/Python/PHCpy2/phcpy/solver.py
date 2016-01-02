"""
The main functionality of PHCpack is its blackbox solver
and the wide variety of start systems and homotopies.
The blackbox solver is exported by the function solve.
An important task of the solver is the mixed-volume computation,
available in the function mixed_volume.
For start systems based on the degrees of the polynomials,
we have the plain total degree, m-homogeneous Bezout numbers,
and general linear-product start systems.
"""

def random_trinomials():
    """
    Returns a system of two trinomials equations for testing.
    A trinomial consists of three monomials in two variables.
    Exponents are uniform between 0 and 5 and coefficients are
    on the complex unit circle.
    """
    from random import randint as r
    exponents = [(r(0, 5), r(0, 5)) for _ in range(0, 6)]
    makemonf = lambda e: 'x^%d*y^%d' % e
    monomials = [makemonf(e) for e in exponents]
    from random import uniform as u
    from math import cos, sin, pi
    angles = [u(0, 2*pi) for _ in range(0, 6)]
    makecff = lambda a: '(' + str(cos(a)) + '%+.14f' % sin(a) + '*i)'
    cff = [makecff(a) for a in angles]
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

def random_system(dim, nbrmon, deg, cff):
    """
    Generates a random polynomial system based on the following:
    dim : number of equations and variables,
    nbrmon : maximum number of monomials per equation,
    deg : upper bound on the degree of the monomials,
    cff : type of coefficients, must be 0, 1, or 2,
    if 0, then random complex numbers on the unit circle,
    if 1, then coefficients are one (or integer multiples of one),
    if 2, then coefficients are floats in [-1,+1].
    """
    from phcpy2c import py2c_syscon_random_system
    from interface import load_standard_system
    py2c_syscon_random_system(dim, nbrmon, deg, cff)
    return load_standard_system()

def standard_solve(pols, silent=False, tasks=0):
    """
    Calls the blackbox solver.  On input in pols is a list of strings.
    By default, the solver will print to screen the computed root counts.
    To make the solver silent, set the flag silent to True.
    The number of tasks for multithreading is given by tasks.
    The solving happens in standard double precision arithmetic.
    """
    from phcpy2c import py2c_syscon_clear_standard_Laurent_system
    from phcpy2c import py2c_syscon_initialize_number_of_standard_Laurentials
    from phcpy2c import py2c_syscon_store_standard_Laurential
    from phcpy2c import py2c_solcon_clear_standard_solutions
    from phcpy2c import py2c_solve_Laurent_system
    from interface import load_standard_solutions
    py2c_syscon_clear_standard_Laurent_system()
    py2c_solcon_clear_standard_solutions()
    dim = len(pols)
    py2c_syscon_initialize_number_of_standard_Laurentials(dim)
    for ind in range(0, dim):
        pol = pols[ind]
        nchar = len(pol)
        py2c_syscon_store_standard_Laurential(nchar, dim, ind+1, pol)
    py2c_solve_Laurent_system(silent, tasks)
    return load_standard_solutions()

def dobldobl_solve(pols, silent=False, tasks=0):
    """
    Calls the blackbox solver.  On input in pols is a list of strings.
    By default, the solver will print to screen the computed root counts.
    To make the solver silent, set the flag silent to True.
    The number of tasks for multithreading is given by tasks.
    The solving happens in double double precision arithmetic.
    """
    from phcpy2c import py2c_syscon_clear_dobldobl_Laurent_system
    from phcpy2c import py2c_syscon_initialize_number_of_dobldobl_Laurentials
    from phcpy2c import py2c_syscon_store_dobldobl_Laurential
    from phcpy2c import py2c_solcon_clear_dobldobl_solutions
    from phcpy2c import py2c_solve_dobldobl_Laurent_system
    from interface import load_dobldobl_solutions
    py2c_syscon_clear_dobldobl_Laurent_system()
    py2c_solcon_clear_dobldobl_solutions()
    dim = len(pols)
    py2c_syscon_initialize_number_of_dobldobl_Laurentials(dim)
    for ind in range(0, dim):
        pol = pols[ind]
        nchar = len(pol)
        py2c_syscon_store_dobldobl_Laurential(nchar, dim, ind+1, pol)
    py2c_solve_dobldobl_Laurent_system(silent, tasks)
    return load_dobldobl_solutions()

def quaddobl_solve(pols, silent=False, tasks=0):
    """
    Calls the blackbox solver.  On input in pols is a list of strings.
    By default, the solver will print to screen the computed root counts.
    To make the solver silent, set the flag silent to True.
    The number of tasks for multithreading is given by tasks.
    The solving happens in quad double precision arithmetic.
    """
    from phcpy2c import py2c_syscon_clear_quaddobl_Laurent_system
    from phcpy2c import py2c_syscon_initialize_number_of_quaddobl_Laurentials
    from phcpy2c import py2c_syscon_store_quaddobl_Laurential
    from phcpy2c import py2c_solcon_clear_quaddobl_solutions
    from phcpy2c import py2c_solve_quaddobl_Laurent_system
    from interface import load_quaddobl_solutions
    py2c_syscon_clear_quaddobl_Laurent_system()
    py2c_solcon_clear_quaddobl_solutions()
    dim = len(pols)
    py2c_syscon_initialize_number_of_quaddobl_Laurentials(dim)
    for ind in range(0, dim):
        pol = pols[ind]
        nchar = len(pol)
        py2c_syscon_store_quaddobl_Laurential(nchar, dim, ind+1, pol)
    py2c_solve_quaddobl_Laurent_system(silent, tasks)
    return load_quaddobl_solutions()

def solve(pols, silent=False, tasks=0, precision='d'):
    """
    Calls the blackbox solver.  On input in pols is a list of strings.
    By default, the solver will print to screen the computed root counts.
    To make the solver silent, set the flag silent to True.
    The number of tasks for multithreading is given by tasks.
    Three levels of precision are supported:
    d  : standard double precision (1.1e-15 or 2^(-53)),
    dd : double double precision (4.9e-32 or 2^(-104)),
    qd : quad double precision (1.2e-63 or 2^(-209)).
    """
    if(precision == 'd'):
        return standard_solve(pols, silent, tasks)
    elif(precision == 'dd'):
        return dobldobl_solve(pols, silent, tasks)
    elif(precision == 'qd'):
        return quaddobl_solve(pols, silent, tasks)
    else:
        print 'wrong level of precision, use d, dd, or qd'

def newton_step(system, solutions, precision='d', decimals=100):
    """
    Applies one Newton step to the solutions of the system.
    For each solution, prints its last line of diagnostics.
    Four levels of precision are supported:
    d  : standard double precision (1.1e-15 or 2^(-53)),
    dd : double double precision (4.9e-32 or 2^(-104)),
    qd : quad double precision (1.2e-63 or 2^(-209)).
    mp : arbitrary precision, where the number of decimal places
    in the working precision is determined by decimals.
    """
    if(precision == 'd'):
        from interface import store_standard_system
        from interface import store_standard_solutions, load_standard_solutions
        store_standard_system(system)
        store_standard_solutions(len(system), solutions)
        from phcpy2c import py2c_standard_Newton_step
        py2c_standard_Newton_step()
        result = load_standard_solutions()
    elif(precision == 'dd'):
        from interface import store_dobldobl_system
        from interface import store_dobldobl_solutions, load_dobldobl_solutions
        store_dobldobl_system(system)
        store_dobldobl_solutions(len(system), solutions)
        from phcpy2c import py2c_dobldobl_Newton_step
        py2c_dobldobl_Newton_step()
        result = load_dobldobl_solutions()
    elif(precision == 'qd'):
        from interface import store_quaddobl_system
        from interface import store_quaddobl_solutions, load_quaddobl_solutions
        store_quaddobl_system(system)
        store_quaddobl_solutions(len(system), solutions)
        from phcpy2c import py2c_quaddobl_Newton_step
        py2c_quaddobl_Newton_step()
        result = load_quaddobl_solutions()
    elif(precision == 'mp'):
        from interface import store_multprec_system
        from interface import store_multprec_solutions, load_multprec_solutions
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

def newton_laurent_step(system, solutions, precision='d', decimals=100):
    """
    Applies one Newton step to the solutions of the Laurent system.
    For each solution, prints its last line of diagnostics.
    Four levels of precision are supported:
    d  : standard double precision (1.1e-15 or 2^(-53)),
    dd : double double precision (4.9e-32 or 2^(-104)),
    qd : quad double precision (1.2e-63 or 2^(-209)).
    mp : arbitrary precision, where the number of decimal places
    in the working precision is determined by decimals.
    """
    if(precision == 'd'):
        from interface import store_standard_laurent_system
        from interface import store_standard_solutions, load_standard_solutions
        store_standard_laurent_system(system)
        store_standard_solutions(len(system), solutions)
        from phcpy2c import py2c_standard_Newton_Laurent_step
        py2c_standard_Newton_Laurent_step()
        result = load_standard_solutions()
    elif(precision == 'dd'):
        from interface import store_dobldobl_laurent_system
        from interface import store_dobldobl_solutions, load_dobldobl_solutions
        store_dobldobl_laurent_system(system)
        store_dobldobl_solutions(len(system), solutions)
        from phcpy2c import py2c_dobldobl_Newton_Laurent_step
        py2c_dobldobl_Newton_Laurent_step()
        result = load_dobldobl_solutions()
    elif(precision == 'qd'):
        from interface import store_quaddobl_laurent_system
        from interface import store_quaddobl_solutions, load_quaddobl_solutions
        store_quaddobl_laurent_system(system)
        store_quaddobl_solutions(len(system), solutions)
        from phcpy2c import py2c_quaddobl_Newton_Laurent_step
        py2c_quaddobl_Newton_Laurent_step()
        result = load_quaddobl_solutions()
    elif(precision == 'mp'):
        from interface import store_multprec_laurent_system
        from interface import store_multprec_solutions, load_multprec_solutions
        store_multprec_laurent_system(system, decimals)
        store_multprec_solutions(len(system), solutions)
        from phcpy2c import py2c_multprec_Newton_Laurent_step
        py2c_multprec_Newton_Laurent_step(decimals)
        result = load_multprec_solutions()
    else:
        print 'wrong argument for precision'
        return None
    for sol in result:
        strsol = sol.split('\n')
        print strsol[-1]
    return result

def newton_steps(system, solutions, accuracy=8, maxsteps=4, maxprec=256):
    """
    Runs a sequence of variable precision Newton steps to approximate
    solutions accurate up to a specified number of decimal places.
    In addition to the system and solutions, there are three parameters:
    accuracy : number of decimal places wanted to be accurate,
    maxsteps : maximum number of Newton steps,
    maxprec : maximum number of decimal places in the precision used
    to estimate the condition numbers.
    """
    from phcpy2c import py2c_varbprec_Newton_Laurent_steps as vmpnewt
    from interface import store_multprec_solutions, load_multprec_solutions
    dim = len(system)
    store_multprec_solutions(dim, solutions)
    pols = ""
    for polynomial in system:
        pols = pols + polynomial
    vmpnewt(dim, accuracy, maxsteps, maxprec, len(pols), pols)
    result = load_multprec_solutions()
    return result

def standard_deflate(system, solutions):
    """
    The deflation method augments the given system with
    derivatives to restore the quadratic convergence of
    Newton's method at isolated singular solutions,
    in standard double precision.
    After application of deflation with default settings,
    the new approximate solutions are returned.
    """
    from phcpy2c import py2c_standard_deflate
    from interface import store_standard_system
    from interface import store_standard_solutions, load_standard_solutions
    store_standard_system(system)
    store_standard_solutions(len(system), solutions)
    py2c_standard_deflate()
    result = load_standard_solutions()
    return result

def dobldobl_deflate(system, solutions):
    """
    The deflation method augments the given system with
    derivatives to restore the quadratic convergence of
    Newton's method at isolated singular solutions,
    in double double precision.
    After application of deflation with default settings,
    the new approximate solutions are returned.
    """
    from phcpy2c import py2c_dobldobl_deflate
    from interface import store_dobldobl_system
    from interface import store_dobldobl_solutions, load_dobldobl_solutions
    store_dobldobl_system(system)
    store_dobldobl_solutions(len(system), solutions)
    py2c_dobldobl_deflate()
    result = load_dobldobl_solutions()
    return result

def quaddobl_deflate(system, solutions):
    """
    The deflation method augments the given system with
    derivatives to restore the quadratic convergence of
    Newton's method at isolated singular solutions,
    in quad double precision.
    After application of deflation with default settings,
    the new approximate solutions are returned.
    """
    from phcpy2c import py2c_quaddobl_deflate
    from interface import store_quaddobl_system
    from interface import store_quaddobl_solutions, load_quaddobl_solutions
    store_quaddobl_system(system)
    store_quaddobl_solutions(len(system), solutions)
    py2c_quaddobl_deflate()
    result = load_quaddobl_solutions()
    return result

def total_degree(pols):
    """
    Given in pols a list of string representations of polynomials,
    returns the product of the degrees of the polynomials,
    the so-called total degree which bounds the number of
    isolated solutions of the polynomial system.
    """
    from phcpy2c import py2c_syscon_total_degree
    from interface import store_standard_system
    store_standard_system(pols)
    return py2c_syscon_total_degree()

def total_degree_start_system(pols):
    """
    Returns the system and solutions of the total degree start system
    for the polynomials represented by the strings in the list pols.
    """
    from phcpy2c import py2c_syscon_number_of_standard_polynomials
    from phcpy2c import py2c_syscon_string_of_symbols
    from phcpy2c import py2c_syscon_degree_of_standard_polynomial
    from interface import store_standard_system
    store_standard_system(pols)
    dim = py2c_syscon_number_of_standard_polynomials()
    svars = py2c_syscon_string_of_symbols()
    nvars = svars.split(' ')
    degrees = [py2c_syscon_degree_of_standard_polynomial(k+1) \
               for k in range(dim)]
    result = []
    for ind in range(dim):
        result.append(nvars[ind]+'^'+str(degrees[ind])+' - 1;')
    return (result, solve(result))

def m_homogeneous_bezout_number(pols):
    """
    Given in pols a list of string representations of polynomials,
    in as many variables as the elements in the list,
    this function applies a heuristic to generate a partition of the
    set of unknowns to exploit the product structure of the system.
    On return are the m-homogeneous Bezout number and the partition
    of the set of unknowns.  If the partition equals the entire
    set of unknowns, then the 1-homogeneous Bezout number equals
    the total degree of the system.
    """
    from phcpy2c import py2c_product_m_homogeneous_Bezout_number
    from interface import store_standard_system
    store_standard_system(pols)
    result = py2c_product_m_homogeneous_Bezout_number()
    return result

def m_partition_bezout_number(pols, partition):
    """
    There are as many m-homogeneous Bezout numbers as there are
    partitions of the set of unknowns of a polynomial system.
    Given in pols the string representations of a polynomial system
    in as many variables as equations, and a string representation of
    a partition of the set of unknowns, this function returns the
    m-homogeneous Bezout number corresponding to the given partition.
    """
    from phcpy2c import py2c_product_m_partition_Bezout_number
    from interface import store_standard_system
    store_standard_system(pols)
    return py2c_product_m_partition_Bezout_number(len(partition), partition)

def m_homogeneous_start_system(pols, partition):
    """
    For an m-homogeneous Bezout number of a polynomial system defined by
    a partition of the set of unknowns, one can define a linear-product
    system that has exactly as many regular solutions as the Bezount number.
    This linear-product system can then be used as start system in a
    homotopy to compute all isolated solutions of any polynomial system
    with the same m-homogeneous structure.
    This function returns a linear-product start system with random
    coefficients and its solutions for the given polynomials in pols
    and the partition.
    """
    from phcpy2c import py2c_product_m_homogeneous_start_system
    from phcpy2c import py2c_product_solve_linear_product_system
    from interface import store_standard_system, load_standard_system
    from interface import load_standard_solutions
    store_standard_system(pols)
    py2c_product_m_homogeneous_start_system(len(partition), partition)
    result = load_standard_system()
    py2c_product_solve_linear_product_system()
    sols = load_standard_solutions()
    return (result, sols)

def linear_product_root_count(pols, silent=False):
    """
    Given in pols a list of string representations of polynomials,
    returns a linear-product root count based on a supporting
    set structure of the polynomials in pols.  This root count is
    an upper bound for the number of isolated solutions.
    """
    from phcpy2c import py2c_product_supporting_set_structure
    from phcpy2c import py2c_product_write_set_structure
    from phcpy2c import py2c_product_linear_product_root_count
    from interface import store_standard_system
    store_standard_system(pols)
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
    from interface import store_standard_system, load_standard_system
    from interface import load_standard_solutions
    store_standard_system(pols)
    py2c_product_supporting_set_structure()
    py2c_product_random_linear_product_system()
    result = load_standard_system()
    if not tosolve:
        return result
    py2c_product_solve_linear_product_system()
    sols = load_standard_solutions()
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
    Incorrectly parsed strings will result in a negative value on return.
    """
    from phcpy2c import py2c_celcon_clear_container
    from phcpy2c import py2c_mixed_volume
    from interface import store_standard_system, store_standard_laurent_system
    py2c_celcon_clear_container()
    if stable:
        fail = store_standard_system(pols)
    else:
        fail = store_standard_laurent_system(pols)
    if(fail != 0):
        return -fail  # a negative number is clearly wrong
    else:
        return py2c_mixed_volume(stable)

def standard_random_coefficient_system(silent=False):
    """
    Runs the polyhedral homotopies and returns a random coefficient
    system based on the contents of the cell container,
    in standard double precision arithmetic.
    For this to work, the mixed_volume function must be called first.
    """
    from phcpy2c import py2c_celcon_standard_random_coefficient_system
    from phcpy2c import py2c_celcon_copy_into_standard_systems_container
    from phcpy2c import py2c_celcon_standard_polyhedral_homotopy
    from phcpy2c import py2c_celcon_number_of_cells
    from phcpy2c import py2c_solcon_clear_standard_solutions
    from phcpy2c import py2c_celcon_solve_standard_start_system
    from phcpy2c import py2c_celcon_track_standard_solution_path
    from phcpy2c import py2c_celcon_copy_target_standard_solution_to_container
    py2c_celcon_standard_random_coefficient_system()
    py2c_celcon_copy_into_standard_systems_container()
    from interface import load_standard_system, load_standard_solutions
    # py2c_syscon_write_system()
    result = load_standard_system()
    # print result
    py2c_celcon_standard_polyhedral_homotopy()
    nbcells = py2c_celcon_number_of_cells()
    py2c_solcon_clear_standard_solutions()
    for cell in range(1, nbcells+1):
        mixvol = py2c_celcon_solve_standard_start_system(cell)
        if not silent:
            print 'system %d has %d solutions' % (cell, mixvol)
        for j in range(1, mixvol+1):
            if not silent:
                print '-> tracking path %d out of %d' % (j, mixvol)
            py2c_celcon_track_standard_solution_path(cell, j, 0)
            py2c_celcon_copy_target_standard_solution_to_container(cell, j)
    sols = load_standard_solutions()
    # print sols
    # newton_step(result, sols)
    return (result, sols)

def dobldobl_random_coefficient_system(silent=False):
    """
    Runs the polyhedral homotopies and returns a random coefficient
    system based on the contents of the cell container,
    in double double precision arithmetic.
    For this to work, the mixed_volume function must be called first.
    """
    from phcpy2c import py2c_celcon_dobldobl_random_coefficient_system
    from phcpy2c import py2c_celcon_copy_into_dobldobl_systems_container
    from phcpy2c import py2c_celcon_dobldobl_polyhedral_homotopy
    from phcpy2c import py2c_celcon_number_of_cells
    from phcpy2c import py2c_solcon_clear_dobldobl_solutions
    from phcpy2c import py2c_celcon_solve_dobldobl_start_system
    from phcpy2c import py2c_celcon_track_dobldobl_solution_path
    from phcpy2c \
        import py2c_celcon_copy_target_dobldobl_solution_to_container
    from interface import load_dobldobl_system, load_dobldobl_solutions
    py2c_celcon_dobldobl_random_coefficient_system()
    py2c_celcon_copy_into_dobldobl_systems_container()
    # py2c_syscon_write_dobldobl_system()
    result = load_dobldobl_system()
    # print result
    py2c_celcon_dobldobl_polyhedral_homotopy()
    nbcells = py2c_celcon_number_of_cells()
    py2c_solcon_clear_dobldobl_solutions()
    for cell in range(1, nbcells+1):
        mixvol = py2c_celcon_solve_dobldobl_start_system(cell)
        if not silent:
            print 'system %d has %d solutions' % (cell, mixvol)
        for j in range(1, mixvol+1):
            if not silent:
                print '-> tracking path %d out of %d' % (j, mixvol)
            py2c_celcon_track_dobldobl_solution_path(cell, j, 0)
            py2c_celcon_copy_target_dobldobl_solution_to_container(cell, j)
    sols = load_dobldobl_solutions()
    # print sols
    # newton_step(result, sols)
    return (result, sols)

def quaddobl_random_coefficient_system(silent=False):
    """
    Runs the polyhedral homotopies and returns a random coefficient
    system based on the contents of the cell container,
    in quad double precision arithmetic.
    For this to work, the mixed_volume function must be called first.
    """
    from phcpy2c import py2c_celcon_quaddobl_random_coefficient_system
    from phcpy2c import py2c_celcon_copy_into_quaddobl_systems_container
    from phcpy2c import py2c_celcon_quaddobl_polyhedral_homotopy
    from phcpy2c import py2c_celcon_number_of_cells
    from phcpy2c import py2c_solcon_clear_quaddobl_solutions
    from phcpy2c import py2c_celcon_solve_quaddobl_start_system
    from phcpy2c import py2c_celcon_track_quaddobl_solution_path
    from phcpy2c \
        import py2c_celcon_copy_target_quaddobl_solution_to_container
    from interface import load_quaddobl_system, load_quaddobl_solutions
    py2c_celcon_quaddobl_random_coefficient_system()
    py2c_celcon_copy_into_quaddobl_systems_container()
    # py2c_syscon_write_dobldobl_system()
    result = load_quaddobl_system()
    # print result
    py2c_celcon_quaddobl_polyhedral_homotopy()
    nbcells = py2c_celcon_number_of_cells()
    py2c_solcon_clear_quaddobl_solutions()
    for cell in range(1, nbcells+1):
        mixvol = py2c_celcon_solve_quaddobl_start_system(cell)
        if not silent:
            print 'system %d has %d solutions' % (cell, mixvol)
        for j in range(1, mixvol+1):
            if not silent:
                print '-> tracking path %d out of %d' % (j, mixvol)
            py2c_celcon_track_quaddobl_solution_path(cell, j, 0)
            py2c_celcon_copy_target_quaddobl_solution_to_container(cell, j)
    sols = load_quaddobl_solutions()
    # print sols
    # newton_step(result, sols)
    return (result, sols)

def random_coefficient_system(silent=False, precision='d'):
    """
    Runs the polyhedral homotopies and returns a random coefficient
    system based on the contents of the cell container.
    For this to work, the mixed_volume function must be called first.
    Three levels of precision are supported:
    d  : standard double precision (1.1e-15 or 2^(-53)),
    dd : double double precision (4.9e-32 or 2^(-104)),
    qd : quad double precision (1.2e-63 or 2^(-209)).
    """
    if(precision == 'd'):
        return standard_random_coefficient_system(silent)
    elif(precision == 'dd'):
        return dobldobl_random_coefficient_system(silent)

def permute_standard_system(pols):
    """
    Permutes the equations in the list of polynomials in pols
    with coefficients in standard double precision,
    along the permutation used in the mixed volume computation.
    """
    from phcpy2c import py2c_celcon_permute_standard_system
    from interface import store_standard_system, load_standard_system
    store_standard_system(pols)
    py2c_celcon_permute_standard_system()
    return load_standard_system()

def permute_dobldobl_system(pols):
    """
    Permutes the equations in the list of polynomials in pols
    with coefficients in double double precision,
    along the permutation used in the mixed volume computation.
    """
    from phcpy2c import py2c_celcon_permute_dobldobl_system
    from interface import store_dobldobl_system, load_dobldobl_system
    store_dobldobl_system(pols)
    py2c_celcon_permute_dobldobl_system()
    return load_dobldobl_system()

def permute_quaddobl_system(pols):
    """
    Permutes the equations in the list of polynomials in pols
    with coefficients in quad double precision,
    along the permutation used in the mixed volume computation.
    """
    from phcpy2c import py2c_celcon_permute_quaddobl_system
    from interface import store_quaddobl_system, load_quaddobl_system
    store_quaddobl_system(pols)
    py2c_celcon_permute_quaddobl_system()
    return load_quaddobl_system()

def standard_usolve(pol, mxi, eps):
    """
    Applies the method of Durand-Kerner (aka Weierstrass)
    to the polynomial in the string pol, in standard double precision
    The maximum number of iterations is in mxi,
    the requirement on the accuracy in eps.
    """
    from phcpy2c import py2c_usolve_standard
    from interface import store_standard_system, load_standard_solutions
    store_standard_system([pol])
    nit = py2c_usolve_standard(mxi, eps)
    rts = load_standard_solutions()
    return (nit, rts)

def dobldobl_usolve(pol, mxi, eps):
    """
    Applies the method of Durand-Kerner (aka Weierstrass)
    to the polynomial in the string pol, in double double precision
    The maximum number of iterations is in mxi,
    the requirement on the accuracy in eps.
    """
    from phcpy2c import py2c_usolve_dobldobl
    from interface import store_dobldobl_system, load_dobldobl_solutions
    store_dobldobl_system([pol])
    nit = py2c_usolve_dobldobl(mxi, eps)
    rts = load_dobldobl_solutions()
    return (nit, rts)

def quaddobl_usolve(pol, mxi, eps):
    """
    Applies the method of Durand-Kerner (aka Weierstrass)
    to the polynomial in the string pol, in quad double precision
    The maximum number of iterations is in mxi,
    the requirement on the accuracy in eps.
    """
    from phcpy2c import py2c_usolve_quaddobl
    from interface import store_quaddobl_system, load_quaddobl_solutions
    store_quaddobl_system([pol])
    nit = py2c_usolve_quaddobl(mxi, eps)
    rts = load_quaddobl_solutions()
    return (nit, rts)

def multprec_usolve(pol, mxi, eps, decimals):
    """
    Applies the method of Durand-Kerner (aka Weierstrass)
    to the polynomial in the string pol, in arbitrary multiprecision,
    the number of decimal places in the precision is in decimals.
    The maximum number of iterations is in mxi,
    the requirement on the accuracy in eps.
    """
    from phcpy2c import py2c_usolve_multprec
    from interface import store_multprec_system, load_multprec_solutions
    store_multprec_system([pol], decimals)
    nit = py2c_usolve_multprec(decimals, mxi, eps)
    rts = load_multprec_solutions()
    return (nit, rts)

def usolve(pol, mxi, eps, precision='d', decimals=100):
    """
    Applies the method of Durand-Kerner (aka Weierstrass)
    to the polynomial in the string pol.
    The maximum number of iterations is in mxi,
    the requirement on the accuracy in eps.
    Four levels of precision are supported:
    d  : standard double precision (1.1e-15 or 2^(-53)),
    dd : double double precision (4.9e-32 or 2^(-104)),
    qd : quad double precision (1.2e-63 or 2^(-209)).
    mp : arbitrary precision, where the number of decimal places
    in the working precision is determined by decimals.
    """
    if(precision == 'd'):
        return standard_usolve(pol, mxi, eps)
    elif(precision == 'dd'):
        return dobldobl_usolve(pol, mxi, eps)
    elif(precision == 'qd'):
        return quaddobl_usolve(pol, mxi, eps)
    else:
        return multprec_usolve(pol, mxi, eps, decimals)

def test_usolve():
    """
    Does a simple sanity check on solving a univariate polynomial
    at various levels of precision.
    """
    sqrt2pol = 'x^2 - 2;'
    (nit, roots) = usolve(sqrt2pol, 10, 1.0e-12)
    print 'number of iterations :', nit
    for root in roots:
        print root
    (nit, roots) = usolve(sqrt2pol, 20, 1.0e-28, 'dd')
    print 'number of iterations :', nit
    for root in roots:
        print root
    (nit, roots) = usolve(sqrt2pol, 30, 1.0e-48, 'qd')
    print 'number of iterations :', nit
    for root in roots:
        print root
    (nit, roots) = usolve(sqrt2pol, 40, 1.0e-88, 'md', 100)
    print 'number of iterations :', nit
    for root in roots:
        print root

def standard_scale_system(pols):
    """
    Applies equation and variable scaling in standard double precision
    to the polynomials in the list pols.
    On return is the list of scaled polynomials and the scaling coefficients.
    """
    from interface import store_standard_system, load_standard_system
    from phcpy2c import py2c_scale_standard_system
    store_standard_system(pols)
    cffs = py2c_scale_standard_system(2)
    spol = load_standard_system()
    return (spol, cffs)

def dobldobl_scale_system(pols):
    """
    Applies equation and variable scaling in double double precision
    to the polynomials in the list pols.
    On return is the list of scaled polynomials and the scaling coefficients.
    """
    from interface import store_dobldobl_system, load_dobldobl_system
    from phcpy2c import py2c_scale_dobldobl_system
    store_dobldobl_system(pols)
    cffs = py2c_scale_dobldobl_system(2)
    spol = load_dobldobl_system()
    return (spol, cffs)

def quaddobl_scale_system(pols):
    """
    Applies equation and variable scaling in quad double precision
    to the polynomials in the list pols.
    On return is the list of scaled polynomials and the scaling coefficients.
    """
    from interface import store_quaddobl_system, load_quaddobl_system
    from phcpy2c import py2c_scale_quaddobl_system
    store_quaddobl_system(pols)
    cffs = py2c_scale_quaddobl_system(2)
    spol = load_quaddobl_system()
    return (spol, cffs)

def standard_scale_solutions(nvar, sols, cffs):
    """
    Scales the solutions in the list sols using the coefficients in cffs,
    using standard double precision arithmetic.
    The number of variables is given in the parameter nvar.
    If the sols are the solution of the polynomials in the output of
    standard_scale_system(pols), then the solutions on return will be
    solutions of the original polynomials in the list pols.
    """
    from interface import store_standard_solutions, load_standard_solutions
    from phcpy2c import py2c_scale_standard_solutions
    store_standard_solutions(nvar, sols)
    py2c_scale_standard_solutions(len(cffs), str(cffs))
    return load_standard_solutions()

def dobldobl_scale_solutions(nvar, sols, cffs):
    """
    Scales the solutions in the list sols using the coefficients in cffs,
    using double double precision arithmetic.
    The number of variables is given in the parameter nvar.
    If the sols are the solution of the polynomials in the output of
    dobldobl_scale_system(pols), then the solutions on return will be
    solutions of the original polynomials in the list pols.
    """
    from interface import store_dobldobl_solutions, load_dobldobl_solutions
    from phcpy2c import py2c_scale_dobldobl_solutions
    store_dobldobl_solutions(nvar, sols)
    py2c_scale_dobldobl_solutions(len(cffs), str(cffs))
    return load_dobldobl_solutions()

def quaddobl_scale_solutions(nvar, sols, cffs):
    """
    Scales the solutions in the list sols using the coefficients in cffs,
    using quad double precision arithmetic.
    The number of variables is given in the parameter nvar.
    If the sols are the solution of the polynomials in the output of
    quaddobl_scale_system(pols), then the solutions on return will be
    solutions of the original polynomials in the list pols.
    """
    from interface import store_quaddobl_solutions, load_quaddobl_solutions
    from phcpy2c import py2c_scale_quaddobl_solutions
    store_quaddobl_solutions(nvar, sols)
    py2c_scale_quaddobl_solutions(len(cffs), str(cffs))
    return load_quaddobl_solutions()

def test_scale():
    """
    Performs a basic test on variable scaling.
    """
    orgp = ['100*x^2 + 10*x + 1;']
    print 'a badly scaled problem :', orgp
    print 'scaling a polynomial in standard double precision'
    (scp, cff) = standard_scale_system(orgp)
    print 'the scaled polynomial', scp
    print 'the scaling coefficients :', cff[0:-1]
    print 'estimated inverse condition number :', cff[-1]
    print 'solving the scaled problem ...'
    scpsols = solve(scp)
    print 'the solutions of the scaled problem :'
    for sol in scpsols:
        print sol
    orgpsols = standard_scale_solutions(len(scp), scpsols, cff[0:-1])
    print 'the solutions of the original problem :'
    for sol in orgpsols:
        print sol
    print 'scaling a polynomial in double double precision'
    (scp, cff) = dobldobl_scale_system(orgp)
    print 'the scaled polynomial', scp
    print 'the scaling coefficients :', cff[0:-1]
    print 'estimated inverse condition number :', cff[-1]
    print 'solving the scaled problem ...'
    scpsols = solve(scp)
    dd_scpsols = newton_step(scp, scpsols, precision='dd')
    print 'the solutions of the scaled problem :'
    for sol in dd_scpsols:
        print sol
    scpsols = dobldobl_scale_solutions(len(scp), dd_scpsols, cff[0:-1])
    print 'the solutions of the original problem :'
    for sol in scpsols:
        print sol
    print 'scaling a polynomial in quad double precision'
    (scp, cff) = quaddobl_scale_system(orgp)
    print 'the scaled polynomial', scp
    print 'the scaling coefficients :', cff[0:-1]
    print 'estimated inverse condition number :', cff[-1]
    print 'solving the scaled problem ...'
    scpsols = solve(scp)
    dd_scpsols = newton_step(scp, scpsols, precision='dd')
    qd_scpsols = newton_step(scp, dd_scpsols, precision='qd')
    print 'the solutions of the scaled problem :'
    for sol in qd_scpsols:
        print sol
    orgpsols = quaddobl_scale_solutions(len(scp), qd_scpsols, cff[0:-1])
    print 'the solutions of the original problem :'
    for sol in orgpsols:
        print sol

def test_dobldobl_polyhedral_homotopy():
    """
    Test polyhedral homotopy in double double precision
    on a small random polynomial system.
    """
    qrt = random_system(2, 3, 3, 1)
    print 'a random system :', qrt
    mixvol = mixed_volume(qrt)
    print 'the mixed volume is', mixvol
    (rqs, rsols) = dobldobl_random_coefficient_system()
    print 'the random coefficient system:'
    for pol in rqs:
        print pol
    print 'found %d solutions' % len(rsols)
    newton_step(rqs, rsols, precision='dd')

def test_quaddobl_polyhedral_homotopy():
    """
    Test polyhedral homotopy in quad double precision
    on a small random polynomial system.
    """
    qrt = random_system(2, 3, 3, 1)
    print 'a random system :', qrt
    mixvol = mixed_volume(qrt)
    print 'the mixed volume is', mixvol
    (rqs, rsols) = quaddobl_random_coefficient_system()
    print 'the random coefficient system:'
    for pol in rqs:
        print pol
    print 'found %d solutions' % len(rsols)
    newton_step(rqs, rsols, precision='qd')

def test_standard_polyhedral_homotopy():
    """
    Test on jumpstarting a polyhedral homotopy
    in standard precision.
    """
    from phcpy2c import py2c_syscon_clear_standard_system
    from phcpy2c import py2c_syscon_clear_standard_Laurent_system
    from trackers import track
    py2c_syscon_clear_standard_system()
    py2c_syscon_clear_standard_Laurent_system()
    qrt = random_trinomials()
    mixvol = mixed_volume(qrt)
    print 'the mixed volume is', mixvol
    (rqs, rsols) = random_coefficient_system()
    print 'found %d solutions' % len(rsols)
    newton_step(rqs, rsols)
    print 'tracking to target...'
    pqr = permute_standard_system(qrt)
    qsols = track(pqr, rqs, rsols)
    newton_step(qrt, qsols)

def test_polyhedral_homotopy(precision='d'):
    """
    Test polyhedral homotopy on small random systems
    for standard double precision (d), double double precision (dd),
    or quad double precision (qd).
    """
    if(precision == 'd'):
        test_standard_polyhedral_homotopy()
    elif(precision == 'dd'):
        test_dobldobl_polyhedral_homotopy()
    elif(precision == 'qd'):
        test_quaddobl_polyhedral_homotopy()
    else:
        print 'precision not recognized'

def test_newton():
    """
    Tests Newton's method on simple polynomial system,
    refining the square root of 2 with increasing precision.
    """
    pols = ['x*y - 1;', 'x^2 - 2;']
    from phcpy.solutions import make_solution
    sol = make_solution(['x', 'y'], [1.414, 0.707])
    sols = [sol]
    print 'start solution :\n', sols[0]
    for k in range(1, 4):
        sols = newton_step(pols, sols)
        print 'at step', k, ':\n', sols[0]
    for k in range(4, 6):
        sols = newton_step(pols, sols, precision='dd')
        print 'at step', k, ':\n', sols[0]
    for k in range(6, 8):
        sols = newton_step(pols, sols, precision='qd')
        print 'at step', k, ':\n', sols[0]
    for k in range(8, 10):
        sols = newton_step(pols, sols, precision='mp')
        print 'at step', k, ':\n', sols[0]

def test_newton_laurent():
    """
    Tests Newton's method on simple Laurent system,
    refining the square root of 2 with increasing precision.
    """
    laurpols = ['x - y^-1;', 'x^2 - 2;']
    from phcpy.solutions import make_solution
    sol = make_solution(['x', 'y'], [1.414, 0.707])
    sols = [sol]
    print 'start solution :\n', sols[0]
    for k in range(1, 4):
        sols = newton_laurent_step(laurpols, sols)
        print 'at step', k, ':\n', sols[0]
    for k in range(4, 6):
        sols = newton_laurent_step(laurpols, sols, precision='dd')
        print 'at step', k, ':\n', sols[0]
    for k in range(6, 8):
        sols = newton_laurent_step(laurpols, sols, precision='qd')
        print 'at step', k, ':\n', sols[0]
    for k in range(8, 10):
        sols = newton_laurent_step(laurpols, sols, precision='mp')
        print 'at step', k, ':\n', sols[0]

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
    result = standard_deflate(pols, sols)
    print 'the solutions after deflation in standard double precision:'
    for sol in result:
        print sol
    result = dobldobl_deflate(pols, sols)
    print 'the solutions after deflation in double double precision:'
    for sol in result:
        print sol
    result = quaddobl_deflate(pols, sols)
    print 'the solutions after deflation in quad double precision:'
    for sol in result:
        print sol

def test():
    """
    Runs test_polyhedral_homotopy, test_solver and test_deflate.
    """
    print '\ntesting polyhedral homotopy...\n'
    test_polyhedral_homotopy()
    print '\ntesting solver...\n'
    test_solver()
    print '\ntesting deflation...\n'
    test_deflate()

if __name__ == "__main__":
    test()
