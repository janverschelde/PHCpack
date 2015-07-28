"""
This module exports routines of PHCpack to manipulate
positive dimensional solution sets of polynomial systems.
"""

def embed(nvar, topdim, pols):
    """
    Given in pols a list of strings that represent
    polynomials in nvar variables, this function
    returns an embedding of pols of dimension topdim.
    The topdim abbreviates the top dimension which
    equals the expected highest dimension of a component
    of the solution set of the system of polynomials.
    """
    from phcpy2c import py2c_syscon_clear_standard_system
    from phcpy2c import py2c_syscon_initialize_number_of_standard_polynomials
    from phcpy2c import py2c_syscon_store_standard_polynomial
    from phcpy2c import py2c_syscon_load_standard_polynomial
    from phcpy2c import py2c_embed_system
    py2c_syscon_clear_standard_system()
    nequ = len(pols)
    if nequ > nvar:
        py2c_syscon_initialize_number_of_standard_polynomials(nequ)
        nbres = nequ
    else:
        py2c_syscon_initialize_number_of_standard_polynomials(nvar)
        nbres = nvar
    for i in range(0, nequ):
        nchar = len(pols[i])
        py2c_syscon_store_standard_polynomial(nchar, nvar, i+1, pols[i])
    py2c_embed_system(topdim)
    # py2c_syscon_write_system()
    result = []
    for i in range(1, nbres+topdim+1):
        result.append(py2c_syscon_load_standard_polynomial(i))
    return result

def witness_set_of_hypersurface(nvar, hpol):
    """
    Given in hpol the string representation of a polynomial
    in nvar variables (ending with a semicolon),
    on return is an embedded system and its solutions
    which represents a witness set for hpol.
    The number of solutions on return should equal
    the degree of the polynomial in hpol.
    """
    from phcpy2c import py2c_witness_set_of_hypersurface
    py2c_witness_set_of_hypersurface(nvar, len(hpol), hpol)
    from interface import load_standard_system, load_standard_solutions
    pols = load_standard_system()
    sols = load_standard_solutions()
    return (pols, sols)

def drop_variable_from_polynomials(pols, svar):
    """
    Removes the variable with symbol in the string svar
    from the list pols of strings that represented
    polynomials in several variables.
    """
    from phcpy2c import py2c_syscon_standard_drop_variable_by_name
    from phcpy2c import py2c_syscon_remove_symbol_name
    from interface import store_standard_system, load_standard_system
    store_standard_system(pols)
    py2c_syscon_standard_drop_variable_by_name(len(svar), svar)
    py2c_syscon_remove_symbol_name(len(svar), svar)
    return load_standard_system()

def drop_coordinate_from_solutions(sols, nbvar, svar):
    """
    Removes the variable with symbol in the string svar
    from the list sols of strings that represent solutions
    in nbvar variables.
    """
    from phcpy2c import py2c_syscon_clear_symbol_table
    from phcpy2c import py2c_solcon_standard_drop_coordinate_by_name
    from phcpy2c import py2c_syscon_remove_symbol_name
    from interface import store_standard_solutions, load_standard_solutions
    py2c_syscon_clear_symbol_table()
    store_standard_solutions(nbvar, sols)
    py2c_solcon_standard_drop_coordinate_by_name(len(svar), svar)
    py2c_syscon_remove_symbol_name(len(svar), svar)
    return load_standard_solutions()

def standard_double_cascade_step(embsys, esols):
    """
    Given in embsys an embedded polynomial system and
    solutions with nonzero slace variables in esols,
    does one step in the homotopy cascade,
    with standard double precision arithmetic.
    The list on return contains witness points on
    lower dimensional solution components.
    """
    from phcpy2c import py2c_copy_container_to_start_system
    from phcpy2c import py2c_copy_container_to_start_solutions
    from phcpy2c import py2c_standard_cascade_homotopy
    from phcpy2c import py2c_solve_by_standard_homotopy_continuation
    from phcpy2c import py2c_solcon_clear_standard_solutions
    from phcpy2c import py2c_copy_target_solutions_to_container
    from interface import store_standard_system
    from interface import store_standard_solutions, load_standard_solutions
    store_standard_system(embsys)
    py2c_copy_container_to_start_system()
    store_standard_solutions(len(embsys), esols)
    py2c_copy_container_to_start_solutions()
    py2c_standard_cascade_homotopy()
    py2c_solve_by_standard_homotopy_continuation(0)
    py2c_solcon_clear_standard_solutions()
    py2c_copy_target_solutions_to_container()
    return load_standard_solutions()

def double_double_cascade_step(embsys, esols):
    """
    Given in embsys an embedded polynomial system and
    solutions with nonzero slace variables in esols,
    does one step in the homotopy cascade,
    with double double precision arithmetic.
    The list on return contains witness points on
    lower dimensional solution components.
    """
    from phcpy2c import py2c_copy_dobldobl_container_to_start_system
    from phcpy2c import py2c_copy_dobldobl_container_to_start_solutions
    from phcpy2c import py2c_dobldobl_cascade_homotopy
    from phcpy2c import py2c_solve_by_dobldobl_homotopy_continuation
    from phcpy2c import py2c_solcon_clear_dobldobl_solutions
    from phcpy2c import py2c_copy_dobldobl_target_solutions_to_container
    from interface import store_dobldobl_system
    from interface import store_dobldobl_solutions, load_dobldobl_solutions
    store_dobldobl_system(embsys)
    py2c_copy_dobldobl_container_to_start_system()
    store_dobldobl_solutions(len(embsys), esols)
    py2c_copy_dobldobl_container_to_start_solutions()
    py2c_dobldobl_cascade_homotopy()
    py2c_solve_by_dobldobl_homotopy_continuation()
    py2c_solcon_clear_dobldobl_solutions()
    py2c_copy_dobldobl_target_solutions_to_container()
    return load_dobldobl_solutions()

def quad_double_cascade_step(embsys, esols):
    """
    Given in embsys an embedded polynomial system and
    solutions with nonzero slace variables in esols,
    does one step in the homotopy cascade,
    with quad double precision arithmetic.
    The list on return contains witness points on
    lower dimensional solution components.
    """
    from phcpy2c import py2c_copy_quaddobl_container_to_start_system
    from phcpy2c import py2c_copy_quaddobl_container_to_start_solutions
    from phcpy2c import py2c_quaddobl_cascade_homotopy
    from phcpy2c import py2c_solve_by_quaddobl_homotopy_continuation
    from phcpy2c import py2c_solcon_clear_quaddobl_solutions
    from phcpy2c import py2c_copy_quaddobl_target_solutions_to_container
    from interface import store_quaddobl_system
    from interface import store_quaddobl_solutions, load_quaddobl_solutions
    store_quaddobl_system(embsys)
    py2c_copy_quaddobl_container_to_start_system()
    store_quaddobl_solutions(len(embsys), esols)
    py2c_copy_quaddobl_container_to_start_solutions()
    py2c_quaddobl_cascade_homotopy()
    py2c_solve_by_quaddobl_homotopy_continuation()
    py2c_solcon_clear_quaddobl_solutions()
    py2c_copy_quaddobl_target_solutions_to_container()
    return load_quaddobl_solutions()

def cascade_step(embsys, esols, precision='d'):
    """
    Given in embsys an embedded polynomial system and
    solutions with nonzero slack variables in esols,
    does one step in the homotopy cascade, with precision
    'd'  : standard double precision (1.1e-15 or 2^(-53)),
    'dd' : double double precision (4.9e-32 or 2^(-104)),
    'qd' : quad double precision (1.2e-63 or 2^(-209)).
    The list on return contains witness points on
    lower dimensional solution components.
    """
    if(precision == 'd'):
        return standard_double_cascade_step(embsys, esols)
    elif(precision == 'dd'):
        return double_double_cascade_step(embsys, esols)
    elif(precision == 'qd'):
        return quad_double_cascade_step(embsys, esols)
    else:
        print 'wrong argument for precision'
        return None

def test_cascade():
    """
    Does one cascade step on simple example.
    In the top embedding we first find the 2-dimensional
    solution set x = 1.  In the cascade step we compute
    the three witness points on the twisted cubic.
    """
    from solver import solve
    pols = ['(x - 1)*(y-x^2);', \
            '(x - 1)*(z-x^3);', \
            '(x^2 - 1)*(y-x^2);' ]
    print pols
    embpols = embed(3, 2, pols)
    print 'the embedded system :'
    print embpols
    raw_input('hit enter to continue...')
    sols = solve(embpols, silent=True)
    for sol in sols:
        print sol
    print 'number of solutions :', len(sols)
    raw_input('hit enter to continue...')
    from solutions import filter_zero_coordinates, filter_regular
    sols0 = filter_zero_coordinates(sols, 'zz1', 1.0e-8, 'select')
    sols1 = filter_zero_coordinates(sols, 'zz1', 1.0e-8, 'remove')
    print 'solutions with zero slack variables :'
    for sol in sols0:
        print sol
    print 'solutions with nonzero slack variables :'
    for sol in sols1:
        print sol
    print len(sols), '=' , len(sols0), '+', len(sols1)
    rs1 = filter_regular(sols1, 1.0e-8, 'select')
    print 'number of nonsolutions :', len(rs1)
    raw_input('hit enter to continue...')
    print '... running cascade step ...'
    s2c = cascade_step(embpols, rs1, precision='d')
    print '... after running the cascade ...'
    for sol in s2c:
        print sol
    wp1 = drop_variable_from_polynomials(embpols, 'zz2')
    print 'the 1-dimensional embedding :'
    for pol in wp1:
        print pol
    ws1 = drop_coordinate_from_solutions(s2c, len(embpols), 'zz2')
    ws1f1 = filter_zero_coordinates(ws1, 'zz1', 1.0e-8, 'select')
    ws1f2 = filter_regular(ws1f1, 1.0e-8, 'select')
    print 'the witness points :'
    for sol in ws1f2:
        print sol

def decomposition(deg):
    """
    Returns the decomposition as a list of labels
    of witness points on the components.
    """
    from phcpy2c import py2c_factor_number_of_components
    from phcpy2c import py2c_factor_witness_points_of_component
    from phcpy2c import py2c_factor_trace_sum_difference
    nbcmp = py2c_factor_number_of_components()
    result = []
    for i in range(1, nbcmp+1):
        compnt = py2c_factor_witness_points_of_component(deg, i)
        tsd = py2c_factor_trace_sum_difference(deg, i, len(compnt), compnt)
        result.append((eval(compnt), tsd))
    return result

def monodromy_breakup(embsys, esols, dim):
    """
    Applies the monodromy breakup algorithm to factor
    the d-dimensional algebraic set represented by the
    embedded system e and its solutions esols.
    """
    from phcpy2c import py2c_factor_set_to_mute
    from phcpy2c import py2c_factor_assign_labels
    from phcpy2c import py2c_factor_initialize_monodromy
    from phcpy2c import py2c_factor_initialize_sampler
    from phcpy2c import py2c_factor_set_trace_slice
    from phcpy2c import py2c_factor_store_gammas
    from phcpy2c import py2c_factor_track_paths
    from phcpy2c import py2c_factor_store_solutions
    from phcpy2c import py2c_factor_restore_solutions
    from phcpy2c import py2c_factor_new_slices
    from phcpy2c import py2c_factor_swap_slices
    from phcpy2c import py2c_factor_permutation_after_loop
    from phcpy2c import py2c_factor_number_of_components
    from phcpy2c import py2c_factor_update_decomposition
    from phcpy2c import py2c_solcon_clear_standard_solutions
    from interface import store_standard_solutions
    print '... applying monodromy factorization ...'
    py2c_factor_set_to_mute()
    deg = len(esols)
    nvar = len(embsys)
    print 'dim =', dim
    store_standard_solutions(nvar, esols)
    py2c_factor_assign_labels(nvar, deg)
    # py2c_solcon_write_solutions()
    py2c_factor_initialize_sampler(dim)
    nbloops = input('give the maximum number of loops : ')
    py2c_factor_initialize_monodromy(nbloops, deg, dim)
    py2c_factor_store_solutions()
    print '... initializing the grid ...'
    for i in range(1, 3):
        py2c_factor_set_trace_slice(i)
        py2c_factor_store_gammas(nvar)
        py2c_factor_track_paths()
        py2c_factor_store_solutions()
        py2c_factor_restore_solutions()
        py2c_factor_swap_slices()
    for i in range(1, nbloops+1):
        print '... starting loop %d ...' % i
        py2c_factor_new_slices(dim, nvar)
        py2c_factor_store_gammas(nvar)
        py2c_factor_track_paths()
        py2c_solcon_clear_standard_solutions()
        py2c_factor_store_gammas(nvar)
        py2c_factor_track_paths()
        py2c_factor_store_solutions()
        sprm = py2c_factor_permutation_after_loop(deg)
        perm = eval(sprm)
        print 'the permutation :', perm
        nb0 = py2c_factor_number_of_components()
        done = py2c_factor_update_decomposition(deg, len(sprm), sprm)
        nb1 = py2c_factor_number_of_components()
        print 'number of factors : %d -> %d' % (nb0, nb1)
        print 'decomposition :', decomposition(deg)
        if(done == 1):
            break
        py2c_factor_restore_solutions()

def test_monodromy():
    """
    Runs a test on applying monodromy loops
    to factor a curve into irreducible components.
    """
    from solver import solve
    pols = ['(x^2 - y)*(x-y);', 'x^3 - z;']
    embsys = embed(3, 1, pols)
    # patch : make sure zz1 is last symbol!
    embsys[0] = 'x - x + y - y + z - z + ' + embsys[0]
    print embsys
    sols = solve(embsys, silent=True)
    # for sol in sols: print sol
    print 'the degree is', len(sols)
    monodromy_breakup(embsys, sols, 1)

def test():
    """
    Runs a test on algebraic sets.
    """
    from phcpy2c import py2c_set_seed
    py2c_set_seed(234798272)
    test_cascade()
    # test_monodromy()

if __name__ == "__main__":
    test()
