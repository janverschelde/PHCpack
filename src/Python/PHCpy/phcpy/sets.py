"""
This module exports routines of PHCpack to manipulate
positive dimensional solution sets of polynomial systems.
"""

def standard_embed(nvar, topdim, pols):
    """
    Given in pols a list of strings representing polynomials in nvar
    variables, with coefficients in standard double precision,
    this function returns an embedding of pols of dimension topdim.
    The topdim is the top dimension which equals the expected highest
    dimension of a component of the solution set of the system of polynomials.
    """
    from phcpy.phcpy2c import py2c_syscon_clear_standard_system
    from phcpy.phcpy2c \
    import py2c_syscon_initialize_number_of_standard_polynomials
    from phcpy.phcpy2c import py2c_syscon_store_standard_polynomial
    from phcpy.phcpy2c import py2c_syscon_load_standard_polynomial
    from phcpy.phcpy2c import py2c_embed_standard_system
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
    py2c_embed_standard_system(topdim)
    result = []
    for i in range(1, nbres+topdim+1):
        result.append(py2c_syscon_load_standard_polynomial(i))
    return result

def dobldobl_embed(nvar, topdim, pols):
    """
    Given in pols a list of strings that represent polynomials in nvar
    variables, with coefficients in double double precision,
    this function returns an embedding of pols of dimension topdim.
    The topdim is the top dimension which equals the expected highest
    dimension of a component of the solution set of the system of polynomials.
    """
    from phcpy.phcpy2c import py2c_syscon_clear_dobldobl_system
    from phcpy.phcpy2c \
    import py2c_syscon_initialize_number_of_dobldobl_polynomials
    from phcpy.phcpy2c import py2c_syscon_store_dobldobl_polynomial
    from phcpy.phcpy2c import py2c_syscon_load_dobldobl_polynomial
    from phcpy.phcpy2c import py2c_embed_dobldobl_system
    py2c_syscon_clear_dobldobl_system()
    nequ = len(pols)
    if nequ > nvar:
        py2c_syscon_initialize_number_of_dobldobl_polynomials(nequ)
        nbres = nequ
    else:
        py2c_syscon_initialize_number_of_dobldobl_polynomials(nvar)
        nbres = nvar
    for i in range(0, nequ):
        nchar = len(pols[i])
        py2c_syscon_store_dobldobl_polynomial(nchar, nvar, i+1, pols[i])
    py2c_embed_dobldobl_system(topdim)
    result = []
    for i in range(1, nbres+topdim+1):
        result.append(py2c_syscon_load_dobldobl_polynomial(i))
    return result

def quaddobl_embed(nvar, topdim, pols):
    """
    Given in pols a list of strings that represent polynomials in nvar
    variables, with coefficients in quad double precision,
    this function returns an embedding of pols of dimension topdim.
    The topdim is the top dimension which equals the expected highest
    dimension of a component of the solution set of the system of polynomials.
    """
    from phcpy.phcpy2c import py2c_syscon_clear_quaddobl_system
    from phcpy.phcpy2c \
    import py2c_syscon_initialize_number_of_quaddobl_polynomials
    from phcpy.phcpy2c import py2c_syscon_store_quaddobl_polynomial
    from phcpy.phcpy2c import py2c_syscon_load_quaddobl_polynomial
    from phcpy.phcpy2c import py2c_embed_quaddobl_system
    py2c_syscon_clear_quaddobl_system()
    nequ = len(pols)
    if nequ > nvar:
        py2c_syscon_initialize_number_of_quaddobl_polynomials(nequ)
        nbres = nequ
    else:
        py2c_syscon_initialize_number_of_quaddobl_polynomials(nvar)
        nbres = nvar
    for i in range(0, nequ):
        nchar = len(pols[i])
        py2c_syscon_store_quaddobl_polynomial(nchar, nvar, i+1, pols[i])
    py2c_embed_quaddobl_system(topdim)
    result = []
    for i in range(1, nbres+topdim+1):
        result.append(py2c_syscon_load_quaddobl_polynomial(i))
    return result

def embed(nvar, topdim, pols, precision='d'):
    """
    Given in pols a list of strings that represent polynomials in nvar
    variables, this function returns an embedding of pols of dimension topdim.
    The topdim is the top dimension which equals the expected highest
    dimension of a component of the solution set of the system of polynomials.
    The default precision of the coefficients is 'd', for standard double
    precision.  For double double and quad double precision, set the value
    of precision to 'dd' or 'qd' respectively.
    """
    if(precision == 'd'):
        return standard_embed(nvar, topdim, pols)
    elif(precision == 'dd'):
        return dobldobl_embed(nvar, topdim, pols)
    elif(precision == 'qd'):
        return quaddobl_embed(nvar, topdim, pols)
    else:
        print 'wrong argument for precision'
        return None

def witness_set_of_hypersurface(nvar, hpol, precision='d'):
    """
    Given in hpol the string representation of a polynomial
    in nvar variables (ending with a semicolon),
    on return is an embedded system and its solutions
    which represents a witness set for hpol.
    The number of solutions on return should equal
    the degree of the polynomial in hpol.
    Three different precisions are supported, by default double ('d'),
    or otherwise double double ('dd') or quad double ('qd').
    """
    if(precision == 'd'):
        from phcpy.phcpy2c import py2c_standard_witset_of_hypersurface
        from phcpy.interface import load_standard_system
        from phcpy.interface import load_standard_solutions
        py2c_standard_witset_of_hypersurface(nvar, len(hpol), hpol)
        return (load_standard_system(), load_standard_solutions())
    elif(precision == 'dd'):
        from phcpy.phcpy2c import py2c_dobldobl_witset_of_hypersurface
        from phcpy.interface import load_dobldobl_system
        from phcpy.interface import load_dobldobl_solutions
        py2c_dobldobl_witset_of_hypersurface(nvar, len(hpol), hpol)
        return (load_dobldobl_system(), load_dobldobl_solutions())
    elif(precision == 'qd'):
        from phcpy.phcpy2c import py2c_quaddobl_witset_of_hypersurface
        from phcpy.interface import load_quaddobl_system
        from phcpy.interface import load_quaddobl_solutions
        py2c_quaddobl_witset_of_hypersurface(nvar, len(hpol), hpol)
        return (load_quaddobl_system(), load_quaddobl_solutions())
    else:
        print 'wrong argument for precision'
        return None

def drop_variable_from_polynomials(pols, svar):
    """
    Removes the variable with symbol in the string svar
    from the list pols of strings that represented
    polynomials in several variables.
    """
    from phcpy.phcpy2c import py2c_syscon_standard_drop_variable_by_name
    from phcpy.phcpy2c import py2c_syscon_remove_symbol_name
    from phcpy.interface import store_standard_system, load_standard_system
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
    from phcpy.phcpy2c import py2c_syscon_clear_symbol_table
    from phcpy.phcpy2c import py2c_solcon_standard_drop_coordinate_by_name
    from phcpy.phcpy2c import py2c_syscon_remove_symbol_name
    from phcpy.interface import store_standard_solutions
    from phcpy.interface import load_standard_solutions
    py2c_syscon_clear_symbol_table()
    store_standard_solutions(nbvar, sols)
    py2c_solcon_standard_drop_coordinate_by_name(len(svar), svar)
    py2c_syscon_remove_symbol_name(len(svar), svar)
    return load_standard_solutions()

def standard_double_cascade_step(embsys, esols, tasks=0):
    """
    Given in embsys an embedded polynomial system and
    solutions with nonzero slace variables in esols,
    does one step in the homotopy cascade,
    with standard double precision arithmetic.
    The list on return contains witness points on
    lower dimensional solution components.
    """
    from phcpy.phcpy2c import py2c_copy_standard_container_to_start_system
    from phcpy.phcpy2c import py2c_copy_standard_container_to_start_solutions
    from phcpy.phcpy2c import py2c_standard_cascade_homotopy
    from phcpy.phcpy2c import py2c_solve_by_standard_homotopy_continuation
    from phcpy.phcpy2c import py2c_solcon_clear_standard_solutions
    from phcpy.phcpy2c import py2c_copy_standard_target_solutions_to_container
    from phcpy.interface import store_standard_system
    from phcpy.interface import store_standard_solutions
    from phcpy.interface import load_standard_solutions
    store_standard_system(embsys)
    py2c_copy_standard_container_to_start_system()
    store_standard_solutions(len(embsys), esols)
    py2c_copy_standard_container_to_start_solutions()
    py2c_standard_cascade_homotopy()
    py2c_solve_by_standard_homotopy_continuation(tasks)
    py2c_solcon_clear_standard_solutions()
    py2c_copy_standard_target_solutions_to_container()
    return load_standard_solutions()

def double_double_cascade_step(embsys, esols, tasks=0):
    """
    Given in embsys an embedded polynomial system and
    solutions with nonzero slace variables in esols,
    does one step in the homotopy cascade,
    with double double precision arithmetic.
    The list on return contains witness points on
    lower dimensional solution components.
    """
    from phcpy.phcpy2c import py2c_copy_dobldobl_container_to_start_system
    from phcpy.phcpy2c import py2c_copy_dobldobl_container_to_start_solutions
    from phcpy.phcpy2c import py2c_dobldobl_cascade_homotopy
    from phcpy.phcpy2c import py2c_solve_by_dobldobl_homotopy_continuation
    from phcpy.phcpy2c import py2c_solcon_clear_dobldobl_solutions
    from phcpy.phcpy2c import py2c_copy_dobldobl_target_solutions_to_container
    from phcpy.interface import store_dobldobl_system
    from phcpy.interface import store_dobldobl_solutions
    from phcpy.interface import load_dobldobl_solutions
    store_dobldobl_system(embsys)
    py2c_copy_dobldobl_container_to_start_system()
    store_dobldobl_solutions(len(embsys), esols)
    py2c_copy_dobldobl_container_to_start_solutions()
    py2c_dobldobl_cascade_homotopy()
    py2c_solve_by_dobldobl_homotopy_continuation(tasks)
    py2c_solcon_clear_dobldobl_solutions()
    py2c_copy_dobldobl_target_solutions_to_container()
    return load_dobldobl_solutions()

def quad_double_cascade_step(embsys, esols, tasks=0):
    """
    Given in embsys an embedded polynomial system and
    solutions with nonzero slace variables in esols,
    does one step in the homotopy cascade,
    with quad double precision arithmetic.
    The list on return contains witness points on
    lower dimensional solution components.
    """
    from phcpy.phcpy2c import py2c_copy_quaddobl_container_to_start_system
    from phcpy.phcpy2c import py2c_copy_quaddobl_container_to_start_solutions
    from phcpy.phcpy2c import py2c_quaddobl_cascade_homotopy
    from phcpy.phcpy2c import py2c_solve_by_quaddobl_homotopy_continuation
    from phcpy.phcpy2c import py2c_solcon_clear_quaddobl_solutions
    from phcpy.phcpy2c import py2c_copy_quaddobl_target_solutions_to_container
    from phcpy.interface import store_quaddobl_system
    from phcpy.interface import store_quaddobl_solutions
    from phcpy.interface import load_quaddobl_solutions
    store_quaddobl_system(embsys)
    py2c_copy_quaddobl_container_to_start_system()
    store_quaddobl_solutions(len(embsys), esols)
    py2c_copy_quaddobl_container_to_start_solutions()
    py2c_quaddobl_cascade_homotopy()
    py2c_solve_by_quaddobl_homotopy_continuation(tasks)
    py2c_solcon_clear_quaddobl_solutions()
    py2c_copy_quaddobl_target_solutions_to_container()
    return load_quaddobl_solutions()

def cascade_step(embsys, esols, precision='d', tasks=0):
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
        return standard_double_cascade_step(embsys, esols, tasks)
    elif(precision == 'dd'):
        return double_double_cascade_step(embsys, esols, tasks)
    elif(precision == 'qd'):
        return quad_double_cascade_step(embsys, esols, tasks)
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
    from phcpy.solver import solve
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
    from phcpy.solutions import filter_zero_coordinates, filter_regular
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
    from phcpy.phcpy2c import py2c_factor_number_of_components
    from phcpy.phcpy2c import py2c_factor_witness_points_of_component
    from phcpy.phcpy2c import py2c_factor_trace_sum_difference
    nbcmp = py2c_factor_number_of_components()
    result = []
    for i in range(1, nbcmp+1):
        compnt = py2c_factor_witness_points_of_component(deg, i)
        tsd = py2c_factor_trace_sum_difference(deg, i, len(compnt), compnt)
        result.append((eval(compnt), tsd))
    return result

def monodromy_breakup(embsys, esols, dim, verbose=True, nbloops=0):
    """
    Applies the monodromy breakup algorithm to factor
    the d-dimensional algebraic set represented by the
    embedded system e and its solutions esols.
    If verbose is False, then no output is written.
    If nbloops equals zero, then the user is prompted to give
    the maximum number of loops.
    """
    from phcpy.phcpy2c import py2c_factor_set_to_mute
    from phcpy.phcpy2c import py2c_factor_assign_labels
    from phcpy.phcpy2c import py2c_factor_initialize_monodromy
    from phcpy.phcpy2c import py2c_factor_initialize_sampler
    from phcpy.phcpy2c import py2c_factor_set_trace_slice
    from phcpy.phcpy2c import py2c_factor_store_gammas
    from phcpy.phcpy2c import py2c_factor_track_paths
    from phcpy.phcpy2c import py2c_factor_store_solutions
    from phcpy.phcpy2c import py2c_factor_restore_solutions
    from phcpy.phcpy2c import py2c_factor_new_slices
    from phcpy.phcpy2c import py2c_factor_swap_slices
    from phcpy.phcpy2c import py2c_factor_permutation_after_loop
    from phcpy.phcpy2c import py2c_factor_number_of_components
    from phcpy.phcpy2c import py2c_factor_update_decomposition
    from phcpy.phcpy2c import py2c_solcon_clear_standard_solutions
    from phcpy.interface import store_standard_solutions
    if(verbose):
        print '... applying monodromy factorization ...'
    py2c_factor_set_to_mute()
    deg = len(esols)
    nvar = len(embsys)
    if(verbose):
        print 'dim =', dim
    store_standard_solutions(nvar, esols)
    py2c_factor_assign_labels(nvar, deg)
    # py2c_solcon_write_solutions()
    py2c_factor_initialize_sampler(dim)
    if(nbloops == 0):
        strnbloops = raw_input('give the maximum number of loops : ')
        nbloops = int(strnbloops)
    py2c_factor_initialize_monodromy(nbloops, deg, dim)
    py2c_factor_store_solutions()
    if(verbose):
        print '... initializing the grid ...'
    for i in range(1, 3):
        py2c_factor_set_trace_slice(i)
        py2c_factor_store_gammas(nvar)
        py2c_factor_track_paths()
        py2c_factor_store_solutions()
        py2c_factor_restore_solutions()
        py2c_factor_swap_slices()
    for i in range(1, nbloops+1):
        if(verbose):
            print '... starting loop %d ...' % i
        py2c_factor_new_slices(dim, nvar)
        py2c_factor_store_gammas(nvar)
        py2c_factor_track_paths()
        py2c_solcon_clear_standard_solutions()
        py2c_factor_store_gammas(nvar)
        py2c_factor_track_paths()
        py2c_factor_store_solutions()
        sprm = py2c_factor_permutation_after_loop(deg)
        if(verbose):
            perm = eval(sprm)
            print 'the permutation :', perm
        nb0 = py2c_factor_number_of_components()
        done = py2c_factor_update_decomposition(deg, len(sprm), sprm)
        nb1 = py2c_factor_number_of_components()
        if(verbose):
            print 'number of factors : %d -> %d' % (nb0, nb1)
            print 'decomposition :', decomposition(deg)
        if(done == 1):
            break
        py2c_factor_restore_solutions()

def factor(dim, witsys, witsols, verbose=True, nbloops=20):
    """
    Applies monodromy to factor an equidimensional algebraic set,
    given as a witness sets, with the embedded polynomials in witsys,
    and corresponding generic points in witsols.
    The dimension of the algebraic set is given in dim.
    """
    monodromy_breakup(witsys, witsols, dim, verbose, nbloops)
    return decomposition(len(witsols))

def test_monodromy():
    """
    Runs a test on applying monodromy loops
    to factor a curve into irreducible components.
    """
    from phcpy.solver import solve
    pols = ['(x^2 - y)*(x-y);', 'x^3 - z;']
    embsys = embed(3, 1, pols)
    # patch : make sure zz1 is last symbol!
    embsys[0] = 'x - x + y - y + z - z + ' + embsys[0]
    print embsys
    sols = solve(embsys, silent=True)
    # for sol in sols: print sol
    print 'the degree is', len(sols)
    monodromy_breakup(embsys, sols, 1)

def test_factor():
    """
    Simple test on the factor method.
    """
    hyp = '(x+1)*(x^2 + y^2 + 1);'
    (wsys, wsols) = witness_set_of_hypersurface(2, hyp)
    fac = factor(1, wsys, wsols)
    print fac

def standard_diagonal_homotopy(dim1, sys1, esols1, dim2, sys2, esols2):
    """
    Defines a diagonal homotopy to intersect the witness sets defined
    by (sys1, esols1) and (sys2, esols2), respectively of dimensions
    dim1 and dim2.  The systems sys1 and sys2 are assumed to be square
    and with as many slack variables as the dimension of the solution sets.
    The data is stored in standard double precision.
    """
    from phcpy.interface import store_standard_system as storesys
    from phcpy.interface import store_standard_solutions as storesols
    from phcpy.phcpy2c import py2c_copy_standard_container_to_target_system
    from phcpy.phcpy2c import py2c_copy_standard_container_to_target_solutions
    from phcpy.phcpy2c import py2c_copy_standard_container_to_start_system
    from phcpy.phcpy2c import py2c_copy_standard_container_to_start_solutions
    from phcpy.phcpy2c import py2c_standard_diagonal_homotopy
    from phcpy.phcpy2c import py2c_syscon_number_of_symbols
    from phcpy.phcpy2c import py2c_syscon_string_of_symbols
    from phcpy.phcpy2c import py2c_diagonal_symbols_doubler
    storesys(sys1)
    symbols = py2c_syscon_string_of_symbols()
    nbsymbs = py2c_syscon_number_of_symbols()
    print 'number of symbols :', nbsymbs
    print 'names of variables :', symbols
    storesols(len(sys1), esols1)
    if(dim1 >= dim2):
        py2c_copy_standard_container_to_target_system()
        py2c_copy_standard_container_to_target_solutions()
    else:
        py2c_copy_standard_container_to_start_system()
        py2c_copy_standard_container_to_start_solutions()
    storesys(sys2)
    storesols(len(sys2), esols2)
    if(dim1 >= dim2):
        py2c_copy_standard_container_to_start_system()
        py2c_copy_standard_container_to_start_solutions()
    else:
        py2c_copy_standard_container_to_target_system()
        py2c_copy_standard_container_to_target_solutions()
    if(dim1 >= dim2):
        py2c_standard_diagonal_homotopy(dim1, dim2)
    else:
        py2c_standard_diagonal_homotopy(dim2, dim1)
    py2c_diagonal_symbols_doubler(nbsymbs-dim1, dim1, len(symbols), symbols)

def dobldobl_diagonal_homotopy(dim1, sys1, esols1, dim2, sys2, esols2):
    """
    Defines a diagonal homotopy to intersect the witness sets defined
    by (sys1, esols1) and (sys2, esols2), respectively of dimensions
    dim1 and dim2.  The systems sys1 and sys2 are assumed to be square
    and with as many slack variables as the dimension of the solution sets.
    The data is stored in double double precision.
    """
    from phcpy.interface import store_dobldobl_system as storesys
    from phcpy.interface import store_dobldobl_solutions as storesols
    from phcpy.phcpy2c import py2c_copy_dobldobl_container_to_target_system
    from phcpy.phcpy2c import py2c_copy_dobldobl_container_to_target_solutions
    from phcpy.phcpy2c import py2c_copy_dobldobl_container_to_start_system
    from phcpy.phcpy2c import py2c_copy_dobldobl_container_to_start_solutions
    from phcpy.phcpy2c import py2c_dobldobl_diagonal_homotopy
    from phcpy.phcpy2c import py2c_syscon_number_of_symbols
    from phcpy.phcpy2c import py2c_syscon_string_of_symbols
    from phcpy.phcpy2c import py2c_diagonal_symbols_doubler
    storesys(sys1)
    symbols = py2c_syscon_string_of_symbols()
    nbsymbs = py2c_syscon_number_of_symbols()
    print 'number of symbols :', nbsymbs
    print 'names of variables :', symbols
    storesols(len(sys1), esols1)
    if(dim1 >= dim2):
        py2c_copy_dobldobl_container_to_target_system()
        py2c_copy_dobldobl_container_to_target_solutions()
    else:
        py2c_copy_dobldobl_container_to_start_system()
        py2c_copy_dobldobl_container_to_start_solutions()
    storesys(sys2)
    storesols(len(sys2), esols2)
    if(dim1 >= dim2):
        py2c_copy_dobldobl_container_to_start_system()
        py2c_copy_dobldobl_container_to_start_solutions()
    else:
        py2c_copy_dobldobl_container_to_target_system()
        py2c_copy_dobldobl_container_to_target_solutions()
    if(dim1 >= dim2):
        py2c_dobldobl_diagonal_homotopy(dim1, dim2)
    else:
        py2c_dobldobl_diagonal_homotopy(dim2, dim1)
    py2c_diagonal_symbols_doubler(nbsymbs-dim1, dim1, len(symbols), symbols)

def quaddobl_diagonal_homotopy(dim1, sys1, esols1, dim2, sys2, esols2):
    """
    Defines a diagonal homotopy to intersect the witness sets defined
    by (sys1, esols1) and (sys2, esols2), respectively of dimensions
    dim1 and dim2.  The systems sys1 and sys2 are assumed to be square
    and with as many slack variables as the dimension of the solution sets.
    The data is stored in quad double precision.
    """
    from phcpy.interface import store_quaddobl_system as storesys
    from phcpy.interface import store_quaddobl_solutions as storesols
    from phcpy.phcpy2c import py2c_copy_quaddobl_container_to_target_system
    from phcpy.phcpy2c import py2c_copy_quaddobl_container_to_target_solutions
    from phcpy.phcpy2c import py2c_copy_quaddobl_container_to_start_system
    from phcpy.phcpy2c import py2c_copy_quaddobl_container_to_start_solutions
    from phcpy.phcpy2c import py2c_quaddobl_diagonal_homotopy
    from phcpy.phcpy2c import py2c_syscon_number_of_symbols
    from phcpy.phcpy2c import py2c_syscon_string_of_symbols
    from phcpy.phcpy2c import py2c_diagonal_symbols_doubler
    storesys(sys1)
    symbols = py2c_syscon_string_of_symbols()
    nbsymbs = py2c_syscon_number_of_symbols()
    print 'number of symbols :', nbsymbs
    print 'names of variables :', symbols
    storesols(len(sys1), esols1)
    if(dim1 >= dim2):
        py2c_copy_quaddobl_container_to_target_system()
        py2c_copy_quaddobl_container_to_target_solutions()
    else:
        py2c_copy_quaddobl_container_to_start_system()
        py2c_copy_quaddobl_container_to_start_solutions()
    storesys(sys2)
    storesols(len(sys2), esols2)
    if(dim1 >= dim2):
        py2c_copy_quaddobl_container_to_start_system()
        py2c_copy_quaddobl_container_to_start_solutions()
    else:
        py2c_copy_quaddobl_container_to_target_system()
        py2c_copy_quaddobl_container_to_target_solutions()
    if(dim1 >= dim2):
        py2c_quaddobl_diagonal_homotopy(dim1, dim2)
    else:
        py2c_quaddobl_diagonal_homotopy(dim2, dim1)
    py2c_diagonal_symbols_doubler(nbsymbs-dim1, dim1, len(symbols), symbols)

def standard_diagonal_cascade_solutions(dim1, dim2):
    """
    Defines the start solutions in the cascade to start the diagonal
    homotopy to intersect a set of dimension dim1 with another set
    of dimension dim2, in standard double precision.  For this to work,
    standard_diagonal_homotopy must have been executed successfully.
    """
    from phcpy.phcpy2c import py2c_standard_diagonal_cascade_solutions
    if(dim1 >= dim2):
        py2c_standard_diagonal_cascade_solutions(dim1, dim2)
    else:
        py2c_standard_diagonal_cascade_solutions(dim2, dim1)

def dobldobl_diagonal_cascade_solutions(dim1, dim2):
    """
    Defines the start solutions in the cascade to start the diagonal
    homotopy to intersect a set of dimension dim1 with another set
    of dimension dim2, in double double precision.  For this to work,
    dobldobl_diagonal_homotopy must have been executed successfully.
    """
    from phcpy.phcpy2c import py2c_dobldobl_diagonal_cascade_solutions
    if(dim1 >= dim2):
        py2c_dobldobl_diagonal_cascade_solutions(dim1, dim2)
    else:
        py2c_dobldobl_diagonal_cascade_solutions(dim2, dim1)

def quaddobl_diagonal_cascade_solutions(dim1, dim2):
    """
    Defines the start solutions in the cascade to start the diagonal
    homotopy to intersect a set of dimension dim1 with another set
    of dimension dim2, in quad double precision.  For this to work,
    quaddobl_diagonal_homotopy must have been executed successfully.
    """
    from phcpy.phcpy2c import py2c_quaddobl_diagonal_cascade_solutions
    if(dim1 >= dim2):
        py2c_quaddobl_diagonal_cascade_solutions(dim1, dim2)
    else:
        py2c_quaddobl_diagonal_cascade_solutions(dim2, dim1)

def standard_start_diagonal_cascade(gamma=0, tasks=0):
    """
    Does the path tracking to start a diagonal cascade in standard double
    precision.  For this to work, the functions standard_diagonal_homotopy
    and standard_diagonal_cascade_solutions must be executed successfully.
    If gamma equals 0 on input, then a random gamma constant is generated,
    otherwise, the given complex gamma will be used in the homotopy.
    Multitasking is available, and activated by the tasks parameter.
    Returns the target (system and its corresponding) solutions.
    """
    from phcpy.phcpy2c import py2c_create_standard_homotopy
    from phcpy.phcpy2c import py2c_create_standard_homotopy_with_gamma
    from phcpy.phcpy2c import py2c_solve_by_standard_homotopy_continuation
    from phcpy.phcpy2c import py2c_solcon_clear_standard_solutions
    from phcpy.phcpy2c import py2c_syscon_clear_standard_system
    from phcpy.phcpy2c import py2c_copy_standard_target_solutions_to_container
    from phcpy.phcpy2c import py2c_copy_standard_target_system_to_container
    from phcpy.interface import load_standard_solutions
    from phcpy.interface import load_standard_system
    if(gamma == 0):
        py2c_create_standard_homotopy()
    else:
        py2c_create_standard_homotopy_with_gamma(gamma.real, gamma.imag)
    py2c_solve_by_standard_homotopy_continuation(tasks)
    py2c_solcon_clear_standard_solutions()
    py2c_syscon_clear_standard_system()
    py2c_copy_standard_target_solutions_to_container()
    from phcpy.phcpy2c import py2c_write_standard_target_system
    # print 'the standard target system :'
    # py2c_write_standard_target_system()
    py2c_copy_standard_target_system_to_container()
    tsys = load_standard_system()
    sols = load_standard_solutions()
    return (tsys, sols)

def dobldobl_start_diagonal_cascade(gamma=0, tasks=0):
    """
    Does the path tracking to start a diagonal cascade in double double
    precision.  For this to work, the functions dobldobl_diagonal_homotopy
    and dobldobl_diagonal_cascade_solutions must be executed successfully.
    If gamma equals 0 on input, then a random gamma constant is generated,
    otherwise, the given complex gamma will be used in the homotopy.
    Multitasking is available, and activated by the tasks parameter.
    Returns the target (system and its corresponding) solutions.
    """
    from phcpy.phcpy2c import py2c_create_dobldobl_homotopy
    from phcpy.phcpy2c import py2c_create_dobldobl_homotopy_with_gamma
    from phcpy.phcpy2c import py2c_solve_by_dobldobl_homotopy_continuation
    from phcpy.phcpy2c import py2c_solcon_clear_dobldobl_solutions
    from phcpy.phcpy2c import py2c_syscon_clear_dobldobl_system
    from phcpy.phcpy2c import py2c_copy_dobldobl_target_solutions_to_container
    from phcpy.phcpy2c import py2c_copy_dobldobl_target_system_to_container
    from phcpy.interface import load_dobldobl_solutions
    from phcpy.interface import load_dobldobl_system
    if(gamma == 0):
        py2c_create_dobldobl_homotopy()
    else:
        py2c_create_dobldobl_homotopy_with_gamma(gamma.real, gamma.imag)
    py2c_solve_by_dobldobl_homotopy_continuation(tasks)
    py2c_solcon_clear_dobldobl_solutions()
    py2c_syscon_clear_dobldobl_system()
    py2c_copy_dobldobl_target_solutions_to_container()
    from phcpy.phcpy2c import py2c_write_dobldobl_target_system
    # print 'the dobldobl target system :'
    # py2c_write_dobldobl_target_system()
    py2c_copy_dobldobl_target_system_to_container()
    tsys = load_dobldobl_system()
    sols = load_dobldobl_solutions()
    return (tsys, sols)

def quaddobl_start_diagonal_cascade(gamma=0, tasks=0):
    """
    Does the path tracking to start a diagonal cascade in quad double
    precision.  For this to work, the functions quaddobl_diagonal_homotopy
    and quaddobl_diagonal_cascade_solutions must be executed successfully.
    If gamma equals 0 on input, then a random gamma constant is generated,
    otherwise, the given complex gamma will be used in the homotopy.
    Multitasking is available, and is activated by the tasks parameter.
    Returns the target (system and its corresponding) solutions.
    """
    from phcpy.phcpy2c import py2c_create_quaddobl_homotopy
    from phcpy.phcpy2c import py2c_create_quaddobl_homotopy_with_gamma
    from phcpy.phcpy2c import py2c_solve_by_quaddobl_homotopy_continuation
    from phcpy.phcpy2c import py2c_solcon_clear_quaddobl_solutions
    from phcpy.phcpy2c import py2c_syscon_clear_quaddobl_system
    from phcpy.phcpy2c import py2c_copy_quaddobl_target_solutions_to_container
    from phcpy.phcpy2c import py2c_copy_quaddobl_target_system_to_container
    from phcpy.interface import load_quaddobl_solutions
    from phcpy.interface import load_quaddobl_system
    if(gamma == 0):
        py2c_create_quaddobl_homotopy()
    else:
        py2c_create_quaddobl_homotopy_with_gamma(gamma.real, gamma.imag)
    py2c_solve_by_quaddobl_homotopy_continuation(tasks)
    py2c_solcon_clear_quaddobl_solutions()
    py2c_syscon_clear_quaddobl_system()
    py2c_copy_quaddobl_target_solutions_to_container()
    from phcpy.phcpy2c import py2c_write_quaddobl_target_system
    # print 'the quaddobl target system :'
    # py2c_write_quaddobl_target_system()
    py2c_copy_quaddobl_target_system_to_container()
    tsys = load_quaddobl_system()
    sols = load_quaddobl_solutions()
    return (tsys, sols)

def standard_diagonal_solver(dim, dm1, sys1, sols1, dm2, sys2, sols2, tasks=0):
    """
    Runs the diagonal homotopies in standard double precision
    to intersect two witness sets stored in (sys1, sols1) and
    (sys2, sols2), of respective dimensions dm1 and dm2.
    The ambient dimension equals dim.
    Multitasking is available, and is activated by the tasks parameter.
    Returns the last system in the cascade and its solutions.
    """
    from phcpy.phcpy2c import py2c_standard_collapse_diagonal
    from phcpy.interface import store_standard_solutions as storesols
    from phcpy.interface import load_standard_solutions as loadsols
    from phcpy.interface import load_standard_system as loadsys
    from phcpy.phcpy2c import py2c_extrinsic_top_diagonal_dimension
    topdim = py2c_extrinsic_top_diagonal_dimension(dim+dm1, dim+dm2, dm1, dm2)
    print 'the top dimension :', topdim
    standard_diagonal_homotopy(dm1, sys1, sols1, dm2, sys2, sols2)
    print 'defining the start solutions'
    standard_diagonal_cascade_solutions(dm1, dm2)
    print 'starting the diagonal cascade'
    (topsys, startsols) = standard_start_diagonal_cascade()
    print 'the system solved in the start of the cascade :'
    for pol in topsys:
        print pol
    print 'the solutions after starting the diagonal cascade :'
    for sol in startsols:
        print sol
    endsols = standard_double_cascade_step(topsys, startsols)
    print 'after running one cascade step :'
    for sol in endsols:
        print sol
    storesols(len(topsys), endsols)
    py2c_standard_collapse_diagonal(topdim - 2*dim, 0)
    result = (loadsys(), loadsols())
    return result

def dobldobl_diagonal_solver(dim, dm1, sys1, sols1, dm2, sys2, sols2, tasks=0):
    """
    Runs the diagonal homotopies in double double precision
    to intersect two witness sets stored in (sys1, sols1) and
    (sys2, sols2), of respective dimensions dm1 and dm2.
    The ambient dimension equals dim.
    Multitasking is available, and is activated by the tasks parameter.
    Returns the last system in the cascade and its solutions.
    """
    from phcpy.phcpy2c import py2c_dobldobl_collapse_diagonal
    from phcpy.interface import store_dobldobl_solutions as storesols
    from phcpy.interface import load_dobldobl_solutions as loadsols
    from phcpy.interface import load_dobldobl_system as loadsys
    from phcpy.phcpy2c import py2c_extrinsic_top_diagonal_dimension
    topdim = py2c_extrinsic_top_diagonal_dimension(dim+dm1, dim+dm2, dm1, dm2)
    print 'the top dimension :', topdim
    dobldobl_diagonal_homotopy(dm1, sys1, sols1, dm2, sys2, sols2)
    print 'defining the start solutions'
    dobldobl_diagonal_cascade_solutions(dm1, dm2)
    print 'starting the diagonal cascade'
    (topsys, startsols) = dobldobl_start_diagonal_cascade()
    print 'the system solved in the start of the cascade :'
    for pol in topsys:
        print pol
    print 'the solutions after starting the diagonal cascade :'
    for sol in startsols:
        print sol
    endsols = double_double_cascade_step(topsys, startsols)
    print 'after running one cascade step :'
    for sol in endsols:
        print sol
    storesols(len(topsys), endsols)
    py2c_dobldobl_collapse_diagonal(topdim - 2*dim, 0)
    result = (loadsys(), loadsols())
    return result

def quaddobl_diagonal_solver(dim, dm1, sys1, sols1, dm2, sys2, sols2, tasks=0):
    """
    Runs the diagonal homotopies in quad double precision
    to intersect two witness sets stored in (sys1, sols1) and
    (sys2, sols2), of respective dimensions dm1 and dm2.
    The ambient dimension equals dim.
    Multitasking is available, and is activated by the tasks parameter.
    Returns the last system in the cascade and its solutions.
    """
    from phcpy.phcpy2c import py2c_quaddobl_collapse_diagonal
    from phcpy.interface import store_quaddobl_solutions as storesols
    from phcpy.interface import load_quaddobl_solutions as loadsols
    from phcpy.interface import load_quaddobl_system as loadsys
    from phcpy.phcpy2c import py2c_extrinsic_top_diagonal_dimension
    topdim = py2c_extrinsic_top_diagonal_dimension(dim+dm1, dim+dm2, dm1, dm2)
    print 'the top dimension :', topdim
    quaddobl_diagonal_homotopy(dm1, sys1, sols1, dm2, sys2, sols2)
    print 'defining the start solutions'
    quaddobl_diagonal_cascade_solutions(dm1, dm2)
    print 'starting the diagonal cascade'
    (topsys, startsols) = quaddobl_start_diagonal_cascade()
    print 'the system solved in the start of the cascade :'
    for pol in topsys:
        print pol
    print 'the solutions after starting the diagonal cascade :'
    for sol in startsols:
        print sol
    endsols = quad_double_cascade_step(topsys, startsols)
    print 'after running one cascade step :'
    for sol in endsols:
        print sol
    storesols(len(topsys), endsols)
    py2c_quaddobl_collapse_diagonal(topdim - 2*dim, 0)
    result = (loadsys(), loadsols())
    return result

def diagonal_solver(dim, dm1, sys1, sols1, dm2, sys2, sols2, tasks=0, prc='d'):
    """
    Runs the diagonal homotopies to intersect two witness sets stored in
    (sys1, sols1) and (sys2, sols2), of respective dimensions dim1 and dim2.
    The ambient dimension equals dim.
    Multitasking is available, and is activated by the tasks parameter.
    The precision is set by the parameter prc, which takes the default
    value 'd' for standard double, 'dd' for double double, or 'qd' for
    quad double precision.
    Returns the last system in the cascade and its solutions.
    """
    if(prc == 'd'):
        return standard_diagonal_solver\
                   (dim, dm1, sys1, sols1, dm2, sys2, sols2, tasks)
    elif(prc == 'dd'):
        return dobldobl_diagonal_solver\
                   (dim, dm1, sys1, sols1, dm2, sys2, sols2, tasks)
    elif(prc == 'qd'):
        return quaddobl_diagonal_solver\
                   (dim, dm1, sys1, sols1, dm2, sys2, sols2, tasks)
    else:
        print 'wrong argument for precision'
        return None

def test_diaghom(precision='d'):
    """
    Test on the diagonal homotopy.
    """
    hyp1 = 'x1*x2;'
    hyp2 = 'x1 - x2;'
    (w1sys, w1sols) = witness_set_of_hypersurface(2, hyp1, precision)
    print 'the witness sets for', hyp1
    for pol in w1sys:
        print pol
    for sol in w1sols:
        print sol
    (w2sys, w2sols) = witness_set_of_hypersurface(2, hyp2, precision)
    print 'the witness sets for', hyp2
    for pol in w2sys:
        print pol
    for sol in w2sols:
        print sol
    (sys, sols) = diagonal_solver\
        (2, 1, w1sys, w1sols, 1, w2sys, w2sols, 0, precision)
    print 'the end system :'
    for pol in sys:
        print pol
    print 'the solutions of the diagonal solver :'
    for sol in sols:
        print sol

def test():
    """
    Runs a test on algebraic sets.
    """
    from phcpy.phcpy2c import py2c_set_seed
    py2c_set_seed(234798272)
    # test_cascade()
    # test_monodromy()
    test_diaghom('d')

if __name__ == "__main__":
    test()
