"""
Given two witness sets for two pure dimensional solution sets,
a diagonal homotopy computes a sets of witness sets for all components
of the intersection of the two pure dimensional solution sets.
"""

def top_diagonal_dimension(kdm, dim1, dim2):
    r"""
    Returns the number of slack variables at the top in the cascade of
    diagonal homotopies to intersect two sets of dimension *dim1* and *dim2*,
    where *dim1* >= *dim2* and *kdm* is the dimension before the embedding.
    Typically, *kdm* is the number of equations in the first witness set
    minus *dim1*.
    """
    if dim1 + dim2 < kdm:
        return dim2
    else:
        return kdm - dim1

def standard_diagonal_homotopy(dim1, sys1, esols1, dim2, sys2, esols2):
    r"""
    Defines a diagonal homotopy to intersect the witness sets defined
    by (*sys1*, *esols1*) and (*sys2*, *esols2*), respectively of dimensions
    *dim1* and *dim2*.  The systems *sys1* and *sys2* are assumed to be square
    and with as many slack variables as the dimension of the solution sets.
    The data is stored in standard double precision.
    """
    from phcpy.interface import store_standard_system as storesys
    from phcpy.interface import store_standard_solutions as storesols
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_target_solutions
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_start_system
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_start_solutions
    from phcpy.phcpy2c3 import py2c_standard_diagonal_homotopy
    from phcpy.phcpy2c3 import py2c_syscon_number_of_symbols
    from phcpy.phcpy2c3 import py2c_syscon_string_of_symbols
    from phcpy.phcpy2c3 import py2c_diagonal_symbols_doubler
    storesys(sys1)
    symbols = py2c_syscon_string_of_symbols()
    nbsymbs = py2c_syscon_number_of_symbols()
    print('number of symbols :', nbsymbs)
    print('names of variables :', symbols)
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
    r"""
    Defines a diagonal homotopy to intersect the witness sets defined
    by (*sys1*, *esols1*) and (*sys2*, *esols2*), respectively of dimensions
    *dim1* and *dim2*.  The systems *sys1* and *sys2* are assumed to be square
    and with as many slack variables as the dimension of the solution sets.
    The data is stored in double double precision.
    """
    from phcpy.interface import store_dobldobl_system as storesys
    from phcpy.interface import store_dobldobl_solutions as storesols
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_target_solutions
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_start_system
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_start_solutions
    from phcpy.phcpy2c3 import py2c_dobldobl_diagonal_homotopy
    from phcpy.phcpy2c3 import py2c_syscon_number_of_symbols
    from phcpy.phcpy2c3 import py2c_syscon_string_of_symbols
    from phcpy.phcpy2c3 import py2c_diagonal_symbols_doubler
    storesys(sys1)
    symbols = py2c_syscon_string_of_symbols()
    nbsymbs = py2c_syscon_number_of_symbols()
    print('number of symbols :', nbsymbs)
    print('names of variables :', symbols)
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
    r"""
    Defines a diagonal homotopy to intersect the witness sets defined
    by (*sys1*, *esols1*) and (*sys2*, *esols2*), respectively of dimensions
    *dim1* and *dim2*.  The systems *sys1* and *sys2* are assumed to be square
    and with as many slack variables as the dimension of the solution sets.
    The data is stored in quad double precision.
    """
    from phcpy.interface import store_quaddobl_system as storesys
    from phcpy.interface import store_quaddobl_solutions as storesols
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_target_solutions
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_start_system
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_start_solutions
    from phcpy.phcpy2c3 import py2c_quaddobl_diagonal_homotopy
    from phcpy.phcpy2c3 import py2c_syscon_number_of_symbols
    from phcpy.phcpy2c3 import py2c_syscon_string_of_symbols
    from phcpy.phcpy2c3 import py2c_diagonal_symbols_doubler
    storesys(sys1)
    symbols = py2c_syscon_string_of_symbols()
    nbsymbs = py2c_syscon_number_of_symbols()
    print('number of symbols :', nbsymbs)
    print('names of variables :', symbols)
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
    homotopy to intersect a set of dimension *dim1* with another set
    of dimension *dim2*, in standard double precision.  For this to work,
    standard_diagonal_homotopy must have been executed successfully.
    """
    from phcpy.phcpy2c3 import py2c_standard_diagonal_cascade_solutions
    if(dim1 >= dim2):
        py2c_standard_diagonal_cascade_solutions(dim1, dim2)
    else:
        py2c_standard_diagonal_cascade_solutions(dim2, dim1)

def dobldobl_diagonal_cascade_solutions(dim1, dim2):
    r"""
    Defines the start solutions in the cascade to start the diagonal
    homotopy to intersect a set of dimension *dim1* with another set
    of dimension *dim2*, in double double precision.  For this to work,
    dobldobl_diagonal_homotopy must have been executed successfully.
    """
    from phcpy.phcpy2c3 import py2c_dobldobl_diagonal_cascade_solutions
    if(dim1 >= dim2):
        py2c_dobldobl_diagonal_cascade_solutions(dim1, dim2)
    else:
        py2c_dobldobl_diagonal_cascade_solutions(dim2, dim1)

def quaddobl_diagonal_cascade_solutions(dim1, dim2):
    r"""
    Defines the start solutions in the cascade to start the diagonal
    homotopy to intersect a set of dimension *dim1* with another set
    of dimension *dim2*, in quad double precision.  For this to work,
    quaddobl_diagonal_homotopy must have been executed successfully.
    """
    from phcpy.phcpy2c3 import py2c_quaddobl_diagonal_cascade_solutions
    if(dim1 >= dim2):
        py2c_quaddobl_diagonal_cascade_solutions(dim1, dim2)
    else:
        py2c_quaddobl_diagonal_cascade_solutions(dim2, dim1)

def standard_start_diagonal_cascade(gamma=0, tasks=0):
    r"""
    Does the path tracking to start a diagonal cascade in standard double
    precision.  For this to work, the functions standard_diagonal_homotopy
    and standard_diagonal_cascade_solutions must be executed successfully.
    If *gamma* equals 0 on input, then a random gamma constant is generated,
    otherwise, the given complex gamma will be used in the homotopy.
    Multitasking is available, and activated by the *tasks* parameter.
    Returns the target (system and its corresponding) solutions.
    """
    from phcpy.phcpy2c3 import py2c_create_standard_homotopy
    from phcpy.phcpy2c3 import py2c_create_standard_homotopy_with_gamma
    from phcpy.phcpy2c3 import py2c_solve_by_standard_homotopy_continuation
    from phcpy.phcpy2c3 import py2c_solcon_clear_standard_solutions
    from phcpy.phcpy2c3 import py2c_syscon_clear_standard_system
    from phcpy.phcpy2c3 import py2c_copy_standard_target_solutions_to_container
    from phcpy.phcpy2c3 import py2c_copy_standard_target_system_to_container
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
    # from phcpy.phcpy2c3 import py2c_write_standard_target_system
    # print 'the standard target system :'
    # py2c_write_standard_target_system()
    py2c_copy_standard_target_system_to_container()
    tsys = load_standard_system()
    sols = load_standard_solutions()
    return (tsys, sols)

def dobldobl_start_diagonal_cascade(gamma=0, tasks=0):
    r"""
    Does the path tracking to start a diagonal cascade in double double
    precision.  For this to work, the functions dobldobl_diagonal_homotopy
    and dobldobl_diagonal_cascade_solutions must be executed successfully.
    If *gamma* equals 0 on input, then a random gamma constant is generated,
    otherwise, the given complex gamma will be used in the homotopy.
    Multitasking is available, and activated by the *tasks* parameter.
    Returns the target (system and its corresponding) solutions.
    """
    from phcpy.phcpy2c3 import py2c_create_dobldobl_homotopy
    from phcpy.phcpy2c3 import py2c_create_dobldobl_homotopy_with_gamma
    from phcpy.phcpy2c3 import py2c_solve_by_dobldobl_homotopy_continuation
    from phcpy.phcpy2c3 import py2c_solcon_clear_dobldobl_solutions
    from phcpy.phcpy2c3 import py2c_syscon_clear_dobldobl_system
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_target_solutions_to_container
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_target_system_to_container
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
    # from phcpy.phcpy2c3 import py2c_write_dobldobl_target_system
    # print 'the dobldobl target system :'
    # py2c_write_dobldobl_target_system()
    py2c_copy_dobldobl_target_system_to_container()
    tsys = load_dobldobl_system()
    sols = load_dobldobl_solutions()
    return (tsys, sols)

def quaddobl_start_diagonal_cascade(gamma=0, tasks=0):
    r"""
    Does the path tracking to start a diagonal cascade in quad double
    precision.  For this to work, the functions quaddobl_diagonal_homotopy
    and quaddobl_diagonal_cascade_solutions must be executed successfully.
    If *gamma* equals 0 on input, then a random gamma constant is generated,
    otherwise, the given complex gamma will be used in the homotopy.
    Multitasking is available, and is activated by the *tasks* parameter.
    Returns the target (system and its corresponding) solutions.
    """
    from phcpy.phcpy2c3 import py2c_create_quaddobl_homotopy
    from phcpy.phcpy2c3 import py2c_create_quaddobl_homotopy_with_gamma
    from phcpy.phcpy2c3 import py2c_solve_by_quaddobl_homotopy_continuation
    from phcpy.phcpy2c3 import py2c_solcon_clear_quaddobl_solutions
    from phcpy.phcpy2c3 import py2c_syscon_clear_quaddobl_system
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_target_solutions_to_container
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_target_system_to_container
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
    # from phcpy.phcpy2c3 import py2c_write_quaddobl_target_system
    # print 'the quaddobl target system :'
    # py2c_write_quaddobl_target_system()
    py2c_copy_quaddobl_target_system_to_container()
    tsys = load_quaddobl_system()
    sols = load_quaddobl_solutions()
    return (tsys, sols)

def standard_diagonal_solver(dim, dm1, sys1, sols1, dm2, sys2, sols2,
    tasks=0, verbose=True):
    r"""
    Runs the diagonal homotopies in standard double precision
    to intersect two witness sets stored in (*sys1*, *sols1*) and
    (*sys2*, *sols2*), of respective dimensions *dm1* and *dm2*.
    The ambient dimension equals *dim*.
    Multitasking is available, and is activated by the *tasks* parameter.
    Returns the last system in the cascade and its solutions.
    If *verbose*, then the solver runs in interactive mode, printing
    intermediate results to screen and prompting the user to continue.
    """
    from phcpy.phcpy2c3 import py2c_standard_collapse_diagonal
    from phcpy.interface import store_standard_solutions as storesols
    from phcpy.interface import load_standard_solutions as loadsols
    from phcpy.interface import load_standard_system as loadsys
    from phcpy.phcpy2c3 import py2c_extrinsic_top_diagonal_dimension
    from phcpy.solutions import filter_vanishing
    from phcpy.sets import drop_coordinate_from_standard_solutions
    from phcpy.sets import drop_variable_from_standard_polynomials
    from phcpy.cascades import standard_double_cascade_step
    topdim = py2c_extrinsic_top_diagonal_dimension(dim+dm1, dim+dm2, dm1, dm2)
    kdm = len(sys1) - dm1
    topdiagdim = top_diagonal_dimension(kdm, dm1, dm2)
    if verbose:
        print('the top dimension :', topdim, 'dim :', dim)
        print('number of slack variables at the top :', topdiagdim)
    standard_diagonal_homotopy(dm1, sys1, sols1, dm2, sys2, sols2)
    if verbose:
        print('defining the start solutions')
    standard_diagonal_cascade_solutions(dm1, dm2)
    if verbose:
        print('starting the diagonal cascade')
    (topsys, startsols) = standard_start_diagonal_cascade()
    if verbose:
        print('the system solved in the start of the cascade :')
        for pol in topsys:
            print(pol)
        print('the solutions after starting the diagonal cascade :')
        for sol in startsols:
            print(sol)
        input('hit enter to continue')
    for k in range(topdiagdim, 0, -1):
        endsols = standard_double_cascade_step(k, topsys, startsols)
        if verbose:
            print('after running cascade step %d :' % k)
            for sol in endsols:
                print(sol)
        endsolsf1 = filter_vanishing(endsols, 1.0e-8)
        if verbose:
            print('computed', len(endsolsf1), 'solutions')
            input('hit enter to continue')
        slack = 'zz' + str(k)
        nbvar = len(topsys)
        endsolsf2 = drop_coordinate_from_standard_solutions\
            (endsolsf1, nbvar, slack)
        if verbose:
            print('after dropping the slack coordinate from the solutions :')
            for sol in endsolsf2:
                print(sol)
            input('hit enter to continue')
        nextsys = drop_variable_from_standard_polynomials(topsys, slack)
        if verbose:
            print('after dropping the variable', slack, 'from the system :')
            for pol in nextsys:
                print(pol)
        (topsys, startsols) = (nextsys[:-1], endsolsf2)
    storesols(len(topsys), startsols)
    # py2c_standard_collapse_diagonal(topdim - 2*dim, 0)
    py2c_standard_collapse_diagonal(0, 0)
    result = (loadsys(), loadsols())
    return result

def dobldobl_diagonal_solver(dim, dm1, sys1, sols1, dm2, sys2, sols2, \
    tasks=0, verbose=True):
    r"""
    Runs the diagonal homotopies in double double precision
    to intersect two witness sets stored in (*sys1*, *sols1*) and
    (*sys2*, *sols2*), of respective dimensions *dm1* and *dm2*.
    The ambient dimension equals *dim*.
    Multitasking is available, and is activated by the *tasks* parameter.
    Returns the last system in the cascade and its solutions.
    If *verbose*, then the solver runs in interactive mode, printing
    intermediate results to screen and prompting the user to continue.
    """
    from phcpy.phcpy2c3 import py2c_dobldobl_collapse_diagonal
    from phcpy.interface import store_dobldobl_solutions as storesols
    from phcpy.interface import load_dobldobl_solutions as loadsols
    from phcpy.interface import load_dobldobl_system as loadsys
    from phcpy.phcpy2c3 import py2c_extrinsic_top_diagonal_dimension
    from phcpy.solutions import filter_vanishing
    from phcpy.sets import drop_coordinate_from_dobldobl_solutions
    from phcpy.sets import drop_variable_from_dobldobl_polynomials
    from phcpy.cascades import double_double_cascade_step
    topdim = py2c_extrinsic_top_diagonal_dimension(dim+dm1, dim+dm2, dm1, dm2)
    kdm = len(sys1) - dm1
    topdiagdim = top_diagonal_dimension(kdm, dm1, dm2)
    if verbose:
        print('the top dimension :', topdim, 'dim :', dim)
        print('number of slack variables at the top :', topdiagdim)
    dobldobl_diagonal_homotopy(dm1, sys1, sols1, dm2, sys2, sols2)
    if verbose:
        print('defining the start solutions')
    dobldobl_diagonal_cascade_solutions(dm1, dm2)
    if verbose:
        print('starting the diagonal cascade')
    (topsys, startsols) = dobldobl_start_diagonal_cascade()
    if verbose:
        print('the system solved in the start of the cascade :')
        for pol in topsys:
            print(pol)
        print('the solutions after starting the diagonal cascade :')
        for sol in startsols:
            print(sol)
        input('hit enter to continue')
    for k in range(topdiagdim, 0, -1):
        endsols = double_double_cascade_step(k, topsys, startsols)
        if verbose:
            print('after running cascade step %d :' % k)
            for sol in endsols:
                print(sol)
        endsolsf1 = filter_vanishing(endsols, 1.0e-8)
        if verbose:
            print('computed', len(endsolsf1), 'solutions')
            input('hit enter to continue')
        slack = 'zz' + str(k)
        nbvar = len(topsys)
        endsolsf2 = drop_coordinate_from_dobldobl_solutions\
            (endsolsf1, nbvar, slack)
        if verbose:
            print('after dropping the slack coordinate from the solutions :')
            for sol in endsolsf2:
                print(sol)
            input('hit enter to continue')
        nextsys = drop_variable_from_dobldobl_polynomials(topsys, slack)
        if verbose:
            print('after dropping the variable', slack, 'from the system :')
            for pol in nextsys:
                print(pol)
        (topsys, startsols) = (nextsys[:-1], endsolsf2)
    storesols(len(topsys), startsols)
    # py2c_dobldobl_collapse_diagonal(topdim - 2*dim, 0)
    py2c_dobldobl_collapse_diagonal(0, 0)
    result = (loadsys(), loadsols())
    return result

def quaddobl_diagonal_solver(dim, dm1, sys1, sols1, dm2, sys2, sols2, \
    tasks=0, verbose=True):
    r"""
    Runs the diagonal homotopies in quad double precision
    to intersect two witness sets stored in (*sys1*, *sols1*) and
    (*sys2*, *sols2*), of respective dimensions *dm1* and *dm2*.
    The ambient dimension equals *dim*.
    Multitasking is available, and is activated by the *tasks* parameter.
    Returns the last system in the cascade and its solutions.
    If *verbose*, then the solver runs in interactive mode, printing
    intermediate results to screen and prompting the user to continue.
    """
    from phcpy.phcpy2c3 import py2c_quaddobl_collapse_diagonal
    from phcpy.interface import store_quaddobl_solutions as storesols
    from phcpy.interface import load_quaddobl_solutions as loadsols
    from phcpy.interface import load_quaddobl_system as loadsys
    from phcpy.phcpy2c3 import py2c_extrinsic_top_diagonal_dimension
    from phcpy.solutions import filter_vanishing
    from phcpy.sets import drop_coordinate_from_quaddobl_solutions
    from phcpy.sets import drop_variable_from_quaddobl_polynomials
    from phcpy.cascades import quad_double_cascade_step
    topdim = py2c_extrinsic_top_diagonal_dimension(dim+dm1, dim+dm2, dm1, dm2)
    kdm = len(sys1) - dm1
    topdiagdim = top_diagonal_dimension(kdm, dm1, dm2)
    if verbose:
        print('the top dimension :', topdim, 'dim :', dim)
        print('number of slack variables at the top :', topdiagdim)
    quaddobl_diagonal_homotopy(dm1, sys1, sols1, dm2, sys2, sols2)
    if verbose:
        print('defining the start solutions')
    quaddobl_diagonal_cascade_solutions(dm1, dm2)
    if verbose:
        print('starting the diagonal cascade')
    (topsys, startsols) = quaddobl_start_diagonal_cascade()
    if verbose:
        for pol in topsys:
            print(pol)
        print('the solutions after starting the diagonal cascade :')
        for sol in startsols:
            print(sol)
        input('hit enter to continue')
    for k in range(topdiagdim, 0, -1):
        endsols = quad_double_cascade_step(k, topsys, startsols)
        if verbose:
            print('after running cascade step %d :' % k)
            for sol in endsols:
                print(sol)
        endsolsf1 = filter_vanishing(endsols, 1.0e-8)
        if verbose:
            print('computed', len(endsolsf1), 'solutions')
            input('hit enter to continue')
        slack = 'zz' + str(k)
        nbvar = len(topsys)
        endsolsf2 = drop_coordinate_from_quaddobl_solutions\
            (endsolsf1, nbvar, slack)
        if verbose:
            print('after dropping the slack coordinate from the solutions :')
            for sol in endsolsf2:
                print(sol)
            input('hit enter to continue')
        nextsys = drop_variable_from_quaddobl_polynomials(topsys, slack)
        if verbose:
            print('after dropping the variable', slack, 'from the system :')
            for pol in nextsys:
                print(pol)
        (topsys, startsols) = (nextsys[:-1], endsolsf2)
    storesols(len(topsys), startsols)
    # py2c_quaddobl_collapse_diagonal(topdim - 2*dim, 0)
    py2c_quaddobl_collapse_diagonal(0, 0)
    result = (loadsys(), loadsols())
    return result

def diagonal_solver(dim, dm1, sys1, sols1, dm2, sys2, sols2, tasks=0, \
    prc='d', verbose=True):
    r"""
    Runs the diagonal homotopies to intersect two witness sets stored in
    (*sys1*, *sols1*) and (*sys2*, *sols2*), of respective dimensions *dim1*
    and *dim2*.  The ambient dimension equals *dim*.
    Multitasking is available, and is activated by the *tasks* parameter.
    The precision is set by the parameter *prc*, which takes the default
    value 'd' for standard double, 'dd' for double double, or 'qd' for
    quad double precision.
    Returns the last system in the cascade and its solutions.
    """
    if(prc == 'd'):
        return standard_diagonal_solver\
                   (dim, dm1, sys1, sols1, dm2, sys2, sols2, tasks, verbose)
    elif(prc == 'dd'):
        return dobldobl_diagonal_solver\
                   (dim, dm1, sys1, sols1, dm2, sys2, sols2, tasks, verbose)
    elif(prc == 'qd'):
        return quaddobl_diagonal_solver\
                   (dim, dm1, sys1, sols1, dm2, sys2, sols2, tasks, verbose)
    else:
        print('wrong argument for precision')
        return None

def test_diaghom(precision='d'):
    """
    Test on the diagonal homotopy.
    """
    from phcpy.sets import witness_set_of_hypersurface
    hyp1 = 'x1*x2;'
    hyp2 = 'x1 - x2;'
    (w1sys, w1sols) = witness_set_of_hypersurface(2, hyp1, precision)
    print('the witness sets for', hyp1)
    for pol in w1sys:
        print(pol)
    for sol in w1sols:
        print(sol)
    (w2sys, w2sols) = witness_set_of_hypersurface(2, hyp2, precision)
    print('the witness sets for', hyp2)
    for pol in w2sys:
        print(pol)
    for sol in w2sols:
        print(sol)
    (sys, sols) = diagonal_solver\
        (2, 1, w1sys, w1sols, 1, w2sys, w2sols, 0, precision)
    print('the end system :')
    for pol in sys:
        print(pol)
    print('the solutions of the diagonal solver :')
    for sol in sols:
        print(sol)

def test():
    """
    Runs a test on algebraic sets.
    """
    from phcpy.phcpy2c3 import py2c_set_seed
    py2c_set_seed(234798272)
    test_diaghom('d')

if __name__ == "__main__":
    test()
