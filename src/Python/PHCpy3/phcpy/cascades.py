"""
A cascade homotopy removes one hyperplane from an embedded system,
taking the solutions with nonzero slack variables to solutions on
lower dimensional components of the solution set of the original system.
"""

def standard_double_cascade_step(dim, embsys, esols, tasks=0):
    r"""
    Given in *embsys* an embedded polynomial system and
    solutions with nonzero slack variables in *esols*, does one step 
    in the homotopy cascade, with standard double precision arithmetic.
    The dimension of the solution set represented by *embsys*
    and *esols* is the value of *dim*.
    The number of tasks in multithreaded path tracking is given by *tasks*.
    The default zero value of *tasks* indicates no multithreading.
    The list on return contains witness points on
    lower dimensional solution components.
    """
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_start_system
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_start_solutions
    from phcpy.phcpy2c3 import py2c_standard_cascade_homotopy
    from phcpy.phcpy2c3 import py2c_solve_by_standard_homotopy_continuation
    from phcpy.phcpy2c3 import py2c_solcon_clear_standard_solutions
    from phcpy.phcpy2c3 import py2c_copy_standard_target_solutions_to_container
    from phcpy.interface import store_standard_witness_set
    from phcpy.interface import load_standard_solutions
    store_standard_witness_set(len(embsys), dim, embsys, esols)
    py2c_copy_standard_container_to_start_system()
    py2c_copy_standard_container_to_start_solutions()
    py2c_standard_cascade_homotopy()
    py2c_solve_by_standard_homotopy_continuation(tasks)
    py2c_solcon_clear_standard_solutions()
    py2c_copy_standard_target_solutions_to_container()
    return load_standard_solutions()

def double_double_cascade_step(dim, embsys, esols, tasks=0):
    r"""
    Given in *embsys* an embedded polynomial system and
    solutions with nonzero slack variables in *esols*, does one step 
    in the homotopy cascade, with double double precision arithmetic.
    The dimension of the solution set represented by *embsys*
    and *esols* is the value of *dim*.
    The number of tasks in multithreaded path tracking is given by *tasks*.
    The default zero value of *tasks* indicates no multithreading.
    The list on return contains witness points on
    lower dimensional solution components.
    """
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_start_system
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_start_solutions
    from phcpy.phcpy2c3 import py2c_dobldobl_cascade_homotopy
    from phcpy.phcpy2c3 import py2c_solve_by_dobldobl_homotopy_continuation
    from phcpy.phcpy2c3 import py2c_solcon_clear_dobldobl_solutions
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_target_solutions_to_container
    from phcpy.interface import store_dobldobl_witness_set
    from phcpy.interface import load_dobldobl_solutions
    store_dobldobl_witness_set(len(embsys), dim, embsys, esols)
    py2c_copy_dobldobl_container_to_start_system()
    py2c_copy_dobldobl_container_to_start_solutions()
    py2c_dobldobl_cascade_homotopy()
    py2c_solve_by_dobldobl_homotopy_continuation(tasks)
    py2c_solcon_clear_dobldobl_solutions()
    py2c_copy_dobldobl_target_solutions_to_container()
    return load_dobldobl_solutions()

def quad_double_cascade_step(dim, embsys, esols, tasks=0):
    r"""
    Given in *embsys* an embedded polynomial system and
    solutions with nonzero slack variables in *esols*, does one step 
    in the homotopy cascade, with quad double precision arithmetic.
    The dimension of the solution set represented by *embsys*
    and *esols* is the value of *dim*.
    The number of tasks in multithreaded path tracking is given by *tasks*.
    The default zero value of *tasks* indicates no multithreading.
    The list on return contains witness points on
    lower dimensional solution components.
    """
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_start_system
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_start_solutions
    from phcpy.phcpy2c3 import py2c_quaddobl_cascade_homotopy
    from phcpy.phcpy2c3 import py2c_solve_by_quaddobl_homotopy_continuation
    from phcpy.phcpy2c3 import py2c_solcon_clear_quaddobl_solutions
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_target_solutions_to_container
    from phcpy.interface import store_quaddobl_witness_set
    from phcpy.interface import load_quaddobl_solutions
    store_quaddobl_witness_set(len(embsys), dim, embsys, esols)
    py2c_copy_quaddobl_container_to_start_system()
    py2c_copy_quaddobl_container_to_start_solutions()
    py2c_quaddobl_cascade_homotopy()
    py2c_solve_by_quaddobl_homotopy_continuation(tasks)
    py2c_solcon_clear_quaddobl_solutions()
    py2c_copy_quaddobl_target_solutions_to_container()
    return load_quaddobl_solutions()

def cascade_step(dim, embsys, esols, precision='d', tasks=0):
    r"""
    Given in *embsys* an embedded polynomial system and
    solutions with nonzero slack variables in *esols*,
    does one step in the homotopy cascade, with *precision*

    *d*: standard double precision (1.1e-15 or 2^(-53)),

    *dd*: double double precision (4.9e-32 or 2^(-104)),

    *qd*: quad double precision (1.2e-63 or 2^(-209)).

    The dimension of the solution set represented by *embsys*
    and *esols* is the value of *dim*.
    The number of tasks in multithreaded path tracking is given by *tasks*.
    The default zero value of *tasks* indicates no multithreading.
    The list on return contains witness points on
    lower dimensional solution components.
    """
    if(precision == 'd'):
        return standard_double_cascade_step(dim, embsys, esols, tasks)
    elif(precision == 'dd'):
        return double_double_cascade_step(dim, embsys, esols, tasks)
    elif(precision == 'qd'):
        return quad_double_cascade_step(dim, embsys, esols, tasks)
    else:
        print('wrong argument for precision')
        return None

def standard_double_laurent_cascade_step(dim, embsys, esols, tasks=0):
    r"""
    Given in *embsys* an embedded Laurent polynomial system and
    solutions with nonzero slack variables in *esols*, does one step
    in the homotopy cascade, with standard double precision arithmetic.
    The dimension of the solution set represented by *embsys*
    and *esols* is the value of *dim*.
    The number of tasks in multithreaded path tracking is given by *tasks*.
    The default zero value of *tasks* indicates no multithreading.
    The list on return contains witness points on
    lower dimensional solution components.
    """
    from phcpy.phcpy2c3 \
    import py2c_copy_standard_Laurent_container_to_start_system
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_start_solutions
    from phcpy.phcpy2c3 import py2c_standard_Laurent_cascade_homotopy
    from phcpy.phcpy2c3 \
    import py2c_solve_by_standard_Laurent_homotopy_continuation
    from phcpy.phcpy2c3 import py2c_solcon_clear_standard_solutions
    from phcpy.phcpy2c3 import py2c_copy_standard_target_solutions_to_container
    from phcpy.interface import store_standard_laurent_witness_set
    from phcpy.interface import load_standard_solutions
    store_standard_laurent_witness_set(len(embsys), dim, embsys, esols)
    py2c_copy_standard_Laurent_container_to_start_system()
    py2c_copy_standard_container_to_start_solutions()
    py2c_standard_Laurent_cascade_homotopy()
    py2c_solve_by_standard_Laurent_homotopy_continuation(tasks)
    py2c_solcon_clear_standard_solutions()
    py2c_copy_standard_target_solutions_to_container()
    return load_standard_solutions()

def double_double_laurent_cascade_step(dim, embsys, esols, tasks=0):
    r"""
    Given in *embsys* an embedded Laurent polynomial system and
    solutions with nonzero slack variables in *esols*, does one step
    in the homotopy cascade, with double double precision arithmetic.
    The dimension of the solution set represented by *embsys*
    and *esols* is the value of *dim*.
    The number of tasks in multithreaded path tracking is given by *tasks*.
    The default zero value of *tasks* indicates no multithreading.
    The list on return contains witness points on
    lower dimensional solution components.
    """
    from phcpy.phcpy2c3 \
    import py2c_copy_dobldobl_Laurent_container_to_start_system
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_start_solutions
    from phcpy.phcpy2c3 import py2c_dobldobl_Laurent_cascade_homotopy
    from phcpy.phcpy2c3 \
    import py2c_solve_by_dobldobl_Laurent_homotopy_continuation
    from phcpy.phcpy2c3 import py2c_solcon_clear_dobldobl_solutions
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_target_solutions_to_container
    from phcpy.interface import store_dobldobl_laurent_witness_set
    from phcpy.interface import load_dobldobl_solutions
    store_dobldobl_laurent_witness_set(len(embsys), dim, embsys, esols)
    py2c_copy_dobldobl_Laurent_container_to_start_system()
    py2c_copy_dobldobl_container_to_start_solutions()
    py2c_dobldobl_Laurent_cascade_homotopy()
    py2c_solve_by_dobldobl_Laurent_homotopy_continuation(tasks)
    py2c_solcon_clear_dobldobl_solutions()
    py2c_copy_dobldobl_target_solutions_to_container()
    return load_dobldobl_solutions()

def quad_double_laurent_cascade_step(dim, embsys, esols, tasks=0):
    r"""
    Given in *embsys* an embedded Laurent polynomial system and
    solutions with nonzero slack variables in *esols*, does one step
    in the homotopy cascade, with quad double precision arithmetic.
    The dimension of the solution set represented by *embsys*
    and *esols* is the value of *dim*.
    The number of tasks in multithreaded path tracking is given by *tasks*.
    The default zero value of *tasks* indicates no multithreading.
    The list on return contains witness points on
    lower dimensional solution components.
    """
    from phcpy.phcpy2c3 \
    import py2c_copy_quaddobl_Laurent_container_to_start_system
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_start_solutions
    from phcpy.phcpy2c3 import py2c_quaddobl_Laurent_cascade_homotopy
    from phcpy.phcpy2c3 \
    import py2c_solve_by_quaddobl_Laurent_homotopy_continuation
    from phcpy.phcpy2c3 import py2c_solcon_clear_quaddobl_solutions
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_target_solutions_to_container
    from phcpy.interface import store_quaddobl_laurent_system
    from phcpy.interface import store_quaddobl_solutions
    from phcpy.interface import load_quaddobl_solutions
    store_quaddobl_laurent_witness_set(len(embsys), dim, embsys, esols)
    py2c_copy_quaddobl_Laurent_container_to_start_system()
    py2c_copy_quaddobl_container_to_start_solutions()
    py2c_quaddobl_Laurent_cascade_homotopy()
    py2c_solve_by_quaddobl_Laurent_homotopy_continuation(tasks)
    py2c_solcon_clear_quaddobl_solutions()
    py2c_copy_quaddobl_target_solutions_to_container()
    return load_quaddobl_solutions()

def laurent_cascade_step(dim, embsys, esols, precision='d', tasks=0):
    r"""
    Given in *embsys* an embedded Laurent polynomial system and
    solutions with nonzero slack variables in *esols*,
    does one step in the homotopy cascade, with precision

    *d*: standard double precision (1.1e-15 or 2^(-53)),

    *dd*: double double precision (4.9e-32 or 2^(-104)),

    *qd*: quad double precision (1.2e-63 or 2^(-209)).

    The dimension of the solution set represented by *embsys*
    and *esols* is the value of *dim*.
    The number of tasks in multithreaded path tracking is given by *tasks*.
    The default zero value of *tasks* indicates no multithreading.
    The list on return contains witness points on
    lower dimensional solution components.
    """
    if(precision == 'd'):
        return standard_double_laurent_cascade_step(dim, embsys, esols, tasks)
    elif(precision == 'dd'):
        return double_double_laurent_cascade_step(dim, embsys, esols, tasks)
    elif(precision == 'qd'):
        return quad_double_laurent_cascade_step(dim, embsys, esols, tasks)
    else:
        print('wrong argument for precision')
        return None

def split_filter(sols, dim, tol, verbose=True):
    r"""
    Given in *sols* is a list of solutions of dimension *dim*,
    which contain a variable with name 'zz' + str(*dim*),
    which is the name of the last slack variable.
    The tolerance *tol* is used to split the list of solution in two.
    On return is a tuple of two lists of solutions (possibly empty).
    The first list of solutions has the last slack variable equal
    to zero (with respect to the tolerance *tol) and the last slack
    variable of each solution in the second list has a magnitude
    larger than *tol*.
    If *verbose*, then the length of each solution list is printed.
    """
    from phcpy.solutions import filter_zero_coordinates as filter
    lastslack = 'zz' + str(dim)
    zerosols = filter(sols, lastslack, tol, 'select')
    nonzsols = filter(sols, lastslack, tol, 'remove')
    if verbose:
        print('number of candidate generic points :', len(zerosols))
        print('number of nonsolutions :', len(nonzsols))
    return (zerosols, nonzsols)

def top_cascade(nvr, dim, pols, tol, nbtasks=0, prc='d', verbose=True):
    r"""
    Constructs an embedding of the polynomials in *pols*,
    with the number of variables in *pols* equal to *nvr*,
    where *dim* is the top dimension of the solution set.
    Applies the blackbox solver to the embedded system.
    The tolerance *tol* is used to split the solution list in the list
    of generic points and the nonsolutions for use in the cascade.
    Returns a tuple with three items:

    1. the embedded system,

    2. the solutions with zero last coordinate w.r.t. *tol*,

    3. the solutions with nonzero last coordinate w.r.t. *tol*.

    The three parameters are

    1. *nbtasks* is the number of tasks, 0 if no multitasking;

    2. the working precision *prc*, 'd' for double, 'dd' for double double,
       or 'qd' for quad double;

    3. if *verbose*, then some output is written to screen.
    """
    from phcpy.sets import embed
    from phcpy.solver import solve
    topemb = embed(nvr, dim, pols)
    if verbose:
        print('solving the embedded system at the top ...')
    topsols = solve(topemb, verbose=verbose, tasks=nbtasks, precision=prc)
    if verbose:
         print('number of solutions found :', len(topsols))
    (sols0, sols1) = split_filter(topsols, dim, tol, verbose=verbose)
    return (topemb, sols0, sols1)

def laurent_top_cascade(nvr, dim, pols, tol, \
    nbtasks=0, prc='d', verbose=True):
    r"""
    Constructs an embedding of the Laurent polynomials in *pols*,
    with the number of variables in *pols* equal to *nvr*,
    where *dim* is the top dimension of the solution set.
    Applies the blackbox solver to the embedded system.
    The tolerance *tol* is used to split the solution list in the list
    of generic points and the nonsolutions for use in the cascade.
    Returns a tuple with three items:

    1. the embedded system,

    2. the solutions with zero last coordinate w.r.t. *tol*,

    3. the solutions with nonzero last coordinate w.r.t. *tol*.

    The three parameters are

    1. *nbtasks* is the number of tasks, 0 if no multitasking;

    2. the working precision *prc*, 'd' for double, 'dd' for double double,
       or 'qd' for quad double;

    3. if *verbose*, then some output is written to screen.
    """
    from phcpy.sets import laurent_embed
    from phcpy.solver import solve
    topemb = laurent_embed(nvr, dim, pols)
    if verbose:
        print('solving the embedded system at the top ...')
    topsols = solve(topemb, verbose=verbose, tasks=nbtasks, precision=prc)
    if verbose:
         print('number of solutions found :', len(topsols))
    (sols0, sols1) = split_filter(topsols, dim, tol, verbose=verbose)
    return (topemb, sols0, sols1)

def cascade_filter(dim, embpols, nonsols, tol, \
    nbtasks=0, prc='d', verbose=True):
    r"""
    Runs one step in the cascade homotopy defined by the embedding of
    polynomials in *embpols*, starting at the solutions in *nonsols*,
    removing the last hyperplane from *embpols* at dimension *dim*.
    The tolerance *tol* is used to split filter the solutions.
    By default, the precision *prc* is double ('d').  Other valid values
    for *prc* are 'dd' (for double double) and 'qd' (for quad double).
    If *verbose*, then some output is written to screen.
    """
    if(prc=='d'):
        from phcpy.sets \
        import drop_variable_from_standard_polynomials as drop1poly
        from phcpy.sets \
        import drop_coordinate_from_standard_solutions as drop1sols
    elif(prc=='dd'):
        from phcpy.sets \
        import drop_variable_from_dobldobl_polynomials as drop1poly
        from phcpy.sets \
        import drop_coordinate_from_dobldobl_solutions as drop1sols
    elif(prc=='qd'):
        from phcpy.sets \
        import drop_variable_from_quaddobl_polynomials as drop1poly
        from phcpy.sets \
        import drop_coordinate_from_quaddobl_solutions as drop1sols
    else:
        print('invalid value for precision as argument for prc')
        return
    if verbose:
        print('running a cascade with %d paths ...' % len(nonsols))
    sols = cascade_step(dim, embpols, nonsols, precision=prc, tasks=nbtasks)
    dimslackvar = 'zz' + str(dim)
    embdown = drop1poly(embpols, dimslackvar)
    solsdrop = drop1sols(sols, len(embpols), dimslackvar)
    if dim <= 1:
        return (embdown[:-1], solsdrop)
    else:
        (sols0, sols1) = split_filter(solsdrop, dim-1, tol, verbose)
        return (embdown[:-1], sols0, sols1) 

def laurent_cascade_filter(dim, embpols, nonsols, tol, \
    nbtasks=0, prc='d', verbose=True):
    r"""
    Runs one step in the cascade homotopy defined by the embedding of
    Laurent polynomials in *embpols*, starting at the solutions in *nonsols*,
    removing the last hyperplane from *embpols* at dimension *dim*.
    The tolerance *tol* is used to split filter the solutions.
    By default, the precision *prc* is double ('d').  Other valid values
    for *prc* are 'dd' (for double double) and 'qd' (for quad double).
    If *verbose*, then some output is written to screen.
    """
    if(prc=='d'):
        from phcpy.sets \
        import drop_variable_from_standard_laurent_polynomials as drop1poly
        from phcpy.sets \
        import drop_coordinate_from_standard_solutions as drop1sols
    elif(prc=='dd'):
        from phcpy.sets \
        import drop_variable_from_dobldobl_laurent_polynomials as drop1poly
        from phcpy.sets \
        import drop_coordinate_from_dobldobl_solutions as drop1sols
    elif(prc=='qd'):
        from phcpy.sets \
        import drop_variable_from_quaddobl_laurent_polynomials as drop1poly
        from phcpy.sets \
        import drop_coordinate_from_quaddobl_solutions as drop1sols
    else:
        print('invalid value for precision as argument for prc')
        return
    if verbose:
        print('running a cascade with %d paths ...' % len(nonsols))
    sols = laurent_cascade_step\
               (dim, embpols, nonsols, precision=prc, tasks=nbtasks)
    dimslackvar = 'zz' + str(dim)
    embdown = drop1poly(embpols, dimslackvar)
    solsdrop = drop1sols(sols, len(embpols), dimslackvar)
    if dim <= 1:
        return (embdown[:-1], solsdrop)
    else:
        (sols0, sols1) = split_filter(solsdrop, dim-1, tol, verbose)
        return (embdown[:-1], sols0, sols1) 

def run_cascade(nvr, dim, pols, islaurent=False, \
    tol=1.0e-8, rcotol=1.0e-6, evatol=1.0e-6, memtol=1.0e-6, \
    tasks=0, prc='d', verbose=True):
    r"""
    Runs a cascade on the polynomials *pols*,
    in the number of variables equal to *nvr*,
    starting at the top dimension *dim*.
    If islaurent, then the polynomials in *pols* may have negative exponents.
    Returns a dictionary with as keys the dimensions
    and as values the tuples with the embedded systems
    and the corresponding generic points.
    Four tolerance parameters have default values on input:
    *tol* is used to decide which slack variables are zero,
    *rcotol* is the tolerance on the estimated inverse condition number,
    *evatol* is the tolerance on the residual to filter junk points,
    *memtol* is the tolerance for the homotopy membership test.
    The number of tasks is given by *tasks* (0 for no multitasking)
    and the default precision is double.  Other supported values
    for *prc* are 'dd' for double double and 'qd' for quad double.
    If *verbose*, then a summary of the filtering is printed.
    """
    from phcpy.sets import ismember_filter, laurent_ismember_filter
    from phcpy.cascades import top_cascade, laurent_top_cascade
    from phcpy.cascades import cascade_filter, laurent_cascade_filter
    lowdim = max(0, nvr-len(pols)) # lower bound on the dimension
    result = {}
    if islaurent:
        (topemb, topsols, nonsols) \
            = laurent_top_cascade(nvr, dim, pols, tol, tasks, prc, verbose)
    else:
        (topemb, topsols, nonsols) \
            = top_cascade(nvr, dim, pols, tol, tasks, prc, verbose)
    result[dim] = (topemb, topsols)
    for idx in range(dim, lowdim, -1):
        emb = result[idx][0]
        if(idx == 1):
            if islaurent:
                (embp, sols) = laurent_cascade_filter\
                    (idx, emb, nonsols, tol, tasks, prc, verbose)
            else:
                (embp, sols) = cascade_filter\
                    (idx, emb, nonsols, tol, tasks, prc, verbose)
        else:
            if islaurent:
                (embp, sols, nonsols) = laurent_cascade_filter\
                    (idx, emb, nonsols, tol, tasks, prc, verbose)
            else:
                (embp, sols, nonsols) = cascade_filter\
                    (idx, emb, nonsols, tol, tasks, prc, verbose)
        dims = list(result.keys())
        dims.sort(reverse=True)
        for dim in dims:
            (epols, esols) = result[dim]
            if islaurent:
                sols = laurent_ismember_filter\
                          (epols, esols, dim, sols, \
                           rcotol, evatol, memtol, False, prc)
            else:
                sols = ismember_filter(epols, esols, dim, sols, \
                                       rcotol, evatol, memtol, False, prc)
            if verbose:
                print('number of generic points after filtering :', len(sols))
        result[idx-1] = (embp, sols)
    return result

def test_run_cascade():
    """
    Runs the cascade on a list of polynomials.
    """
    testpols = ['(x1-1)*(x1-2)*(x1-3)*(x1-4);', \
                '(x1-1)*(x2-1)*(x2-2)*(x2-3);', \
                '(x1-1)*(x1-2)*(x3-1)*(x3-2);', \
                '(x1-1)*(x2-1)*(x3-1)*(x4-1);']
    ans = input('Verbose ? (y/n) ')
    vrb = (ans == 'y')
    deco = run_cascade(nvr=4, dim=3, pols=testpols, verbose=vrb)
    dims = list(deco.keys())
    dims.sort(reverse=True)
    for dim in dims:
        (epols, esols) = deco[dim]
        print('#generic points at dimension', dim, ':', len(esols))

def test_cascade():
    """
    Does one cascade step on simple example.
    In the top embedding we first find the 2-dimensional
    solution set x = 1.  In the cascade step we compute
    the three witness points on the twisted cubic.
    """
    from phcpy.solver import solve
    from phcpy.sets import embed
    from phcpy.sets import drop_variable_from_standard_polynomials
    from phcpy.sets import drop_coordinate_from_standard_solutions
    pols = ['(x - 1)*(y-x^2);', \
            '(x - 1)*(z-x^3);', \
            '(x^2 - 1)*(y-x^2);' ]
    print(pols)
    (embpols, sols0, sols1) = top_cascade(3, 2, pols, 1.0e-8)
    print('the embedded system :')
    print(embpols)
    print('the generic points :')
    for sol in sols0:
        print(sol)
    input('hit enter to continue...')
    print('solutions with nonzero slack variables :')
    for sol in sols1:
        print(sol)
    input('hit enter to continue...')
    print('... running cascade step ...')
    (wp1, ws0, ws1) = cascade_filter(2, embpols, sols1, 1.0e-8)
    print('the 1-dimensional embedding :')
    for pol in wp1:
        print(pol)
    print('the candidate witness points :')
    for sol in ws0:
        print(sol)

def test():
    """
    Fixes a seed for the random number generators
    before running the test on the cascade homotopies.
    """
    from phcpy.phcpy2c3 import py2c_set_seed
    py2c_set_seed(234798272)
    # test_cascade()
    test_run_cascade()

if __name__ == "__main__":
    test()
