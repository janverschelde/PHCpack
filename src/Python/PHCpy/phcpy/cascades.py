r"""
A cascade homotopy removes one hyperplane from an embedded system,
taking the solutions with nonzero slack variables to solutions on
lower dimensional components of the solution set of the original system.
"""
from ctypes import c_int32, c_double, pointer
from phcpy.version import get_phcfun
from phcpy.solutions import filter_zero_coordinates
from phcpy.homotopies import copy_double_system_into_start
from phcpy.homotopies import copy_double_solutions_into_start
from phcpy.homotopies import copy_double_double_system_into_start
from phcpy.homotopies import copy_double_double_solutions_into_start
from phcpy.homotopies import copy_quad_double_system_into_start
from phcpy.homotopies import copy_quad_double_solutions_into_start
from phcpy.homotopies import get_double_target_solutions
from phcpy.homotopies import get_double_double_target_solutions
from phcpy.homotopies import get_quad_double_target_solutions
from phcpy.trackers import do_double_track
from phcpy.trackers import do_double_double_track
from phcpy.trackers import do_quad_double_track
from phcpy.solver import solve
from phcpy.sets import set_double_witness_set
from phcpy.sets import set_double_double_witness_set
from phcpy.sets import set_quad_double_witness_set
from phcpy.sets import double_embed, double_double_embed, quad_double_embed
from phcpy.sets import drop_variable_from_double_polynomials
from phcpy.sets import drop_variable_from_double_double_polynomials
from phcpy.sets import drop_variable_from_quad_double_polynomials
from phcpy.sets import drop_coordinate_from_double_solutions
from phcpy.sets import drop_coordinate_from_double_double_solutions
from phcpy.sets import drop_coordinate_from_quad_double_solutions

def set_double_cascade_homotopy(vrblvl=0):
    """
    Defines the cascade homotopy in double precision,
    called as a function in the double cascade step.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_double_cascade_homotopy ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_cascade_homotopy calls phc', end='')
    retval = phc(164, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_double_cascade_homotopy(vrblvl=0):
    """
    Defines the cascade homotopy in double double precision,
    called as a function in the double double cascade step.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_double_double_cascade_homotopy ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_double_cascade_homotopy calls phc', end='')
    retval = phc(178, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_quad_double_cascade_homotopy(vrblvl=0):
    """
    Defines the cascade homotopy in quad double precision,
    called as a function in the quad double cascade step.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_quad_double_cascade_homotopy ...')
    phc = get_phcfun()
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_quad_double_cascade_homotopy calls phc', end='')
    retval = phc(188, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_cascade_step(dim, embsys, esols, tasks=0, vrblvl=0):
    r"""
    Given in *embsys* an embedded polynomial system and
    solutions with nonzero slack variables in *esols*, does one step 
    in the homotopy cascade, with double precision arithmetic.
    The dimension of the solution set represented by *embsys*
    and *esols* is the value of *dim*.
    The number of tasks in multithreaded path tracking is given by *tasks*.
    The default zero value of *tasks* indicates no multithreading.
    The list on return contains witness points on
    lower dimensional solution components.
    """
    if vrblvl > 0:
        print('in double_cascade_step, dim :', dim, end='')
        print(', tasks :', tasks)
        print('the polynomials :')
        for pol in embsys:
            print(pol)
        print('the solutions :')
        for (idx, sol) in enumerate(esols):
            print('Solution', idx+1, ':')
            print(sol)
    set_double_witness_set(len(embsys), dim, embsys, esols, vrblvl)
    copy_double_system_into_start(vrblvl)
    copy_double_solutions_into_start(vrblvl)
    set_double_cascade_homotopy(vrblvl)
    do_double_track(tasks, vrblvl)
    sols = get_double_target_solutions(vrblvl)
    return sols

def double_double_cascade_step(dim, embsys, esols, tasks=0, vrblvl=0):
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
    if vrblvl > 0:
        print('in double_double_cascade_step, dim :', dim, end='')
        print(', tasks :', tasks)
        print('the polynomials :')
        for pol in embsys:
            print(pol)
        print('the solutions :')
        for (idx, sol) in enumerate(esols):
            print('Solution', idx+1, ':')
            print(sol)
    set_double_double_witness_set(len(embsys), dim, embsys, esols, vrblvl)
    copy_double_double_system_into_start(vrblvl)
    copy_double_double_solutions_into_start(vrblvl)
    set_double_double_cascade_homotopy(vrblvl)
    do_double_double_track(tasks, vrblvl)
    sols = get_double_double_target_solutions(vrblvl)
    return sols

def quad_double_cascade_step(dim, embsys, esols, tasks=0, vrblvl=0):
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
    if vrblvl > 0:
        print('in quad_double_cascade_step, dim :', dim, end='')
        print(', tasks :', tasks)
        print('the polynomials :')
        for pol in embsys:
            print(pol)
        print('the solutions :')
        for (idx, sol) in enumerate(esols):
            print('Solution', idx+1, ':')
            print(sol)
    set_quad_double_witness_set(len(embsys), dim, embsys, esols, vrblvl)
    copy_quad_double_system_into_start(vrblvl)
    copy_quad_double_solutions_into_start(vrblvl)
    set_quad_double_cascade_homotopy(vrblvl)
    do_quad_double_track(tasks, vrblvl)
    sols = get_quad_double_target_solutions(vrblvl)
    return sols

def split_filter(sols, dim, tol, vrblvl=0):
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
    If *vrblvl* is nonzero, then the length of each solution list is printed.
    """
    if vrblvl > 0:
        print('in split_filter, dim :', dim, end='')
        print(', tol :', tol)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    lastslack = 'zz' + str(dim)
    zerosols = filter_zero_coordinates(sols, lastslack, tol, 'select')
    nonzsols = filter_zero_coordinates(sols, lastslack, tol, 'remove')
    if vrblvl > 0:
        print('number of candidate generic points :', len(zerosols))
        print('number of nonsolutions :', len(nonzsols))
    return (zerosols, nonzsols)

def double_top_cascade(nvr, dim, pols, tol=1.0e-6, nbtasks=0, vrblvl=0):
    r"""
    Constructs an embedding of the polynomials in *pols*,
    with the number of variables in *pols* equal to *nvr*,
    where *dim* is the top dimension of the solution set.
    Applies the blackbox solver to the embedded system to compute
    in double precision the generic points and the nonsolutions 
    for use in the cascade.  Returns a tuple with three items:

    1. the embedded system,

    2. the solutions with zero last coordinate w.r.t. *tol*,

    3. the solutions with nonzero last coordinate w.r.t. *tol*.

    The three parameters on input are

    1. *tol* is the tolerance to decide whether a number is zero or not,
       used to split the solution list of the embedded system;

    2. *nbtasks* is the number of tasks, 0 if no multitasking,

    3. if *verbose*, then some output is written to screen.
    """
    topemb = double_embed(nvr, dim, pols, vrblvl)
    if vrblvl > 0:
        print('solving the embedded system at the top ...')
    topsols = solve(topemb, tasks=nbtasks, verbose_level=vrblvl)
    if vrblvl > 0:
        print('number of solutions found :', len(topsols))
    (sols0, sols1) = split_filter(topsols, dim, tol, vrblvl)
    return (topemb, sols0, sols1)

def double_double_top_cascade(nvr, dim, pols, \
    tol=1.0e-6, nbtasks=0, vrblvl=0):
    r"""
    Constructs an embedding of the polynomials in *pols*,
    with the number of variables in *pols* equal to *nvr*,
    where *dim* is the top dimension of the solution set.
    Applies the blackbox solver to the embedded system to compute
    in double double precision the generic points and the nonsolutions 
    for use in the cascade.  Returns a tuple with three items:

    1. the embedded system,

    2. the solutions with zero last coordinate w.r.t. *tol*,

    3. the solutions with nonzero last coordinate w.r.t. *tol*.

    The three parameters on input are

    1. *tol* is the tolerance to decide whether a number is zero or not,
       used to split the solution list of the embedded system;

    2. *nbtasks* is the number of tasks, 0 if no multitasking,

    3. if *verbose*, then some output is written to screen.
    """
    topemb = double_double_embed(nvr, dim, pols, vrblvl)
    if vrblvl > 0:
        print('solving the embedded system at the top ...')
    topsols = solve(topemb, tasks=nbtasks, precision='dd', \
        verbose_level=vrblvl)
    if vrblvl > 0:
        print('number of solutions found :', len(topsols))
    (sols0, sols1) = split_filter(topsols, dim, tol, vrblvl)
    return (topemb, sols0, sols1)

def quad_double_top_cascade(nvr, dim, pols, \
    tol=1.0e-6, nbtasks=0, vrblvl=0):
    r"""
    Constructs an embedding of the polynomials in *pols*,
    with the number of variables in *pols* equal to *nvr*,
    where *dim* is the top dimension of the solution set.
    Applies the blackbox solver to the embedded system to compute
    in quad double precision the generic points and the nonsolutions 
    for use in the cascade.  Returns a tuple with three items:

    1. the embedded system,

    2. the solutions with zero last coordinate w.r.t. *tol*,

    3. the solutions with nonzero last coordinate w.r.t. *tol*.

    The three parameters on input are

    1. *tol* is the tolerance to decide whether a number is zero or not,
       used to split the solution list of the embedded system;

    2. *nbtasks* is the number of tasks, 0 if no multitasking,

    3. if *verbose*, then some output is written to screen.
    """
    topemb = quad_double_embed(nvr, dim, pols, vrblvl)
    if vrblvl > 0:
        print('solving the embedded system at the top ...')
    topsols = solve(topemb, tasks=nbtasks, precision='qd', \
        verbose_level=vrblvl)
    if vrblvl > 0:
        print('number of solutions found :', len(topsols))
    (sols0, sols1) = split_filter(topsols, dim, tol, vrblvl)
    return (topemb, sols0, sols1)

def double_cascade_filter(dim, embpols, nonsols, tol, \
    nbtasks=0, vrblvl=0):
    r"""
    Runs one step in the cascade homotopy defined by the embedding of
    polynomials in *embpols*, starting at the solutions in *nonsols*,
    removing the last hyperplane from *embpols* at dimension *dim*.
    The tolerance *tol* is used to split filter the solutions.
    Computations happen in double precision.
    If *verbose*, then some output is written to screen.
    """
    if vrblvl > 0:
        print('in double_cascade_filter, dim :', dim)
        print('the polynomials :')
        for pol in embpols:
            print(pol)
        print('the nonsolutions :')
        for (idx, sol) in enumerate(nonsols):
            print('Solution', idx+1, ':')
            print(sol)
    if vrblvl > 0:
        print('running a cascade with %d paths ...' % len(nonsols))
    sols = double_cascade_step(dim, embpols, nonsols, \
        tasks=nbtasks, vrblvl=vrblvl)
    dimslackvar = 'zz' + str(dim)
    embdown = drop_variable_from_double_polynomials(embpols, \
        dimslackvar, vrblvl)
    solsdrop = drop_coordinate_from_double_solutions(sols, len(embpols), \
        dimslackvar, vrblvl)
    if dim <= 1:
        return (embdown[:-1], solsdrop, None)
    (sols0, sols1) = split_filter(solsdrop, dim-1, tol, vrblvl)
    return (embdown[:-1], sols0, sols1) 

def double_double_cascade_filter(dim, embpols, nonsols, tol, \
    nbtasks=0, vrblvl=0):
    r"""
    Runs one step in the cascade homotopy defined by the embedding of
    polynomials in *embpols*, starting at the solutions in *nonsols*,
    removing the last hyperplane from *embpols* at dimension *dim*.
    Computations happen in double double precision.
    The tolerance *tol* is used to split filter the solutions.
    If *verbose*, then some output is written to screen.
    """
    if vrblvl > 0:
        print('in double_double_cascade_filter, dim :', dim)
        print('the polynomials :')
        for pol in embpols:
            print(pol)
        print('the nonsolutions :')
        for (idx, sol) in enumerate(nonsols):
            print('Solution', idx+1, ':')
            print(sol)
    if vrblvl > 0:
        print('running a cascade with %d paths ...' % len(nonsols))
    sols = double_double_cascade_step(dim, embpols, nonsols, \
        tasks=nbtasks, vrblvl=vrblvl)
    dimslackvar = 'zz' + str(dim)
    embdown = drop_variable_from_double_double_polynomials(embpols, \
        dimslackvar, vrblvl)
    solsdrop = drop_coordinate_from_double_double_solutions(sols, \
        len(embpols), dimslackvar, vrblvl)
    if dim <= 1:
        return (embdown[:-1], solsdrop, None)
    (sols0, sols1) = split_filter(solsdrop, dim-1, tol, vrblvl)
    return (embdown[:-1], sols0, sols1) 

def quad_double_cascade_filter(dim, embpols, nonsols, tol, \
    nbtasks=0, vrblvl=0):
    r"""
    Runs one step in the cascade homotopy defined by the embedding of
    polynomials in *embpols*, starting at the solutions in *nonsols*,
    removing the last hyperplane from *embpols* at dimension *dim*.
    Computations happen in quad double precision.
    The tolerance *tol* is used to split filter the solutions.
    If *verbose*, then some output is written to screen.
    """
    if vrblvl > 0:
        print('in quad_double_cascade_filter, dim :', dim)
        print('the polynomials :')
        for pol in embpols:
            print(pol)
        print('the nonsolutions :')
        for (idx, sol) in enumerate(nonsols):
            print('Solution', idx+1, ':')
            print(sol)
    if vrblvl > 0:
        print('running a cascade with %d paths ...' % len(nonsols))
    sols = quad_double_cascade_step(dim, embpols, nonsols, \
        tasks=nbtasks, vrblvl=vrblvl)
    dimslackvar = 'zz' + str(dim)
    embdown = drop_variable_from_quad_double_polynomials(embpols, \
        dimslackvar, vrblvl)
    solsdrop = drop_coordinate_from_quad_double_solutions(sols, \
        len(embpols), dimslackvar, vrblvl)
    if dim <= 1:
        return (embdown[:-1], solsdrop, None)
    (sols0, sols1) = split_filter(solsdrop, dim-1, tol, vrblvl)
    return (embdown[:-1], sols0, sols1) 

def test_double_cascade(vrblvl=0):
    """
    Tests one cascade step in double precision.
    In the top embedding we first find the 2-dimensional
    solution set x = 1.  In the cascade step we compute
    the candidate witness points on the twisted cubic.
    """
    if vrblvl > 0:
        print('in test_double_cascade ...')
    pols = ['(x - 1)*(y-x^2);', \
            '(x - 1)*(z-x^3);', \
            '(x^2 - 1)*(y-x^2);' ]
    if vrblvl > 0:
        for pol in pols:
            print(pol)
    (embpols, sols0, sols1) = double_top_cascade(3, 2, pols, 1.0e-8, \
        vrblvl=vrblvl)
    if vrblvl > 0:
        print('the embedded system :')
        for pol in embpols:
            print(pol)
        print('the generic points :')
        for (idx, sol) in enumerate(sols0):
            print('Solution', idx+1, ':')
            print(sol)
    fail = int(len(sols0) != 1)
    if vrblvl > 0:
        if fail != 0:
            print('Failure: expected one generic point!')
        else:
            print('As expected, got one generic point.')
    # input('hit enter to continue...')
    if vrblvl > 0:
        print('solutions with nonzero slack variables :')
        for (idx, sol) in enumerate(sols1):
            print('Solution', idx+1, ':')
            print(sol)
    # input('hit enter to continue...')
    if vrblvl > 0:
        print('... running cascade step ...')
    (wp1, ws0, ws1) = double_cascade_filter(2, embpols, sols1, 1.0e-8, \
        vrblvl=vrblvl)
    if vrblvl > 0:
        print('the 1-dimensional embedding :')
        for pol in wp1:
            print(pol)
        print('the candidate witness points :')
        for (idx, sol) in enumerate(ws0):
            print('Solution', idx+1, ':')
            print(sol)
    if len(ws0) != 4:
        if vrblvl > 0:
            print('Failure: expected four candidate generic points.')
        fail = fail + 1
    else:
        if vrblvl > 0:
            print('As expected, got four candidate generic points.')
    return fail

def test_double_double_cascade(vrblvl=0):
    """
    Tests one cascade step in double double precision.
    In the top embedding we first find the 2-dimensional
    solution set x = 1.  In the cascade step we compute
    the candidate witness points on the twisted cubic.
    """
    if vrblvl > 0:
        print('in test_double_double_cascade ...')
    pols = ['(x - 1)*(y-x^2);', \
            '(x - 1)*(z-x^3);', \
            '(x^2 - 1)*(y-x^2);' ]
    if vrblvl > 0:
        for pol in pols:
            print(pol)
    (embpols, sols0, sols1) = double_double_top_cascade(3, 2, pols, \
        1.0e-8, vrblvl=vrblvl)
    if vrblvl > 0:
        print('the embedded system :')
        for pol in embpols:
            print(pol)
        print('the generic points :')
        for (idx, sol) in enumerate(sols0):
            print('Solution', idx+1, ':')
            print(sol)
    fail = int(len(sols0) != 1)
    if vrblvl > 0:
        if fail != 0:
            print('Failure: expected one generic point!')
        else:
            print('As expected, got one generic point.')
    # input('hit enter to continue...')
    if vrblvl > 0:
        print('solutions with nonzero slack variables :')
        for (idx, sol) in enumerate(sols1):
            print('Solution', idx+1, ':')
            print(sol)
    # input('hit enter to continue...')
    if vrblvl > 0:
        print('... running cascade step ...')
    (wp1, ws0, ws1) = double_double_cascade_filter(2, embpols, sols1, \
        1.0e-8, nbtasks=1, vrblvl=vrblvl)
    if vrblvl > 0:
        print('the 1-dimensional embedding :')
        for pol in wp1:
            print(pol)
        print('the candidate witness points :')
        for (idx, sol) in enumerate(ws0):
            print('Solution', idx+1, ':')
            print(sol)
    if len(ws0) != 4:
        if vrblvl > 0:
            print('Failure: expected four candidate generic points.')
        fail = fail + 1
    else:
        print('As expected, got four candidate generic points.')
    return fail

def test_quad_double_cascade(vrblvl=0):
    """
    Tests one cascade step in quad double precision.
    In the top embedding we first find the 2-dimensional
    solution set x = 1.  In the cascade step we compute
    the candidate witness points on the twisted cubic.
    """
    if vrblvl > 0:
        print('in test_quad_double_cascade ...')
    pols = ['(x - 1)*(y-x^2);', \
            '(x - 1)*(z-x^3);', \
            '(x^2 - 1)*(y-x^2);' ]
    if vrblvl > 0:
        for pol in pols:
            print(pol)
    (embpols, sols0, sols1) = quad_double_top_cascade(3, 2, pols, \
        1.0e-8, vrblvl=vrblvl)
    if vrblvl > 0:
        print('the embedded system :')
        for pol in embpols:
            print(pol)
        print('the generic points :')
        for (idx, sol) in enumerate(sols0):
            print('Solution', idx+1, ':')
            print(sol)
    fail = int(len(sols0) != 1)
    if vrblvl > 0:
        if fail != 0:
            print('Failure: expected one generic point!')
        else:
            print('As expected, got one generic point.')
    # input('hit enter to continue...')
    if vrblvl > 0:
        print('solutions with nonzero slack variables :')
        for (idx, sol) in enumerate(sols1):
            print('Solution', idx+1, ':')
            print(sol)
    # input('hit enter to continue...')
    if vrblvl > 0:
        print('... running cascade step ...')
    (wp1, ws0, ws1) = quad_double_cascade_filter(2, embpols, sols1, \
        1.0e-8, nbtasks=1, vrblvl=vrblvl)
    if vrblvl > 0:
        print('the 1-dimensional embedding :')
        for pol in wp1:
            print(pol)
        print('the candidate witness points :')
        for (idx, sol) in enumerate(ws0):
            print('Solution', idx+1, ':')
            print(sol)
    if len(ws0) != 4:
        if vrblvl > 0:
            print('Failure: expected four candidate generic points.')
        fail = fail + 1
    else:
        print('As expected, got four candidate generic points.')
    return fail

def main():
    """
    Runs some tests.
    """
    lvl = 10
    fail = test_double_cascade(lvl)
    fail = fail + test_double_double_cascade(lvl)
    fail = fail + test_quad_double_cascade(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()
