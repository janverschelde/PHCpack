"""
Given two witness sets for two pure dimensional solution sets,
a diagonal homotopy computes a sets of witness sets for all components
of the intersection of the two pure dimensional solution sets.
"""
from ctypes import c_int32, c_double, pointer
from phcpy.version import get_phcfun, str2int4a
from phcpy.polynomials import number_of_symbols, string_of_symbols
from phcpy.polynomials import set_double_system
from phcpy.polynomials import get_double_system
from phcpy.polynomials import clear_double_system
from phcpy.polynomials import set_double_double_system
from phcpy.polynomials import get_double_double_system
from phcpy.polynomials import clear_double_double_system
from phcpy.polynomials import set_quad_double_system
from phcpy.polynomials import get_quad_double_system
from phcpy.polynomials import clear_quad_double_system
from phcpy.solutions import set_double_solutions
from phcpy.solutions import get_double_solutions
from phcpy.solutions import clear_double_solutions
from phcpy.solutions import set_double_double_solutions
from phcpy.solutions import get_double_double_solutions
from phcpy.solutions import clear_double_double_solutions
from phcpy.solutions import set_quad_double_solutions
from phcpy.solutions import get_quad_double_solutions
from phcpy.solutions import clear_quad_double_solutions
from phcpy.solutions import filter_vanishing
from phcpy.homotopies import copy_double_system_into_start
from phcpy.homotopies import copy_double_solutions_into_start
from phcpy.homotopies import copy_double_system_into_target
from phcpy.homotopies import copy_double_system_from_target
from phcpy.homotopies import copy_double_solutions_into_target
from phcpy.homotopies import copy_double_double_system_into_start
from phcpy.homotopies import copy_double_double_solutions_into_start
from phcpy.homotopies import copy_double_double_system_into_target
from phcpy.homotopies import copy_double_double_system_from_target
from phcpy.homotopies import copy_double_double_solutions_into_target
from phcpy.homotopies import copy_quad_double_system_into_start
from phcpy.homotopies import copy_quad_double_solutions_into_start
from phcpy.homotopies import copy_quad_double_system_into_target
from phcpy.homotopies import copy_quad_double_system_from_target
from phcpy.homotopies import copy_quad_double_solutions_into_target
from phcpy.homotopies import set_double_homotopy
from phcpy.homotopies import set_double_double_homotopy
from phcpy.homotopies import set_quad_double_homotopy
from phcpy.homotopies import get_double_target_solutions
from phcpy.homotopies import get_double_double_target_solutions
from phcpy.homotopies import get_quad_double_target_solutions
from phcpy.trackers import do_double_track
from phcpy.trackers import do_double_double_track
from phcpy.trackers import do_quad_double_track
from phcpy.sets import double_hypersurface_set
from phcpy.sets import double_double_hypersurface_set
from phcpy.sets import quad_double_hypersurface_set
from phcpy.sets import drop_coordinate_from_double_solutions
from phcpy.sets import drop_variable_from_double_polynomials
from phcpy.sets import drop_coordinate_from_double_double_solutions
from phcpy.sets import drop_variable_from_double_double_polynomials
from phcpy.sets import drop_coordinate_from_quad_double_solutions
from phcpy.sets import drop_variable_from_quad_double_polynomials
from phcpy.cascades import double_cascade_step
from phcpy.cascades import double_double_cascade_step
from phcpy.cascades import quad_double_cascade_step

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

def set_double_diagonal_homotopy(dim1, dim2, vrblvl=0):
    """
    Defines a diagonal homotopy to intersect two solution sets of
    dimensions dim1 and dim2 respectively, where dim1 >= dim2.
    The systems that define the witness sets must have been set
    already in double precision before the call to this function.
    """
    if vrblvl > 0:
        print('in set_double_diagonal_homotopy, dim1 :', dim1, end='')
        print(', dim2 :', dim2)
    phc = get_phcfun()
    aaa = pointer(c_int32(dim1))
    bbb = pointer(c_int32(dim2))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_diagonal_homotopy calls phc', end='')
    retval = phc(165, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_double_diagonal_homotopy(dim1, dim2, vrblvl=0):
    """
    Defines a diagonal homotopy to intersect two solution sets of
    dimensions dim1 and dim2 respectively, where dim1 >= dim2.
    The systems that define the witness sets must have been set
    already in double double precision before the call to this function.
    """
    if vrblvl > 0:
        print('in set_double_double_diagonal_homotopy, dim1 :', \
            dim1, end='')
        print(', dim2 :', dim2)
    phc = get_phcfun()
    aaa = pointer(c_int32(dim1))
    bbb = pointer(c_int32(dim2))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_double_diagonal_homotopy calls phc', end='')
    retval = phc(289, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_quad_double_diagonal_homotopy(dim1, dim2, vrblvl=0):
    """
    Defines a diagonal homotopy to intersect two solution sets of
    dimensions dim1 and dim2 respectively, where dim1 >= dim2.
    The systems that define the witness sets must have been set
    already in quad double precision before the call to this function.
    """
    if vrblvl > 0:
        print('in set_quad_double_diagonal_homotopy, dim1 :', \
            dim1, end='')
        print(', dim2 :', dim2)
    phc = get_phcfun()
    aaa = pointer(c_int32(dim1))
    bbb = pointer(c_int32(dim2))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_quad_double_diagonal_homotopy calls phc', end='')
    retval = phc(290, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def diagonal_symbols_doubler(nbr, dim, symbols, vrblvl=0):
    """
    Doubles the number of symbols in the symbol table to enable the
    writing of the target system to string properly when starting the
    cascade of a diagonal homotopy in extrinsic coordinates.
    On input are nbr, the ambient dimension = #variables before the embedding,
    dim is the number of slack variables, or the dimension of the first set,
    and in symbols are the symbols for the first witness set.
    This function takes the symbols in s and combines those symbols with
    those in the current symbol table for the second witness set stored
    in the standard systems container.  On return, the symbol table
    contains then all symbols to write the top system in the cascade
    to start the diagonal homotopy.
    """
    if vrblvl > 0:
        print('in diagonal_symbols_doubler, nbr :', nbr, end='')
        print(', dim :', dim)
        print('the symbols :', symbols)
    symseq = ' '.join(symbols)
    if vrblvl > 0:
        print('the sequence of symbols :', symseq)
    phc = get_phcfun()
    pars = (c_int32 * 3)()
    pars[0] = nbr
    pars[1] = dim
    pars[2] = len(symseq)
    apars = pointer(pars)
    bsymseq = str2int4a(symseq)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> diagonal_symbols_doubler calls phc', end='')
    retval = phc(230, apars, bsymseq, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_diagonal_homotopy(dim1, sys1, esols1, \
    dim2, sys2, esols2, vrblvl=0):
    r"""
    Defines a diagonal homotopy to intersect the witness sets defined
    by (*sys1*, *esols1*) and (*sys2*, *esols2*), respectively of dimensions
    *dim1* and *dim2*.  The systems *sys1* and *sys2* are assumed to be square
    and with as many slack variables as the dimension of the solution sets.
    The data is stored in double precision.
    """
    if vrblvl > 0:
        print('in double_diagonal_homotopy, dim1 :', dim1, end='')
        print(', dim2 :', dim2)
        print('polynomials in the first set :')
        for pol in sys1:
            print(pol)
        print('generic points in the first set :')
        for (idx, sol) in enumerate(esols1):
            print('Solution', idx+1, ':')
            print(sol)
    nbsymbs = number_of_symbols(sys1, vrblvl)
    symbols = string_of_symbols(maxlen=1024, vrblvl=vrblvl)
    if vrblvl > 0:
        print('number of symbols :', nbsymbs)
        print('names of variables :', symbols)
    set_double_system(len(sys1), sys1, vrblvl)
    set_double_solutions(len(sys1), esols1, vrblvl)
    if(dim1 >= dim2):
        copy_double_system_into_target(vrblvl)
        copy_double_solutions_into_target(vrblvl)
    else:
        copy_double_system_into_start(vrblvl)
        copy_double_solutions_into_start(vrblvl)
    set_double_system(len(sys2), sys2, vrblvl)
    set_double_solutions(len(sys2), esols2, vrblvl)
    if(dim1 >= dim2):
        copy_double_system_into_start(vrblvl)
        copy_double_solutions_into_start(vrblvl)
    else:
        copy_double_system_into_target(vrblvl)
        copy_double_solutions_into_target(vrblvl)
    if(dim1 >= dim2):
        set_double_diagonal_homotopy(dim1, dim2, vrblvl)
    else:
        set_double_diagonal_homotopy(dim2, dim1, vrblvl)
    return diagonal_symbols_doubler(nbsymbs-dim1, dim1, symbols, vrblvl)

def double_double_diagonal_homotopy(dim1, sys1, esols1, \
    dim2, sys2, esols2, vrblvl=0):
    r"""
    Defines a diagonal homotopy to intersect the witness sets defined
    by (*sys1*, *esols1*) and (*sys2*, *esols2*), respectively of dimensions
    *dim1* and *dim2*.  The systems *sys1* and *sys2* are assumed to be square
    and with as many slack variables as the dimension of the solution sets.
    The data is stored in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_diagonal_homotopy, dim1 :', dim1, end='')
        print(', dim2 :', dim2)
        print('polynomials in the first set :')
        for pol in sys1:
            print(pol)
        print('generic points in the first set :')
        for (idx, sol) in enumerate(esols1):
            print('Solution', idx+1, ':')
            print(sol)
    nbsymbs = number_of_symbols(sys1, vrblvl)
    symbols = string_of_symbols(maxlen=1024, vrblvl=vrblvl)
    if vrblvl > 0:
        print('number of symbols :', nbsymbs)
        print('names of variables :', symbols)
    set_double_double_system(len(sys1), sys1, vrblvl)
    set_double_double_solutions(len(sys1), esols1, vrblvl)
    if(dim1 >= dim2):
        copy_double_double_system_into_target(vrblvl)
        copy_double_double_solutions_into_target(vrblvl)
    else:
        copy_double_double_system_into_start(vrblvl)
        copy_double_double_solutions_into_start(vrblvl)
    set_double_double_system(len(sys2), sys2, vrblvl)
    set_double_double_solutions(len(sys2), esols2, vrblvl)
    if(dim1 >= dim2):
        copy_double_double_system_into_start(vrblvl)
        copy_double_double_solutions_into_start(vrblvl)
    else:
        copy_double_double_system_into_target(vrblvl)
        copy_double_double_solutions_into_target(vrblvl)
    if(dim1 >= dim2):
        set_double_double_diagonal_homotopy(dim1, dim2, vrblvl)
    else:
        set_double_double_diagonal_homotopy(dim2, dim1, vrblvl)
    return diagonal_symbols_doubler(nbsymbs-dim1, dim1, symbols, vrblvl)

def quad_double_diagonal_homotopy(dim1, sys1, esols1, \
    dim2, sys2, esols2, vrblvl=0):
    r"""
    Defines a diagonal homotopy to intersect the witness sets defined
    by (*sys1*, *esols1*) and (*sys2*, *esols2*), respectively of dimensions
    *dim1* and *dim2*.  The systems *sys1* and *sys2* are assumed to be square
    and with as many slack variables as the dimension of the solution sets.
    The data is stored in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_diagonal_homotopy, dim1 :', dim1, end='')
        print(', dim2 :', dim2)
        print('polynomials in the first set :')
        for pol in sys1:
            print(pol)
        print('generic points in the first set :')
        for (idx, sol) in enumerate(esols1):
            print('Solution', idx+1, ':')
            print(sol)
    nbsymbs = number_of_symbols(sys1, vrblvl)
    symbols = string_of_symbols(maxlen=1024, vrblvl=vrblvl)
    if vrblvl > 0:
        print('number of symbols :', nbsymbs)
        print('names of variables :', symbols)
    set_quad_double_system(len(sys1), sys1, vrblvl)
    set_quad_double_solutions(len(sys1), esols1, vrblvl)
    if(dim1 >= dim2):
        copy_quad_double_system_into_target(vrblvl)
        copy_quad_double_solutions_into_target(vrblvl)
    else:
        copy_quad_double_system_into_start(vrblvl)
        copy_quad_double_solutions_into_start(vrblvl)
    set_quad_double_system(len(sys2), sys2, vrblvl)
    set_quad_double_solutions(len(sys2), esols2, vrblvl)
    if(dim1 >= dim2):
        copy_quad_double_system_into_start(vrblvl)
        copy_quad_double_solutions_into_start(vrblvl)
    else:
        copy_quad_double_system_into_target(vrblvl)
        copy_quad_double_solutions_into_target(vrblvl)
    if(dim1 >= dim2):
        set_quad_double_diagonal_homotopy(dim1, dim2, vrblvl)
    else:
        set_quad_double_diagonal_homotopy(dim2, dim1, vrblvl)
    return diagonal_symbols_doubler(nbsymbs-dim1, dim1, symbols, vrblvl)

def double_diagonal_cascade_solutions(dim1, dim2, vrblvl=0):
    """
    Defines the start solutions in the cascade to start the diagonal
    homotopy to intersect a set of dimension *dim1* with another set
    of dimension *dim2*, in double precision.  For this to work,
    double_diagonal_homotopy must have been executed successfully.
    """
    if vrblvl > 0:
        print('in double_diagonal_cascade_solutions, dim1 :', dim1, end='')
        print(', dim2 :', dim2)
    phc = get_phcfun()
    if dim1 >= dim2: 
        adim1 = pointer(c_int32(dim1))
        bdim2 = pointer(c_int32(dim2))
    else:
        adim1 = pointer(c_int32(dim2))
        bdim2 = pointer(c_int32(dim1))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_diagonal_cascade_solutions calls phc', end='')
    retval = phc(271, adim1, bdim2, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_double_diagonal_cascade_solutions(dim1, dim2, vrblvl=0):
    """
    Defines the start solutions in the cascade to start the diagonal
    homotopy to intersect a set of dimension *dim1* with another set
    of dimension *dim2*, in double double precision.  For this to work,
    double_double_diagonal_homotopy must have been executed successfully.
    """
    if vrblvl > 0:
        print('in double_double_diagonal_cascade_solutions, dim1 :', \
            dim1, end='')
        print(', dim2 :', dim2)
    phc = get_phcfun()
    if dim1 >= dim2: 
        adim1 = pointer(c_int32(dim1))
        bdim2 = pointer(c_int32(dim2))
    else:
        adim1 = pointer(c_int32(dim2))
        bdim2 = pointer(c_int32(dim1))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_double_diagonal_cascade_solutions calls phc', end='')
    retval = phc(297, adim1, bdim2, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def quad_double_diagonal_cascade_solutions(dim1, dim2, vrblvl=0):
    """
    Defines the start solutions in the cascade to start the diagonal
    homotopy to intersect a set of dimension *dim1* with another set
    of dimension *dim2*, in quad double precision.  For this to work,
    quad_double_diagonal_homotopy must have been executed successfully.
    """
    if vrblvl > 0:
        print('in quad_double_diagonal_cascade_solutions, dim1 :', \
            dim1, end='')
        print(', dim2 :', dim2)
    phc = get_phcfun()
    if dim1 >= dim2: 
        adim1 = pointer(c_int32(dim1))
        bdim2 = pointer(c_int32(dim2))
    else:
        adim1 = pointer(c_int32(dim2))
        bdim2 = pointer(c_int32(dim1))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> quad_double_diagonal_cascade_solutions calls phc', end='')
    retval = phc(298, adim1, bdim2, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_double_diagonal_cascade_solutions(dim1, dim2, vrblvl=0):
    """
    Defines the start solutions in the cascade to start the diagonal
    homotopy to intersect a set of dimension *dim1* with another set
    of dimension *dim2*, in double double precision.  For this to work,
    double_double_diagonal_homotopy must have been executed successfully.
    """
    if vrblvl > 0:
        print('in double_double_diagonal_cascade_solutions, dim1 :', \
            dim1, end='')
        print(', dim2 :', dim2)
    phc = get_phcfun()
    if dim1 >= dim2: 
        adim1 = pointer(c_int32(dim1))
        bdim2 = pointer(c_int32(dim2))
    else:
        adim1 = pointer(c_int32(dim2))
        bdim2 = pointer(c_int32(dim1))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_double_diagonal_cascade_solutions calls phc', end='')
    retval = phc(297, adim1, bdim2, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def quad_double_diagonal_cascade_solutions(dim1, dim2, vrblvl=0):
    """
    Defines the start solutions in the cascade to start the diagonal
    homotopy to intersect a set of dimension *dim1* with another set
    of dimension *dim2*, in quad double precision.  For this to work,
    quad_double_diagonal_homotopy must have been executed successfully.
    """
    if vrblvl > 0:
        print('in quad_double_diagonal_cascade_solutions, dim1 :', \
            dim1, end='')
        print(', dim2 :', dim2)
    phc = get_phcfun()
    if dim1 >= dim2: 
        adim1 = pointer(c_int32(dim1))
        bdim2 = pointer(c_int32(dim2))
    else:
        adim1 = pointer(c_int32(dim2))
        bdim2 = pointer(c_int32(dim1))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> quad_double_diagonal_cascade_solutions calls phc', end='')
    retval = phc(298, adim1, bdim2, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_start_diagonal_cascade(gamma=0, tasks=0, vrblvl=0):
    r"""
    Does the path tracking to start a diagonal cascade in double
    precision.  For this to work, the functions double_diagonal_homotopy
    and double_diagonal_cascade_solutions must be executed successfully.
    If *gamma* equals 0 on input, then a random gamma constant is generated,
    otherwise, the given complex gamma will be used in the homotopy.
    Multitasking is available, and activated by the *tasks* parameter.
    Returns the target (system and its corresponding) solutions.
    """
    if vrblvl > 0:
        print('in double_start_diagonal_cascade, gamma :', gamma, end='')
        print(', tasks :', tasks)
    usedgamma = set_double_homotopy(gamma, vrblvl=vrblvl-1)
    do_double_track(tasks, vrblvl=vrblvl-1)
    clear_double_solutions(vrblvl-1)
    clear_double_system(vrblvl-1)
    copy_double_system_from_target(vrblvl-1)
    tsys = get_double_system(vrblvl-1)
    sols = get_double_target_solutions(vrblvl-1)
    return (tsys, sols)

def double_double_start_diagonal_cascade(gamma=0, tasks=0, vrblvl=0):
    r"""
    Does the path tracking to start a diagonal cascade in double double
    precision.  For this to work, the double_double_diagonal_homotopy and
    double_double_diagonal_cascade_solutions must be executed successfully.
    If *gamma* equals 0 on input, then a random gamma constant is generated,
    otherwise, the given complex gamma will be used in the homotopy.
    Multitasking is available, and activated by the *tasks* parameter.
    Returns the target (system and its corresponding) solutions.
    """
    if vrblvl > 0:
        print('in double_double_start_diagonal_cascade, gamma :', \
            gamma, end='')
        print(', tasks :', tasks)
    usedgamma = set_double_double_homotopy(gamma, vrblvl=vrblvl-1)
    do_double_double_track(tasks, vrblvl=vrblvl-1)
    clear_double_double_solutions(vrblvl-1)
    clear_double_double_system(vrblvl-1)
    copy_double_double_system_from_target(vrblvl-1)
    tsys = get_double_double_system(vrblvl-1)
    sols = get_double_double_target_solutions(vrblvl-1)
    return (tsys, sols)

def quad_double_start_diagonal_cascade(gamma=0, tasks=0, vrblvl=0):
    r"""
    Does the path tracking to start a diagonal cascade in quad double
    precision.  For this to work, the quad_double_diagonal_homotopy and
    quad_double_diagonal_cascade_solutions must be executed successfully.
    If *gamma* equals 0 on input, then a random gamma constant is generated,
    otherwise, the given complex gamma will be used in the homotopy.
    Multitasking is available, and activated by the *tasks* parameter.
    Returns the target (system and its corresponding) solutions.
    """
    if vrblvl > 0:
        print('in quad_double_start_diagonal_cascade, gamma :', \
            gamma, end='')
        print(', tasks :', tasks)
    usedgamma = set_quad_double_homotopy(gamma, vrblvl=vrblvl-1)
    do_quad_double_track(tasks, vrblvl=vrblvl-1)
    clear_quad_double_solutions(vrblvl-1)
    clear_quad_double_system(vrblvl-1)
    copy_quad_double_system_from_target(vrblvl-1)
    tsys = get_quad_double_system(vrblvl-1)
    sols = get_quad_double_target_solutions(vrblvl-1)
    return (tsys, sols)

def extrinsic_top_diagonal_dimension(ambdim1, ambdim2, setdim1, setdim2, \
    vrblvl=0):
    """
    Returns the dimension of the start and target system to
    start the extrinsic cascade to intersect two witness sets,
    respectively of dimensions setdim1 and setdim2,
    with ambient dimensions respectively equal to ambdim1 and ambdim2.
    """
    if vrblvl > 0:
        print('in extrinsic_top_diagonal_dimension')
        print('ambdim1 :', ambdim1, end='')
        print(', ambdim2 :', ambdim2, end='')
        print(', setdim1 :', setdim1, end='')
        print(', setdim2 :', setdim2)
    phc = get_phcfun()
    alpha = (c_int32 * 2)()
    alpha[0] = ambdim1
    alpha[1] = ambdim2
    palpha = pointer(alpha)
    beta = (c_int32 * 2)()
    beta[0] = setdim1
    beta[1] = setdim2
    pbeta = pointer(beta)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> extrinsic_top_diagonal_dimension calls phc', end='')
    retval = phc(168, palpha, pbeta, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the top diagonal dimension :', palpha[0][0])
    return palpha[0][0]

def double_collapse_diagonal(ksl, dim, vrblvl=0):
    """
    Eliminates the extrinsic diagonal for the system and solutions
    in double precision.
    The number of slack variables in the current embedding is ksl, and
    dim is the number of slack variables to add to the final embedding. 
    """
    if vrblvl > 0:
        print('in double_collapse_diagonal, ksl :', ksl, end='')
        print(', dim :', dim)
    phc = get_phcfun()
    alpha = (c_int32 * 2)()
    alpha[0] = ksl
    alpha[1] = dim
    palpha = pointer(alpha)
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_collapse_diagonal calls phc', end='')
    retval = phc(170, palpha, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_double_collapse_diagonal(ksl, dim, vrblvl=0):
    """
    Eliminates the extrinsic diagonal for the system and solutions
    in double double precision.
    The number of slack variables in the current embedding is ksl, and
    dim is the number of slack variables to add to the final embedding. 
    """
    if vrblvl > 0:
        print('in double_double_collapse_diagonal, ksl :', ksl, end='')
        print(', dim :', dim)
    phc = get_phcfun()
    alpha = (c_int32 * 2)()
    alpha[0] = ksl
    alpha[1] = dim
    palpha = pointer(alpha)
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_double_collapse_diagonal calls phc', end='')
    retval = phc(299, palpha, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def quad_double_collapse_diagonal(ksl, dim, vrblvl=0):
    """
    Eliminates the extrinsic diagonal for the system and solutions
    in quad double precision.
    The number of slack variables in the current embedding is ksl, and
    dim is the number of slack variables to add to the final embedding. 
    """
    if vrblvl > 0:
        print('in quad_double_collapse_diagonal, ksl :', ksl, end='')
        print(', dim :', dim)
    phc = get_phcfun()
    alpha = (c_int32 * 2)()
    alpha[0] = ksl
    alpha[1] = dim
    palpha = pointer(alpha)
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> quad_double_collapse_diagonal calls phc', end='')
    retval = phc(312, palpha, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_diagonal_solver(dim, dm1, sys1, sols1, dm2, sys2, sols2,
    tasks=0, vrblvl=0):
    r"""
    Runs the diagonal homotopies in double precision
    to intersect two witness sets stored in (*sys1*, *sols1*) and
    (*sys2*, *sols2*), of respective dimensions *dm1* and *dm2*.
    The ambient dimension equals *dim*.
    Multitasking is available, and is activated by the *tasks* parameter.
    Returns the last system in the cascade and its solutions.
    """
    if vrblvl > 0:
        print('in double_diagonal_solver, dim :', dim, end='')
        print(', dm1 :', dm1, ', dm2 :', dm2)
        print('the first system :')
        for pol in sys1:
            print(pol)
        print('generic points in the first set :')
        for (idx, sol) in enumerate(sols1):
            print('Solution', idx+1, ':')
            print(sol)
        print('the second system :')
        for pol in sys2:
            print(pol)
        print('generic points in the second set :')
        for (idx, sol) in enumerate(sols2):
            print('Solution', idx+1, ':')
            print(sol)
    topdim = extrinsic_top_diagonal_dimension(dim+dm1, dim+dm2, \
        dm1, dm2, vrblvl)
    if vrblvl > 0:
        print('the top dimension :', topdim)
    kdm = len(sys1) - dm1
    topdiagdim = top_diagonal_dimension(kdm, dm1, dm2)
    if vrblvl > 0:
        print('number of slack variables at the top :', topdiagdim)
    double_diagonal_homotopy(dm1, sys1, sols1, dm2, sys2, sols2, vrblvl)
    double_diagonal_cascade_solutions(dm1, dm2, vrblvl)
    (topsys, startsols) = double_start_diagonal_cascade(tasks=tasks, \
        vrblvl=vrblvl)
    if vrblvl > 0:
        print('the system solved in the start of the cascade :')
        for pol in topsys:
            print(pol)
        print('the solutions after starting the diagonal cascade :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    for k in range(topdiagdim, 0, -1):
        endsols = double_cascade_step(k, topsys, startsols, tasks=tasks, \
            vrblvl=vrblvl)
        if vrblvl > 0:
            print('after running cascade step %d :' % k)
            for (idx, sol) in enumerate(endsols):
                print('Solution', idx+1, ':')
                print(sol)
        endsolsf1 = filter_vanishing(endsols, 1.0e-8)
        if vrblvl > 0: 
            print('computed', len(endsolsf1), 'solutions')
        slack = 'zz' + str(k)
        nbvar = len(topsys)
        endsolsf2 = drop_coordinate_from_double_solutions\
            (endsolsf1, nbvar, slack, vrblvl-1)
        if vrblvl > 0:
            print('after dropping the slack coordinate from the solutions :')
            for (idx, sol) in enumerate(endsolsf2):
                print('Solution', idx+1, ':')
                print(sol)
        nextsys = drop_variable_from_double_polynomials\
            (topsys, slack, vrblvl-1)
        if vrblvl > 0:
            print('after dropping the variable', slack, 'from the system :')
            for pol in nextsys:
                print(pol)
        (topsys, startsols) = (nextsys[:-1], endsolsf2)
    double_collapse_diagonal(0, 0, vrblvl)
    rsys = get_double_system(vrblvl-1)
    sols = get_double_solutions(vrblvl-1)
    return (rsys, sols)

def double_double_diagonal_solver(dim, dm1, sys1, sols1, dm2, sys2, sols2,
    tasks=0, vrblvl=0):
    r"""
    Runs the diagonal homotopies in double double precision
    to intersect two witness sets stored in (*sys1*, *sols1*) and
    (*sys2*, *sols2*), of respective dimensions *dm1* and *dm2*.
    The ambient dimension equals *dim*.
    Multitasking is available, and is activated by the *tasks* parameter.
    Returns the last system in the cascade and its solutions.
    """
    if vrblvl > 0:
        print('in double_double_diagonal_solver, dim :', dim, end='')
        print(', dm1 :', dm1, ', dm2 :', dm2)
        print('the first system :')
        for pol in sys1:
            print(pol)
        print('generic points in the first set :')
        for (idx, sol) in enumerate(sols1):
            print('Solution', idx+1, ':')
            print(sol)
        print('the second system :')
        for pol in sys2:
            print(pol)
        print('generic points in the second set :')
        for (idx, sol) in enumerate(sols2):
            print('Solution', idx+1, ':')
            print(sol)
    topdim = extrinsic_top_diagonal_dimension(dim+dm1, dim+dm2, \
        dm1, dm2, vrblvl)
    if vrblvl > 0:
        print('the top dimension :', topdim)
    kdm = len(sys1) - dm1
    topdiagdim = top_diagonal_dimension(kdm, dm1, dm2)
    if vrblvl > 0:
        print('number of slack variables at the top :', topdiagdim)
    double_double_diagonal_homotopy(dm1, sys1, sols1, dm2, sys2, sols2, vrblvl)
    double_double_diagonal_cascade_solutions(dm1, dm2, vrblvl)
    (topsys, startsols) = double_double_start_diagonal_cascade(tasks=tasks, \
        vrblvl=vrblvl)
    if vrblvl > 0:
        print('the system solved in the start of the cascade :')
        for pol in topsys:
            print(pol)
        print('the solutions after starting the diagonal cascade :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    for k in range(topdiagdim, 0, -1):
        endsols = double_double_cascade_step(k, topsys, startsols, \
            tasks=tasks, vrblvl=vrblvl)
        if vrblvl > 0:
            print('after running cascade step %d :' % k)
            for (idx, sol) in enumerate(endsols):
                print('Solution', idx+1, ':')
                print(sol)
        endsolsf1 = filter_vanishing(endsols, 1.0e-8)
        if vrblvl > 0: 
            print('computed', len(endsolsf1), 'solutions')
        slack = 'zz' + str(k)
        nbvar = len(topsys)
        endsolsf2 = drop_coordinate_from_double_double_solutions\
            (endsolsf1, nbvar, slack, vrblvl-1)
        if vrblvl > 0:
            print('after dropping the slack coordinate from the solutions :')
            for (idx, sol) in enumerate(endsolsf2):
                print('Solution', idx+1, ':')
                print(sol)
        nextsys = drop_variable_from_double_double_polynomials\
            (topsys, slack, vrblvl-1)
        if vrblvl > 0:
            print('after dropping the variable', slack, 'from the system :')
            for pol in nextsys:
                print(pol)
        (topsys, startsols) = (nextsys[:-1], endsolsf2)
    double_double_collapse_diagonal(0, 0, vrblvl)
    rsys = get_double_double_system(vrblvl-1)
    sols = get_double_double_solutions(vrblvl-1)
    return (rsys, sols)

def quad_double_diagonal_solver(dim, dm1, sys1, sols1, dm2, sys2, sols2,
    tasks=0, vrblvl=0):
    r"""
    Runs the diagonal homotopies in quad double precision
    to intersect two witness sets stored in (*sys1*, *sols1*) and
    (*sys2*, *sols2*), of respective dimensions *dm1* and *dm2*.
    The ambient dimension equals *dim*.
    Multitasking is available, and is activated by the *tasks* parameter.
    Returns the last system in the cascade and its solutions.
    """
    if vrblvl > 0:
        print('in quad_double_diagonal_solver, dim :', dim, end='')
        print(', dm1 :', dm1, ', dm2 :', dm2)
        print('the first system :')
        for pol in sys1:
            print(pol)
        print('generic points in the first set :')
        for (idx, sol) in enumerate(sols1):
            print('Solution', idx+1, ':')
            print(sol)
        print('the second system :')
        for pol in sys2:
            print(pol)
        print('generic points in the second set :')
        for (idx, sol) in enumerate(sols2):
            print('Solution', idx+1, ':')
            print(sol)
    topdim = extrinsic_top_diagonal_dimension(dim+dm1, dim+dm2, \
        dm1, dm2, vrblvl)
    if vrblvl > 0:
        print('the top dimension :', topdim)
    kdm = len(sys1) - dm1
    topdiagdim = top_diagonal_dimension(kdm, dm1, dm2)
    if vrblvl > 0:
        print('number of slack variables at the top :', topdiagdim)
    quad_double_diagonal_homotopy(dm1, sys1, sols1, dm2, sys2, sols2, vrblvl)
    quad_double_diagonal_cascade_solutions(dm1, dm2, vrblvl)
    (topsys, startsols) = quad_double_start_diagonal_cascade(tasks=tasks, \
        vrblvl=vrblvl)
    if vrblvl > 0:
        print('the system solved in the start of the cascade :')
        for pol in topsys:
            print(pol)
        print('the solutions after starting the diagonal cascade :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    for k in range(topdiagdim, 0, -1):
        endsols = quad_double_cascade_step(k, topsys, startsols, \
            tasks=tasks, vrblvl=vrblvl)
        if vrblvl > 0:
            print('after running cascade step %d :' % k)
            for (idx, sol) in enumerate(endsols):
                print('Solution', idx+1, ':')
                print(sol)
        endsolsf1 = filter_vanishing(endsols, 1.0e-8)
        if vrblvl > 0: 
            print('computed', len(endsolsf1), 'solutions')
        slack = 'zz' + str(k)
        nbvar = len(topsys)
        endsolsf2 = drop_coordinate_from_quad_double_solutions\
            (endsolsf1, nbvar, slack, vrblvl-1)
        if vrblvl > 0:
            print('after dropping the slack coordinate from the solutions :')
            for (idx, sol) in enumerate(endsolsf2):
                print('Solution', idx+1, ':')
                print(sol)
        nextsys = drop_variable_from_quad_double_polynomials\
            (topsys, slack, vrblvl-1)
        if vrblvl > 0:
            print('after dropping the variable', slack, 'from the system :')
            for pol in nextsys:
                print(pol)
        (topsys, startsols) = (nextsys[:-1], endsolsf2)
    quad_double_collapse_diagonal(0, 0, vrblvl)
    rsys = get_quad_double_system(vrblvl-1)
    sols = get_quad_double_solutions(vrblvl-1)
    return (rsys, sols)

def test_double_diagonal_homotopy(vrblvl=0):
    """
    Tests the diagonal homotopy in double precision.
    """
    hyp1 = 'x1*x2;'
    hyp2 = 'x1 - x2;'
    (w1sys, w1sols) = double_hypersurface_set(2, hyp1, vrblvl-1)
    if vrblvl > 0:
        print('the witness sets for', hyp1)
        for pol in w1sys:
            print(pol)
        for (idx, sol) in enumerate(w1sols):
            print('Solution', idx+1, ':')
            print(sol)
    (w2sys, w2sols) = double_hypersurface_set(2, hyp2, vrblvl-1)
    if vrblvl > 0:
        print('the witness sets for', hyp2)
        for pol in w2sys:
            print(pol)
        for (idx, sol) in enumerate(w2sols):
            print('Solution', idx+1, ':')
            print(sol)
    (psys, sols) = double_diagonal_solver\
        (2, 1, w1sys, w1sols, 1, w2sys, w2sols, 0, vrblvl)
    print('the end system :')
    for pol in psys:
        print(pol)
    print('the solutions of the diagonal solver :')
    for (idx, sol) in enumerate(sols):
        print('Solution', idx+1, ':')
        print(sol)
    return int(len(sols) != 2)

def test_double_double_diagonal_homotopy(vrblvl=0):
    """
    Tests the diagonal homotopy in double double precision.
    """
    hyp1 = 'x1*x2;'
    hyp2 = 'x1 - x2;'
    (w1sys, w1sols) = double_double_hypersurface_set(2, hyp1, vrblvl-1)
    if vrblvl > 0:
        print('the witness sets for', hyp1)
        for pol in w1sys:
            print(pol)
        for (idx, sol) in enumerate(w1sols):
            print('Solution', idx+1, ':')
            print(sol)
    (w2sys, w2sols) = double_double_hypersurface_set(2, hyp2, vrblvl-1)
    if vrblvl > 0:
        print('the witness sets for', hyp2)
        for pol in w2sys:
            print(pol)
        for (idx, sol) in enumerate(w2sols):
            print('Solution', idx+1, ':')
            print(sol)
    (psys, sols) = double_double_diagonal_solver\
        (2, 1, w1sys, w1sols, 1, w2sys, w2sols, tasks=1, vrblvl=vrblvl)
    print('the end system :')
    for pol in psys:
        print(pol)
    print('the solutions of the diagonal solver :')
    for (idx, sol) in enumerate(sols):
        print('Solution', idx+1, ':')
        print(sol)
    return int(len(sols) != 2)

def test_quad_double_diagonal_homotopy(vrblvl=0):
    """
    Tests the diagonal homotopy in quad double precision.
    """
    hyp1 = 'x1*x2;'
    hyp2 = 'x1 - x2;'
    (w1sys, w1sols) = quad_double_hypersurface_set(2, hyp1, vrblvl-1)
    if vrblvl > 0:
        print('the witness sets for', hyp1)
        for pol in w1sys:
            print(pol)
        for (idx, sol) in enumerate(w1sols):
            print('Solution', idx+1, ':')
            print(sol)
    (w2sys, w2sols) = quad_double_hypersurface_set(2, hyp2, vrblvl-1)
    if vrblvl > 0:
        print('the witness sets for', hyp2)
        for pol in w2sys:
            print(pol)
        for (idx, sol) in enumerate(w2sols):
            print('Solution', idx+1, ':')
            print(sol)
    (psys, sols) = quad_double_diagonal_solver\
        (2, 1, w1sys, w1sols, 1, w2sys, w2sols, tasks=1, vrblvl=vrblvl)
    print('the end system :')
    for pol in psys:
        print(pol)
    print('the solutions of the diagonal solver :')
    for (idx, sol) in enumerate(sols):
        print('Solution', idx+1, ':')
        print(sol)
    return int(len(sols) != 2)

def main():
    """
    Runs some tests.
    """
    lvl = 1
    fail = test_double_diagonal_homotopy(lvl)
    fail = fail + test_double_double_diagonal_homotopy(lvl)
    fail = fail + test_quad_double_diagonal_homotopy(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()
