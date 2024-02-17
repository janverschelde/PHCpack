"""
This module exports routines of PHCpack to manipulate
positive dimensional solution sets of polynomial systems.

A numerical data structure to represent a pure dimensional 
solution set is a witness set, which consists of two parts:

1. the original polynomial system, augmented with as many
   random linear equations as the dimension of the set; and

2. isolated solutions of the augmented system,
   as many as the degree of the solution set.

The solutions of the augmented system are called witness points
and are generic points on the algebraic set.

The embed functions add slack variables and random hyperplanes.
The number of slack variables equals the number of random hyperplanes,
which in turn equals the dimension of the solution set.
The drop functions remove the added slack variables from the polynomials
and the coordinates of the solutions.
Given a witness set and a point, a homotopy membership determines whether
the point belongs to the solution set represented by the witness set.
"""
from ctypes import c_int32, c_double, pointer
from phcpy.version import get_phcfun, str2int4a
from phcpy.polynomials import number_of_symbols, string_of_symbols
from phcpy.polynomials import clear_symbol_table
from phcpy.polynomials import set_double_system, get_double_system
from phcpy.polynomials import set_double_double_system
from phcpy.polynomials import get_double_double_system
from phcpy.polynomials import set_quad_double_system
from phcpy.polynomials import get_quad_double_system
from phcpy.polynomials import set_double_laurent_system
from phcpy.polynomials import get_double_laurent_system
from phcpy.polynomials import set_double_double_laurent_system
from phcpy.polynomials import get_double_double_laurent_system
from phcpy.polynomials import set_quad_double_laurent_system
from phcpy.polynomials import get_quad_double_laurent_system
from phcpy.solutions import set_double_solutions
from phcpy.solutions import get_double_solutions
from phcpy.solutions import clear_double_solutions
from phcpy.solutions import set_double_double_solutions
from phcpy.solutions import get_double_double_solutions
from phcpy.solutions import clear_double_double_solutions
from phcpy.solutions import set_quad_double_solutions
from phcpy.solutions import get_quad_double_solutions
from phcpy.solutions import clear_quad_double_solutions
from phcpy.solver import solve

def double_embed(nvr, topdim, pols, vrblvl=0):
    r"""
    Given in *pols* a list of strings representing polynomials in *nvr*
    variables, with coefficients in double precision,
    this function returns an embedding of *pols* of dimension *topdim*.
    The *topdim* is the top dimension which equals the expected highest
    dimension of a component of the solution set of the system of polynomials.
    If *vrblvl* is larger than 0, then extra information is printed.
    """
    if vrblvl > 0:
        print('in double_embed, nvr :', nvr, end='')
        print(', topdim :', topdim)
        print('the polynomials :')
        for pol in pols:
            print(pol)
    set_double_system(nvr, pols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    btopdim = pointer(c_int32(topdim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_embed calls phc', end='')
    retval = phc(66, btopdim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    result = get_double_system(vrblvl-1)
    return result

def double_double_embed(nvr, topdim, pols, vrblvl=0):
    r"""
    Given in *pols* a list of strings representing polynomials in *nvr*
    variables, with coefficients in double double precision,
    this function returns an embedding of *pols* of dimension *topdim*.
    The *topdim* is the top dimension which equals the expected highest
    dimension of a component of the solution set of the system of polynomials.
    If *vrblvl* is larger than 0, then extra information is printed.
    """
    if vrblvl > 0:
        print('in double_double_embed, nvr :', nvr, end='')
        print(', topdim :', topdim)
        print('the polynomials :')
        for pol in pols:
            print(pol)
    set_double_double_system(nvr, pols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    btopdim = pointer(c_int32(topdim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_embed calls phc', end='')
    retval = phc(129, btopdim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    result = get_double_double_system(vrblvl-1)
    return result

def quad_double_embed(nvr, topdim, pols, vrblvl=0):
    r"""
    Given in *pols* a list of strings representing polynomials in *nvr*
    variables, with coefficients in quad double precision,
    this function returns an embedding of *pols* of dimension *topdim*.
    The *topdim* is the top dimension which equals the expected highest
    dimension of a component of the solution set of the system of polynomials.
    If *vrblvl* is larger than 0, then extra information is printed.
    """
    if vrblvl > 0:
        print('in quad_double_embed, nvr :', nvr, end='')
        print(', topdim :', topdim)
        print('the polynomials :')
        for pol in pols:
            print(pol)
    set_quad_double_system(nvr, pols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    btopdim = pointer(c_int32(topdim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> quad_double_embed calls phc', end='')
    retval = phc(260, btopdim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    result = get_quad_double_system(vrblvl-1)
    return result

def double_laurent_embed(nvr, topdim, pols, vrblvl=0):
    r"""
    Given in *pols* a list of strings representing Laurent polynomials 
    in *nvr* variables, with coefficients in double precision,
    this function returns an embedding of *pols* of dimension *topdim*.
    The *topdim* is the top dimension which equals the expected highest
    dimension of a component of the solution set of the system of polynomials.
    If *vrblvl* is larger than 0, then the names of the procedures
    called in the running of the blackbox solver will be listed.
    """
    if vrblvl > 0:
        print('in double_laurent_embed, nvr :', nvr, end='')
        print(', topdim :', topdim)
        print('the polynomials :')
        for pol in pols:
            print(pol)
    set_double_laurent_system(nvr, pols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    btopdim = pointer(c_int32(topdim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_laurent_embed calls phc', end='')
    retval = phc(625, btopdim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    result = get_double_laurent_system(vrblvl-1)
    return result

def double_double_laurent_embed(nvr, topdim, pols, vrblvl=0):
    r"""
    Given in *pols* a list of strings representing Laurent polynomials 
    in *nvr* variables, with coefficients in double double precision,
    this function returns an embedding of *pols* of dimension *topdim*.
    The *topdim* is the top dimension which equals the expected highest
    dimension of a component of the solution set of the system of polynomials.
    If *vrblvl* is larger than 0, then the names of the procedures
    called in the running of the blackbox solver will be listed.
    """
    if vrblvl > 0:
        print('in double_double_laurent_embed, nvr :', nvr, end='')
        print(', topdim :', topdim)
        print('the polynomials :')
        for pol in pols:
            print(pol)
    set_double_double_laurent_system(nvr, pols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    btopdim = pointer(c_int32(topdim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_laurent_embed calls phc', end='')
    retval = phc(626, btopdim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    result = get_double_double_laurent_system(vrblvl-1)
    return result

def quad_double_laurent_embed(nvr, topdim, pols, vrblvl=0):
    r"""
    Given in *pols* a list of strings representing Laurent polynomials 
    in *nvr* variables, with coefficients in quad double precision,
    this function returns an embedding of *pols* of dimension *topdim*.
    The *topdim* is the top dimension which equals the expected highest
    dimension of a component of the solution set of the system of polynomials.
    If *vrblvl* is larger than 0, then the names of the procedures
    called in the running of the blackbox solver will be listed.
    """
    if vrblvl > 0:
        print('in quad_double_laurent_embed, nvr :', nvr, end='')
        print(', topdim :', topdim)
        print('the polynomials :')
        for pol in pols:
            print(pol)
    set_quad_double_laurent_system(nvr, pols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    btopdim = pointer(c_int32(topdim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_laurent_embed calls phc', end='')
    retval = phc(627, btopdim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    result = get_quad_double_laurent_system(vrblvl)
    return result

def set_double_witness_set(nvr, dim, pols, sols, vrblvl=0):
    r"""
    Given in *nvr* is the total number of variables in the list of
    polynomials in *pols* and its list of solutions in *sols*.
    The coefficients in the polynomials and the coordinates of the
    solutions will be parsed and stored in double precision.
    The parameter *dim* equals the number of slack variables used in
    the embedding of *pols* and *sols*.  This *dim* also equals the
    dimension of the solution set represented by the witness set
    given by the lists *pols* and *sols*.
    The symbols for the slack variables are swapped to the end of the
    symbol table in both the polynomials and the solutions.
    """
    if vrblvl > 0:
        print('in set_double_witness_set, nvr :', nvr, end='')
        print(', dim :', dim)
        print('the polynomials :')
        for pol in pols:
            print(pol)
        print('the solutions :')
        for sol in sols:
            print(sol)
    set_double_system(nvr, pols, vrblvl-1)
    set_double_solutions(nvr, sols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    anvr = pointer(c_int32(nvr))
    bdim = pointer(c_int32(dim))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_witness_set calls phc', end='')
    retval = phc(816, anvr, bdim, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_double_witness_set(nvr, dim, pols, sols, vrblvl=0):
    r"""
    Given in *nvr* is the total number of variables in the list of
    polynomials in *pols* and its list of solutions in *sols*.
    The coefficients in the polynomials and the coordinates of the
    solutions will be parsed and stored in double double precision.
    The parameter *dim* equals the number of slack variables used in
    the embedding of *pols* and *sols*.  This *dim* also equals the
    dimension of the solution set represented by the witness set
    given by the lists *pols* and *sols*.
    The symbols for the slack variables are swapped to the end of the
    symbol table in both the polynomials and the solutions.
    """
    if vrblvl > 0:
        print('in set_double_double_witness_set, nvr :', nvr, end='')
        print(', dim :', dim)
        print('the polynomials :')
        for pol in pols:
            print(pol)
        print('the solutions :')
        for sol in sols:
            print(sol)
    set_double_double_system(nvr, pols, vrblvl-1)
    set_double_double_solutions(nvr, sols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    anvr = pointer(c_int32(nvr))
    bdim = pointer(c_int32(dim))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_double_witness_set calls phc', end='')
    retval = phc(817, anvr, bdim, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_quad_double_witness_set(nvr, dim, pols, sols, vrblvl=0):
    r"""
    Given in *nvr* is the total number of variables in the list of
    polynomials in *pols* and its list of solutions in *sols*.
    The coefficients in the polynomials and the coordinates of the
    solutions will be parsed and stored in quad double precision.
    The parameter *dim* equals the number of slack variables used in
    the embedding of *pols* and *sols*.  This *dim* also equals the
    dimension of the solution set represented by the witness set
    given by the lists *pols* and *sols*.
    The symbols for the slack variables are swapped to the end of the
    symbol table in both the polynomials and the solutions.
    """
    if vrblvl > 0:
        print('in set_quad_double_witness_set, nvr :', nvr, end='')
        print(', dim :', dim)
        print('the polynomials :')
        for pol in pols:
            print(pol)
        print('the solutions :')
        for sol in sols:
            print(sol)
    set_quad_double_system(nvr, pols, vrblvl-1)
    set_quad_double_solutions(nvr, sols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    anvr = pointer(c_int32(nvr))
    bdim = pointer(c_int32(dim))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_quad_double_witness_set calls phc', end='')
    retval = phc(818, anvr, bdim, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_laurent_witness_set(nvr, dim, pols, sols, vrblvl=0):
    r"""
    Given in *nvr* is the total number of variables in the list of
    Laurent polynomials in *pols* and its list of solutions in *sols*.
    The coefficients in the Laurent polynomials and the coordinates of
    the solutions will be parsed and stored in double precision.
    The parameter *dim* equals the number of slack variables used in
    the embedding of *pols* and *sols*.  This *dim* also equals the
    dimension of the solution set represented by the witness set
    given by the lists *pols* and *sols*.
    The symbols for the slack variables are swapped to the end of the
    symbol table in both the polynomials and the solutions.
    """
    if vrblvl > 0:
        print('in set_double_laurent_witness_set, nvr :', nvr, end='')
        print(', dim :', dim)
        print('the polynomials :')
        for pol in pols:
            print(pol)
        print('the solutions :')
        for sol in sols:
            print(sol)
    set_double_laurent_system(nvr, pols, vrblvl-1)
    set_double_solutions(nvr, sols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    anvr = pointer(c_int32(nvr))
    bdim = pointer(c_int32(dim))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_laurent_witness_set calls phc', end='')
    retval = phc(819, anvr, bdim, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_double_laurent_witness_set(nvr, dim, pols, sols, vrblvl=0):
    r"""
    Given in *nvr* is the total number of variables in the list of
    Laurent polynomials in *pols* and its list of solutions in *sols*.
    The coefficients in the Laurent polynomials and the coordinates of the
    solutions will be parsed and stored in double double precision.
    The parameter *dim* equals the number of slack variables used in
    the embedding of *pols* and *sols*.  This *dim* also equals the
    dimension of the solution set represented by the witness set
    given by the lists *pols* and *sols*.
    The symbols for the slack variables are swapped to the end of the
    symbol table in both the polynomials and the solutions.
    """
    if vrblvl > 0:
        print('in set_double_double_laurent_witness_set, nvr :', nvr, end='')
        print(', dim :', dim)
        print('the polynomials :')
        for pol in pols:
            print(pol)
        print('the solutions :')
        for sol in sols:
            print(sol)
    set_double_double_laurent_system(nvr, pols, vrblvl-1)
    set_double_double_solutions(nvr, sols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    anvr = pointer(c_int32(nvr))
    bdim = pointer(c_int32(dim))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_double_laurent_witness_set calls phc', end='')
    retval = phc(820, anvr, bdim, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_quad_double_laurent_witness_set(nvr, dim, pols, sols, vrblvl=0):
    r"""
    Given in *nvr* is the total number of variables in the list of
    Laurent polynomials in *pols* and its list of solutions in *sols*.
    The coefficients in the Laurent polynomials and the coordinates of the
    solutions will be parsed and stored in quad double precision.
    The parameter *dim* equals the number of slack variables used in
    the embedding of *pols* and *sols*.  This *dim* also equals the
    dimension of the solution set represented by the witness set
    given by the lists *pols* and *sols*.
    The symbols for the slack variables are swapped to the end of the
    symbol table in both the polynomials and the solutions.
    """
    if vrblvl > 0:
        print('in set_quad_double_laurent_witness_set, nvr :', nvr, end='')
        print(', dim :', dim)
        print('the polynomials :')
        for pol in pols:
            print(pol)
        print('the solutions :')
        for sol in sols:
            print(sol)
    set_quad_double_system(nvr, pols, vrblvl-1)
    set_quad_double_solutions(nvr, sols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    anvr = pointer(c_int32(nvr))
    bdim = pointer(c_int32(dim))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_quad_double_laurent_witness_set calls phc', end='')
    retval = phc(821, anvr, bdim, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def drop_variable_from_table(name, vrblvl=0):
    """
    Drops a variable with the given name from the symbol table.
    The verbose level is given by the value of vrblvl.
    """
    if vrblvl > 0:
        print('in drop_variable_from_table, name :', name)
    phc = get_phcfun(vrblvl-1)
    anbc = pointer(c_int32(len(name)))
    bvar = str2int4a(name, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> drop_variable_from_table calls phc', end='')
    retval = phc(296, anbc, bvar, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def drop_variable_from_double_polynomials(pols, svar, vrblvl=0):
    r"""
    Removes the variable with symbol in the string *svar*
    from the list *pols* of strings that represent polynomials
    in several variables, with coefficients in double precision.
    The system in *pols* must be square.
    """
    if vrblvl > 0:
        print('in drop_variable_from_double_polynomials, svar :', svar)
        print('the polynomials :')
        for pol in pols:
            print(pol)
    set_double_system(len(pols), pols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    anbc = pointer(c_int32(len(svar)))
    bvar = str2int4a(svar, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> drop_variable_from_double_polynomials calls phc', end='')
    retval = phc(309, anbc, bvar, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    drop_variable_from_table(svar, vrblvl-1)
    return get_double_system(vrblvl-1)

def drop_variable_from_double_double_polynomials(pols, svar, vrblvl=0):
    r"""
    Removes the variable with symbol in the string *svar*
    from the list *pols* of strings that represent polynomials
    in several variables, with coefficients in double double precision.
    The system in *pols* must be square.
    """
    if vrblvl > 0:
        print('in drop_variable_from_double_double_polynomials, svar :', svar)
        print('the polynomials :')
        for pol in pols:
            print(pol)
    set_double_double_system(len(pols), pols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    anbc = pointer(c_int32(len(svar)))
    bvar = str2int4a(svar, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> drop_variable_from_double_double_polynomials calls phc', \
            end='')
    retval = phc(310, anbc, bvar, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    drop_variable_from_table(svar, vrblvl-1)
    return get_double_double_system(vrblvl-1)

def drop_variable_from_quad_double_polynomials(pols, svar, vrblvl=0):
    r"""
    Removes the variable with symbol in the string *svar*
    from the list *pols* of strings that represent polynomials
    in several variables, with coefficients in quad double precision.
    The system in *pols* must be square.
    """
    if vrblvl > 0:
        print('in drop_variable_from_quad_double_polynomials, svar :', svar)
        print('the polynomials :')
        for pol in pols:
            print(pol)
    set_quad_double_system(len(pols), pols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    anbc = pointer(c_int32(len(svar)))
    bvar = str2int4a(svar, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> drop_variable_from_quad_double_polynomials calls phc', \
            end='')
    retval = phc(311, anbc, bvar, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    drop_variable_from_table(svar, vrblvl-1)
    return get_quad_double_system(vrblvl-1)

def drop_variable_from_double_laurent_polynomials(pols, svar, vrblvl=0):
    r"""
    Removes the variable with symbol in the string *svar*
    from the list *pols* of strings that represent Laurent polynomials
    in several variables, with coefficients in double precision.
    The system in *pols* must be square.
    """
    if vrblvl > 0:
        print('in drop_variable_from_double_laurent_polynomials,', end='')
        print(', svar :', svar)
        print('the polynomials :')
        for pol in pols:
            print(pol)
    set_double_laurent_system(len(pols), pols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    anbc = pointer(c_int32(len(svar)))
    bvar = str2int4a(svar, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> drop_variable_from_double_laurent_polynomials calls phc', \
            end='')
    retval = phc(831, anbc, bvar, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    drop_variable_from_table(svar, vrblvl-1)
    return get_double_laurent_system(vrblvl-1)

def drop_variable_from_double_double_laurent_polynomials(pols, svar, \
    vrblvl=0):
    r"""
    Removes the variable with symbol in the string *svar*
    from the list *pols* of strings that represent Laurent polynomials
    in several variables, with coefficients in double double precision.
    The system in *pols* must be square.
    """
    if vrblvl > 0:
        print('in drop_variable_from_double_double_laurent_polynomials,', \
            end='')
        print(', svar :', svar)
        print('the polynomials :')
        for pol in pols:
            print(pol)
    set_double_double_laurent_system(len(pols), pols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    anbc = pointer(c_int32(len(svar)))
    bvar = str2int4a(svar, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print( \
        '-> drop_variable_from_double_double_laurent_polynomials calls phc', \
            end='')
    retval = phc(832, anbc, bvar, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    drop_variable_from_table(svar, vrblvl-1)
    return get_double_double_laurent_system(vrblvl-1)

def drop_variable_from_quad_double_laurent_polynomials(pols, svar, \
    vrblvl=0):
    r"""
    Removes the variable with symbol in the string *svar*
    from the list *pols* of strings that represent Laurent polynomials
    in several variables, with coefficients in quad double precision.
    The system in *pols* must be square.
    """
    if vrblvl > 0:
        print('in drop_variable_from_quad_double_laurent_polynomials,', \
            end='')
        print(', svar :', svar)
        print('the polynomials :')
        for pol in pols:
            print(pol)
    set_quad_double_laurent_system(len(pols), pols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    anbc = pointer(c_int32(len(svar)))
    bvar = str2int4a(svar, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print( \
        '-> drop_variable_from_quad_double_laurent_polynomials calls phc', \
            end='')
    retval = phc(833, anbc, bvar, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    drop_variable_from_table(svar, vrblvl-1)
    return get_quad_double_laurent_system(vrblvl-1)

def drop_coordinate_from_double_solutions(sols, nvr, svar, vrblvl=0):
    r"""
    Removes the variable with symbol in the string *svar*
    from the list *sols* of strings that represent solutions
    in *nvr* variables, in double precision.
    """
    if vrblvl > 0:
        print('in drop_coordinate_from_double_solutions', end='')
        print(', nvr :', nvr, end='')
        print(', svar :', svar)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    clear_symbol_table(vrblvl-1)
    clear_double_solutions(vrblvl-1)
    set_double_solutions(nvr, sols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    anbc = pointer(c_int32(len(svar)))
    bvar = str2int4a(svar, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> drop_coordinate_from_double_solutions calls phc', end='')
    retval = phc(146, anbc, bvar, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    drop_variable_from_table(svar, vrblvl-1)
    return get_double_solutions(vrblvl-1)

def drop_coordinate_from_double_double_solutions(sols, nvr, svar, vrblvl=0):
    r"""
    Removes the variable with symbol in the string *svar*
    from the list *sols* of strings that represent solutions
    in *nvr* variables, in double double precision.
    """
    if vrblvl > 0:
        print('in drop_coordinate_from_double_double_solutions', end='')
        print(', nvr :', nvr, end='')
        print(', svar :', svar)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    clear_symbol_table(vrblvl-1)
    clear_double_double_solutions(vrblvl-1)
    set_double_double_solutions(nvr, sols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    anbc = pointer(c_int32(len(svar)))
    bvar = str2int4a(svar, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> drop_coordinate_from_double_double_solutions calls phc', \
            end='')
    retval = phc(349, anbc, bvar, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    drop_variable_from_table(svar, vrblvl-1)
    return get_double_double_solutions(vrblvl-1)

def drop_coordinate_from_quad_double_solutions(sols, nvr, svar, vrblvl=0):
    r"""
    Removes the variable with symbol in the string *svar*
    from the list *sols* of strings that represent solutions
    in *nvr* variables, in quad double precision.
    """
    if vrblvl > 0:
        print('in drop_coordinate_from_quad_double_solutions', end='')
        print(', nvr :', nvr, end='')
        print(', svar :', svar)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    clear_symbol_table(vrblvl-1)
    clear_quad_double_solutions(vrblvl-1)
    set_quad_double_solutions(nvr, sols, vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    anbc = pointer(c_int32(len(svar)))
    bvar = str2int4a(svar, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> drop_coordinate_from_quad_double_solutions calls phc', \
            end='')
    retval = phc(399, anbc, bvar, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    drop_variable_from_table(svar, vrblvl-1)
    return get_quad_double_solutions(vrblvl-1)

def double_membertest(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, tasks=0, vrblvl=0):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedding polynomial
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list *point*,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    """
    if vrblvl > 0:
        print('in double_membertest, dim :', dim)
        print('the polynomials :')
        for pol in wsys:
            print(pol)
        print('the generic points in the witness set :')
        for (idx, sol) in enumerate(gpts):
            print('Generic point', idx+1, ':')
            print(sol)
        print('the point :')
        print(point)
    nvr = number_of_symbols(wsys)
    set_double_witness_set(nvr, dim, wsys, gpts, vrblvl)
    nvr = len(point)//2
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(vrblvl))
    bdims = (c_int32 * 3)()
    bdims[0] = nvr
    bdims[1] = dim
    bdims[2] = tasks
    dims = pointer(bdims)
    size = 2 + 2*nvr
    tolstpt = (c_double * size)()
    tolstpt[0] = evatol
    tolstpt[1] = memtol
    for (idx, crd) in enumerate(point):
        tolstpt[2+idx] = c_double(crd)
    ctoltpt = pointer(tolstpt)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_membertest calls phc', end='')
    retval = phc(537, aaa, dims, ctoltpt, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    onsys = aaa[0]
    vals = dims[:2]
    onset = vals[0][0]
    if vrblvl > 0:
        print('onsys :', onsys, '  onset :', onset)
    return onset == 1

def double_double_membertest(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, tasks=0, vrblvl=0):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedding polynomial
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list *point*,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in double double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    """
    if vrblvl > 0:
        print('in double_double_membertest, dim :', dim)
        print('the polynomials :')
        for pol in wsys:
            print(pol)
        print('the generic points in the witness set :')
        for (idx, sol) in enumerate(gpts):
            print('Generic point', idx+1, ':')
            print(sol)
        print('the point :')
        print(point)
    nvr = number_of_symbols(wsys)
    set_double_double_witness_set(nvr, dim, wsys, gpts, vrblvl)
    nvr = len(point)//4
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(vrblvl))
    bdims = (c_int32 * 3)()
    bdims[0] = nvr
    bdims[1] = dim
    bdims[2] = tasks
    dims = pointer(bdims)
    size = 2 + 4*nvr
    tolstpt = (c_double * size)()
    tolstpt[0] = evatol
    tolstpt[1] = memtol
    for (idx, crd) in enumerate(point):
        tolstpt[2+idx] = c_double(crd)
    ctoltpt = pointer(tolstpt)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_membertest calls phc', end='')
    retval = phc(538, aaa, dims, ctoltpt, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    onsys = aaa[0]
    vals = dims[:2]
    onset = vals[0][0]
    if vrblvl > 0:
        print('onsys :', onsys, '  onset :', onset)
    return onset == 1

def quad_double_membertest(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, tasks=0, vrblvl=0):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedding polynomial
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list *point*,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in quad double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    """
    if vrblvl > 0:
        print('in quad_double_membertest, dim :', dim)
        print('the polynomials :')
        for pol in wsys:
            print(pol)
        print('the generic points in the witness set :')
        for (idx, sol) in enumerate(gpts):
            print('Generic point', idx+1, ':')
            print(sol)
        print('the point :')
        print(point)
    nvr = number_of_symbols(wsys)
    set_quad_double_witness_set(nvr, dim, wsys, gpts, vrblvl)
    nvr = len(point)//8
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(vrblvl))
    bdims = (c_int32 * 3)()
    bdims[0] = nvr
    bdims[1] = dim
    bdims[2] = tasks
    dims = pointer(bdims)
    size = 2 + 8*nvr
    tolstpt = (c_double * size)()
    tolstpt[0] = evatol
    tolstpt[1] = memtol
    for (idx, crd) in enumerate(point):
        tolstpt[2+idx] = c_double(crd)
    ctoltpt = pointer(tolstpt)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_membertest calls phc', end='')
    retval = phc(539, aaa, dims, ctoltpt, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    onsys = aaa[0]
    vals = dims[:2]
    onset = vals[0][0]
    if vrblvl > 0:
        print('onsys :', onsys, '  onset :', onset)
    return onset == 1

def double_laurent_membertest(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, tasks=0, vrblvl=0):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedded Laurent
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list *point*,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    """
    if vrblvl > 0:
        print('in double_laurent_membertest, dim :', dim)
        print('the polynomials :')
        for pol in wsys:
            print(pol)
        print('the generic points in the witness set :')
        for (idx, sol) in enumerate(gpts):
            print('Generic point', idx+1, ':')
            print(sol)
        print('the point :')
        print(point)
    nvr = number_of_symbols(wsys)
    set_double_laurent_witness_set(nvr, dim, wsys, gpts, vrblvl)
    nvr = len(point)//2
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(vrblvl))
    bdims = (c_int32 * 3)()
    bdims[0] = nvr
    bdims[1] = dim
    bdims[2] = tasks
    dims = pointer(bdims)
    size = 2 + 2*nvr
    tolstpt = (c_double * size)()
    tolstpt[0] = evatol
    tolstpt[1] = memtol
    for (idx, crd) in enumerate(point):
        tolstpt[2+idx] = c_double(crd)
    ctoltpt = pointer(tolstpt)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_laurent_membertest calls phc', end='')
    retval = phc(795, aaa, dims, ctoltpt, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    onsys = aaa[0]
    vals = dims[:2]
    onset = vals[0][0]
    if vrblvl > 0:
        print('onsys :', onsys, '  onset :', onset)
    return onset == 1

def double_double_laurent_membertest(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, tasks=0, vrblvl=0):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedded Laurent
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list *point*,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in double double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    """
    if vrblvl > 0:
        print('in double_double_laurent_membertest, dim :', dim)
        print('the polynomials :')
        for pol in wsys:
            print(pol)
        print('the generic points in the witness set :')
        for (idx, sol) in enumerate(gpts):
            print('Generic point', idx+1, ':')
            print(sol)
        print('the point :')
        print(point)
    nvr = number_of_symbols(wsys)
    set_double_double_laurent_witness_set(nvr, dim, wsys, gpts, vrblvl)
    nvr = len(point)//4
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(vrblvl))
    bdims = (c_int32 * 3)()
    bdims[0] = nvr
    bdims[1] = dim
    bdims[2] = tasks
    dims = pointer(bdims)
    size = 2 + 4*nvr
    tolstpt = (c_double * size)()
    tolstpt[0] = evatol
    tolstpt[1] = memtol
    for (idx, crd) in enumerate(point):
        tolstpt[2+idx] = c_double(crd)
    ctoltpt = pointer(tolstpt)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_laurent_membertest calls phc', end='')
    retval = phc(796, aaa, dims, ctoltpt, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    onsys = aaa[0]
    vals = dims[:2]
    onset = vals[0][0]
    if vrblvl > 0:
        print('onsys :', onsys, '  onset :', onset)
    return onset == 1

def quad_double_laurent_membertest(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, tasks=0, vrblvl=0):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedded Laurent
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list *point*,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in quad double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    """
    if vrblvl > 0:
        print('in quad_double_laurent_membertest, dim :', dim)
        print('the polynomials :')
        for pol in wsys:
            print(pol)
        print('the generic points in the witness set :')
        for (idx, sol) in enumerate(gpts):
            print('Generic point', idx+1, ':')
            print(sol)
        print('the point :')
        print(point)
    nvr = number_of_symbols(wsys)
    set_quad_double_laurent_witness_set(nvr, dim, wsys, gpts, vrblvl)
    nvr = len(point)//8
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(vrblvl))
    bdims = (c_int32 * 3)()
    bdims[0] = nvr
    bdims[1] = dim
    bdims[2] = tasks
    dims = pointer(bdims)
    size = 2 + 8*nvr
    tolstpt = (c_double * size)()
    tolstpt[0] = evatol
    tolstpt[1] = memtol
    for (idx, crd) in enumerate(point):
        tolstpt[2+idx] = c_double(crd)
    ctoltpt = pointer(tolstpt)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_laurent_membertest calls phc', end='')
    retval = phc(797, aaa, dims, ctoltpt, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    onsys = aaa[0]
    vals = dims[:2]
    onset = vals[0][0]
    if vrblvl > 0:
        print('onsys :', onsys, '  onset :', onset)
    return onset == 1

def double_hypersurface_set(nvr, hpol, vrblvl=0):
    r"""
    Given in *hpol* the string representation of a polynomial
    in *nvar* variables (ending with a semicolon),
    on return is an embedded system and its solutions
    which represents a witness set for *hpol*.
    The number of solutions on return should equal
    the degree of the polynomial in *hpol*.
    """
    clear_symbol_table(vrblvl-1)
    if vrblvl > 0:
        print('in double_hypersurface_set, nvr :', nvr)
        print('the polynomial :', hpol)
    phc = get_phcfun(vrblvl-1)
    nbc = (c_int32 * 2)()
    nbc[0] = nvr
    nbc[1] = len(hpol)
    anbc = pointer(nbc)
    bpol = str2int4a(hpol, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_hypersurface_set calls phc', end='')
    retval = phc(270, anbc, bpol, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    wsys = get_double_system(vrblvl-1)
    sols = get_double_solutions(vrblvl-1)
    return (wsys, sols)

def double_double_hypersurface_set(nvr, hpol, vrblvl=0):
    r"""
    Given in *hpol* the string representation of a polynomial
    in *nvar* variables (ending with a semicolon),
    on return is an embedded system and its solutions
    which represents a witness set for *hpol*.
    The number of solutions on return should equal
    the degree of the polynomial in *hpol*.
    """
    if vrblvl > 0:
        print('in double_double_hypersurface_set, nvr :', nvr)
        print('the polynomial :', hpol)
    clear_symbol_table(vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    nbc = (c_int32 * 2)()
    nbc[0] = nvr
    nbc[1] = len(hpol)
    anbc = pointer(nbc)
    bpol = str2int4a(hpol, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_hypersurface_set calls phc', end='')
    retval = phc(259, anbc, bpol, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    wsys = get_double_double_system(vrblvl-1)
    sols = get_double_double_solutions(vrblvl-1)
    return (wsys, sols)

def quad_double_hypersurface_set(nvr, hpol, vrblvl=0):
    r"""
    Given in *hpol* the string representation of a polynomial
    in *nvar* variables (ending with a semicolon),
    on return is an embedded system and its solutions
    which represents a witness set for *hpol*.
    The number of solutions on return should equal
    the degree of the polynomial in *hpol*.
    """
    if vrblvl > 0:
        print('in quad_double_hypersurface_set, nvr :', nvr)
        print('the polynomial :', hpol)
    clear_symbol_table(vrblvl-1)
    phc = get_phcfun(vrblvl-1)
    nbc = (c_int32 * 2)()
    nbc[0] = nvr
    nbc[1] = len(hpol)
    anbc = pointer(nbc)
    bpol = str2int4a(hpol, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_hypersurface_set calls phc', end='')
    retval = phc(269, anbc, bpol, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    wsys = get_quad_double_system(vrblvl-1)
    sols = get_quad_double_solutions(vrblvl-1)
    return (wsys, sols)

def test_double_twisted(vrblvl=0):
    """
    Tests the computation of a witness set for the twisted cubic
    in double precision.  As the degree of the twisted cubic is three,
    the number of witness points must equal three as well.
    """
    if vrblvl > 0:
        print('in test_double_twisted ...')
    twisted = ['x^2 - y;', 'x^3 - z;']
    twistede1 = double_embed(3, 1, twisted, vrblvl)
    # trick to order the variables in the solutions
    twistede1[0] = 'x + y + z - x - y - z + ' + twistede1[0]
    if vrblvl > 0:
        print('the embedded system :')
        for pol in twistede1:
            print(pol)
    sols = solve(twistede1, vrblvl=vrblvl)
    if vrblvl > 0:
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    return int(len(sols) != 3)

def test_double_double_twisted(vrblvl=0):
    """
    Tests the computation of a witness set for the twisted cubic
    in double double precision.  As the degree of the twisted cubic
    is three, the number of witness points must equal three as well.
    """
    if vrblvl > 0:
        print('in test_double_double_twisted ...')
    twisted = ['x^2 - y;', 'x^3 - z;']
    twistede1 = double_double_embed(3, 1, twisted, vrblvl)
    # trick to order the variables in the solutions
    twistede1[0] = 'x + y + z - x - y - z + ' + twistede1[0]
    if vrblvl > 0:
        print('the embedded system :')
        for pol in twistede1:
            print(pol)
    sols = solve(twistede1, tasks=1, precision='dd', vrblvl=vrblvl)
    if vrblvl > 0:
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    return int(len(sols) != 3)

def test_quad_double_twisted(vrblvl=0):
    """
    Tests the computation of a witness set for the twisted cubic
    in quad double precision.  As the degree of the twisted cubic
    is three, the number of witness points must equal three as well.
    """
    if vrblvl > 0:
        print('in test_quad_double_twisted ...')
    twisted = ['x^2 - y;', 'x^3 - z;']
    twistede1 = quad_double_embed(3, 1, twisted, vrblvl)
    # trick to order the variables in the solutions
    twistede1[0] = 'x + y + z - x - y - z + ' + twistede1[0]
    if vrblvl > 0:
        print('the embedded system :')
        for pol in twistede1:
            print(pol)
    sols = solve(twistede1, tasks=1, precision='qd', vrblvl=vrblvl)
    if vrblvl > 0:
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    return int(len(sols) != 3)

def test_double_laurent_twisted(vrblvl=0):
    """
    Tests the computation of a witness set for the twisted cubic
    defined as a laurent system in double precision.
    As the degree of the twisted cubic is three,
    the number of witness points must equal three as well.
    """
    if vrblvl > 0:
        print('in test_double_laurent_twisted ...')
    twisted = ['x^2*y^(-1) - 1;', 'x^3*z^(-1) - 1;']
    twistede1 = double_laurent_embed(3, 1, twisted, vrblvl)
    # trick to order the variables in the solutions
    twistede1[0] = 'x + y + z - x - y - z + ' + twistede1[0]
    if vrblvl > 0:
        print('the embedded system :')
        for pol in twistede1:
            print(pol)
    # second trick to avoid single monomial equation
    twistede1[2] = twistede1[2][:-1] + twistede1[3]
    sols = solve(twistede1, vrblvl=vrblvl)
    if vrblvl > 0:
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    return int(len(sols) != 3)

def test_double_double_laurent_twisted(vrblvl=0):
    """
    Tests the computation of a witness set for the twisted cubic
    defined as a laurent system in double double precision.
    As the degree of the twisted cubic is three,
    the number of witness points must equal three as well.
    """
    if vrblvl > 0:
        print('in test_double_double_laurent_twisted ...')
    twisted = ['x^2*y^(-1) - 1;', 'x^3*z^(-1) - 1;']
    twistede1 = double_double_laurent_embed(3, 1, twisted, vrblvl)
    # trick to order the variables in the solutions
    twistede1[0] = 'x + y + z - x - y - z + ' + twistede1[0]
    if vrblvl > 0:
        print('the embedded system :')
        for pol in twistede1:
            print(pol)
    # second trick to avoid single monomial equation
    twistede1[2] = twistede1[2][:-1] + twistede1[3]
    sols = solve(twistede1, tasks=1, precision='dd', vrblvl=vrblvl)
    if vrblvl > 0:
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    return int(len(sols) != 3)

def test_quad_double_laurent_twisted(vrblvl=0):
    """
    Tests the computation of a witness set for the twisted cubic
    defined as a laurent system in quad double precision.
    As the degree of the twisted cubic is three,
    the number of witness points must equal three as well.
    """
    if vrblvl > 0:
        print('in test_quad_double_laurent_twisted ...')
    twisted = ['x^2*y^(-1) - 1;', 'x^3*z^(-1) - 1;']
    twistede1 = quad_double_laurent_embed(3, 1, twisted, vrblvl)
    # trick to order the variables in the solutions
    twistede1[0] = 'x + y + z - x - y - z + ' + twistede1[0]
    if vrblvl > 0:
        print('the embedded system :')
        for pol in twistede1:
            print(pol)
    # second trick to avoid single monomial equation
    twistede1[2] = twistede1[2][:-1] + twistede1[3]
    sols = solve(twistede1, tasks=1, precision='qd', vrblvl=vrblvl)
    if vrblvl > 0:
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    return int(len(sols) != 3)

def test_double_member(vrblvl=0):
    """
    Tests the membership in double precision.
    """
    if vrblvl > 0:
        print('in test_double_member ...')
    twisted = ['x^2 - y;', 'x^3 - z;']
    twistede1 = double_embed(3, 1, twisted, vrblvl=vrblvl)
    twistede1[0] = 'x + y + z - x - y - z + ' + twistede1[0]
    sols = solve(twistede1, vrblvl=vrblvl)
    inpoint = [1, 0, 1, 0, 1, 0]
    ismbr = double_membertest(twistede1, sols, 1, inpoint, vrblvl=vrblvl)
    if vrblvl > 0:
        print('on point in the set :', ismbr)
    fail = int(ismbr is not True)
    outpoint = [1, 0, 1, 0, 2, 0]
    ismbr = double_membertest(twistede1, sols, 1, outpoint, vrblvl=vrblvl)
    if vrblvl > 0:
        print('on point not in the set :', ismbr)
    fail = fail + int(ismbr is not False)
    return fail

def test_double_double_member(vrblvl=0):
    """
    Tests the membership in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_member ...')
    twisted = ['x^2 - y;', 'x^3 - z;']
    twistede1 = double_double_embed(3, 1, twisted, vrblvl=vrblvl)
    twistede1[0] = 'x + y + z - x - y - z + ' + twistede1[0]
    sols = solve(twistede1, tasks=1, precision='dd', vrblvl=vrblvl)
    inpoint = [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]
    ismbr = double_double_membertest(twistede1, sols, 1, inpoint, \
        vrblvl=vrblvl)
    if vrblvl > 0:
        print('on point in the set :', ismbr)
    fail = int(ismbr is not True)
    outpoint = [1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0]
    ismbr = double_double_membertest(twistede1, sols, 1, outpoint, \
        vrblvl=vrblvl)
    if vrblvl > 0:
        print('on point not in the set :', ismbr)
    fail = fail + int(ismbr is not False)
    return fail

def test_quad_double_member(vrblvl=0):
    """
    Tests the membership in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_member ...')
    twisted = ['x^2 - y;', 'x^3 - z;']
    twistede1 = quad_double_embed(3, 1, twisted, vrblvl=vrblvl)
    twistede1[0] = 'x + y + z - x - y - z + ' + twistede1[0]
    sols = solve(twistede1, tasks=1, precision='qd', vrblvl=vrblvl)
    inpoint = [1, 0, 0, 0, 0, 0, 0, 0, \
               1, 0, 0, 0, 0, 0, 0, 0, \
               1, 0, 0, 0, 0, 0, 0, 0 ]
    ismbr = quad_double_membertest(twistede1, sols, 1, inpoint, \
        vrblvl=vrblvl)
    if vrblvl > 0:
        print('on point in the set :', ismbr)
    fail = int(ismbr is not True)
    outpoint = [1, 0, 0, 0, 0, 0, 0, 0, \
                1, 0, 0, 0, 0, 0, 0, 0, \
                2, 0, 0, 0, 0, 0, 0, 0]
    ismbr = quad_double_membertest(twistede1, sols, 1, outpoint, \
        vrblvl=vrblvl)
    if vrblvl > 0:
        print('on point not in the set :', ismbr)
    fail = fail + int(ismbr is not False)
    return fail

def test_double_drop(vrblvl=0):
    """
    Tests the removal of a slack variable in double precision.
    """
    if vrblvl > 0:
        print('in test_double_drop ...')
    twisted = ['x^2 - y;', 'x^3 - z;']
    twistede1 = double_embed(3, 1, twisted, vrblvl=vrblvl)
    # trick needed to move zz1 to the end
    twistede1[0] = 'x + y + z - x - y - z + ' + twistede1[0]
    if vrblvl > 0:
        print('the embedded system :')
        for pol in twistede1:
            print(pol)
    clear_double_solutions(vrblvl)
    sols = solve(twistede1, vrblvl=vrblvl)
    if vrblvl > 0:
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
        print('the symbols :', string_of_symbols())
    name = 'zz1'
    dropped = drop_variable_from_double_polynomials(twistede1, name, vrblvl)
    if vrblvl > 0:
        print('the polynomials after the drop :')
        for pol in dropped:
            print(pol)
        print('the symbols :', string_of_symbols())
    dropsols = drop_coordinate_from_double_solutions(sols, 4, name, vrblvl)
    if vrblvl > 0:
        print('the solutions after the drop :')
        for (idx, sol) in enumerate(dropsols):
            print('Solution', idx+1, ':')
            print(sol)
    return 0

def test_double_double_drop(vrblvl=0):
    """
    Tests the removal of a slack variable in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_drop ...')
    twisted = ['x^2 - y;', 'x^3 - z;']
    twistede1 = double_double_embed(3, 1, twisted, vrblvl=vrblvl)
    # trick needed to move zz1 to the end
    twistede1[0] = 'x + y + z - x - y - z + ' + twistede1[0]
    if vrblvl > 0:
        print('the embedded system :')
        for pol in twistede1:
            print(pol)
    clear_double_double_solutions(vrblvl)
    sols = solve(twistede1, tasks=1, precision='dd', vrblvl=vrblvl)
    if vrblvl > 0:
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
        print('the symbols :', string_of_symbols())
    name = 'zz1'
    dropped = drop_variable_from_double_double_polynomials(twistede1, \
        name, vrblvl)
    if vrblvl > 0:
        print('the polynomials after the drop :')
        for pol in dropped:
            print(pol)
        print('the symbols :', string_of_symbols())
    dropsols = drop_coordinate_from_double_double_solutions(sols, 4, \
        name, vrblvl)
    if vrblvl > 0:
        print('the solutions after the drop :')
        for (idx, sol) in enumerate(dropsols):
            print('Solution', idx+1, ':')
            print(sol)
    return 0

def test_quad_double_drop(vrblvl=0):
    """
    Tests the removal of a slack variable in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_drop ...')
    twisted = ['x^2 - y;', 'x^3 - z;']
    twistede1 = quad_double_embed(3, 1, twisted, vrblvl=vrblvl)
    # trick needed to move zz1 to the end
    twistede1[0] = 'x + y + z - x - y - z + ' + twistede1[0]
    if vrblvl > 0:
        print('the embedded system :')
        for pol in twistede1:
            print(pol)
    clear_quad_double_solutions()
    sols = solve(twistede1, precision='qd', vrblvl=vrblvl)
    if vrblvl > 0:
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
        print('the symbols :', string_of_symbols())
    name = 'zz1'
    dropped = drop_variable_from_quad_double_polynomials(twistede1, \
        name, vrblvl)
    if vrblvl > 0:
        print('after the drop :')
        for pol in dropped:
            print(pol)
        print('the symbols :', string_of_symbols())
    dropsols = drop_coordinate_from_quad_double_solutions(sols, 4, \
        name, vrblvl)
    if vrblvl > 0:
        print('the solutions after the drop :')
        for (idx, sol) in enumerate(dropsols):
            print('Solution', idx+1, ':')
            print(sol)
    return 0

def test_double_hypersurface_set(vrblvl=0):
    """
    Tests the construction of a witness set of a hypersurface
    in double precision.
    """
    if vrblvl > 0:
        print('in test_double_hypersurface_set ...')
    hyp = 'x*y*z + 3*x^2 - 1;'
    (pols, sols) = double_hypersurface_set(3, hyp, vrblvl)
    if vrblvl > 0:
        print('the polynomials :')
        for pol in pols:
            print(pol)
        print('the generic points :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    return int(len(sols) != 3)

def test_double_double_hypersurface_set(vrblvl=0):
    """
    Tests the construction of a witness set of a hypersurface
    in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_hypersurface_set ...')
    hyp = 'x*y*z + 3*x^2 - 1;'
    (pols, sols) = double_double_hypersurface_set(3, hyp, vrblvl)
    if vrblvl > 0:
        print('the polynomials :')
        for pol in pols:
            print(pol)
        print('the generic points :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    return int(len(sols) != 3)

def test_quad_double_hypersurface_set(vrblvl=0):
    """
    Tests the construction of a witness set of a hypersurface
    in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_hypersurface_set ...')
    hyp = 'x*y*z + 3*x^2 - 1;'
    (pols, sols) = quad_double_hypersurface_set(3, hyp, vrblvl)
    if vrblvl > 0:
        print('the polynomials :')
        for pol in pols:
            print(pol)
        print('the generic points :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    return int(len(sols) != 3)

def main():
    """
    Runs some tests.
    """
    lvl = 1
    fail = test_double_twisted(lvl)
    fail = fail + test_double_double_twisted(lvl)
    fail = fail + test_quad_double_twisted(lvl)
    fail = fail + test_double_laurent_twisted(lvl)
    fail = fail + test_double_double_laurent_twisted(lvl)
    fail = fail + test_quad_double_laurent_twisted(lvl)
    fail = fail + test_double_member(lvl)
    fail = fail + test_double_double_member(lvl)
    fail = fail + test_quad_double_member(lvl)
    fail = fail + test_double_drop(lvl)
    fail = fail + test_double_double_drop(lvl)
    fail = fail + test_quad_double_drop(lvl)
    fail = fail + test_double_hypersurface_set(lvl)
    fail = fail + test_double_double_hypersurface_set(lvl)
    fail = fail + test_quad_double_hypersurface_set(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()
