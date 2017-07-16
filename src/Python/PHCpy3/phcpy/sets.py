"""
This module exports routines of PHCpack to manipulate
positive dimensional solution sets of polynomial systems.
"""

def standard_embed(nvar, topdim, pols):
    r"""
    Given in *pols* a list of strings representing polynomials in *nvar*
    variables, with coefficients in standard double precision,
    this function returns an embedding of *pols* of dimension *topdim*.
    The *topdim* is the top dimension which equals the expected highest
    dimension of a component of the solution set of the system of polynomials.
    """
    from phcpy.phcpy2c3 import py2c_embed_standard_system
    from phcpy.interface import store_standard_system, load_standard_system
    store_standard_system(pols, nbvar=nvar)
    py2c_embed_standard_system(topdim)
    return load_standard_system()

def dobldobl_embed(nvar, topdim, pols):
    r"""
    Given in *pols* a list of strings that represent polynomials in *nvar*
    variables, with coefficients in double double precision,
    this function returns an embedding of *pols* of dimension *topdim*.
    The *topdim* is the top dimension which equals the expected highest
    dimension of a component of the solution set of the system of polynomials.
    """
    from phcpy.phcpy2c3 import py2c_embed_dobldobl_system
    from phcpy.interface import store_dobldobl_system, load_dobldobl_system
    store_dobldobl_system(pols, nbvar=nvar)
    py2c_embed_dobldobl_system(topdim)
    return load_dobldobl_system()

def quaddobl_embed(nvar, topdim, pols):
    r"""
    Given in *pols* a list of strings that represent polynomials in *nvar*
    variables, with coefficients in quad double precision,
    this function returns an embedding of *pols* of dimension topdim.
    The *topdim* is the top dimension which equals the expected highest
    dimension of a component of the solution set of the system of polynomials.
    """
    from phcpy.phcpy2c3 import py2c_embed_quaddobl_system
    from phcpy.interface import store_quaddobl_system, load_quaddobl_system
    store_quaddobl_system(pols, nbvar=nvar)
    py2c_embed_quaddobl_system(topdim)
    return load_quaddobl_system()

def standard_laurent_embed(nvar, topdim, pols):
    r"""
    Given in *pols* a list of strings representing Laurent polynomials 
    in *nvar* variables, with coefficients in standard double precision,
    this function returns an embedding of *pols* of dimension *topdim*.
    The *topdim* is the top dimension which equals the expected highest
    dimension of a component of the solution set of the system of polynomials.
    """
    from phcpy.phcpy2c3 import py2c_embed_standard_Laurent_system
    from phcpy.interface import store_standard_laurent_system
    from phcpy.interface import load_standard_laurent_system
    store_standard_laurent_system(pols, nbvar=nvar)
    py2c_embed_standard_Laurent_system(topdim)
    return load_standard_laurent_system()

def dobldobl_laurent_embed(nvar, topdim, pols):
    r"""
    Given in *pols* a list of strings that represent Laurent polynomials 
    in *nvar* variables, with coefficients in double double precision,
    this function returns an embedding of *pols* of dimension *topdim*.
    The *topdim* is the top dimension which equals the expected highest
    dimension of a component of the solution set of the system of polynomials.
    """
    from phcpy.phcpy2c3 import py2c_embed_dobldobl_Laurent_system
    from phcpy.interface import store_dobldobl_laurent_system
    from phcpy.interface import load_dobldobl_laurent_system
    store_dobldobl_laurent_system(pols, nbvar=nvar)
    py2c_embed_dobldobl_Laurent_system(topdim)
    return load_dobldobl_laurent_system()

def quaddobl_laurent_embed(nvar, topdim, pols):
    r"""
    Given in *pols* a list of strings that represent Laurent polynomials 
    in *nvar* variables, with coefficients in quad double precision,
    this function returns an embedding of *pols* of dimension topdim.
    The *topdim* is the top dimension which equals the expected highest
    dimension of a component of the solution set of the system of polynomials.
    """
    from phcpy.phcpy2c3 import py2c_embed_quaddobl_Laurent_system
    from phcpy.interface import store_quaddobl_laurent_system
    from phcpy.interface import load_quaddobl_laurent_system
    store_quaddobl_laurent_system(pols, nbvar=nvar)
    py2c_embed_quaddobl_Laurent_system(topdim)
    return load_quaddobl_laurent_system()

def embed(nvar, topdim, pols, precision='d'):
    r"""
    Given in *pols* a list of strings that represent polynomials in *nvar*
    variables, this function returns an embedding of *pols* 
    of dimension *topdim*.
    The *topdim* is the top dimension which equals the expected highest
    dimension of a component of the solution set of the system of polynomials.
    The default *precision* of the coefficients is 'd', for standard double
    precision.  For double double and quad double precision, set the value
    of *precision* to 'dd' or 'qd' respectively.
    """
    if(precision == 'd'):
        return standard_embed(nvar, topdim, pols)
    elif(precision == 'dd'):
        return dobldobl_embed(nvar, topdim, pols)
    elif(precision == 'qd'):
        return quaddobl_embed(nvar, topdim, pols)
    else:
        print('wrong argument for precision')
        return None

def laurent_embed(nvar, topdim, pols, precision='d'):
    r"""
    Given in *pols* a list of strings that represent Laurent polynomials
    in *nvar* variables, this function returns an embedding of *pols* 
    of dimension *topdim*.
    The *topdim* is the top dimension which equals the expected highest
    dimension of a component of the solution set of the system of polynomials.
    The default *precision* of the coefficients is 'd', for standard double
    precision.  For double double and quad double precision, set the value
    of *precision* to 'dd' or 'qd' respectively.
    """
    if(precision == 'd'):
        return standard_laurent_embed(nvar, topdim, pols)
    elif(precision == 'dd'):
        return dobldobl_laurent_embed(nvar, topdim, pols)
    elif(precision == 'qd'):
        return quaddobl_laurent_embed(nvar, topdim, pols)
    else:
        print('wrong argument for precision')
        return None

def witness_set_of_hypersurface(nvar, hpol, precision='d'):
    r"""
    Given in *hpol* the string representation of a polynomial
    in *nvar* variables (ending with a semicolon),
    on return is an embedded system and its solutions
    which represents a witness set for *hpol*.
    The number of solutions on return should equal
    the degree of the polynomial in *hpol*.
    Three different precisions are supported, by default double ('d'),
    or otherwise double double ('dd') or quad double ('qd').
    """
    if(precision == 'd'):
        from phcpy.phcpy2c3 import py2c_standard_witset_of_hypersurface
        from phcpy.interface import load_standard_system
        from phcpy.interface import load_standard_solutions
        py2c_standard_witset_of_hypersurface(nvar, len(hpol), hpol)
        return (load_standard_system(), load_standard_solutions())
    elif(precision == 'dd'):
        from phcpy.phcpy2c3 import py2c_dobldobl_witset_of_hypersurface
        from phcpy.interface import load_dobldobl_system
        from phcpy.interface import load_dobldobl_solutions
        py2c_dobldobl_witset_of_hypersurface(nvar, len(hpol), hpol)
        return (load_dobldobl_system(), load_dobldobl_solutions())
    elif(precision == 'qd'):
        from phcpy.phcpy2c3 import py2c_quaddobl_witset_of_hypersurface
        from phcpy.interface import load_quaddobl_system
        from phcpy.interface import load_quaddobl_solutions
        py2c_quaddobl_witset_of_hypersurface(nvar, len(hpol), hpol)
        return (load_quaddobl_system(), load_quaddobl_solutions())
    else:
        print('wrong argument for precision')
        return None

def drop_variable_from_standard_polynomials(pols, svar):
    r"""
    Removes the variable with symbol in the string *svar*
    from the list *pols* of strings that represent polynomials
    in several variables, with coefficients in standard double precision.
    Note that the system in *pols* must be square.
    """
    from phcpy.phcpy2c3 import py2c_syscon_standard_drop_variable_by_name
    from phcpy.phcpy2c3 import py2c_syscon_remove_symbol_name
    from phcpy.interface import store_standard_system, load_standard_system
    store_standard_system(pols)
    py2c_syscon_standard_drop_variable_by_name(len(svar), svar)
    py2c_syscon_remove_symbol_name(len(svar), svar)
    return load_standard_system()

def drop_variable_from_dobldobl_polynomials(pols, svar):
    r"""
    Removes the variable with symbol in the string *svar*
    from the list *pols* of strings that represent polynomials
    in several variables, with coefficients in double double precision.
    Note that the system in *pols* must be square.
    """
    from phcpy.phcpy2c3 import py2c_syscon_dobldobl_drop_variable_by_name
    from phcpy.phcpy2c3 import py2c_syscon_remove_symbol_name
    from phcpy.interface import store_dobldobl_system, load_dobldobl_system
    store_dobldobl_system(pols)
    py2c_syscon_dobldobl_drop_variable_by_name(len(svar), svar)
    py2c_syscon_remove_symbol_name(len(svar), svar)
    return load_dobldobl_system()

def drop_variable_from_quaddobl_polynomials(pols, svar):
    r"""
    Removes the variable with symbol in the string *svar*
    from the list *pols* of strings that represent polynomials
    in several variables, with coefficients in quad double precision.
    Note that the system in *pols* must be square.
    """
    from phcpy.phcpy2c3 import py2c_syscon_quaddobl_drop_variable_by_name
    from phcpy.phcpy2c3 import py2c_syscon_remove_symbol_name
    from phcpy.interface import store_quaddobl_system, load_quaddobl_system
    store_quaddobl_system(pols)
    py2c_syscon_quaddobl_drop_variable_by_name(len(svar), svar)
    py2c_syscon_remove_symbol_name(len(svar), svar)
    return load_quaddobl_system()

def drop_coordinate_from_standard_solutions(sols, nbvar, svar):
    r"""
    Removes the variable with symbol in the string *svar*
    from the list *sols* of strings that represent solutions
    in *nbvar* variables, in standard double precision.
    """
    from phcpy.phcpy2c3 import py2c_syscon_clear_symbol_table
    from phcpy.phcpy2c3 import py2c_solcon_standard_drop_coordinate_by_name
    from phcpy.phcpy2c3 import py2c_syscon_remove_symbol_name
    from phcpy.interface import store_standard_solutions
    from phcpy.interface import load_standard_solutions
    py2c_syscon_clear_symbol_table()
    store_standard_solutions(nbvar, sols)
    py2c_solcon_standard_drop_coordinate_by_name(len(svar), svar)
    py2c_syscon_remove_symbol_name(len(svar), svar)
    return load_standard_solutions()

def drop_coordinate_from_dobldobl_solutions(sols, nbvar, svar):
    r"""
    Removes the variable with symbol in the string *svar*
    from the list *sols* of strings that represent solutions
    in *nbvar* variables, in double double precision.
    """
    from phcpy.phcpy2c3 import py2c_syscon_clear_symbol_table
    from phcpy.phcpy2c3 import py2c_solcon_dobldobl_drop_coordinate_by_name
    from phcpy.phcpy2c3 import py2c_syscon_remove_symbol_name
    from phcpy.interface import store_dobldobl_solutions
    from phcpy.interface import load_dobldobl_solutions
    py2c_syscon_clear_symbol_table()
    store_dobldobl_solutions(nbvar, sols)
    py2c_solcon_dobldobl_drop_coordinate_by_name(len(svar), svar)
    py2c_syscon_remove_symbol_name(len(svar), svar)
    return load_dobldobl_solutions()

def drop_coordinate_from_quaddobl_solutions(sols, nbvar, svar):
    r"""
    Removes the variable with symbol in the string *svar*
    from the list *sols* of strings that represent solutions
    in *nbvar* variables, in quad double precision.
    """
    from phcpy.phcpy2c3 import py2c_syscon_clear_symbol_table
    from phcpy.phcpy2c3 import py2c_solcon_quaddobl_drop_coordinate_by_name
    from phcpy.phcpy2c3 import py2c_syscon_remove_symbol_name
    from phcpy.interface import store_quaddobl_solutions
    from phcpy.interface import load_quaddobl_solutions
    py2c_syscon_clear_symbol_table()
    store_quaddobl_solutions(nbvar, sols)
    py2c_solcon_quaddobl_drop_coordinate_by_name(len(svar), svar)
    py2c_syscon_remove_symbol_name(len(svar), svar)
    return load_quaddobl_solutions()

def standard_membertest(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, verbose=True):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedding polynomial
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list *point*,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True.
    Calculations happen in standard double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    """
    from phcpy.interface import store_standard_system as storesys
    from phcpy.interface import store_standard_solutions as storesols
    from phcpy.phcpy2c3 import py2c_witset_standard_membertest as membtest
    storesys(wsys)
    storesols(len(wsys), gpts)
    nvr = len(point)//2
    strpt = str(point)
    nbc = len(strpt)
    result = membtest(int(verbose), nvr, dim, nbc, evatol, memtol, strpt)
    return (result[2] == 1)

def dobldobl_membertest(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, verbose=True):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedding polynomial
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list *point*,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True.
    Calculations happen in double double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    """
    from phcpy.interface import store_dobldobl_system as storesys
    from phcpy.interface import store_dobldobl_solutions as storesols
    from phcpy.phcpy2c3 import py2c_witset_dobldobl_membertest as membtest
    storesys(wsys)
    storesols(len(wsys), gpts)
    nvr = len(point)//4
    strpt = str(point)
    nbc = len(strpt)
    result = membtest(int(verbose), nvr, dim, nbc, evatol, memtol, strpt)
    return (result[2] == 1)

def quaddobl_membertest(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, verbose=True):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedding polynomial
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list point,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True.
    Calculations happen in quad double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    """
    from phcpy.interface import store_quaddobl_system as storesys
    from phcpy.interface import store_quaddobl_solutions as storesols
    from phcpy.phcpy2c3 import py2c_witset_quaddobl_membertest as membtest
    storesys(wsys)
    storesols(len(wsys), gpts)
    nvr = len(point)//8
    strpt = str(point)
    nbc = len(strpt)
    result = membtest(int(verbose), nvr, dim, nbc, evatol, memtol, strpt)
    return (result[2] == 1)

def standard_laurent_membertest(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, verbose=True):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedding Laurent
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list *point*,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True.
    Calculations happen in standard double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    """
    from phcpy.interface import store_standard_laurent_system as storesys
    from phcpy.interface import store_standard_solutions as storesols
    from phcpy.phcpy2c3 \
    import py2c_witset_standard_Laurent_membertest as membtest
    storesys(wsys)
    storesols(len(wsys), gpts)
    nvr = len(point)//2
    strpt = str(point)
    nbc = len(strpt)
    result = membtest(int(verbose), nvr, dim, nbc, evatol, memtol, strpt)
    return (result[2] == 1)

def dobldobl_laurent_membertest(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, verbose=True):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedding Laurent
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list *point*,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True.
    Calculations happen in double double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    """
    from phcpy.interface import store_dobldobl_laurent_system as storesys
    from phcpy.interface import store_dobldobl_solutions as storesols
    from phcpy.phcpy2c3 \
    import py2c_witset_dobldobl_Laurent_membertest as membtest
    storesys(wsys)
    storesols(len(wsys), gpts)
    nvr = len(point)//4
    strpt = str(point)
    nbc = len(strpt)
    result = membtest(int(verbose), nvr, dim, nbc, evatol, memtol, strpt)
    return (result[2] == 1)

def quaddobl_laurent_membertest(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, verbose=True):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedding Laurent
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list point,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True.
    Calculations happen in quad double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    """
    from phcpy.interface import store_quaddobl_laurent_system as storesys
    from phcpy.interface import store_quaddobl_solutions as storesols
    from phcpy.phcpy2c3 \
    import py2c_witset_quaddobl_Laurent_membertest as membtest
    storesys(wsys)
    storesols(len(wsys), gpts)
    nvr = len(point)//8
    strpt = str(point)
    nbc = len(strpt)
    result = membtest(int(verbose), nvr, dim, nbc, evatol, memtol, strpt)
    return (result[2] == 1)

def membertest(wsys, gpts, dim, point, evatol=1.0e-6, memtol=1.0e-6, \
    verbose=True, precision='d'):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedding polynomial
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list *point*,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True, and the
    working *precision* is double 'd'.  Other levels of *precision* are
    double double precision 'dd' and quad double precision 'qd'.
    There are two tolerances: *evatol* is the tolerance on the residual
    of the evaluation of the polynomial equations at the test point.
    If the residual of the evalution is not less than *evatol*,
    then the membertest returns False.  Otherwise, the homotopy
    membership test is called and the *memtol* is used to compare
    the coordinates of the point with the newly computed generic points.
    If there is a match between the coordinates within the given
    tolerance *memtol*, then True is returned.
    """
    if(precision == 'd'):
        return standard_membertest(wsys, gpts, dim, point, \
                                   evatol, memtol, verbose)
    elif(precision == 'dd'):
        return dobldobl_membertest(wsys, gpts, dim, point, \
                                   evatol, memtol, verbose)
    elif(precision == 'qd'):
        return quaddobl_membertest(wsys, gpts, dim, point, \
                                   evatol, memtol, verbose)
    else:
        print('wrong argument for precision')
        return None

def laurent_membertest(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, verbose=True, precision='d'):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedding Laurent
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list *point*,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True, and the
    working *precision* is double 'd'.  Other levels of *precision* are
    double double precision 'dd' and quad double precision 'qd'.
    There are two tolerances: *evatol* is the tolerance on the residual
    of the evaluation of the Laurent system at the test point.
    If the residual of the evalution is not less than *evatol*,
    then the membertest returns False.  Otherwise, the homotopy
    membership test is called and the *memtol* is used to compare
    the coordinates of the point with the newly computed generic points.
    If there is a match between the coordinates within the given
    tolerance *memtol*, then True is returned.
    """
    if(precision == 'd'):
        return standard_laurent_membertest\
                   (wsys, gpts, dim, point, evatol, memtol, verbose)
    elif(precision == 'dd'):
        return dobldobl_membertest\
                   (wsys, gpts, dim, point, evatol, memtol, verbose)
    elif(precision == 'qd'):
        return quaddobl_membertest\
                   (wsys, gpts, dim, point, evatol, memtol, verbose)
    else:
        print('wrong argument for precision')
        return None

def is_slackvar(var):
    r"""
    Given in *var* is a string with a variable name.
    Returns True if the variable name starts with 'zz',
    followed by a number.  Returns False otherwise.
    """
    if len(var) <= 2: # too short to be a slack variable
        return False
    elif var[:2] != 'zz': # a slack variable starts with zz
        return False
    else:
        rest = var[2:]
        return rest.isdigit()

def is_signed(pol):
    r"""
    Given in *pol* is the string representation of a polynomial.
    Returns True if the first nonspace character in the string *pol*
    is either '+' or '-'.  Returns False otherwise.
    """
    wrk = pol.lstrip()
    if len(wrk) == 0: # all characters in pol are spaces
        return False
    else:
        return (wrk[0] == '+' or wrk[0] == '-')

def is_member(wsys, gpts, dim, solpt, evatol=1.0e-6, memtol=1.0e-6, \
    verbose=True, precision='d'):
    r"""
    This function wraps the membertest where the point is a solution,
    given in *solpt*.  All other parameters have the same meaning as
    in the function membertest.
    """
    from phcpy.solutions import strsol2dict, variables
    from phcpy.solver import names_of_variables
    soldc = strsol2dict(solpt)
    svars = variables(soldc)
    pvars = names_of_variables(wsys)
    if verbose:
        print('variables in solution :', svars)
        print('variables in system   :', pvars)
    (ordvar, solval) = ('', [])
    for var in svars:
        if not is_slackvar(var): # skip the slack variables
            if var not in pvars:
                print('Error: mismatch of variables names.')
                return None
            ordvar = ordvar + '+' + var + '-' + var
            solval.append(soldc[var].real)
            if precision == 'dd':
                solval.append(0.0)
            elif precision == 'qd':
                solval.append(0.0)
                solval.append(0.0)
                solval.append(0.0)
            solval.append(soldc[var].imag)
            if precision == 'dd':
                solval.append(0.0)
            elif precision == 'qd':
                solval.append(0.0)
                solval.append(0.0)
                solval.append(0.0)
    if verbose:
        print('order of the variables :', ordvar)
        print('values in the point :', solval)
    wpol0 = wsys[0]
    if is_signed(wpol0):
        newpol0 = ordvar + wpol0
    else:
        newpol0 = ordvar + '+' + wpol0
    newsys = [newpol0] + wsys[1:]
    return membertest(newsys, gpts, dim, solval, evatol, memtol, \
                      verbose, precision)

def test_member(prc='d'):
    """
    To test the membership, we take the twisted cubic.
    """
    twisted = ['x^2 - y;', 'x^3 - z;']
    twiste1 = embed(3, 1, twisted)
    twiste1[0] = 'x + y + z - x - y - z + ' + twiste1[0]
    from phcpy.solver import solve
    twtsols = solve(twiste1, precision=prc)
    for sol in twtsols:
        print(sol)
    if(prc == 'd'):
        inpoint = [1, 0, 1, 0, 1, 0]
        outpoint = [1, 0, 1, 0, 2, 0]
    elif(prc == 'dd'):
        inpoint = [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]
        outpoint = [1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0]
    elif(prc == 'qd'):
        inpoint = [1, 0, 0, 0, 0, 0, 0, 0, \
                   1, 0, 0, 0, 0, 0, 0, 0, \
                   1, 0, 0, 0, 0, 0, 0, 0]
        outpoint = [1, 0, 0, 0, 0, 0, 0, 0, \
                    1, 0, 0, 0, 0, 0, 0, 0, \
                    2, 0, 0, 0, 0, 0, 0, 0, ]
    print(membertest(twiste1, twtsols, 1, inpoint, precision=prc))
    print(membertest(twiste1, twtsols, 1, outpoint, precision=prc))

def test_ismember(prc='d'):
    """
    To test the membertest wrapper, we take the twisted cubic again.
    """
    from phcpy.solver import solve
    from phcpy.solutions import make_solution
    twisted = ['x^2 - y;', 'x^3 - z;']
    twiste1 = embed(3, 1, twisted)
    twtsols = solve(twiste1, precision=prc)
    for sol in twtsols:
        print(sol)
    print(is_member(twiste1, twtsols, 1, twtsols[0], precision=prc))
    outsol = make_solution(['x', 'y', 'z'], [1, 2, 3])
    print(is_member(twiste1, twtsols, 1, outsol, precision=prc))

def test():
    """
    Runs a test on algebraic sets.
    """
    from phcpy.phcpy2c3 import py2c_set_seed
    py2c_set_seed(234798272)
    test_ismember('d')

if __name__ == "__main__":
    test()
