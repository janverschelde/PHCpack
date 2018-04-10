"""
This module exports routines of PHCpack to manipulate
positive dimensional solution sets of polynomial systems.
The embed functions add slack variables and random hyperplanes.
The number of slack variables equals the number of random hyperplanes,
which in turn equals the dimension of the solution set.
The drop functions remove the added slack variables from the polynomials
and the coordinates of the solutions.
Given a witness set and a point, a homotopy membership determines whether
the point belongs to the solution set represented by the witness set.
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

def witness_set_of_laurent_hypersurface(nvar, hpol, precision='d'):
    r"""
    Given in *hpol* the string representation of a laurent polynomial
    in *nvar* variables (ending with a semicolon),
    on return is an embedded system and its solutions
    which represents a witness set for *hpol*.
    The number of solutions on return may differ from the actual degree
    of hpol if the polynomial represented by hpol has negative exponents.
    Three different precisions are supported, by default double ('d'),
    or otherwise double double ('dd') or quad double ('qd').
    """
    if(precision == 'd'):
        from phcpy.phcpy2c3 import py2c_standard_witset_of_Laurent_hypersurface
        from phcpy.interface import load_standard_laurent_system
        from phcpy.interface import load_standard_solutions
        py2c_standard_witset_of_Laurent_hypersurface(nvar, len(hpol), hpol)
        return (load_standard_laurent_system(), load_standard_solutions())
    elif(precision == 'dd'):
        from phcpy.phcpy2c3 import py2c_dobldobl_witset_of_Laurent_hypersurface
        from phcpy.interface import load_dobldobl_laurent_system
        from phcpy.interface import load_dobldobl_solutions
        py2c_dobldobl_witset_of_Laurent_hypersurface(nvar, len(hpol), hpol)
        return (load_dobldobl_laurent_system(), load_dobldobl_solutions())
    elif(precision == 'qd'):
        from phcpy.phcpy2c3 import py2c_quaddobl_witset_of_Laurent_hypersurface
        from phcpy.interface import load_quaddobl_laurent_system
        from phcpy.interface import load_quaddobl_solutions
        py2c_quaddobl_witset_of_Laurent_hypersurface(nvar, len(hpol), hpol)
        return (load_quaddobl_laurent_system(), load_quaddobl_solutions())
    else:
        print('wrong argument for precision')

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

def drop_variable_from_standard_laurent_polynomials(pols, svar):
    r"""
    Removes the variable with symbol in the string *svar*
    from the list *pols* of strings that represent Laurent polynomials
    in several variables, with coefficients in standard double precision.
    Note that the system in *pols* must be square.
    """
    from phcpy.phcpy2c3 \
    import py2c_syscon_standard_Laurent_drop_variable_by_name
    from phcpy.phcpy2c3 import py2c_syscon_remove_symbol_name
    from phcpy.interface import store_standard_laurent_system
    from phcpy.interface import load_standard_laurent_system
    store_standard_laurent_system(pols)
    py2c_syscon_standard_Laurent_drop_variable_by_name(len(svar), svar)
    py2c_syscon_remove_symbol_name(len(svar), svar)
    return load_standard_laurent_system()

def drop_variable_from_dobldobl_laurent_polynomials(pols, svar):
    r"""
    Removes the variable with symbol in the string *svar*
    from the list *pols* of strings that represent Laurent polynomials
    in several variables, with coefficients in double double precision.
    Note that the system in *pols* must be square.
    """
    from phcpy.phcpy2c3 \
    import py2c_syscon_dobldobl_Laurent_drop_variable_by_name
    from phcpy.phcpy2c3 import py2c_syscon_remove_symbol_name
    from phcpy.interface import store_dobldobl_laurent_system
    from phcpy.interface import load_dobldobl_laurent_system
    store_dobldobl_laurent_system(pols)
    py2c_syscon_dobldobl_Laurent_drop_variable_by_name(len(svar), svar)
    py2c_syscon_remove_symbol_name(len(svar), svar)
    return load_dobldobl_laurent_system()

def drop_variable_from_quaddobl_laurent_polynomials(pols, svar):
    r"""
    Removes the variable with symbol in the string *svar*
    from the list *pols* of strings that represent Laurent polynomials
    in several variables, with coefficients in quad double precision.
    Note that the system in *pols* must be square.
    """
    from phcpy.phcpy2c3 \
    import py2c_syscon_quaddobl_Laurent_drop_variable_by_name
    from phcpy.phcpy2c3 import py2c_syscon_remove_symbol_name
    from phcpy.interface import store_quaddobl_laurent_system
    from phcpy.interface import load_quaddobl_laurent_system
    store_quaddobl_laurent_system(pols)
    py2c_syscon_quaddobl_Laurent_drop_variable_by_name(len(svar), svar)
    py2c_syscon_remove_symbol_name(len(svar), svar)
    return load_quaddobl_laurent_system()

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
    evatol=1.0e-6, memtol=1.0e-6, verbose=True, tasks=0):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedding polynomial
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list *point*,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in standard double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    """
    from phcpy.interface import store_standard_witness_set
    from phcpy.phcpy2c3 import py2c_witset_standard_membertest as membtest
    store_standard_witness_set(len(wsys), dim, wsys, gpts)
    nvr = len(point)//2
    spt = str(point)
    nbc = len(spt)
    result = membtest(int(verbose), tasks, nvr, dim, nbc, evatol, memtol, spt)
    return (result[2] == 1)

def dobldobl_membertest(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, verbose=True, tasks=0):
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
    from phcpy.interface import store_dobldobl_witness_set
    from phcpy.phcpy2c3 import py2c_witset_dobldobl_membertest as membtest
    store_dobldobl_witness_set(len(wsys), dim, wsys, gpts)
    nvr = len(point)//4
    spt = str(point)
    nbc = len(spt)
    result = membtest(int(verbose), tasks, nvr, dim, nbc, evatol, memtol, spt)
    return (result[2] == 1)

def quaddobl_membertest(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, verbose=True, tasks=0):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedding polynomial
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list point,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in quad double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    """
    from phcpy.interface import store_quaddobl_witness_set
    from phcpy.phcpy2c3 import py2c_witset_quaddobl_membertest as membtest
    store_quaddobl_witness_set(len(wsys), dim, wsys, gpts)
    nvr = len(point)//8
    spt = str(point)
    nbc = len(spt)
    result = membtest(int(verbose), tasks, nvr, dim, nbc, evatol, memtol, spt)
    return (result[2] == 1)

def standard_laurent_membertest(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, verbose=True, tasks=0):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedding Laurent
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list *point*,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in standard double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    """
    from phcpy.interface import store_standard_laurent_witness_set
    from phcpy.phcpy2c3 \
    import py2c_witset_standard_Laurent_membertest as membtest
    store_standard_laurent_witness_set(len(wsys), dim, wsys, gpts)
    nvr = len(point)//2
    spt = str(point)
    nbc = len(spt)
    result = membtest(int(verbose), tasks, nvr, dim, nbc, evatol, memtol, spt)
    return (result[2] == 1)

def dobldobl_laurent_membertest(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, verbose=True, tasks=0):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedding Laurent
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
    from phcpy.interface import store_dobldobl_laurent_witness_set
    from phcpy.phcpy2c3 \
    import py2c_witset_dobldobl_Laurent_membertest as membtest
    store_dobldobl_laurent_witness_set(len(wsys), dim, wsys, gpts)
    nvr = len(point)//4
    spt = str(point)
    nbc = len(spt)
    result = membtest(int(verbose), tasks, nvr, dim, nbc, evatol, memtol, spt)
    return (result[2] == 1)

def quaddobl_laurent_membertest(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, verbose=True, tasks=0):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedding Laurent
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list point,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in quad double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    """
    from phcpy.interface import store_quaddobl_laurent_witness_set
    from phcpy.phcpy2c3 \
    import py2c_witset_quaddobl_Laurent_membertest as membtest
    store_quaddobl_laurent_witness_set(len(wsys), dim, wsys, gpts)
    nvr = len(point)//8
    spt = str(point)
    nbc = len(spt)
    result = membtest(int(verbose), tasks, nvr, dim, nbc, evatol, memtol, spt)
    return (result[2] == 1)

def membertest(wsys, gpts, dim, point, evatol=1.0e-6, memtol=1.0e-6, \
    verbose=True, precision='d', tasks=0):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedding polynomial
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list *point*,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    The default working *precision* is double 'd'.  Other levels of precision
    are double double precision 'dd' and quad double precision 'qd'.
    There are two tolerances: *evatol* is the tolerance on the residual
    of the evaluation of the polynomial equations at the test point.
    If the residual of the evaluation is not less than *evatol*,
    then the membertest returns False.  Otherwise, the homotopy
    membership test is called and the *memtol* is used to compare
    the coordinates of the point with the newly computed generic points.
    If there is a match between the coordinates within the given
    tolerance *memtol*, then True is returned.
    """
    if(precision == 'd'):
        return standard_membertest(wsys, gpts, dim, point, \
                                   evatol, memtol, verbose, tasks)
    elif(precision == 'dd'):
        return dobldobl_membertest(wsys, gpts, dim, point, \
                                   evatol, memtol, verbose, tasks)
    elif(precision == 'qd'):
        return quaddobl_membertest(wsys, gpts, dim, point, \
                                   evatol, memtol, verbose, tasks)
    else:
        print('wrong argument for precision')
        return None

def laurent_membertest(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, verbose=True, precision='d', tasks=0):
    r"""
    Applies the homotopy membership test for a point to belong to
    a witness set of dimension *dim*, given by an embedding Laurent
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the list *point*,
    as a list of doubles, with the real and imaginary part of each
    coordinate of the point.  By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    The default working *precision* is double 'd'.  Other levels of precision
    are double double precision 'dd' and quad double precision 'qd'.
    There are two tolerances: *evatol* is the tolerance on the residual
    of the evaluation of the Laurent system at the test point.
    If the residual of the evaluation is not less than *evatol*,
    then the membertest returns False.  Otherwise, the homotopy
    membership test is called and the *memtol* is used to compare
    the coordinates of the point with the newly computed generic points.
    If there is a match between the coordinates within the given
    tolerance *memtol*, then True is returned.
    """
    if(precision == 'd'):
        return standard_laurent_membertest\
                   (wsys, gpts, dim, point, evatol, memtol, verbose, tasks)
    elif(precision == 'dd'):
        return dobldobl_membertest\
                   (wsys, gpts, dim, point, evatol, memtol, verbose, tasks)
    elif(precision == 'qd'):
        return quaddobl_membertest\
                   (wsys, gpts, dim, point, evatol, memtol, verbose, tasks)
    else:
        print('wrong argument for precision')
        return None

def standard_ismember(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, verbose=True, tasks=0):
    r"""
    Applies the homotopy membership test for a *point* to belong to
    a witness set of dimension *dim*, given by an embedding polynomial
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the string *point*,
    which is the string representation of a solution in PHCpack format,
    with symbols of the variables before the values of the coordinates.
    By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in standard double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    Returns a tuple of two booleans.  The first boolean is True if the
    point satisfies the equations, otherwise it is False.  The second
    boolean is True if the point belongs to the witness set,
    otherwise, the second boolean is False.
    """
    from phcpy.interface import store_standard_witness_set
    from phcpy.phcpy2c3 import py2c_witset_standard_ismember as membtest
    store_standard_witness_set(len(wsys), dim, wsys, gpts)
    nbc = len(point)
    nvr = len(wsys) - dim # test point should not have slack variables
    if verbose:
        print('calling standard_ismember with test point :')
        print(point)
    result = membtest(int(verbose), tasks, nvr, dim, nbc, evatol, memtol, point)
    (fail, onpolsys, inwitset) = result
    return (onpolsys, inwitset)

def dobldobl_ismember(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, verbose=True, tasks=0):
    r"""
    Applies the homotopy membership test for a *point* to belong to
    a witness set of dimension *dim*, given by an embedding polynomial
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the string *point*,
    which is the string representation of a solution in PHCpack format,
    with symbols of the variables before the values of the coordinates.
    By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in double double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    Returns a tuple of two booleans.  The first boolean is True if the
    point satisfies the equations, otherwise it is False.  The second
    boolean is True if the point belongs to the witness set,
    otherwise, the second boolean is False.
    """
    from phcpy.interface import store_dobldobl_witness_set
    from phcpy.phcpy2c3 import py2c_witset_dobldobl_ismember as membtest
    store_dobldobl_witness_set(len(wsys), dim, wsys, gpts)
    nbc = len(point)
    nvr = len(wsys) - dim # test point should not have slack variables
    if verbose:
        print('calling dobldobl_ismember with test point :')
        print(point)
    result = membtest(int(verbose), tasks, nvr, dim, nbc, evatol, memtol, point)
    (fail, onpolsys, inwitset) = result
    return (onpolsys, inwitset)

def quaddobl_ismember(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, verbose=True, tasks=0):
    r"""
    Applies the homotopy membership test for a *point* to belong to
    a witness set of dimension *dim*, given by an embedding polynomial
    system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the string *point*,
    which is the string representation of a solution in PHCpack format,
    with symbols of the variables before the values of the coordinates.
    By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in quad double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    Returns a tuple of two booleans.  The first boolean is True if the
    point satisfies the equations, otherwise it is False.  The second
    boolean is True if the point belongs to the witness set,
    otherwise, the second boolean is False.
    """
    from phcpy.interface import store_quaddobl_witness_set
    from phcpy.phcpy2c3 import py2c_witset_quaddobl_ismember as membtest
    store_quaddobl_witness_set(len(wsys), dim, wsys, gpts)
    nbc = len(point)
    nvr = len(wsys) - dim # test point should not have slack variables
    if verbose:
        print('calling quaddobl_ismember with test point :')
        print(point)
    result = membtest(int(verbose), tasks, nvr, dim, nbc, evatol, memtol, point)
    (fail, onpolsys, inwitset) = result
    return (onpolsys, inwitset)

def standard_laurent_ismember(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, verbose=True, tasks=0):
    r"""
    Applies the homotopy membership test for a *point* to belong to
    a witness set of dimension *dim*, given by an embedding polynomial
    Laurent system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the string *point*,
    which is the string representation of a solution in PHCpack format,
    with symbols of the variables before the values of the coordinates.
    By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in standard double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    Returns a tuple of two booleans.  The first boolean is True if the
    point satisfies the equations, otherwise it is False.  The second
    boolean is True if the point belongs to the witness set,
    otherwise, the second boolean is False.
    """
    from phcpy.interface import store_standard_laurent_witness_set
    from phcpy.phcpy2c3 \
    import py2c_witset_standard_Laurent_ismember as membtest
    store_standard_laurent_witness_set(len(wsys), dim, wsys, gpts)
    nbc = len(point)
    nvr = len(wsys) - dim # test point should not have slack variables
    if verbose:
        print('calling standard_ismember with test point :')
        print(point)
    result = membtest(int(verbose), tasks, nvr, dim, nbc, evatol, memtol, point)
    (fail, onpolsys, inwitset) = result
    return (onpolsys, inwitset)

def dobldobl_laurent_ismember(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, verbose=True, tasks=0):
    r"""
    Applies the homotopy membership test for a *point* to belong to
    a witness set of dimension *dim*, given by an embedding polynomial
    Laurent system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the string *point*,
    which is the string representation of a solution in PHCpack format,
    with symbols of the variables before the values of the coordinates.
    By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in double double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    Returns a tuple of two booleans.  The first boolean is True if the
    point satisfies the equations, otherwise it is False.  The second
    boolean is True if the point belongs to the witness set,
    otherwise, the second boolean is False.
    """
    from phcpy.interface import store_dobldobl_laurent_witness_set
    from phcpy.phcpy2c3 \
    import py2c_witset_dobldobl_Laurent_ismember as membtest
    store_dobldobl_laurent_witness_set(len(wsys), dim, wsys, gpts)
    nbc = len(point)
    nvr = len(wsys) - dim # test point should not have slack variables
    if verbose:
        print('calling dobldobl_ismember with test point :')
        print(point)
    result = membtest(int(verbose), tasks, nvr, dim, nbc, evatol, memtol, point)
    (fail, onpolsys, inwitset) = result
    return (onpolsys, inwitset)

def quaddobl_laurent_ismember(wsys, gpts, dim, point, \
    evatol=1.0e-6, memtol=1.0e-6, verbose=True, tasks=0):
    r"""
    Applies the homotopy membership test for a *point* to belong to
    a witness set of dimension *dim*, given by an embedding polynomial
    Laurent system in *wsys*, with corresponding generic points in *gpts*.
    The coordinates of the test point are given in the string *point*,
    which is the string representation of a solution in PHCpack format,
    with symbols of the variables before the values of the coordinates.
    By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in quad double precision.
    The default values for the evaluation (*evatol*) and the membership
    (*memtol*) allow for singular values at the end points of the paths
    in the homotopy membership test.
    Returns a tuple of two booleans.  The first boolean is True if the
    point satisfies the equations, otherwise it is False.  The second
    boolean is True if the point belongs to the witness set,
    otherwise, the second boolean is False.
    """
    from phcpy.interface import store_quaddobl_laurent_witness_set
    from phcpy.phcpy2c3 \
    import py2c_witset_quaddobl_Laurent_ismember as membtest
    store_quaddobl_laurent_witness_set(len(wsys), dim, wsys, gpts)
    nbc = len(point)
    nvr = len(wsys) - dim # test point should not have slack variables
    if verbose:
        print('calling quaddobl_ismember with test point :')
        print(point)
    result = membtest(int(verbose), tasks, nvr, dim, nbc, evatol, memtol, point)
    (fail, onpolsys, inwitset) = result
    return (onpolsys, inwitset)

def standard_ismember_filter(wsys, gpts, dim, points, \
    rcotol=1.0e-6, evatol=1.0e-6, memtol=1.0e-6, verbose=True, tasks=0):
    r"""
    Given in *wsys* and *gpts* is a witness set of dimension *dim*,
    where *wsys* is an embedded polynomial system,
    and in *points* a list of strings.  The strings represent points
    as solutions in PHCpack format.  The homotopy membership test is
    applied to each point in the list *points*.  The list on return
    contains the points that do NOT belong to the witness set.
    Points that belong to the witness set are considered junk.
    By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in standard double precision.
    The parameter *rcotol* is used to bypass the homotopy membership test,
    for points with their estimated inverse condition number larger than
    *rcotol* will be considered isolated and not in the witness set.
    """
    from phcpy.solutions import diagnostics
    result = []
    for point in points:
        rco = diagnostics(point)[1]
        if rco > rcotol:
            (isgood, ismember) = (True, False)
        else:
            tst = standard_ismember(wsys, gpts, dim, point, \
                                    evatol, memtol, verbose, tasks)
            (isgood, ismember) = tst
        if isgood and not ismember:
            result.append(point)
    return result

def dobldobl_ismember_filter(wsys, gpts, dim, points, \
    rcotol=1.0e-6, evatol=1.0e-6, memtol=1.0e-6, verbose=True, tasks=0):
    r"""
    Given in *wsys* and *gpts* is a witness set of dimension *dim*,
    where *wsys* is an embedded polynomial system,
    and in *points* a list of strings.  The strings represent points
    as solutions in PHCpack format.  The homotopy membership test is
    applied to each point in the list *points*.  The list on return
    contains the points that do NOT belong to the witness set.
    Points that belong to the witness set are considered junk.
    By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in double double precision.
    The parameter *rcotol* is used to bypass the homotopy membership test,
    for points with their estimated inverse condition number larger than
    *rcotol* will be considered isolated and not in the witness set.
    """
    from phcpy.solutions import diagnostics
    result = []
    for point in points:
        rco = diagnostics(point)[1]
        if rco > rcotol:
            (isgood, ismember) = (True, False)
        else:
            tst = dobldobl_ismember(wsys, gpts, dim, point, \
                                    evatol, memtol, verbose, tasks)
            (isgood, ismember) = tst
        if isgood and not ismember:
            result.append(point)
    return result

def quaddobl_ismember_filter(wsys, gpts, dim, points, \
    rcotol=1.0e-6, evatol=1.0e-6, memtol=1.0e-6, verbose=True, tasks=0):
    r"""
    Given in *wsys* and *gpts* is a witness set of dimension *dim*,
    where *wsys* is an embedded polynomial system,
    and in *points* a list of strings.  The strings represent points
    as solutions in PHCpack format.  The homotopy membership test is
    applied to each point in the list *points*.  The list on return
    contains the points that do NOT belong to the witness set.
    Points that belong to the witness set are considered junk.
    By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in quad double precision.
    The parameter *rcotol* is used to bypass the homotopy membership test,
    for points with their estimated inverse condition number larger than
    *rcotol* will be considered isolated and not in the witness set.
    """
    from phcpy.solutions import diagnostics
    result = []
    for point in points:
        rco = diagnostics(point)[1]
        if rco > rcotol:
            (isgood, ismember) = (True, False)
        else:
            tst = quaddobl_ismember(wsys, gpts, dim, point, \
                                    evatol, memtol, verbose, tasks)
            (isgood, ismember) = tst
        if isgood and not ismember:
            result.append(point)
    return result

def standard_laurent_ismember_filter(wsys, gpts, dim, points, \
    rcotol=1.0e-6, evatol=1.0e-6, memtol=1.0e-6, verbose=True, tasks=0):
    r"""
    Given in *wsys* and *gpts* is a witness set of dimension *dim*,
    where *wsys* is an embedded Laurent polynomial system,
    and in *points* a list of strings.  The strings represent points
    as solutions in PHCpack format.  The homotopy membership test is
    applied to each point in the list *points*.  The list on return
    contains the points that do NOT belong to the witness set.
    Points that belong to the witness set are considered junk.
    By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in standard double precision.
    The parameter *rcotol* is used to bypass the homotopy membership test,
    for points with their estimated inverse condition number larger than
    *rcotol* will be considered isolated and not in the witness set.
    """
    from phcpy.solutions import diagnostics
    result = []
    for point in points:
        rco = diagnostics(point)[1]
        if rco > rcotol:
            (isgood, ismember) = (True, False)
        else:
            tst = standard_laurent_ismember(wsys, gpts, dim, point, \
                                            evatol, memtol, verbose, tasks)
            (isgood, ismember) = tst
        if isgood and not ismember:
            result.append(point)
    return result

def dobldobl_laurent_ismember_filter(wsys, gpts, dim, points, \
    rcotol=1.0e-6, evatol=1.0e-6, memtol=1.0e-6, verbose=True, tasks=0):
    r"""
    Given in *wsys* and *gpts* is a witness set of dimension *dim*,
    where *wsys* is an embedded Laurent polynomial system,
    and in *points* a list of strings.  The strings represent points
    as solutions in PHCpack format.  The homotopy membership test is
    applied to each point in the list *points*.  The list on return
    contains the points that do NOT belong to the witness set.
    Points that belong to the witness set are considered junk.
    By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in double double precision.
    The parameter *rcotol* is used to bypass the homotopy membership test,
    for points with their estimated inverse condition number larger than
    *rcotol* will be considered isolated and not in the witness set.
    """
    from phcpy.solutions import diagnostics
    result = []
    for point in points:
        rco = diagnostics(points)[1]
        if rco > rcotol:
            (isgood, ismember) = (True, False)
        else:
            tst = dobldobl_laurent_ismember(wsys, gpts, dim, point, \
                                            evatol, memtol, verbose, tasks)
            (isgood, ismember) = tst
        if isgood and not ismember:
            result.append(point)
    return result

def quaddobl_laurent_ismember_filter(wsys, gpts, dim, points, \
    rcotol=1.0e-6, evatol=1.0e-6, memtol=1.0e-6, verbose=True, tasks=0):
    r"""
    Given in *wsys* and *gpts* is a witness set of dimension *dim*,
    where *wsys* is an embedded Laurent polynomial system,
    and in *points* a list of strings.  The strings represent points
    as solutions in PHCpack format.  The homotopy membership test is
    applied to each point in the list *points*.  The list on return
    contains the points that do NOT belong to the witness set.
    Points that belong to the witness set are considered junk.
    By default, *verbose* is True.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    Calculations happen in quad double precision.
    The parameter *rcotol* is used to bypass the homotopy membership test,
    for points with their estimated inverse condition number larger than
    *rcotol* will be considered isolated and not in the witness set.
    """
    from phcpy.solutions import diagnostics
    result = []
    for point in points:
        rco = diagnostics(point)[1]
        if rco > rcotol:
            (isgood, ismember) = (True, False)
        else:
            tst = quaddobl_laurent_ismember(wsys, gpts, dim, point, \
                                            evatol, memtol, verbose, tasks)
            (isgood, ismember) = tst
        if isgood and not ismember:
            result.append(point)
    return result

def ismember_filter(wsys, gpts, dim, points, \
    rcotol=1.0e-6, evatol=1.0e-6, memtol=1.0e-6, \
    verbose=True, precision='d', tasks=0):
    r"""
    Filters *points* so the list on return contains only those points
    which do not belong to the witness set of dimension *dim*,
    given by an embedded polynomial system in *wsys*,
    with corresponding generic points in *gpts*.
    The list *points* is a list of strings.  Each string is the
    symbolic string representation of a solution.
    By default, *verbose* is True, and the precision is double 'd'.
    Other levels of precision are double double precision 'dd' 
    and quad double precision 'qd'.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    The parameter *rcotol* is used to bypass the homotopy membership test,
    for points with their estimated inverse condition number larger than
    *rcotol* will be considered isolated and not in the witness set.
    The homotopy membership test has two tolerances: *evatol* and *memtol*.
    The *evatol* is the tolerance on the residual
    of the evaluation of the polynomial equations at the test point.
    If the residual of the evaluation is not less than *evatol*,
    then the point is not a member.  Otherwise, the homotopy
    membership test is called and the *memtol* is used to compare
    the coordinates of the point with the newly computed generic points.
    If there is a match between the coordinates within the given
    tolerance *memtol*, then the points is a member and filtered out.
    """
    if(precision == 'd'):
        return standard_ismember_filter\
                   (wsys, gpts, dim, points, rcotol, evatol, memtol, \
                    verbose, tasks)
    elif(precision == 'dd'):
        return dobldobl_ismember_filter\
                   (wsys, gpts, dim, points, rcotol, evatol, memtol, \
                    verbose, tasks)
    elif(precision == 'qd'):
        return quaddobl_ismember_filter\
                   (wsys, gpts, dim, points, rcotol, evatol, memtol, \
                    verbose, tasks)
    else:
        print('wrong argument for precision')
        return points

def laurent_ismember_filter(wsys, gpts, dim, points, \
    rcotol=1.0e-6, evatol=1.0e-6, memtol=1.0e-6, \
    verbose=True, precision='d', tasks=0):
    r"""
    Filters *points* so the list on return contains only those points
    which do not belong to the witness set of dimension *dim*,
    given by an embedded Laurent polynomial system in *wsys*,
    with corresponding generic points in *gpts*.
    The list *points* is a list of strings.  Each string is the
    symbolic string representation of a solution.
    By default, *verbose* is True, and the precision is double 'd'.
    Other levels of precision are double double precision 'dd' 
    and quad double precision 'qd'.
    The number of threads is given in *tasks*.  If *tasks* is zero,
    then no multithreading is applied in the homotopy membership test.
    The parameter *rcotol* is used to bypass the homotopy membership test,
    for points with their estimated inverse condition number larger than
    *rcotol* will be considered isolated and not in the witness set.
    The homotopy membership test has two tolerances: *evatol* and *memtol*.
    The *evatol* is the tolerance on the residual
    of the evaluation of the polynomial equations at the test point.
    If the residual of the evaluation is not less than *evatol*,
    then the point is not a member.  Otherwise, the homotopy
    membership test is called and the *memtol* is used to compare
    the coordinates of the point with the newly computed generic points.
    If there is a match between the coordinates within the given
    tolerance *memtol*, then the point is a member and filtered out.
    """
    if(precision == 'd'):
        return standard_laurent_ismember_filter\
                   (wsys, gpts, dim, points, rcotol, evatol, memtol,
                    verbose, tasks)
    elif(precision == 'dd'):
        return dobldobl_laurent_ismember_filter\
                   (wsys, gpts, dim, points, rcotol, evatol, memtol,
                    verbose, tasks)
    elif(precision == 'qd'):
        return quaddobl_laurent_ismember_filter\
                   (wsys, gpts, dim, points, rcotol, evatol, memtol,
                    verbose, tasks)
    else:
        print('wrong argument for precision')
        return points

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
    verbose=True, precision='d', tasks=0):
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

def test_twisted_ismember(prc='d', laurent=True):
    """
    To test the membertest wrapper, we take the twisted cubic again.
    The test point is given as a solution in PHCpack format.
    """
    from phcpy.solver import solve
    from phcpy.solutions import make_solution, strsol2dict
    twisted = ['x^2 - y;', 'x^3 - z;']
    twiste1 = embed(3, 1, twisted)
    twtsols = solve(twiste1, precision=prc)
    for sol in twtsols:
        print(sol)
    input('hit enter to continue')
    # print(is_member(twiste1, twtsols, 1, twtsols[0], precision=prc))
    genpt = twtsols[0]
    dicpt = strsol2dict(genpt)
    coord = [dicpt['x'], dicpt['y'], dicpt['z']]
    tstpt = make_solution(['x', 'y', 'z'], coord)
    if laurent:
        if prc == 'd':
            print(standard_laurent_ismember\
                      (twiste1, twtsols, 1, tstpt, tasks=3))
        elif prc == 'dd':
            print(dobldobl_laurent_ismember\
                      (twiste1, twtsols, 1, tstpt, tasks=3))
        elif prc == 'qd':
            print(quaddobl_laurent_ismember\
                      (twiste1, twtsols, 1, tstpt, tasks=3))
        else:
            print( 'wrong level of precision')
    else:
        if prc == 'd':
            print(standard_ismember(twiste1, twtsols, 1, tstpt, tasks=3))
        elif prc == 'dd':
            print(dobldobl_ismember(twiste1, twtsols, 1, tstpt, tasks=3))
        elif prc == 'qd':
            print(quaddobl_ismember(twiste1, twtsols, 1, tstpt, tasks=3))
        else:
            print( 'wrong level of precision')
    input('hit enter to continue')
    outsol = make_solution(['x', 'y', 'z'], [1, 2, 3])
    # print(is_member(twiste1, twtsols, 1, outsol, precision=prc))
    if laurent:
        if prc == 'd':
            print(standard_laurent_ismember\
                      (twiste1, twtsols, 1, outsol, tasks=3))
        elif prc == 'dd':
            print(dobldobl_laurent_ismember\
                      (twiste1, twtsols, 1, outsol, tasks=3))
        elif prc == 'qd':
            print(quaddobl_laurent_ismember\
                      (twiste1, twtsols, 1, outsol, tasks=3))
        else:
            print('wrong level of precision')
    else:
        if prc == 'd':
            print(standard_ismember(twiste1, twtsols, 1, outsol, tasks=3))
        elif prc == 'dd':
            print(dobldobl_ismember(twiste1, twtsols, 1, outsol, tasks=3))
        elif prc == 'qd':
            print(quaddobl_ismember(twiste1, twtsols, 1, outsol, tasks=3))
        else:
            print('wrong level of precision')

def test():
    """
    Runs a test on algebraic sets.
    """
    from phcpy.phcpy2c3 import py2c_set_seed
    py2c_set_seed(234798272)
    test_twisted_ismember(prc='qd', laurent=False)

if __name__ == "__main__":
    test()
