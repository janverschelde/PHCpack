"""
This module provides the data interface to PHCpack.
The user gives as input string representations of polynomials or solutions
to the interface functions which store the data.
"""

def store_standard_tableau(poltab, verbose=False):
    r"""
    Stores the polynomial system given in the list of lists *poltab* in
    the container for systems with coefficients in standard double precision.
    Every polynomial in the system is represented by a list of tuples.
    A monomial is represented by a 2-tuple:

    1. the coefficient of the monomial is a complex number,

    2. the exponent are a tuple of natural numbers.

    For example, the system x^2 - y = 0, x^3 - z = 0 is represented as
    [[((1+0j), (2, 0, 0)), ((-1+0j), (0, 1, 0))], \
    [((1+0j), (3, 0, 0)), ((-1+0j), (0, 0, 1))]].
    """
    from phcpy.phcpy2c3 import py2c_tabform_store_standard_tableau as store
    neq = len(poltab)                               # number of equations
    if verbose:
        print('number of equations :', neq)
    if neq > 0:
        nvr = len(poltab[0][0][1])                  # number of variables
        if verbose:
            print('number of variables :', nvr)
        nbt = [len(pol) for pol in poltab]          # number of terms
        if verbose:
            print('number of terms :', nbt)
        cff = [[c for (c, e) in p] for p in poltab] # coefficients
        if verbose:
            print('the coefficients on input :', cff)
        flatcff = sum(cff, [])                      # flatten list of lists
        frimcff = sum([[x.real, x.imag] for x in flatcff], [])
        if verbose:
            print('flat list of coefficients :', frimcff)
        xps = [[e for (c, e) in p] for p in poltab] # exponents
        if verbose:
            print('the exponents on input :', xps)
        txps = sum(xps, [])                         # list of tuples
        lxps = [list(x) for x in txps]              # list of lists
        flatxps = sum(lxps, [])                     # flatten list of lists
        if verbose:
            print('flat list of exponents :', flatxps)
        strnbt = str(nbt)
        strcff = str(frimcff)
        strxps = str(flatxps)
        fail = store(neq, nvr, len(strnbt), strnbt, \
                     len(strcff), strcff, len(strxps), strxps, int(verbose))

def load_standard_tableau(verbose=False):
    """
    Returns the tableau form of the system stored in the container for
    double precision coefficients.
    """
    from phcpy.phcpy2c3 import py2c_tabform_load_standard_tableau as load
    neq, nvr, nbterms, coefficients, exponents = load(int(verbose))
    nbt = eval(nbterms)
    cff = eval(coefficients)
    xps = eval(exponents)
    xpsidx = 0
    cffidx = 0
    tableau = []
    for i in range(neq):
        pol = []
        for j in range(nbt[i]):
            trm = []
            for k in range(nvr):
                trm.append(xps[xpsidx])
                xpsidx = xpsidx + 1
            cfftrm = complex(cff[cffidx], cff[cffidx+1])
            cffidx = cffidx + 2
            pol.append((cfftrm ,tuple(trm)))
        tableau.append(pol)
    return (neq, nvr, tableau)

def store_standard_system(polsys, **nbvar):
    r"""
    Stores the polynomials represented by the list of strings in *polsys* into
    the container for systems with coefficients in standard double precision.
    The number of variables is an optional argument given in *nbvar*.
    If *nbvar* is omitted, then the system is assumed to be square.
    Otherwise, suppose the number of variables equals 2 and pols is the list
    of polynomials, then the call **store_standard_system(pols, nbvar=2)**
    will store the polynomials in pols in the standard systems container.
    """
    from phcpy.phcpy2c3 import py2c_syscon_clear_standard_system
    from phcpy.phcpy2c3 \
        import py2c_syscon_initialize_number_of_standard_polynomials
    from phcpy.phcpy2c3 import py2c_syscon_store_standard_polynomial
    py2c_syscon_clear_standard_system()
    dim = len(polsys)
    fail = 0
    py2c_syscon_initialize_number_of_standard_polynomials(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        if(len(nbvar) == 0):
            fail = py2c_syscon_store_standard_polynomial(nchar, dim, cnt+1, pol)
        else:
            nvr = list(nbvar.values())[0]
            fail = py2c_syscon_store_standard_polynomial(nchar, nvr, cnt+1, pol)
        if(fail != 0):
            break
    return fail

def store_dobldobl_system(polsys, **nbvar):
    r"""
    Stores the polynomials represented by the list of strings in *polsys*
    into the systems container for double double arithmetic.
    The number of variables is an optional argument given in *nbvar*.
    If *nbvar* is omitted, then the system is assumed to be square.
    Otherwise, suppose the number of variables equals 2 and pols is the list
    of polynomials, then the call **store_dobldobl_system(pols, nbvar=2)**
    will store the polynomials in pols in the dobldobl systems container.
    """
    from phcpy.phcpy2c3 import py2c_syscon_clear_dobldobl_system
    from phcpy.phcpy2c3 \
        import py2c_syscon_initialize_number_of_dobldobl_polynomials
    from phcpy.phcpy2c3 import py2c_syscon_store_dobldobl_polynomial
    py2c_syscon_clear_dobldobl_system()
    dim = len(polsys)
    fail = 0
    py2c_syscon_initialize_number_of_dobldobl_polynomials(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        if(len(nbvar) == 0):
            fail = py2c_syscon_store_dobldobl_polynomial(nchar, dim, cnt+1, pol)
        else:
            nvr = list(nbvar.values())[0]
            fail = py2c_syscon_store_dobldobl_polynomial(nchar, nvr, cnt+1, pol)
        if(fail != 0):
            break
    return fail

def store_quaddobl_system(polsys, **nbvar):
    r"""
    Stores the polynomials represented by the list of strings in *polsys*
    into the systems container for quad double arithmetic.
    The number of variables is an optional argument given in *nbvar*.
    If *nbvar* is omitted, then the system is assumed to be square.
    Otherwise, suppose the number of variables equals 2 and pols is the list
    of polynomials, then the call **store_quaddobl_system(pols, nbvar=2)**
    will store the polynomials in pols in the quaddobl systems container.
    """
    from phcpy.phcpy2c3 import py2c_syscon_clear_quaddobl_system
    from phcpy.phcpy2c3 \
        import py2c_syscon_initialize_number_of_quaddobl_polynomials
    from phcpy.phcpy2c3 import py2c_syscon_store_quaddobl_polynomial
    py2c_syscon_clear_quaddobl_system()
    dim = len(polsys)
    fail = 0
    py2c_syscon_initialize_number_of_quaddobl_polynomials(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        if(len(nbvar) == 0):
            fail = py2c_syscon_store_quaddobl_polynomial(nchar, dim, cnt+1, pol)
        else:
            nvr = list(nbvar.values())[0]
            fail = py2c_syscon_store_quaddobl_polynomial(nchar, nvr, cnt+1, pol)
        if(fail != 0):
            break
    return fail

def store_multprec_system(polsys, decimals, **nbvar):
    r"""
    Stores the polynomials represented by the list of strings in *polsys*
    into the systems container for multiprecision arithmetic.
    The parameter *decimals* equals the number of decimal places
    in the working precision for the parsing of the strings in *polsys*.
    The number of variables is an optional argument given in *nbvar*.
    If *nbvar* is omitted, then the system is assumed to be square.
    Otherwise, suppose the number of variables equals 2 and pols is the list of
    polynomials, then the call **store_multprec_system(pols, nbvar=2)**
    will store the polynomials in pols in the multiprecision systems container.
    """
    from phcpy.phcpy2c3 import py2c_syscon_clear_multprec_system
    from phcpy.phcpy2c3 \
        import py2c_syscon_initialize_number_of_multprec_polynomials
    from phcpy.phcpy2c3 import py2c_syscon_store_multprec_polynomial
    py2c_syscon_clear_multprec_system()
    dim = len(polsys)
    fail = 0
    py2c_syscon_initialize_number_of_multprec_polynomials(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        if(len(nbvar) == 0):
            fail = py2c_syscon_store_multprec_polynomial\
                       (nchar, dim, cnt+1, decimals, pol)
        else:
            nvr = list(nbvar.values())[0]
            fail = py2c_syscon_store_multprec_polynomial\
                       (nchar, nvr, cnt+1, decimals, pol)
        if(fail != 0):
            break
    return fail

def load_standard_system():
    """
    Returns the polynomials stored in the system container
    for standard double precision arithmetic.
    """
    from phcpy.phcpy2c3 import py2c_syscon_number_of_standard_polynomials
    from phcpy.phcpy2c3 import py2c_syscon_load_standard_polynomial
    dim = py2c_syscon_number_of_standard_polynomials()
    result = []
    for ind in range(1, dim+1):
        result.append(py2c_syscon_load_standard_polynomial(ind))
    return result

def load_dobldobl_system():
    """
    Returns the polynomials stored in the system container
    with double double complex coefficients.
    """
    from phcpy.phcpy2c3 import py2c_syscon_number_of_dobldobl_polynomials
    from phcpy.phcpy2c3 import py2c_syscon_load_dobldobl_polynomial
    dim = py2c_syscon_number_of_dobldobl_polynomials()
    result = []
    for ind in range(1, dim+1):
        result.append(py2c_syscon_load_dobldobl_polynomial(ind))
    return result

def load_quaddobl_system():
    """
    Returns the polynomials stored in the system container
    with quad double complex coefficients.
    """
    from phcpy.phcpy2c3 import py2c_syscon_number_of_quaddobl_polynomials
    from phcpy.phcpy2c3 import py2c_syscon_load_quaddobl_polynomial
    dim = py2c_syscon_number_of_quaddobl_polynomials()
    result = []
    for ind in range(1, dim+1):
        result.append(py2c_syscon_load_quaddobl_polynomial(ind))
    return result

def load_multprec_system():
    """
    Returns the polynomials stored in the system container
    with arbitrary multiprecision complex coefficients.
    """
    from phcpy.phcpy2c3 import py2c_syscon_number_of_multprec_polynomials
    from phcpy.phcpy2c3 import py2c_syscon_load_multprec_polynomial
    dim = py2c_syscon_number_of_multprec_polynomials()
    result = []
    for ind in range(1, dim+1):
        result.append(py2c_syscon_load_multprec_polynomial(ind))
    return result

def store_standard_laurent_system(polsys, **nbvar):
    r"""
    Stores the Laurent polynomials represented by the list of
    strings in *polsys* into the container for systems
    with coefficients in standard double precision.
    If *nbvar* is omitted, then the system is assumed to be square.
    Otherwise, suppose the number of variables equals 2 and pols is the list
    of polynomials, then **store_standard_laurent_system(pols, nbvar=2)**
    stores the polynomials in pols in the standard Laurent systems container.
    """
    from phcpy.phcpy2c3 import py2c_syscon_clear_standard_Laurent_system
    from phcpy.phcpy2c3 \
        import py2c_syscon_initialize_number_of_standard_Laurentials
    from phcpy.phcpy2c3 import py2c_syscon_store_standard_Laurential
    py2c_syscon_clear_standard_Laurent_system()
    dim = len(polsys)
    fail = 0
    py2c_syscon_initialize_number_of_standard_Laurentials(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        if len(nbvar) == 0:
            fail = py2c_syscon_store_standard_Laurential(nchar, dim, cnt+1, pol)
        else:
            nvr = list(nbvar.values())[0]
            fail = py2c_syscon_store_standard_Laurential(nchar, nvr, cnt+1, pol)
        if(fail != 0):
            break
    return fail

def store_dobldobl_laurent_system(polsys, **nbvar):
    r"""
    Stores the Laurent polynomials represented by the list of
    strings in *polsys* into the container for systems
    with coefficients in double double precision.
    If *nbvar* is omitted, then the system is assumed to be square.
    Otherwise, suppose the number of variables equals 2 and pols is the list
    of polynomials, then **store_dobldobl_laurent_system(pols, nbvar=2)**
    stores the polynomials in pols in the dobldobl Laurent systems container.
    """
    from phcpy.phcpy2c3 import py2c_syscon_clear_dobldobl_Laurent_system
    from phcpy.phcpy2c3 \
        import py2c_syscon_initialize_number_of_dobldobl_Laurentials
    from phcpy.phcpy2c3 import py2c_syscon_store_dobldobl_Laurential
    py2c_syscon_clear_dobldobl_Laurent_system()
    dim = len(polsys)
    fail = 0
    py2c_syscon_initialize_number_of_dobldobl_Laurentials(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        if len(nbvar) == 0:
            fail = py2c_syscon_store_dobldobl_Laurential(nchar, dim, cnt+1, pol)
        else:
            nvr = list(nbvar.values())[0]
            fail = py2c_syscon_store_dobldobl_Laurential(nchar, nvr, cnt+1, pol)
        if(fail != 0):
            break
    return fail

def store_quaddobl_laurent_system(polsys, **nbvar):
    r"""
    Stores the Laurent polynomials represented by the list
    of strings in *polsys* into the container for systems
    with coefficients in quad double precision.
    If *nbvar* is omitted, then the system is assumed to be square.
    Otherwise, suppose the number of variables equals 2 and pols is the list
    of polynomials, then **store_quaddobl_laurent_system(pols, nbvar=2)**
    stores the polynomials in pols in the quaddobl Laurent systems container.
    """
    from phcpy.phcpy2c3 import py2c_syscon_clear_quaddobl_Laurent_system
    from phcpy.phcpy2c3 \
        import py2c_syscon_initialize_number_of_quaddobl_Laurentials
    from phcpy.phcpy2c3 import py2c_syscon_store_quaddobl_Laurential
    py2c_syscon_clear_quaddobl_Laurent_system()
    dim = len(polsys)
    fail = 0
    py2c_syscon_initialize_number_of_quaddobl_Laurentials(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        if len(nbvar) == 0:
            fail = py2c_syscon_store_quaddobl_Laurential(nchar, dim, cnt+1, pol)
        else:
            nvr = list(nbvar.values())[0]
            fail = py2c_syscon_store_quaddobl_Laurential(nchar, nvr, cnt+1, pol)
        if(fail != 0):
            break
    return fail

def store_multprec_laurent_system(polsys, decimals, **nbvar):
    r"""
    Stores the Laurent polynomials represented by the list
    of strings in *polsys* into the container for systems
    with coefficients in multiprecision.
    The parameter *decimals* equals the number of decimal places
    in the working precision for the parsing of the strings in *polsys*.
    If *nbvar* is omitted, then the system is assumed to be square.
    Otherwise, suppose the number of variables equals 2 and pols is the list
    of polynomials, then **store_multprec_laurent_system(pols, nbvar=2)**
    stores the polynomials in pols in the multprec Laurent systems container.
    """
    from phcpy.phcpy2c3 import py2c_syscon_clear_multprec_Laurent_system
    from phcpy.phcpy2c3 \
        import py2c_syscon_initialize_number_of_multprec_Laurentials
    from phcpy.phcpy2c3 import py2c_syscon_store_multprec_Laurential
    py2c_syscon_clear_multprec_Laurent_system()
    dim = len(polsys)
    fail = 0
    py2c_syscon_initialize_number_of_multprec_Laurentials(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        if len(nbvar) == 0:
            fail = py2c_syscon_store_multprec_Laurential\
                       (nchar, dim, cnt+1, decimals, pol)
        else:
            nvr = list(nbvar.values())[0]
            fail = py2c_syscon_store_multprec_Laurential\
                       (nchar, nvr, cnt+1, decimals, pol)
        if(fail != 0):
            break
    return fail

def load_standard_laurent_system():
    """
    Returns the Laurent polynomials stored in the system container
    for standard double precision arithmetic.
    """
    from phcpy.phcpy2c3 import py2c_syscon_number_of_standard_Laurentials
    from phcpy.phcpy2c3 import py2c_syscon_load_standard_Laurential
    dim = py2c_syscon_number_of_standard_Laurentials()
    result = []
    for ind in range(1, dim+1):
        result.append(py2c_syscon_load_standard_Laurential(ind))
    return result

def load_dobldobl_laurent_system():
    """
    Returns the Laurent polynomials stored in the system container
    with double double complex coefficients.
    """
    from phcpy.phcpy2c3 import py2c_syscon_number_of_dobldobl_Laurentials
    from phcpy.phcpy2c3 import py2c_syscon_load_dobldobl_Laurential
    dim = py2c_syscon_number_of_dobldobl_Laurentials()
    result = []
    for ind in range(1, dim+1):
        result.append(py2c_syscon_load_dobldobl_Laurential(ind))
    return result

def load_quaddobl_laurent_system():
    """
    Returns the Laurent polynomials stored in the system container
    with quad double complex coefficients.
    """
    from phcpy.phcpy2c3 import py2c_syscon_number_of_quaddobl_Laurentials
    from phcpy.phcpy2c3 import py2c_syscon_load_quaddobl_Laurential
    dim = py2c_syscon_number_of_quaddobl_Laurentials()
    result = []
    for ind in range(1, dim+1):
        result.append(py2c_syscon_load_quaddobl_Laurential(ind))
    return result

def load_multprec_laurent_system():
    """
    Returns the Laurent polynomials stored in the system container
    with multiprecision complex coefficients.
    """
    from phcpy.phcpy2c3 import py2c_syscon_number_of_multprec_Laurentials
    from phcpy.phcpy2c3 import py2c_syscon_load_multprec_Laurential
    dim = py2c_syscon_number_of_multprec_Laurentials()
    result = []
    for ind in range(1, dim+1):
        result.append(py2c_syscon_load_multprec_Laurential(ind))
    return result

def read_standard_system(filename):
    r"""
    Opens the *filename* for reading a polynomial system
    with coefficients in standard double precision.
    Returns the list of polynomials in the system
    or *None* if something went wrong.
    """
    from phcpy.phcpy2c3 import py2c_syscon_clear_symbol_table
    from phcpy.phcpy2c3 import py2c_read_standard_target_system_from_file
    from phcpy.phcpy2c3 import py2c_copy_standard_target_system_to_container
    py2c_syscon_clear_symbol_table()
    lnf = len(filename)
    fail = py2c_read_standard_target_system_from_file(lnf,filename)
    if(fail != 0):
        return None
    else:
        py2c_copy_standard_target_system_to_container()
        return load_standard_system()

def read_dobldobl_system(filename):
    r"""
    Opens the *filename* for reading a polynomial system
    with coefficients in double double precision.
    Returns the list of polynomials in the system
    or *None* if something went wrong.
    """
    from phcpy.phcpy2c3 import py2c_syscon_clear_symbol_table
    from phcpy.phcpy2c3 import py2c_read_dobldobl_target_system_from_file
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_target_system_to_container
    py2c_syscon_clear_symbol_table()
    lnf = len(filename)
    fail = py2c_read_dobldobl_target_system_from_file(lnf,filename)
    if(fail != 0):
        return []
    else:
        py2c_copy_dobldobl_target_system_to_container()
        return load_dobldobl_system()

def read_quaddobl_system(filename):
    r"""
    Opens the *filename* for reading a polynomial system
    with coefficients in quad double precision.
    Returns the list of polynomials in the system
    or *None* if something went wrong.
    """
    from phcpy.phcpy2c3 import py2c_syscon_clear_symbol_table
    from phcpy.phcpy2c3 import py2c_read_quaddobl_target_system_from_file
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_target_system_to_container
    py2c_syscon_clear_symbol_table()
    lnf = len(filename)
    fail = py2c_read_quaddobl_target_system_from_file(lnf,filename)
    if(fail != 0):
        return []
    else:
        py2c_copy_quaddobl_target_system_to_container()
        return load_quaddobl_system()

def read_standard_system_and_solutions(filename):
    r"""
    Opens the *filename* for reading a polynomial system
    with coefficients in standard double precision,
    and its corresponding list of solutions.
    Returns *None* if the reading went wrong, or otherwise
    returns a tuple with first the list of polynomials
    and second the list of solutions.
    """
    from phcpy.phcpy2c3 import py2c_syscon_clear_symbol_table
    from phcpy.phcpy2c3 import py2c_read_standard_start_system_from_file
    from phcpy.phcpy2c3 import py2c_copy_start_system_to_container
    from phcpy.phcpy2c3 import py2c_copy_start_solutions_to_container
    py2c_syscon_clear_symbol_table()
    lnf = len(filename)
    fail = py2c_read_standard_start_system_from_file(lnf,filename)
    if(fail != 0):
        return None
    else:
        py2c_copy_start_system_to_container()
        py2c_copy_start_solutions_to_container()
        return (load_standard_system(), load_standard_solutions())

def read_dobldobl_system_and_solutions(filename):
    r"""
    Opens the *filename* for reading a polynomial system
    with coefficients in double double precision,
    and its corresponding list of solutions.
    Returns *None* if the reading went wrong, or otherwise
    returns a tuple with first the list of polynomials
    and second the list of solutions.
    """
    from phcpy.phcpy2c3 import py2c_syscon_clear_symbol_table
    from phcpy.phcpy2c3 import py2c_read_dobldobl_start_system_from_file
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_start_system_to_container
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_start_solutions_to_container
    py2c_syscon_clear_symbol_table()
    lnf = len(filename)
    fail = py2c_read_dobldobl_start_system_from_file(lnf,filename)
    if(fail != 0):
        return None
    else:
        py2c_copy_dobldobl_start_system_to_container()
        py2c_copy_dobldobl_start_solutions_to_container()
        return (load_dobldobl_system(), load_dobldobl_solutions())

def read_quaddobl_system_and_solutions(filename):
    r"""
    Opens the *filename* for reading a polynomial system
    with coefficients in quad double precision,
    and its corresponding list of solutions.
    Returns *None* if the reading went wrong, or otherwise
    returns a tuple with first the list of polynomials
    and second the list of solutions.
    """
    from phcpy.phcpy2c3 import py2c_syscon_clear_symbol_table
    from phcpy.phcpy2c3 import py2c_read_quaddobl_start_system_from_file
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_start_system_to_container
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_start_solutions_to_container
    py2c_syscon_clear_symbol_table()
    lnf = len(filename)
    fail = py2c_read_quaddobl_start_system_from_file(lnf,filename)
    if(fail != 0):
        return None
    else:
        py2c_copy_quaddobl_start_system_to_container()
        py2c_copy_quaddobl_start_solutions_to_container()
        return (load_quaddobl_system(), load_quaddobl_solutions())

def store_standard_solutions(nvar, sols):
    r"""
    Stores the solutions in the list *sols*, represented as strings
    in PHCpack format into the container for solutions
    with standard double precision.
    The number *nvar* equals the number of variables.
    """
    from phcpy.phcpy2c3 import py2c_solcon_clear_standard_solutions
    from phcpy.phcpy2c3 import py2c_solcon_append_standard_solution_string
    py2c_solcon_clear_standard_solutions()
    fail = 0
    for ind in range(0, len(sols)):
        fail = py2c_solcon_append_standard_solution_string\
                   (nvar, len(sols[ind]), sols[ind])
        if(fail != 0):
            # break
            print('Solution at position', ind, 'is not appended.')
    return fail

def store_dobldobl_solutions(nvar, sols):
    r"""
    Stores the solutions in the list *sols*, represented as strings
    in PHCpack format into the solution container for processing
    with complex double double arithmetic.
    The number *nvar* equals the number of variables.
    """
    from phcpy.phcpy2c3 import py2c_solcon_clear_dobldobl_solutions
    from phcpy.phcpy2c3 import py2c_solcon_append_dobldobl_solution_string
    py2c_solcon_clear_dobldobl_solutions()
    fail = 0
    for ind in range(0, len(sols)):
        fail = py2c_solcon_append_dobldobl_solution_string\
                   (nvar, len(sols[ind]), sols[ind])
        if(fail != 0):
            # break
            print('Solution at position', ind, 'is not appended.')
    return fail

def store_quaddobl_solutions(nvar, sols):
    r"""
    Stores the solutions in the list *sols*, represented as strings
    in PHCpack format into the solution container for processing
    with complex quad double arithmetic.
    The number *nvar* equals the number of variables.
    """
    from phcpy.phcpy2c3 import py2c_solcon_clear_quaddobl_solutions
    from phcpy.phcpy2c3 import py2c_solcon_append_quaddobl_solution_string
    py2c_solcon_clear_quaddobl_solutions()
    fail = 0
    for ind in range(0, len(sols)):
        fail = py2c_solcon_append_quaddobl_solution_string\
                   (nvar, len(sols[ind]), sols[ind])
        if(fail != 0):
            # break
            print('Solution at position', ind, 'is not appended.')
    return fail

def store_multprec_solutions(nvar, sols):
    r"""
    Stores the solutions in the list *sols*, represented as strings
    in PHCpack format into the solution container for processing
    with complex multiprecision arithmetic.
    The number *nvar* equals the number of variables.
    """
    from phcpy.phcpy2c3 import py2c_solcon_clear_multprec_solutions
    from phcpy.phcpy2c3 import py2c_solcon_append_multprec_solution_string
    py2c_solcon_clear_multprec_solutions()
    fail = 0
    for ind in range(0, len(sols)):
        fail = py2c_solcon_append_multprec_solution_string\
                   (nvar, len(sols[ind]), sols[ind])
        if(fail != 0):
            # break
            print('Solution at position', ind, 'is not appended.')
    return fail

def load_standard_solutions():
    """
    Returns the list of solutions stored in the container
    for solutions with standard double precision.
    """
    from phcpy.phcpy2c3 import py2c_solcon_retrieve_next_standard_initialize
    from phcpy.phcpy2c3 import py2c_solcon_move_current_standard_to_next
    from phcpy.phcpy2c3 \
        import py2c_solcon_length_current_standard_solution_string
    from phcpy.phcpy2c3 \
        import py2c_solcon_write_current_standard_solution_string
    result = []
    py2c_solcon_retrieve_next_standard_initialize()
    while True:
        lns = py2c_solcon_length_current_standard_solution_string()
        if(lns == 0):
            break
        sol = py2c_solcon_write_current_standard_solution_string(lns)
        result.append(sol)
        ind = py2c_solcon_move_current_standard_to_next()
        if(ind == 0):
            break
    return result

def load_dobldobl_solutions():
    """
    Returns the list of solutions stored in the container
    for complex double double solutions.
    """
    from phcpy.phcpy2c3 import py2c_solcon_retrieve_next_dobldobl_initialize
    from phcpy.phcpy2c3 import py2c_solcon_move_current_dobldobl_to_next
    from phcpy.phcpy2c3 \
        import py2c_solcon_length_current_dobldobl_solution_string
    from phcpy.phcpy2c3 \
        import py2c_solcon_write_current_dobldobl_solution_string
    result = []
    py2c_solcon_retrieve_next_dobldobl_initialize()
    while True:
        lns = py2c_solcon_length_current_dobldobl_solution_string()
        if(lns == 0):
            break
        sol = py2c_solcon_write_current_dobldobl_solution_string(lns)
        result.append(sol)
        ind = py2c_solcon_move_current_dobldobl_to_next()
        if(ind == 0):
            break
    return result

def load_quaddobl_solutions():
    """
    Returns the list of solutions stored in the container
    for complex quad double solutions.
    """
    from phcpy.phcpy2c3 import py2c_solcon_retrieve_next_quaddobl_initialize
    from phcpy.phcpy2c3 import py2c_solcon_move_current_quaddobl_to_next
    from phcpy.phcpy2c3 \
        import py2c_solcon_length_current_quaddobl_solution_string
    from phcpy.phcpy2c3 \
        import py2c_solcon_write_current_quaddobl_solution_string
    result = []
    py2c_solcon_retrieve_next_quaddobl_initialize()
    while True:
        lns = py2c_solcon_length_current_quaddobl_solution_string()
        if(lns == 0):
            break
        sol = py2c_solcon_write_current_quaddobl_solution_string(lns)
        result.append(sol)
        ind = py2c_solcon_move_current_quaddobl_to_next()
        if(ind == 0):
            break
    return result

def load_multprec_solutions():
    """
    Returns the list of solutions stored in the container
    for complex multiprecision solutions.
    """
    from phcpy.phcpy2c3 import py2c_solcon_retrieve_next_multprec_initialize
    from phcpy.phcpy2c3 import py2c_solcon_move_current_multprec_to_next
    from phcpy.phcpy2c3 \
        import py2c_solcon_length_current_multprec_solution_string
    from phcpy.phcpy2c3 \
        import py2c_solcon_write_current_multprec_solution_string
    result = []
    py2c_solcon_retrieve_next_multprec_initialize()
    while True:
        lns = py2c_solcon_length_current_multprec_solution_string()
        if(lns == 0):
            break
        sol = py2c_solcon_write_current_multprec_solution_string(lns)
        result.append(sol)
        ind = py2c_solcon_move_current_multprec_to_next()
        if(ind == 0):
            break
    return result

def store_standard_witness_set(nbvar, dim, pols, sols):
    r"""
    Given in *nbar* is the total number of variables in the list of
    polynomials in *pols* and its list of solutions in *sols*.
    The coefficients in the polynomials and the coordinates of the
    solutions will be parsed and stored in standard double precision.
    The parameter *dim* equals the number of slack variables used in
    the embedding of *pols* and *sols*.  This *dim* also equals the
    dimension of the solution set represented by the witness set
    given by the lists *pols* and *sols*.
    The symbols for the slack variables are swapped to the end of the
    symbol table in both the polynomials and the solutions.
    """
    from phcpy.phcpy2c3 import py2c_swap_symbols_for_standard_witness_set
    store_standard_system(pols)
    store_standard_solutions(len(pols), sols)
    py2c_swap_symbols_for_standard_witness_set(nbvar, dim)

def store_dobldobl_witness_set(nbvar, dim, pols, sols):
    r"""
    Given in *nbar* is the total number of variables in the list of
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
    from phcpy.phcpy2c3 import py2c_swap_symbols_for_dobldobl_witness_set
    store_dobldobl_system(pols)
    store_dobldobl_solutions(len(pols), sols)
    py2c_swap_symbols_for_dobldobl_witness_set(nbvar, dim)

def store_quaddobl_witness_set(nbvar, dim, pols, sols):
    r"""
    Given in *nbar* is the total number of variables in the list of
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
    from phcpy.phcpy2c3 import py2c_swap_symbols_for_quaddobl_witness_set
    store_quaddobl_system(pols)
    store_quaddobl_solutions(len(pols), sols)
    py2c_swap_symbols_for_quaddobl_witness_set(nbvar, dim)

def store_standard_laurent_witness_set(nbvar, dim, pols, sols):
    r"""
    Given in *nbar* is the total number of variables in the list of
    Laurent polynomials in *pols* and its list of solutions in *sols*.
    The coefficients in the Laurent polynomials and the coordinates of the
    solutions will be parsed and stored in standard double precision.
    The parameter *dim* equals the number of slack variables used in
    the embedding of *pols* and *sols*.  This *dim* also equals the
    dimension of the solution set represented by the witness set
    given by the lists *pols* and *sols*.
    The symbols for the slack variables are swapped to the end of the
    symbol table in both the Laurent polynomials and the solutions.
    """
    from phcpy.phcpy2c3 \
    import py2c_swap_symbols_for_standard_Laurent_witness_set
    store_standard_laurent_system(pols)
    store_standard_solutions(len(pols), sols)
    py2c_swap_symbols_for_standard_Laurent_witness_set(nbvar, dim)

def store_dobldobl_laurent_witness_set(nbvar, dim, pols, sols):
    r"""
    Given in *nbar* is the total number of variables in the list of
    Laurent polynomials in *pols* and its list of solutions in *sols*.
    The coefficients in the Laurent polynomials and the coordinates of the
    solutions will be parsed and stored in double double precision.
    The parameter *dim* equals the number of slack variables used in
    the embedding of *pols* and *sols*.  This *dim* also equals the
    dimension of the solution set represented by the witness set
    given by the lists *pols* and *sols*.
    The symbols for the slack variables are swapped to the end of the
    symbol table in both the Laurent polynomials and the solutions.
    """
    from phcpy.phcpy2c3 \
    import py2c_swap_symbols_for_dobldobl_Laurent_witness_set
    store_dobldobl_laurent_system(pols)
    store_dobldobl_solutions(len(pols), sols)
    py2c_swap_symbols_for_dobldobl_Laurent_witness_set(nbvar, dim)

def store_quaddobl_laurent_witness_set(nbvar, dim, pols, sols):
    r"""
    Given in *nbar* is the total number of variables in the list of
    Laurent polynomials in *pols* and its list of solutions in *sols*.
    The coefficients in the Laurent polynomials and the coordinates of the
    solutions will be parsed and stored in quad double precision.
    The parameter *dim* equals the number of slack variables used in
    the embedding of *pols* and *sols*.  This *dim* also equals the
    dimension of the solution set represented by the witness set
    given by the lists *pols* and *sols*.
    The symbols for the slack variables are swapped to the end of the
    symbol table in both the Laurent polynomials and the solutions.
    """
    from phcpy.phcpy2c3 \
    import py2c_swap_symbols_for_quaddobl_Laurent_witness_set
    store_quaddobl_laurent_system(pols)
    store_quaddobl_solutions(len(pols), sols)
    py2c_swap_symbols_for_quaddobl_Laurent_witness_set(nbvar, dim)

def read_standard_solutions(filename):
    """
    Returns the list of solutions stored on file with the given file name.
    The solutions are parsed in standard double precision.
    """
    from phcpy.phcpy2c3 import py2c_solcon_read_standard_solutions_from_file
    namelen = len(filename)
    py2c_solcon_read_standard_solutions_from_file(namelen, filename)
    return load_standard_solutions()

def read_dobldobl_solutions(filename):
    """
    Returns the list of solutions stored on file with the given file name.
    The solutions are parsed in double double precision.
    """
    from phcpy.phcpy2c3 import py2c_solcon_read_dobldobl_solutions_from_file
    namelen = len(filename)
    py2c_solcon_read_dobldobl_solutions_from_file(namelen, filename)
    return load_dobldobl_solutions()

def read_quaddobl_solutions(filename):
    """
    Returns the list of solutions stored on file with the given file name.
    The solutions are parsed in quad double precision.
    """
    from phcpy.phcpy2c3 import py2c_solcon_read_quaddobl_solutions_from_file
    namelen = len(filename)
    py2c_solcon_read_quaddobl_solutions_from_file(namelen, filename)
    return load_quaddobl_solutions()

def test(prc='d', laurent=False):
    """
    Tests the storing of a witness set for the twisted cubic.
    The embedding induces the order x, y, zz1, z on the variables.
    After storing the witness set, the order is x, y, z, zz1,
    in both the system and solutions.
    The default precision prc is double 'd'.
    Other supported precisions are double double 'dd' and quad double 'qd'.
    """
    print('testing the storing of a witness set for the twisted cubic')
    e0 = 'x^2 - y + (9.22092060529474E-01 + 3.86970582743067E-01*i)*zz1;'
    e1 = 'x^3 - z + (6.47649182027044E-01-7.61938670117025E-01*i)*zz1;'
    e2 = 'zz1;'
    e3 = ' + (-7.98273302152795E-01-6.02295388551227E-01*i)*x' + \
         ' + (-9.99580934345184E-01 + 2.89474643727148E-02*i)*y' + \
         ' + (9.36299591763384E-01 + 3.51202326962281E-01*i)*z' + \
         ' + (8.34685978078231E-01-5.50726173338063E-01*i)*zz1' + \
         ' + (-6.91333588932132E-01 + 7.22535721479719E-01*i);'
    pols = [e0, e1, e2, e3]
    s0 = \
"""
t :  1.00000000000000E+00   0.00000000000000E+00
m : 1
the solution for t :
 x :  1.68397283871286E+00  -3.40775544483616E-01
 y :  2.71963654980455E+00  -1.14771352201599E+00
 zz1 : -2.76305054422513E-32  -1.43165074087146E-32
 z :  4.18868138066541E+00  -2.85950402375559E+00
== err :  1.663E-15 = rco :  1.823E-02 = res :  2.109E-15 =
"""
    s1 = \
"""
t :  1.00000000000000E+00   0.00000000000000E+00
m : 1
the solution for t :
 x : -7.91596946923866E-01  -6.16514124459453E-01
 y :  2.46536060721179E-01   9.76061397315087E-01
 zz1 :  0.00000000000000E+00   0.00000000000000E+00
 z :  4.06598444810859E-01  -9.24640185748066E-01
== err :  2.293E-16 = rco :  8.842E-02 = res :  2.776E-16 =
"""
    s2 = \
"""
t :  1.00000000000000E+00   0.00000000000000E+00
m : 1
the solution for t :
 x :  3.33649121255065E-02   5.79131019739152E-01
 y : -3.34279520662967E-01   3.86453111655036E-02
 zz1 : -2.66882752423716E-33   9.35041489367531E-33
 z : -3.35339052956912E-02  -1.92302242268359E-01
== err :  1.570E-16 = rco :  8.857E-02 = res :  1.457E-16 =
"""
    sols = [s0[1:-1], s1[1:-1], s2[1:-1]]
    if laurent:
        if prc == 'd':
            store_standard_laurent_system(pols)
            pols = load_standard_laurent_system()
        elif prc == 'dd':
            store_dobldobl_laurent_system(pols)
            pols = load_dobldobl_laurent_system()
        elif prc == 'qd':
            store_quaddobl_laurent_system(pols)
            pols = load_quaddobl_laurent_system()
        else:
            print('invalid precision level')
    else:
        if prc == 'd':
            store_standard_system(pols)
            pols = load_standard_system()
        elif prc == 'dd':
            store_dobldobl_system(pols)
            pols = load_dobldobl_system()
        elif prc == 'qd':
            store_quaddobl_system(pols)
            pols = load_quaddobl_system()
        else:
            print('invalid precision level')
    print('the embedded polynomials :')
    for pol in pols:
        print(pol)
    print('the generic points :')
    for sol in sols:
        print(sol)
    input('hit enter to continue')
    if laurent:
        if prc == 'd':
            store_standard_laurent_witness_set(4, 1, pols, sols)
            storedpols = load_standard_laurent_system()
        elif prc == 'dd':
            store_dobldobl_laurent_witness_set(4, 1, pols, sols)
            storedpols = load_dobldobl_laurent_system()
        elif prc == 'qd':
            store_quaddobl_laurent_witness_set(4, 1, pols, sols)
            storedpols = load_quaddobl_laurent_system()
        else:
            print('invalid precision level')
    else:
        if prc == 'd':
            store_standard_witness_set(4, 1, pols, sols)
            storedpols = load_standard_system()
        elif prc == 'dd':
            store_dobldobl_witness_set(4, 1, pols, sols)
            storedpols = load_dobldobl_system()
        elif prc == 'qd':
            store_quaddobl_witness_set(4, 1, pols, sols)
            storedpols = load_quaddobl_system()
        else:
            print('invalid precision level')
    for pol in storedpols:
        print(pol)
    if prc == 'd':
        storedsols = load_standard_solutions()
    elif prc == 'dd':
        storedsols = load_dobldobl_solutions()
    elif prc == 'qd':
        storedsols = load_quaddobl_solutions()
    else:
        print('invalid precision level')
    for sol in storedsols:
        print(sol)

def test_tableau():
    """
    Tests on storing and loading of a tableau.
    """
    from phcpy.phcpy2c3 import py2c_syscon_random_system
    from phcpy.phcpy2c3 import py2c_syscon_clear_standard_system
    dim, nbrmon, deg, cff = 2, 4, 3, 0
    py2c_syscon_random_system(dim, nbrmon, deg, cff)
    pols = load_standard_system()
    print('a random polynomial system :\n', pols)
    neq, nvr, ptab = load_standard_tableau()
    print('its tableau form :\n', ptab)
    py2c_syscon_clear_standard_system()
    store_standard_tableau(ptab)
    storedpols = load_standard_system()
    print('after clearing the container :\n', storedpols)

if __name__ == "__main__":
    test('d')
    test('d', True)
    test('dd')
    test('dd', True)
    test('qd')
    test('qd', True)
