"""
This module provides the data interface to PHCpack.
The user gives as input string representations of polynomials or solutions
to the interface functions which store the data.
"""

def store_standard_system(polsys, **nbvar):
    """
    Stores the polynomials represented by the list of strings in polsys into
    the container for systems with coefficients in standard double precision.
    The number of variables is an optional argument given in nbvar.
    If nbvar is omitted, then the system is assumed to be square.
    Otherwise, suppose the number of variables equals 2 and pols is the list
    of polynomials, then the call store_standard_system(pols, nbvar=2)
    will store the polynomials in pols in the standard systems container.
    """
    from phcpy.phcpy2c2 import py2c_syscon_clear_standard_system
    from phcpy.phcpy2c2 \
        import py2c_syscon_initialize_number_of_standard_polynomials
    from phcpy.phcpy2c2 import py2c_syscon_store_standard_polynomial
    py2c_syscon_clear_standard_system()
    dim = len(polsys)
    py2c_syscon_initialize_number_of_standard_polynomials(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        if(len(nbvar) == 0):
            fail = py2c_syscon_store_standard_polynomial(nchar, dim, cnt+1, pol)
        else:
            nvr = nbvar.values()[0]
            fail = py2c_syscon_store_standard_polynomial(nchar, nvr, cnt+1, pol)
        if(fail != 0):
            break
    return fail

def store_dobldobl_system(polsys, **nbvar):
    """
    Stores the polynomials represented by the list of strings in polsys
    into the systems container for double double arithmetic.
    The number of variables is an optional argument given in nbvar.
    If nbvar is omitted, then the system is assumed to be square.
    Otherwise, suppose the number of variables equals 2 and pols is the list
    of polynomials, then the call store_dobldobl_system(pols, nbvar=2)
    will store the polynomials in pols in the dobldobl systems container.
    """
    from phcpy.phcpy2c2 import py2c_syscon_clear_dobldobl_system
    from phcpy.phcpy2c2 \
        import py2c_syscon_initialize_number_of_dobldobl_polynomials
    from phcpy.phcpy2c2 import py2c_syscon_store_dobldobl_polynomial
    py2c_syscon_clear_dobldobl_system()
    dim = len(polsys)
    py2c_syscon_initialize_number_of_dobldobl_polynomials(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        if(len(nbvar) == 0):
            fail = py2c_syscon_store_dobldobl_polynomial(nchar, dim, cnt+1, pol)
        else:
            nvr = nbvar.values()[0]
            fail = py2c_syscon_store_dobldobl_polynomial(nchar, nvr, cnt+1, pol)
        if(fail != 0):
            break
    return fail

def store_quaddobl_system(polsys, **nbvar):
    """
    Stores the polynomials represented by the list of strings in polsys
    into the systems container for quad double arithmetic.
    The number of variables is an optional argument given in nbvar.
    If nbvar is omitted, then the system is assumed to be square.
    Otherwise, suppose the number of variables equals 2 and pols is the list
    of polynomials, then the call store_quaddobl_system(pols, nbvar=2)
    will store the polynomials in pols in the quaddobl systems container.
    """
    from phcpy.phcpy2c2 import py2c_syscon_clear_quaddobl_system
    from phcpy.phcpy2c2 \
        import py2c_syscon_initialize_number_of_quaddobl_polynomials
    from phcpy.phcpy2c2 import py2c_syscon_store_quaddobl_polynomial
    py2c_syscon_clear_quaddobl_system()
    dim = len(polsys)
    py2c_syscon_initialize_number_of_quaddobl_polynomials(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        if(len(nbvar) == 0):
            fail = py2c_syscon_store_quaddobl_polynomial(nchar, dim, cnt+1, pol)
        else:
            nvr = nbvar.values()[0]
            fail = py2c_syscon_store_quaddobl_polynomial(nchar, nvr, cnt+1, pol)
        if(fail != 0):
            break
    return fail

def store_multprec_system(polsys, decimals, **nbvar):
    """
    Stores the polynomials represented by the list of strings in polsys
    into the systems container for multiprecision arithmetic.
    The parameter decimals equals the number of decimal places
    in the working precision for the parsing of the strings in polsys.
    The number of variables is an optional argument given in nbvar.
    If nbvar is omitted, then the system is assumed to be square.
    Otherwise, suppose the number of variables equals 2 and pols is the list of
    polynomials, then the call store_multprec_system(pols, nbvar=2)
    will store the polynomials in pols in the multiprecision systems container.
    """
    from phcpy.phcpy2c2 import py2c_syscon_clear_multprec_system
    from phcpy.phcpy2c2 \
        import py2c_syscon_initialize_number_of_multprec_polynomials
    from phcpy.phcpy2c2 import py2c_syscon_store_multprec_polynomial
    py2c_syscon_clear_multprec_system()
    dim = len(polsys)
    py2c_syscon_initialize_number_of_multprec_polynomials(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        if(len(nbvar) == 0):
            fail = py2c_syscon_store_multprec_polynomial\
                       (nchar, dim, cnt+1, decimals, pol)
        else:
            nvr = nbvar.values()[0]
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
    from phcpy.phcpy2c2 import py2c_syscon_number_of_standard_polynomials
    from phcpy.phcpy2c2 import py2c_syscon_load_standard_polynomial
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
    from phcpy.phcpy2c2 import py2c_syscon_number_of_dobldobl_polynomials
    from phcpy.phcpy2c2 import py2c_syscon_load_dobldobl_polynomial
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
    from phcpy.phcpy2c2 import py2c_syscon_number_of_quaddobl_polynomials
    from phcpy.phcpy2c2 import py2c_syscon_load_quaddobl_polynomial
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
    from phcpy.phcpy2c2 import py2c_syscon_number_of_multprec_polynomials
    from phcpy.phcpy2c2 import py2c_syscon_load_multprec_polynomial
    dim = py2c_syscon_number_of_multprec_polynomials()
    result = []
    for ind in range(1, dim+1):
        result.append(py2c_syscon_load_multprec_polynomial(ind))
    return result

def store_standard_laurent_system(polsys):
    """
    Stores the Laurent polynomials represented by the list of
    strings in polsys into the container for systems
    with coefficients in standard double precision.
    """
    from phcpy.phcpy2c2 import py2c_syscon_clear_standard_Laurent_system
    from phcpy.phcpy2c2 \
        import py2c_syscon_initialize_number_of_standard_Laurentials
    from phcpy.phcpy2c2 import py2c_syscon_store_standard_Laurential
    py2c_syscon_clear_standard_Laurent_system()
    dim = len(polsys)
    py2c_syscon_initialize_number_of_standard_Laurentials(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        fail = py2c_syscon_store_standard_Laurential(nchar, dim, cnt+1, pol)
        if(fail != 0):
            break
    return fail

def store_dobldobl_laurent_system(polsys):
    """
    Stores the Laurent polynomials represented by the list of
    strings in polsys into the container for systems
    with coefficients in double double precision.
    """
    from phcpy.phcpy2c2 import py2c_syscon_clear_dobldobl_Laurent_system
    from phcpy.phcpy2c2 \
        import py2c_syscon_initialize_number_of_dobldobl_Laurentials
    from phcpy.phcpy2c2 import py2c_syscon_store_dobldobl_Laurential
    py2c_syscon_clear_dobldobl_Laurent_system()
    dim = len(polsys)
    py2c_syscon_initialize_number_of_dobldobl_Laurentials(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        fail = py2c_syscon_store_dobldobl_Laurential(nchar, dim, cnt+1, pol)
        if(fail != 0):
            break
    return fail

def store_quaddobl_laurent_system(polsys):
    """
    Stores the Laurent polynomials represented by the list
    of strings in polsys into the container for systems
    with coefficients in quad double precision.
    """
    from phcpy.phcpy2c2 import py2c_syscon_clear_quaddobl_Laurent_system
    from phcpy.phcpy2c2 \
        import py2c_syscon_initialize_number_of_quaddobl_Laurentials
    from phcpy.phcpy2c2 import py2c_syscon_store_quaddobl_Laurential
    py2c_syscon_clear_quaddobl_Laurent_system()
    dim = len(polsys)
    py2c_syscon_initialize_number_of_quaddobl_Laurentials(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        fail = py2c_syscon_store_quaddobl_Laurential(nchar, dim, cnt+1, pol)
        if(fail != 0):
            break
    return fail

def store_multprec_laurent_system(polsys, decimals):
    """
    Stores the Laurent polynomials represented by the list
    of strings in polsys into the container for systems
    with coefficients in multiprecision.
    The parameter decimals equals the number of decimal places
    in the working precision for the parsing of the strings in polsys.
    """
    from phcpy.phcpy2c2 import py2c_syscon_clear_multprec_Laurent_system
    from phcpy.phcpy2c2 \
        import py2c_syscon_initialize_number_of_multprec_Laurentials
    from phcpy.phcpy2c2 import py2c_syscon_store_multprec_Laurential
    py2c_syscon_clear_multprec_Laurent_system()
    dim = len(polsys)
    py2c_syscon_initialize_number_of_multprec_Laurentials(dim)
    for cnt in range(0, dim):
        pol = polsys[cnt]
        nchar = len(pol)
        fail = py2c_syscon_store_multprec_Laurential\
                   (nchar, dim, cnt+1, decimals, pol)
        if(fail != 0):
            break
    return fail

def load_standard_laurent_system():
    """
    Returns the Laurent polynomials stored in the system container
    for standard double precision arithmetic.
    """
    from phcpy.phcpy2c2 import py2c_syscon_number_of_standard_Laurentials
    from phcpy.phcpy2c2 import py2c_syscon_load_standard_Laurential
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
    from phcpy.phcpy2c2 import py2c_syscon_number_of_dobldobl_Laurentials
    from phcpy.phcpy2c2 import py2c_syscon_load_dobldobl_Laurential
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
    from phcpy.phcpy2c2 import py2c_syscon_number_of_quaddobl_Laurentials
    from phcpy.phcpy2c2 import py2c_syscon_load_quaddobl_Laurential
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
    from phcpy.phcpy2c2 import py2c_syscon_number_of_multprec_Laurentials
    from phcpy.phcpy2c2 import py2c_syscon_load_multprec_Laurential
    dim = py2c_syscon_number_of_multprec_Laurentials()
    result = []
    for ind in range(1, dim+1):
        result.append(py2c_syscon_load_multprec_Laurential(ind))
    return result

def store_standard_solutions(nvar, sols):
    """
    Stores the solutions in the list sols, represented as strings
    in PHCpack format into the container for solutions
    with standard double precision.
    The number nvar equals the number of variables.
    """
    from phcpy.phcpy2c2 import py2c_solcon_clear_standard_solutions
    from phcpy.phcpy2c2 import py2c_solcon_append_standard_solution_string
    py2c_solcon_clear_standard_solutions()
    for ind in range(0, len(sols)):
        fail = py2c_solcon_append_standard_solution_string\
                   (nvar, len(sols[ind]), sols[ind])
        if(fail != 0):
            break
    return fail

def store_dobldobl_solutions(nvar, sols):
    """
    Stores the solutions in the list sols, represented as strings
    in PHCpack format into the solution container for processing
    with complex double double arithmetic.
    The number nvar equals the number of variables.
    """
    from phcpy.phcpy2c2 import py2c_solcon_clear_dobldobl_solutions
    from phcpy.phcpy2c2 import py2c_solcon_append_dobldobl_solution_string
    py2c_solcon_clear_dobldobl_solutions()
    for ind in range(0, len(sols)):
        fail = py2c_solcon_append_dobldobl_solution_string\
                   (nvar, len(sols[ind]), sols[ind])
        if(fail != 0):
            break
    return fail

def store_quaddobl_solutions(nvar, sols):
    """
    Stores the solutions in the list sols, represented as strings
    in PHCpack format into the solution container for processing
    with complex quad double arithmetic.
    The number n equals the number of variables.
    """
    from phcpy.phcpy2c2 import py2c_solcon_clear_quaddobl_solutions
    from phcpy.phcpy2c2 import py2c_solcon_append_quaddobl_solution_string
    py2c_solcon_clear_quaddobl_solutions()
    for ind in range(0, len(sols)):
        fail = py2c_solcon_append_quaddobl_solution_string\
                   (nvar, len(sols[ind]), sols[ind])
        if(fail != 0):
            break
    return fail

def store_multprec_solutions(nvar, sols):
    """
    Stores the solutions in the list sols, represented as strings
    in PHCpack format into the solution container for processing
    with complex multiprecision arithmetic.
    The number n equals the number of variables.
    """
    from phcpy.phcpy2c2 import py2c_solcon_clear_multprec_solutions
    from phcpy.phcpy2c2 import py2c_solcon_append_multprec_solution_string
    py2c_solcon_clear_multprec_solutions()
    for ind in range(0, len(sols)):
        fail = py2c_solcon_append_multprec_solution_string\
                   (nvar, len(sols[ind]), sols[ind])
        if(fail != 0):
            break
    return fail

def load_standard_solutions():
    """
    Returns the list of solutions stored in the container
    for solutions with standard double precision.
    """
    from phcpy.phcpy2c2 import py2c_solcon_retrieve_next_standard_initialize
    from phcpy.phcpy2c2 import py2c_solcon_move_current_standard_to_next
    from phcpy.phcpy2c2 \
        import py2c_solcon_length_current_standard_solution_string
    from phcpy.phcpy2c2 \
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
    from phcpy.phcpy2c2 import py2c_solcon_retrieve_next_dobldobl_initialize
    from phcpy.phcpy2c2 import py2c_solcon_move_current_dobldobl_to_next
    from phcpy.phcpy2c2 \
        import py2c_solcon_length_current_dobldobl_solution_string
    from phcpy.phcpy2c2 \
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
    from phcpy.phcpy2c2 import py2c_solcon_retrieve_next_quaddobl_initialize
    from phcpy.phcpy2c2 import py2c_solcon_move_current_quaddobl_to_next
    from phcpy.phcpy2c2 \
        import py2c_solcon_length_current_quaddobl_solution_string
    from phcpy.phcpy2c2 \
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
    from phcpy.phcpy2c2 import py2c_solcon_retrieve_next_multprec_initialize
    from phcpy.phcpy2c2 import py2c_solcon_move_current_multprec_to_next
    from phcpy.phcpy2c2 \
        import py2c_solcon_length_current_multprec_solution_string
    from phcpy.phcpy2c2 \
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
