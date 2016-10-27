"""
This module allows to work with monomial maps, defined by binomial systems.
"""

def is_binomial_system(silent=True):
    r"""
    Returns True if the system stored in the Laurent systems
    container is a binomial system, returns False otherwise.
    if not *silent*, then the number of terms in each Laurent
    polynomial is written to screen.
    """
    from phcpy.phcpy2c3 import py2c_syscon_number_of_standard_Laurentials
    from phcpy.phcpy2c3 import py2c_syscon_number_of_Laurent_terms
    nbequ = py2c_syscon_number_of_standard_Laurentials()
    if not silent:
        print('checking if binomial system ...')
        print('  number of Laurent polynomials :', nbequ)
    for i in range(1, nbequ+1):
        nbterms = py2c_syscon_number_of_Laurent_terms(i)
        if not silent:
            print('  -> number of terms in polynomial', i, ':', nbterms)
        if(nbterms != 2):
            if not silent:
                print('  the system is not a binomial system')
            return False
    if not silent:
        print('  the system is a binomial system')
    return True

def store_laurent_system(nbvar, pols):
    r"""
    Given in *pols* a list of string representing Laurent polynomials
    into the systems container.  The number of variables equals *nbvar*.
    """
    from phcpy.phcpy2c3 import py2c_syscon_clear_standard_Laurent_system
    from phcpy.phcpy2c3 \
        import py2c_syscon_initialize_number_of_standard_Laurentials
    from phcpy.phcpy2c3 import py2c_syscon_store_standard_Laurential
    py2c_syscon_clear_standard_Laurent_system()
    nbequ = len(pols)
    py2c_syscon_initialize_number_of_standard_Laurentials(nbequ)
    for ind in range(0, nbequ):
        pol = pols[ind]
        # print 'storing' , pol
        nbchar = len(pol)
        py2c_syscon_store_standard_Laurential(nbchar, nbvar, ind+1, pol)

def monomial_map_strings(dim, ind, nbvar):
    r"""
    Returns the list of strings representing the components of
    the monomial map of dimension *dim*, with index *ind*, and
    where the number of variables equals *nbvar*.
    """
    from phcpy.phcpy2c3 import py2c_mapcon_coefficients_of_map
    from phcpy.phcpy2c3 import py2c_mapcon_exponents_of_map
    from phcpy.phcpy2c3 import py2c_syscon_string_of_symbols
    sym = py2c_syscon_string_of_symbols()
    symvars = sym.split(' ')
    # py2c_mapcon_write_maps may have inserted t-varables
    thevars = symvars[-nbvar:]
    # print '-> the variables', symvars, thevars
    cff = py2c_mapcon_coefficients_of_map(dim, ind, nbvar)
    # print '-> coefficients of map', ind, ':', cff
    exp = py2c_mapcon_exponents_of_map(dim, ind, nbvar)
    # print '-> exponents of map', ind, ':', exp
    result = []
    exp_ind = 0
    for i in range(0, nbvar):
        str_var = thevars[i]
        if(cff[i] == 0):
            str_var = str_var + ' - 0'
            exp_ind = exp_ind + dim
        else:
            str_var = str_var + ' - ' + str(cff[i])
            for j in range(0, dim):
                pwr = exp[exp_ind]
                if(pwr != 0):
                    str_var = str_var + '*' + 't' + str(j+1) + '**' + str(pwr)
                exp_ind = exp_ind + 1
        result.append(str_var)
    return result

def write_monomial_map(dim, ind, nbvar):
    r"""
    Write the monomial map of dimension *dim* and of index *ind*,
    with number of variables equal to *nbvar*.
    """
    str_map = monomial_map_strings(dim, ind, nbvar)
    print(str_map)
    for str_var in str_map:
        print(str_var)

def monomial_map_solutions(nbvar, with_degree=True):
    r"""
    Returns the list of lists of strings,
    each list of strings representing a monomial map
    stored in the container.
    The number of variables equals *nbvar*.
    """
    from phcpy.phcpy2c3 import py2c_mapcon_top_dimension
    from phcpy.phcpy2c3 import py2c_mapcon_number_of_maps
    from phcpy.phcpy2c3 import py2c_mapcon_degree_of_map
    result = []
    topdim = py2c_mapcon_top_dimension()
    for dim in range(topdim, -1, -1):
        nbmaps = py2c_mapcon_number_of_maps(dim)
        for ind in range(1, nbmaps+1):
            str_map = monomial_map_strings(dim, ind, nbvar)
            str_map.append('dimension = %d' % dim)
            if with_degree:
                degmap = py2c_mapcon_degree_of_map(dim, ind)
                str_map.append('degree = %d' % degmap)
            result.append(str_map)
    return result

def write_monomial_maps(nbvar):
    r"""
    Writes the maps stored in the container.
    The number of variables is given in *nbvar*.
    """
    from phcpy.phcpy2c3 import py2c_mapcon_top_dimension
    from phcpy.phcpy2c3 import py2c_mapcon_number_of_maps
    from phcpy.phcpy2c3 import py2c_mapcon_degree_of_map
    topdim = py2c_mapcon_top_dimension()
    print('the top dimension :', topdim)
    for dim in range(topdim, -1, -1):
        nbmaps = py2c_mapcon_number_of_maps(dim)
        print('number of maps of dimension', dim, ':', nbmaps)
        for ind in range(1, nbmaps+1):
            print('monomial map', ind, 'of dimension', dim, ':')
            write_monomial_map(dim, ind, nbvar)
            degmap = py2c_mapcon_degree_of_map(dim, ind)
            print('degree of map', ind, ':', degmap)

def solve_binomials(nbvar, pols, silent=True, puretopdim=False):
    r"""
    If the system given in *pols* as a list of strings in as many
    variables as the value of *nbvar* is a binomial system
    (that is: it has exactly two monomials with a nonzero coefficient
    in every equation), then this function will return monomial maps
    to represent the solution sets.
    By default, *silent* is True and no additional output is written.
    If only the expected pure top dimensional solution sets are of interest,
    then switch the default *puretopdim* to True for faster results.
    The expected top dimension equals the number of variables minus
    the number of equations.
    """
    from phcpy.phcpy2c3 import py2c_syscon_write_standard_Laurent_system
    from phcpy.phcpy2c3 import py2c_mapcon_solve_system
    from phcpy.phcpy2c3 import py2c_mapcon_write_maps
    from phcpy.phcpy2c3 import py2c_mapcon_clear_maps
    store_laurent_system(nbvar, pols)
    if not silent:
        print('the polynomials on input :')
        py2c_syscon_write_standard_Laurent_system()
    isbin = is_binomial_system(silent)
    result = []
    if isbin:
        py2c_mapcon_solve_system(puretopdim)
        if not silent:
            py2c_mapcon_write_maps()
            write_monomial_maps(nbvar)
        result = monomial_map_solutions(nbvar)
        py2c_mapcon_clear_maps()
    return result

def test():
    """
    Solves a binomial system which has the x-axis, the yz-plane,
    and the twisted cubic as solution components.
    The yz-plane is a solution set of the unexpected dimension 2.
    """
    twisted2 = ['x**2*y - z*x;', 'x**2*z - y**2*x;']
    maps = solve_binomials(3, twisted2, silent=False)
    for solmap in maps:
        print(solmap)
    # note that the expected top dimension is 1,
    # the output below contains one curve on a 2-dimensional set
    print('looking only for expected pure dimensional sets ...')
    maps = solve_binomials(3, twisted2, puretopdim=True)
    for solmap in maps:
        print(solmap)

if __name__ == "__main__":
    test()
