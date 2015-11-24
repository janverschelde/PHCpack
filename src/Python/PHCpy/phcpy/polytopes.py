"""
This module exports routines of PHCpack to work with Newton polytopes.
"""

def random_points(dim, nbr, low, upp):
    """
    Generates a list of random integer points.
    Returns a list of nbr points in dimension dim,
    with integer coordinates in the range low..upp.
    """
    from random import randint
    result = []
    for _ in range(nbr):
        coords = [randint(low, upp) for _ in range(dim)]
        point = tuple(coords)
        result.append(point)
    return result

def support(nvr, pol):
    """
    The support of a multivariate polynomial is a set of exponents
    of the monomials that appear with nonzero coefficient.
    Given in nvr the number of variables and in pol a string
    representation of a polynomial in nvr variables,
    returns the support of the polynomial as a list of tuples.
    """
    from phcpy.phcpy2c import py2c_syscon_clear_standard_Laurent_system
    from phcpy.phcpy2c \
    import py2c_syscon_initialize_number_of_standard_Laurentials
    from phcpy.phcpy2c import py2c_syscon_store_standard_Laurential
    py2c_syscon_clear_standard_Laurent_system()
    py2c_syscon_initialize_number_of_standard_Laurentials(nvr)
    nchar = len(pol)
    fail = py2c_syscon_store_standard_Laurential(nchar, nvr, 1, pol)
    if(fail != 0):
        return fail
    else:
        # from phcpy.phcpy2c import py2c_syscon_number_of_Laurent_terms
        from phcpy.phcpy2c import py2c_giftwrap_support_size
        from phcpy.phcpy2c import py2c_giftwrap_support_string
        from phcpy.phcpy2c import py2c_giftwrap_clear_support_string
        # ntm = py2c_syscon_number_of_Laurent_terms(1)
        # print 'number of terms : %d' % ntm
        size = py2c_giftwrap_support_size()
        # print 'size of support :', size
        supp = py2c_giftwrap_support_string(size)
        # print 'support string :', supp
        py2c_giftwrap_clear_support_string()
        result = eval(supp)
        return result

def initial_form(pols, normal):
    """
    Returns the initial form of the polynomials in pols
    with respect to the inner normal with coordinates in normal.
    """
    from phcpy.phcpy2c import py2c_syscon_clear_standard_Laurent_system
    from phcpy.phcpy2c \
    import py2c_syscon_initialize_number_of_standard_Laurentials
    from phcpy.phcpy2c import py2c_syscon_store_standard_Laurential
    from phcpy.phcpy2c import py2c_syscon_load_standard_Laurential
    from phcpy.phcpy2c import py2c_giftwrap_initial_form
    py2c_syscon_clear_standard_Laurent_system()
    dim = max(len(pols), len(normal))
    py2c_syscon_initialize_number_of_standard_Laurentials(dim)
    for i in range(len(pols)):
        nchar = len(pols[i])
        fail = py2c_syscon_store_standard_Laurential(nchar, dim, i+1, pols[i])
    strnrm = str(tuple(normal))
    fail = py2c_giftwrap_initial_form(len(normal), len(strnrm), strnrm)
    result = []
    for i in range(len(pols)):
        result.append(py2c_syscon_load_standard_Laurential(i+1))
    return result

def initial_support(points, normal):
    """
    Returns the list of elements in points that make the minimal inner
    product with the given normal, as the second element in a tuple.
    The first element is the value of that minimal inner product.
    Every tuple in points must have the same length as normal.
    """
    inprod = lambda x, y: sum([a*b for (a, b) in zip(x, y)])
    vals = [inprod(point, normal) for point in points]
    minv = min(vals)
    idx = [ind for ind in range(len(vals)) if vals[ind] == minv]
    result = [points[ind] for ind in idx]
    return (minv, result)

def convex_hull_checkin(dim, points):
    """
    Checks whether the input arguments satisfy the requirements:
    points is a list of tuples that each contain as many integer
    numbers as the value of dim.
    Returns True if the requirements are satisfied,
    returns False otherwise.
    """
    if not isinstance(points, list):
        print 'the argument points is not a list'
        return False
    else:
        tup = [isinstance(x, tuple) for x in points]
        if(sum(tup) != len(points)):
            print 'not every element in points is a tuple'
        else:
            for point in points:
                if(len(point) != dim):
                    print 'the point', point, 'is not of length', dim
                else:
                    coord = [isinstance(x, int) for x in point]
                    if(sum(coord) != dim):
                        print point, 'contains non integer values'
    return True

def planar_convex_hull(points, checkin=True):
    """
    The convex hull of a point configuration in the plane
    consists of an ordered list of vertex points, ordered
    such that any two consecutive points span an edge,
    with the list of corresponding inner normals.
    If checkin (by default), the type of the input is checked.
    """
    if checkin:
        if not convex_hull_checkin(2, points):
            return None
    from phcpy.phcpy2c import py2c_giftwrap_planar
    strpoints = str(points)
    strhull = py2c_giftwrap_planar(len(strpoints), strpoints)
    hull = eval(strhull)
    return hull

def convex_hull(dim, points, checkin=True):
    """
    Returns the list of facets of the convex hull of the points,
    given in points.  The dimension of the ambient space is in dim.
    If checkin (by default), the type of the input is checked.
    """
    if checkin:
        if not convex_hull_checkin(dim, points):
            return None
    from phcpy.phcpy2c import py2c_giftwrap_convex_hull
    from phcpy.phcpy2c import py2c_giftwrap_number_of_facets
    from phcpy.phcpy2c import py2c_giftwrap_retrieve_facet
    strpoints = str(points)
    fail = py2c_giftwrap_convex_hull(len(strpoints), strpoints)
    nbrfacets = py2c_giftwrap_number_of_facets(dim)
    print 'computed', nbrfacets, 'facets'
    result = []
    for k in range(nbrfacets):
        strfacet = py2c_giftwrap_retrieve_facet(dim, k)
        facet = eval(strfacet)
        result.append(facet)
    if(dim == 3):
        from phcpy.phcpy2c import py2c_giftwrap_clear_3d_facets
        fail = py2c_giftwrap_clear_3d_facets()
    if(dim == 4):
        from phcpy.phcpy2c import py2c_giftwrap_clear_4d_facets
        fail = py2c_giftwrap_clear_4d_facets()
    return result

def mixed_volume(mixture, points):
    """
    Returns the mixed volume of the tuple in points.
    Both mixture and points have the same length.
    The list mixture counts the number of times each support
    in points should be counted.
    For example, to compute the volume of a three dimensional polytope,
    the mixture is [3].  In general, the mixture determines the powers
    of the unknowns in the Minkowski polynomial of which the computed
    mixed volume is its coefficient.
    """
    from phcpy.phcpy2c import py2c_celcon_initialize_supports as init
    from phcpy.phcpy2c import py2c_celcon_set_type_of_mixture as setmix
    from phcpy.phcpy2c import py2c_celcon_append_lifted_point as applft
    from phcpy.phcpy2c import py2c_celcon_mixed_volume_of_supports as mixvol
    nbr = len(mixture)
    init(nbr)
    setmix(nbr, str(mixture))
    for k in range(nbr):
        for point in points[k]:
            lpt = list(point)
            lpt.append(0)
            applft(len(lpt), k+1, str(lpt))
    return mixvol()

def test_planar_hull():
    """
    Generates a random point configuration in the plane
    and then computes its convex hull.
    """
    pts = random_points(2, 7, -9, 9)
    print 'the points :', pts
    (vertices, normals) = planar_convex_hull(pts)
    print 'the vertices :', vertices
    print 'inner normals :', normals

def test_convex_hull():
    """
    Generates a random point configuration in 3-space
    and then computes its convex hull.
    """
    pts = random_points(3, 10, -9, 9)
    print 'the points :', pts
    facets = convex_hull(3, pts)
    print 'the facets :'
    for facet in facets:
        print facet

def test_mixed_volume():
    """
    Runs some simple tests on mixed volume computation.
    """
    simple = [(0, 0, 0), (1, 0, 0), (0, 2, 0), (0, 0, 3)]
    mixvol = mixed_volume([3], tuple([simple]))
    print simple, 'has volume', mixvol
    points = random_points(3, 10, -9, 9)
    mixvol = mixed_volume([3], tuple([points]))
    print points, 'has volume', mixvol
    mixvol = mixed_volume([2, 1], (points, simple))
    print 'the mixed volume of the two is', mixvol
    cyclic5 = ([(1, 0, 0, 0, 0), (0, 1, 0, 0, 0), (0, 0, 1, 0, 0), \
                (0, 0, 0, 1, 0), (0, 0, 0, 0, 1)],
               [(1, 1, 0, 0, 0), (0, 1, 1, 0, 0), (0, 0, 1, 1, 0), \
                (0, 0, 0, 1, 1), (1, 0, 0, 0, 1)],
               [(1, 1, 1, 0, 0), (0, 1, 1, 1, 0), (0, 0, 1, 1, 1), \
                (1, 0, 0, 1, 1), (1, 1, 0, 0, 1)],
               [(1, 1, 1, 1, 0), (0, 1, 1, 1, 1), (1, 0, 1, 1, 1), \
                (1, 1, 0, 1, 1), (1, 1, 1, 0, 1)],
               [(1, 1, 1, 1, 1), (0, 0, 0, 0, 0)])
    mixvol = mixed_volume([1, 1, 1, 1, 1], cyclic5)
    print 'the mixed volume of the cyclic 5-roots problem is', mixvol

if __name__ == "__main__":
    test_planar_hull()
    test_convex_hull()
