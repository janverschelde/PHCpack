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
    for i in range(nbr):
        coords = [randint(low, upp) for k in range(dim)]
        point = tuple(coords)
        result.append(point)
    return result

def planar_convex_hull(points):
    """
    The convex hull of a point configuration in the plane
    consists of an ordered list of vertex points, ordered
    such that any two consecutive points span an edge,
    with the list of corresponding inner normals.
    """
    from phcpy2c import py2c_giftwrap_planar
    strpoints = str(points)
    strhull = py2c_giftwrap_planar(len(strpoints), strpoints)
    hull = eval(strhull)
    return hull

def convex_hull(dim, points):
    """
    Returns the list of facets of the convex hull of the points,
    given in points.  The dimension of the ambient space is in dim.
    """
    from phcpy2c import py2c_giftwrap_convex_hull
    from phcpy2c import py2c_giftwrap_number_of_facets
    from phcpy2c import py2c_giftwrap_retrieve_facet
    from phcpy2c import py2c_giftwrap_clear_3d_facets
    from phcpy2c import py2c_giftwrap_clear_3d_facets
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
        from phcpy2c import py2c_giftwrap_clear_3d_facets
        fail = py2c_giftwrap_clear_3d_facets()
    if(dim == 4):
        from phcpy2c import py2c_giftwrap_clear_4d_facets
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
