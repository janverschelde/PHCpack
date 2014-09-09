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

if __name__ == "__main__":
    test_planar_hull()
    test_convex_hull()
