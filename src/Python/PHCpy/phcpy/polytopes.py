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

if __name__ == "__main__":
    test_planar_hull()
