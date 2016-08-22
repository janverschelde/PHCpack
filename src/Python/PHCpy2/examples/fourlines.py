"""
Given four lines in general position,
there are two lines which meet all four given lines.
With Pieri homotopies we can solve this Schubert problem.
For the verification of the intersection conditions, numpy is used.
The plots are made with matplotlib.
"""
from numpy import zeros, array, concatenate, matrix
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def indices(name):
    """
    For the string name in the format xij
    return (i, j) as two integer indices.
    """
    return (int(name[1]), int(name[2]))

def solution_plane(rows, cols, sol):
    """
    Returns a sympy matrix with as many rows
    as the value of rows and with as many columns
    as the value of columns, using the string
    represention of a solution in sol.
    """
    from phcpy.solutions import coordinates
    result = zeros((rows, cols), dtype=complex)
    for k in range(cols):
        result[k][k] = 1
    (vars, vals) = coordinates(sol)
    for (name, value) in zip(vars, vals):
        i, j = indices(name)
        result[i-1][j-1] = value
    return result

def verify_determinants(inps, sols, verbose=True):
    """
    Verifies the intersection conditions with determinants,
    concatenating the planes in inps with those in the sols.
    Both inps and sols are lists of numpy arrays.
    Returns the sum of the absolute values of all determinants.
    If verbose, then for all solutions in sols, the computed
    determinants are printed to screen.
    """
    from numpy.linalg import det
    checksum = 0
    for sol in sols:
        if verbose:
            print('checking solution\n', sol)
        for plane in inps:
            cat = concatenate([plane, sol], axis=-1)
            mat = matrix(cat)
            dcm = det(mat)
            if verbose:
                print('the determinant :', dcm)
            checksum = checksum + abs(dcm)
    return checksum

def solve_general(mdim, pdim, qdeg):
    """
    Solves a general instance of Pieri problem, computing the
    p-plane producing curves of degree qdeg which meet a number
    of general m-planes at general interpolation points,
    where p = pdim and m = mdim on input.
    For the problem of computing the two lines which meet
    four general lines, mdim = 2, pdim = 2, and qdeg = 0.
    Returns a tuple with four lists.
    The first two lists contain matrices with the input planes
    and the solution planes respectively.
    The third list is the list of polynomials solved
    and the last list is the solution list.
    """
    from phcpy.schubert import random_complex_matrix
    from phcpy.schubert import run_pieri_homotopies
    dim = mdim*pdim + qdeg*(mdim+pdim)
    ranplanes = [random_complex_matrix(mdim+pdim, mdim) for _ in range(0, dim)]
    (pols, sols) = run_pieri_homotopies(mdim, pdim, qdeg, ranplanes, \
        verbose=False)
    inplanes = [array(plane) for plane in ranplanes]
    outplanes = [solution_plane(mdim+pdim, pdim, sol) for sol in sols]
    return (inplanes, outplanes, pols, sols)

def solve_real(mdim, pdim, start, sols):
    """
    Solves a real instance of Pieri problem, for input planes
    of dimension mdim osculating a rational normal curve.
    On return are the planes of dimension pdim.
    """
    from phcpy.schubert import real_osculating_planes
    from phcpy.schubert import make_pieri_system
    from phcpy.trackers import track
    oscplanes = real_osculating_planes(mdim, pdim, 0)
    target = make_pieri_system(mdim, pdim, 0, oscplanes, False)
    rtsols = track(target, start, sols)
    inplanes = [array(plane) for plane in oscplanes]
    outplanes = [solution_plane(mdim+pdim, pdim, sol) for sol in rtsols]
    return (inplanes, outplanes, target, rtsols)

def input_generators(plane):
    """
    Given in plane is a numpy matrix, with in its columns
    the coordinates of the points which span a line, in 4-space.
    The first coordinate must not be zero.
    Returns the affine representation of the line,
    after dividing each generator by its first coordinate.
    """
    pone = list(plane[:,0])
    ptwo = list(plane[:,1])
    aone = [x/pone[0] for x in pone]
    atwo = [x/ptwo[0] for x in ptwo]
    return (aone[1:], atwo[1:])

def output_generators(plane):
    """
    Given in plane is a numpy matrix, with in its columns
    the coordinates of the points which span a line, in 4-space.
    The solution planes follow the localization pattern
    1, *, *, 0 for the first point and 0, 1, *, * for
    the second point, which means that the second point
    in standard projective coordinates lies at infinity.
    For the second generator, the sum of the points is taken.
    The imaginary part of each coordinate is omitted.
    """
    pone = list(plane[:,0])
    ptwo = list(plane[:,1])
    aone = [x.real for x in pone]
    atwo = [x.real + y.real for (x, y) in zip(pone, ptwo)]
    return (aone[1:], atwo[1:])

def boxrange(inlines, outlines):
    """
    Returns a list of three lists with the [min, max]
    values of each coordinate of each generator in the lists
    inlines and outlines.
    The ranges are adjusted for the particular real case.
    """
    fst = inlines[0][0]
    result = {'xmin': fst[0], 'xmax': fst[0], \
              'ymin': fst[1], 'ymax': fst[1], \
              'zmin': fst[2], 'zmax': fst[2]} 
    pts = [x for (x, y) in inlines] + [y for (x, y) in inlines] \
        + [x for (x, y) in outlines] + [y for (x, y) in outlines]
    print('the points :\n', pts)
    for point in pts:
        result['xmin'] = min(result['xmin'], point[0])
        result['ymin'] = min(result['ymin'], point[1])
        result['zmin'] = min(result['zmin'], point[2])
        result['xmax'] = max(result['xmax'], point[0])
        result['ymax'] = max(result['ymax'], point[1])
        result['zmax'] = max(result['zmax'], point[2])
    return ((result['xmin']+3, result['xmax']-3), \
            (result['ymin']+8, result['ymax']-11), \
            (result['zmin']+3, result['zmax']-5))

def inbox(point, lims):
    """
    Returns true if the coordinates of the point
    are in the box defined by the 3-tuple lims
    which contain the minima and maxima for the coordinates.
    """
    tol = 1.0e-8 # this is essential for roundoff
    (xlim, ylim, zlim) = lims
    if point[0] < xlim[0] - tol:
        return False
    elif point[0] > xlim[1] + tol:
        return False
    elif point[1] < ylim[0] - tol:
        return False
    elif point[1] > ylim[1] + tol:
        return False
    elif point[2] < zlim[0] - tol:
        return False
    elif point[2] > zlim[1] + tol:
        return False
    else:
        return True

def equal(pt1, pt2):
    """
    Returns true if the all coordinates of pt1 and pt2
    match up to a tolerance of 1.0e-10.
    """
    tol = 1.0e-8
    if abs(pt1[0] - pt2[0]) > tol:
        return False
    elif abs(pt1[1] - pt2[1]) > tol:
        return False
    elif abs(pt1[2] - pt2[2]) > tol:
        return False
    return True

def isin(points, pnt):
    """
    Returns true if pnt belongs to the list points.
    """
    if len(points) == 0:
        return False
    else:
        for point in points:
            if equal(point, pnt):
                return True
        return False;

def plot_line(axs, line, lims, color):
    """
    Plots the line defined as a tuple of two points,
    using the axis object in axs.
    The 3-tuple lims contains three lists with limits [min, max]
    for the x, y, and z coordinates.
    """
    (fst, snd) = line
    axs.set_xlabel('x')
    axs.set_ylabel('y')
    axs.set_zlabel('z')
    axs.set_xlim(lims[0])
    axs.set_ylim(lims[1])
    axs.set_zlim(lims[2])
    dir = (fst[0] - snd[0], fst[1] - snd[1], fst[2] - snd[2])
    result = []
    for k in range(3):
        fac = (lims[k][1]-fst[k])/dir[k]
        pnt = (fst[0] + fac*dir[0], fst[1] + fac*dir[1], fst[2] + fac*dir[2])
        if inbox(pnt, lims):
            if not isin(result, pnt): result.append(pnt)
    for k in range(3):
        fac = (lims[k][0]-fst[k])/dir[k]
        pnt = (fst[0] + fac*dir[0], fst[1] + fac*dir[1], fst[2] + fac*dir[2])
        if inbox(pnt, lims):
            if not isin(result, pnt): result.append(pnt)
    (one, two) = (result[0], result[1])
    # axs.plot([fst[0], snd[0]], [fst[1], snd[1]], [fst[2], snd[2]], 'bo')
    # axs.plot([one[0], two[0]], [one[1], two[1]], [one[2], two[2]], 'ro')
    axs.plot([one[0], two[0]], [one[1], two[1]], [one[2], two[2]], color)

def plot_lines(inlines, outlines, points, lims):
    """
    Generates coordinates of the points in a random line
    and then plots this line.  The intersection points are
    in the list points and limits for the bounding box in lims
    """
    fig = plt.figure()
    axs = fig.add_subplot(111, projection='3d')
    for line in inlines:
        plot_line(axs, line, lims, 'b')
    for line in outlines:
        plot_line(axs, line, lims, 'r')
    for point in points:
        axs.plot([point[0]], [point[1]], [point[2]], 'ro')
    axs.view_init(azim=5, elev=20)
    plt.show()

def intersection_point(apl, bpl, check=True):
    """
    Given in apl the two points that define a line
    and in bpl the two points that define another line,
    returns the intersection point.
    If check, then additional tests are done
    and the outcome of the tests is written to screen.
    """
    from numpy.linalg import solve
    (apt, bpt) = apl
    (cpt, dpt) = bpl
    mat = array([[apt[0], bpt[0], -cpt[0]], \
                 [apt[1], bpt[1], -cpt[1]], \
                 [apt[2], bpt[2], -cpt[2]]])
    rhs = array([[dpt[0]], [dpt[1]], [dpt[2]]])
    sol = solve(mat, rhs)
    cff = list(sol[:,0])
    csm = cff[0] + cff[1]
    result = ((cff[0]*apt[0] + cff[1]*bpt[0])/csm, \
              (cff[0]*apt[1] + cff[1]*bpt[1])/csm, \
              (cff[0]*apt[2] + cff[1]*bpt[2])/csm)
    if check:
        csm = cff[2] + 1.0
        verify = ((cff[2]*cpt[0] + dpt[0])/csm, \
                  (cff[2]*cpt[1] + dpt[1])/csm, \
                  (cff[2]*cpt[2] + dpt[2])/csm)
        print('the solution :\n', result)
        print('the solution verified :\n', verify)
        res = matrix(rhs) - matrix(mat)*matrix(sol)
        print('the residual :\n', res)
    return result

def intersection_points(ipl, opl):
    """
    Returns the list of intersection points between
    the input planes in ipl and the output planes in opl.
    """
    result = []
    for inplane in ipl:
        for outplane in opl:
            result.append(intersection_point(inplane, outplane))
    return result

def show_planes(ipl, opl):
    """
    Shows the input and the output planes.
    """
    (inlines, outlines) = ([], [])
    for plane in ipl:
        inlines.append(input_generators(plane))
    for plane in opl:
        outlines.append(output_generators(plane))
    print('The generators of the input lines :')
    for line in inlines:
        print(line)
    print('The generators of the output lines :')
    for line in outlines:
        print(line)
    brg = boxrange(inlines, outlines)
    print('the range:', brg)
    intpts = intersection_points(inlines, outlines)
    print('the intersection points :')
    for point in intpts:
        print(point)
    plot_lines(inlines, outlines, intpts, brg)

def main():
    """
    We start with the formalism of the root count,
    solve a general configuration and then a special problem.
    """
    from phcpy.schubert import pieri_root_count
    (mdim, pdim, deg) = (2, 2, 0)
    pcnt = pieri_root_count(mdim, pdim, deg, False)
    print('The Pieri root count :', pcnt)
    print('Solving a general case ...')
    (inp, otp, pols, sols) = solve_general(mdim, pdim, deg)
    print('The input planes :')
    for plane in inp:
        print(plane)
    print('The solution planes :')
    for plane in otp:
        print(plane)
    check = verify_determinants(inp, otp)
    print('Sum of absolute values of determinants :', check)
    input('Hit enter to continue.')
    from random import seed
    seed(400)
    (oscp, otp2, pols2, sols2) = solve_real(mdim, pdim, pols, sols)
    print('The input planes :')
    for plane in oscp:
        print(plane)
    print('The solution planes :')
    for plane in otp2:
        print(plane)
    check = verify_determinants(oscp, otp2)
    print('Sum of absolute values of determinants :', check)
    show_planes(oscp, otp2)

if __name__ == "__main__":
    main()
