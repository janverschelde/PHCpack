"""
Sage 7.4 code to generate the polynomial system for the problem
of all lines tangent to four given spheres.
This scripts runs typing sage tangents4spheres.sage at the command prompt,
provided the python interpreter of sage is extended with phcpy.
"""
def tangent_system(centers, radii, verbose=True):
    """
    Returns a tuple with the list of polynomials
    and variables for the common tangents to four
    spheres with given centers and squares of radii.
    """
    x0, x1, x2 = var('x0, x1, x2')
    t = (x0, x1, x2) 
    vt = vector(t)   # tangent vector
    normt = vt.dot_product(vt) - 1
    if verbose:
        print 'normalized tangent :'
        print normt
    x3, x4, x5 = var('x3, x4, x5')
    m = (x3, x4, x5)
    vm = vector(m)   # moment vector
    momt = vt.dot_product(vm)
    if verbose:
        print 'moment is perpendicular to tangent :'
        print momt
    eqs = [normt, momt]
    for (ctr, rad) in zip(centers, radii):
        if verbose:
            print 'the center :', ctr
        vc = vector(ctr)
        left = vm - vc.cross_product(vt)
        equ = left.dot_product(left) - rad**2
        if verbose:
            print equ
        eqs.append(equ)
    vrs = [x0, x1, x2, x3, x4, x5]
    return eqs, vrs

def quadruples():
    """
    If each pair of the four given spheres touch each other, then the
    tangent lines connect the points where the spheres touch each other.
    In that case, the tangent lines then have multiplicity four.
    This function returns the centers and the radii of the four
    spheres as a tuple of two lists.
    """
    c0 = (+1, +1, +1)
    c1 = (+1, -1, -1)
    c2 = (-1, +1, -1)
    c3 = (-1, -1, +1)
    ctr = [c0, c1, c2, c3]
    rad = [sqrt(2) for _ in range(4)]
    return (ctr, rad)

def quadruples_perturbed():
    """
    If each pair of the four given spheres touch each other, then the
    tangent lines connect the points where the spheres touch each other.
    In that case, the tangent lines then have multiplicity four.
    This function returns the centers and the radii of the four
    spheres as a tuple of two lists.
    """
    c0 = (+1, +1, +1)
    c1 = (+1, -1, -1)
    c2 = (-1, +1, -1)
    c3 = (-1, -1, +1)
    ctr = [c0, c1, c2, c3]
    rad = [sqrt(2.01) for _ in range(4)]
    return (ctr, rad)

def doubles():
    """
    Returns the centers and the radii for a case where the solutions 
    occur with multiplicity two.  The reference for this case is
    the paper by Frank Sottile and Thorsten Theobald:
    "Line problems in nonlinear computational geometry"
    In Computational Geometry - Twenty Years Later, pages 411-432,
    edited by J.E. Goodman, J. Pach, and R. Pollack, AMS, 2008.
    """
    c0 = (2, 2, 0)
    c1 = (2, 0, 2)
    c2 = (0, 2, 2)
    c3 = (0, 0, 0)
    ctr = [c0, c1, c2, c3]
    rad = [3/2 for _ in range(4)]
    return (ctr, rad)

def tostrings(eqs):
    """
    Returns the string representations of the equations
    for input to the solve of phcpy.
    """
    result = []
    for equ in eqs:
        pol = str(equ) + ';'
        result.append(pol)
    return result

def real_coordinates(vars, sols):
    """
    Given in vars the string representations of the variables
    and in sols the list of solutions in string format, returns
    the real parts of the coordinates of the solutions as tuples,
    collected in a list.
    """
    from phcpy.solutions import coordinates
    result = []
    for sol in sols:
        (names, values) = coordinates(sol)
        crd = []
        for var in vars:
            idx = names.index(var)
            nbr = RDF(values[idx].real)
            val = float(nbr.n(digits=12)) # remove round off errors
            crd.append((val if abs(val) > 1.e-12 else 0.0))
        result.append(tuple(crd))
    return result

def iszero(vec, tol):
    """
    Returns True if every component of the vector vec
    has a smaller magnitude than the tolerance tol.
    """
    for x in vec:
        if abs(x) > tol:
            return False
    return True

def tangent_lines(solpts, verbose=True):
    """
    Given the list of solution points in solpts, returns the list of 
    tuples which represent the lines.  Each tuple contains a point on 
    the line and the tangent vector.  Because the tangent t is normalized,
    a point on the line is computed as t x m, where m is the moment vector.
    """
    result = []
    for point in solpts:
        if verbose:
            print point
        tan = vector(point[0:3])
        mom = vector(point[3:6])
        pnt = tan.cross_product(mom) # solves m = p x t
        result.append((pnt, tan))
    return result

def plot_quadruple_spheres(ctr, rad):
    """
    Returns the spheres with centers in ctr and radii in rad
    as golden balls, for the multiplicity four tangent lines.
    """
    x, y, z = var('x, y, z')
    eqs = [(x - c[0])^2 + (y - c[1])^2 + (z - c[2])^2 - r^2 \
        for (c, r) in zip(ctr, rad)]
    xr = (x, -3, 3)
    yr = (y, -3, 3)
    zr = (z, -3, 3)
    balls = [implicit_plot3d(equ, xr, yr, zr, color='gold') for equ in eqs]
    return sum(balls)

def plot_double_spheres(ctr, rad):
    """
    Returns the spheres with centers in ctr and radii in rad
    as golden balls, for the multiplicity two tangent lines.
    """
    x, y, z = var('x, y, z')
    eqs = [(x - c[0])^2 + (y - c[1])^2 + (z - c[2])^2 - r^2 \
        for (c, r) in zip(ctr, rad)]
    xr = (x, -2, 4)
    yr = (y, -2, 4)
    zr = (z, -2, 4)
    balls = [implicit_plot3d(equ, xr, yr, zr, color='gold') for equ in eqs]
    return sum(balls)

def insigned(vectors, tangent):
    """
    Returns True if -tangent occurs in the list vectors,
    returns False otherwise.
    """
    for vec in vectors:
        if ((abs(vec[0] + tangent[0]) < 1.0e-12) and
            (abs(vec[1] + tangent[1]) < 1.0e-12) and
            (abs(vec[2] + tangent[2]) < 1.0e-12)):
            return True
    return False

def filter(points, tangents):
    """
    Every tangent vector t occurs as a pairs: (+t, -t), as +t and -t.
    Given the list of points and corresponding tangent vectors,
    returns a tuple of two lists: a list of points and a list of tangents,
    where every tangent vector occurs only once.
    Therefore, the lists on return should be half the length
    of the lists on input.
    """
    pts = []
    tgs = []
    for (pnt, tan) in zip(points, tangents):
        if not insigned(tgs, tan):
            pts.append(pnt)
            tgs.append(tan)
    return (pts, tgs)

def plot_quadruple_tangents(lines, verbose=True):
    """
    Given in lines the coordinates of the points and the tangents,
    returns the plot of the lines in the range [-3, +3] for all three
    coordinates, in the case of the tangents of multiplicity four.
    """
    points = [pnt for (pnt, tan) in lines]
    tangents = [tan for (pnt, tan) in lines]
    if verbose:
        for (pnt, tan) in zip(points, tangents):
            print 'point :', pnt
            print 'tangent :', tan
    (filpts, filtan) = filter(points, tangents)
    if verbose:
        print 'after filtering :'
        for (pnt, tan) in zip(filpts, filtan):
            print 'point :', pnt
            print 'tangent :', tan
    a1 = vector(RR, filpts[0]) + 3*vector(RR, filtan[0])
    a2 = vector(RR, filpts[1]) + 3*vector(RR, filtan[1])
    a3 = vector(RR, filpts[2]) + 3*vector(RR, filtan[2])
    b1 = vector(RR, filpts[0]) - 3*vector(RR, filtan[0])
    b2 = vector(RR, filpts[1]) - 3*vector(RR, filtan[1])
    b3 = vector(RR, filpts[2]) - 3*vector(RR, filtan[2])
    L1 = line([a1, b1], thickness=10, color='red')
    L2 = line([a2, b2], thickness=10, color='blue')
    L3 = line([a3, b3], thickness=10, color='green')
    fig = L1 + L2 + L3
    return fig

def plot_twelve_tangents(lines, verbose=True):
    """
    Given in lines the coordinates of the points and the tangents,
    returns the plot of the lines in the range [-3, +3] for all three
    coordinates, in the case of the twelve tangents of the
    perturbed case of the multiplicity four case.
    """
    points = [pnt for (pnt, tan) in lines]
    tangents = [tan for (pnt, tan) in lines]
    if verbose:
        for (pnt, tan) in zip(points, tangents):
            print 'point :', pnt
            print 'tangent :', tan
    (filpts, filtan) = filter(points, tangents)
    if verbose:
        print 'after filtering :'
        for (pnt, tan) in zip(filpts, filtan):
            print 'point :', pnt
            print 'tangent :', tan
    first = True
    for (pnt, tan) in zip(filpts, filtan):
        apt = vector(RR, pnt) + 5*vector(RR, tan)
        bpt = vector(RR, pnt) - 5*vector(RR, tan)
        abL = line([apt, bpt], thickness=3, color='red')
        fig = (abL if first else fig + abL)
        first = False
    return fig

def plot_double_tangents(lines, verbose=True):
    """
    Given in lines the coordinates of the points and the tangents,
    plots the lines in the range [-2, 4] for the x, y, z coordinates,
    for the case of the tangents of multiplicity two.
    """
    points = [pnt for (pnt, tan) in lines]
    tangents = [tan for (pnt, tan) in lines]
    if verbose:
        for (pnt, tan) in zip(points, tangents):
            print 'point :', pnt
            print 'tangent :', tan
    (filpts, filtan) = filter(points, tangents)
    if verbose:
        print 'after filtering :'
        for (pnt, tan) in zip(filpts, filtan):
            print 'point :', pnt
            print 'tangent :', tan
    a1 = vector(RR, filpts[0]) + 5*vector(RR, filtan[0])
    a2 = vector(RR, filpts[1]) + 5*vector(RR, filtan[1])
    a3 = vector(RR, filpts[2]) + 5*vector(RR, filtan[2])
    a4 = vector(RR, filpts[3]) + 5*vector(RR, filtan[3])
    a5 = vector(RR, filpts[4]) + 5*vector(RR, filtan[4])
    a6 = vector(RR, filpts[5]) + 5*vector(RR, filtan[5])
    b1 = vector(RR, filpts[0]) - 5*vector(RR, filtan[0])
    b2 = vector(RR, filpts[1]) - 5*vector(RR, filtan[1])
    b3 = vector(RR, filpts[2]) - 5*vector(RR, filtan[2])
    b4 = vector(RR, filpts[3]) - 5*vector(RR, filtan[3])
    b5 = vector(RR, filpts[4]) - 3*vector(RR, filtan[4])
    b6 = vector(RR, filpts[5]) - 5*vector(RR, filtan[5])
    L0 = line([a1, b1], thickness=10, color='red')
    L1 = line([a2, b2], thickness=10, color='blue')
    L2 = line([a3, b3], thickness=10, color='green')
    L3 = line([a4, b4], thickness=10, color='black')
    L4 = line([a5, b5], thickness=10, color='orange')
    L5 = line([a6, b6], thickness=10, color='purple')
    fig = L0 + L1 + L2 + L3 + L4 + L5
    return fig

def main():
    """
    Solves the problem of computing all tangent lines
    to four given spheres.
    """
    ans = raw_input('Double lines ? (y/n) ')
    dbl = (ans == 'y')
    if dbl:
        (ctr, rad) = doubles()
    else:
        ans = raw_input('Perturbed ? (y/n) ')
        prb = (ans == 'y')
        if prb:
            (ctr, rad) = quadruples_perturbed()
        else:
            (ctr, rad) = quadruples()
    ans = raw_input('Verbose ? (y/n) ')
    verbose = (ans == 'y')
    eqs, vrs = tangent_system(ctr, rad, verbose)
    if verbose:
        print 'the polynomial system :'
        for equ in eqs:
            print equ
    pols = tostrings(eqs)
    if verbose:
        for pol in pols:
            print pol
        print 'calling the solver of phcpy :'
    from phcpy.solver import solve
    sols = solve(pols, verbose=False)
    if verbose:
        print 'the solutions :'
        for sol in sols:
            print sol
    vars = ['x0', 'x1', 'x2', 'x3', 'x4', 'x5']
    pts = real_coordinates(vars, sols)
    if verbose:
        print 'the coordinates of the solution points :'
        for point in pts:
            print point
    lines = tangent_lines(pts, verbose)
    if verbose:
        print 'the tangent lines :'
        for line in lines:
            print line
    if dbl:
        fig1 = plot_double_spheres(ctr, rad)
        fig2 = plot_double_tangents(lines, verbose)
        thefig = fig1 + fig2
    else:
        fig1 = plot_quadruple_spheres(ctr, rad)
        if prb:
            fig2 = plot_twelve_tangents(lines, verbose)
            thefig = (fig1 + fig2).rotateZ(-pi/12)
        else:
            fig2 = plot_quadruple_tangents(lines, verbose)
            thefig = fig1 + fig2
    thefig.save('tangents.png')

main()
