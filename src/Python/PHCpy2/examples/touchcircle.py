"""
This script computes the tangents lines to a circle
via a witness set intersection, with phcpy in Python 2.
With matplotlib a plot is made of the tangent lines.
"""
import matplotlib.pyplot as plt

def polynomials(a, b, r):
    """
    Returns string representations of two polynomials:
    1) a circle with radius r centered at (a, b);
    2) a line through the origin with slope s.
    """
    crc = '(x - %.15e)^2 + (y - %.15e)^2 - %.15e;' % (a, b, r**2)
    lin = 'y - s*x;'
    return [crc, lin]

def special_solutions(pols, slope):
    """
    Given in pols the polynomials for the line intersecting a circle
    for a general slope s, solves the problem for a special numerical
    value for s, given in the slope.
    The value of the slope is added as last coordinate to each solution.
    """
    from phcpy.solver import solve
    from phcpy.solutions import make_solution, coordinates
    special = []
    for pol in pols:
        rpl = pol.replace('s', '(' + str(slope) + ')')
        special.append(rpl)
    sols = solve(special, verbose=False)
    result = []
    for sol in sols:
        (vars, vals) = coordinates(sol)
        vars.append('s')
        vals.append(slope)
        extsol = make_solution(vars, vals)
        result.append(extsol)
    return result

def make_witness_set(pols, verbose=True):
    """
    We have two equations in three variables in pols
    and therefore we expect a one dimensional set.
    """
    from phcpy.sets import embed
    from phcpy.solver import solve
    embpols = embed(3, 1, pols)
    if verbose:
        print 'the embedded system :'
        for pol in embpols:
            print pol
    embsols = solve(embpols, verbose=verbose)
    if verbose:
        print 'the witness points :'
        for sol in embsols:
            print sol
    return (embpols, embsols)

def membership_test(witsys, witsols, verbose=True):
    """
    Given a witness sets in the tuple witset,
    runs a membership test on a solution.
    """
    from phcpy.solutions import make_solution
    from phcpy.sets import is_member
    point = make_solution(['x', 'y', 's'], [2, 2, 1])
    print 'testing the point\n', point
    ismb = is_member(witsys, witsols, 1, point, verbose=False)
    if ismb:
        print 'the point is a member'
    else:
        print 'the point is NOT a member'
    return ismb

def circle_line_set():
    """
    Generates the system and its witness set.
    As a sanity check, a membership test is done.
    Returns the witness set for a fixed circle 
    intersect with a one parameter family of lines 
    as a tuple of polynomials and solutions
    """
    syst = polynomials(3, 2, 1)
    for pol in syst:
        print pol
    spsols = special_solutions(syst, 1)
    for sol in spsols:
        print sol
    raw_input('hit enter to continue')
    (embsyst, embsols) = make_witness_set(syst, False)
    print 'the polynomials in the witness set:'
    for pol in embsyst:
        print pol
    print 'the solutions :'
    for sol in embsols:
        print sol
    print 'degree of the set :', len(embsols)
    membership_test(embsyst, embsols)
    return (embsyst, embsols)

def random_complex():
    """
    Returns a random complex number on the unit circle.
    """
    from math import cos, sin, pi
    from random import uniform
    theta = uniform(0, 2*pi)
    return complex(cos(theta), sin(theta))

def random_hyperplane(vars):
    """
    Returns a linear equation in the variables in
    the list vars, with random complex coefficients.
    """
    cf0 = str(random_complex())
    tf0 = cf0.replace('j', '*i')
    result = tf0
    for var in vars:
        cff = str(random_complex())
        tcf = cff.replace('j', '*i')
        result = result + '+' + tcf + '*' + var
    return result + ';'

def jacobian(a, b):
    """
    For the circle centered at (a, b),
    returns the equations which define the points
    where the Jacobian matrix is singular,
    as a random linear combination of the columns.
    Random complex coefficients are generated to
    scale the multiplier variables.
    """
    eq1 = '2*(x-%.15e)*L1 + 2*(y-%.15e)*L2;' % (a, b)
    eq2 = '-s*L1 + L2;'
    eq3 = random_hyperplane(['L1', 'L2'])
    return [eq1, eq2, eq3]

def extend_solutions(sols):
    """
    To each solution in sols, adds L1 and L2 with values 1,
    and zz2 and zz3 with values zero.
    """
    from phcpy.solutions import make_solution, coordinates
    result = []
    for sol in sols:
        (vars, vals) = coordinates(sol)
        vars = vars + ['L1', 'L2', 'zz2', 'zz3']
        vals = vals + [1, 1, 0, 0]
        extsol = make_solution(vars, vals)
        result.append(extsol)
    return result

def extend(pols, sols, verbose=True):
    """
    Extends the witness set with two free variables
    L1 and L2, addition two linear equations,
    and two slack variables zz2 and zz3.
    """
    vars = ['zz2', 'zz3']
    eq1 = 'zz2;'
    eq2 = 'zz3;'
    eq3 = 'L1 - 1;'
    eq4 = 'L2 - 1;'
    extpols = pols[:-1] + [eq1, eq2, eq3, eq4, pols[-1]]
    extsols = extend_solutions(sols)
    if verbose:
        print 'the extended polynomials :'
        for pol in extpols:
            print pol
        raw_input('hit enter to continue')
        print 'the extended solutions :'
        for sol in extsols:
            print sol
        raw_input('hit enter to continue')
    return (extpols, extsols)

def witset(pols, verbose=True):
    """
    We have three equations in pols in five variables:
    x, y, s, L1, and L2.
    Therefore we expect a two dimensional set.
    """
    from phcpy.sets import embed
    from phcpy.solver import solve
    embpols = embed(5, 2, pols)
    if verbose:
        print 'the embedded system :'
        for pol in embpols:
            print pol
    embsols = solve(embpols, verbose=verbose)
    if verbose:
	print 'the witness points :'
        for sol in embsols:
            print sol
    return (embpols, embsols)

def insert_symbols(pol):
    """
    To the string pol, adds the sequence of symbols.
    """
    q = pol.lstrip()
    if q[0] == '+' or q[0] == '-':
        smb = 'x - x + y - y + s + L1 - L1 + L2 - L2 - s '
    else:
        smb = 'x - x + y - y + s + L1 - L1 + L2 - L2 - s + '
    return smb + pol

def intersect(dim, w1d, w2d, ws1, ws2):
    """
    Applies the diagonal homotopy to intersect two witness sets
    w1 and w2 of dimensions w1d and w2d in a space of dimension dim.
    """
    from phcpy.diagonal import diagonal_solver as diagsolve
    w1eqs, w1sols = ws1
    w2eqs, w2sols = ws2
    nw1eq0 = insert_symbols(w1eqs[0])
    nw2eq0 = insert_symbols(w2eqs[0])
    nw1eqs = [nw1eq0] + w1eqs[1:]
    nw2eqs = [nw2eq0] + w2eqs[1:]
    print 'number of equations in first witness set :', len(w1eqs)
    print 'number of equations in second witness set :', len(w2eqs)
    ans = raw_input('Verbose version of diagonal solve ? (y/n) ')
    otp = (ans == 'y')
    result = diagsolve(dim, w1d, nw1eqs, w1sols, w2d, nw2eqs, w2sols, \
        verbose=otp)
    (eqs, sols) = result
    print 'the equations :'
    for pol in eqs:
        print pol
    raw_input('hit enter to continue')
    print 'the solutions :'
    for sol in sols:
        print sol
    print 'computed', len(sols), 'solutions'
    raw_input('hit enter to continue')
    return result

def singular_locus_set():
    """
    Generates a witness set for the singular locus
    of the algebraic set of a fixed circle intersected
    with a one parameter family of lines.
    """
    syst = jacobian(3, 2)
    for pol in syst:
        print pol
    (embsyst, embsols) = witset(syst, False)
    print 'the polynomials in the witness set:'
    for pol in embsyst:
        print pol
    raw_input('hit enter to continue')
    print 'the solutions :'
    for sol in embsols:
        print sol
    print 'degree of the singular locus set :', len(embsols)
    raw_input('hit enter to continue')
    return (embsyst, embsols)

def coordinates_and_slopes(sol):
    """
    Given a solution, return the 3-tuple with the x and y coordinates
    of the tangent point and the slope of the tangent line.
    The real parts of the coordinates are selected.
    """
    from phcpy.solutions import strsol2dict
    sdc = strsol2dict(sol)
    return (sdc['x'].real, sdc['y'].real, sdc['s'].real)

def plotgeneral(fig):
    """
    Plots the circle with radius 1, centered at (3, 2),
    along with a line with slope one through (0, 0).
    """
    axs = fig.gca()
    center = (3, 2)
    radius = 1
    circle = plt.Circle(center, radius, edgecolor='blue', \
       facecolor='none')
    axs.add_artist(circle)
    plt.plot([0, 4], [0, 4], 'r') 
    plt.plot([2, 3], [2, 3], 'go')
    plt.axis([0, 5, 0, 4])

def plotspecial(fig, sol1, sol2):
    """
    Plots the circle with radius 1, centered at (3, 2),
    along with its two lines tangent through (0, 0).
    The two solutions sol1 and sol2 define 3-tuples of
    x and y coordinates and a slope.
    """
    (xp1, yp1, s1) = sol1
    (xp2, yp2, s2) = sol2
    axs = fig.gca()
    center = (3, 2)
    radius = 1
    circle = plt.Circle(center, radius, edgecolor='blue', \
       facecolor='none')
    axs.add_artist(circle)
    y1 = 5*s1
    y2 = 5*s2
    plt.plot([0, 5], [0, y1], 'r') 
    plt.plot([0, 5], [0, y2], 'r') 
    plt.plot([xp1, xp2], [yp1, yp2], 'go')
    plt.axis([0, 5, 0, 4])

def show_solutions(sol1, sol2):
    """
    Makes a matplotlib figure with a general configuration
    and the tangent lines defined by sol1 and sol2.
    """
    fig = plt.figure()
    fig.add_subplot(121, aspect='equal')
    plotgeneral(fig)
    fig.add_subplot(122, aspect='equal')
    plotspecial(fig, sol1, sol2)
    fig.show()
    raw_input('hit enter to exit')

def main():
    """
    Defines two witness sets and intersects them.
    """
    (witpols, witsols) = circle_line_set()
    raw_input('hit enter to continue')
    witset1 = extend(witpols, witsols)
    witset2 = singular_locus_set()
    intwitset = intersect(5, 3, 2, witset1, witset2)
    (eqs, sols) = intwitset
    print 'the solutions :'
    for sol in sols:
        print sol
    sol1 = coordinates_and_slopes(sols[0])
    sol2 = coordinates_and_slopes(sols[1])
    show_solutions(sol1, sol2)

if __name__ == "__main__":
    main()
