"""
This script computes the tangents lines to a circle
via a witness set intersection, with phcpy in Python 3.
With matplotlib a plot is made of the tangent lines.
"""
from math import cos, sin, pi
from random import uniform
import matplotlib.pyplot as plt
from phcpy.solver import solve
from phcpy.solutions import make_solution, coordinates, strsol2dict
from phcpy.sets import double_embed
from phcpy.diagonal import double_diagonal_solve

def polynomials(apt, bpt, rad):
    """
    Returns string representations of two polynomials:
    1) a circle with radius r centered at (apt, bpt);
    2) a line through the origin with slope s.
    """
    crc = f'(x - {apt:.15e})^2 + (y - {bpt:.15e})^2 - {rad**2:.15e};'
    lin = 'y - s*x;'
    return [crc, lin]

def special_solutions(pols, slope):
    """
    Given in pols the polynomials for the line intersecting a circle
    for a general slope s, solves the problem for a special numerical
    value for s, given in the slope.
    The value of the slope is added as last coordinate to each solution.
    """
    special = []
    for pol in pols:
        rpl = pol.replace('s', '(' + str(slope) + ')')
        special.append(rpl)
    sols = solve(special)
    result = []
    for sol in sols:
        (names, values) = coordinates(sol)
        names.append('s')
        values.append(slope)
        extsol = make_solution(names, values)
        result.append(extsol)
    return result

def make_witness_set(pols, verbose=True):
    """
    We have two equations in three variables in pols
    and therefore we expect a one dimensional set.
    """
    embpols = double_embed(3, 1, pols)
    if verbose:
        print('the embedded system :')
        for pol in embpols:
            print(pol)
    embsols = solve(embpols)
    if verbose:
        print('the witness points :')
        for sol in embsols:
            print(sol)
    return (embpols, embsols)

def circle_line_set():
    """
    Generates the system and its witness set.
    Returns the witness set for a fixed circle 
    intersect with a one parameter family of lines 
    as a tuple of polynomials and solutions
    """
    syst = polynomials(3, 2, 1)
    for pol in syst:
        print(pol)
    spsols = special_solutions(syst, 1)
    for sol in spsols:
        print(sol)
    input('hit enter to continue')
    (embsyst, embsols) = make_witness_set(syst, False)
    print('the polynomials in the witness set:')
    for pol in embsyst:
        print(pol)
    print('the solutions :')
    for sol in embsols:
        print(sol)
    print('degree of the set :', len(embsols))
    return (embsyst, embsols)

def random_complex():
    """
    Returns a random complex number on the unit circle.
    """
    theta = uniform(0, 2*pi)
    return complex(cos(theta), sin(theta))

def random_hyperplane(variables):
    """
    Returns a linear equation in the variables,
    with random complex coefficients.
    """
    cf0 = str(random_complex())
    tf0 = cf0.replace('j', '*i')
    result = tf0
    for var in variables:
        cff = str(random_complex())
        tcf = cff.replace('j', '*i')
        result = result + '+' + tcf + '*' + var
    return result + ';'

def jacobian(apt, bpt):
    """
    Returns the equations which define the points
    where the Jacobian matrix is singular,
    as a random linear combination of the columns.
    Random complex coefficients are generated to
    scale the multiplier variables.
    """
    eq1 = f'2*(x-{apt:.15e})*L1 + 2*(y-{bpt:.15e})*L2;'
    eq2 = '-s*L1 + L2;'
    eq3 = random_hyperplane(['L1', 'L2'])
    print('eq3 = ', eq3)
    return [eq1, eq2, eq3]

def extend_solutions(sols):
    """
    To each solution in sols, adds L1 and L2 with values 1,
    and zz2 and zz3 with values zero.
    """
    result = []
    for sol in sols:
        (names, values) = coordinates(sol)
        names = names + ['L1', 'L2', 'zz2', 'zz3']
        values = values + [1, 1, 0, 0]
        extsol = make_solution(names, values)
        result.append(extsol)
    return result

def extend(pols, sols, verbose=True):
    """
    Extends the witness set with two free variables
    L1 and L2, addition two linear equations,
    and two slack variables zz2 and zz3.
    """
    eq1 = 'zz2;'
    eq2 = 'zz3;'
    eq3 = 'L1 - 1;'
    eq4 = 'L2 - 1;'
    extpols = pols[:-1] + [eq1, eq2, eq3, eq4, pols[-1]]
    extsols = extend_solutions(sols)
    if verbose:
        print('the extended polynomials :')
        for pol in extpols:
            print(pol)
        input('hit enter to continue')
        print('the extended solutions :')
        for sol in extsols:
            print(sol)
        input('hit enter to continue')
    return (extpols, extsols)

def witset(pols, verbose=True):
    """
    We have three equations in pols in five variables:
    x, y, s, L1, and L2.
    Therefore we expect a two dimensional set.
    """
    embpols = double_embed(5, 2, pols)
    if verbose:
        print('the embedded system :')
        for pol in embpols:
            print(pol)
    embsols = solve(embpols)
    if verbose:
        print('the witness points :')
        for sol in embsols:
            print(sol)
    return (embpols, embsols)

def insert_symbols(pol):
    """
    To the string pol, adds the sequence of symbols.
    """
    qpol = pol.lstrip()
    if qpol[0] == '+' or qpol[0] == '-':
        smb = 'x - x + y - y + s + L1 - L1 + L2 - L2 - s '
    else:
        smb = 'x - x + y - y + s + L1 - L1 + L2 - L2 - s + '
    return smb + pol

def intersect(dim, w1d, w2d, ws1, ws2):
    """
    Applies the diagonal homotopy to intersect two witness sets
    w1 and w2 of dimensions w1d and w2d in a space of dimension dim.
    """
    w1eqs, w1sols = ws1
    w2eqs, w2sols = ws2
    nw1eq0 = insert_symbols(w1eqs[0])
    nw2eq0 = insert_symbols(w2eqs[0])
    nw1eqs = [nw1eq0] + w1eqs[1:]
    nw2eqs = [nw2eq0] + w2eqs[1:]
    print('number of equations in first witness set :', len(w1eqs))
    print('number of equations in second witness set :', len(w2eqs))
    ans = input('Verbose version of diagonal solve ? (y/n) ')
    otp = int(ans == 'y')
    result = double_diagonal_solve(dim, w1d, nw1eqs, w1sols, w2d, nw2eqs, w2sols, \
        vrblvl=otp)
    (eqs, sols) = result
    print('the equations :')
    for pol in eqs:
        print(pol)
    input('hit enter to continue')
    print('the solutions :')
    for sol in sols:
        print(sol)
    print('computed', len(sols), 'solutions')
    input('hit enter to continue')
    return result

def singular_locus_set():
    """
    Generates a witness set for the singular locus
    of the algebraic set of a fixed circle intersected
    with a one parameter family of lines.
    """
    syst = jacobian(3, 2)
    for pol in syst:
        print(pol)
    (embsyst, embsols) = witset(syst, False)
    print('the polynomials in the witness set :')
    for pol in embsyst:
        print(pol)
    input('hit enter to continue')
    print('the solutions :')
    for sol in embsols:
        print(sol)
    print('degree of the singular locus set :', len(embsols))
    input('hit enter to continue')
    return (embsyst, embsols)

def coordinates_and_slopes(sol):
    """
    Given a solution, return the 3-tuple with the x and y coordinates
    of the tangent point and the slope of the tangent line.
    The real parts of the coordinates are selected.
    """
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
    (xp1, yp1, slope1) = sol1
    (xp2, yp2, slope2) = sol2
    axs = fig.gca()
    center = (3, 2)
    radius = 1
    circle = plt.Circle(center, radius, edgecolor='blue', \
       facecolor='none')
    axs.add_artist(circle)
    y1r = 5*slope1
    y2r = 5*slope2
    plt.plot([0, 5], [0, y1r], 'r')
    plt.plot([0, 5], [0, y2r], 'r')
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
    input('hit enter to exit')

def main():
    """
    Defines two witness sets and intersects them.
    """
    (witpols, witsols) = circle_line_set()
    input('hit enter to continue')
    witset1 = extend(witpols, witsols)
    witset2 = singular_locus_set()
    intwitset = intersect(5, 3, 2, witset1, witset2)
    (eqs, sols) = intwitset
    print('the equations :')
    for equ in eqs:
        print(equ)
    print('the solutions :')
    for (idx, sol) in enumerate(sols):
        print('Solution', idx+1, ':')
        print(sol)
    sol1 = coordinates_and_slopes(sols[0])
    sol2 = coordinates_and_slopes(sols[1])
    show_solutions(sol1, sol2)

if __name__ == "__main__":
    main()
