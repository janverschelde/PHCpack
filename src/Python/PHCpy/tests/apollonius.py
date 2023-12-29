"""
The circle problem of Apollonius has the following input/output specification:
Given three circles, find all circles that are tangent to the given circles.
Without loss of generality, we take the first circle to be the unit circle,
centered at (0, 0) and with radius 1.  The origin of the second circle lies
on the first coordinate axis, so its center has coordinates (c2x, 0) and
radius r2.  The third circle has center (c3x, c3y) and radius r3.
So there are five parameters in this problem: c2x, r2, c3x, c3y, and r3.
Values for the five parameters are defined by the first five equations.
The next three equations determine the center (x, y) and the radius r
of the circle which touches the three given circles.
The condition on the center of the touching circle is that its distance
to the center of the given circle is either the difference or the sum of
the radii of both circles.  So we arrive at eight polynomial systems.
This Python 3 script solves the systems with phcpy and then prints those
solution circles that have a positive radius.
With matplotlib, the given circles and the solutions are plotted.
"""
def polynomials(c2x, r2, c3x, c3y, r3):
    """
    On input are the five parameters of the circle problem of Apollonius:
    c2x : the x-coordinate of the center of the second circle,
    r2 : the radius of the second circle,
    c3x : the x-coordinate of the center of the third circle,
    c3y : the y-coordinate of the center of the third circle,
    r3 : the radius of the third circle.
    Returns a list of lists.  Each list contains a polynomial system.
    Solutions to each polynomial system define center (x, y) and radius r
    of a circle touching three given circles.
    """
    e1m = 'x^2 + y^2 - (r-1)^2;'
    e1p = 'x^2 + y^2 - (r+1)^2;'
    e2m = '(x-%.15f)^2 + y^2 - (r-%.15f)^2;' % (c2x, r2)
    e2p = '(x-%.15f)^2 + y^2 - (r+%.15f)^2;' % (c2x, r2)
    e3m = '(x-%.15f)^2 + (y-%.15f)^2 - (r-%.15f)^2;' % (c3x, c3y, r3)
    e3p = '(x-%.15f)^2 + (y-%.15f)^2 - (r+%.15f)^2;' % (c3x, c3y, r3)
    eqs0 = [e1m,e2m,e3m]
    eqs1 = [e1m,e2m,e3p]
    eqs2 = [e1m,e2p,e3m]
    eqs3 = [e1m,e2p,e3p]
    eqs4 = [e1p,e2m,e3m]
    eqs5 = [e1p,e2m,e3p]
    eqs6 = [e1p,e2p,e3m]
    eqs7 = [e1p,e2p,e3p]
    return [eqs0,eqs1,eqs2,eqs3,eqs4,eqs5,eqs6,eqs7]

def solve4circles(syst, verbose=True):
    """
    Given in syst is a list of polynomial systems.
    Returns a list of tuples.  Each tuple in the list of return
    consists of the coordinates of the center and the radius of
    a circle touching the three given circles.
    """
    from phcpy.solver import solve
    from phcpy.solutions import strsol2dict, is_real
    (circle, eqscnt) = (0, 0)
    result = []
    for eqs in syst:
        eqscnt = eqscnt + 1
        if verbose:
            print('solving system', eqscnt, ':')
            for pol in eqs:
                print(pol)
        sols = solve(eqs)
        if verbose:
            print('system', eqscnt, 'has', len(sols), 'solutions')
        for sol in sols:
            if is_real(sol, 1.0e-8):
                soldic = strsol2dict(sol)
                if soldic['r'].real > 0:
                    circle = circle + 1
                    ctr = (soldic['x'].real, soldic['y'].real)
                    rad = soldic['r'].real
                    result.append((ctr, rad))
                    if verbose:
                        print('solution circle', circle)
                        print('center =', ctr)
                        print('radius =', rad)
    return result

def makecircles(plt, data, color, disk=False, verbose=True):
    """
    Gives in data a list of centers and radii,
    returns a list of matplotlib objects, using the color as edgecolor,
    The pyplot in matplotlib is passed as plt.
    """
    result = []
    for (center, radius) in data:
        if verbose:
            print('circle with center', center, 'and radius', radius)
        if disk:
            crc = plt.Circle(center, radius, edgecolor=color, \
                 facecolor=color) #, linewidth=2)
        else:
            crc = plt.Circle(center, radius, edgecolor=color, \
                 facecolor='none') # , linewidth=2)
        result.append(crc)
    return result

def plotcircles(plt, circles, xa, xb, ya, yb):
    """
    Given on input in circles a list of matplotlib objects,
    the circles are rendered in a figure.
    The pyplot in matplotlib is passed as plt.
    The range of the axis in the matplotlib plot is defined
    by the four numbers xa, xb, ya, and yb.
    """
    fig = plt.figure()
    fig.add_subplot(111, aspect='equal')
    axs = plt.gcf().gca()
    for circle in circles:
        axs.add_artist(circle)
    plt.axis([xa, xb, ya, yb])
    fig.show() # plt.show()
    ans = input('hit enter to continue')

def solve_general_problem():
    """
    Solves a general configuration of three circles.
    Defines a list of polynomial systems, solves each system
    and extract those real solutions with positive radius.
    The given circles are plotted as blue disks,
    while the eight solution circles are plotted in red.
    """
    syst = polynomials(2, 2.0/3, 1, 1, 1.0/3)
    sols = solve4circles(syst)
    print('the solution list :')
    print(sols)
    ans = input('Continue with matplotlib ? (y/n) ')
    if ans == 'y':
        import matplotlib.pyplot as plt
        crcdata = [((0, 0), 1), ((2, 0), 2.0/3), ((1, 1), 1.0/3)]
        incircles = makecircles(plt, crcdata, disk=True, color='blue')
        outcircles = makecircles(plt, sols, color='red')
        plotcircles(plt, incircles + outcircles, -3, 11, -3, 8)

def solve_special_problem():
    """
    Solves a special configuration of three circles,
    where the three circles are mutually touching each other.
    Defines a list of polynomial systems, solves each system
    and extract those real solutions with positive radius.
    The given circles are plotted as blue disks,
    while the solution circles are plotted in red.
    """
    from math import sqrt
    height = sqrt(3)
    syst = polynomials(2, 1, 1, height, 1)
    sols = solve4circles(syst)
    print('the solution list :')
    print(sols)
    ans = input('Continue with matplotlib ? (y/n) ')
    if ans == 'y':
        import matplotlib.pyplot as plt
        crcdata = [((0, 0), 1), ((2, 0), 1), ((1, height), 1)]
        incircles = makecircles(plt, crcdata, disk=True, color='blue')
        outcircles = makecircles(plt, sols, color='red')
        plotcircles(plt, incircles + outcircles, -2, 4, -2, 3)

def solve_perturbed_problem():
    """
    Solves a small perturbation of a special configuration of three circles,
    where the three circles are mutually touching each other.
    Defines a list of polynomial systems, solves each system
    and extract those real solutions with positive radius.
    The given circles are plotted as blue disks,
    while the solution circles are plotted in red.
    """
    from math import sqrt
    height = sqrt(3)
    syst = polynomials(2.05, 1, 1.025, height+0.025, 1)
    sols = solve4circles(syst)
    print('the solution list :')
    print(sols)
    ans = input('Continue with matplotlib ? (y/n) ')
    if ans == 'y':
        import matplotlib.pyplot as plt
        crcdata = [((0, 0), 1), ((2.05, 0), 1), ((1.025, height+0.025), 1)]
        incircles = makecircles(plt, crcdata, disk=True, color='blue')
        outcircles = makecircles(plt, sols, color='red')
        plotcircles(plt, incircles + outcircles, -2, 4, -2, 4)

def main():
    """
    Solves a general and a special instance of the circle problem
    of Apollonius.
    """
    import phcpy
    print('solving a general instance of the Apollonius circle problem')
    solve_general_problem()
    print('solving a special instance of the Apollonius circle problem')
    solve_special_problem()
    print('solving a perturbed instance of the Apollonius circle problem')
    solve_perturbed_problem()

if __name__=="__main__":
    main()
