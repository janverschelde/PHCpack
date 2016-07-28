"""
The circle problem of Appolonius has the following input/output specification:
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
This Python 2 script solves the systems with phcpy and then prints those
solution circles that have a positive radius.
With matplotlib, the given circles and the solutions are plotted.
"""
def polynomials():
    """
    Returns a list of lists.  Each list contains a polynomial system.
    Solutions to each polynomial system define center (x, y) and radius r
    of a circle touching three given circles.
    """
    e1 = 'c2x - 2;'
    e2 = 'r2 - 2/3;'
    e3 = 'c3x - 1;'
    e4 = 'c3y - 1;'
    e5 = 'r3 - 1/3;'
    e6m = 'x^2 + y^2 - (r-1)^2;'
    e6p = 'x^2 + y^2 - (r+1)^2;'
    e7m = '(x-c2x)^2 + y^2 - (r-r2)^2;'
    e7p = '(x-c2x)^2 + y^2 - (r+r2)^2;'
    e8m = '(x-c3x)^2 + (y-c3y)^2 - (r-r3)^2;'
    e8p = '(x-c3x)^2 + (y-c3y)^2 - (r+r3)^2;'
    eqs0 = [e1,e2,e3,e4,e5,e6m,e7m,e8m]
    eqs1 = [e1,e2,e3,e4,e5,e6m,e7m,e8p]
    eqs2 = [e1,e2,e3,e4,e5,e6m,e7p,e8m]
    eqs3 = [e1,e2,e3,e4,e5,e6m,e7p,e8p]
    eqs4 = [e1,e2,e3,e4,e5,e6p,e7m,e8m]
    eqs5 = [e1,e2,e3,e4,e5,e6p,e7m,e8p]
    eqs6 = [e1,e2,e3,e4,e5,e6p,e7p,e8m]
    eqs7 = [e1,e2,e3,e4,e5,e6p,e7p,e8p]
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
        sols = solve(eqs, silent=True)
        eqscnt = eqscnt + 1
        if verbose:
            print 'system', eqscnt, 'has', len(sols), 'solutions'
        for sol in sols:
            if is_real(sol, 1.0e-8):
                soldic = strsol2dict(sol)
                if soldic['r'].real > 0:
                    circle = circle + 1
                    ctr = (soldic['x'].real, soldic['y'].real)
                    rad = soldic['r'].real
                    result.append((ctr, rad))
                    if verbose:
                        print 'solution circle', circle
                        print 'center =', ctr
                        print 'radius =', rad
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
            print 'circle with center', center, 'and radius', radius
        if disk:
            crc = plt.Circle(center, radius, edgecolor=color, \
                 facecolor=color) #, linewidth=2)
        else:
            crc = plt.Circle(center, radius, edgecolor=color, \
                 facecolor='none') # , linewidth=2)
        result.append(crc)
    return result

def plotcircles(plt, circles):
    """
    Given on input in circles a list of matplotlib objects,
    the circles are rendered in a figure.
    The pyplot in matplotlib is passed as plt.
    """
    fig = plt.figure()
    fig.add_subplot(111, aspect='equal')
    axs = plt.gcf().gca()
    for circle in circles:
        axs.add_artist(circle)
    plt.axis([-3, 11, -3, 8])
    fig.show()
    ans = raw_input('hit enter to exit')

def main():
    """
    Defines a list of polynomial systems, solves each system
    and extract those real solutions with positive radius.
    The given circles are plotted as blue disks,
    while the eight solution circles are plotted in red.
    """
    syst = polynomials()
    sols = solve4circles(syst)
    print 'the solution list :'
    print sols
    ans = raw_input('Continue with matplotlib ? (y/n) ')
    if ans == 'y':
        import matplotlib.pyplot as plt
        crcdata = [((0, 0), 1), ((2, 0), 2.0/3), ((1, 1), 1.0/3)]
        incircles = makecircles(plt, crcdata, disk=True, color='blue')
        outcircles = makecircles(plt, sols, color='red')
        plotcircles(plt, incircles + outcircles)

if __name__=="__main__":
    main()
