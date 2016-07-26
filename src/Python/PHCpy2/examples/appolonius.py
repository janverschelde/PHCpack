"""
The circle problem of Appolonius has the following input/output specification:
given three circles, find all circles that are tangent to the given circles.
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
syst = [eqs0,eqs1,eqs2,eqs3,eqs4,eqs5,eqs6,eqs7]
from phcpy.solver import solve
from phcpy.solutions import strsol2dict
(circle, eqscnt) = (0, 0)
for eqs in syst:
    sols = solve(eqs, silent=True)
    eqscnt = eqscnt + 1
    print 'system', eqscnt, 'has', len(sols), 'solutions'
    for sol in sols:
        soldic = strsol2dict(sol)
        if soldic['r'].real > 0:
            circle = circle + 1
            print 'solution circle', circle
            print 'x =', soldic['x']
            print 'y =', soldic['y']
            print 'r =', soldic['r']
