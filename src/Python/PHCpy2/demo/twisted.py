"""
Making a witness set of the twisted cubic.
The twisted cubic is a curve in 3-space of degree 3,
defined as the intersection of a quadratic and a cubic cylinder.
Adding a random plane to the original system leads to three
generic points on the cubic.  These three points, jointly
with the embedded system form a witness set for the twisted cubic.
"""
twisted = ['x^2 - y;', 'x^3 - z;']
print 'polynomials that define the twisted cubic :'
for pol in twisted:
    print pol
from phcpy.sets import embed, cascade_step
from phcpy.solver import solve
et = embed(3, 1, twisted)
print 'polynomials in the embedded system :'
for pol in et:
    print pol
etsols = solve(et)
print 'the witness points :'
for point in etsols:
    print point
