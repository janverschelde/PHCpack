"""
Illustration of a mixed volume computation, polyhedral
homotopies and the path tracking in double precision.
First we construct a random coefficient start system,
computing the mixed volume and running polyhedral homotopies.
Then we use the constructed start system g and corresponding
start solutions sols to solve the original system f.
Finally, we verify the solutions with one Newton step
on each computed solution, in double double precision.
"""
f = ['x^3*y + x*y^2 - 8;', 'x*y - 3;']
print 'the polynomials in the target system :'
for pol in f:
    print pol
from phcpy.solver import mixed_volume as mv
from phcpy.solver import random_coefficient_system as rcs
print 'the mixed volume :', mv(f)
(g, gsols) = rcs()
print 'the polynomials in the start system :'
for pol in g:
    print pol
print 'number of start solutions :', len(gsols)
from phcpy.trackers import standard_double_track as track
fsols = track(f, g, gsols)
from phcpy.solver import newton_step
for sol in fsols:
    print sol
print 'one Newton step in double double precision...'
nsols = newton_step(f, fsols, precision='dd')
for sol in nsols:
    print sol
