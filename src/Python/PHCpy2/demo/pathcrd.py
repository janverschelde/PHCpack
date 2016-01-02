"""
Get the x and y coordinates of a solution path.
First we take two quadrics and construct a start system
based on the total degree applying the theorem of Bezout.
To extract the coordinates of the points on a solution path,
we convert the string representation into a dictionary.
"""
f = ['x*y + y^2 - 8;', 'x^2 - 3*y + 4;']
print 'the polynomials in the target system :'
for pol in f:
    print pol
from phcpy.solver import total_degree_start_system as tds
(g, gsols) = tds(f)
print 'the polynomials in the start system :'
for pol in g:
    print pol
print 'number of start solutions :', len(gsols)
from phcpy.trackers import initialize_standard_tracker
from phcpy.trackers import initialize_standard_solution
from phcpy.trackers import next_standard_solution as nxtsol
initialize_standard_tracker(f, g)
initialize_standard_solution(len(g), gsols[0])
points = [nxtsol() for k in range(10)]
from phcpy.solutions import strsol2dict
dicpts = [strsol2dict(sol) for sol in points]
coords = [(sol['x'].real, sol['y'].real) for sol in dicpts]
print 'the real parts of the first ten points on the path :'
for point in coords:
    print point
