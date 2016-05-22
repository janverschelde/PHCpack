"""
Illustration of the witness set computation of the cyclic 4-roots system.
"""
from phcpy.families import cyclic
c4 = cyclic(4)
from phcpy.sets import embed
c4e1 = embed(4, 1, c4)
print 'the embedded cyclic 4-roots problem :'
for pol in c4e1:
    print pol
from phcpy.solver import solve
sols = solve(c4e1)
print 'computed', len(sols), 'solutions'
from phcpy.solutions import filter_zero_coordinates as filter
genpts = filter(sols, 'zz1', 1.0e-8, 'select')
print 'generic points :'
for sol in genpts:
    print sol
from phcpy.sets import membertest
sdpoint = [-1, 0, -1, 0, 1, 0, 1, 0]
print 'testing in standard double precision ...'
print membertest(c4e1, genpts, 1, sdpoint, verbose=True, precision='d')
raw_input('*** hit enter to continue ***')
print 'testing in double double precision ...'
ddpoint = [-1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]
print membertest(c4e1, genpts, 1, ddpoint, verbose=True, precision='dd')
raw_input('*** hit enter to continue ***')
print 'testing in quad double precision ...'
ddpoint = [-1, 0, 0, 0, 0, 0, 0, 0, \
           -1, 0, 0, 0, 0, 0, 0, 0, \
            1, 0, 0, 0, 0, 0, 0, 0, \
            1, 0, 0, 0, 0, 0, 0, 0]
print membertest(c4e1, genpts, 1, ddpoint, verbose=True, precision='qd')
