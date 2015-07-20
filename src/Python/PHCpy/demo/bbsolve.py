"""
Demonstration of the multitasked blackbox solver.
Run this at the command line as time python d1.py.
"""
from phcpy.families import cyclic
from phcpy.solver import solve
c7 = cyclic(7)
s = solve(c7, tasks=4)
print 'number of solutions :', len(s)
