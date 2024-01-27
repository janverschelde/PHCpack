"""
Counter example of Bertrand Haas, against the Koushnirenko conjecture,
executed on one core and on many cores.

Bertrand Haas: A simple counterexample to Kouchnirenko's conjecture.
Beitraege zur Algebra und Geometrie/Contributions to Algebra and Geometry, 
volume 43, number 1, pages 1 to 8, 2002.

Bertrand constructed a polynomial system with three monomials in every
equation and with five positive real roots.  Kouchnirenko's conjecture
predicted there could only four positive real roots for such a system.
"""
from datetime import datetime
from phcpy.dimension import get_core_count
from phcpy.solver import solve

H = [ 'x**108 + 1.1*y**54 - 1.1*y;', 
      'y**108 + 1.1*x**54 - 1.1*x;' ]
print('Solving on one core ...')
wstart = datetime.now()
sols = solve(H)
wstop = datetime.now()
print('  Number of solutions :', len(sols))
print('start time :', wstart)
print(' stop time :', wstop)
print('   elapsed :', wstop - wstart)
nbcores = get_core_count()
print('Solving on', nbcores, 'cores ...')
wstart = datetime.now()
sols = solve(H, tasks=nbcores)
wstop = datetime.now()
print('  Number of solutions :', len(sols))
print('start time :', wstart)
print(' stop time :', wstop)
print('   elapsed :', wstop - wstart)
