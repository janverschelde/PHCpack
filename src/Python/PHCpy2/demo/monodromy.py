"""
As an extension of the witness set for the twisted cubic,
we consider the decomposition of the twisted cubic and line.
The decomposition is computed via monodromy on the witness
set for the one dimensional solution set.
"""
# f = ['(x^2 - y)*(y-0.5);', '(x^3 - z)*(z-0.5);']
f = ['x^2 - y;', 'x^3 - z;']
print 'polynomials that define an algebraic set:'
for pol in f:
    print pol
from phcpy.sets import embed, cascade_step
from phcpy.solver import solve
ef = embed(3, 1, f)
# patch : make sure zz1 is last symbol!
ef[0] = 'x - x + y - y + z - z + ' + ef[0]
efsols = solve(ef)
print 'degree of the set :', len(efsols)
from phcpy.sets import monodromy_breakup
monodromy_breakup(ef, efsols, 1)
