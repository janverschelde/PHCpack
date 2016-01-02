"""
Illustration of equation and variable scaling.
On a simple example with coefficients of different magnitudes,
we see that after scaling, the coefficients differ less in magnitude.
Morever, the local condition numbers of the solution improve from
10^4 to 10, thus gaining about three decimal places in accuracy.
"""
from phcpy.solver import solve
from phcpy.solver import standard_scale_system as scale_system
from phcpy.solver import standard_scale_solutions as scale_solutions
from phcpy.phcpy2c import py2c_set_seed as set_seed
set_seed(2015)
p = ['0.000001*x^2 + 0.000004*y^2 - 4;', '0.000002*y^2 - 0.001*x;']
psols = solve(p, silent=True) 
print 'the first solution of the original system :'
print psols[0]
(q, c) = scale_system(p)
print 'the scaling coefficients :', c
print 'the scaled system :'
for pol in q:
    print pol
qsols = solve(q, silent=True)
print 'the corresponding solution of the scaled system :'
print qsols[3]
ssols = scale_solutions(len(q), qsols, c)
print 'the corresponding solution in original coordinates :'
print ssols[3]
