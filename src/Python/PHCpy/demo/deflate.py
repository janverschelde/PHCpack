"""
We illustrate deflation on the example of Griewank and Osborne.
The point (0, 0) is an irregular singularity and Newton fails
to converge even if we start already quite close to (0, 0).
Deflation solves this problem.
"""
p = ['(29/16)*x^3 - 2*x*y;', 'x^2 - y;']
from phcpy.solutions import make_solution
sol = make_solution(['x', 'y'],[1.0e-6, 1.0e-6])
print 'the initial solution :'
print sol
from phcpy.solver import newton_step
from phcpy.solutions import diagnostics
sols = [sol]
for k in range(5):
    sols = newton_step(p, sols)
    (err, rco, res) = diagnostics(sols[0])
    print 'forward error :', err
    print 'estimate for inverse condition :', rco
    print 'backward error :', res
print 'the solution after five Newton steps :'
print sols[0]
from phcpy.solver import standard_deflate
sols = standard_deflate(p, sols)
print 'after deflation :'
print sols[0]
