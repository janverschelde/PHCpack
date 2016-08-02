the design of a 4-bar mechanism
===============================

Given two pivot points and five precision points for the coupler,
the design problem asks to determine the lengths of the bars that
allow the coupler to pass through the given precision points.

This chapter presents a *use case* for phcpy.
The equations are generated with ``sympy``.

The system is taken from a paper by A.P. Morgan and C.W. Wampler
on **Solving a Planar Four-Bar Design Using Continuation**, published in
the *Journal of Mechanical Design*, volume 112, pages 544-550, 1990.

a general configuration
-----------------------

For a general configuration, we generate 5 points,
with coordinates uniformly distributed in the interval :math:`[-1, +1]`.
With ``sympy``, the points are stored in object of the type ``Matrix``
for the computations in the formulation of the equations.

::

   from sympy.matrices import Matrix
   from random import uniform as u
   d0 = Matrix(2, 1, lambda i,j: u(-1,+1))
   d1 = Matrix(2, 1, lambda i,j: u(-1,+1))
   d2 = Matrix(2, 1, lambda i,j: u(-1,+1))
   d3 = Matrix(2, 1, lambda i,j: u(-1,+1))
   d4 = Matrix(2, 1, lambda i,j: u(-1,+1))

The rotation matrices involve cosines and sines of angles.

::

   c1, s1 = var('c1, s1')
   c2, s2 = var('c2, s2')
   c3, s3 = var('c3, s3')
   c4, s4 = var('c4, s4')
   R1 = Matrix([[c1, -s1], [s1, c1]])
   R2 = Matrix([[c2, -s2], [s2, c2]])
   R3 = Matrix([[c3, -s3], [s3, c3]])
   R4 = Matrix([[c4, -s4], [s4, c4]])

Then the first four equations reflect the identity
:math:`cos^2(t) + sin^(t) = 1` for any angle :math:`t`.

::
 
   p1 = 'c1^2 + s1^2 - 1;'
   p2 = 'c2^2 + s2^2 - 1;'
   p3 = 'c3^2 + s3^2 - 1;'
   p4 = 'c4^2 + s4^2 - 1;'

Two more additional unknowns are the end points of the first bar,
which is connected to the origin.

::

   x1, x2 = var('x1, x2')
   X = Matrix([[x1], [x2]])
   c1x = 0.5*(d1.transpose()*d1 - d0.transpose()*d0)
   c2x = 0.5*(d2.transpose()*d2 - d0.transpose()*d0)
   c3x = 0.5*(d3.transpose()*d3 - d0.transpose()*d0)
   c4x = 0.5*(d3.transpose()*d4 - d0.transpose()*d0)
   e1x = (d1.transpose()*R1 - d0.transpose())*X + c1x
   e2x = (d2.transpose()*R2 - d0.transpose())*X + c2x
   e3x = (d3.transpose()*R3 - d0.transpose())*X + c3x
   e4x = (d4.transpose()*R4 - d0.transpose())*X + c4x

For the equations on ``x1`` and ``x2`` we conver to strings:

::

   s1 = str(e1x[0]) + ';'
   s2 = str(e2x[0]) + ';'
   s3 = str(e3x[0]) + ';'
   s4 = str(e4x[0]) + ';'

The third group of equations on Y involve the pivot ``a``.

::

   a = Matrix([[1], [0]])
   y1, y2 = var('y1, y2')
   Y = Matrix([[y1], [y2]])
   c1y = c1x - a.transpose()*(d1 - d0)
   c2y = c2x - a.transpose()*(d2 - d0)
   c3y = c3x - a.transpose()*(d3 - d0)
   c4y = c4x - a.transpose()*(d4 - d0)
   e1y = ((d1.transpose() - a.transpose())*R1 \
        - (d0.transpose() - a.transpose()))*Y + c1y
   e2y = ((d2.transpose() - a.transpose())*R2 \
        - (d0.transpose() - a.transpose()))*Y + c2y
   e3y = ((d3.transpose() - a.transpose())*R3 \
        - (d0.transpose() - a.transpose()))*Y + c3y
   e4y = ((d4.transpose() - a.transpose())*R4 \
        - (d0.transpose() - a.transpose()))*Y + c4y

The string representations are defined as follows:

::

   s5 = str(e1y[0]) + ';'
   s6 = str(e2y[0]) + ';'
   s7 = str(e3y[0]) + ';'
   s8 = str(e4y[0]) + ';'

Then we have the polynomial system in the list:

::

   equ = [p1, p2, p3, p4, s1, s2, s3, s4, s5, s6, s7, s8]
   print('the polynomial system :')
   for pol in equ:
       print(pol)

Then, at last, we run the blackbox solver:

::

   from phcpy.solver import solve
   sols = solve(equ)
   print('the solutions :')
   for sol in sols:
       print(sol)
   print('computed', len(sols), 'solutions')

For any general choice of precision points,
the number of solutions should always be the same, that is: 36.
