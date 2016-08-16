the design of a 4-bar mechanism
===============================

Given two pivot points and five precision points for the coupler,
the design problem asks to determine the lengths of the bars that
allow the coupler to pass through the given precision points.

The equations are generated with version 1.0 of ``sympy``
and the plots are made with version 1.5.1 of ``matplotlib``.

The system is taken from a paper by A.P. Morgan and C.W. Wampler
on **Solving a Planar Four-Bar Design Using Continuation**, published in
the *Journal of Mechanical Design*, volume 112, pages 544-550, 1990.
In :numref:`fig4barline`, the precision points are
taken from Problem 7 in the paper.

.. _fig4barline:

.. figure:: ./fbarline.png
    :align: center

    A straight line design of a 4-bar mechanism.

The first plot in :numref:`fig4barline`,
at the top leftmost corner shows the five precision points,
labeled with the numbers 0, 1, 2, 3, and 4.
The two white triangles in each plot represent the fixed pivots.
The next five plots show one position of the 4-bar mechanism.
Each position passes through one of the prescribed precision points.
The rotation angles and the coordinates for ``x`` and ``y`` for the
initial position are obtained as solutions of a polynomial system.

For the formulation of the equations we follow the notation of the problem 
statement in the second section of the paper by Morgan and Wampler.
The first pivot point is fixed at the origin and the coordinates of the
other fixed pivot point are in :math:`a = (a_1, a_2)`.
The coordinates of the five precision points are denoted by
:math:`d_0, d_1, d_2, d_3,` and :math:`d_4`.
All vectors are column vectors and the superscript :math:`T` denotes
the transpose.  The planar rotation matrices are defined by

.. math::

   R_j = \left[
           \begin{array}{rr}
              c_j & -s_j \\
              s_j & c_j
           \end{array}
         \right], \quad j = 1, 2, 3, 4,

where :math:`c_j` and :math:`s_j` are respectively the cosines
and sines of the rotation angles.  The first four equations express
the relationship between cosines and sines in the identities

.. math::

   c_j^2 + s_j^2 - 1 = 0, \quad j = 1, 2, 3, 4.

The second group of equations involves the vector :math:`x = (x_1, x_2)`
of variables.  
The first bar in the mechanism is between the pivots.
The line segment between the first pivot at :math:`(0, 0)`
and :math:`x` represents the second bar in the 4-bar mechanism.

.. math::

   \left[ d^T_j R_j - d^T_0 \right] x
   + \frac{1}{2} \left[ d^T_j d_j - d^T_0 d_0 \right] = 0,
   \quad j = 1, 2, 3, 4.

The third bar in the mechanism is represented between :math:`x`
and :math:`y = (y_1, y_2)` and the fourth bar connects :math:`y`
and the second pivot at :math:`a`.
The third group of equations involving :math:`y` is defined by

.. math::

  \left[ \left( d^T_j - a^T \right) R_j
       - \left( d^T_0 - a^T \right) \right] y
  + \left[ \frac{1}{2} \left( d^T_j d_j - d^T_0 d_0 \right)
           - a^T \left( d_j - d_0 \right) \right], \quad j = 1, 2, 3, 4.

So we end up with a system of twelve equations in twelve unknowns:
:math:`c_1, s_1, c_2, s_2, c_3, s_3, c_4, s_4, x_1, x_2, y_1, y_2`
and ten parameters, the coordinates of the precision points
:math:`d_0, d_1, d_2, d_3`, and :math:`d_4`.
The coordinates of the second fixed pivot point :math:`a`
are typically set to be :math:`(1, 0)`.

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
   print 'the polynomial system :'
   for pol in equ:
       print pol

Then, at last, we run the blackbox solver:

::

   from phcpy.solver import solve
   sols = solve(equ)
   print 'the solutions :'
   for sol in sols:
       print sol
   print 'computed', len(sols), 'solutions'

For any general choice of precision points,
the number of solutions should always be the same, that is: 36.

visualization of a straight line design
---------------------------------------

Of special interest are those 4-bar mechanisms where the five
precision points are on a line, as such mechanisms can be applied
to translate circular into linear motion or otherwise.

The coordinates of the following five precision points are
copied from Problem 7 of the paper by Morgan and Wampler:

::

    pt0 = Matrix([[ 0.50], [ 1.06]])
    pt1 = Matrix([[-0.83], [-0.27]])
    pt2 = Matrix([[-0.34], [ 0.22]])
    pt3 = Matrix([[-0.13], [ 0.43]])
    pt4 = Matrix([[ 0.22], [ 0.78]])

These are the coordinates shown in :numref:`fig4barline`.
There are 33 solutions to the polynomial system formulated in the
same fashion as in the previous section.  
Of those 33 solutions, 15 are real.
Only real solutions can lead to valid designs.
Not every real solution leads to a valid design.
One condition is that the four angles computed from the cosine
and sine coordinates must be ordered, so that the precision points
are reached the same order as they are listed in the input.

In :numref:`fig4barcoupler`,
the *coupler curve* for the straight line mechanism is shown.
This coupler curve is traced by the tip of the triangle moved
by the 4-bar mechanism.
The ``matplotlib`` code is available in the script ``fourbar.py``
in the ``examples`` folder of the source code for the Python 2
and Python 3 distributions.

.. _fig4barcoupler:

.. figure:: ./fbarcoupler.png
    :align: center

    The coupler curve of a straight line 4-bar mechanism.
