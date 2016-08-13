tangent lines to a circle
=========================

Given a fixed circle in the plane, compute the lines tangent to
the circle and passing through the origin.
In :numref:`figtouchcircle` we see a general line through the origin
and two lines touching the circle.

.. _figtouchcircle:

.. figure:: ./touchcircle.png
    :align: center

    A general line through a circle and two lines touching a circle.

The tangent lines are special: at the points where the lines touch
the circle, we have a double solution, a solution of multiplicity two.
One method is to consider the one parameter family of lines through
the origin and intersect this family with the polynomials which express
the singularity condition.

lines through the origin intersecting a circle
----------------------------------------------

The polynomials which express all lines through the origin
intersecting a fixed circle, fixed by its center and radius,
are returned by a function.

::

   def polynomials(a, b, r):
       """
       Returns string representations of two polynomials:
       1) a circle with radius r centered at (a, b);
       2) a line through the origin with slope s.
       """
       crc = '(x - %.15e)^2 + (y - %.15e)^2 - %.15e;' % (a, b, r**2)
       lin = 'y - s*x;'
       return [crc, lin]

There are two equations, one for the circle and one for the line.
The variables are two coordinates ``x``, ``y``, and the slope ``s``.
When given two equations in three variables we expect a one
dimensional solution set.  To represent this space curve,
we intersect the curve with a general hyperplane and compute
the points on the curve and on the hyperplane.

The code snippet below defines the problem for a circle
centered at ``(3, 2)`` with radius one.  The ``embed`` function
returns the original polynomials with one general hyperplane added
and also one slack variable.  
The blackbox solver computes the generic points.

::

   pols = polynomials(3, 2, 1)
   from phcpy.sets import embed
   from phcpy.solver import solve
   embpols = embed(3, 1, pols)
   embsols = solve(embpols)

The tuple ``(embpols, embsols)`` is a numerical representation for
the set of all lines through the origin intersecting a fixed circle.

As a sanity check, consider a point on the set of all lines as
in the left of :numref:`figtouchcircle`.
Such a point is for instance the line with slope one.
The coordinates for the intersection points,
as can be seen from :numref:`figtouchcircle` are ``(2, 2)`` and ``(3, 3)``.  
In the code below, the intersection point ``(2, 2)`` is joined with
the slope ``1`` in a solution string, called ``point``.

::

   from phcpy.solutions import make_solution
   point = make_solution(['x', 'y', 's'], [2, 2, 1])
   ismb = is_member(embpols, embsols, 1, point)

The call to ``is_member`` returns a boolean value,
which should be ``True`` for this point.
