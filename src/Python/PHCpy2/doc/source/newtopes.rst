Newton polytopes, monomial maps, and power series
=================================================

The Newton polytopes of the polynomial system provide important
information about the structure of the solution sets.
The module **phcpy.polytopes** provides an interface to the convex hull
methods of PHCpack.  It also provides a directer interface to the
mixed volume calculator, directer in the sense that the user can enter
the supports directly, without having to formulate a polynomial system.

Systems that have exactly two monomials with nonzero coefficient
in every equation are called binomial systems.
Although such binomial systems are very particular,
because of their sparse structure, they can be solved much faster.
The module **phcpy.maps** provides a Python interface to the
solvers of binomial systems.

The classical arithmetic can be extended to the field of truncated
power series.  In this field, Newton's method computes power series
solutions of polynomial systems.
The module **phcpy.series** exports functions to compute power series
with Newton's method.  
Power series are input to Padé approximants,
which lead to more accurate path trackers,
exported by the module **phcpy.curves**.

convex hulls of lattice polytopes
---------------------------------

The session below illustration the calculation of the convex hull
of a configuration of seven points in the plane.
The points are generated at random, with coordinates between -9 and +9.

::

   >>> from phcpy.polytopes import random_points as rp
   >>> from phcpy.polytopes import planar_convex_hull as pch
   >>> points = rp(2,7,-9,9)
   >>> points
   [(9, 8), (5, 6), (6, 0), (2, -5), (4, -1), (9, -4), (-1, -6)]
   >>> (vertices, normals) = pch(points)
   >>> vertices
   [(9, 8), (5, 6), (-1, -6), (9, -4)]
   >>> normals
   [(1, -2), (2, -1), (-1, 5), (-1, 0)]

The output of the convex hull method consists of a tuple of two lists.
The first list is the list of vertices.  For this particular example,
seven points were given on input, and only four of those points are 
corners of the convex hull.  The list of vertices is ordered cyclically:
two consecutive vertices span an edge of the polygon and the last and
first vertex also span an edge as a polygon has exactly as many vertices
as edges.  The second list in the output is the list of inner normals,
which are vectors perpendicular to the edges.  
Taking the inner product of the normal with the points that span an edge
yields the same value for each point on the edge, and that value is minimal
for all points in the polygon.  For the example above 
for the inner normal (1, -2) and the two points (9, 8) and (5, 6),
we have

.. math::

    9 - 2 \times 8 = 5 - 2 \times 6 = -7 

as the edge lies on the edge of the half plane defined by
the inequality

.. math::

    x_1 - 2 x_2 \geq -7

which holds for all points in the polygon spanned by the points
in the example.  The inner normals define the half planes that cut
out the polygon.

For a convex hull of a point configuration in 3-space, 
consider the example in the session below:

::

   >>> from phcpy.polytopes import random_points as rp
   >>> points = rp(3,10,-9,9)
   >>> for point in points: print point
   ... 
   (5, 9, -5)
   (0, 0, 1)
   (-3, -4, -1)
   (-9, -3, -3)
   (-5, 3, -8)
   (-4, 3, 7)
   (2, -3, 8)
   (9, 3, -9)
   (7, 4, -2)
   (1, -8, 1)
   >>> from phcpy.polytopes import convex_hull as ch
   >>> facets = ch(3, points)
   computed 12 facets
   >>> for facet in facets: print facet
   ... 
   (-597, [90, -65, -6], [4, 5, 6], [1, 2, 3])
   (-84, [1, 11, 14], [4, 8, 5], [7, 11, 0])
   (-281, [30, -49, -2], [5, 1, 6], [11, 4, 0])
   (-51, [6, 5, -6], [4, 6, 7], [0, 4, 5])
   (-203, [-22, -27, -30], [6, 1, 7], [2, 9, 3])
   (-48, [5, 6, -5], [4, 7, 10], [3, 6, 7])
   (-684, [-127, 66, -29], [10, 7, 8], [5, 8, 7])
   (-150, [1, 22, 25], [4, 10, 8], [5, 6, 1])
   (-315, [-59, 15, -19], [7, 9, 8], [9, 10, 6])
   (-265, [-29, -35, -39], [7, 1, 9], [4, 10, 8])
   (-165, [-19, -10, -4], [1, 8, 9], [11, 8, 9])
   (-429, [3, -26, 42], [5, 8, 1], [1, 10, 2])

The output of the ``convex_hull`` function returns a list of facets.
Each facet is represented as a tuple of four items.
The first number is the value of the inner product of the vector
perpendicular to the facet, given by the list in the second item
of the tuple.  So the first two items in the tuple define the
half space defined by the facet.  For the first facet, we have
the inequality defined by the number -597 and the vector [90, -65, -6]:

.. math::

   90 x_1 - 65 x_2 - 6 x_3 \geq -597

which holds for all points \ :math:`(x_1, x_2, x_3)` in the convex hull.  
The equality \ :math:`90 x_1 - 65 x_2 - 6 x_3 = -597` holds
for all points that lie on the first facet in the list of facets above.
The third item in the representation of a facet is the list of numbers
to the points that span the facet.  In the example above, the first
facet is spanned by the points 4, 5, 6 in the input list points.
Note that the counting of the points starts at one and not at zero.
The last item in the representation of a facet is the list of 
facets that are adjacent to the facet.  For the first facet,
facets 1, 2, and 3 are adjacent to it.  The counting of the facets
starts at zero, so the first facet has label zero.

From the list of facets we can extract all vertex points.
If we continue with the session from above:

::

   >>> vertices = []
   >>> for facet in facets:
   ...     for point in facet[2]:
   ...         if not point in vertices:
   ...             vertices.append(point)
   ... 
   >>> vertices
   [4, 5, 6, 8, 1, 7, 10, 9]
   >>> len(vertices)
   8

We have 8 vertices and 12 facets.  The points the span the facets are
ordered cyclically so that two consecutive points span an edge and the
last and first point span also an edge.  Every edge lies in the intersection
of exactly two facets.  Edges of adjacent facets are ordered in opposite
order.  For example, facet 0 is spanned by [4, 5, 6] and its adjacent
facet 1 is spanned by [4, 8, 5], with the edge shared between both of
them oriented from 4 to 5 in facet 0 and from 5 to 4 in facet 1.

As the points in the configuration were generated sufficiently at
random, the polytope is simplicial: every facet is spanned by exactly
3 points and has exactly 3 edges.  As every edge is shared by exactly
two facets we count every edge twice if we multiply the number of facets
by three, so we have 36/2 = 18 edges.

mixed volumes
-------------

The mixed volume of a tuple of Newton polytopes
if defined as the coefficient in the expansion of the volume
of a linear combination of Newton polytopes.
For example, for a 3-tuple of Newton polytopes:

.. math::

    \begin{array}{rcl}
      vol(\lambda_1 P_1 + \lambda_2 P_2 + \lambda_3 P_3)  
      & = & V(P_1, P_1, P_1) \lambda_1^3 \\
      & + & V(P_1, P_1, P_2) \lambda_1^2 \lambda_2 \\
      & + & V(P_1, P_2, P_2) \lambda_1 \lambda_2^2 \\
      & + & V(P_1, P_2, P_3) \lambda_1 \lambda_2 \lambda_3 \\
      & + & V(P_2, P_2, P_2) \lambda_2^3 \\
      & + & V(P_2, P_2, P_3) \lambda_2^2 \lambda_3 \\
      & + & V(P_2, P_3, P_3) \lambda_2 \lambda_3^2 \\
      & + & V(P_3, P_3, P_3) \lambda_3^3
    \end{array}

where \ :math:`vol(\cdot)` is the volume function
and \ :math:`V(\cdot)` is the mixed volume.
For the tuple \ :math:`(P_1, P_2, P_3)`, its mixed volume
is \ :math:`V(P_1,P_2,P_3)` in the expansion above.

The function ``mixed_volume`` expects two arguments.
The first argument is the list of exponents of
the \ :math:`\lambda` variables in the volume expansion formula.
The second argument of ``mixed_volume`` is a tuple of Newton polytopes.
The session below illustrates the computation of the volume of one
single polytope.

::

   >>> from phcpy.polytopes import random_points as rp
   >>> from phcpy.polytopes import mixed_volume as mv
   >>> p1 = rp(3, 5, -9, 9)
   >>> print p1
   [(3, 7, -3), (-1, 0, 8), (-6, -6, 8), (-6, 9, 4), (-3, 4, -7)]
   >>> mv([3], [p1])
   2107

The volume is normalized, so the standard unit simplex has volume one.
To compute mixed volumes of two polytopes, we continue the session,
generating another polytope:

::

   >>> p2 = rp(3, 5, -9, 9)
   >>> mv([2, 1],(p1, p2))
   3910
   >>> mv([1, 2],(p1, p2))
   3961

The ``mixed_volume`` function executes and Ada translation of
MixedVol, ACM TOMS Algorithm 846 of 2005.
This algorithm generates random floating point values to lift the
points in the supports.  
The function ``integer_mixed_cells`` allows the user to specify integer 
lifting values as the last coordinate of the points in the supports.

solving binomial systems
------------------------

The irreducible components of
positive dimensional solution sets of binomial systems
have coordinates that can be represented by maps of monomials 
in free independent variables.  In this representation, there
are as many free variables as the dimension of the solution set.
The module ``maps`` exports a solver for binomial systems.

In the example below, we consider a simple system
of two binomials in three variables:

::

   >>> f = [ 'x**2*y - z*x;', 'x**2*z - y**2*x;' ]
   >>> from phcpy.maps import solve_binomials
   >>> maps = solve_binomials(3,f)
   >>> for map in maps: print map

In the printed maps, we recognize the twisted cubic,
the x-axis, and the yz-plane as the three solution sets.

power series solutions
----------------------

Newton's method applies also to systems where the coefficients are
truncated power series.
The module ``series`` exports functions to compute power series solutions
in double, double double, and quad double precision.
The function ``test()`` of the ``series`` module provides an example.

As example, we consider the Viviani curve and intersect the curve
with a moving plane.  The parameter ``s`` defines the movement of
the plane ``y = 0`` to the plane ``y = 1``, as in the setup below:

::

   >>> vivplane = ['(1-s)*y + s*(y-1);',
   ... 'x^2 + y^2 + z^2 - 4;',
   ... '(x-1)^2 + y^2 - 1;']
   >>> vivs0 = vivplane + ['s;']
   >>> from phcpy.solver import solve
   >>> sols = solve(vivs0, silent=True)
   >>> print sols[0]
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    s :  0.00000000000000E+00   0.00000000000000E+00
    y :  0.00000000000000E+00   0.00000000000000E+00
    x :  0.00000000000000E+00   0.00000000000000E+00
    z :  2.00000000000000E+00   0.00000000000000E+00
   == err :  0.000E+00 = rco :  3.186E-01 = res :  0.000E+00 =
   >>> 

It is important that the parameter ``s`` is the first symbol
in the polynomials in the input (in the list ``vivplane`` above)
for Newton's method to compute series solutions.
In the session below, the output is formatted with continuation symbols.

::

   >>> from phcpy.series import standard_newton_series
   >>> sersols = standard_newton_series(vivplane, sols, verbose=False)
   >>> sersols[0]
   ['s;', '3.12500000000000E-02*s^8 + 6.25000000000000E-02*s^6 \
   + 1.25000000000000E-01*s^4 + 5.00000000000000E-01*s^2;', \
   ' - 2.07519531250000E-02*s^8 - 4.10156250000000E-02*s^6 \
   - 7.81250000000000E-02*s^4 - 2.50000000000000E-01*s^2 + 2;']
   >>>

Starting at the solution for ``s = 0``, the series solution
allows to predict the solution as the plane moves away from ``y = 0``
towards ``y = 1``.

approximating algebraic curves
------------------------------

Power series are the input to algorithms to construct rational
approximations, also called Padé approximants.

A session can start with the tuning of the homotopy
continuation parameters, as indicated below:

::

   >>> from phcpy.curves import tune_homotopy_continuation_parameters as tune
   >>> tune()
   Values of the HOMOTOPY CONTINUATION PARAMETERS :
    1. gamma : (-0.797398052335-0.603453681844j)
    2. degree of numerator of Pade approximant    : 4
    3. degree of denominator of Pade approximant  : 4
    4. maximum step size                          : 0.1
    5. minimum step size                          : 1e-06
    6. multiplication factor of the series step   : 0.5
    7. multiplication factor of the pole radius   : 0.5
    8. tolerance on the residual of the predictor : 0.001
    9. tolerance on the residual of the corrector : 1e-8
   10. tolerance on zero series coefficients      : 1e-12
   11. maximum number of corrector steps          : 4
   12. maximum steps on a path                    : 1000
   To change a value, give an index (0 to exit) : 0

The function ``tune()`` enters an interactive loop.
In this loop the user can alter all values for the twelve parameters.
In the session above, the user entered zero to exit the loop
without modifications to the default values of the parameters.
The parameters can be queried and altered separately with the
``get_`` and ``set_`` functions of the module curves.

To illustrate the application of the trackers,
we consider a small example of the Katsura family of systems.

::

   >>> from phcpy.families import katsura
   >>> k3 = katsura(3)
   >>> from phcpy.solver import total_degree_start_system as startsys
   >>> (k3q, k3qsols) = startsys(k3)
   >>> len(k3qsols)
   8

The total degree start system for this problem of three quadrics
has eight solutions.  With the target system in ``k3``, 
the start system in ``k3q``, and the start solutions in ``k3qsols``,
we can launch the path tracker in standard double precision as follows:

::

   >>> from phcpy.curves import standard_track as track
   >>> k3sols = track(k3, k3q, k3qsols, "/tmp/out", True)

The solutions at the end of the paths are assigned to ``k3sols``.
The output will be written to the file with name ``/tmp/out``
and the verbose flag is set to ``True``.
The verbose flag is optional.  If no output to file is needed,
then also the file name may be omitted.
If an emtpy string is given as the value for the file name
and ``True`` for the verbose option, then all extra output
is written to screen.

To print the last 40 lines of ``/tmp/out`` one can do the following:

::

   >>> file = open("/tmp/out")
   >>> lines = file.readlines()
   >>> for k in range(-40,0): print lines[k][:-1]

The module exports functions to run the path trackers step by step.
In each step, the user can retrieve the power series,
the Padé approximants, its poles, the step size and solution.

To test this step-by-step path tracking, do

::

   >>> from phcpy.curves import test_next_track as test
   >>> test()

The default precision is double, providing 'dd' or 'qd' to the argument
of ``test()`` sets the precision respectively to double double
and quad double.
