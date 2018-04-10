positive dimensional solution sets
==================================

The modules **sets**, **cascades**, **factor**, and **diagonal** provide 
some functionality to work with positive dimensional solution sets.
In particular, the **solve** function of the **factor** module computes
a numerical irreducible decomposition of the solution set of a polynomial
system.  Also polynomials that have variables raised to negative powers,
so-called *Laurent polynomials* are supported.

For isolated solutions, the main outcome of the numerical solver is 
a list of points, given as tuples of values for the coordinates.
For positive dimensional solutions, with numerical homotopy continuation
methods we can compute a *numerical irreducible decomposition*
of the solution set.  Such a decomposition has two layers:

1. For every dimension of the solution set,
   we have as many :index:`generic points` as the degree 
   of the solution set of that dimension.

2. For every dimension of the solution set,
   those generic points that belong to the same irreducible factor
   are stored in the same list.

A *generic point* on a *d*-dimensional solution set is computed
as a solution of the given system of polynomial equations, augmented
with *d* linear equations with randomly generated complex coefficients.

Generic points occur as solutions in the data structure that is
called a witness set.  With embeddings and cascades, we define
homotopies in a top down calculation of a numerical irreducible
decomposition.  Diagonal homotopies define a bottom up construction
of a :index:`numerical irreducible decomposition`.
The application of monodromy loops leads to a factorization of a pure
dimensional solution set into irreducible components.
Even without having explicit equations for the irreducible factors,
with a homotopy membership test we can determine whether any given
point belongs to any given factor in the decomposition.

witness sets
------------

A *witness set* is a data structure to represent a positive dimensional
solution set, which is stored as a tuple of two items:

1. An *embedding* of the polynomial equations that define the solution set,
   augmented with as many generic linear equations as the dimension of 
   the solution set.
   To every linear equation corresponds one *slack variable*.

2. Witness points are solutions in the intersection of the original
   polynomial equations and the generic linear equations.
   For generic coefficients of the added linear equations,
   we obtain *generic points* on the solution set.
   The number of witness points equals the degree of the solution set.

In the example below we consider the twisted cubic:

::

   >>> twisted = ['x^2 - y;', 'x^3 - z;']
   >>> from phcpy.sets import embed
   >>> e = embed(3,1,twisted)
   >>> e[0]
   'x^2 - y + (-8.23538851649530E-01-5.67259869745581E-01*i)*zz1;'
   >>> e[1]
   'x^3 - z + (9.35464826338593E-01-3.53419805165623E-01*i)*zz1;'

The last equation of the embedded system is a linear equation
with randomly generated complex coefficient.
The ``zz1`` denotes the slack variable.
Continuing the session:

::

   >>> terms = e[-1].split(')*')
   >>> for t in terms: print(t)
   ... 
    + (-8.85038627286137E-01 + 4.65517591731472E-01*i
   x + (-2.12324313395875E-02 + 9.99774566519578E-01*i
   y + (-9.52478263619098E-01-3.04606561539880E-01*i
   z + (-9.59619713608467E-01 + 2.81300560351385E-01*i
   zz1+(-3.24025444378001E-01 + 9.46048366308847E-01*i);
   >>> from phcpy.solver import solve
   >>> s = solve(e, silent=True)
   >>> len(s)
   3
   >>> for sol in s: print(sol)
   ... 
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -4.06510360753325E-01   5.24436731997959E-01
    y : -1.09783212468901E-01  -4.26377930233570E-01
    zz1 : -4.27642353614751E-50  -2.73691106313441E-48
    z :  2.68236261633139E-01   1.15752697061077E-01
   == err :  3.693E-16 = rco :  1.041E-01 = res :  1.804E-16 =
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x :  7.97989261058868E-01   1.37578286563110E+00
    y : -1.25599163259885E+00   2.19571990464483E+00
    zz1 :  1.44884468793274E-32  -5.03715049801216E-33
    z : -4.02310165732919E+00   2.41891366942435E-02
   == err :  1.240E-15 = rco :  1.463E-02 = res :  2.120E-15 =
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -1.07164436617733E-01  -9.41488516596475E-01
    y : -8.74916410407434E-01   2.01788172926252E-01
    zz1 :  0.00000000000000E+00   0.00000000000000E+00
    z :  2.83741171803972E-01   8.02099237512644E-01
   == err :  9.857E-17 = rco :  8.220E-02 = res :  1.110E-16 =

The variable ``zz1`` is an artificial slack variable.
Adding the slack variable via an embedding is a general technique
to make overdetermined polynomial systems *square*,
that is: having as many equations as unknowns.
Only solutions with zero slack variables matter.

There are four homotopies which involve witness sets.

1. Given a witness set and a point, a *homotopy membership test*
   decides whether the point lies on the solution set represented
   by the witness set.

2. Given all solutions with nonzero values for the slack variables
   of an embedded system, a *cascade homotopy* takes those solutions
   as the start points of solution paths leading to generic points
   on lower dimensional solution sets.

3. Given a witness set, a *monodromy homotopy* separates the
   generic points in the witness set according to the irreducible
   factors of the solution set.

4. Given two witness sets, a *diagonal homotopy* computes witness set
   representations for all components of the intersection of the two
   given witness sets.

homotopy membership test
------------------------

Given a witness set and a point, with a homotopy we can decide
whether the point belongs to the algebraic set represented by
the given witness set.  We illustrate this membership test on
the cyclic 4-roots problem.  First we compute a witness set.

::

   >>> from phcpy.families import cyclic
   >>> c4 = cyclic(4)
   >>> from phcpy.sets import embed
   >>> c4e1 = embed(4, 1, c4)
   >>> from phcpy.solver import solve
   >>> sols = solve(c4e1)
   >>> from phcpy.solutions import filter_zero_coordinates as filter
   >>> genpts = filter(sols, 'zz1', 1.0e-8, 'select')
   >>> for sol in genpts:
   ...     print(sol)

Because there are four solutions that satisfy the original cyclic 4-roots
problem and a hyperplane with randomly generated coefficients,
there is a one dimensional solution set of cyclic 4-roots.

The function ``membertest`` takes as input the witness set,
represented by the polynomials in ``c4e1`` 
and the generic points in ``genpts``, and a point.
The point is given as a list of doubles, with the real and imaginary
parts of all coordinates.  The point ``(1, -1, 1, -1)`` is thus
given as the list ``[1, 0, -1, 0, 1, 0, -1, 0]``.  The four extra
zeroes are the zero imaginary parts of the four coordinates.

::

   >>> point = [1, 0, -1, 0, 1, 0, -1, 0]
   >>> from phcpy.sets import membertest
   >>> membertest(c4e1, genpts, 1, point)
   residual is  4.00000000000000E+00
     point does not lie on the component, as residual >  1.000E-06
   False

The function ``membertest`` returns ``False`` as the residual of
the evaluation of the point at the equations does not satisfy the
default tolerance. 

Testing the point ``(-1, -1, 0, 0)`` proceeds as follows.
The ``...`` below stands for omitted output.

::

   >>> point = [-1, 0, -1, 0, 1, 0, 1, 0]
   >>> membertest(c4e1, genpts, 1, point)
   residual is  0.00000000000000E+00
     point satisfies the equations, as residual <=  1.000E-06
   ===========================================================================
   == 1 =  #step :  48 #fail :  0 #iter :  61 = regular solution ==
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1                  Length of path :  3.20529115248889E-02
   the solution for t : 
    x0 : -9.99999989106364E-01   7.42569305742753E-09
    x1 : -1.00000001089364E+00  -7.42569325135928E-09
    x2 :  1.00000000001044E+00   7.19523989696616E-12
    x3 :  9.99999999989557E-01  -7.19544395567596E-12
    zz1 : -3.89859685757465E-16  -1.56365136224798E-16
   == err :  2.235E-08 = rco :  2.239E-09 = res :  1.079E-15 ==
   ...
       match with generic point 1, as difference is  1.074E-08 <=  1.000E-06
     Point lies on a solution component.
   True

The point passes the residual test.  The test continues
with the computation of new generic points for a hyperplane
that passes through the test point.  If the test point is among
the new generic points, then the test point belongs to the
positive dimensional solution set represented by the witness set.
For this example we see that the point ``(-1, -1, 1, 1)`` is
a singular point on the curve, as can be seen from the estimate
for the inverse condition number, ``rco :  2.239E-09``.
The default tolerance of ``1.0e-6`` is high enough in this case
for the point to satisfy the membership test.

If the tolerance ``1.0e-6`` is deemed too sloppy,
then we can allow for a stronger tolerance and execute
the homotopy membership test in double double precision.
More zeroes must be inserted in the test point for the second part
(the least significant double) in the double double representation
for the real and imaginary parts of the coordinates:

::

   >>> ddpoint = [-1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]

Instead of ``1.0e-6``, the new tolerance is ``1.0e-12``:

::

   >>> membertest(c4e1, genpts, 1, ddpoint, memtol=1.e-12, precision='dd')
   residual is 0.00000000000000000000000000000000E+00
     point satisfies the equations, as residual <=  1.000E-06
   ...
   == 3 =  #step :  43 #fail :  9 #iter : 119 = singular solution
   t : 1.00000000000000000000000000000000E+00  0.00000000000000000000000000000000E+00
   m : 1                  Length of path :  2.58366452243257E+01
   the solution for t : 
    x0 : -1.00000000000000000000460097514793E+00  -5.03546515557825836758428944198701E-21
    x1 : -9.99999999999999999995399044465921E-01  5.03544121195596798412885036927768E-21
    x2 : 9.99999999999999999999995602319048E-01  -4.84887236589678161814305345115216E-24
    x3 : 1.00000000000000000000000441729480E+00  4.87281598848246759732685155434612E-24
    zz1 : 3.21588798303478472454372240227796E-34  -4.61270682182747444272903697352050E-34
   == err :   1.525E-13 = rco :   1.343E-14 = res :   8.687E-26 ==
   ...
       match with generic point 3, as difference is 4.183E-17 <=  1.000E-12
     Point lies on a solution component.
   True

In double double precision, the condition number estimate for the
inverse condition number drops to ``1.343E-14`` (see the ``rco`` field).

To perform the membership test in quad double precision,
invoke ``membertest`` with ``precision='qd'``.

For solution sets of large degree, the homotopy membership test will
run faster in its multitasked version.  To run the membership test
with 8 tasks, add ``tasks=8`` as last argument of the call to the function. 

cascade of homotopies
---------------------

With a cascade of homotopies, we separate generic points on one
equidimensional component from another equidimensional component
of the solution set.  A cascade starts at the top dimension.
We consider an illustrative example:

::

   >>> pols = ['(x^2 + y^2 + z^2 - 1)*(y - x^2)*(x - 0.5);', \
               '(x^2 + y^2 + z^2 - 1)*(z - x^3)*(y - 0.5);', \
               '(x^2 + y^2 + z^2 - 1)*(z - x*y)*(z - 0.5);']

The polynomials in ``pols`` are defined in factored form
so for this illustrative example we may read of the equidimensional
components of the solution set, which contain the two dimensional
sphere, the one dimensional twisted cubic, and the isolated point
``(0.5, 0.5, 0.5)``.

To initialize the cascade, we must have solved an embedded polynomial system.
With ``embed(3, 2, pols)`` we make an embedding of the 3-dimensional
system in ``pols`` adding two linear equations with random complex
coefficients.  Two slack variables ``zz1`` and ``zz2`` are added to
make this overdetermined system square.

::

   >>> from phcpy.sets import embed
   >>> topemb = embed(3, 2, pols)
   >>> from phcpy.solver import solve
   >>> topsols = solve(topemb, silent=True)

The list ``topsols`` contains two types of solutions:
those with nonzero values for the slack variables, and
those with zero slack variables, which thus satisfy the original
equations in ``pols`` and the two added linear equations with random
complex coefficients.  The solutions with zero values for the slack
variables define generic points on the two dimensional solution set.
We filter the solutions, as follows:

::

   >>> from phcpy.solutions import filter_zero_coordinates as filter
   >>> topsols0 = filter(topsols, 'zz2', 1.0e-8, 'select')
   >>> topsols1 = filter(topsols, 'zz2', 1.0e-8, 'remove')
   >>> print('generic points on the two dimensional surface :')
   >>> for sol in topsols0:
   ...     print(sol)

The solutions with nonzero values for the slack variables are
called *nonsolutions*.  These solutions are regular and serve
as start solutions in a cascade to compute generic points on 
the lower dimensional components of the solution set.

::

   >>> from phcpy.cascades import cascade_step
   >>> lvl1sols = cascade_step(2, topemb, topsols1)

After the filtering, we must drop variables, coordinates,
and hyperplane for the next level in the cascade.

::

   >>> from phcpy.sets import drop_variable_from_polynomials as drop1poly
   >>> from phcpy.sets import drop_coordinate_from_solutions as drop1sols
   >>> lvl1emb = drop1poly(topemb, 'zz2')
   >>> lvl1emb = lvl1emb[:-1]  # dropping the last polynomial
   >>> lvl1solsdrop = drop1sols(lvl1sols, len(topemb), 'zz2')
   >>> lvl1sols0 = filter(lvl1solsdrop, 'zz1', 1.0e-8, 'select') 
   >>> lvl1sols1 = filter(lvl1solsdrop, 'zz1', 1.0e-8, 'remove') 

Among the solutions at the end of the paths defined by the cascade
homotopy are solutions that belong to the two dimensional sphere.
These solutions are singular and we filter then away based on
threshold for the estimate of the inverse condition number.

::

   >>> from phcpy.solutions import filter_regular as regfilt
   >>> reglvl1sols0 = regfilt(lvl1sols0, 1.0e-8, 'select')
   >>> for sol in reglvl1sols0:
   ...     print(sol)

To find the isolated solutions, another cascade homotopy is applied,
tracking the paths starting at the nonsolutions at the end of the
previous cascade.

::

   >>> lvl2sols = cascade_step(1, lvl1emb, lvl1sols1)
   >>> lvl2solsdrop = drop1sols(lvl2sols, len(lvl1emb), 'zz1')
   >>> for sol in reglvl2solsdrop:
   ...     print(sol)

To perform the filtering of the solutions properly,
we apply a membership test, defined in the **sets** module.

The function ``run_cascade()`` takes on input the number of variables
in the polynomials and the top dimenion of the solution set.
Starting at the top dimension, a witness set representation for
each pure dimensional component of the solution set is computed.

factoring into irreducibles
---------------------------

A witness set consists of two parts.
The first part of a witness set is a polynomial system with as many added
linear equations with random coefficients as the dimension.
The number of slack variables (variables that start with the name ``zz``)
equals the dimension of the witness set.
The second part of a witness set is a list of solutions of the first part.
Because the added linear equations have random coefficients,
the solutions are generic points on the positive dimensional algebraic set.

Given a witness set, applying monodromy loops those points in a witness set
that lie on the same irreducible factor are joined.
The application of monodromy is a probabilistic method with unknown
probability of failure because it relies on the unknown distribution
of the singular solutions.

Below is a simple example, given already in factored form:

::

    >>> p = '(x+1)*(x^2 + y^2 + 1);'

To construct a witness set we import
``witness_set_of_hypersurface`` from ``phcpy.sets``:

::

   >>> from phcpy.sets import witness_set_of_hypersurface as wh
   >>> (w, s) = wh(2, p)
   >>> len(s)

Because the degree of ``p`` is three,
we see ``3`` as the outcome of ``len(s)``.

::

   >>> from phcpy.factor import factor
   >>> f = factor(1, w, s)
   >>> f

The result in ``f`` is a a list of tuples:

::

   [([1, 2], 8.537360146292391e-15), ([3], 2.1316282072803006e-14)]

The factorization joined the first two solutions of `s` 
as they represent the quadratic factor.
A generic point for the linear factor is in the second tuple.
The second floating point number in each tuple is the residual
obtained via the linear trace test, used as stop criterion in
the running of monodromy loops.

For polynomials of higher degrees, double double or even quad double
could be required to obtain accurate results.
The following two commands illustrate how to apply monodromy
respectively in double double and quad double precision:

::

    >>> f = factor(1, w, s, precision='dd')
    >>> f = factor(1, w, s, precision='qd')

The witness set ``(w, s)`` should also have been computed in
double double and quad double precision.

The function ``decompose()`` takes the output of the ``run_cascade()``
function of the *cascades* module and factors every witness set for
the pure dimensional components into irreducible factors.
The functions ``run_cascade()`` and ``decompose()`` lead to a
numerical irreducible decomposition of the solution set.

numerical irreducible decomposition
-----------------------------------

Consider the polynomials defined by the list ``pols`` as follows:

::

    >>> pols = ['(x1-1)*(x1-2)*(x1-3)*(x1-4);', \
                '(x1-1)*(x2-1)*(x2-2)*(x2-3);', \
                '(x1-1)*(x1-2)*(x3-1)*(x3-2);', \
                '(x1-1)*(x2-1)*(x3-1)*(x4-1);']
    >>> from phcpy.factor import solve, write_decomposition
    >>> deco = solve(4, 3, pols, verbose=False)
    >>> write_decomposition(deco)

We see the common factor ``x1-1`` which defines a three dimensional
solution plane.  The factor ``x1-2`` leads to a two dimensional
solution plane, with the additional factor ``x2-1``.
Furthermore, the system in ``pols`` has twelve lines as solutions
and four isolated solution points.

The first argument ``4`` in the ``solve(4, 3, pols, verbose=False)`` is the
number of variables in the polynomials in ``pols``.
The second argument ``3`` equals the top dimension of the solution set.
The ``write_decomposition()`` confirms there is one three dimensional
linear component, one two dimensional linear component, twelve lines,
and four isolated solutions.

diagonal homotopies
-------------------

Given two witness sets, with diagonal homotopies we can compute 
generic points on the intersection of the algebraic sets represented
by the witness sets, and thus obtain a witness set of the intersection.
This section illustrates the intersection of the unit sphere with
a cylinder.  This intersection defines a quartic curve.

We start with equations for the unit sphere and a cylinder:

::

   >>> sph = 'x^2 + y^2 + z^2 - 1;'
   >>> cyl = 'x^2 + y - y + (z - 0.5)^2 - 1;'

Observe the ``+ y - y`` line in the assignment to ``cyl``.
With this trick we initialize the symbol table for the witness set
computation, ensuring that ``y`` is present.

Next, we compute a witness sets for the sphere and the cylinder:

::

   >>> from phcpy.sets import witness_set_of_hypersurface as witsurf
   >>> sphwit = witsurf(3, sph)
   >>> spheqs, sphpts = sphwit
   >>> cylwit = witsurf(3, cyl)
   >>> cyleqs, cylpts = cylwit

Once we have two witness sets, we call the ``diagonal_solver``
method to compute a witness set for the intersection:

::

   >>> from phcpy.diagonal import diagonal_solver as diagsolve
   >>> quawit = diagsolve(3, 2, spheqs, sphpts, 2, cyleqs, cylpts)
   >>> quaeqs, quapts = quawit
   >>> for pol in quaeqs:
   ...     print(pol)
   >>> for sol in quapts:
   ...     print(sol)
