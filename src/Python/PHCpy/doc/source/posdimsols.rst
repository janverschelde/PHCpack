positive dimensional solution sets
==================================

The module sets.py provides some functionality of PHCpack
to work with positive dimensional solution sets.

witness sets
------------

A witness set is a data structure to represent a positive dimensional
solution set.  A witness set consists of an embedding of the polynomial
equations that define the solution set, augmented with as many generic
linear equations as the dimension of the solution set.
Witness points are solutions in the intersection of the original
polynomial equations and the generic linear equations.
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
with randomly generated complex coefficient.  Continuing the session:

::

   >>> terms = e[-1].split(')*')
   >>> for t in terms: print t
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
   >>> for sol in s: print sol
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

cascade of homotopies
---------------------

With a cascade of homotopies, we separate generic points on one
equidimensional component from another equidimensional component
of the solution set.  A cascade starts at the top dimension.
We consider an illustrative example:

::

   pols = ['(x^2 + y^2 + z^2 - 1)*(y - x^2)*(x - 0.5);', \
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

   from phcpy.sets import embed
   topemb = embed(3, 2, pols)
   from phcpy.solver import solve
   topsols = solve(topemb, silent=True)

The list ``topsols`` contains two types of solutions:
those with nonzero values for the slack variables, and
those with zero slack variables, which thus satisfy the original
equations in ``pols`` and the two added linear equations with random
complex coefficients.  The solutions with zero values for the slack
variables define generic points on the two dimensional solution set.
We filter the solutions, as follows:

::


   from phcpy.solutions import filter_zero_coordinates as filter
   topsols0 = filter(topsols, 'zz2', 1.0e-8, 'select')
   topsols1 = filter(topsols, 'zz2', 1.0e-8, 'remove')
   print 'generic points on the two dimensional surface :'
   for sol in topsols0:
       print sol

The solutions with nonzero values for the slack variables are
called *nonsolutions*.  These solutions are regular and serve
as start solutions in a cascade to compute generic points on 
the lower dimensional components of the solution set.

::

   from phcpy.sets import cascade_step
   lvl1sols = cascade_step(topemb, topsols1)

After the filtering, we must drop variables, coordinates,
and hyperplane for the next level in the cascade.

::

   from phcpy.sets import drop_variable_from_polynomials as drop1poly
   from phcpy.sets import drop_coordinate_from_solutions as drop1sols
   lvl1emb = drop1poly(topemb, 'zz2')
   lvl1emb = lvl1emb[:-1]  # dropping the last polynomial
   lvl1solsdrop = drop1sols(lvl1sols, len(topemb), 'zz2')
   lvl1sols0 = filter(lvl1solsdrop, 'zz1', 1.0e-8, 'select') 
   lvl1sols1 = filter(lvl1solsdrop, 'zz1', 1.0e-8, 'remove') 

Among the solutions at the end of the paths defined by the cascade
homotopy are solutions that belong to the two dimensional sphere.
These solutions are singular and we filter then away based on
threshold for the estimate of the inverse condition number.

::

   from phcpy.solutions import filter_regular as regfilt
   reglvl1sols0 = regfilt(lvl1sols0, 1.0e-8, 'select')
   for sol in reglvl1sols0:
       print sol

To find the isolated solutions, another cascade homotopy is applied,
tracking the paths starting at the nonsolutions at the end of the
previous cascade.

::

   lvl2sols = cascade_step(lvl1emb, lvl1sols1)
   lvl2solsdrop = drop1sols(lvl2sols, len(lvl1emb), 'zz1')
   for sol in reglvl2solsdrop:
       print sol

To perform the filtering of the solutions properly, we apply
a membership test.

diagonal homotopies
-------------------

Given two witness sets, with diagonal homotopies we can compute 
generic points on the intersection of the algebraic sets represented
by the witness sets, and thus obtain a witness set of the intersection.
This section illustrates the intersection of the unit sphere with
a cylinder.  This intersection defines a quartic curve.

We start with equations for the unit sphere and a cylinder:

::

   sph = 'x^2 + y^2 + z^2 - 1;'
   cyl = 'x^2 + y - y + (z - 0.5)^2 - 1;'

Observe the ``+ y - y`` line in the assignment to ``cyl``.
With this trick we initialize the symbol table for the witness set
computation, ensuring that ``y`` is present.

Next, we compute a witness sets for the sphere and the cylinder:

::

   from phcpy.sets import witness_set_of_hypersurface as witsurf
   sphwit = witsurf(3, sph)
   spheqs, sphpts = sphwit
   cylwit = witsurf(3, cyl)
   cyleqs, cylpts = cylwit

Once we have two witness sets, we call the ``diagonal_solver``
method to compute a witness set for the intersection:

::

   from phcpy.sets import diagonal_solver as diagsolve
   quawit = diagsolve(3, 2, spheqs, sphpts, 2, cyleqs, cylpts)
   quaeqs, quapts = quawit
   for pol in quaeqs:
       print pol
   for sol in quapts:
       print sol

functions in the module
-----------------------

.. automodule:: sets
   :members:
