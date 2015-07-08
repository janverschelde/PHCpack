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

functions in the module
-----------------------

.. automodule:: sets
   :members:
