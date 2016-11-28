some interesting examples and families
======================================

PHCpack has been tested on many examples of polynomial systems
taken from the research literature.
The module examples exports some of those examples.
Running **python examples.py** at the command prompt
performs a regression test, solving all examples.

Polynomial systems often occur in families and are defined
for any number of equations and variables.

interactive regression testing
------------------------------

An interactive use of examples.py at the Python prompt can go as follows:

::

   >>> from phcpy.examples import noon3
   >>> f = noon3()
   >>> for p in f: print p
   ... 
   x1*x2^2 + x1*x3^2 - 1.1*x1 + 1;
   x2*x1^2 + x2*x3^2 - 1.1*x2 + 1;
   x3*x1^2 + x3*x2^2 - 1.1*x3 + 1;
   >>> 

The functions in examples.py returns the polynomials as lists of strings.
If we want to solve the system defined by f, we continue the above session as

::

   >>> from phcpy.solver import solve
   >>> s = solve(f,silent=True)
   >>> len(s)
   21
   >>> print s[0]
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x1 : -1.65123467890611E-01  -7.61734168646636E-01
    x2 :  8.98653694263692E-01  -3.48820047576431E-01
    x3 :  8.98653694263692E-01  -3.48820047576431E-01
   == err :  3.034E-16 = rco :  2.761E-01 = res :  5.974E-16 =
   >>> 

The example session continues in the description of the module solutions.

the cyclic n-roots problem
--------------------------

One such noteworthy family is the cyclic n-roots problem:

::

   >>> from phcpy.families import cyclic
   >>> c4 = cyclic(4)
   >>> for p in c4: print p
   ... 
   x0 + x1 + x2 + x3;
   x0*x1 + x1*x2 + x2*x3 + x3*x0;
   x0*x1*x2 + x1*x2*x3 + x2*x3*x0 + x3*x0*x1;
   x0*x1*x2*x3 - 1;
   >>> 
