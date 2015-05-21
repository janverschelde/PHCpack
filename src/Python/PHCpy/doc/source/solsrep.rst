representations of solutions of polynomial systems 
==================================================

Solutions of phcpy.solve are returned as lists of PHCpack
solution strings.  The solutions module contains functions to
parse a PHCpack solution string into a dictionary.

The solutions module exports operations 

1. to parse strings in the PHCpack solution format into dictionaries;

2. to evaluate these dictionaries into polynomials substituting the
   values for the variables into the strings representing the polynomials.

The main test in the module solutions is the solution of a small
trinomial system and the evaluation of the computed solutions
at the trinomial system.

The information of a solution as a dictionary contains the following:

1. `t` : value of the continuation parameter

   `m` : multiplicity of the solution

2. symbols for the variables are keys in the dictionary,
   the corresponding values are complex floating-point numbers

3. `err` : magnitude of the last correction term of Newton's method
   (forward error)

   `rco` : estimate for the inverse of the condition number of
   the Jacobian matrix at the solution

   `res` : magnitude of the residual (backward error)

The triplet (err, rco, res) measures the numerical quality of the solution.
The residual `res` is normally interpreted as an estimate of the backward
error: by how much should we change the original problem such that the
approximate solution becomes an exact solution of the changed problem.
The estimate `rco` gives a (sometimes too pessimistic) bound on the
number of correct decimal places in the approximate solution.
In particular: `abs(log(rco, 10))` bounds the number of lost decimal
places in the approximate solution.
For example, if `rco` equals `1.0E-8`, then the last 8 decimal places
in the coordinates of the solution could be wrong.

For a solution of the example ``noon3`` from the module examples,
we convert the PHCpack format solution string to a dictionary as follows:

::

   >>> print s[0]
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x1 : -1.65123467890611E-01  -7.61734168646636E-01
    x2 :  8.98653694263692E-01  -3.48820047576431E-01
    x3 :  8.98653694263692E-01  -3.48820047576431E-01
   == err :  3.034E-16 = rco :  2.761E-01 = res :  5.974E-16 =
   >>> from phcpy.solutions import strsol2dict
   >>> d = strsol2dict(s[0])
   >>> d.keys()
   ['err', 'res', 'm', 'rco', 't', 'x2', 'x3', 'x1']
   >>> d['x1']
   (-0.165123467890611-0.761734168646636j)
   >>> 

Note that the values of the dictionary d are evaluated strings,
parsed into Python objects.

By plain substitution of the values of the dictionary representation
of the solution into the string representation of the polynomial system
we can verify that the coordinates of the solution evaluate to numbers
close to the numerical working precision:

::

   >>> from phcpy.solutions import evaluate
   >>> e = evaluate(f,d)
   >>> for x in e: print x
   ... 
   (1.11022302463e-15+4.4408920985e-16j)
   (7.77156117238e-16+9.99200722163e-16j)
   (7.77156117238e-16+9.99200722163e-16j)
   >>> 

A more elaborate verification of the solution is provided by
the function **newton_step** of the module ``solver`` of phcpy.

functions in the module
-----------------------

The documentation strings of the functions
exported by the module ``solutions`` are listed below.
The script **test()** runs when typing **python solutions.py**
at the command prompt.

.. automodule:: solutions
   :members:
