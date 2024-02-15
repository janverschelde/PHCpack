Equation and Variable Scaling
=============================

A polynomial system may have coefficients which vary widely in magnitude.
This may lead to inaccurate solutions.   Equation scaling multiplies 
all coefficients in the same equation by the same constant.  
Variable scaling replaces the original variables by new variables times
some constant.  Both equation and variable scaling can reduce 
the variability of the magnitudes of the coefficients of the system.

In this section, the following functions are applied:

::

    from phcpy.solutions import verify
    from phcpy.solver import solve
    from phcpy.scaling import double_scale_system
    from phcpy.scaling import double_scale_solutions

solving without scaling
-----------------------

Consider the polynomials below.

::

    p = ['0.000001*x^2 + 0.000004*y^2 - 4;', '0.000002*y^2 - 0.001*x;']

Observe the variation in magnitude between the coefficients
of the polynomials in ``p``.  First, the system is solved,
without scaling the coefficients.

::

    psols = solve(p)
    for (idx, sol) in enumerate(psols):
        print('Solution', idx+1, ':')
        print(sol)

Then the solutions are

::

    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  1.23606797749980E+03  -4.34288187660045E-10
     y : -7.86151377757428E+02   2.63023405572965E-10
    == err :  2.598E-03 = rco :  4.053E-04 = res :  8.002E-10 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  1.23606797749981E+03   1.58919085844426E-11
     y :  7.86151377757436E+02   9.62482268770893E-12
    == err :  5.733E-04 = rco :  4.053E-04 = res :  6.661E-11 =
    Solution 3 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x : -3.23606797750004E+03  -1.73527763271358E-11
     y :  4.58877485485452E-12  -1.27201964951414E+03
    == err :  1.310E-03 = rco :  2.761E-04 = res :  1.855E-10 =
    Solution 4 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x : -3.23606797749979E+03   5.48229606785362E-14
     y :  1.44974035789780E-14   1.27201964951407E+03
    == err :  6.525E-05 = rco :  2.761E-04 = res :  4.063E-14 =

Observe the magnitude of the values for the coordinates of the solutions
and the estimates for the inverse condition numbers in ``rco``.
The forward errors ``err`` are rather large and the residual ``res``
not so close to the machine precision.

The function ``verify`` evaluates the polynomials in the system 
at the solutions and returns the sum of the absolute values.

::

    verify(p, psols)

yields the number

::

    4.3339840153217635e-12

With scaling this error can be reduced.

solving after scaling
---------------------

Equation and variable scaling is applied in double precision.

::

    (q, c) = double_scale_system(p)

Consider now the coefficients of the scaled system:

::

    for pol in q:
        print(pol)

and here are the scaled polynomials:

::

    9.99999999999994E-01*x^2 + 9.99999999999996E-01*y^2 - 1.00000000000000E+00;
      " + 9.99999999999996E-01*y^2 - 9.99999999999997E-01*x;


and scaling coefficients in ``c`` are

::

    [3.30102999566398,
     0.0,
     2.999999999999999,
     0.0,
     -0.6020599913279623,
     0.0,
     -0.30102999566398136,
     0.0,
     0.04028876017114153,
     0.0]

We solve the scaled system:

::

    qsols = solve(q)
    for (idx, sol) in enumerate(qsols):
        print('Solution', idx+1, ':')
        print(sol)

which gives the solutions

::

    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x : -1.61803398874990E+00   2.94545607917864E-90
     y :  1.47272803958932E-90   1.27201964951407E+00
    == err :  1.475E-16 = rco :  2.268E-01 = res :  6.661E-16 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x : -1.61803398874990E+00   2.94545607917864E-90
     y : -1.47272803958932E-90  -1.27201964951407E+00
    == err :  1.475E-16 = rco :  2.268E-01 = res :  6.661E-16 =
    Solution 3 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  6.18033988749897E-01  -2.92604772168262E-98
     y :  7.86151377757425E-01   0.00000000000000E+00
    == err :  8.868E-17 = rco :  4.601E-01 = res :  1.110E-16 =
    Solution 4 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  6.18033988749897E-01  -2.92604772168262E-98
     y : -7.86151377757425E-01   0.00000000000000E+00
    == err :  8.868E-17 = rco :  4.601E-01 = res :  1.110E-16 =

All solutions of the scaled system are well conditioned 
with small forward and backward errors.

The scaling coefficients in ``c`` are used to bring the solutions 
of the scaled problem to the original coordinates.

::

    ssols = double_scale_solutions(len(q), qsols, c)

The output of

::

    for (idx, sol) in enumerate(ssols):
        print('Solution', idx+1, ':')
        print(sol)

is

::

    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x : -3.23606797749979E+03   5.89091215835726E-87
     y :  1.47272803958932E-87   1.27201964951407E+03
    == err :  1.475E-16 = rco :  2.268E-01 = res :  6.661E-16 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x : -3.23606797749979E+03   5.89091215835726E-87
     y : -1.47272803958932E-87  -1.27201964951407E+03
    == err :  1.475E-16 = rco :  2.268E-01 = res :  6.661E-16 =
    Solution 3 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  1.23606797749979E+03  -5.85209544336522E-95
     y :  7.86151377757424E+02   0.00000000000000E+00
    == err :  8.868E-17 = rco :  4.601E-01 = res :  1.110E-16 =
    Solution 4 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  1.23606797749979E+03  -5.85209544336522E-95
     y : -7.86151377757424E+02   0.00000000000000E+00
    == err :  8.868E-17 = rco :  4.601E-01 = res :  1.110E-16 =

Let us look at the sum of the backward errors.
Executing

::

    verify(p, ssols)

yields the error 

::

    4.4853010194856324e-14

Observe that the error is about one hundred times 
smaller than without scaling."
