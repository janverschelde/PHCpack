Newton's Method, Deflation and Multiplicity
===========================================

The instructions in the code cells in this chapter
require the following functions to be imported:

::

    from phcpy.solutions import make_solution, strsol2dict
    from phcpy.deflation import double_newton_step
    from phcpy.deflation import double_double_newton_step
    from phcpy.deflation import quad_double_newton_step
    from phcpy.deflation import double_deflate
    from phcpy.deflation import double_multiplicity

Newton's method
---------------

Let us start with a simple system of two polynomials in two unknowns.

::

    pols = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']

At a regular solution, Newton's method doubles the accuracy in each step.
For the example, we start at a solution where the coordinates are
given with only three decimal places.

::

    sols = [make_solution(['x', 'y'], [1.23, -0.786])]
    print(sols[0])

Here is then the solution with three decimal places of accuracy:

::

    t : 0.000000000000000E+00 0.000000000000000E+00
    m : 1
    the solution for t :
     x : 1.230000000000000E+00  0.0
     y : -7.860000000000000E-01  0.0
    == err : 0.000E+00 = rco : 1.000E+00 = res : 0.000E+00 =

The function ``double_newton_step`` expects a list of solution strings
as its second argument.

::

    sols = double_newton_step(pols, sols)
    print(sols[0])

The outcome of one step with Newton's method in double precision is

::

    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 0
    the solution for t :
     x :  1.23607623318386E+00   0.00000000000000E+00
     y : -7.86154018188250E-01   0.00000000000000E+00
    == err :  6.230E-03 = rco :  1.998E-01 = res :  3.706E-05 =

Observe the values for ``err`` (forward error),
``rco`` (estimated inverse of the condition number), and
``res`` (the residual or backward error).

The multiplicity field ``m`` turned ``0`` because the default 
tolerance was not reached and the solution could not be counted 
as a proper solution.  Let us reset the multiplicity, 
making a new solution with the values from the ``sols[0]``.
The string representation of the solution is converted
into a dictionary first:

::

    sold = strsol2dict(sols[0])
    sold

which yields

::

    {'t': 0j,
     'm': 0,
     'err': 0.00623,
     'rco': 0.1998,
     'res': 3.706e-05,
     'x': (1.23607623318386+0j),
     'y': (-0.78615401818825+0j)}

and then we make a new solution:

::

    sol = make_solution(['x', 'y'], [sold['x'], sold['y']])
    print(sol)

As input for the next Newton step we will use

::

    t : 0.000000000000000E+00 0.000000000000000E+00
    m : 1
    the solution for t :
     x : 1.236076233183860E+00  0.000000000000000E+00
     y : -7.861540181882500E-01  0.000000000000000E+00
    == err : 0.000E+00 = rco : 1.000E+00 = res : 0.000E+00 =

Calling ``double_newton_step`` for the second time

::

    sols = [sol]
    sols = double_newton_step(pols, sols)
    print(sols[0])

gives

::

    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 0
    the solution for t :
     x :  1.23606797751503E+00   0.00000000000000E+00
     y : -7.86151377766704E-01   0.00000000000000E+00
    == err :  1.090E-05 = rco :  1.998E-01 = res :  1.100E-10 =

Observe that the value of ``res`` dropped from magnitude ``1.0e-5``
down to ``1.0e-10``, corresponding to the well conditioning of the root.
However, the multiplicy field is still zero because the estimate 
for the forward error is still too high.

The third application of the Newton step

::

    sold = strsol2dict(sols[0])
    sol = make_solution(['x', 'y'], [sold['x'], sold['y']])
    sols = [sol]
    sols = double_newton_step(pols, sols)
    print(sols[0])

then yields the solution

::

    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  1.23606797749979E+00   0.00000000000000E+00
     y : -7.86151377757423E-01   0.00000000000000E+00
    == err :  2.452E-11 = rco :  1.998E-01 = res :  4.441E-16 =

The value for the residual ``res`` is very close to machine precision
and the solution is considered a proper regular solution.

We can double the precision to double double,
and continue with the previous approximate solution:

::

    sols = double_double_newton_step(pols, sols)
    print(sols[0])

what is printed is below

::

    t : 0.00000000000000000000000000000000E+00      0.00000000000000000000000000000000E+00    
    m : 1
    the solution for t :
     x : 1.23606797749978969640917366873130E+00      0.00000000000000000000000000000000E+00    
     y : -7.86151377757423286069558585843026E-01     0.00000000000000000000000000000000E+00    
    == err :  3.036E-16 = rco :  1.998E-01 = res :  3.944E-31 =

and doing this once more

::

    sols = double_double_newton_step(pols, sols)\n",
    print(sols[0])"

shows

::

    t : 0.00000000000000000000000000000000E+00      0.00000000000000000000000000000000E+00    
    m : 1
    the solution for t :
     x : 1.23606797749978969640917366873128E+00      0.00000000000000000000000000000000E+00    
     y : -7.86151377757423286069558585842966E-01     0.00000000000000000000000000000000E+00    
    == err :  6.272E-32 = rco :  1.998E-01 = res :  0.000E+00 =

and then on to quad double precision:

::

    sols = quad_double_newton_step(pols, sols)
    print(sols[0])

which then gives the solution:

::

    t : 0.0000000000000000000000000000000000000000000000000000000000000000E+00      0.0000000000000000000000000000000000000000000000000000000000000000E+00    
    m : 1
    the solution for t :
     x : 1.2360679774997896964091736687312762354406183596115257242708972454E+00      0.0000000000000000000000000000000000000000000000000000000000000000E+00    
     y : -7.8615137775742328606955858584295892952312205783772323766490197015E-01     0.0000000000000000000000000000000000000000000000000000000000000000E+00    
    == err :  7.070E-33 = rco :  1.998E-01 = res :  2.101E-64 =

deflation
---------

At an isolated singular solution, *deflation* is a method to restore
the quadratic convergence of Newton's method.

::

    pols = ['(29/16)*x^3 - 2*x*y;', 'x^2 - y;']

Obviously, ``(0, 0)`` is a root of the ``pols`` but even if we start
very close to ``(0, 0)``, Newton's method converges extremely slowly.

::

    sols = [make_solution(['x', 'y'],[1.0e-6, 1.0e-6])]
    print(sols[0])

The initial solution is then

::

    t : 0.000000000000000E+00 0.000000000000000E+00
    m : 1
    the solution for t :
     x : 1.000000000000000E-06  0.0
     y : 1.000000000000000E-06  0.0
    == err : 0.000E+00 = rco : 1.000E+00 = res : 0.000E+00 =

Running one step with Newton's method:

::

    sols = double_newton_step(pols, sols)
    print(sols[0])

gives

::

    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 0
    the solution for t :
     x :  9.99999906191101E-07   0.00000000000000E+00
     y :  9.99999812409806E-13   0.00000000000000E+00
    == err :  1.000E-06 = rco :  5.625E-13 = res :  1.875E-19 =

Observe the tiny value for ``rco``, the estimate of the inverse
of the condition number.

The next step with Newton's method

::

    sols = double_newton_step(pols, sols)
    print(sols[0])

shows

::

     t :  0.00000000000000E+00   0.00000000000000E+00
     m : 0
     the solution for t :\n",
      x :  6.66666604160106E-07   0.00000000000000E+00
      y :  3.33333270859482E-13   0.00000000000000E+00
     == err :  3.333E-07 = rco :  2.778E-14 = res :  1.111E-13 =

which does not brings us much closer to ``(0, 0)``.

Now we apply deflation in double precison:

::

    solsd = double_deflate(pols, sols)
    print(solsd[0])

which then yields

::

    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  9.46532112069346E-24   4.09228221015004E-24
     y :  1.02357542351685E-24  -2.03442589046821E-24
    == err :  5.292E-12 = rco :  5.608E-03 = res :  1.885E-15 =

Deflation also works on systems with more equations than unknowns.

::

    pols = ['x^2;', 'x*y;', 'y^2;']

Once again, ``(0, 0)`` is the obvious root of ``pols``.

::

    sols = [make_solution(['x', 'y'], [1.0e-6, 1.0e-6])]
    print(sols[0]

So, we start at

::

    t : 0.000000000000000E+00 0.000000000000000E+00
    m : 1
    the solution for t :
     x : 1.000000000000000E-06  0.0
     y : 1.000000000000000E-06  0.0
    == err : 0.000E+00 = rco : 1.000E+00 = res : 0.000E+00 =

and apply deflation

::

    sols = double_deflate(pols, sols, tolrnk=1.0e-8)
    print(sols[0])

which shows

::

    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  1.25000000000000E-07   0.00000000000000E+00
     y :  1.25000000000000E-07   0.00000000000000E+00
    == err :  1.250E-07 = rco :  8.165E-01 = res :  1.562E-14 =

This is not an improvement, because the tolerance on the numerical
rank, given by ``tolrnk=1.0e-8`` is too severe.
If we relax this tolerance to ``1.0e-4``, 

::

    sols = double_deflate(pols, sols, tolrnk=1.0e-4)
    print(sols[0])

then we obtain

::

    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x : -5.87747175411144E-39  -5.14278778484751E-39
     y :  2.93873587705572E-39  -1.83670992315982E-40
    == err :  2.757E-23 = rco :  4.082E-01 = res :  1.984E-38 =

and the coordinates of this numerical approximation are
well below the double precision.

multiplicity structure
----------------------

The multiplicity can be computed *locally* starting at the solution.
Consider the following example:

::

    pols = [ 'x**2+y-3;', 'x+0.125*y**2-1.5;']

The solution ``(1, 2)`` is a multiple root of ``pols``.

Let us prepare the input:

::

    sol = make_solution(['x', 'y'], [1, 2])
    print(sol)

so, we work with

::

    t : 0.000000000000000E+00 0.000000000000000E+00
    m : 1
    the solution for t :
     x : 1.000000000000000E+00  0.0
     y : 2.000000000000000E+00  0.0
    == err : 0.000E+00 = rco : 1.000E+00 = res : 0.000E+00 =

as input to ``double_multiplicity``

::

    multiplicity, hilbert_function = double_multiplicity(pols, sol)

The value of ``multiplicity`` equals ``3``
which is confirmed by the content of ``hilbert_function``:

::

    [1, 1, 1, 0, 0, 0]

Thus, the multiplicity of ``(1, 2)`` as a root of ``pols`` equals three.
