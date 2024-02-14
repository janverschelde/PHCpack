Arc Length Parameter Continuation
=================================

With increment and fix continuation, the continuation parameter :math:`t`
is fixed during the corrector stage.  
With *arc length parameter continuation*, the parameter :math:`t` is 
variable during the correction stage.  This leads to a numerically 
effective algorithm to compute quadratic turning points.

In this section, the following functions are used:

::

    from phcpy.solutions import make_solution
    from phcpy.sweepers import double_real_sweep
    from phcpy.sweepers import double_complex_sweep

computing a real quadratic turning point
----------------------------------------

Arc length parameter continuation is applied in a real sweep
which ends at a quadratic turning point.  
The homotopy is defined below.

::

    circle = ['x^2 + y^2 - 1;', 'y*(1-s) + (y-2)*s;']

At ``s`` equal to zero, ``y`` is zero and we have :math:`\pm 1` as 
the solutions for ``x``.  The two start solutions are defined by

::

    first = make_solution(['x', 'y', 's'], [1.0, 0.0, 0.0])
    second = make_solution(['x', 'y', 's'], [-1.0, 0.0, 0.0])
    startsols = [first, second]

With the start solutions defined, 
we then start the sweep in the real parameter space.

::

    newsols = double_real_sweep(circle, startsols)

and print the solutions with

::

    for (idx, sol) in enumerate(newsols):
        print('Solution', idx+1, ':')
        print(sol)

The solutions at the end of the two paths are

::

    Solution 1 :
    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x : -2.46519032881566E-32   0.00000000000000E+00
     y :  1.00000000000000E+00   0.00000000000000E+00
     s :  5.00000000000000E-01   0.00000000000000E+00
    == err :  0.000E+00 = rco :  1.000E+00 = res :  0.000E+00 =
    Solution 2 :
    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  2.46519032881566E-32   0.00000000000000E+00
     y :  1.00000000000000E+00   0.00000000000000E+00
     s :  5.00000000000000E-01   0.00000000000000E+00
    == err :  0.000E+00 = rco :  1.000E+00 = res :  0.000E+00 =

The homotopy stopped at ``s`` equal to ``0.5`` for at that value 
of ``s`` the solution for ``(x, y)`` is the double point ``(0, 1)``.

At a turning point, the real paths turn back in real space.  
If moving forward, the solution points can continue, 
but then as a pair of complex conjugated paths,
with nonzero imaginary parts.

complex parameter homotopy continuation
---------------------------------------

Sweeping the parameter space with a convex linear combination
of the parameters, no singularities are encountered.

::

    circle = ['x^2 + y^2 - 1;']

The ``circle`` defines a natural parameter homotopy,
where ``y`` is the parameter, starting at zero.
For ``y`` equal to zero, the corresponding values for ``x`` in 
the start solutions are :math:`\pm 1`, defined below:

::

    first = make_solution(['x', 'y'], [1.0, 0.0])
    second = make_solution(['x', 'y'], [-1.0, 0.0])
    startsols = [first, second]

In a complex sweep on a natural parameter homotopy, 
we have to define which variable will be the continuation parameter
and we have to provide start and target values,
giving real and imaginary parts.

::

    par = ['y']
    start = [0, 0]
    target = [2, 0]

Then we run a complex sweep in double precision as:

::

    newsols = double_complex_sweep(circle, startsols, 2, par, start, target)

The output of

::

    for (idx, sol) in enumerate(newsols):
        print('Solution', idx+1, ':')
        print(sol)

is

::

    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  0.00000000000000E+00  -1.73205080756888E+00
     y :  2.00000000000000E+00   0.00000000000000E+00
    == err :  1.282E-16 = rco :  1.000E+00 = res :  4.441E-16 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  0.00000000000000E+00   1.73205080756888E+00
     y :  2.00000000000000E+00   0.00000000000000E+00
    == err :  1.282E-16 = rco :  1.000E+00 = res :  4.441E-16 =

At the target value for ``y``, we arrived 
at a complex conjugated pair of solutions.
