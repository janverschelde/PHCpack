Reproducible Runs
=================

The correctness of a polynomial homotopy relies on the choice 
of random complex constants.  Except for an algebraic set of 
bad choices of complex constants, the solution paths are free 
of singularities, except at the end of the paths, 
when the system that is solved has singular solutions.

For correctness, it is important that the random constants 
are generated *after* the user has provided the input.

solving the cyclic 5-roots system twice
---------------------------------------

The random choice of constants makes that the solutions are computed 
in a random order, as the experiment in this section shows.

::

    from phcpy.families import cyclic

The cyclic 5-roots problem belongs to an academic family
of benchmark problems.  The code

::

    c5 = cyclic(5)
    for pol in c5:
        print(pol)

prints the polynomials:

::

    x0 + x1 + x2 + x3 + x4;
    x0*x1 + x1*x2 + x2*x3 + x3*x4 + x4*x0;
    x0*x1*x2 + x1*x2*x3 + x2*x3*x4 + x3*x4*x0 + x4*x0*x1;
    x0*x1*x2*x3 + x1*x2*x3*x4 + x2*x3*x4*x0 + x3*x4*x0*x1 + x4*x0*x1*x2;
    x0*x1*x2*x3*x4 - 1;

Now let us solve the system twice.
In each run we print the first solution.

::

    from phcpy.solver import solve

The outcome of the first run

::
  
    s1 = solve(c5)
    print(s1[0])

is

::

    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x0 : -8.09016994374948E-01   5.87785252292473E-01
     x1 :  2.11803398874989E+00  -1.53884176858763E+00
     x2 :  3.09016994374947E-01  -2.24513988289793E-01
     x3 : -8.09016994374947E-01   5.87785252292473E-01
     x4 : -8.09016994374947E-01   5.87785252292473E-01
    == err :  1.152E-15 = rco :  1.034E-01 = res :  1.665E-15 =

and the second run

::

    s2 = solve(c5)
    print(s2[0])

shows

::

    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x0 : -1.18033988749895E-01   3.63271264002680E-01
     x1 :  3.09016994374947E-01  -9.51056516295154E-01
     x2 :  3.09016994374947E-01  -9.51056516295154E-01
     x3 :  3.09016994374947E-01  -9.51056516295154E-01
     x4 : -8.09016994374947E-01   2.48989828488278E+00
    == err :  6.682E-16 = rco :  4.527E-02 = res :  1.874E-15 =

The cyclic 5-roots problem has 70 different solutions and 
with high likelihood, the first solution will be different in each run.

fixing the seed
---------------

Fixing the seed of the random number generators 
makes the solver deterministic.
Consider the following experiment.

::

    from phcpy.dimension import set_seed

We set the seed to the value ``2024`` and print the first
solution again after solving the system:

::

    set_seed(2024)
    s3 = solve(c5)
    print(s3[0])

which shows 

::

    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x0 :  1.00000000000000E+00  -8.03364927637306E-17
     x1 :  3.09016994374947E-01   9.51056516295154E-01
     x2 : -8.09016994374947E-01   5.87785252292473E-01
     x3 : -8.09016994374947E-01  -5.87785252292473E-01
     x4 :  3.09016994374947E-01  -9.51056516295154E-01
    == err :  6.315E-16 = rco :  2.393E-01 = res :  4.631E-16 =

Let us do it again:

::

    set_seed(2024)
    s4 = solve(c5)
    print(s4[0])

which then shows

::

    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x0 :  1.00000000000000E+00  -8.03364927637306E-17
     x1 :  3.09016994374947E-01   9.51056516295154E-01
     x2 : -8.09016994374947E-01   5.87785252292473E-01
     x3 : -8.09016994374947E-01  -5.87785252292473E-01
     x4 :  3.09016994374947E-01  -9.51056516295154E-01
    == err :  6.315E-16 = rco :  2.393E-01 = res :  4.631E-16 =

And of course, the point is that we see twice the same first solution.
