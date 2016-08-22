a blackbox solver for isolated solutions
========================================

The package phcpy depends on the shared object file phcpy2c.so.
The module **phcpy.solver**
exports the blackbox solver of PHCpack, a fast mixed volume
calculator, and several methods to predict the number of isolated
solutions of a polynomial system.  
The `test_solver()` function of the module generates two trinomials 
(a polynomial with three monomials)
with randomly generated complex coefficients.

By default, the input polynomial systems are expected to be *square*,
that is: the number of polynomials in the input list equals the number
of variables in the polynomials.  The blackbox solver then returns
a list of numerical approximations to the isolated solutions of the
input polynomial system.  Some capabilities of PHCpack to deal with
positive dimensional solution sets are exported by 
the module **phcpy.sets**.

solving random trinomials
-------------------------

Polynomials and solutions are represented as strings.
Below is an illustration of a session with the blackbox solver
on a system of two random trinomials, polynomials with three
monomials with random complex coefficients.

::

   >>> from phcpy.solver import random_trinomials
   >>> f = random_trinomials()

To see what the polynomials in f look like, 
let us print its terms (splitting on the `+(` 
to avoid long swirling lines in the output):

::
   
   >>> terms = f[0].split('+(')
   >>> for t in terms: print t
   ...
   -0.991457094247+0.13043324065066*i)*x^0*y^5
   -0.0509953121395-0.99869889262970*i)*x^5*y^3
   0.232308584664+0.97264213433887*i)*x^4*y^4;

To solve the system defined by f, we call the blackbox solver:

::

   >>> from phcpy.solver import solve
   >>> s = solve(f,silent=True)
   >>> len(s)
   15
   >>> print s[2]
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x :  7.10290847804173E-01  -4.69841154290980E-01
    y : -1.79580717006598E-01  -7.16541556066137E-01
   == err :  1.986E-16 = rco :  2.676E-01 = res :  1.232E-16 =

The *solve* command returned a list of 15 strings in s,
each string represents a solution that makes the polynomials in f vanish.
The module **phcpy.solutions** (documented in the next section)
offers a function to evaluate the solutions
in the polynomials given as strings.

By default, the option ``silent`` is set to ``False`` and the solver
prints the computed root counts.  Other options of the solver are

1. **tasks**: the number of tasks for multithreaded path tracking.
   Solving sufficiently large systems on 4 processor cores may
   result in a speedup of close to a factor 4 if ``tasks=4`` is
   given as input argument of ``solve.``

2. **precision**: by default the ``precision`` is set to ``d`` for
   standard hardware double precision.  While this precision may suffice,
   the blackbox solver supports two additional levels of precision:
   ``dd`` for double double preicsion (about 32 decimal places), and
   ``qd`` for quad double precision (about 64 decimal places).
   Given ``precision=dd`` as extra input parameter to ``solve``
   is likely to yield more accurate results, at an extra cost,
   which may be compensated by multithreading.

3. **checkin**: by default this flag is set to ``True`` to check
   whether the system given on input has as many polynomials as
   variables.  The current version of the blackbox solver accepts
   only square systems.  See the section on positive dimensional
   solution sets for functions that deal with overdetermined or
   underdetermined polynomial systems.

representations of solutions of polynomial systems 
--------------------------------------------------

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

Note that the values of the dictionary d are evaluated strings,
parsed into Python objects.

By plain substitution of the values of the dictionary representation
of the solution into the string representation of the polynomial system
we can verify that the coordinates of the solution evaluate to numbers
close to the numerical working precision:

::

   >>> from phcpy.solutions import evaluate
   >>> e = evaluate(f, d)
   >>> for x in e: print x
   ... 
   (1.11022302463e-15+4.4408920985e-16j)
   (7.77156117238e-16+9.99200722163e-16j)
   (7.77156117238e-16+9.99200722163e-16j)

A more elaborate verification of the solution is provided by
the function **newton_step** of the module ``solver`` of phcpy.

The module exports function to filter regular solutions, solutions
with zero coordinates or real solutions.  The filtering of real
solutions is illustrated in the session below.
We first define one real solution and another with a coordinate
that has a nonzero imaginary part.

::

   >>> from phcpy.solutions import make_solution
   >>> s0 = make_solution(['x', 'y'], [1, complex(0, 1)])
   >>> print s0
   t : 0.0 0.0
   m : 1
   the solution for t :
    x : 1.000000000000000E+00  0.0
    y : 0.000000000000000E+00  1.000000000000000E+00
   == err : 0.0 = rco : 1.0 = res : 0.0 ==
   >>> s1 = make_solution(['x', 'y'], [2, 3])
   >>> print s1
   t : 0.0 0.0
   m : 1
   the solution for t :
    x : 2.000000000000000E+00  0.0
    y : 3.000000000000000E+00  0.0
   == err : 0.0 = rco : 1.0 = res : 0.0 ==

The filtering of real solutions (with respect to a given tolerance)
is provided by the functions ``is_real`` (on one solution)
and ``filter_real`` (on a list of solutions).

::

   >>> from phcpy.solutions import is_real, filter_real
   >>> is_real(s0, 1.0e-8)
   False
   >>> is_real(s1, 1.0e-8)
   True
   >>> realsols = filter_real([s0, s1], 1.0e-8, 'select')
   >>> for sol in realsols: print sol
   ... 
   t : 0.0 0.0
   m : 1
   the solution for t :
    x : 2.000000000000000E+00  0.0
    y : 3.000000000000000E+00  0.0
   == err : 0.0 = rco : 1.0 = res : 0.0 ==

The functions ``filter_regular`` and ``filter_zero_coordinates``
operate in a manner similar as ``filter_real.``

reproducible runs with fixed seeds
----------------------------------

The solver in PHCpack generates different random numbers with each run,
which may very well cause the solutions to appear in a different order
after a second application of solve on the same system.
To prevent this behaviour (to check reproducibility for example),
we can fix the seed of the random number generators in PHCpack,
as follows:

::

   >>> from phcpy.phcpy2c import py2c_set_seed
   >>> py2c_set_seed(2013)
   0

The above session continues as

::

   >>> from phcpy.phcpy2c import py2c_get_seed
   >>> py2c_get_seed()
   2013

To reproduce a computation, we can thus request the seed that was used
(with ``py2c_get_seed``) and then restart the session setting the seed
to what was used before (with ``py2c_set_seed``).

root counting methods
---------------------

The performance of the solver is very sensitive to how accurately
we can predict the number of solutions.  For dense polynomial systems,
looking at the highest degrees of the polynomials in the system suffices,
whereas for sparse polynomial systems, computing the mixed volume of
the Newton polytopes of the polynomials yields much better results.
Below is a simple example, illustrating the bounds based on the
degrees and the mixed volume:

::

   >>> f = ['x^3*y^2 + x*y^2 + x^2;', 'x^5 + x^2*y^3 + y^2;']
   >>> from phcpy.solver import total_degree
   >>> total_degree(f)
   25
   >>> from phcpy.solver import m_homogeneous_bezout_number as mbz
   >>> mbz(f)
   (19, '{ x }{ y }')
   >>> from phcpy.solver import linear_product_root_count as lrc
   >>> lrc(f)
   a supporting set structure :
        { x }{ x }{ x }{ y }{ y }
        { x }{ x }{ x y }{ x y }{ x y }
   the root count : 19
   19
   >>> from phcpy.solver import mixed_volume
   >>> mixed_volume(f, stable=True)
   (14, 18)

The mixed volume is a generically sharp root count for the number of 
isolated solutions with all coordinates different from zero. 
The term *generically sharp* means: except for systems with coefficients 
in a specific collection of algebraic sets, the root count is an exact count.
The stable mixed volume counts all affine solutions, 
that is: also the solutions with zero coordinates.
For the example above, we may expect at most 14 isolated solutions 
with all coordinates different from zero, 
and, also considering solutions with zero coordinates, 
at most 18 isolated solutions, counted with multiplicities.

For every root count, total degree, m-homogeneous Bezout number,
linear-product root count, and mixed volume, there is a corresponding
method to construct a polynomial system with exactly as many regular
solutions at the root count, which can then be used as a start system
in a homotopy to compute all isolated solutions of the polynomial system 
for which the root count was computed.
Examples of the methods to construct start systems in phcpy
are illustrated in the documentation for the module **phcpy.trackers**.

Newton's method and deflation
-----------------------------

Newton's method fails when the Jacobian matrix is singular
(or close to singular) at a solution.  Below is a session
on the example of A. Griewank and M. R. Osborne, in their paper
*Analysis of Newton's method at irregular singularities,*
published in *SIAM J. Numer. Anal.* 20(4): 747-773, 1983.
The origin (0,0) is an irregular singularity: Newton's method
fails no matter how close the initial guess is taken.
With deflation we can restore the quadratic convergence
of Newton's method:

::

   >>> p = ['(29/16)*x^3 - 2*x*y;', 'x^2 - y;']
   >>> from phcpy.solutions import make_solution
   >>> s = make_solution(['x','y'],[1.0e-6,1.0e-6])
   >>> print s
   t : 0.0 0.0
   m : 1
   the solution for t :
    x : 1.000000000000000E-06  0.0
    y : 1.000000000000000E-06  0.0
   == err : 0.0 = rco : 1.0 = res : 0.0 ==
   >>> from phcpy.solver import newton_step
   >>> s2 = newton_step(p,[s])
   == err :  1.000E-06 = rco :  5.625E-13 = res :  1.875E-19 =
   >>> print s2[0]
   t :  0.00000000000000E+00   0.00000000000000E+00
   m : 0
   the solution for t :
    x :  9.99999906191101E-07   0.00000000000000E+00
    y :  9.99999812409806E-13   0.00000000000000E+00
   == err :  1.000E-06 = rco :  5.625E-13 = res :  1.875E-19 =
   >>> s3 = newton_step(p,s2)
   == err :  3.333E-07 = rco :  2.778E-14 = res :  1.111E-13 =
   >>> print s3[0]
   t :  0.00000000000000E+00   0.00000000000000E+00
   m : 0
   the solution for t :
    x :  6.66666604160106E-07   0.00000000000000E+00
    y :  3.33333270859482E-13   0.00000000000000E+00
   == err :  3.333E-07 = rco :  2.778E-14 = res :  1.111E-13 =
   >>> from phcpy.solver import deflate
   >>> sd = deflate(p,[s])
   >>> print sd[0]
   t :  0.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -4.55355758042535E-25   2.75154683741089E-26
    y :  1.57904709676279E-25  -8.86785799319512E-26
   == err :  5.192E-13 = rco :  5.314E-03 = res :  1.388E-16 =

The decision to deflate or not depend on the tolerance to
decide the numerical rank.  Consider the following session:

::

   from phcpy.solutions import make_solution
   from phcpy.solver import standard_deflate
   sol = make_solution(['x', 'y'], [1.0e-6, 1.0e-6])
   print sol
   pols = ['x**2;', 'x*y;', 'y**2;']
   sols = standard_deflate(pols, [sol], tolrnk=1.0e-8)
   print sols[0]
   sols = standard_deflate(pols, [sol], tolrnk=1.0e-4)
   print sols[0]

The default value for ``tolrnk`` equals ``1.0e-6``.
If we do not want to deflate that soon, we can lower the tolerance
to ``1.0e-8`` and in that case, there is no deflation when the
approximation is still as far as ``1.0e-6`` from the exact solution.
Increasing the value for the tolerance to ``1.0e-4`` leads to the
deflation at the approximation for the solution.

equation and variable scaling
-----------------------------

Another source of numerical difficulties are systems
that have extreme values as coefficients.
With equation and variable scaling we solve an optimization problem
to find coordinate transformations that lead to better values for
the coefficients.  The common sense approach to scaling is 
described in Chapter 5 of the book of Alexander Morgan on
*Solving Polynomial Systems Using Continuation for Engineering
and Scientific Problems*, volume 57 in the SIAM Classics in
Applied Mathematics, 2009.  We consider a simple example.

::

   >>> from phcpy.solver import solve
   >>> p = ['0.000001*x^2 + 0.000004*y^2 - 4;', '0.000002*y^2 - 0.001*x;']
   >>> psols = solve(p, silent=True)
   >>> print psols[0]
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -3.23606797749979E+03   8.71618409420601E-19
    y :  2.30490982555757E-19   1.27201964951407E+03
   == err :  2.853E-07 = rco :  2.761E-04 = res :  9.095E-13 =

Observe the rather large values of the coordinates in the first solution
and the estimate for the inverse condition number.
We scale the system as follows:

::

   >>> from phcpy.solver import standard_scale_system as scalesys
   >>> from phcpy.solver import standard_scale_solutions as scalesols
   >>> (q, c) = scalesys(p)
   >>> q[0]
   'x^2 + 9.99999999999998E-01*y^2 - 1.00000000000000E+00;'
   >>> q[1]
   'y^2 - 1.00000000000000E+00*x;'

The coefficients in the scaled system look indeed a lot nicer.
In the parameter ``c`` returned along with the scaled system
are the scaling coefficients, which we need to bring the solutions
of the scaled system into the original coordinates.

::

   >>> qsols = solve(q, silent=True)
   >>> ssols = scalesols(len(q), qsols, c)
   >>> for sol in ssols: print sol
   ... 
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -3.23606797749978E+03  -1.98276706040285E-115
    y :  0.00000000000000E+00  -1.27201964951407E+03
   == err :  1.746E-16 = rco :  2.268E-01 = res :  2.220E-16 =
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -3.23606797749978E+03  -1.98276706040285E-115
    y :  0.00000000000000E+00   1.27201964951407E+03
   == err :  1.746E-16 = rco :  2.268E-01 = res :  2.220E-16 =
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x :  1.23606797749979E+03   0.00000000000000E+00
    y :  7.86151377757423E+02   0.00000000000000E+00
   == err :  4.061E-17 = rco :  4.601E-01 = res :  5.551E-17 =
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x :  1.23606797749979E+03   0.00000000000000E+00
    y : -7.86151377757423E+02   7.38638289422858E-124
   == err :  4.061E-17 = rco :  4.601E-01 = res :  5.551E-17 =

The estimates of the condition numbers in ``ssols`` are for
the scaled problem.  With scaling, the condition numbers were
reduced from 10^4 to 10.  For more extreme values of the
coefficients, we may have to perform the scaling in higher precision,
such as available in the functions
``dobldobl_scale_system`` and ``quaddobl_scale_system``,
respectively with double double and quad double arithmetic.
