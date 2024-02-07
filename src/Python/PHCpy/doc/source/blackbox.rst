The Blackbox Solver
===================

A *blackbox* solver runs with default settings and values 
of the tolerances, executing a selection of algorithms 
which has shown to work for a large collection of problems.  
As to what the *solver* computes, two different types of solvers 
are available

1. to approximate all isolated solutions,
   typically of systems with as many equations as unknows; and

2. to compute a numerical irreducible decomposition of the entire 
   solution set, which includes the isolated solution points,
   but also all positive dimensional sets (curves and surfaces), 
   factored into irreducible components.

Because the output of the two types is so vastly different and 
because the complexity of a numerical irreducible decomposition 
is much higher than computing only the isolated solutions,
the user must decide in advance which solver function to call.

approximating all isolated solutions
------------------------------------

The input to the solver is a list of strings,
with symbolic representations of the polynomials in the system.
The ``;`` at the end signifies the ``= 0`` of the equations.

::

    polynomials = ['x^3 + 2*x*y - x^2;', 'x + y - x^3;']

To call the blackbox solver, we import the ``solve`` function 
from the ``solver`` module.

::

    from phcpy.solver import solve
    solutions = solve(polynomials)

The output of ``solve`` is also a list of strings.

::

    for (idx, sol) in enumerate(solutions):
        print('Solution', idx+1, ':')
        print(sol)

has as output

::

    Solution 1 :
    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 2
    the solution for t :
     x :  0.00000000000000E+00   0.00000000000000E+00
     y :  0.00000000000000E+00   0.00000000000000E+00
    == err :  7.124E-17 = rco :  0.000E+00 = res :  0.000E+00 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  1.00000000000000E+00   0.00000000000000E+00
     y :  0.00000000000000E+00   0.00000000000000E+00
    == err :  2.100E-74 = rco :  5.348E-01 = res :  0.000E+00 =
    Solution 3 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x : -1.50000000000000E+00   0.00000000000000E+00
     y : -1.87500000000000E+00   3.30053725596684E-238
    == err :  6.747E-80 = rco :  1.220E-01 = res : 1.073E-237 =

How did this actually work?  
We can ask to see more output of the solver, 
giving a value to the verbose level parameter, 
the ``vrblvl`` as argument to the solver.

::

    solutions = solve(polynomials, vrblvl=1)

shows the following

::

    in solve, tasks : 0, mvfocus : 0, precision : d
    the polynomials on input :
    x^3 + 2*x*y - x^2;
    x + y - x^3;
    nbr : 4, roco :
    total degree : 9
    2-homogeneous Bezout number : 6
      with with partition : { x }{ y }
    general linear-product Bezout number : 5
      based on the set structure :
         { x }{ x y }{ x }
         { x y }{ x }{ x }
    mixed volume : 2
    stable mixed volume : 4

The ``roco`` is an abbreviation of *root count*.
A root count is an upper bound on the number of solutions.
In the example, the two input polynomials are cubics.
Therefore, the largest number of isolated solutions equals nine,
the product of the degrees of the polynomials.
However, not all monomials of degree three or less that could 
appear with nonzero coefficient are present.  
The numbers ``6`` and ``5`` are better bounds computed on the degrees.
Mixed volumes are generically sharp bounds for the number 
of isolated solutions with all coordinates different from zero.
As in the example, there is one double root at ``(0, 0)``,
which is counted by the stable mixed volume.

If the coefficients of the polynomials are sufficiently independent
from each other, then the number of isolated solutions counted with
multiplicity will match one of the computed root counts.

options of the solve function
-----------------------------

The blackbox solver runs with default values,
which can be changed when calling the ``solve`` function.

In addition to the list of polynomials, 
there are six parameters with default values:
 
1. ``tasks`` equals the number of the threads and 
   is by default set to zero.
   The solutions can approximated independently to each other 
   using *p* threads could in the best case speed up the solver
   by a factor of *p*.

2. ``mvfocus`` is a flag to apply only mixed volumes as root counts,
   if set to one.
   By default, this flag is set to zero and the solver will compute
   all bounds based on the degrees which may be too time consuming 
   for sparse polynomial systems.

3. ``precision`` is by default double (``d``).
   Other values are ``dd`` for double double and ``qd`` for quad double.
   While running in higher precision leads to more accurate results,
   the computational overhead can be significant.
   The overhead may be compensated (in part) by multithreading.

4. ``checkin`` is the option to run some basic checks on the input.
   By default ``True``, setting it to ``False`` is okay
   if the input polynomials are automatically generated
   in the correct way.

5. ``dictionary_output`` is by default ``False`` and a list of strings
   is returned.  If ``True``, then the output is a list of dictionaries,
   often convenient for processing.

6. ``vrblvl`` is the verbose level parameter to make the blackbox 
   solver less black when set to higher values.

Of course, then there is always the ``help(solve)``:

::

    help(solve)

which shows the documentation string of the function.

When changing the default values, consider the following.

1. The number of threads should never be set to a value higher
   than the number of available cores on the system.
   To find out the number of available cores, do the following:

   ::

       from phcpy.dimension import get_core_count
       get_core_count()

   So, use up to the value returned 
   by ``phcpy.dimension.get_core_count()`` as the value 
   to assign to the parameter ``tasks``.

2. The focus on mixed volumes (in the option ``mvfocus``) is
   automatically applied when the polynomials have negative exponents.

3. When computing in higher precision, keep in mind that also 
   the coefficients of the polynomials then must also be evaluated
   in higher precision.  Consider ``1/3`` in double precision:

   ::

       1/3

   which evaluates to

   ::

       0.3333333333333333

   which is of course not equal to the rational number ``1/3``.

4. One of the checks done by default (``checkin=True``) is 
   whether the number of polynomials in the list equals
   the number of unknowns.  At this stage, if the syntax of 
   the polynomial is incorrect, an error message will be printed as well.

5. If ``dictionary_output`` is wanted after a run,
   then it can be computed afterwards, the ``solve()``
   should not be called again,
   but can be computed with the ``strsol2dict()``
   of the ``solutions`` module.  For example:

   ::

      from phcpy.solutions import strsol2dict
      strsol2dict(solutions[0])

   shows the output

   ::

       {'t': 0j, 'm': 2, 'err': 1.17e-16, 'rco': 0.0, 'res': 0.0, 'x': 0j, 'y': 0j}

6. Higher values of the verbose level ``vrblvl`` are mainly meant
   for debugging purposes as it should procedures are executed. 
   As the solving of a polynomial system could take a long time,
   the user can see which procedures are currently running
   if the solver appears to be stuck.
