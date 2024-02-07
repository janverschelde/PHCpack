Numerical Solutions
===================

Solutions of are numerical and returned as lists of PHCpack
solution strings.  The solutions module contains functions to
parse a PHCpack solution string into a dictionary.

The solutions module exports operations

1. to parse strings in the PHCpack solution format into dictionaries;

2. to evaluate these dictionaries into polynomials substituting the
   values for the variables into the strings representing the polynomials.

Another useful operation is the ``verify`` function, 
to evaluate the polynomials at solutions.

attributes of numerical solutions
---------------------------------

The information of a solution as a dictionary contains the following:

1. `t` : value of the continuation parameter

2. `m` : multiplicity of the solution

3. symbols for the variables are keys in the dictionary,
   the corresponding values are complex floating-point numbers

4. `err` : magnitude of the last correction term of Newton's method
   (forward error)

   `rco` : estimate for the inverse of the condition number of
   the Jacobian matrix at the solution

   `res` : magnitude of the residual (backward error)

The triplet (`err`, `rco`, `res`) measures
the numerical quality of the solution.
The residual `res` is normally interpreted as an estimate of the backward
error: by how much should we change the original problem such that the
approximate solution becomes an exact solution of the changed problem?
The estimate `rco` gives a (sometimes too pessimistic) bound on the
number of correct decimal places in the approximate solution.
In particular: ``abs(log(rco, 10))`` bounds the number of lost decimal
places in the approximate solution.
For example, if `rco` equals ``1.0E-8``, then the last 8 decimal places
in the coordinates of the solution could be wrong.

The best numerically conditioned linear systems arise when the
normals to the coefficient vectors of the linear equations are
perpendicular to each other, as in the next session:

::

    from phcpy.solver import solve
    p = ['x + y - 1;', 'x - y - 1;']
    s = solve(p)
    print(s[0])

with the output of the print statement:

::

    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  1.00000000000000E+00   0.00000000000000E+00
     y :  0.00000000000000E+00  -0.00000000000000E+00
    == err :  2.220E-16 = rco :  5.000E-01 = res :  0.000E+00 =

The value of `rco` is ``5.0E-1`` which implies that the
condition number is bounded by 2, as `rco` is an estimate
for the inverse of the condition number.
Roundoff errors are doubled at most.

At the opposite end of the best numerically conditioned linear systems
are those where the the normals to the coefficient vectors of the
linear equations are almost parallel to each other,
as illustrated in the next example:

::

    p = ['x + y - 1;', 'x + 0.999*y - 1;']
    s = solve(p)
    print(s[0])

and the output of the above code cell is

::

    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  1.00000000000000E+00   0.00000000000000E+00
     y :  0.00000000000000E+00  -0.00000000000000E+00
    == err :  2.220E-16 = rco :  2.501E-04 = res :  0.000E+00 =

The reported estimate for the inverse of the condition number
`rco` is ``2.5E-4``, which implies that the condition number is
estimated at 4,000.  Thus for this example, roundoff errors
may magnify thousandfold.  In the next example, the condition
number becomes a 10-digit number:

::

    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  1.00000000000000E+00   0.00000000000000E+00
     y :  0.00000000000000E+00  -0.00000000000000E+00
    == err :  2.220E-16 = rco :  2.500E-10 = res :  0.000E+00 =

which is produced by the code in the cell below:

::

    p = ['x + y - 1;', 'x + 0.999999999*y - 1;']
    s = solve(p)
    print(s[0])

Observe that the actual value of the solution remains ``(1,0)``,
which on the one hand indicates that the condition number is
a pessimistic bound on the accuracy of the solution.
But on the other hand, (1,0) may give the false security that
the solution is right, because the problem on input is very close
to a linear system which has infinitely many solutions
(the line ``x + y - 1 = 0``) and not the isolated point ``(1,0)``.

For a solution of the example ``noon`` from the module ``families``,
we convert the PHCpack format solution string to a dictionary as follows:

::

    from phcpy.families import noon
    s = solve(noon(3))
    print(s[0])

which shows

::

    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x1 : -6.77804511269800E-01   5.27500584353303E-01
     x2 :  1.35560902253960E+00   2.32882178444166E-17
     x3 : -6.77804511269800E-01  -5.27500584353303E-01
    == err :  1.601E-16 = rco :  2.303E-01 = res :  4.996E-16 =

To turn this string into a dictionary, do the following:

::

    from phcpy.solutions import strsol2dict
    d = strsol2dict(s[0])
    d.keys()

which shows

::

    dict_keys(['t', 'm', 'err', 'rco', 'res', 'x1', 'x2', 'x3'])

To select the value of the ``x1`` coordinate, which is

::

    (-0.6778045112698+0.527500584353303j)

then just do ``d['x1']``.

Observe that the values of the dictionary ``d`` are evaluated strings,
parsed into Python objects.

By plain substitution of the values of the dictionary representation
of the solution into the string representation of the polynomial system
we can verify that the coordinates of the solution evaluate to numbers
close to the numerical working precision:

::

    from phcpy.solutions import evaluate
    e = evaluate(noon(3), d)
    for x in e: print(x)

shows

::

    (-2.886579864025407e-15+6.661338147750939e-16j)
    (-4.440892098500626e-16-1.475351643981535e-17j)
    (-2.886579864025407e-15-6.661338147750939e-16j)

The ``evaluate`` is applied in the ``verify`` which computes 
the sum of all evaluated polynomials, in absolute value,
summed over all solutions.

::

    from phcpy.solutions import verify
    err = verify(noon(3), s)

The number ``err`` can be abbreviated into ``2.2e-13`` 
which is close enough to zero.

filtering solution lists
------------------------

The module exports function to filter regular solutions, solutions
with zero coordinates or real solutions.  The filtering of real
solutions is illustrated in the session below.
We first define one real solution and another with a coordinate
that has a nonzero imaginary part.

::

    from phcpy.solutions import make_solution
    s0 = make_solution(['x', 'y'], [complex(1, 0), complex(0, 2)])
    print(s0)

shows

::

    t : 0.000000000000000E+00 0.000000000000000E+00
    m : 1\n",
    the solution for t :\n",
     x : 1.000000000000000E+00  0.000000000000000E+00
     y : 0.000000000000000E+00  2.000000000000000E+00
    == err : 0.000E+00 = rco : 1.000E+00 = res : 0.000E+00 =

and the output 

::

    t : 0.000000000000000E+00 0.000000000000000E+00
    m : 1
    the solution for t :
     x : 2.000000000000000E+00  0.0
     y : 3.000000000000000E+00  0.0
    == err : 0.000E+00 = rco : 1.000E+00 = res : 0.000E+00 =

is produced by the the statements

::

    s1 = make_solution(['x', 'y'], [float(2), float(3)])
    print(s1)

The filtering of real solutions (with respect to a given tolerance)
is provided by the functions ``is_real`` (on one solution)
and ``filter_real`` (on a list of solutions).

Observe the tolerance ``1.0e-8`` as the second argument
in the application of the ``is_real`` function.

::

    from phcpy.solutions import is_real, filter_real
    is_real(s0, 1.0e-8)

with respect to the tolerance ``1.0e-8``, ``is_real``
returns ``False``, as ``s0`` is not a real solution.
For ``s1``, ``is_real(s1, 1.0e-8)`` returns ``True``.

Putting ``[s0, s1]`` into a list, to illustrate the
selection of the real solutions, with

::

    realsols = filter_real([s0, s1], 1.0e-8, 'select')
    for sol in realsols: print(sol)

shows then indeed

::

    t : 0.000000000000000E+00 0.000000000000000E+00
    m : 1
    the solution for t :
     x : 2.000000000000000E+00  0.0
     y : 3.000000000000000E+00  0.0
    == err : 0.000E+00 = rco : 1.000E+00 = res : 0.000E+00 =

The functions ``filter_regular`` and ``filter_zero_coordinates``
to filter the regular solutions and those solutions with
zero coordinates respectively
operate in a manner similar as ``filter_real.``

Another application of ``make_solution`` is to turn the solution
at the end of path (with value 1.0 for ``t``) to a solution which
can serve at the start of another path (with value 0.0 for ``t``).
This is illustrated in the session below.
We start by solving a simple system.

::

    p = ['x**2 - 3*y + 1;', 'x*y - 3;']
    s = solve(p)
    print(s[0])

which shows

::

    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x : -9.60087560673590E-01   1.94043922153735E+00
     y : -6.14512082773443E-01  -1.24199437256077E+00
    == err :  3.317E-16 = rco :  2.770E-01 = res :  4.441E-16 =

Then we import the functions ``coordinates`` and ``make_solution``
of the module ``solutions``.

::

    from phcpy.solutions import coordinates
    (names, values) = coordinates(s[0])
    names"

shows 

::

    ['x', 'y']

and

:: 

    values

shows

::

    [(-0.96008756067359+1.94043922153735j), (-0.614512082773443-1.24199437256077j)]

With the ``names`` and the ``value`` 
we can reconstruct the solution string.

::

    s0 = make_solution(names, values)
    print(s0)

with output

::

    t : 0.000000000000000E+00 0.000000000000000E+00
    m : 1
    the solution for t :
     x : -9.600875606735900E-01  1.940439221537350E+00
     y : -6.145120827734430E-01  -1.241994372560770E+00
    == err : 0.000E+00 = rco : 1.000E+00 = res : 0.000E+00 =

Observe that also the diagnostics are set to the defaults.
