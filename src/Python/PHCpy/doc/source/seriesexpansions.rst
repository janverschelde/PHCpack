Power Series Expansions
=======================

A polynomial homotopy defines algebraic curves.  With Newton's method, 
power series expansions of the algebraic can be computed.
We distinguish between

1. Taylor series that start at a regular point; and

2. power series starting at leading terms of a series.

The functions used in this section are imported below.

::

    from phcpy.solutions import make_solution
    from phcpy.series import double_newton_at_point
    from phcpy.series import double_newton_at_series
    from phcpy.series import double_pade_approximants

Taylor series and Pade approximants
-----------------------------------

The function

.. math::

   f(z) = \sqrt{\frac{1 + z/2}{1 + 2 z}}

is a solution :math:`x(s)` of the homotopy

.. math::

   (1-s)(x^2 - 1) + s(3x^2 - 3/2) = 0.

which is defined by the ``pol`` below:

::

   pol = ['(x^2 - 1)*(1-s) + (3*x^2 - 3/2)*s;']

At ``s`` equal to zero, the values for ``x`` are :math:`\pm 1`.
The start solutions are defined as

::

    variables = ['x', 's']
    sol1 = make_solution(variables, [1, 0])
    sol2 = make_solution(variables, [-1, 0])
    sols = [sol1, sol2]

It is always good to double check by printing the solutions ``sols``:

::

    for (idx, sol) in enumerate(sols):
        print('Solution', idx+1, ':')
        print(sol)

which shows

::

    Solution 1 :
    t : 0.000000000000000E+00 0.000000000000000E+00
    m : 1
    the solution for t :
     x : 1.000000000000000E+00  0.0
     s : 0.000000000000000E+00  0.0
    == err : 0.000E+00 = rco : 1.000E+00 = res : 0.000E+00 =
    Solution 2 :
    t : 0.000000000000000E+00 0.000000000000000E+00
    m : 1
    the solution for t :
     x : -1.000000000000000E+00  0.0
     s : 0.000000000000000E+00  0.0
    == err : 0.000E+00 = rco : 1.000E+00 = res : 0.000E+00 =

In calling ``double_newton_at_point`` we must declare the index 
of the variable which serves as the parameter.
As the parameter ``s`` is the second parameter, we set ``idx = 2``.
Other parameters are the degree of the series, 
by default ``maxdeg=4`` and the number of Newton steps,
by default ``nbr=4``.

::

    srs = double_newton_at_point(pol, sols, idx=2)
    srs

which returns

::

    [[' + 3.69287109375000E+00*s^4 - 2.08593750000000E+00*s^3 + 1.21875000000000E+00*s^2 - 7.50000000000000E-01*s + 1;'],
     [' - 3.69287109375000E+00*s^4 + 2.08593750000000E+00*s^3 - 1.21875000000000E+00*s^2 + 7.50000000000000E-01*s - 1;']]

The power series are polynomials which should be read
from right to left.

Pade approximants are rational expressions which agree 
with the power series expansions.

::

    pad = double_pade_approximants(pol, sols, idx=2)
    pad

which shows

::

    [['(9.53124999999999E-01*s^2 + 2.12500000000000E+00*s + 1)/(1.89062500000000E+00*s^2 + 2.87500000000000E+00*s + 1)'],
     ['( - 9.53124999999999E-01*s^2 - 2.12500000000000E+00*s - 1)/(1.89062500000000E+00*s^2 + 2.87500000000000E+00*s + 1)']]

To evaluate the first Pade approximant at ``0.1`` as the value
for ``s``, we first replace the string ``"s"`` by the string ``"0.1"``.

::

    p0 = pad[0][0].replace('s', '0.1')
    p0

which leads to the string

::

    '(9.53124999999999E-01*0.1^2 + 2.12500000000000E+00*0.1 + 1)/(1.89062500000000E+00*0.1^2 + 2.87500000000000E+00*0.1 + 1)'

Observe the ``^`` which must be replaced by ``**``, as done by

::

    p1 = p0.replace('^', '**')
    p1

and then we have the expression in the string ``p1``:

::

    '(9.53124999999999E-01*0.1**2 + 2.12500000000000E+00*0.1 + 1)/(1.89062500000000E+00*0.1**2 + 2.87500000000000E+00*0.1 + 1)'

Its numerical value, obtained via 

::

    ep1 = eval(p1)
    ep1

shows

::

    0.9354144241119484

Let us compare this now to the function

.. math::

   f(z) = \sqrt{\frac{1 + z/2}{1 + 2 z}}

evaluated at :math:`z = 0.1`.  The statements

::

    from math import sqrt
    ef1 = sqrt((1 + 0.1/2)/(1+2*0.1))
    ef1

produce the value

::

    0.9354143466934854


and we obtain ``7.741846297371069e-08`` as the outcome of
the statement ``abs(ep1 - ef1)``.
The error shows we have about 8 decimal places correct 
for the value of :math:`f(z)` at :math:`z = 0.1`.

expansions starting at series
-----------------------------

Starting the power series expansions at a series allows to start 
at a singular solution, as illustrated by the Viviani curve,
defines as the intersection of a sphere with a quadratic cylinder.

::

    pols = [ '2*t^2 - x;', \
             'x^2 + y^2 + z^2 - 4;' , \
             '(x-1)^2 + y^2 - 1;']

The series starts at :math:`x = 2 t^2`.

::

    lser = [ '2*t^2;', '2*t;', '2;']

Then Newton's method is executed, here using the default ``idx=1``
as ``t`` is the first parameter, asking for a series truncated
at degree 12, using no more than 8 iterations.

::

    nser = double_newton_at_series(pols, lser, maxdeg=12, nbr=8)

To print the expansions in ``nser``, the names
of the variables are used:

::

    variables = ['x', 'y', 'z']
    for (var, pol) in zip(variables, nser):
        print(var, '=', pol)

::

    x = 2*t^2;
    y =  - 5.46875000000000E-02*t^11 - 7.81250000000000E-02*t^9 - 1.25000000000000E-01*t^7 - 2.50000000000000E-01*t^5 - t^3 + 2*t;
    z =  - 4.10156250000000E-02*t^12 - 5.46875000000000E-02*t^10 - 7.81250000000000E-02*t^8 - 1.25000000000000E-01*t^6 - 2.50000000000000E-01*t^4 - t^2 + 2;

The coefficients of the power series expansions indicate how fast 
the solutions change once we move away from the singularity.

The example below compares the series expansions at two solutions
for the problem of Apollonius."

::

   pols = [ 'x1^2 + 3*x2^2 - r^2 - 2*r - 1;', \
            'x1^2 + 3*x2^2 - r^2 - 4*x1 - 2*r + 3;', \
          '3*t^2 + x1^2 - 6*t*x2 + 3*x2^2 - r^2 + 6*t - 2*x1 - 6*x2 + 2*r + 3;']

The ``pols`` define an instance where the input circles are mutually
touching each other.  Once we start moving the input circles apart
from each other, how fast do the touching circles grow?
The developments start at the following terms:

::

    lser1 = ['1;', '1 + 0.536*t;', '1 + 0.904*t;']
    lser2 = ['1;', '1 + 7.464*t;', '1 + 11.196*t;']

We run Newton's method twice:

::

    nser1 = double_newton_at_series(pols, lser1, idx=4, nbr=7)
    nser2 = double_newton_at_series(pols, lser2, idx=4, nbr=7)

The statements

::

    variables = ['x', 'y', 'z']
    print('the first solution series :')
    for (var, pol) in zip(variables, nser1):
        print(var, '=', pol)
    print('the second solution series :')
    for (var, pol) in zip(variables, nser2):
        print(var, '=', pol)

have as output

::

    the first solution series :
    x =  - 4.03896783473158E-28*t^4 - 4.62223186652937E-33*t + 1;
    y =  - 7.73216421430735E-03*t^4 + 7.73216421430735E-03*t^3 - 1.66604983954048E-02*t^2 + 5.35898384862246E-01*t + 1;
    z =  - 1.33925012716462E-02*t^3 + 2.88568297002609E-02*t^2 + 8.03847577293368E-01*t + 1;
    the second solution series :
    x = 1.00974195868290E-28*t^3 + 1.97215226305253E-31*t + 1;
    y =  - 2.90992267835785E+02*t^4 + 2.90992267835785E+02*t^3 + 4.50166604983953E+01*t^2 + 7.46410161513775E+00*t + 1.00000000000000E+00;
    z = 5.04013392501271E+02*t^3 + 7.79711431702996E+01*t^2 + 1.11961524227066E+01*t + 1.00000000000000E+00;

Observe the difference in magnitudes of the coefficients
of the series expansions, indicating that one solution will
change more than the other, as we move away from the singularity.
