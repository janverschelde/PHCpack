Root Counts and Start Systems
=============================

A formal root count can be viewed as a count on the number of solutions 
of a system with a particular structure.  
The structure can be determined by the degrees or the Newton polytopes.

Each of the four root counts below is illustrated by a proper example.

total degree
------------

The *total degree* is the product of the degrees of all polynomials
in the system.

::

   from phcpy.starters import total_degree
   from phcpy.starters import total_degree_start_system

One family of polynomial systems for which the total degree equals
the number of solutions is the Katsura benchmark problem.

::

   from phcpy.families import katsura
   k3 = katsura(3)
   for pol in k3:
       print(pol)

Our example consists of one linear and three quadrics.

::

   u3 + u2 + u1 + u0 + u1 + u2 + u3 - 1;
   u3*u3 + u2*u2 + u1*u1 + u0*u0 + u1*u1 + u2*u2 + u3*u3 - u0;
   u3*u3 + u2*u3 + u1*u2 + u0*u1 + u1*u0 + u2*u1 + u3*u2 - u1;
   u3*u3 + u2 + u1*u3 + u0*u2 + u1*u1 + u2*u0 + u3*u1 - u2;

The total degree is ``8``, computed with

::

   degk3 = total_degree(k3)

and the corresponding start system has also ``8`` solutions,
as can be seen from the output of the code cell below.

::

   q, qsols = total_degree_start_system(k3)
   len(qsols)

The polynomials in the start system ``q`` are shown as

::

   u3^1 - 1;
   u2^2 - 1;
   u1^2 - 1;
   u0^2 - 1;

as output of

::

   for pol in q:
       print(pol)

multihomogeneous Bezout numbers
-------------------------------

Most polynomial systems arising in application have far fewer
solutions than the total degree.

::

   from phcpy.starters import m_homogeneous_bezout_number
   from phcpy.starters import m_homogeneous_start_system"

We illustrate m-homogenous Bezout numbers
with an application from game theory.

::

   from phcpy.families import generic_nash_system
   game4two = generic_nash_system(4)
   for pol in game4two:
       print(pol)"

The (omitted) output of the code cell above shows four cubics.
The output of

::

   mbn = m_homogeneous_bezout_number(game4two)
   mbn

is the tuple

::

   (9, '{ p2 }{ p3 }{ p4 }{ p1 }')

The tuple on return contains first the root count, 
and then the partition of the variables. 
Observe the difference with the total degree, which is ``81``.

To construct a start system with ``9`` solutions,
respecting the structure of the ``game4two`` system, we do:

::

   q, qsols = m_homogeneous_start_system(game4two, partition=mbn[1])
   len(qsols)

and we see ``9`` as the output of ``len(qsols)``.

linear-product start systems
----------------------------

The multihomogeneous Bezout numbers are computed based 
on a partition of the unknowns.
This partition is the same for all polynomials in the system.
A sharper bound can be obtained if this restriction is relaxed.

::

   from phcpy.starters import linear_product_root_count
   from phcpy.starters import random_linear_product_system

The example we use to illustrate this root count
consists of cubics.

::

    from phcpy.families import noon
    n3 = noon(3)
    for pol in n3:
        print(pol)

In dimension three, the system is then

::

    x1*(x2^2 + x3^2) - 1.1*x1 + 1;
    x2*(x1^2 + x3^2) - 1.1*x2 + 1;
    x3*(x1^2 + x2^2) - 1.1*x3 + 1;

The linear-product root count for this system
is computed by the instructions in the cell below:

::

    lprc = linear_product_root_count(n3)
    lprc

which gives as output the tuple

::

    (21,
    '{ x1 }{ x2 x3 }{ x2 x3 };{ x1 x3 }{ x1 x3 }{ x2 };{ x1 x2 }{ x1 x2 }{ x3 };')

The tuple on return contains first the upper bound on the number 
of solutions, followed by the set structure used to compute this bound.
Every set corresponds to a linear equation 
in the linear-product start system.

The start system is constructed and solved by the following:

::

    q, qsols = random_linear_product_system(n3, lprc[1])
    len(qsols)

and ``len(qsols)`` returns ``21``.

mixed volumes
-------------

For sparse polynomial systems, that are systems with relatively few
monomials appearing with nonzero coefficients, the mixed volume of
the Newton polytopes provides a much sharper bound than any of
the degree bounds.

::

   from phcpy.volumes import mixed_volume
   from phcpy.volumes import make_random_coefficient_system

The cyclic n-roots problem illustrates the need to
apply mixed volumes very well.

::

   from phcpy.families import cyclic
   c5 = cyclic(5)
   for pol in c5:
       print(pol)

which shows

::

   x0 + x1 + x2 + x3 + x4;
   x0*x1 + x1*x2 + x2*x3 + x3*x4 + x4*x0;
   x0*x1*x2 + x1*x2*x3 + x2*x3*x4 + x3*x4*x0 + x4*x0*x1;
   x0*x1*x2*x3 + x1*x2*x3*x4 + x2*x3*x4*x0 + x3*x4*x0*x1 + x4*x0*x1*x2;
   x0*x1*x2*x3*x4 - 1;

The mixed volume equals ``70`` for this system,
and is computed via

::

   mv = mixed_volume(c5)

A *random coefficient system* has the same monomial structure
as the input system, but random coefficients.
The start system is made and solved with the code below

::

   vol, q, qsols = make_random_coefficient_system(c5)
   len(qsols)

and we see ``70`` as the output of ``len(qsols)``.

The mixed volume does not count the solutions with zero coordinates.
To count all solutions, also those with zero coordinates, 
the stable mixed volume should be computed.

::

   from phcpy.volumes import stable_mixed_volume

Consider the following example.

::

   pols = ['x^3 + x*y^2 - x^2;', 'x + y - y^3;']

The command

::

   stable_mixed_volume(pols)

returns a tuple of two integers.
The first number in the tuple is the mixed volume, 
the second number is the stable mixed volume,
which takes into account the solutions with zero coordinates.
