Numerical Polynomials
=====================

If the input is wrong, then the output will be wrong as well.

symbols
-------

The solver computes with complex numbers and the imaginary unit
can be denoted by ``i`` or ``I``.  For example:

::

    p = 'x - i;'

which represents the polynomial ``x - i``.
While ``i`` and ``I`` are symbols, they should not be used as 
the names of variables, because ``i`` and ``I`` are interpreted as numbers.

The polynomials can be written in factored form.  For example:

::

    p = '(y - 2)*(x + y - 1)^2;'

To enter this polynomial, first, the number of variables must be set.

::

    from phcpy.dimension import set_double_dimension, get_double_dimension
    set_double_dimension(2)
    get_double_dimension()

which returns ``2``, confirming the number of variables.

::

    from phcpy.polynomials import set_double_polynomial, get_double_polynomial
    set_double_polynomial(1, 2, p)
    get_double_polynomial(1)

shows

::

    'y^3 + 2*y^2*x + y*x^2 - 4*y^2 - 6*y*x - 2*x^2 + 5*y + 4*x - 2;'

Observe that the variable ``y`` appears first in the polynomial.

::

    from phcpy.polynomials import string_of_symbols
    string_of_symbols()

confirms the list of symbols used as variable names:

::

    ['y', 'x']

If that is inconvenient, the simple trick is to add 
the null polynomial ``x - x`` in front of the string
representation to ensure ``x`` comes first.

::

    q = 'x - x + ' + p
    q

defines the string ``'x - x + (y - 2)*(x + y - 1)^2;'``

To continue, the polynomial that has been set must be cleared.

::

    from phcpy.polynomials import clear_double_system
    clear_double_system()

::

    set_double_dimension(2)
    set_double_polynomial(1, 2, q)
    get_double_polynomial(1)

which shows

::

    'x^2*y + 2*x*y^2 + y^3 - 2*x^2 - 6*x*y - 4*y^2 + 4*x + 5*y - 2;'

Polynomials added to the system will follow the same order of symbols.

::

    r = 'y**2 + 4*x + y - 1;'
    set_double_polynomial(2, 2, r)
    get_double_polynomial(2)

confirms

::

    'y^2 + 4*x + y - 1;'

Consider the entire system of polynomials:

::

    from phcpy.polynomials import get_double_system
    get_double_system()

returns the list

::

    ['x^2*y + 2*x*y^2 + y^3 - 2*x^2 - 6*x*y - 4*y^2 + 4*x + 5*y - 2;',
     'y^2 + 4*x + y - 1;']

numbers
-------

The letters ``e`` and ``E`` appear in the scientific notation 
of floating-point numbers.  Therefore, these letters ``e`` and ``E``
should not be used as the names of variables.

Even as coefficients can be entered as exact rational numbers, 
they are evaluated to floating-point numbers.

::

    clear_double_system() # restart

::

    p = 'x - 1/3;'

Entering ``p`` goes then as follows:

::

    set_double_dimension(1)
    set_double_polynomial(1, 1, p)
    get_double_polynomial(1)

which then shows the decimal floating-point expansion of ``1/3``:

::

    ' + x - 3.33333333333333E-01;'

Let us double the precision, 
first importing the corresponding double double functions:

::

    from phcpy.dimension import set_double_double_dimension
    from phcpy.polynomials import set_double_double_polynomial

The statements

::

    from phcpy.polynomials import get_double_double_polynomial
    set_double_double_dimension(1)
    set_double_double_polynomial(1, 1, p)
    get_double_double_polynomial(1)"

then show

::

    ' + x-3.33333333333333333333333333333324E-1;'

Laurent polynomials
-------------------

Laurent polynomials have integer numbers as expeonents
and therefore admit negative values.
Some problems benefit from representations as Laurent polynomials.
Consider for example the cyclic 5-roots problem:

::

   from phcpy.families import cyclic
   for pol in cyclic(5): print(pol)

which shows

::


   x0 + x1 + x2 + x3 + x4;
   x0*x1 + x1*x2 + x2*x3 + x3*x4 + x4*x0;
   x0*x1*x2 + x1*x2*x3 + x2*x3*x4 + x3*x4*x0 + x4*x0*x1;
   x0*x1*x2*x3 + x1*x2*x3*x4 + x2*x3*x4*x0 + x3*x4*x0*x1 + x4*x0*x1*x2;
   x0*x1*x2*x3*x4 - 1;

The cyclic n-roots problem reformulated as a Laurent system,
appears in a paper of Uffe Haagerup, available as

::

   from phcpy.families import cyclic_reformulated
   for pol in cyclic_reformulated(5): print(pol)

which shows

::


   x1 + x1^-1*x2 + x2^-1*x3 + x3^-1*x4 + x4^-1;
   x2 + x1^-1*x3 + x2^-1*x4 + x3^-1 + x4^-1*x1;
   x3 + x1^-1*x4 + x2^-1 + x3^-1*x1 + x4^-1*x2;
   x4 + x1^-1 + x2^-1*x1 + x3^-1*x2 + x4^-1*x3;

Compared to the original cyclic 5-root problem,
the Laurent version has a more regular structure.
