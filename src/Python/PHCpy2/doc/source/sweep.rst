sweep homotopies
================

A *sweep homotopy* is a family of polynomial systems with a least one
natural parameter and one artificial parameter.
As the artificial parameter moves from zero to one, the natural parameter
changes from a given start value to another given target value.
By arc length continuation, the solution paths are tracked from the
given start values for the parameters to the target values.

a simple example: sweeping the circle
-------------------------------------

We consider the unit circle \ :math:`x^2 + y^2 - 1 = 0`,
intersected by a horizontal line, at the start equal to \ :math:`y = 0`.
In a Python session, we could define the sweep homotopy that takes
the line from \ :math:`y = 0` to \ :math:`y = 2`.

::

   >>> circle = ['x^2 + y^2 - 1;', 'y*(1-s) + (y-2)*s;']

For ``s = 0`` there are two solutions,
with values for ``x`` and ``y`` in the 
tuples \ :math:`(1, 0)` and \ :math:`(-1,0)`.

Geometrically, as the horizontal line moves up, the two solutions
(the intersection points on the circle and the line), move closer
to each other to join at a *quadratic turning point*.
At the left of the picture below we see the line transversally intersecting
the circle at a perfect right angle.  At the right, the two distinct
solutions have merged into one point where the line is tangent to the circle.

.. image:: ./circles.png

The tracking of solution paths in a real sweep homotopy will stop
at the first singular point it encounters.  
The continuation of the code with the definition of ``circle``
to launch this path tracking is listed below:

::

   >>> from phcpy.solutions import make_solution as makesol
   >>> first = makesol(['x', 'y', 's'], [1, 0, 0])
   >>> second = makesol(['x', 'y', 's'], [-1, 0, 0])
   >>> startsols = [first, second]
   >>> from phcpy.sweepers import standard_real_sweep as sweep
   >>> newsols = sweep(circle, startsols)
   >>> print newsols[0]

and then we see as output of the ``print`` statement:

::

   t :  0.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -2.46519032881566E-32   0.00000000000000E+00
    y :  1.00000000000000E+00   0.00000000000000E+00
    s :  5.00000000000000E-01   0.00000000000000E+00
   == err :  0.000E+00 = rco :  1.000E+00 = res :  0.000E+00 =

The sweep stopped where ``s`` is equal to 0.5,
with corresponding values for ``x`` and ``y`` in the 
tuple \ :math:`(0, 1)`.

real versus complex sweeps
--------------------------

In a *complex sweep*, an addition random gamma constant is generated
in the convex-linear combination between the sets of start and target
values for the parameters.  If the solutions for the start values of
the parameters are regular, then the application of the *gamma trick*
applies for problems where the parameter space is convex.
This means that, if the problem formulation makes sense for convex
combinations of the parameters, then the solution paths will remain
regular, except for finitely many bad choices of the random gamma constant,
and except perhaps at the very end of the paths, when the target values
for the parameters lead to polynomial systems with singular solutions.

Conducting a complex sweep on the circle can be done as follows:

::

    >>> circle = ['x^2 + y^2 - 1;']
    >>> from phcpy.solutions import make_solution as makesol
    >>> first = makesol(['x', 'y'], [1, 0])
    >>> second = makesol(['x', 'y'], [-1, 0])
    >>> startsols = [first, second]
    >>> par = ['y']
    >>> start = [0, 0] 
    >>> target = [2, 0]
    >>> from phcpy.sweepers import standard_complex_sweep as sweep
    >>> newsols = sweep(circle, startsols, 2, par, start, target)

The setup of the homotopy defines ``y`` as the parameter 
(in the list ``['y']`` assigned to ``par``).
The parameter \ :math:`y` will move from the complex zero \ :math:`0 + 0 I` 
(given by the list ``[0, 0]`` assigned to ``start``) 
to \ :math:`2 + 0 I`
(given by the list ``[2, 0]`` assigned to ``target``).
The corresponding start solutions for \ :math:`y = 0` are stored
in the tuples \ :math:`(1,0)` and \ :math:`(-1,0)`.
Then, at the end of the sweep, we will find two complex conjugated solutions.

::

    >>> print newsols[0]
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x : -6.27554627321419E-26   1.73205080756888E+00
     y :  2.00000000000000E+00   0.00000000000000E+00
    == err :  6.642E-13 = rco :  1.000E+00 = res :  4.441E-16 =

What is now the difference between real versus complex?
The *real sweep* stopped at the singular solution \ :math:`(0,1)`
while the *complex sweep* hopped over this singularity because
of complex random gamma constant in the convex combination between
the start and target values of the parameters.

functions in the module
-----------------------
   
The documentation strings of the functions
exported by the module ``sweepers`` of the package phcpy are listed below.

.. automodule:: sweepers
   :members:
