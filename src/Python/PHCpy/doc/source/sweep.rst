sweeping homotopies
===================

A sweep homotopy is a family of polynomial systems with a least one
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

Geometrically, as the horizontal line moves up, the two solutions,
the intersection points on the circle and the line, move closer
to each other to join at a *quadratic turning point*.

.. image:: ./circles.png

functions in the module
-----------------------
   
The documentation strings of the functions
exported by the module ``sweepers`` of the package phcpy are listed below.

.. automodule:: sweepers
   :members:
