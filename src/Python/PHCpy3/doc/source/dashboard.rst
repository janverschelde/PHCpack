a graphical user interface
==========================

While phcpy is oriented towards command line interactive computing
and running scripts, Tkinter provides tools for building
a graphic user interface.  The goal of the module **dashboard**
is to develop a graphical user interface to the methods of phcpy.

solving with a click of a button
--------------------------------

A very basic graphical user interface to the blackbox solver consists
of two text widgets: one for the input and another for the output;
one button for the user to call the blackbox solver, and then two labels
to document the functionality of the text widgets.

A screen shot of a basic interface to the blackbox solver is shown
in the figure with caption :ref:`figlaunchsolver`

.. _figlaunchsolver:

.. figure:: ./launchsolver.png
    :align: center

    Solving cyclic 5-roots with a click of a button.

The code to launch this GUI is as follows.

::

    >>> from phcpy.dashboard import launchsolver
    >>> from phcpy.families import cyclic
    >>> launchsolver(cyclic(5))

If called without arguments, as ``launchsolver()``,
then the input text widget is empty and the user must enter
the polynomials in the system.

scrolling a list of solutions
-----------------------------

The blackbox solver **solve** of **phcpy.solver** returns a list of strings.
In the command line mode, we can print the solutions, one after the other.
The ``scrollsols`` function launches a simple interface to scroll through
the list of solutions, by clicking on ``previous`` or ``next`` buttons.
The session below illustrates the scrolling through the solutions of
the cyclic 5-roots problem.

::

    >>> from phcpy.families import cyclic
    >>> c5 = cyclic(5)
    >>> from phcpy.solver import solve
    >>> sols = solve(c5, silent=True)
    >>> from phcpy.dashboard import scrollsols
    >>> scrollsols(sols)

making a coordinate plot of solutions
-------------------------------------

Solutions have coordinates in the complex plane.
As in the case of the cyclic 5-roots problem,
a plot of one of the coordinates in the complex plane
reveals the pattern of the distribution in the roots,
see the figure with caption :ref:`figc5x0plot`

.. _figc5x0plot:

.. figure:: ./c5x0plot.png
    :align: center

    The first coordinate of cyclic 5-roots in the complex plane.


The plot appears in a canvas widget, in the GUI launched
by the function ``plotcoordinate(sols, idx)`` where ``sols``
is the list of solutions and ``idx`` an index to a coordinate
of the solutions.

functions in the module dashboard
---------------------------------

.. automodule:: dashboard
   :members:
