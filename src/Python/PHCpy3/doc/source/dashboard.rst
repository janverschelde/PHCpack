a graphical user interface
==========================

While phcpy is oriented towards command line interactive computing
and running scripts, the Tkinter module provides tools for building
a graphic user interface.  The goal of the module **dashboard**
is to develop a graphical user interface to the methods of phcpy.

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

A very basic graphical user interface to the blackbox solver consists
of two text widgets: one for the input and another for the output;
one button for the user to call the blackbox solver, and then two labels
to document the functionality of the text widgets.

A screen shot of a basic interface to the blackbox solver is shown below.

.. image:: ./launchsolver.png

The code to launch this GUI is as follows.

::

    >>> from phcpy.dashboard import launchsolver
    >>> from phcpy.families import cyclic
    >>> launchsolver(cyclic(5))

If called without arguments, as ``launchsolver()``,
then the input text widget is empty and the user must enter
the polynomials in the system.

functions in the module dashboard
---------------------------------

.. automodule:: dashboard
   :members:
