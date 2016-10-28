****************
Reference Manual
****************

This chapter contains the documentation of the modules in the package,
mostly automatically generated from the documentation strings in the
module and of the functions exported by each module.
The order of the sections in this chapter follows the order of
the previous chapters.  The first section on the function in the
solver module corresponds with chapter 3, on a blackbox solver.
The section headings correspond to the earlier chapter headings.

a blackbox solver for isolated solutions
========================================

The two most essential modules to solve polynomial systems
are the ``solver`` module, which exports the blackbox solver,
and the ``solutions`` module, to parse the computed solutions.

functions in the module solver
------------------------------

The documentation strings of the functions
exported by the module ``solver`` of the package phcpy are listed below.

.. automodule:: solver
   :members:

functions in the module solutions
---------------------------------

The documentation strings of the functions
exported by the module ``solutions`` are listed below.
The script **test()** runs when typing **python solutions.py**
at the command prompt.

.. automodule:: solutions
   :members:

path trackers and sweep homotopies
==================================

functions in the module trackers
--------------------------------
   
The documentation strings of the functions
exported by the module ``trackers`` are listed below.

.. automodule:: trackers
   :members:

functions in the module sweepers
--------------------------------
   
The documentation strings of the functions
exported by the module ``sweepers`` are listed below.

.. automodule:: sweepers
   :members:

functions in the module tuning
------------------------------

The documentation strings of the functions
exported by the module ``tuning`` are listed below.

.. automodule:: tuning
   :members:

functions in the module tropisms
--------------------------------

The module ``tropisms`` provides access to the numerically computed
tropisms via a polyhedral end game. 
The functions exported by the module ``tropisms`` are listed below.

.. automodule:: tropisms
   :members:

positive dimensional solution sets
==================================

functions in the module sets
----------------------------

.. automodule:: sets
   :members:

functions in the module examples
--------------------------------

The documentation strings of the functions that
return the polynomials of the example systems as strings of characters
are listed below.
The regression test is exported by the function **test()**
of the module ``examples``.

.. automodule:: examples
   :members:

functions in the module families
--------------------------------

.. automodule:: families
   :members:

numerical Schubert calculus
===========================

functions in the module schubert
--------------------------------

.. automodule:: schubert
   :members:

Newton polytopes and monomial maps
==================================

functions in the module polytopes
---------------------------------

.. automodule:: polytopes
   :members:

functions in the module maps
----------------------------

.. automodule:: maps
   :members:

a graphical user interface
==========================

functions in the module dashboard
---------------------------------

.. automodule:: dashboard
   :members:

the module phcpy.phcpy2c2
=========================

Almost all computations in phcpy are done with compiled code,
provided in one object.

functions in the module interface
---------------------------------

Polynomial systems and solutions are passed through as strings.
The coefficients of the polynomials and coordinates of the solutions
are evaluated in standard double, double double, quad double precision,
or arbitrary multiprecision.

.. automodule:: interface
   :members:

functions in the module phcpy2c2
--------------------------------

The module ``phcpy2c2`` wraps the C functions in the C interface to PHCpack.
The C interface to PHCpack was developed in the application of message passing
(MPI) to run the path trackers on distributed memory multiprocessing computers.
All functions documented below have their counterpart in C 
that are therefore then also directly accessible from C programs.

.. automodule:: phcpy2c2
   :members:
