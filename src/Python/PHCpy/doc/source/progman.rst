*****************
Programmer Manual
*****************

The code is spread over more than twenty modules.
The directed acyclic graph shows the dependencies of
the first seven modules:

::

    version
       |
       +----------------------------> solutions         
       |                                  |
       +--> dimension                     |
       |        |                         |
       +--------+--> polynomials          |
       |        |         |               |
       +------------------+--> volumes    |
       |        |         |       |       |
       +--------+---------+-------+-------+--> solver 
       |        |                                |
       +--------+--------------------------------+--> examples

At the root is the version module.  If version.py works,
then the interfacing with libPHCpack works.
Four other modules are needed for the blackbox solver.

In addition to functions, each module has test functions,
starting with ``test_``.  Each test function returns the
number of failed tests.  The ``main()`` of each module
runs all tests.

the version module
==================

The PHCpack version string is retrieved via this module,
along with some utility functions for the basic type conversions.

.. automodule:: version
   :members:

a blackbox solver for isolated solutions
========================================

The three most essential modules to solve polynomial systems
are the ``solver`` module, which exports the blackbox solver,
the ``polynomials`` module, to define the input,
and the ``solutions`` module, to parse the computed solutions.

functions in the module dimension
---------------------------------

In addition to the number of polynomials in the system,
the dimension module controls some other important numbers 
for each run: the seed of the random number generators
and the number of available cores.

.. automodule:: dimension
   :members:

functions in the module polynomials
-----------------------------------

The module polynomials exports the definition of a class
to represent systems of polynomials.

.. automodule:: polynomials
   :members:

functions in the module solutions
---------------------------------

The documentation strings of the functions
exported by the module ``solutions`` are listed below.

.. automodule:: solutions
   :members:

functions in the module volumes
-------------------------------

The mixed volume of the tuple of Newton polytopes of
a polynomial systems is a generically sharp root count
on the number of isolated solutions with nonzero coordinates.
For all affine solutions (also counting solutions with zero
coordinates), the stable mixed volume provides a generically
sharp root count.

.. automodule:: volumes
   :members:

functions in the module solver
------------------------------

The documentation strings of the functions
exported by the module ``solver`` of the package phcpy are listed below.

.. automodule:: solver
   :members:

functions in the module examples
--------------------------------

To demonstrate the capabilities of the solver
and the relevance of polynomial systems to various fields
of science and engineering, several examples from the literature
are provided by the ``examples`` module.

.. automodule:: examples
   :members:

functions in the module families
--------------------------------

Polynomial system are often defined for any dimension
that is: for any number of equations and number of variables.
Such families of polynomial systems are important benchmark systems,
for example: the cyclic n-roots system.

.. automodule:: families
   :members:

homotopy methods and path tracking algorithms
=============================================

A homotopy method constructs a start system for
an artificial parameter homotopy.  Solution paths
start at the known solutions of the start system
and end at the solutions of the system that has
to be solved, called the target system.

Path tracking algorithms compute approximations for
the points on the solution paths, 
using predictor-corrector methods.
The step size control can be done either *aposteriori*,
based on the performance of the corrector; 
or *apriori*, using a predictor which estimates the
distance to the nearest path and the convergence radius
of the power series expansions of the solution curves.
For natural parameter homotopies, arc length parameter
continuation with aposteriori step size control is available.

functions in the module homotopies
----------------------------------

.. automodule:: homotopies
   :members:

functions in the module starters
--------------------------------

.. automodule:: starters
   :members:

functions in the module trackers
--------------------------------

.. automodule:: trackers
   :members:

functions in the module tropisms
--------------------------------

.. automodule:: tropisms
   :members:

functions in the module sweepers
--------------------------------

.. automodule:: sweepers
   :members:

functions in the module series
------------------------------

.. automodule:: series
   :members:

functions in the module curves
------------------------------

.. automodule:: curves
   :members:

functions in the module deflation
---------------------------------

By random choices of complex constants in an artificial parameter homotopy,
singularities along a path are avoided.  At the end of solution paths,
we may encounter isolated singular solutions.
The method of deflation is often effective at accurately computing
isolated singular solutions.

.. automodule:: deflation
   :members:

homotopies for enumerative geometry
===================================

Homotopy methods to two types of problems in enumerative geometry
are provided, applying Pieri and Littlewood-Richardson root counts.

functions in the module schubert
--------------------------------

.. automodule:: schubert
   :members:

numerical irreducible decomposition
===================================

The second blackbox solver computes a numerical irreducible decomposition
of the solution set of a polynomial system.

functions in the module sets
----------------------------

.. automodule:: sets
   :members:

functions in the module cascades
--------------------------------

.. automodule:: cascades
   :members:

functions in the module diagonal
--------------------------------

.. automodule:: diagonal
   :members:

functions in the module factor
------------------------------

.. automodule:: factor
   :members:

functions in the module decomposition
-------------------------------------

.. automodule:: decomposition
   :members:

functions in the module binomials
---------------------------------

.. automodule:: binomials
   :members:
