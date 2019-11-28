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

the module polynomials
----------------------

The module polynomials exports the definition of a class
to represent systems of polynomials.
The class ``Polynomials`` is intended as the starting point
of the object-oriented interface to the solver.

.. automodule:: polynomials
   :members:

path trackers and sweep homotopies
==================================

In the path tracking, we distinguish between paths defined by

1. homotopies with sufficiently random complex numbers,
   which then are free of singular solutions, except perhaps at the end; and

2. homotopies with real coefficients and real parameters,
   which most likely may contain singularities.

Functions to track solution paths
defined by complex artificial parameter homotopies 
are exported by the module ``trackers`` 
while the module ``sweepers``
exports path trackers for real natural parameter homotopies.
The module ``tuning`` helps to manage the tolerances of
the function to track the solution paths
in the ``trackers`` module.
Access to polyhedral end games is provided by the ``tropisms`` module.

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

Numerical representations of positive dimensional solution sets
are called witness sets and are
computed by the functions exported by the module sets.
Cascades of homotopies compute generic points on each component
of all dimensions.
In a numerical irreducible decomposition, all equidimensional solution
sets are factored into irreducible components.

functions in the module sets
----------------------------

.. automodule:: sets
   :members:

functions in the module cascades
--------------------------------

.. automodule:: cascades
   :members:

functions in the module factor
------------------------------

.. automodule:: factor
   :members:

functions in the module diagonal
--------------------------------

.. automodule:: diagonal
   :members:

some interesting families and examples
======================================

One of the motivations for phcpy was to perform regression tests
on the blackbox solver.

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

The module schubert exports Pieri homotopies
and Littlewood-Richardson homotopies to solve Schubert problems.

functions in the module schubert
--------------------------------

.. automodule:: schubert
   :members:

Newton polytopes, monomial maps, and power series
=================================================

The Newton polytope of a polynomial is spanned by the exponents 
of monomials which occur with nonzero coefficient in the polynomial.

functions in the module polytopes
---------------------------------

Given a polynomial, its support is the set of exponents of monomials
which occur with nonzero coefficient.  The convex hull of the support
of a polynomial is the Newton polytope of the polynomial.
For a polynomial system, the mixed volume of the Newton polytopes of
the polynomials in the systems gives a generically sharp upper bound
on the number of isolated solutions (not in coordinate planes)
of the polynomial system.

.. automodule:: polytopes
   :members:

functions in the module maps
----------------------------

A binomial system is a system where every equation has exactly two
monomials with nonzero coefficient. 
The solution set of a binomial system is a set of monomial maps.

.. automodule:: maps
   :members:

functions in the module series
------------------------------

Newton's method over the field of truncated power series
computes series expansions for solution curves.

.. automodule:: series
   :members:

functions in the module curves
------------------------------

Power series are input to Padé approximants,
which provide accurate predictors to approximate algebraic curves.

.. automodule:: curves
   :members:

a graphical user interface
==========================

With Tkinter we can develop a graphical user interface.
The module exports some possible development for a GUI
to solve polynomial systems and to process solutions.

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
