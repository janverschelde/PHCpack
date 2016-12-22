.. PHCpack documentation master file, created by
   sphinx-quickstart on Sun Jan 27 13:05:16 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

****************
Reference Manual
****************

The code is written in the following languages:
Ada, C, C++ (NVIDIA CUDA), and Python.

Organization of the Ada code
============================

The code in the first release was written in Ada 83.
Release 2 featured a new mathematical library,
rebuilt using Ada 95 concepts and offering multi-precision arithmetic.

There are four major layers in the code:

1. ``Math_Lib``: linear algebra, representations of polynomials,
   Newton polytopes, and power series 

2. ``Deformations``: Newton's method, path trackers, end games, 
   solutions and homotopies, deflation

3. ``Root_Counts``: root counting methods and constructions of homotopies,
   linear-product start start systems based on Bezout bounds,
   mixed volumes and polyhedral homotopies

4. ``Components``: witness sets, cascades of homotopies, monodromy, 
   diagonal homotopies, to compute a numerical irreducible decomposition

There are five other parts, called ``System``, ``Schubert``, ``CtoPHC``,
``PHCtoC``, and ``Tasking``.  The top down perspective starts at the
folder ``Main``.

The Ada sources are organized in a tree of directories:

::

   Ada                       : Ada source code of PHC
    |-- System               : 0. OS dependencies, e.g.: timing package
    |-- Math_Lib             : 1. general mathematical library
    |      |-- Numbers       : 1.1. number representations
    |      |-- QD            : 1.2. quad doubles
    |      |-- Vectors       : 1.3. vectors and vectors of vectors
    |      |-- Matrices      : 1.4. matrices and linear-system solvers
    |      |-- Divisors      : 1.5. common divisors, integer linear algebra
    |      |-- Reduction     : 1.6. row reduction, numerical linear algebra
    |      |-- Polynomials   : 1.7. multivariate polynomial systems
    |      |-- Functions     : 1.8. evaluation and differentiation
    |      |-- Supports      : 1.9. support sets and linear programming
    |      |-- Circuits      : 1.A. circuits for algorithmic differentation
    |      |-- Series        : 1.B. manipulating truncated series
    |-- Deformations         : 2. homotopies, Newton's method & path trackers
    |      |-- Solutions     : 2.1. solutions of systems and homotopies
    |      |-- Homotopy      : 2.2. homotopies, scaling and reduction
    |      |-- Newton        : 2.3. root refining and modified Newton's method
    |      |-- Curves        : 2.4. univariate solving & plane algebraic curves
    |      |-- End_Games     : 2.5. extrapolation end games with Puiseux series
    |      |-- Trackers      : 2.6. path-tracking routines
    |      |-- Sweep         : 2.7. sweeping for singularities
    |      |-- Continuation  : 2.8. drivers and data management
    |-- Root_Counts          : 3. root counts and homotopy construction
    |      |-- Product       : 3.1. linear-product start systems
    |      |-- Binomials     : 3.2. solvers for binomial and simplicial systems
    |      |-- Implift       : 3.3. implicit lifting
    |      |-- Stalift       : 3.4. static lifting
    |      |-- Dynlift       : 3.5. dynamic lifting
    |      |-- Symmetry      : 3.6. exploitation of symmetry relations
    |      |-- MixedVol      : 3.7. translation of ACM TOMS Algorithm 846
    |      |-- Puiseux       : 3.8. Puiseux series for curves
    |-- Schubert             : 4. numerical Schubert calculus
    |      |-- SAGBI         : 4.1. SAGBI homotopies
    |      |-- Pieri         : 4.2. deformations based on Pieri's rule
    |      |-- Induction     : 4.3. Schubert induction
    |-- Components           : 5. numerical irreducible decomposition
    |      |-- Samplers      : 5.1. computing witness points
    |      |-- Interpolators : 5.2. finding equations for components
    |      |-- Factorization : 5.3. factorization into irreducible components
    |      |-- Decomposition : 5.4. sequence of homotopies to filter and factor
    |      |-- Solver        : 5.5. incremental equation by equation solver
    |      |-- Tropical      : 5.6. tropical view on witness sets
    |-- CtoPHC               : 6. interface from C to phc
    |      |-- Funky         : 6.1. functional interface, C -> Ada -> C 
    |      |-- State         : 6.2. state machine gateway, C <-> Ada
    |-- PHCtoC               : 7. GPU acceleration via a C interface
    |-- Tasking              : 8. multitasking
    |-- Main                 : 9. main dispatcher
