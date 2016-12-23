.. PHCpack documentation master file, created by
   sphinx-quickstart on Sun Jan 27 13:05:16 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

****************
Reference Manual
****************

The code is written in the following languages:
Ada, C, C++ (NVIDIA CUDA), and Python.
The following description documents the organization and
design decisions which led to the current state of the code.

The main executable ``phc`` compiles on Linux, MacOS X,
and Windows computers.  Shared memory parallelism works
on all three operating systems.
The message passing with MPI has not been tested on Windows.
The development of the accelerator code with NVIDIA CUDA 
was done on Linux computers.

The Python code was developed and tested on Linux and MacOS X,
not on Windows.

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

Every directory contains a collection of test procedures.

Organization of the C and C++ code
==================================

C code can be called from within Ada, as is the case
with the realization of the feedback laws in the output
placement problem, as defined in the ``Feedback`` directory.
A C (or C++) function may call Ada code, as was done in
the message passing code in the ``MPI`` directory.

Via the options of the main executable ``phc`` the user
navigates through menus and the data is stored on files.
The C interface defines a state machine with persistent objects.
As an example for the state machine metaphor,
consider a vending machine for snacks.  The user deposits coins,
makes a selection, and then retrieves the snacks.
The solution of a polynomial system via the C library happens
in the same manner.  The user enters the polynomials, either
from file or via their string representations, 
selects some algorithms, and then retrieves the solutions,
either from file, or in strings.

The Main Gateway Function
-------------------------

The directory ``Lib`` defines the C interface libraries.
In analogy with the single main executable ``phc``,
there is only one interface function which serves at the main gateway 
exporting the Ada functionality to the C and C++ programmers.

The header files in the definitions of the prototypes of the
library functions typically start with the following declarations:

::

   #ifdef compilewgpp
   extern "C" void adainit( void );
   extern "C" int _ada_use_c2phc ( int task, int *a, int *b, double *c );
   extern "C" void adafinal( void );
   #else
   extern void adainit( void );
   extern int _ada_use_c2phc ( int task, int *a, int *b, double *c );
   extern void adafinal( void );
   #endif

The ``adainit`` and ``adafinal`` are defined by the gnu-ada compiler.
They are required when the main program is not written in Ada.
Before the first call of the Ada code, ``adainit`` must be executed
and ``adafinal`` is required after the last call, before termination
of the program.

Persistent Objects
------------------

The C (or C++) can pass data via files or strings.
The definition of the data structures for the polynomials
and solution lists should not be duplicated in C (or C++).
Unless an explicit deallocation job is performed,
the objects remain in memory after a call to the Ada code.

The blackbox solver is exported by the C program ``phc_solve``.
The version which prompts the user for input and output files
starts as follows:

::

   int input_output_on_files ( int precision )
   {
      int fail,rc,nbtasks;

      if(precision == 0)
      {
         fail = syscon_read_standard_system();
         printf("\nThe system in the container : \n");
         fail = syscon_write_standard_system();
         printf("\nGive the number of tasks : "); scanf("%d",&nbtasks);
         fail = solve_system(&rc,nbtasks);
         printf("\nThe root count : %d\n",rc);
         printf("\nThe solutions :\n");
         fail = solcon_write_standard_solutions();
      }

The ``precision`` equal to zero is the default
standard double precision.  Other precisions that are supported
are double double and quad double precision.
If the number of tasks in ``nbtasks`` is a positive integer,
then the shared multicore version of the path trackers is executed.
The code below illustrates the use of persistent objects:
after the call to ``solve_system``, the solutions remain in main
memory even though only the value of the root count is returned
in ``rc``.  The solutions are printed with the call to
``solcon_write_standard_solutions()``.

Message Passing
===============

The shared memory parallelism is based on the tasking mechanism
defined by the Ada language and implemented by the gnu-ada compiler.
This section describes the distributed memory parallelism with
message passing, using the MPI library.  

The tracking of all solution paths is a pleasingly parallel computation
as the paths can be tracked independently from each other.
Some paths are more difficult to track than others and may require
more time, so dynamic load balancing in a manager/worker paradigm
often gives close to optimal speedups.
The setup suggested by :numref:`figprograminversion1`
is one wherein the manager solves the start system and
then distributes the start solutions to the worker nodes.

.. _figprograminversion1:

.. figure:: ./figprograminversion1.png
    :align: center

    A homotopy solver first solves the start system
    and then tracks all paths from start to target.

The setup in :numref:`figprograminversion1` leads to a top down control
in which the manager dictates the actions of the workers.
A more flexible setup is suggested in :numref:`figprograminversion2`:
start solutions are computed or retrieved when needed by the workers.

.. _figprograminversion2:

.. figure:: ./figprograminversion2.png
    :align: center

    The path tracker in a homotopy solver 
    calls for the next solution of the start system.

The advantage of the inverted control in
:numref:`figprograminversion2` over the more conventional setup in
:numref:`figprograminversion1` is the immediate availability of
solutions of the target system.
Moreover, the inverted control in :numref:`figprograminversion2`
does not require to store all start solutions.
For large polynomial systems, the number of start solutions may be 
too large to store in the main memory of one node.

GPU Acceleration
================

The acceleration with graphics processing units is code with
the NVIDIA compiler.

The Python Package phcpy
========================

The package phcpy provides a scripting interface.
For its functionality phcpy depends mainly on the C interface
and that was done on purpose: as the Python package grows,
so does the C interface.

The scripting interface to PHCpack has its own documentation.
