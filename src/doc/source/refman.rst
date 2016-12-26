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
The following sections describe the functionality defined
in each of the directories.

System: OS Dependencies such as Timing
--------------------------------------

The ``System`` directory defines operations that may have different
definitions on different operation systems.  One such operation is
to compute the elapsed CPU time of a computation.
The timer for Ada on Unix like operation systems was originally
developed by Dave Emory of the MITRE corporation.
Not everything in this timing package could be mapped to Windows,
in particular the resource usage report for Unix.
While the interface of the timing package is the same for all operating
systems, the implementation differs for Windows

When multithreaded runs on multicore processors, the elapsed CPU time
is most often not a good time measurement and one comes interested in
the wall clock time.  The end of the output contains the start and end
date of the computation.  With the ``Ada.Calendar``, the time stamping
is defined in a portable, operating system independent manner.

The directory system contains several very useful utilities,
such as procedures to prompt the user for a yes or no answer,
or for a selection between various alternatives.
While restricting the user selection, the prompting procedures
allow to retry in case of type errors.
Similar user friendly guards are defined when the user gives
the name of an existing file for output.  Before overwriting
the existing file, the user is prompted to confirm.
When reading a file, the user is allowed to retry in case the
given name of the file does not match an existing file.

The handling of the command line options is also defined in this
directory.  Thanks to the ``Ada.Command_Line``, this definition
is operating system independent.

The package ``machines`` wraps some system calls.
One such system call is to get the process identification number (pid).
This pid is used to seed the random number generators.

The Mathematical Library
------------------------

The mathematical library defines code that is not specific
to polynomial homotopy continuation, but nevertheless necessary.
To make PHCpack self contained, the code does not require the
installation of outside libraries.  Although there are eleven
subdirectories, there are three main parts:

1. number representations, general multiprecision and quad doubles;

2. linear algebra with integers and floating-point numbers;

3. polynomials, polynomial functions, series, and Newton polytopes.

The input to a polynomial system solver is a list of polynomials in
several variables.  This input consists of exact data, such as the
integer exponents in the monomials, and approximate data, such as
the floating-point coefficients of the monomials.
Solving a polynomial system with homotopy continuation is therefore
always a hybrid computation, involving exact and approximate data.
While the machine arithmetic may still suffice for many applications,
the increasing available computational power has led to the formulation
of large problems for which software defined multiprecision arithmetic
is required.  The linear algebra operations are defined over exact
number rings and over arbitrary precision floating-point numbers.

The next subsections contain more detailed descriptions of each
subdirectory of the mathematical library.
The following three paragraphs briefly summarize the eleven 
subdirectories in the three main parts.

The number representations are defined in the subdirectory ``Numbers``
and the QD library of Y. Hida, X. S. Li, and D. H. Bailey is integrated
in the subdirectory ``QD``.

The linear algebra data structures are defined in the subdirectories
``Vectors`` and ``Matrices``.  The ``Divisors`` subdirectory relies
on the greatest common divisor algorithm to define the Hermite and
Smith normal forms to solve linear systems over the integer numbers.
The linear system solvers of numerical linear algebra are provided
in the subdirectory ``Reduction``.

The third main part of the mathematical library consists in the
remaining five of the eleven subdirectories.  Multivariate polynomials
over various number rings in the subdirectory ``Polynomials``.
The subdirectory ``Functions`` contains definitions of 
nested Horner schemes to efficiently evaluate dense polynomials.
The support of a polynomial is the set of exponents of the monomials
which appear with nonzero coefficients.  Basic linear programming
and tools to work with polytopes are provided in the subdirectory
``Supports``.  The subdirectory ``Circuits`` defines arithmetic
circuits to evaluate and differentiate polynomials via the reverse
mode of algorithmic differentiation.  Truncated power series define
a field (that is: dividing two series gives again a series)
and the arithmetic to manipulate power series is exported by the
packages in the subdirectory ``Series``.

Deforming Polynomial Systems
----------------------------

A homotopy is a family of polynomial systems defined by one parameter.
The parameter may be introduced in an artificial manner, such as
the parameter :math:`t` in the classical homotopy

.. math::

   h({\bf x}, t) = (1 - t) g({\bf x}) + t f({\bf x}) = {\bf 0}.

The homotopy :math:`h({\bf x}, t)` connects the system
:math:`g({\bf x}) = {\bf 0}` (the so-called *start system*) to the system
:math:`f({\bf x}) = {\bf 0}` (the so-called *target system*),
as :math:`h({\bf x}, 0) = g({\bf x})`
and :math:`h({\bf x}, 1) = f({\bf x})`.
The solutions :math:`{\bf x}(t)` to the homotopy are solution paths,
starting at :math:`t=0` at the solutions of the start system
and ended at :math:`t=1` at the solutions of the target system.

The code was developed mainly for constructing artificial-parameter
homotopies, but there is some still limited support for polynomial
homotopies with natural parameters.  Artificial-parameter homotopies
can be constructed so that singular solutions occur only at the end
of the paths.  For natural-parameter homotopies, the detection and
accurate computation of singularities along the paths becomes an
important topic.

There are eight subdirectories in the ``Deformations`` directory.

Homotopy Construction via Root Counting Methods
-----------------------------------------------

At first, it seems counter intuitive to construct a polynomial homotopy
to solve an unknown system by counting its roots.
But consider the degeneration of two planar quadrics into lines.
Each quadric degenerates to a pair of lines.  How many solutions
could we get intersection two pairs of lines in general position?
Indeed, four, computed as two by two.  Observe that in this simple
argument we have no information about the particular representation
of the quadrics.  To get to this root count, we assumed only that
the lines after degeneration were generic enough and the count
involved only the degrees of the polynomials.

Of critical importance for the performance of a polynomial homotopy
is the accuracy of the root count.  If the root count is a too large
upper bound for the number of solutions of the system that will be
solved, then too many solution paths will diverge to infinity,
representing a very wasteful computation.

We can construct homotopies based on the degree information alone
or rely on the Newton polytopes.
Sparse polynomial systems are systems where relatively few monomials
appear with nonzero coefficient, relative to the degrees of the
polynomials in the system.  
For sparse system, the information of the Newton polytopes provides
a much sharper root count than the ones provided by the degrees.

The are eight subdirecties in the ``Root_Counts`` directory.

Numerical Schubert Calculus
---------------------------

The classical problem in Schubert calculus asks for the number
of lines which meet four given general lines in 3-space.
With polynomial homotopies, we not only count, but also compute
the actual number of solutions to a Schubert problem.

The problem of four lines is a special case of a Pieri problem:
compute all *p*-planes which meet :math:`m \times p` given *m*-planes 
in a space of dimension :math:`m + p`.  If the given *m*-planes are 
sufficiently generic, then all solution *p*-planes are isolated and
finite in number.  Pieri homotopies solve the output pole placement
problem in linear systems control.

There are three subdirectories to the ``Schubert`` directory,
each exporting a different type of homotopy to solve Schubert problems.

The subdirectory ``SAGBI`` applies the concept of
subalgebra analog to Groebner basis for ideals
with polyhedral homotopies to solve Pieri problems.

Pieri homotopies are defined in the subdirectory ``Pieri``.

The subdirectory ``Induction`` implements a geometric
Littlewood-Richardson rule to solve general Schubert problems.

Numerical Irreducible Decomposition
-----------------------------------

Two important characteristics of a pure dimensional solution set of 
a polynomial system are its dimension and its degree.
The dimension of a solution set equals the number of general linear equations
we need to add to the polynomial system so the intersection of the solution
set of the system with the hyperplanes consists of isolated points.
The degree of a solution set then equals the number of isolated points
we find after intersecting the solution set with as many general hyperplanes
as the dimension of the set.
These two characteristics are encoded in the *witness set*
representation of a pure dimensional solution set.
Given a polynomial system, a numerical irreducible decomposition
of its solution set provides a witness set for each irreducible
components, over all dimensions.

The decomposition can be computed in a top down fashion,
with cascades of homotopies, starting a the top dimension.
The bottom up computation applies diagonal homotopies.
Systems can be solved equation-by-equation or subsystem-by-subsystem.

Three types of factorization methods are implemented.
Interpolation with multivariate polynomials of increasing degrees 
is a local procedure.  The second method runs monodromy loops to
connect generic points on the same irreducible component,
using the linear trace test as stop criterion.  
Thirdly, we can apply the linear trace test combinatorially,
which often works very well for components of modest degrees.

The are six subdirectories of the ``Components`` directory.
The ``Samplers`` subdirectory contains the definitions of the data
structures to store witness sets.  The multivariate interpolation
algorithms are implemented in the ``Interpolators`` subdirectory.
The subdirectory ``Factorization`` provides monodromy factorization
and the linear trace test.  Cascades of homotopies and diagonal
homotopies are implemented in the subdirectory ``Decomposition``.
The ``Solver`` subdirectory provides an equation-by-equation solver.
Finally, the ``Tropical`` subdirectory offers code to generalize 
the polyhedral homotopies from isolated solutions to the computation
of representations of positive dimensional solution sets.

Calling Ada Code From C
-----------------------

The directory ``CtoPHC`` has two subdirectories, ``Funky`` and ``State``,
which define two different types of interfacing the Ada code with C.
The first type is a functional interface, the second type is an interface
which operates as a state machine.

In a functional interface, the main C program calls an Ada function,
which then calls a C function to process the results computed by the
Ada function.  This interface was developed for the application of
the Pieri homotopies to compute output feedback laws for linear systems
control.  This type of interface is direct and efficient.
Its main application is in the ``Feedback`` folder which defines C
functions to compute realizations of the computed feedback laws.

The goal of the state interfce in the subdirectory ``State`` is to
export all functionality of the Ada code to the C (and C++) programmer.
The subdirectory ``State`` contains the definition of the
``use_c2phc`` function, which defines more than 700 jobs.
The implementation of this function relies on various container
packages which hold the persistent objects, mainly polynomial systems
and solution lists.

Calling C Code From Ada
-----------------------

The directory ``PHCtoC`` was set up to call the GPU code via a C interface.
In its current state it defines the wrappers to call the accelerated
path trackers with algorithmic differentiation.
Its main goal is to define the extension modules for calling the
accelerated path trackers from the Python package phcpy.

Multitasking
------------

The Ada tasking mechanisms allows to define shared memory parallel
programs at a high level.  Tasks in Ada are mapped to kernel threads.
There are two main applications defined in the ``Tasking`` directory.

Given a queue of path tracking jobs, the tasks are arranged in
a work crew model to execute all jobs.  Dynamic load balancing
is achieved as tasks, when done with their current job, grab the
next job from the queue.  Synchronization overhead is minimal,
as only the movement of the current pointer in the job queue
happens in a critical section.
This parallel work crew path tracking scheme is implemented for
regular homotopies and polyhedral homotopies.

Another application of multitasking is pipelining.
Polyhedral homotopies start at initial form systems computed by
the mixed cells.  For large polynomial systems, the computation
of the mixed volume could be a bottleneck for the parallel execution.
A pipelined multitasked implementation of the polyhedral homotopies
combines the tracking of all paths with the mixed cell computation
as follows.  One task computes the mixed cells and appends the
mixed cells to the job queue.  Other tasks take the mixed cells
as the jobs to solve the random coefficient system.
As soon as one mixed cells is available in the queue,
the path tracking can start.

The Main Program
----------------

The directory ``Main`` contains the main program,
called ``dispatch`` because its main function is to dispatch
the options given at the command line to the specific procedures.

The code for the blackbox solver (invoked by ``phc -b``)
is defined by the packages ``black_box_solvers``
and ``black_box_root_counters``.

A very specific solver is defined by the file ``use_phc.adb``,
mainly as an example how the code could be customized for one
particular application.  The code is below:

::

   with text_io;                            use text_io;
   with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
   with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
   with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
   with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
   with PHCpack;

   procedure use_phc is

     infile,outfile : file_type;        -- input and output file
     p,q : Link_to_Poly_Sys;            -- target and start system
     mixed_volume : natural32;          -- root count is mixed volume
     sols : Solution_List;              -- list of solutions
   
   begin
     Open(infile,in_file,"test.in");
     get(infile,p);
     Create(outfile,out_file,"test.out");
     put(outfile,p.all);
     q := new Poly_Sys(p'range);
     PHCpack.Static_Lifting(outfile,p.all,mixed_volume,q.all,sols);
     PHCpack.Artificial_Parameter_Continuation(outfile,p.all,q.all,sols);
     PHCpack.Refine_Roots(outfile,p.all,sols);
   end use_phc;

Numbers, Linear Algebra, Polynomials and Polytopes
==================================================

In this section we take a closer look at the ``Math_Lib`` directory,
which defines the basic mathematical data structures and operations.

Numbers
-------

The machine numbers are divided in two categories: integer and float.
For the integer types, we distinguish between the 32-bit and 64-bit
versions, between natural and integer numbers.  The following types are
defined: ``natural32``, ``natural64``, ``integer32``, and ``integer64``.
For the float types, we have single precision and double precision,
defined respectively as ``single_float`` and ``double_float``.
The renaming of the hardware number types ensures the independence
of pre-defined number types.

For polynomial system solving, our default field is the field of
complex numbers.  The real and imaginary part of a complex number
are floating-point coefficients.  The homotopy algorithms depend
on the choice of random constants.  Random number generators are
defined.  The default seed for the random number generators is the
process identification number.  For reproducible runs, the user can
set the seed to a fixed number.

Multiprecision numbers are implemented as arrays of machine integers.
Elementary school algorithms defined the arithmetic.
The implementation of the floating-point multiprecision numbers
is directly based on the multiprecision integer numbers,
for the fraction and the exponent part of the multiprecision float.
The precision of each multiprecision number can be adjusted when needed,
which is an advantage.  Mixed-precision arithmetical operations are
supported.  The disadvantage imposed by this flexibility is the
frequent memory allocation and deallocation, which makes this type of
arbitrary multiprecision arithmetic unsuitable for shared memory parallelism.

The directory ``Numbers`` contains definitions of abstract rings, domains,
and fields.  These abstract classes are useful to define composite
generic types.  Multiprecision complex numbers are defined via the
instantiation of a generic complex numbers package.

Quad Doubles
------------

The directory ``QD`` provides the double double and quad double arithmetic,
based on the QDlib package of Y. Hida, X. S. Li, and D. H. Bailey.

Compared to arbitrary multiprecision arithmetic, double double and quad
double numbers exploit the floating-point hardware and have a simple
memory management.  While arbitrary multiprecision numbers are allocated
via the heap, the two doubles of a double double and the four doubles
of a quad double use the stack.  Thus the QD library is very well suited
for shared memory parallelism.  Another advantage is the predictable
cost overhead.  Working with double doubles has a similar cost overhead
as working with complex numbers.  Computations with double doubles are about
five to eight times slower compared to computations in double precision.
With quad doubles, computations that took seconds in double precision
can turn into minutes.

The code in QDlib was hand translated into Ada.
The directory contains the original C versions for comparison
and verification of correctness.

Vectors and Matrices
--------------------

The directories ``Vectors`` and ``Matrices`` contain the definitions
of respectively all vector and matrix types.
In both directories, generic packages are defined, which allow to
specify the ring of numbers (natural32, integer32, natural64, integer64)
or the number fields (double, double double, quad double, or arbitrary
multiprecision).  Input and output for all types is provided.

Although both ``Vectors`` and ``Matrices`` are basic data structures,
random number generators are provided, to generate vectors and matrices
of random numbers.  The test procedures check the basic arithmetical
operations.

The directory ``Vectors`` defines vectors of vectors and 
vectors of matrices are defined in the directory ``Matrices``.

Linear Systems with Integer Coefficients
----------------------------------------

The problem considered in the directory ``Divisors``
is the manipulation of matrices with integer coefficients.

With the greatest common divisor we can define unimodular coordinate
transformations to compute an upper triangular form of a matrix with
integer coefficients.  Such form is call the Hermite normal form.
The diagonalization process results in the Smith normal form.

Even if the input matrices have small integer coefficients,
the size of the integers in the unimodular coordinate transformations
can outgrow the size of the hardware integers.
Therefore, multiprecision versions of the normal forms are provided.

This integer linear algebra is applied in the computation of the
volumes of the mixed cells of subdivisions of Newton polytopes.

Linear Systems with Floating-Point Coefficients
-----------------------------------------------

The directory ``Reduction`` contains several matrix factorizations
as common in numerical linear algebra.

The LU factorization is based on the ``lufac``, ``lufco``,
and ``lusolve`` of the F77 LINPACK libary.
The Fortran77 code was translated into Ada and extended with versions 
for double double, quad double, and arbitrary multiprecision;
both for real and complex number types.

To solve overdetermined linear systems in the least squares sense,
packages are provided for the QR decomposition.  
Also the Singular Value Decomposition (SVD) is implemented,
for all precisions, and for real and complex number types.

To implement a variable precision Newton's method, there are
variable precision linear system solvers.
Given the desired accuracy,
the variable precision linear system solver sets the working
precision based on a condition number estimate.

Polynomials in Several Variables
--------------------------------

Multivariable polynomials and polynomial systems are defined
in the directory ``Polynomials``.  In addition to ordinary polynomials,
polynomials with integer exponents, so-called Laurent polynomials,
are defined as well.  In solving Laurent polynomials, solutions
with zero coordinates are excluded.

Nested Horner Forms for Evaluation
----------------------------------

Because the evaluation and differentiation of polynomials can be
just as expensive as solving a linear system in the application of
Newton's method, the distributed list of terms in a polynomial is
converted into a nested Horner form, for efficient evaluation.
The directory ``Functions`` provides specific data structures
to construct and evaluate the nested Horner forms.

Support Sets and Linear Programming
-----------------------------------

Given a list of vectors with integer coefficients,
via linear programming we can extract from the list those points
which are vertex points of the polytope spanned by the points
in the list.  Another application of linear programming is
the computation of all k-dimensional faces of the polytope.
The directory ``Supports`` provides the primitive operations
for the volume computations in the polyhedral root counts.

Circuits for Algorithmic Differentiation
----------------------------------------

The directory ``Circuits`` contains implementations of the algorithms
which evaluate and differentiate polynomials in several variables using
the reverse mode of algorithmic differentiation.

Truncated Power Series
----------------------

Similar to Taylor series approximations for general functions,
we can approximate roots of polynomials in a parameter by series.
The directory ``Series`` defines truncated power series with
complex numbers as coefficients.  Composite types are vectors,
matrices, and polynomials where the coefficients are series.

The division of two truncated power series is computed via
the solution of a triangular linear system.
So we can have a field and we can solve linear systems over
this field of truncated power series.  However to work efficiently,
instead of working with vectors and matrices of power series,
we apply linearization and consider series where the coefficients
are vectors and matrices.

The directory exports packages to solve linear systems where
the coefficient matrix is a power series of matrix coefficients.
We can solve such linear systems with LU factorization, or
for overdetermined problems we solve in the least squares sense,
either with a QR or an SVD decomposition.
To solve Hermite-Laurent interpolation problems,
a lower triangular echelon form is provided.

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

The acceleration with Graphics Processing Units (GPUs) is coded with
the NVIDIA compiler.  GPUs are designed for data parallel applications.  
Their execution model is single instruction multiple data: 
the same instruction is executed on many different data elements.  
Unlike shared memory parallelism with threads on multicore processors, 
to fully occupy a GPU, one must launch ten thousands of threads.

Polynomial homotopy continuation methods can take advantage of GPUs
by the evaluation and differentiation of polynomials as required in
the frequent application of Newton's method.  The reverse mode of
algorithmic differentiation applied to the monomials with appear
with a nonzero coefficient in the polynomials provides sufficient
parallelism and a granularity fine enough for the data parallel
execution model.  The same arithmetic circuits to evaluate and
differentiate monomials are applied to different solutions when
tracking many solution paths.  For the tracking of one path in
large enough dimension, different threads collaborate in the
evaluation and differentiation algorithms.

To introduce the evaluation and differentiation algorithms
consider :numref:`figcirceval4` and :numref:`figcircdiff4`
to compute the product of four variables and its gradient.
Observe that results from the evaluation can be recycled in
the computation of all partial derivatives.

.. _figcirceval4:

.. figure:: ./figcirceval4.png
    :align: center

    An arithmetic circuit to evaluate the product of four variables
    :math:`x_1`, :math:`x_2`, :math:`x_3`, and :math:`x_4`.

.. _figcircdiff4:

.. figure:: ./figcircdiff4.png
    :align: center

    An arithmetic circuit to compute the gradient of
    the product :math:`x_1 x_2 x_3 x_4`.

The computation of the gradient of :math:`x_1 x_2 \cdots x_8` is
illustrated in :numref:`figcircdiff8`.

.. _figcircdiff8:

.. figure:: ./figcircdiff8.png
    :align: center

    An arithmetic circuit to compute the gradient of the product
    of eight variables
    :math:`x_1`, :math:`x_2`, :math:`\ldots`, and :math:`x_8`.

The Python Package phcpy
========================

The package phcpy provides a scripting interface.
For its functionality phcpy depends mainly on the C interface
and that was done on purpose: as the Python package grows,
so does the C interface.

There are several other scripting interfaces to PHCpack:
to the computer algebra system Maple (PHCmaple), 
PHClab for MATLAB and Octave, and for Macaulay2: PHCpack.m2.
These other interfaces rely only on the executable version of the program.

Another major difference between phcpy and other scripting
interface is the scope of exported functionality.
The main goal of phcpy is to export all functionality of ``phc``
to the Python programmer.  The development of phcpy can be viewed
as a modernization of the PHCpack code, bringing it into 
Python's growing computational ecosystem.

The scripting interface to PHCpack has its own documentation.
