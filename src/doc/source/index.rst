.. PHCpack documentation master file, created by
   sphinx-quickstart on Sun Jan 27 13:05:16 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PHCpack's documentation!
===================================

PHCpack implements a collection of algorithms
to solve polynomial systems by homotopy continuation methods.

On input is a sequence of polynomials in several variables,
on output are the solutions to the polynomial system given on input.
The computational complexity of this problem is #P-hard because of
the exponential growth of the number of solutions as the number of
input polynomials increases.  For example, ten polynomials of degree two
may intersect in 1,024 isolated points (that is two to the power ten).
Twenty quadratic polynomials may lead to 1,048,576 solutions
(that is 1,024 times 1,024).  So it is not too difficult to write
down small input sequences that lead to a huge output.

Even as the computation of the total number of solutions may take
a long time, numerical homotopy continuation methods have the advantage
that they compute one solution after the other.  A homotopy is a family
of polynomial systems, connecting the system we want to solve with an
easier to solve system, which is called the start system.
Numerical continuation methods track the solution paths starting at
known solutions of an easier system to the system we want to solve.
We say that a homotopy is optimal if every path leads to a solution.

PHCpack offers optimal homotopies for systems that resemble 
linear-product structures, for geometric problems in enumerative geometry,
and for sparse polynomial systems with sufficiently generic choices 
of the coefficients.
While mathematically this sounds all good, most systems arising in
practical applications have their own peculiar structure
and so most homotopies will lead to diverging solution paths.
In general, a polynomial system may have solution sets of many
different dimensions, which renders the solving process challenging
but a the same time still very interesting.

Version 1.0 of PHCpack was archived as Algorithm 795
by ACM Transactions on Mathematical Software.  
PHCpack is open source and free software 
which gives any user the same rights as in free speech.

A short tutorial
================

This section provides a quick getting started guide.

Downloading and installing
--------------------------
Executable versions of the program for various machine architectures
and operating systems are available via
<http://www.math.uic.edu/~jan/download.html>.

If you are frugal with disk space,
it may be a good idea to download on the /tmp.
On Linux, type in the following commands to unzip
and untar the downloaded file:

::

   gunzip lin_phcv2p.tar.gz
   tar xpf lin_phcv2p.tar

If all went well, typing /tmp/phc at the command prompt should bring up 
the welcome message and the screen with available options.

Input formats
-------------

A lot of examples are contained in the database of Demo systems,
which can be downloaded in zipped and tarred format from the above
web site.

The input file starts with the number of equations and (optionally,
but necessary in case of an unequal number) the number of unknowns.
For example, the polynomial system of Bertrand Haas (which provided
a counterexample for the conjecture of Koushnirenko) is represented
as follows

::

   2
     x**108 + 1.1*y**54 - 1.1*y;
     y**108 + 1.1*x**54 - 1.1*x;

For future use, we save this system in the file haas.
Observe that every polynomial terminates with a semicolon.
The exponentiation may also be denoted by a hat instead of
a double asterix. 

Symbols that are forbidden to denote names of variables are
i and I, because they both represent the square root of -1.
Also forbidden are e and E because they are used in
the scientific notation of floating-point numbers, like
0.01 = 1.0e-2 = 1.0E-2.

The equations defining the adjacent 2-by-2 minors of
a general 2-by-4 matrix are represented as

::

   3 8
    x11*x22 - x21*x12;
    x12*x23 - x22*x13;
    x13*x24 - x23*x14;

thus as 3 polynomials in the 8 undeterminates of a general
2-by-4 matrix.  We save this file as adjmin4.

The program also accepts polynomials in factored form,
for example,

::

   5
    (a-1)*(b-5)*(c-9)*(d-13) - 21;
    (a-2)*(b-6)*(c-10)*(f-17) - 22;
    (a-3)*(b-7)*(d-14)*(f-18) - 23;
    (a-4)*(c-11)*(d-15)*(f-19) - 24;
    (b-8)*(c-12)*(d-16)*(f-20) - 25;

is a valid input file for phc.
Note that we replaced the logical e variable by f.
We save this input in the file with name multilin.

A very simple Maple interface
-----------------------------

The software is developed for command line interactions.
Because there is no interpreter provided with PHCpack,
there are interfaces to computer algebra systems like Maple.

From the web site mentioned above we can download the Maple procedure
run_phc and an example worksheet on how to use this procedure.
The Maple procedure requires only two arguments: the path name ending
in the name of the executable version of the program, and a list of
polynomials.  This procedure sets up the input file for phc, calls
the blackbox solver and returns the list of approximate solutions.
This list is returned in Maple format.

Other interfaces are PHClab (for Octave and MATLAB),
phc.py (for Sage), and PHCpack.m2 (for Macaulay 2).

Calling the blackbox solver
---------------------------

The blackbox solver works reasonably well to approximate all isolated
solutions of a polynomial system.  On the system we saved earlier in
the file multilin, we invoke the blackbox solver typing
at the command prompt

::

    /tmp/phc -b multilin multilin.phc

The output of the solver will be sent to the file multilin.phc.
In case the input file did not yet contain any solutions, 
the solution list will be appended to the input file.

We now explain the format of the solutions, for example, the last
solution in the list occurs in the following format:

::

   solution 44 :    start residual :  1.887E-14   #iterations : 1   success
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    a :  5.50304308029581E+00  -6.13068078142107E-44
    b :  8.32523889626848E+00  -5.18918337570284E-45
    c :  1.01021324864917E+01  -1.29182202179944E-45
    d :  1.42724963260133E+01   1.38159270467025E-44
    f :  4.34451307203401E+01  -6.26380413553193E-43
   == err :  3.829E-12 = rco :  3.749E-03 = res :  2.730E-14 = real regular ==

This is the actual output of the root refiner.  As the residual
at the end of the solution path and at the start of the root refinement
is already 1.887E-14, one iteration of
Newton's method suffices to confirm the quality of the root.

The next line in the output indicates that we reached the end of
the path, at t=1, properly.  The multiplicity of the root is one,
as indicated by m = 1.  Then we see the values for the five variables,
as pairs of two floating-point numbers: the real and imaginary part of
each value.  The last line summarizes the numerical quality of the root.
The value for err is the magnitude of the last correction term
used in Newton's method.  The number for rco is an estimate for
the inverse condition number of the root.  Here this means that we are
guaranteed to have all decimal places correct, except for the last three
decimal places.  The last number represents the residual, the magnitude
of the vector evaluated at the root.

Running the program in full mode
--------------------------------

If we just type in /tmp/phc without any option, we run the program
in full mode and will pass through all the main menus.
A nice application is the verification of the counterexample of Bertrand
Haas.  We type in haas when the program asks us for the name of
the input file.  As the output may be rather large, we better save the
output file on /tmp.  As we run through all the menus, for this system,
a good choice is given by the default, so we can type in 0 to answer
every question.  At the very end, for the output format, it may be good
to type in 1 instead of 0, so we can see the progress of the program as
it adds solution after solution to the output file.

If we look at the output file for the system in multilin,
then we see that the mixed volume equals the 4-homogeneous Bezout
number.  Since polyhedral methods (e.g. to compute the mixed volume)
are computationally more expensive than the solvers based on product
homotopies, we can solve the same problem faster.
If we run the program on the system in multilin in full mode,
we can construct a multi-homogeneous homotopy as follows.
At the menu for Root Counts and Method to Construct Start Systems,
we type in 1 to select a multi-homogeneous Bezout number.
Since there are only 52 possible partitions of a set of four unknowns,
it does not take that long for the program to try all 52 partitions
and to retain that partition that yields the lowest Bezout number.
Once we have this partition, we leave the root counting menu with 0,
and construct a linear-product system typing 2 in the menu to construct
m-homogeneous start systems.  We can save the start system in the file
multilin\_start (only used for backup).
Now we continue just as before.

Running the program in toolbox mode
-----------------------------------

To avoid the preconditioning stage (scaling and reduction) we can
compute root counts and construct start systems via the option -r,
thus calling the program as phc -r.  One important submenu is
the mixed-volume computation, invoked via phc -m.

Once we created an appropriate start system, we can call the path
trackers via the option -p.  Calling the program as phc -p
is useful if we have to solve a slightly modified problem.  For instance,
suppose we change the coefficients of the system in multilin,
then we can still use multilin_start to solve the system with
modified coefficients, using the -p option.  In this way we use
a cheater's homotopy, performing a kind of coefficient-parameter
polynomial continuation.

Dealing with components of solutions
------------------------------------

Consider the system of adjacent minors, we previously saved 
as adjmin4.  We first must construct a suitable embedding
to get to a system with as many equations as unknowns.
We call phc -c and type 5 as top dimension.  The system
the program produces is saved as adjmin4e5.  The blackbox
solver has no difficulty to solve this problem and appends the
witness points to the file adjmin4e5.  To compute the
irreducible decomposition, we may use the monodromy breakup
algorithm, selecting 2 from the menu that comes up when we
can the program with the option -f.

Options of the executable phc
=============================

For many small to moderate size problems,
the most convenient way to get an answer of phc is to
call the blackbox solver with the option -b.

phc -0 : random numbers with fixed seed for repeatable runs    
-----------------------------------------------------------

Many homotopy algorithms generate random constants.
With each run, the current time is used to generate
another seed for the random number generator, leading
to different random constants in each run.
As different random values give different random start systems,
this may cause differences in the solution paths and fluctuations
in the executation time.  Another notable effect of generating a
different random constant each time is that the order of the
solutions in the list may differ.  Although the same solutions
should be found with each run, a solution that appears first
in one run may turn out last in another run.

With the option -0, a fixed seed is used in each run.
This option can be combined with the blackbox solver (phc -b),
e.g.: phc -b -0 or phc -0 -b.

Since version 2.3.89, the option -0 is extended so the user may
give the digits of the seed to be used.  For example, 
calling phc as phc -089348224 will initialize the random
number generator with the seed 89348224.
Just calling phc as phc -0 will still result in using the same
fixed seed as before in each run.

phc -a : Solving polynomial systems equation-by-equation       
--------------------------------------------------------

The equation-by-equation solver applies the diagonal homotopies
to intersect the solution set of the current set of polynomials
with the next polynomial equation.  The user has the possibility
to shuffle the order of input polynomials.

phc -b : Batch or black-box processing                         
--------------------------------------

The blackbox solver operates in four stages:

1. Preprocessing: scaling (phc -s), 
   handle special cases such as binomial systems.

2. Counting the roots and constructing a start system (phc -r).

3. Track the solution paths from the solutions of the start system
   to the solutions of the target system (phc -p).

4. Verify whether all end points of the solution paths are distinct,
   apply Newton's method with deflation on singular solutions (phc -v).

Through the options -s, -r, -p, and -v, 
the user can go through the stages separately.

phc -c : Irreducible decomposition for solution components     
----------------------------------------------------------

The menu structure for a numerical irreducible decomposition 
consists of three parts:

1. Embedding and running cascade to compute witness points on
   all positive dimensional components of a polynomial system.
   See the option -f to factor a witness set for an pure
   dimensional solution set into irreducible components.

2. Intersecting witness sets with extrinsic and intrinsic
   diagonal homotopies.
   See the option -w for a wrapper to the diagonal homotopies.

3. Compute monomial maps for the irreducible decomposition
   of binomial systems.

phc -d : Linear and nonlinear Reduction w.r.t. the total degree
---------------------------------------------------------------

Degree bounds for the number of isolated solution often overshoot
the actual number of solution because of relationships between the
coefficients.  Consider for example the intersection of two circles.
A simple linear reduction of the coefficient matrix gives 
an equivalent polynomial system (having the same number of affine
solutions) but with lower degrees.  Reducing polynomials to
introduce more sparsity may also benefit polyhedral methods.

phc -e : SAGBI/Pieri homotopies to intersect linear subspaces  
-------------------------------------------------------------

There are three types of homotopies in numerical Schubert calculus:

1. SAGBI homotopies solve hypersurface intersection conditions
   the extrinsic way.

2. Pieri homotopies are intrinsically geometric and are better able
   to solve more general problems in enumerate geometry.

3. Littlewood-Richardson homotopies resolve general Schubert
   intersection conditions.

phc -f : Factor pure dimensional solution set into irreducibles
---------------------------------------------------------------

The -f is the f of factor and filter.

The first basic filter allows for example to extract the real
solutions from a given list.  The second filter implements the
homotopy membership test to decide whether a point belongs to
a witness set.  Given on input a witness set and a point,
this filter runs a homotopy to decide if the point belongs
to the positive dimensional solution set represented by
the given witness set.

The factorization method take on input a witness set for
a pure dimensional solution set of a polynomial system.
For small degrees, a combinatorial factorization method
will be fast.  The second factorization methods applies
monodromy loops, using the linear trace test as a stop criterion.

The last option in the menu of phc -f gives access to a
tropical method to detect a common factor of two Laurent polynomials.

phc -k : realization of dynamic output feedback placing poles  
-------------------------------------------------------------

The homotopies in numerical Schubert calculus (see the option -e)
solve the output pole placement problem in linear systems control.
The option -k applies the Pieri homotopies to compute feedback laws
for plants defined by (A,B,C) matrices.

phc -l : Witness Set for Hypersurface cutting with Random Line 
--------------------------------------------------------------

A hypersurface defined by a polynomial in several variables is
cut with one general line.  The number of points on the hypersurface
and the general line equal the degree of the hypersurface.
This collection of points on the intersection of a hypersurface
and a general line form a witness set.

For example, if the file ``sphere`` contains

::

   1 3
   x^2 + y^2 + z^2 - 1;

then typing at the command prompt

::

   phc -l sphere sphere.out

results in the creation of the file ``sphere_w2`` which contains
a witness set of dimension two for the unit sphere.
The output file ``sphere.out`` contains diagnostics about the computation.

phc -m : Mixed-Volume Computation via lift+prune and MixedVol  
-------------------------------------------------------------

The options of phc -m are a subset of the options of phc -r.
The focus on phc -m is on mixed volumes.  For polynomial systems
with as many equations as unknowns, the mixed volume of the Newton
polytopes gives a generically sharp upper bound on the number of
isolated solutions with coordinates different from zero.

The menu of phc -m offers 5 different algorithms:

0. Static lifting: a lifting function is applied to the points in the
   support sets of the polynomials in the system and the lower hull
   defines the mixed cells.
   The users can specify the lifting values interactively.
   Liftings that do not lead to cells that are fine mixed
   are subdivided with a random lifting.

1. Implicit lifting: based on a recursive formula used in Bernshtein's
   original proof that the mixed volumes bounds the number of isolated
   solutions with nonzero coordinates.

2. Dynamic lifting: points are added one after the other in an
   incremental construction of a mixed cell configuration.
   An implementation of the Cayley trick gives the Minkowski polynomial.

3. Symmetric lifting: many systems have Newton polytopes that are
   invariant to permutation symmetry.  Even if the original system
   is not symmetric, the construction of the start system could
   benefit from the exploitation of this permutation symmetry.

4. The MixedVol Algorithm is a specific implementation of the
   static lifting method, applying a floating random lifting function.

   The code offered with this option is a translation of software
   described in the paper by Tangan Gao, T. Y. Li, Mengnien Wu:
   *Algorithm 846: MixedVol: a software package for mixed-volume 
   computation.*
   ACM Transactions on Mathematical Software, 31(4):555-560, 2005;
   distributed under the terms of the GNU General Public License as
   published by the Free Software Foundation.

   Stable mixed volumes count all affine solutions
   (not only those with nonzero coordinates) 
   and lead to polyhedral homotopies that compute all affine solutions.

phc -p : Polynomial Continuation by a homotopy in one parameter
---------------------------------------------------------------

The user of phc -p is prompted for a target system and a start system
with start solutions.  This option is useful for solving several systems
with the same structure but with different coefficients.

With phc -p, the user has full access to all numerical tolerances
that define how close the numerical approximations have to stay
along a solution path.  Another application of phc -p is to rerun
a selection of solution paths.

In addition to the artificial parameter increment-and-fix continuation,
there is support for complex parameter continuation
and real pseudo arc length path tracking with detection of singularities
using the determinant of the Jacobian along the solution path.

phc -q : Tracking Solution Paths with incremental read/write   
------------------------------------------------------------

The jumpstarting method for a polynomial homotopy
does not require the computation of all solutions of the start system
and neither does it keep the complete solution list in memory.

phc -r : Root counting and Construction of start systems       
--------------------------------------------------------

Methods to bound the number of isolated solutions of a polynomial system
fall in two classes:

1. Bounds based on the highest degrees of polynomials and variable groupings.

2. Bounds based on the Newton polytopes of the polynomials in the system.
   See the documentation for phc -m.

phc -s : Equation and variable Scaling on system and solutions 
--------------------------------------------------------------

We distinguish two types of scaling:

1. Equation scaling: multiplying every coefficient in the same equation
   by the same constant.

2. Variable scaling: multiplying variables with constants.

phc -t : Tasking for tracking paths using multiple threads     
----------------------------------------------------------

The option -t allows the user to take advantage
of multicore processors.
For example, typing at the command prompt.

::

   phc -b -t4 cyclic7 /tmp/cyclic7.out

makes that the blackbox solver uses 4 threads to solve the system.
If there are at least 4 computational cores available,
then the solver may finish its computations up to 4 times faster
than a sequential run.

phc -v : Verification, refinement and purification of solutions
---------------------------------------------------------------

While solution paths do in general not become singular or diverge,
at the end of the paths, solutions may turn out to be singular
and/or at infinity.  With phc -v one can do the following tasks:

1. Perform a basic verification of the solutions based on Newton's method
   and weed out spurious solutions.

2. Apply Newton's method with multiprecision arithmetic.

3. For isolated singular solutions, the deflation method may recondition
   the solutions and restore quadratic convergence.

phc -w : Witness Set Intersection using Diagonal Homotopies    
-----------------------------------------------------------

This option wraps the diagonal homotopies to intersect two witness sets,
see the option -c for more choices in the algorithms.

For example, to intersect the unit sphere 
(see the making of ``sphere_w2`` with phc -l) with a cylinder
to form a quartic curve, we first make a witness set for a cylinder,
putting in the file ``cylinder`` the two lines:

::

   1 3
   x^2 + y - y + (z - 0.5)^2 - 1; 

Please note the introduction of the symbol ``y``
even though the symbol does not appear in the equation of a cylinder
about the y-axis.  But to intersect this cylinder with the unit sphere
the symbols of both witness sets must match.
After executing ``phc -l cylinder cylinder.out`` we get the witness
set ``cylinder_w2`` and then we intersect with phc -w:

::

   phc -w sphere_w2 cylinder_w2 quartic

The file ``quartic`` contains diagnostics of the computation.
Four general points on the quartic solution curve of the intersection
of the sphere and the cylinder are in the file ``quartic_w1``
which represents a witness set.

phc -x : convert solutions from PHCpack into Python dictionary 
--------------------------------------------------------------

To work with solution lists in Python scripts, the script phc -x
convert a solution list in PHCpack format to a list of dictionaries.
Given a Python list of dictionaries, phc -x returns a list of
solutions in PHCpack format.  For example:

::

   phc -x cyclic5 /tmp/cyclic5.dic
   phc -x /tmp/cyclic5.dic

The first phc -x writes to the file /tmp/cyclic5.dic a list of
dictionaries, ready for processing by a Python script.
If no output file is given as second argument, then the output
is written to screen.  The second phc -x writes a solution list
to PHCpack format, because a list of dictionaries is given on input.

phc -y : sample points from an algebraic set, given witness set
---------------------------------------------------------------

Given in ``sphere_w2`` a witness set for the unit sphere
(made with phc -l, see above), we can make a new witness set
with phc -y, typing at the command prompt:

::

   phc -y sphere_w2 new_sphere

and answering two questions with parameter settings
(type 0 for the defaults).  The output file ``new_sphere``
contains diagnostics of the run and a new witness set is
in the file ``new_sphere_w2``.

phc -z : strip phc output solution lists into Maple format     
----------------------------------------------------------

The phc -z commands converts solution lists in PHCpack format
into Maple lists and converts Maple lists into solutions lists 
in PHCpack format.  For example:

::

   phc -z cyclic5 /tmp/cyclic5.mpl
   phc -z /tmp/cyclic5.mpl

If the file ``cyclic5`` contains the solutions of the cyclic 5-roots
problem in PHCpack format, then the first command makes the file 
/tmp/cyclic5.mpl which can be parsed by Maple.  The next command
has no second argument for output file and the output is written
directly to screen, converting the solutions in Maple format to
solution lists is PHCpack format.

An Application Programming Interface to PHCpack
===============================================

Because code development on PHCpack has taken a very long time,
looking at the code may be a bit too overwhelming at first.
A good starting point could be the Python interface
and in particular phcpy, with documentation at
<http://www.math.uic.edu/~jan/phcpy_doc_html/index.html>.

References
==========

1. T. Gao, T. Y. Li, M. Wu:
   **Algorithm 846: MixedVol: a software package for mixed-volume 
   computation.**
   *ACM Transactions on Mathematical Software*, 31(4):555-560, 2005.

2. E. Gross, S. Petrovic, and J. Verschelde: **PHCpack in Macaulay2.**
   *The Journal of Software for Algebra and Geometry: Macaulay2*,
   5:20-25, 2013.

3. Y. Guan and J. Verschelde: 
   **PHClab: A MATLAB/Octave interface to PHCpack.**
   In *IMA Volume 148: Software for Algebraic Geometry*,
   edited by M. E. Stillman, N. Takayama, and J. Verschelde,
   pages 15-32, Springer-Verlag, 2008. 

4. Y. Hida, X.S. Li, and D.H. Bailey:
   **Algorithms for quad-double precision floating point arithmetic.**
   In *15th IEEE Symposium on Computer Arithmetic (Arith-15 2001)*,
   11-17 June 2001, Vail, CO, USA, pages 155-162.
   IEEE Computer Society, 2001.
   Shortened version of Technical Report LBNL-46996.

5. A. Leykin and J. Verschelde: 
   **PHCmaple: A Maple Interface to the Numerical Homotopy Algorithms
   in PHCpack.**
   In the *Proceedings of the Tenth International Conference 
   on Applications of Computer Algebra (ACA'2004)*,
   edited by Q. N. Tran, pages 139-147, 2004.

6. A. Leykin and J. Verschelde: 
   **Interfacing with the Numerical Homotopy Algorithms in PHCpack.**
   In *proceedings of ICMS 2006, LNCS 4151*,
   edited by A. Iglesias and N. Takayama,
   pages 354-360, Springer-Verlag, 2006. 

7. A. J. Sommese, J. Verschelde, and C. W. Wampler.
   **Numerical irreducible decomposition using PHCpack.**
   In *Algebra, Geometry, and Software Systems*, 
   edited by M. Joswig and N. Takayama,
   pages 109-130. Springer-Verlag, 2003.

8. J. Verschelde:
   **Algorithm 795: PHCpack: A general-purpose solver for polynomial
   systems by homotopy continuation.**
   *ACM Transactions on Mathematical Software*, 25(2):251--276, 1999.

9. J. Verschelde:
   **Polynomial homotopy continuation with PHCpack.**
   *ACM Communications in Computer Algebra*, 44(4):217-220, 2010.

10. J. Verschelde:
    **Modernizing PHCpack through phcpy.**
    In Proceedings of the 6th European Conference on Python in Science
    (EuroSciPy 2013), edited by Pierre de Buyl and Nelle Varoquaux,
    pages 71-76, 2014, available at
    <http://arxiv.org/abs/1310.0056>.

Acknowledgments
===============

This material is based upon work supported by the 
National Science Foundation under Grants No. 9804846, 0105739, 0134611,
0410036, 0713018, and 1115777.
Any opinions, findings, and conclusions or recommendations expressed 
in this material are those of the author(s) and do not necessarily 
reflect the views of the National Science Foundation. 
