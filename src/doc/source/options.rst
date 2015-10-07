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

To reproduce a run with any seed (when the option -0 was not used),
we can look at the output file, for the line

::

   Seed used in random number generators : 407.

which appears at the end of the output of phc -b.
Running phc -b -0407 on the same input file as before
will generate the same sequences of random numbers
and thus the same output.

A homotopy is a family of polynomial systems.
In an artificial parameter homotopy there is one parameter \ :math:`t` 
usually going from 0 to 1.  
If we want to solve the system \ :math:`f({\bf x}) = {\bf 0}`
and we have a system \ :math:`g({\bf x}) = {\bf 0}`
then a typical homotopy is \ :math:`h({\bf x},t) = {\bf 0}` as below.

.. math::

   h({\bf x},t) = \gamma (1-t) g({\bf x}) + t f({\bf x}) = {\bf 0}.

The \ :math:`\gamma` is a complex constant on the unit circle
in the complex plane, generated uniformly at random.
If all solutions of \ :math:`g({\bf x}) = {\bf 0}` are isolated
and of multiplicity one, then only for a finite number of complex values 
for \ :math:`t` the homotopy \ :math:`h({\bf x},t) = {\bf 0}` has
singular solutions.  
But since we consider \ :math:`t \in [0,1]` and since the values
for \ :math:`t` for which \ :math:`h({\bf x},t) = {\bf 0}` are complex,
the interval \ :math:`[0,1(` will with probability one not contain
any of the bad complex values for \ :math:`t` and therefore no
singular solutions for \ :math:`t \in [0,1(` will occur.

Note that, for this gamma trick to work, the order of the operations matters.
We first give the program the system \ :math:`f({\bf x}) = {\bf 0}`
and then either also give the system \ :math:`g({\bf x}) = {\bf 0}`
or let the program generate a suitable \ :math:`g({\bf x}) = {\bf 0}`.
Only then, in the construction of the homotopy will a random number
generator determine a constant \ :math:`\gamma`.
If the \ :math:`\gamma` is predetermined, then it is possible to
construct an input system \ :math:`f({\bf x}) = {\bf 0}` and
a start system \ :math:`g({\bf x}) = {\bf 0}` for which there
are bad values for \ :math:`t \in [0,1(`.
But reverting the usual order of operations is a bit similar to guessing
the outcome of a coin toss after the coin toss and not before the coin toss.
Therefore phc -0 should be used only for debugging purposes.

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

New since version 2.4.02 are the options -b2 and -b4 to run the
blackbox solver respectively in double double and quad double precision,
for example as

::

   phc -b2 cyclic7 /tmp/c7out2
   phc -b4 cyclic7 /tmp/c7out4

The most computational intensive stage in the solver is in the
path tracking.  Shared memory multitasked path trackers are
available in the path trackers for both the polyhedral homotopies to solve
a random coefficient system and for the
artificial-parameter homotopy towards the target system.
See the documentation for the option phc -t below.

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

Numerical Schubert calculus is the development of numerical
homotopy algorithms to solve Schubert problems.  
A classical problem in Schubert calculus is the problem of four lines:
Given four lines in three dimensional space, find all lines that meet
the four given lines in a point.  If the lines are in general position,
then there are exactly two lines that satisfy the problem specification.
Numerical homotopy continuation methods deform a given generic problem
into special position, solve the problem in special position, and then
deform the solutions from the special into the generic problem.

As Schubert calculus is a classical topic in algebraic geometry,
what seems less well known is that Schubert calculus offers a solution
to the output pole placement problem in linear systems control.
The option phc -k offers one particular interface dedicated to the
Pieri homotopies to solve the output pole placement problem.
A related problem that can be solved with Schubert calculus is the 
completion of a matrix so that the completed matrix has a prescribed 
set of eigenvalues.

There are three types of homotopies in numerical Schubert calculus:

1. SAGBI homotopies solve hypersurface intersection conditions
   the extrinsic way.

2. Pieri homotopies are intrinsically geometric and are better able
   to solve more general problems in enumerate geometry.

3. Littlewood-Richardson homotopies resolve general Schubert
   intersection conditions.

The earliest instances of SAGBI and Pieri homotopies were already
available in version 2.0 of PHCpack.  
Since version 2.3.95, a more complete implementation of the 
Littlewood-Richardson homotopies is available.

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

The mixed volume of a tuple of Newton polytopes
if defined as the coefficient in the expansion of the volume
of a linear combination of Newton polytopes.
For example, for a 3-tuple of Newton polytopes:

.. math::

    \begin{array}{rcl}
      vol(\lambda_1 P_1 + \lambda_2 P_2 + \lambda_3 P_3)  
      & = & V(P_1, P_1, P_1) \lambda_1^3 \\
      & + & V(P_1, P_1, P_2) \lambda_1^2 \lambda_2 \\
      & + & V(P_1, P_2, P_2) \lambda_1 \lambda_2^2 \\
      & + & V(P_1, P_2, P_3) \lambda_1 \lambda_2 \lambda_3 \\
      & + & V(P_2, P_2, P_2) \lambda_2^3 \\
      & + & V(P_2, P_2, P_3) \lambda_2^2 \lambda_3 \\
      & + & V(P_2, P_3, P_3) \lambda_2 \lambda_3^2 \\
      & + & V(P_3, P_3, P_3) \lambda_3^3
    \end{array}

where \ :math:`vol(\cdot)` is the volume function
and \ :math:`V(\cdot)` is the mixed volume.
For the tuple \ :math:`(P_1, P_2, P_3)`, its mixed volume
is \ :math:`V(P_1,P_2,P_3)` in the expansion above.

The polynomial above can be called the *Minkowski polynomial*
and with the Cayley trick we can compute all its coefficients.
This is implemented with the dynamic lifting algorithm.

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

phc -o : writes the contents of the symbol table of an input system
-------------------------------------------------------------------

Running phc -o with as input argument a polynomial system
writes the symbols for the variables in the order in which they
are stored internally after parsing the system.
For example, if the file /tmp/ex1 contains the lines

::

   2
     y + x + 1;
     x*y - 1;

then running phc -o at the command prompt as

::

   phc -o /tmp/ex1 /tmp/ex1.out

makes the file /tmp/ex1.out which contains the line

::

   y x

because in the formulation of the polynomial system,
the variable with name y occurred before the variable with name x.
Consequently, the order of the coordinates of the solutions will
then also be stored in the same order as of the occurrence of the
variable names.  
If a particular order of variables would be inconvenient,
then a trick to force an order on the variables is to insert
a simple polynomial that simplifies to zero.  For example,
a modification of the file /tmp/ex1 could be

::

   2
    x + y - x - y +
    y + x + 1;
    x*y - 1;

and the first four monomials x + y - x - y will initialize the
symbol table with the names x and y, in that order.

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

Chapter 5 of the book of Alexander Morgan on
*Solving Polynomial Systems Using Continuation for Engineering
and Scientific Problems* (volume 57 in the SIAM Classics in
Applied Mathematics, 2009)
describes the setup of an optimization problem to compute coordinate 
transformations that lead to better values of the coefficients.

If the file ``/tmp/example`` contains the following lines

::

   2
    0.000001*x^2 + 0.000004*y^2 - 4;
    0.000002*y^2 - 0.001*x;

then a session with ``phc -s`` (at the command prompt) 
to scale the system goes as follows.

::

   $ phc -s
   Welcome to PHC (Polynomial Homotopy Continuation) v2.3.99 31 Jul 2015
   Equation/variable Scaling on polynomial system and solution list.

   MENU for the precision of the scalers :
     0. standard double precision;
     1. double double precision;
     2. quad double precision.
   Type 0, 1, or 2 to select the precision : 0

   Is the system on a file ? (y/n/i=info) y 

   Reading the name of the input file.
   Give a string of characters : /tmp/example

   Reading the name of the output file.
   Give a string of characters : /tmp/example.out

   MENU for Scaling Polynomial Systems :
     1 : Equation Scaling : divide by average coefficient      
     2 : Variable Scaling : change of variables, as z = (2^c)*x
     3 : Solution Scaling : back to original coordinates       
   Type 1, 2, or 3 to select scaling, or i for info : 2
     Reducing the variability of coefficients ? (y/n) y
     The inverse condition is  4.029E-02.

   Do you want the scaled system on separate file ? (y/n) y
   Reading the name of the output file.
   Give a string of characters : /tmp/scaled

   $ 

Then the contents of the file ``/tmp/scaled`` is

::

   2
   x^2+ 9.99999999999998E-01*y^2-1.00000000000000E+00;
   y^2-1.00000000000000E+00*x;

   SCALING COEFFICIENTS :

   10
   3.30102999566398E+00   0.00000000000000E+00
   3.00000000000000E+00   0.00000000000000E+00
   -6.02059991327962E-01   0.00000000000000E+00
   -3.01029995663981E-01   0.00000000000000E+00

We see that the coefficients of the scaled system are much nicer
than the coefficients of the original problem.
The scaling coefficients are needed to transform the solutions
of the scaled system into the coordinates of the original problem.
To transform the solutions, choose the third option of the second
menu of ``phc -s``.

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

With the time command, we can compare the wall clock time between
a sequential run and a run with 16 tasks:

::

   time phc -b cyclic7 /tmp/cyc7t1

   real    0m10.256s
   user    0m10.202s
   sys     0m0.009s

   time phc -b -t16 cyclic7 /tmp/cyc7t16

   real    0m0.851s
   user    0m11.149s
   sys     0m0.009s

The speedup on the wall clock time is about 12,
obtained as 10.256/0.851.

The relationship with double double and quad double precision
is interesting, consider the following sequence of runs:

::

   time phc -b cyclic7 /tmp/c7out1

   real    0m9.337s
   user    0m9.292s
   sys     0m0.014s

   time phc -b -t16 cyclic7 /tmp/c7out2

   real    0m0.923s
   user    0m13.034s
   sys     0m0.010s

With 16 tasks we get about a tenfold speedup,
but what if we ask to double the precision?

::

   time phc -b2 -t16 cyclic7 /tmp/c7out3

   real    0m4.107s
   user    0m59.164s
   sys     0m0.018s

We see that with 16 tasks in double precision, the elapsed time
equals 4.107 seconds, whereas the time without tasking was 9.337 seconds.
This means that with 16 tasks, for this example, we can double the
working precision and still finish the computation is less than half
of the time without tasking.

For quad double precision, more than 16 tasks are needed to offset
the overhead caused by the quad double arithmetic:

::

   time phc -b4 -t16 cyclic7 /tmp/c7out4

   real    0m53.865s
   user    11m56.630s
   sys     0m0.248s


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

4. Based on condition number estimates the working precision is set
   to meet the wanted number of accurate decimal places in the solutions
   when applying Newton's method.

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
