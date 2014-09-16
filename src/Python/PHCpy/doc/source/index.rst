.. phcpy documentation master file, created by
   sphinx-quickstart on Thu Nov 29 13:01:46 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to phcpy's documentation!
=================================

This documentation describes a collection of Python modules
to compute solutions of polynomial systems using PHCpack.
While PHCpack has been under development for over twenty years,
phcpy is still working through its proof-of-concept stage.
In its present state, working with phcpy will require persistence
and plenty of patience.

The main executable phc (polynomial homotopy continuation)
defined by the source code in PHCpack is a menu driven
and file oriented program.
The Python interface defined by phcpy replaces the files
with persistent objects allowing the user to work with
scripts or in interactive sessions.
The computationally intensive tasks such as path tracking
and mixed volume computations are executed as compiled code
so there will not be a loss of efficiency.

Both phcpy and PHCpack are free and open source software packages;
you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software Foundation;
either version 2 of the License, or (at your option) any later version.  

One application of phcpy is to run regression tests.
The Python interface phcpy to PHCpack is a programmer's interface.
The long-term goal is to make PHCpack more versatile,
at least for those programmers familiar 
with the Python scripting language.

The source for PHCpack can be downloaded from
<http://www.math.uic.edu/~jan/download.html>
and is also available at
<https://github.com/janverschelde/PHCpack>.
For the installation from source, the gnu-ada compiler 
(available for free at <http://libre.adacore.com/>) is needed.  
Also your Python interpreter must most likely have been built with gcc.
In addition to the Python interpreter, the file Python.h of
of the developer's libraries will need to exist on your system.
Otherwise, PHCpack and phcpy are self contained
and do not require the installation of other software.

The computation of the mixed volume in phcpy calls MixedVol
(ACM TOMS Algorithm 846 of T. Gao, T.Y. Li, M. Wu) 
as it is integrated in PHCpack.
For double double and quad double arithmetic, PHCpack incorporates
the QD library of Y. Hida, X.S. Li, and D.H. Bailey.
See the References section for pointers to the literature.

For Red Hat Linux 64-bit and some Mac OS X versions,
the distribution provides the shared object file phcpy2c.so.
For phcpy to work, the file phcpy2c.so must be present at the
same location as the Python modules.

This document arose as an exercise in exploring restructured text
and Sphinx.  Even with a wonderful tool like Sphinx, this documentation
is (just as phcpy) very much a work in progress...

a blackbox solver for isolated solutions
----------------------------------------

The package phcpy depends on the shared object file phcpy2c.so.
The module **phcpy.solver**
exports the blackbox solver of PHCpack, a fast mixed volume
calculator, and several methods to predict the number of isolated
solutions of a polynomial system.  
The `test_solver()` function of the module generates two trinomials 
(a polynomial with three monomials)
with randomly generated complex coefficients.

By default, the input polynomial systems are expected to be *square*,
that is: the number of polynomials in the input list equals the number
of variables in the polynomials.  The blackbox solver then returns
a list of numerical approximations to the isolated solutions of the
input polynomial system.  Some capabilities of PHCpack to deal with
positive dimensional solution sets are exported by 
the module **phcpy.sets**.

Polynomials and solutions are represented as strings.
Below is an illustration of a session with the blackbox solver
on a system of two random trinomials, polynomials with three
monomials with random complex coefficients.

::

   >>> from phcpy.solver import random_trinomials
   >>> f = random_trinomials()

To see what the polynomials in f look like, 
let us print its terms (splitting on the `+(` 
to avoid long swirling lines in the output):

::
   
   >>> terms = f[0].split('+(')
   >>> for t in terms: print t
   ...
   -0.991457094247+0.13043324065066*i)*x^0*y^5
   -0.0509953121395-0.99869889262970*i)*x^5*y^3
   0.232308584664+0.97264213433887*i)*x^4*y^4;
   >>>

To solve the system defined by f, we call the blackbox solver:

::

   >>> from phcpy.solver import solve
   >>> s = solve(f,silent=True)
   >>> len(s)
   15
   >>> print s[2]
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x :  7.10290847804173E-01  -4.69841154290980E-01
    y : -1.79580717006598E-01  -7.16541556066137E-01
   == err :  1.986E-16 = rco :  2.676E-01 = res :  1.232E-16 =
   >>>

The *solve* command returned a list of 15 strings in s,
each string represents a solution that makes the polynomials in f vanish.
The module **phcpy.solutions** (documented below)
offers a function to evaluate the solutions
in the polynomials given as strings.

The solver in PHCpack generates different random numbers with each run,
which may very well cause the solutions to appear in a different order
after a second application of solve on the same system.
To prevent this behaviour (to check reproducibility for example),
we can fix the seed of the random number generators in PHCpack,
as follows:

::

   >>> from phcpy.phcpy2c import py2c_set_seed
   >>> py2c_set_seed(2013)
   0
   >>>

The performance of the solver is very sensitive to how accurately
we can predict the number of solutions.  For dense polynomial systems,
looking at the highest degrees of the polynomials in the system suffices,
whereas for sparse polynomial systems, computing the mixed volume of
the Newton polytopes of the polynomials yields much better results.
Below is a simple example, illustrating the bounds based on the
degrees and the mixed volume:

::

   >>> f = ['x^3*y^2 + x*y^2 + x^2;', 'x^5 + x^2*y^3 + y^2;']
   >>> from phcpy.solver import total_degree
   >>> total_degree(f)
   25
   >>> from phcpy.solver import m_homogeneous_bezout_number as mbz
   >>> mbz(f)
   (19, '{ x }{ y }')
   >>> from phcpy.solver import linear_product_root_count as lrc
   >>> lrc(f)
   a supporting set structure :
        { x }{ x }{ x }{ y }{ y }
        { x }{ x }{ x y }{ x y }{ x y }
   the root count : 19
   19
   >>> from phcpy.solver import mixed_volume
   >>> mixed_volume(f, stable=True)
   (14, 18)
   >>>

The mixed volume is a generically sharp root count for the number of 
isolated solutions with all coordinates different from zero. 
The term *generically sharp* means: except for systems with coefficients 
in a specific collection of algebraic sets, the root count is an exact count.
The stable mixed volume counts all affine solutions, 
that is: also the solutions with zero coordinates.
For the example above, we may expect at most 14 isolated solutions 
with all coordinates different from zero, 
and, also considering solutions with zero coordinates, 
at most 18 isolated solutions, counted with multiplicities.

For every root count, total degree, m-homogeneous Bezout number,
linear-product root count, and mixed volume, there is a corresponding
method to construct a polynomial system with exactly as many regular
solutions at the root count, which can then be used as a start system
in a homotopy to compute all isolated solutions of the polynomial system 
for which the root count was computed.
Examples of the methods to construct start systems in phcpy
are illustrated in the documentation for the module **phcpy.trackers**.

The other functions exported by **phcpy.solver** concern the movement
of data between Python and PHCpack.  The `store_` methods parse strings
representing polynomials and solutions into the corresponding internal
data structures.  The corresponding `load_` methods take the internal
data structures for polynomials and solutions, stored in containers,
and show their corresponding representations as Python strings.
For example, consider the session

::

   >>> from phcpy.solver import store_standard_system, load_standard_system
   >>> store_standard_system(['x^2 - 1/3;'])
   >>> load_standard_system()
   ['x^2 - 3.33333333333333E-01;']

The session above illustrates the parsing of a system one could use
to approximate the square root of 1/3.  With standard double precision,
the 1/3 is approximated to about 15 decimal places.

Newton's method fails when the Jacobian matrix is singular
(or close to singular) at a solution.  Below is a session
on the example of A. Griewank and M. R. Osborne, in their paper
*Analysis of Newton's method at irregular singularities,*
published in *SIAM J. Numer. Anal.* 20(4): 747-773, 1983.
The origin (0,0) is an irregular singularity: Newton's method
fails no matter how close the initial guess is taken.
With deflation we can restore the quadratic convergence
of Newton's method:

::

   >>> p = ['(29/16)*x^3 - 2*x*y;', 'x^2 - y;']
   >>> from phcpy.solutions import make_solution
   >>> s = make_solution(['x','y'],[1.0e-6,1.0e-6])
   >>> print s
   t : 0.0 0.0
   m : 1
   the solution for t :
    x : 1.000000000000000E-06  0.0
    y : 1.000000000000000E-06  0.0
   == err : 0.0 = rco : 1.0 = res : 0.0 ==
   >>> from phcpy.solver import newton_step
   >>> s2 = newton_step(p,[s])
   == err :  1.000E-06 = rco :  5.625E-13 = res :  1.875E-19 =
   >>> print s2[0]
   t :  0.00000000000000E+00   0.00000000000000E+00
   m : 0
   the solution for t :
    x :  9.99999906191101E-07   0.00000000000000E+00
    y :  9.99999812409806E-13   0.00000000000000E+00
   == err :  1.000E-06 = rco :  5.625E-13 = res :  1.875E-19 =
   >>> s3 = newton_step(p,s2)
   == err :  3.333E-07 = rco :  2.778E-14 = res :  1.111E-13 =
   >>> print s3[0]
   t :  0.00000000000000E+00   0.00000000000000E+00
   m : 0
   the solution for t :
    x :  6.66666604160106E-07   0.00000000000000E+00
    y :  3.33333270859482E-13   0.00000000000000E+00
   == err :  3.333E-07 = rco :  2.778E-14 = res :  1.111E-13 =
   >>> from phcpy.solver import deflate
   >>> sd = deflate(p,[s])
   >>> print sd[0]
   t :  0.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -4.55355758042535E-25   2.75154683741089E-26
    y :  1.57904709676279E-25  -8.86785799319512E-26
   == err :  5.192E-13 = rco :  5.314E-03 = res :  1.388E-16 =
   >>>

The documentation strings of the functions
exported by the module ``solver`` of the package phcpy are listed below.

.. automodule:: solver
   :members:

representations of solutions of polynomial systems 
--------------------------------------------------

Solutions of phcpy.solve are returned as lists of PHCpack
solution strings.  The solutions module contains functions to
parse a PHCpack solution string into a dictionary.

The solutions module exports operations 

1. to parse strings in the PHCpack solution format into dictionaries;

2. to evaluate these dictionaries into polynomials substituting the
   values for the variables into the strings representing the polynomials.

The main test in the module solutions is the solution of a small
trinomial system and the evaluation of the computed solutions
at the trinomial system.

The information of a solution as a dictionary contains the following:

1. `t` : value of the continuation parameter

   `m` : multiplicity of the solution

2. symbols for the variables are keys in the dictionary,
   the corresponding values are complex floating-point numbers

3. `err` : magnitude of the last correction term of Newton's method
   (forward error)

   `rco` : estimate for the inverse of the condition number of
   the Jacobian matrix at the solution

   `res` : magnitude of the residual (backward error)

The triplet (err, rco, res) measures the numerical quality of the solution.
The residual `res` is normally interpreted as an estimate of the backward
error: by how much should we change the original problem such that the
approximate solution becomes an exact solution of the changed problem.
The estimate `rco` gives a (sometimes too pessimistic) bound on the
number of correct decimal places in the approximate solution.
In particular: `abs(log(rco, 10))` bounds the number of lost decimal
places in the approximate solution.
For example, if `rco` equals `1.0E-8`, then the last 8 decimal places
in the coordinates of the solution could be wrong.

For a solution of the example ``noon3`` from the module examples,
we convert the PHCpack format solution string to a dictionary as follows:

::

   >>> print s[0]
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x1 : -1.65123467890611E-01  -7.61734168646636E-01
    x2 :  8.98653694263692E-01  -3.48820047576431E-01
    x3 :  8.98653694263692E-01  -3.48820047576431E-01
   == err :  3.034E-16 = rco :  2.761E-01 = res :  5.974E-16 =
   >>> from phcpy.solutions import strsol2dict
   >>> d = strsol2dict(s[0])
   >>> d.keys()
   ['err', 'res', 'm', 'rco', 't', 'x2', 'x3', 'x1']
   >>> d['x1']
   (-0.165123467890611-0.761734168646636j)
   >>> 

Note that the values of the dictionary d are evaluated strings,
parsed into Python objects.

By plain substitution of the values of the dictionary representation
of the solution into the string representation of the polynomial system
we can verify that the coordinates of the solution evaluate to numbers
close to the numerical working precision:

::

   >>> from phcpy.solutions import evaluate
   >>> e = evaluate(f,d)
   >>> for x in e: print x
   ... 
   (1.11022302463e-15+4.4408920985e-16j)
   (7.77156117238e-16+9.99200722163e-16j)
   (7.77156117238e-16+9.99200722163e-16j)
   >>> 

A more elaborate verification of the solution is provided by
the function **newton_step** of the module ``solver`` of phcpy.

The documentation strings of the functions
exported by the module ``solutions`` are listed below.
The script **test()** runs when typing **python solutions.py**
at the command prompt.

.. automodule:: solutions
   :members:

the path trackers
-----------------

Homotopy continuation methods are applied to solve a polynomial system.
The module **phcpy.trackers** exports the path trackers of PHCpack.

The example session below illustrates the computation of the intersection
of an ellipse with a parabola.  A homotopy method based on the total degree
replaces the two given quadratic equations for the ellipse and the parabola
by a configuration of lines that has exactly as many solutions as the
expected number of intersection points.  The homotopy connects the given
system with the equations of the simpler configuration, which define the
start system.  Continuation methods track the paths starting at the solutions
of the start system to the solutions of the target system.

::

   >>> from phcpy.solver import total_degree
   >>> from phcpy.solver import total_degree_start_system
   >>> from phcpy.trackers import track
   >>> p = ['x^2 + 4*y^2 - 4;','2*y^2 - x;']
   >>> d = total_degree(p)
   >>> d
   4
   >>> (q,qsols) = total_degree_start_system(p)
   >>> len(qsols)
   4
   >>> q
   ['x^2 - 1;', 'y^2 - 1;']
   >>> s = track(p,q,qsols)
   >>> len(s)
   4
   >>> for sol in s: print sol
   ... 
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x :  1.23606797749979E+00   0.00000000000000E+00
    y :  7.86151377757423E-01   0.00000000000000E+00
   == err :  1.309E-16 = rco :  1.998E-01 = res :  4.441E-16 =
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x :  1.23606797749979E+00   0.00000000000000E+00
    y : -7.86151377757423E-01   0.00000000000000E+00
   == err :  1.309E-16 = rco :  1.998E-01 = res :  4.441E-16 =
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -3.23606797749979E+00   0.00000000000000E+00
    y :  0.00000000000000E+00   1.27201964951407E+00
   == err :  1.505E-36 = rco :  1.079E-01 = res :  0.000E+00 =
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -3.23606797749979E+00   0.00000000000000E+00
    y :  0.00000000000000E+00  -1.27201964951407E+00
   == err :  1.505E-36 = rco :  1.079E-01 = res :  0.000E+00 =
   >>>

As expected when we intersect two quadratic equations,
we find four intersection points.  The coordinates of
the solutions are complex numbers, listed as two consecutive
floating-point numbers in scientific notation.
The two consecutive numbers approximate the real and imaginary part
of the complex number.  In the four solutions above, observe that
two solutions are real and two solutions are complex conjugate.

Note that the start system ``q`` in ``['x^2 - 1;', 'y^2 - 1;']``
has four real solutions, while the system ``p`` we solve had two
complex conjugate solutions.  If we connect ``p`` to ``q`` 
with a real homotopy, then at some point along the path, 
two real solutions have to turn into a pair of complex conjugate solutions.
Multiplying the start system with a random complex constant,
we avoid the singularities along the solution paths.
The side effect of this multiplication is that different constants
will results in different orders of the solutions at the end.
For example:

::

   >>> from phcpy.solver import total_degree_start_system
   >>> from phcpy.trackers import track
   >>> p = ['x^2 + 4*y^2 - 4;','2*y^2 - x;']
   >>> (q,qsols) = total_degree_start_system(p)
   >>> s1 = track(p,q,[qsols[2]])
   >>> print s1[0]
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x :  1.23606797749979E+00   0.00000000000000E+00
    y :  7.86151377757423E-01   0.00000000000000E+00
   == err :  1.383E-16 = rco :  1.998E-01 = res :  2.220E-16 =
   >>> s2 = track(p,q,[qsols[2]])
   >>> print s2[0]
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -3.23606797749979E+00   0.00000000000000E+00
    y :  0.00000000000000E+00   1.27201964951407E+00
   == err :  4.815E-35 = rco :  1.079E-01 = res :  0.000E+00 =
   >>>

To avoid this side effect, *track* accepts a complex value
as its last argument for the so-called gamma constant.
As a continuation of the session from above:

::

   >>> s3 = track(p,q,[qsols[2]],gamma=complex(0.824372806319,0.56604723848934))
   >>> print s3[0]
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -3.23606797749979E+00   0.00000000000000E+00
    y :  0.00000000000000E+00   1.27201964951407E+00
   == err :  0.000E+00 = rco :  1.079E-01 = res :  0.000E+00 =
   >>>

If we track all solution paths one after the other,
each time calling track with the same value for gamma,
then all solutions will be found.

The ``track`` function follows a solution path till the end.
Often it could be useful to view all intermediate solutions
computed along a path.  The functions
``next_standard_solution()``,
``next_dobldobl_solution()``,
``next_quaddobl_solution()``, and
``next_multprec_solution()``,
implement generators for the path trackers in standard double,
double double, quad double precision, and arbitrary multiprecision 
respectively.  
With these ``next_`` functions, the user not only gets all solutions
along a path, but also receives control of the order of execution.
Before the application of ``next_``, one must initialize the homotopy
with target and start system and give an initial start solution.
The session below illustrates the use of this generator:

::

   >>> from phcpy.solver import total_degree_start_system
   >>> p = ['x**2 + 4*x**2 - 4;', '2*y**2 - x;']
   >>> (q,s) = total_degree_start_system(p)
   >>> from phcpy.trackers import initialize_standard_tracker
   >>> from phcpy.trackers import initialize_standard_solution
   >>> from phcpy.trackers import next_standard_solution
   >>> initialize_standard_tracker(p,q)
   >>> initialize_standard_solution(len(p),s[0])
   >>> s1 = next_standard_solution()
   >>> print s1
   t :  1.00000000000000E-01   0.00000000000000E+00
   m : 1
   the solution for t :
    x :  9.96338438384030E-01   4.70831004481527E-03
    y :  9.96408320626402E-01   4.95310952563875E-03
   == err :  2.375E-05 = rco :  1.000E+00 = res :  3.619E-10 =
   >>> print next_standard_solution()
   't :  2.00000000000000E-01   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  9.80919860804043E-01   1.78496473654540E-02
     y :  9.81218221286503E-01   2.32056259678926E-02
   == err :  1.671E-08 = rco :  1.000E+00 = res :  1.424E-16 ='
   >>> print next_standard_solution()
   t :  3.00000000000000E-01   0.00000000000000E+00
   m : 1
   the solution for t :
    x :  9.51909891692765E-01   2.71534790078036E-02
    y :  9.42895891640611E-01   5.51080014180090E-02
   == err :  4.812E-09 = rco :  1.000E+00 = res :  1.665E-16 =
   >>>

In the session above, we see the solutions ``s1`` for t = 0.1,
and two other solutions for consecutive values 0.2 and 0.3 for t.
If we continue the session from above with the second start solution
in ``s[1]``, we can select the first 11 points along the path
and view all values for ``x`` of the solutions:

::

   >>> initialize_standard_solution(len(p),s[1])
   >>> points = [next_standard_solution() for i in range(11)]
   >>> from phcpy.solutions import strsol2dict
   >>> dicpts = [strsol2dict(sol) for sol in points]
   >>> xvals = [sol['x'] for sol in dicpts]
   >>> for x in xvals: print x
   ... 
   (0.996338438384+0.00470831004482j)
   (0.980919860804+0.0178496473655j)
   (0.951909891693+0.0271534790078j)
   (0.924234166108+0.0231054530961j)
   (0.908102639672+0.0141598112703j)
   (0.90039366434+0.00726313574566j)
   (0.896843555845+0.00320608226584j)
   (0.895239133202+0.00112430968375j)
   (0.894586634218+0.000224845127444j)
   (0.894427191-2.20881053462e-28j)
   (0.894427191+0j)
   >>>

We see that the last two values differ little from each other
because we arrived at the end of the path.  
To test whether at the end of a path, it suffices to check
whether the value for t equals one.

The image below plots the real parts of the four paths.
Three of the paths converge to the triple solution (1,2).

.. image:: ./plotpaths.png

The code used to make the plot (using matplotlib) is below:

::

   p = ['x^2 + y - 3;', 'x + 0.125*y^2 - 1.5;']
   print 'constructing a total degree start system ...'
   from phcpy.solver import total_degree_start_system as tds
   q, qsols = tds(p)
   print 'number of start solutions :', len(qsols)
   from phcpy.trackers import initialize_standard_tracker
   from phcpy.trackers import initialize_standard_solution
   from phcpy.trackers import next_standard_solution
   initialize_standard_tracker(p, q, False)
   from phcpy.solutions import strsol2dict
   import matplotlib.pyplot as plt
   plt.ion()
   fig = plt.figure()
   for k in range(len(qsols)):
       if(k == 0):
           axs = fig.add_subplot(221)
       elif(k == 1):
           axs = fig.add_subplot(222)
       elif(k == 2):
           axs = fig.add_subplot(223)
       elif(k == 3):
           axs = fig.add_subplot(224)
       startsol = qsols[k]
       initialize_standard_solution(len(p),startsol)
       dictsol = strsol2dict(startsol)
       xpoints =  [dictsol['x']]
       ypoints =  [dictsol['y']]
       for k in range(300):
           ns = next_standard_solution()
           dictsol = strsol2dict(ns)
           xpoints.append(dictsol['x'])
           ypoints.append(dictsol['y'])
           tval = eval(dictsol['t'].lstrip().split(' ')[0])
           if(tval == 1.0):
               break
       print ns
       xre = [point.real for point in xpoints]
       yre = [point.real for point in ypoints]
       axs.set_xlim(min(xre)-0.3, max(xre)+0.3)
       axs.set_ylim(min(yre)-0.3, max(yre)+0.3)
       dots, = axs.plot(xre,yre,'r-')
       fig.canvas.draw()
   fig.canvas.draw()
   ans = raw_input('hit return to exit')

With *False* in
*initialize_standard_tracker(p, q, False)*
the option to generate a fixed gamma constant was turned off,
so rerunning the same code will generate other random constants
and produce different plots.

Below is an interactive session to illustrate the solving 
with polyhedral homotopies.

::

   >>> p = ['x^3*y^2 - 3*x^3 + 7;','x*y^3 + 6*y^3 - 9;']
   >>> from phcpy.solver import mixed_volume
   >>> mixed_volume(p)
   11
   >>> from phcpy.solver import random_coefficient_system
   >>> (q,qsols) = random_coefficient_system(silent=True)
   >>> len(qsols)
   11
   >>> from phcpy.trackers import track
   >>> psols = track(p,q,qsols)
   >>> len(psols)
   11
   >>> print psols[4]
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -7.33932797408386E-01  -9.84310202527377E-01
    y : -6.56632351304388E-01   9.90969278772793E-01
   == err :  1.938E-16 = rco :  5.402E-01 = res :  2.102E-15 =
   >>>

We can apply one Newton step with higher precision to improve
the accuracy of the solutions.  Doubling the precision:

::

   >>> psols_dd = newton_step(p,psols,precision='dd')
   == err :  4.380E-15 = rco :  4.624E-01 = res :  4.239E-28 =
   == err :  5.190E-15 = rco :  3.266E-03 = res :  3.342E-27 =
   == err :  5.051E-15 = rco :  1.978E-02 = res :  2.727E-28 =
   == err :  4.306E-15 = rco :  3.778E-01 = res :  3.768E-28 =
   == err :  5.015E-16 = rco :  5.402E-01 = res :  3.525E-30 =
   == err :  5.015E-16 = rco :  5.402E-01 = res :  3.525E-30 =
   == err :  4.187E-15 = rco :  6.236E-01 = res :  4.236E-28 =
   == err :  4.187E-15 = rco :  6.236E-01 = res :  4.236E-28 =
   == err :  4.611E-15 = rco :  5.158E-01 = res :  2.719E-28 =
   == err :  4.306E-15 = rco :  3.778E-01 = res :  3.768E-28 =
   == err :  4.611E-15 = rco :  5.158E-01 = res :  2.719E-28 =
   >>>

We see that the residual (the parameter *res*) drops for every solution.

Below is an illustration of the use of linear-product start systems:

::

   >>> p = ['x*y^3 + y - 2;', 'x^3*y + x - 8;']
   >>> from phcpy.solver import linear_product_root_count
   >>> r = linear_product_root_count(p)
   a supporting set structure :
        { x }{ y }{ y }{ y }
        { x }{ x }{ x }{ y }
   the root count : 10
   >>> from phcpy.solver import random_linear_product_system
   >>> (q,qsols) = random_linear_product_system(p)
   >>> len(qsols)
   10
   >>> from phcpy.trackers import track
   >>> psols = track(p,q,qsols)
   >>> len(psols)
   10
   >>> from phcpy.solver import newton_step
   >>> psols_dd = newton_step(p,psols,precision='dd')
   == err :  6.197E-15 = rco :  1.606E-01 = res :  2.268E-28 =
   == err :  6.197E-15 = rco :  1.606E-01 = res :  1.446E-28 =
   == err :  2.453E-15 = rco :  2.699E-01 = res :  7.116E-29 =
   == err :  5.269E-15 = rco :  2.918E-01 = res :  1.374E-28 =
   == err :  2.453E-15 = rco :  2.699E-01 = res :  7.116E-29 =
   == err :  4.108E-15 = rco :  2.707E-01 = res :  9.348E-29 =
   == err :  5.855E+30 = rco :  1.078E-92 = res :  7.123E+93 =
   == err :  2.332E-15 = rco :  2.877E-01 = res :  2.931E-29 =
   == err :  5.269E-15 = rco :  2.918E-01 = res :  1.374E-28 =
   == err :  6.753E+29 = rco :  5.037E-91 = res :  2.547E+90 =
   >>> 

Looking at the values for *err* and *res* we see huge values
for two solutions which are spurious.

Last but certainly not least, consider the application of multitasking
to path tracking.  On the benchmark problem of cyclic 7-roots:

::

   $ time python trackcyclic7.py
   number of start solutions : 924
   starting the path tracking with 1 task(s) ...
   tracked 924 solution paths

   real    0m7.147s
   user    0m7.126s
   sys     0m0.016s
   $ time python trackcyclic7.py 2
   number of start solutions : 924
   starting the path tracking with 2 task(s) ...
   tracked 924 solution paths

   real    0m3.927s
   user    0m7.640s
   sys     0m0.017s
   $ 

Observe that the wall clock time (the time following the *real*),
is cut almost in half when 2 tasks are used.
The script is below:

::

   from sys import argv
   if(len(argv) == 1):
       nbtasks = 1
   else:
       nbtasks = eval(argv[1])
   from phcpy.phcpy2c import py2c_read_standard_target_system_from_file
   from phcpy.phcpy2c import py2c_read_standard_start_system_from_file
   from phcpy.phcpy2c import py2c_copy_target_system_to_container
   from phcpy.phcpy2c import py2c_copy_start_system_to_container
   from phcpy.phcpy2c import py2c_copy_start_solutions_to_container
   from phcpy.phcpy2c import py2c_solcon_number_of_solutions
   from phcpy.solver import load_standard_system, load_standard_solutions
   from phcpy.trackers import standard_double_track
   cyclic7 = '/Users/jan/PHCv2/Demo/cyclic7'
   cyclic7q = '/Users/jan/PHCv2/Demo/cyclic7q'
   fail = py2c_read_standard_target_system_from_file(len(cyclic7),cyclic7)
   fail = py2c_copy_target_system_to_container()
   target = load_standard_system()
   fail = py2c_read_standard_start_system_from_file(len(cyclic7q),cyclic7q)
   fail = py2c_copy_start_system_to_container()
   start = load_standard_system()
   fail = py2c_copy_start_solutions_to_container()
   sols = load_standard_solutions()
   print 'number of start solutions :', py2c_solcon_number_of_solutions()
   print 'starting the path tracking with', nbtasks, 'task(s) ...'
   endsols = standard_double_track(target, start, sols, 0, nbtasks)
   print 'tracked', len(endsols), 'solution paths'
   
The documentation strings of the functions
exported by the module ``trackers`` of the package phcpy are listed below.

.. automodule:: trackers
   :members:

some interesting examples
-------------------------

PHCpack has been tested on many examples of polynomial systems
taken from the research literature.
The module examples exports some of those examples.
Running **python examples.py** at the command prompt
performs a regression test, solving all examples.

An interactive use of examples.py at the Python prompt can go as follows:

::

   >>> from phcpy.examples import noon3
   >>> f = noon3()
   >>> for p in f: print p
   ... 
   x1*x2^2 + x1*x3^2 - 1.1*x1 + 1;
   x2*x1^2 + x2*x3^2 - 1.1*x2 + 1;
   x3*x1^2 + x3*x2^2 - 1.1*x3 + 1;
   >>> 

The functions in examples.py returns the polynomials as lists of strings.
If we want to solve the system defined by f, we continue the above session as

::

   >>> from phcpy.solver import solve
   >>> s = solve(f,silent=True)
   >>> len(s)
   21
   >>> print s[0]
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x1 : -1.65123467890611E-01  -7.61734168646636E-01
    x2 :  8.98653694263692E-01  -3.48820047576431E-01
    x3 :  8.98653694263692E-01  -3.48820047576431E-01
   == err :  3.034E-16 = rco :  2.761E-01 = res :  5.974E-16 =
   >>> 

The example session continues in the description of the module solutions.

The documentation strings of the functions that
return the polynomials of the example systems as strings of characters
are listed below.
The regression test is exported by the function **test()**
of the module ``examples``.

.. automodule:: examples
   :members:

some families
-------------

Polynomial systems often occur in families and are defined
for any number of equations and variables.
One such noteworthy family is the cyclic n-roots problem:

::

   >>> from phcpy.families import cyclic
   >>> c4 = cyclic(4)
   >>> for p in c4: print p
   ... 
   x0 + x1 + x2 + x3;
   x0*x1 + x1*x2 + x2*x3 + x3*x0;
   x0*x1*x2 + x1*x2*x3 + x2*x3*x0 + x3*x0*x1;
   x0*x1*x2*x3 - 1;
   >>> 

.. automodule:: families
   :members:

numerical Schubert calculus
---------------------------

The module schubert.py exports the hypersurface and quantum
Pieri homotopies to solve the following Schubert problem:
Given a sequence of generic m-planes and a corresponding
sequence of interpolation points, compute all maps of
degree q that meet the given m-planes nontrivially at
the interpolation points.

The Pieri homotopies illustrates the homotopy approach.

1. Based on the dimension of the input problem, there is a formal
   root count on the number of solutions, a root count that is 
   exact for sufficiently generic instances of the input; and
   an upper bound for the number of isolated solution in all cases.

2. For sufficiently generic instances of the input, the performance
   of homotopies is optimal in the sense that every solution path
   defined by the homotopies ends at an actual solution of the problem.

The methods exported by the schubert module do the following:

1. Compute the formal root count for any m, p, and q.
   This calculation goes fast and gives an impression on the
   hardness of the problem.

2. Generate random instances of the problem for any (m,p,q).

3. Compute all solutions with the Pieri homotopies.

4. Verify the solutions with solutions.

5. Generate a instance of the problem known to be fully real.

The session below runs the Pieri homotopies to compute all linear maps
that produce 2-planes meeting 8 given 2-planes at random interpolation points:

::

   >>> from phcpy.schubert import *
   >>> (m,p,q) = (2,2,1)
   >>> n = m*p + q*(m+p)
   >>> r = Pieri_root_count(m,p,q)
   Pieri root count for (2, 2, 1) is 8
   the localization poset :
   n = 0 : ([3 4],[3 4],1)([2 5],[2 5],1)
   n = 1 : 
   n = 2 : ([2 4],[3 5],2)
   n = 3 : 
   n = 4 : ([2 3],[3 6],2)([2 3],[4 5],2)([1 4],[3 6],2)([1 4],[4 5],2)
   n = 5 : 
   n = 6 : ([1 3],[4 6],8)
   n = 7 : 
   n = 8 : ([1 2],[4 7],8)
   >>> L = [random_complex_matrix(m+p,m) for k in range(n)]
   >>> points = random_complex_matrix(n,1)
   >>> (f,fsols) = run_Pieri_homotopies(m,p,q,L,points)

The function **test()** of the module ``schubert``
runs an interactive session to solve instances that
are fully real (in case q = 0).

With the Littlewood-Richardson homotopies we can solve general Schubert problems.
The input to a Schubert problem is a sequence of n-by-n matrices and a
corresponding list of intersection conditions, represented by brackets.
For example, the bracket [2, 4, 6] imposes on a 3-plane in 6-space that it
meets nontrivially the space spanned by the first two columns of the corresponding
matrix in a line and that it meets the space spanned by the first four columns of
the corresponding matrix in a 2-plane.

For a generic sequence of input matrices, there are exactly two 3-planes
in 6-space that satisfy the conditions imposes by the three brackets
[2,4,6], [2,4,6], and [2,4,6], as computed in the session below.

::

   >>> from phcpy.schubert import resolve_schubert_conditions as rsc
   >>> brackets = [[2,4,6],[2,4,6],[2,4,6]]
   >>> rsc(6,3,brackets)
     the dimension of the planes : 3
     the number of conditions : 3
   [2 4 6] and [2 4 6] form a happy configuration.
   +1[2 4 6]*[2 4 6] = +1[2 3 4]+2[1 3 5]+1[1 2 6]
   [2 3 4] and [2 4 6] are not happy and will not create any children.
   [1 3 5] and [2 4 6] are happy and will create children...
   [1 2 6] and [2 4 6] are not happy and will not create any children.
   The new formal equations : 
   +2[1 3 5]*[2 4 6] = +2[1 2 3]
   All formal equations in the intersection poset :
   +1[2 4 6]*[2 4 6] = +1[2 3 4]+2[1 3 5]+1[1 2 6]
   +2[1 3 5]*[2 4 6] = +2[1 2 3]
   The intersection condition resolved :
   +1[2 4 6]*[2 4 6]*[2 4 6]
    = (+1[2 3 4]+2[1 3 5]+1[1 2 6])*[2 4 6]
    = +2[1 2 3]
   2
   >>>

To compute the 2-planes, we run the Littlewood-Richardson homotopies,
continuing the session from above:

::

   >>> from phcpy.schubert import littlewood_richardson_homotopies as lrh
   >>> (count, flags, sys, sols) = lrh(6,3,brackets,verbose=False)
   >>> count
   2
   >>> for sol in sols: print sol
   ... 
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x11 : -1.95764646993258E-01   1.07253045769427E+00
    x32 :  2.69552376387238E-01  -4.99588315456159E-01
    x53 : -9.21562255223665E-01  -9.28437273121748E-01
   == err :  3.559E-16 = rco :  4.125E-02 = res :  7.772E-16 =
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x11 : -5.85142692165828E-01  -9.14161314393362E-02
    x32 : -5.16006715209336E-01   3.41609194636644E-01
    x53 : -6.60253695726872E-02  -1.15262273262567E+00
   == err :  2.706E-13 = rco :  9.880E-02 = res :  4.219E-15 =
   >>>  len(sys)
   13
   >>>

The Littlewood-Richardson homotopies computed two solutions of a system
of 13 equations in 3 unknowns.

.. automodule:: schubert
   :members:

positive dimensional solution sets
----------------------------------

The module sets.py provides some functionality of PHCpack
to work with positive dimensional solution sets.

A witness set is a data structure to represent a positive dimensional
solution set.  A witness set consists of an embedding of the polynomial
equations that define the solution set, augmented with as many generic
linear equations as the dimenion of the solution set.
Witness points are solutions in the intersection of the original
polynomial equations and the generic linear equations.
The number of witness points equals the degree of the solution set.

In the example below we consider the twisted cubic:

::

   >>> twisted = ['x^2 - y;', 'x^3 - z;']
   >>> from phcpy.sets import embed
   >>> e = embed(3,1,twisted)
   >>> e[0]
   'x^2 - y + (-8.23538851649530E-01-5.67259869745581E-01*i)*zz1;'
   >>> e[1]
   'x^3 - z + (9.35464826338593E-01-3.53419805165623E-01*i)*zz1;'

The last equation of the embedded system is a linear equation
with randomly generated complex coefficient.  Continuing the session:

::

   >>> terms = e[-1].split(')*')
   >>> for t in terms: print t
   ... 
    + (-8.85038627286137E-01 + 4.65517591731472E-01*i
   x + (-2.12324313395875E-02 + 9.99774566519578E-01*i
   y + (-9.52478263619098E-01-3.04606561539880E-01*i
   z + (-9.59619713608467E-01 + 2.81300560351385E-01*i
   zz1+(-3.24025444378001E-01 + 9.46048366308847E-01*i);
   >>> from phcpy.solver import solve
   >>> s = solve(e, silent=True)
   >>> len(s)
   3
   >>> for sol in s: print sol
   ... 
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -4.06510360753325E-01   5.24436731997959E-01
    y : -1.09783212468901E-01  -4.26377930233570E-01
    zz1 : -4.27642353614751E-50  -2.73691106313441E-48
    z :  2.68236261633139E-01   1.15752697061077E-01
   == err :  3.693E-16 = rco :  1.041E-01 = res :  1.804E-16 =
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x :  7.97989261058868E-01   1.37578286563110E+00
    y : -1.25599163259885E+00   2.19571990464483E+00
    zz1 :  1.44884468793274E-32  -5.03715049801216E-33
    z : -4.02310165732919E+00   2.41891366942435E-02
   == err :  1.240E-15 = rco :  1.463E-02 = res :  2.120E-15 =
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -1.07164436617733E-01  -9.41488516596475E-01
    y : -8.74916410407434E-01   2.01788172926252E-01
    zz1 :  0.00000000000000E+00   0.00000000000000E+00
    z :  2.83741171803972E-01   8.02099237512644E-01
   == err :  9.857E-17 = rco :  8.220E-02 = res :  1.110E-16 =

The variable ``zz1`` is an artificial slack variable.
Adding the slack variable via an embedding is a general technique
to make overdetermined polynomial systems *square*,
that is: having as many equations as unknowns.
Only solutions with zero slack variables matter.

.. automodule:: sets
   :members:

monomial maps
-------------

Systems that have exactly two monomials with nonzero coefficient
in every equation are called binomial systems.
Although such binomial systems are very particular,
because of their sparse structure, they can be solved much faster.

The irreducible components of
positive dimensional solution sets of binomial systems
have coordinates that can be represented by maps of monomials 
in free independent variables.  In this representation, there
are as many free variables as the dimension of the solution set.
The module ``maps`` exports a solver for binomial systems.

In the example below, we consider a simple system
of two binomials in three variables:

::

   >>> f = [ 'x**2*y - z*x;', 'x**2*z - y**2*x;' ]
   >>> from phcpy.maps import binomial_solver
   >>> from phcpy.maps import solve_binomials
   >>> maps = solve_binomials(3,f)
   >>> for map in maps: print map
   ... 
   ['x - 0', 'y - (1+0j)*t1**1', 'z - (1+0j)*t2**1', 'dimension = 2', 'degree = 1']
   ['x - (1+0j)*t1**1', 'y - (1+0j)*t1**2', 'z - (1+0j)*t1**3', 'dimension = 1', 'degree = 3']
   ['x - (1+0j)*t1**1', 'y - 0', 'z - 0', 'dimension = 1', 'degree = 1']
   >>> 

In the output above we recognize the twisted cubic,
the x-axis, and the yz-plane as the three solution sets.

.. automodule:: maps
   :members:

Newton polytopes
----------------

The Newton polytopes of the polynomial system provide important
information about the structure of the solution sets.
The module ``polytopes`` provides an interface to the convex hull
methods of PHCpack.  It also provides a directer interface to the
mixed volume calculator, directer in the sense that the user can enter
the supports directly, without having to formulate a polynomial system.

The session below illustration the calculation of the convex hull
of a configuration of seven points in the plane.
The points are generated at random, with coordinates between -9 and +9.

::

   >>> from phcpy.polytopes import random_points as rp
   >>> from phcpy.polytopes import planar_convex_hull as pch
   >>> points = rp(2,7,-9,9)
   >>> points
   [(9, 8), (5, 6), (6, 0), (2, -5), (4, -1), (9, -4), (-1, -6)]
   >>> (vertices, normals) = pch(points)
   >>> vertices
   [(9, 8), (5, 6), (-1, -6), (9, -4)]
   >>> normals
   [(1, -2), (2, -1), (-1, 5), (-1, 0)]
   >>> 

The output of the convex hull method consists of a tuple of two lists.
The first list is the list of vertices.  For this particular example,
seven points were given on input, and only four of those points are 
corners of the convex hull.  The list of vertices is ordered cyclically:
two consecutive vertices span an edge of the polygon and the last and
first vertex also span an edge as a polygon has exactly as many vertices
as edges.  The second list in the output is the list of inner normals,
which are vectors perpendicular to the edges.  
Taking the inner product of the normal with the points that span an edge
yields the same value for each point on the edge, and that value is minimal
for all points in the polygon.  For the example above 
for the inner normal (1, -2) and the two points (9, 8) and (5, 6),
we have

.. math::

    9 - 2 \times 8 = 5 - 2 \times 6 = -7 

as the edge lies on the edge of the half plane defined by
the inequality

.. math::

    x_1 - 2 x_2 \geq -7

which holds for all points in the polygon spanned by the points
in the example.  The inner normals define the half planes that cut
out the polygon.

For a convex hull of a point configuration in 3-space, 
consider the example in the session below:

::

   >>> from phcpy.polytopes import random_points as rp
   >>> points = rp(3,10,-9,9)
   >>> for point in points: print point
   ... 
   (5, 9, -5)
   (0, 0, 1)
   (-3, -4, -1)
   (-9, -3, -3)
   (-5, 3, -8)
   (-4, 3, 7)
   (2, -3, 8)
   (9, 3, -9)
   (7, 4, -2)
   (1, -8, 1)
   >>> from phcpy.polytopes import convex_hull as ch
   >>> facets = ch(3, points)
   computed 12 facets
   >>> for facet in facets: print facet
   ... 
   (-597, [90, -65, -6], [4, 5, 6], [1, 2, 3])
   (-84, [1, 11, 14], [4, 8, 5], [7, 11, 0])
   (-281, [30, -49, -2], [5, 1, 6], [11, 4, 0])
   (-51, [6, 5, -6], [4, 6, 7], [0, 4, 5])
   (-203, [-22, -27, -30], [6, 1, 7], [2, 9, 3])
   (-48, [5, 6, -5], [4, 7, 10], [3, 6, 7])
   (-684, [-127, 66, -29], [10, 7, 8], [5, 8, 7])
   (-150, [1, 22, 25], [4, 10, 8], [5, 6, 1])
   (-315, [-59, 15, -19], [7, 9, 8], [9, 10, 6])
   (-265, [-29, -35, -39], [7, 1, 9], [4, 10, 8])
   (-165, [-19, -10, -4], [1, 8, 9], [11, 8, 9])
   (-429, [3, -26, 42], [5, 8, 1], [1, 10, 2])
   >>> 

The output of the ``convex_hull`` function returns a list of facets.
Each facet is represented as a tuple of four items.
The first number is the value of the inner product of the vector
perpendicular to the facet, given by the list in the second item
of the tuple.  So the first two items in the tuple define the
half space defined by the facet.  For the first facet, we have
the inequality defined by the number -597 and the vector [90, -65, -6]:

.. math::

   90 x_1 - 65 x_2 - 6 x_3 \geq -597

which holds for all points \ :math:`(x_1, x_2, x_3)` in the convex hull.  
The equality \ :math:`90 x_1 - 65 x_2 - 6 x_3 = -597` holds
for all points that lie on the first facet in the list of facets above.
The third item in the representation of a facet is the list of numbers
to the points that span the facet.  In the example above, the first
facet is spanned by the points 4, 5, 6 in the input list points.
Note that the counting of the points starts at one and not at zero.
The last item in the representation of a facet is the list of 
facets that are adjacent to the facet.  For the first facet,
facets 1, 2, and 3 are adjacent to it.  The counting of the facets
starts at zero, so the first facet has label zero.

From the list of facets we can extract all vertex points.
If we continue with the session from above:

::

   >>> vertices = []
   >>> for facet in facets:
   ...     for point in facet[2]:
   ...         if not point in vertices:
   ...             vertices.append(point)
   ... 
   >>> vertices
   [4, 5, 6, 8, 1, 7, 10, 9]
   >>> len(vertices)
   8

We have 8 vertices and 12 facets.  The points the span the facets are
ordered cyclically so that two consecutive points span an edge and the
last and first point span also an edge.  Every edge lies in the intersection
of exactly two facets.  Edges of adjacent facets are ordered in opposite
order.  For example, facet 0 is spanned by [4, 5, 6] and its adjacent
facet 1 is spanned by [4, 8, 5], with the edge shared between both of
them oriented from 4 to 5 in facet 0 and from 5 to 4 in facet 1.

As the points in the configuration were generated sufficiently at
random, the polytope is simplicial: every facet is spanned by exactly
3 points and has exactly 3 edges.  As every edge is shared by exactly
two facets we count every edge twice if we multiply the number of facets
by three, so we have 36/2 = 18 edges.

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

The function ``mixed_volume`` expects two arguments.
The first argument is the list of exponents of
the \ :math:`\lambda` variables in the volume expansion formula.
The second argument of ``mixed_volume`` is a tuple of Newton polytopes.
The session below illustrates the computation of the volume of one
single polytope.

::

   >>> from phcpy.polytopes import random_points as rp
   >>> from phcpy.polytopes import mixed_volume as mv
   >>> p1 = rp(3, 5, -9, 9)
   >>> p1
   [(3, 7, -3), (-1, 0, 8), (-6, -6, 8), (-6, 9, 4), (-3, 4, -7)]
   >>> tp1 = tuple([p1])
   >>> mv([3], tp1)
   2107

The volume is normalized, so the standard unit simplex has volume one.
To compute mixed volumes of two polytopes, we continue the session,
generating another polytope:

::

   >>> p2 = rp(3, 5, -9, 9)
   >>> mv([2, 1],(p1, p2))
   3910
   >>> mv([1, 2],(p1, p2))
   3961
   >>> 



.. automodule:: polytopes
   :members:

the module phcwulf
------------------

The file phcwulf defines a simple client/server interaction
to solve many random trinomials.

.. automodule:: phcwulf
   :members:

The module phcpy.phcpy2c
------------------------

The Python scripts in the package phcpy call the wrappers
for the C interface to PHCpack.  Below is the list of all
functions exported by the shared object file phcpy2c.so.
The source code provides more detailed documentation.

A basic application of the primitive operations in phcpy2c
is an interactive reading of a polynomial system.
Assume the file example at /tmp/ contains a polynomial system,
then we can do the following:

::

   >>> from phcpy.phcpy2c import py2c_syscon_read_system as readsys
   >>> from phcpy.phcpy2c import py2c_syscon_write_system as writesys
   >>> readsys()

   Reading a polynomial system...
   Is the system on a file ? (y/n/i=info) y

   Reading the name of the input file.
   Give a string of characters : /tmp/example
   0
   >>> writesys()
    2
   x^2+4*y^2-4;
   2*y^2-x;
   0
   >>> from phcpy.phcpy2c import py2c_solve_system as solve
   >>> solve()

   ROOT COUNTS :

   total degree : 4
   general linear-product Bezout number : 4
     based on the set structure :
        { x y }{ x y }
        { x y }{ y }
   mixed volume : 4
   stable mixed volume : 4
   4
   >>> from phcpy.phcpy2c import py2c_solcon_write_solutions as writesols
   >>> writesols()
   4 2
   ===========================================================================
   solution 1 :
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x :  1.23606797749979E+00  -9.91383530201425E-119
    y :  7.86151377757423E-01   4.95691765100713E-119
   == err :  1.567E-16 = rco :  3.067E-01 = res :  3.331E-16 ==
   solution 2 :
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x :  1.23606797749979E+00  -9.91383530201425E-119
    y : -7.86151377757423E-01  -4.95691765100713E-119
   == err :  1.567E-16 = rco :  3.067E-01 = res :  3.331E-16 ==
   solution 3 :
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -3.23606797749979E+00   3.17242729664456E-117
    y : -1.58621364832228E-117   1.27201964951407E+00
   == err :  3.703E-16 = rco :  1.515E-01 = res :  4.441E-16 ==
   solution 4 :
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -3.23606797749979E+00   3.17242729664456E-117
    y :  1.58621364832228E-117  -1.27201964951407E+00
   == err :  3.703E-16 = rco :  1.515E-01 = res :  4.441E-16 ==
   0
   >>> 

With these primitive operations in phcpy2c we can bypass the writing
and the parsing to strings.

.. automodule:: phcpy2c
   :members:

Project History
===============

The Python interface to PHCpack got to a first start when
Kathy Piret met William Stein at the software for algebraic geometry
workshop at the IMA in the Fall of 2006.  
The first version of this interface is described
in the 2008 PhD Thesis of Kathy Piret.

The implementation of the Python bindings depend on the C interface
to PHCpack, developed for use with message passing on distributed
memory computers.

Version 0.0.1 originated at lecture 40 of MCS 507 in the Fall of 2012,
as an illustration of Sphinx.  In Spring of 2013, version 0.5.0 was
presented at a graduate computational algebraic geometry seminar.
Version 0.1.0 was prepared for presentation at EuroSciPy 2013 (August 2013).
Improvements using pylint led to version 0.1.1
and the module maps was added in version 0.1.2.
Version 0.1.4 added path trackers with a generator
so all solutions along a path are returned to the user.
Multicore path tracking was added in version 0.1.7.

The paper **Modernizing PHCpack through phcpy**
written for the EuroSciPy 2013 proceedings 
and available at <http://arxiv.org/abs/1310.0056>
describes the design of phcpy.

References
==========

1. T. Gao, T. Y. Li, M. Wu:
   **Algorithm 846: MixedVol: a software package for mixed-volume 
   computation.**
   *ACM Transactions on Mathematical Software*, 31(4):555-560, 2005.

2. Y. Hida, X.S. Li, and D.H. Bailey:
   **Algorithms for quad-double precision floating point arithmetic.**
   In *15th IEEE Symposium on Computer Arithmetic (Arith-15 2001)*,
   11-17 June 2001, Vail, CO, USA, pages 155-162.
   IEEE Computer Society, 2001.
   Shortened version of Technical Report LBNL-46996.

3. A. Leykin and J. Verschelde.
   **Interfacing with the numerical homotopy algorithms in PHCpack.**
   In N. Takayama and A. Iglesias, editors, *Proceedings of ICMS 2006*,
   volume 4151 of *Lecture Notes in Computer Science*,
   pages 354--360. Springer-Verlag, 2006.

4. K. Piret:
   **Computing Critical Points of Polynomial Systems 
   using PHCpack and Python.**
   PhD Thesis, University of Illinois at Chicago, 2008.

5. A. J. Sommese, J. Verschelde, and C. W. Wampler.
   **Numerical irreducible decomposition using PHCpack.**
   In *Algebra, Geometry, and Software Systems*, 
   edited by M. Joswig and N. Takayama,
   pages 109-130. Springer-Verlag, 2003.

6. J. Verschelde:
   **Algorithm 795: PHCpack: A general-purpose solver for polynomial
   systems by homotopy continuation.**
   *ACM Transactions on Mathematical Software*, 25(2):251--276, 1999.

7. J. Verschelde:
   **Modernizing PHCpack through phcpy.**
   In Proceedings of the 6th European Conference on Python in Science
   (EuroSciPy 2013), edited by Pierre de Buyl and Nelle Varoquaux,
   pages 71-76, 2014, available at
   <http://arxiv.org/abs/1310.0056>.

Acknowledgments
===============

This material is based upon work supported by the 
National Science Foundation under Grant 1115777.
Any opinions, findings, and conclusions or recommendations expressed 
in this material are those of the author(s) and do not necessarily 
reflect the views of the National Science Foundation. 

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

