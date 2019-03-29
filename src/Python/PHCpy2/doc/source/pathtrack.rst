path trackers and sweep homotopies
==================================

Homotopy continuation methods are applied to solve a polynomial system.
The module **phcpy.trackers** exports the path trackers of PHCpack.
The functions in this module track paths defined by artificial-parameter
homotopies, of the form

.. math::

    h(x,t) = \gamma (1-t) g(x) + t f(x) = 0,

where \ :math:`\gamma` is a randomly generated complex constant.
The artificial parameter \ :math:`t` goes from zero to one,
from the known solutions of the start system \ :math:`g(x) = 0`
to the solutions of the target system \ :math:`f(x) = 0`.

The module **phcpy.sweepers** exports the sweep homotopies.
A sweep homotopy is a natural parameter homotopy.
Its application is to track solution paths from one set of values
for the parameters to another set of values for the parameters.

The tracking of solution paths defined by an artificial-parameter homotopy
apply an increment-and-fix method: the continuation parameter \ :math:`t`
is incremented by the predictor and remains fixed in the corrector.
The tracking of solution paths defined by a sweep homotopy apply
arc length parameter continuation.

a simple example
----------------

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
   >>> p = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
   >>> d = total_degree(p)
   >>> d
   4
   >>> (q, qsols) = total_degree_start_system(p)
   >>> len(qsols)
   4
   >>> q
   ['x^2 - 1;', 'y^2 - 1;']
   >>> s = track(p, q, qsols)
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
   >>> p = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
   >>> (q,qsols) = total_degree_start_system(p)
   >>> s1 = track(p, q, [qsols[2]])
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

fixing the gamma constant
-------------------------

To avoid this side effect, *track* accepts a complex value
as its last argument for the so-called gamma constant.
As a continuation of the session from above:

::

   >>> s3 = track(p, q, [qsols[2]], gamma=complex(0.824372806319,0.56604723848934))
   >>> print s3[0]
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -3.23606797749979E+00   0.00000000000000E+00
    y :  0.00000000000000E+00   1.27201964951407E+00
   == err :  0.000E+00 = rco :  1.079E-01 = res :  0.000E+00 =

If we track all solution paths one after the other,
each time calling track with the same value for gamma,
then all solutions will be found.

give the next solution on a path
--------------------------------

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
   >>> (q, s) = total_degree_start_system(p)
   >>> from phcpy.trackers import initialize_standard_tracker
   >>> from phcpy.trackers import initialize_standard_solution
   >>> from phcpy.trackers import next_standard_solution
   >>> initialize_standard_tracker(p, q)
   >>> initialize_standard_solution(len(p), s[0])
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

In the session above, we see the solutions ``s1`` for t = 0.1,
and two other solutions for consecutive values 0.2 and 0.3 for t.
If we continue the session from above with the second start solution
in ``s[1]``, we can select the first 11 points along the path
and view all values for ``x`` of the solutions:

::

   >>> initialize_standard_solution(len(p), s[1])
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

We see that the last two values differ little from each other
because we arrived at the end of the path.  
To test whether at the end of a path, it suffices to check
whether the value for t equals one.

The real parts of the four paths are shown in :numref:`figplotpaths`.
Three of the paths converge to the triple solution (1,2).

.. _figplotpaths:

.. figure:: ./plotpaths.png
    :align: center

    The real parts of four solution paths.

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
           tval = dictsol['t'].real
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

With ``False`` in
``initialize_standard_tracker(p, q, False)``
the option to generate a fixed gamma constant is turned off,
so rerunning the same code will generate other random constants
and produce different plots.

To set the value of the gamma constant to a specific value,
e.g.: ``(-0.853618933016554-0.52089799116111j)``, do the following.

::

    gamma = (-0.853618933016554-0.52089799116111j)
    initialize_standard_tracker(target, start, False, gamma.real, gamma.imag)

solving with polyhedral homotopies
----------------------------------

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

Newton's method at higher precision
-----------------------------------

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

Looking at the values for *err* and *res* we see huge values
for two solutions which are spurious.

multitasked path tracking
-------------------------

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

GPU accelerated path tracking
-----------------------------

The script below illustrates the call to the GPU accelerated path
trackers.  As input, the location of a random coefficient system
(as obtained via the polyhedral homotopies) is needed.
With this random coefficient system, we have an artificial-parameter
homotopy that defines 35,940 solution paths to solve the
cyclic 10-roots problem.

::

   GPU = 1 # use the GPU
   DIR = '/home/jan/Problems/GPUdata/MultiPath' # location of systems
   from phcpy.phcpy2c \
   import py2c_read_standard_target_system_from_file as read_target
   from phcpy.phcpy2c \
   import py2c_read_standard_start_system_from_file as read_start
   cyc10tarfile = DIR + '/cyclic10.target'
   cyc10stafile = DIR + '/cyclic10.start'
   fail = read_target(len(cyc10tarfile), cyc10tarfile)
   from phcpy.interface import load_standard_system as loadsys
   from phcpy.interface import load_standard_solutions as loadsols
   cyc10 = loadsys()
   print 'the cyclic 10-roots problem :'
   for pol in cyc10:
       print pol
   fail = read_start(len(cyc10stafile), cyc10stafile)
   cyc10q = loadsys()
   print 'a start system for the cyclic 10-roots problem :'
   for pol in cyc10q:
       print pol
   cyc10qsols = loadsols()
   print 'number of start solutions :', len(cyc10qsols)
   print 'the first solution :'
   print cyc10qsols[0]
   print 'calling the path tracker...'
   if(GPU == 0):
       from phcpy.trackers import ade_double_track
       cyc10sols = ade_double_track(cyc10,cyc10q,cyc10qsols,verbose=0)
   else:
       from phcpy.trackers import gpu_double_track
       cyc10sols = gpu_double_track(cyc10,cyc10q,cyc10qsols,verbose=0)
   print 'number of solutions :', len(cyc10sols)
   for sol in cyc10sols:
       print sol

sweep homotopies
----------------

A *sweep homotopy* is a family of polynomial systems with a least one
natural parameter and one artificial parameter.
As the artificial parameter moves from zero to one, the natural parameter
changes from a given start value to another given target value.
By arc length continuation, the solution paths are tracked from the
given start values for the parameters to the target values.

Consider a simple example: sweeping the circle.
We consider the unit circle \ :math:`x^2 + y^2 - 1 = 0`,
intersected by a horizontal line, at the start equal to \ :math:`y = 0`.
In a Python session, we could define the sweep homotopy that takes
the line from \ :math:`y = 0` to \ :math:`y = 2`.

::

   >>> circle = ['x^2 + y^2 - 1;', 'y*(1-s) + (y-2)*s;']

For ``s = 0`` there are two solutions,
with values for ``x`` and ``y`` in the 
tuples \ :math:`(1, 0)` and \ :math:`(-1,0)`.

Geometrically, as the horizontal line moves up, the two solutions
(the intersection points on the circle and the line), move closer
to each other to join at a *quadratic turning point*,
shown in :numref:`figcircleline`.
At the left picture of :numref:`figcircleline` we see the line transversally 
intersecting the circle at a perfect right angle.  
At the right picture of :numref:`figcircleline`, the two distinct
solutions have merged into one point where the line is tangent to the circle.

.. _figcircleline:

.. figure:: ./circles.png
    :align: center

    Two complex conjugated solutions meet at a quadratic turning point.

The tracking of solution paths in a real sweep homotopy will stop
at the first singular point it encounters.  
The continuation of the code with the definition of ``circle``
to launch this path tracking is listed below:

::

   >>> from phcpy.solutions import make_solution as makesol
   >>> first = makesol(['x', 'y', 's'], [1, 0, 0])
   >>> second = makesol(['x', 'y', 's'], [-1, 0, 0])
   >>> startsols = [first, second]
   >>> from phcpy.sweepers import standard_real_sweep as sweep
   >>> newsols = sweep(circle, startsols)
   >>> print newsols[0]

and then we see as output of the ``print`` statement:

::

   t :  0.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -2.46519032881566E-32   0.00000000000000E+00
    y :  1.00000000000000E+00   0.00000000000000E+00
    s :  5.00000000000000E-01   0.00000000000000E+00
   == err :  0.000E+00 = rco :  1.000E+00 = res :  0.000E+00 =

The sweep stopped where ``s`` is equal to 0.5,
with corresponding values for ``x`` and ``y`` in the 
tuple \ :math:`(0, 1)`.

real versus complex sweeps
--------------------------

In a *complex sweep*, an addition random gamma constant is generated
in the convex-linear combination between the sets of start and target
values for the parameters.  If the solutions for the start values of
the parameters are regular, then the application of the *gamma trick*
applies for problems where the parameter space is convex.
This means that, if the problem formulation makes sense for convex
combinations of the parameters, then the solution paths will remain
regular, except for finitely many bad choices of the random gamma constant,
and except perhaps at the very end of the paths, when the target values
for the parameters lead to polynomial systems with singular solutions.

Conducting a complex sweep on the circle can be done as follows:

::

    >>> circle = ['x^2 + y^2 - 1;']
    >>> from phcpy.solutions import make_solution as makesol
    >>> first = makesol(['x', 'y'], [1, 0])
    >>> second = makesol(['x', 'y'], [-1, 0])
    >>> startsols = [first, second]
    >>> par = ['y']
    >>> start = [0, 0] 
    >>> target = [2, 0]
    >>> from phcpy.sweepers import standard_complex_sweep as sweep
    >>> newsols = sweep(circle, startsols, 2, par, start, target)

The setup of the homotopy defines ``y`` as the parameter 
(in the list ``['y']`` assigned to ``par``).
The parameter \ :math:`y` will move from the complex zero \ :math:`0 + 0 I` 
(given by the list ``[0, 0]`` assigned to ``start``) 
to \ :math:`2 + 0 I`
(given by the list ``[2, 0]`` assigned to ``target``).
The corresponding start solutions for \ :math:`y = 0` are stored
in the tuples \ :math:`(1,0)` and \ :math:`(-1,0)`.
Then, at the end of the sweep, we will find two complex conjugated solutions.

::

    >>> print newsols[0]
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x : -6.27554627321419E-26   1.73205080756888E+00
     y :  2.00000000000000E+00   0.00000000000000E+00
    == err :  6.642E-13 = rco :  1.000E+00 = res :  4.441E-16 =

What is now the difference between real versus complex?
The *real sweep* stopped at the singular solution \ :math:`(0,1)`
while the *complex sweep* hopped over this singularity because
of complex random gamma constant in the convex combination between
the start and target values of the parameters.

tuning parameters, settings, and tolerances
-------------------------------------------

The default values of the numerical parameters were set based
on computational experiences on a large, representative collection
of polynomial systems.  The module ``tuning`` provides functions
to adjust the parameters, settings, and tolerances.
The function ``tune_track_parameters`` gives access to the tuning 
as in ``phc -p``, via an interactive menu.
The other functions in the module allow to get the values and to
set the values of each parameter, setting, or tolerance.

a polyhedral end game
---------------------

In case the mixed volume is not a sharp root count,
there are paths diverging to points with coordinates equal to zero,
or diverging to infinity. 
The directions of those diverging paths coincide with
the leading exponents of the Puiseux series expansions of the points
with coordinates equal to zero and/or at infinity.
In particular, positive leading exponents occur with coordinates
going to zero, while for a coordinate at infinity, the corresponding
leading exponent will be negative.

To activate the polyhedral end game, the extrapolation order needs
to be nonzero.  We can set this order as follows:

::

   >>> from phcpy.tuning import order_endgame_extrapolator_set as set
   >>> set(4)
   0

The ``0`` on return is the failure code, which should equal zero
if all went well.  To double check, we can get the value of the
order of the extrapolator in the end game:

::

   >>> from phcpy.tuning import order_endgame_extrapolator_get as get
   >>> get()
   4

Let us run a polyhedral end game on a very simple example.

::

   >>> f = ['x + y^3 - 1;', 'x + y^3 + 1;']
   >>> from phcpy.solver import mixed_volume as mv
   >>> from phcpy.solver import random_coefficient_system as rcs
   >>> mv(f)
   4
   >>> (g, gsols) = rcs(f)
   >>> len(gsols)
   4

Although the mixed volume equals four (and we have four start solutions
in ``gsols`` of the start system ``g``), we can see that ``f`` has no
solutions, and all four paths will diverge to infinity.

::

   >>> from phcpy.trackers import standard_double_track as track
   >>> sols = track(f, g, gsols)
   >>> from phcpy.tropisms import standard_retrieve as retrieve
   >>> (w, d, e) = retrieve(len(sols), len(f))
   >>> w
   [3, 3, 3, 3]

We see that the winding numbers of the four paths are all equal to 3
and the numerically computed tropisms are approximations of (-1, -1/3),
or (-3, -1) when presented in normal form.
