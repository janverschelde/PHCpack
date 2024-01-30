********
Tutorial
********

Via some interesting use cases, several features are introduced.

1. If one can solve polynomial systems,
   then one can solve many problems.

2. The use cases illustrate the formulation of the polynomial systems
   via the packages in Python's computational ecosystem,
   in particular sympy and numpy.

3. The interpretation of the output often happens via plots.

For the plots, Jupyter notebooks are available in the ``tests`` folder.
The development of the use cases happened via short Python scripts.

A Counter Example to Koushnirenko's Conjecture
==============================================

This chapter illustrates the
counter example of Bertrand Haas, against the Koushnirenko conjecture,
executed on one core and on many cores.
For the mathematical background, consult:

Bertrand Haas: **A simple counterexample to Kouchnirenko's conjecture.**,
*Beitraege zur Algebra und Geometrie/Contributions to Algebra and Geometry*, 
volume 43, number 1, pages 1 to 8, 2002.

Bertrand constructed a polynomial system with three monomials in every
equation and with five positive real roots.  Kouchnirenko's conjecture
predicted there could only four positive real roots for such a system.

To get the proper wall clock time, we have to be mindful 
that the Python code calls the compiled functions in the PHCpack library.  
Therefore, the Python timers will not give accurate timings.  
Instead, we have to rely on the actual date and time, 
from the package `datetime` in Python."

::

    from datetime import datetime

For the plot, the implicit plotting of `sympy` will be used.

::

    from sympy import plot_implicit, symbols, Eq

From `phcpy` we import the following functions:

::

    from phcpy.dimension import get_core_count
    from phcpy.solver import solve
    from phcpy.solutions import filter_real

Solving the System on Many Cores
--------------------------------

The example of Bertrand Haas is defined as

::

    H = [ 'x**108 + 1.1*y**54 - 1.1*y;',
          'y**108 + 1.1*x**54 - 1.1*x;' ]

According to the theorem of Bézout, we may expect a number of 
complex solutions equals to the product of the degrees of the polynomials.
The square of 108 equals 11664.
As the solver computes all complex solutions,
executing the following code block takes some time ...

::

    print('Solving on one core ...')
    wstart = datetime.now()
    sols = solve(H)
    wstop = datetime.now()
    print('  Number of solutions :', len(sols))
    print('start time :', wstart)
    print(' stop time :', wstop)
    print('   elapsed :', wstop - wstart)

The output of the above code cell is

::

     Solving on one core ...
       Number of solutions : 11664
     start time : 2024-01-28 11:57:53.061707
      stop time : 2024-01-28 11:57:59.223344
        elapsed : 0:00:06.161637

We can significantly speed up this computation if the computer has many cores.

::

    nbcores = get_core_count()
    print('Solving on', nbcores, 'cores ...')
    wstart = datetime.now()
    sols = solve(H, tasks=nbcores)
    wstop = datetime.now()
    print('  Number of solutions :', len(sols))
    print('start time :', wstart)
    print(' stop time :', wstop)
    print('   elapsed :', wstop - wstart)

The output of the above code cell is

::

    Solving on 32 cores ...
      Number of solutions : 11664
    start time : 2024-01-28 11:58:07.241324
     stop time : 2024-01-28 11:58:08.747874
       elapsed : 0:00:01.506550

Compared the `elapsed :` above with the previous one.

Extracting the Real Roots
-------------------------

Rather than eyeballing all 11,664 complex solutions ourselves, 
we ask to filter the real solutions.

::

    realsols = filter_real(sols, tol=1.0e-8, oper='select')

The code cell below prints the solutions in `realsols`:

::

    for (idx, sol) in enumerate(realsols):
        print('Solution', idx+1, ':')
        print(sol)

The output is

::

     Solution 1 :
     t :  0.00000000000000E+00   0.00000000000000E+00
     m : 1
     the solution for t :
      x :  0.00000000000000E+00   0.00000000000000E+00
      y :  0.00000000000000E+00   0.00000000000000E+00
     == err :  0.000E+00 = rco :  1.000E+00 = res :  0.000E+00 =
     Solution 2 :
     t :  1.00000000000000E+00   0.00000000000000E+00
     m : 1
     the solution for t :
      x :  9.91489402484465E-01  -2.94004118110142E-49
      y :  9.91489402484465E-01   2.96676882820234E-49
     == err :  7.704E-17 = rco :  8.274E-02 = res :  5.773E-15 =
     Solution 3 :
     t :  1.00000000000000E+00   0.00000000000000E+00
     m : 1
     the solution for t :
      x :  9.99997917489999E-01   9.52445049970774E-46
      y :  9.19904793199125E-01  -1.72639970804817E-42
     == err :  2.708E-16 = rco :  1.601E-03 = res :  8.677E-15 =
     Solution 4 :
     t :  1.00000000000000E+00   0.00000000000000E+00
     m : 1
     the solution for t :
      x :  9.99986016402972E-01   1.19248391761152E-37
      y :  9.36266084294562E-01  -2.97164971887874E-34
     == err :  2.240E-15 = rco :  2.070E-03 = res :  3.610E-15 =
     Solution 5 :
     t :  1.00000000000000E+00   0.00000000000000E+00
     m : 1
     the solution for t :
      x :  9.19904793199125E-01  -2.01786978862774E-41
      y :  9.99997917489999E-01   1.02032044433651E-44
     == err :  4.384E-16 = rco :  1.601E-03 = res :  8.564E-15 =
     Solution 6 :
     t :  1.00000000000000E+00   0.00000000000000E+00
     m : 1
     the solution for t :
      x :  9.36266084294562E-01  -4.40471624223194E-48
      y :  9.99986016402972E-01   3.13214614463929E-51
     == err :  9.256E-17 = rco :  2.070E-03 = res :  3.730E-15 =

We observe (0, 0) and five additional real positive roots.
According to the Koushnirenko conjecture, we would expect
no more than four real positive roots.

Plotting the Curves
-------------------

In converting the strings in the polynomial system `H` 
we have to remove the trailing semicolon

::

    x, y = symbols('x y')
    p0 = eval(H[0][:-1])
    p1 = eval(H[1][:-1])

Without knowing the precise location of the intersection points, 
the curves are hard to plot.
The code below produces the plot :numref:`haasfig1`.

::

    plot0 = plot_implicit(Eq(p0, 0), (x, 0.93, 1.01), (y, 0.93, 1.01),
                          line_color='black', depth=1,
        markers=[{'args': [[0.99148, 0.93626, 0.99998],
                           [0.99148, 0.99998, 0.93626], 'bo']}],
                          axis_center=(0.93, 0.93), show=False)

    plot1 = plot_implicit(Eq(p1, 0), (x, 0.93, 1.01), (y, 0.93, 1.01),
                          line_color='red', depth=1,
                          axis_center=(0.93, 0.93), show=False)

    plot0.append(plot1[0])
    plot0.show()

.. _haasfig1:

.. figure:: ./haasfig1.png
   :align: center
    
   Three positive roots of the counterexample.

Let us zoom in to another root ...
:numref:`haasfig2` is made executing the code below:

::

    plot3 = plot_implicit(Eq(p0, 0), (x, 0.9, 1.0025), (y, 0.99, 1.0025),
                          line_color='black', depth=1,
        markers=[{'args': [[0.99148, 0.93626, 0.91990],
                           [0.99148, 0.99998, 0.99999], 'bo']}],
                          axis_center=(0.9, 0.99), show=False)

    plot4 = plot_implicit(Eq(p1, 0), (x, 0.9, 1.0025), (y, 0.99, 1.0025),
                          line_color='red', depth=1,
                          axis_center=(0.9, 0.99), show=False)

    plot3.append(plot4[0])
    plot3.show()

.. _haasfig2:

.. figure:: ./haasfig2.png
   :align: center
    
   Other positive roots of the counterexample.

The third plot (in :numref:`haasfig3`)
is produced by the following code:

::

    plot5 = plot_implicit(Eq(p0, 0), (x, 0.9, 1.0005), (y, 0.999, 1.0005),
                          line_color='black', depth=1,
        markers=[{'args': [[0.93626, 0.91990],
                           [0.99998, 0.99999], 'bo']}],
                          axis_center=(0.9, 0.999), show=False)

    plot6 = plot_implicit(Eq(p1, 0), (x, 0.9, 1.0005), (y, 0.999, 1.0005),
                          line_color='red', depth=1,
                          axis_center=(0.9, 0.999), show=False)

    plot5.append(plot6[0])
    plot5.show()

.. _haasfig3:

.. figure:: ./haasfig3.png
   :align: center
    
   Two close positive roots of the counterexample.

A 4-Bar Mechanism
=================

The equations to design a 4-bar mechanism are defined with sympy.

The system appears in a paper by A.P. Morgan and C.W. Wampler on
**Solving a Planar Four-Bar Design Using Continuation**, published in
*the Journal of Mechanical Design*, volume 112, pages 544-550, 1990.

The solutions for a straight-line confirmation are shown with matplotlib.
Random numbers will be generated.

::

    from math import sqrt
    from random import uniform

From `sympy` we import the following:

::

    from sympy import var
    from sympy.matrices import Matrix

For the plotting, we import `pyplot` of `matplotlib`.

::

    import matplotlib.pyplot as plt

And then, last an not least, the blackbox solver
of `phcpy` is imported.

::

    from phcpy.solver import solve

As `phcpy` is an API, the problem is solved
via a sequence of functions.

solving a polynomial system
---------------------------

The system of polynomial equations is formulated
by the function in the code cell below.

::

    def polynomials(d0, d1, d2, d3, d4, a):
        """
        Given in d0, d1, d2, d3, d4 are the coordinates of
        the precision points, given as Matrix objects.
        Also the coordinates of the pivot in a are stored in a Matrix.
        Returns the system of polynomials to design the 4-bar
        mechanism with a coupler passing through the precision points.
        """
        # the four rotation matrices
        c1, s1, c2, s2 = var('c1, s1, c2, s2')
        c3, s3, c4, s4 = var('c3, s3, c4, s4')
        R1 = Matrix([[c1, -s1], [s1, c1]])
        R2 = Matrix([[c2, -s2], [s2, c2]])
        R3 = Matrix([[c3, -s3], [s3, c3]])
        R4 = Matrix([[c4, -s4], [s4, c4]])
        # the first four equations reflecting cos^2(t) + sin^(t) = 1
        p1, p2 = 'c1^2 + s1^2 - 1;', 'c2^2 + s2^2 - 1;'
        p3, p4 = 'c3^2 + s3^2 - 1;', 'c4^2 + s4^2 - 1;'
        # the second four equations on X
        x1, x2 = var('x1, x2')
        X = Matrix([[x1], [x2]])
        c1x = 0.5*(d1.transpose()*d1 - d0.transpose()*d0)
        c2x = 0.5*(d2.transpose()*d2 - d0.transpose()*d0)
        c3x = 0.5*(d3.transpose()*d3 - d0.transpose()*d0)
        c4x = 0.5*(d4.transpose()*d4 - d0.transpose()*d0)
        e1x = (d1.transpose()*R1 - d0.transpose())*X + c1x
        e2x = (d2.transpose()*R2 - d0.transpose())*X + c2x
        e3x = (d3.transpose()*R3 - d0.transpose())*X + c3x
        e4x = (d4.transpose()*R4 - d0.transpose())*X + c4x
        s1, s2 = str(e1x[0]) + ';', str(e2x[0]) + ';'
        s3, s4 = str(e3x[0]) + ';', str(e4x[0]) + ';'
        # the third group of equations on Y
        y1, y2 = var('y1, y2')
        Y = Matrix([[y1], [y2]])
        c1y = c1x - a.transpose()*(d1 - d0)
        c2y = c2x - a.transpose()*(d2 - d0)
        c3y = c3x - a.transpose()*(d3 - d0)
        c4y = c4x - a.transpose()*(d4 - d0)
        e1y = ((d1.transpose() - a.transpose())*R1 \
             - (d0.transpose() - a.transpose()))*Y + c1y
        e2y = ((d2.transpose() - a.transpose())*R2 \
             - (d0.transpose() - a.transpose()))*Y + c2y
        e3y = ((d3.transpose() - a.transpose())*R3 \
             - (d0.transpose() - a.transpose()))*Y + c3y
        e4y = ((d4.transpose() - a.transpose())*R4 \
             - (d0.transpose() - a.transpose()))*Y + c4y
        s5, s6 = str(e1y[0]) + ';', str(e2y[0]) + ';'
        s7, s8 = str(e3y[0]) + ';', str(e4y[0]) + ';'
        return [p1, p2, p3, p4, s1, s2, s3, s4, s5, s6, s7, s8]

Let us generate random points and define the polynomial system.


::

    pt0 = Matrix(2, 1, lambda i,j: uniform(-1,+1))
    pt1 = Matrix(2, 1, lambda i,j: uniform(-1,+1))
    pt2 = Matrix(2, 1, lambda i,j: uniform(-1,+1))
    pt3 = Matrix(2, 1, lambda i,j: uniform(-1,+1))
    pt4 = Matrix(2, 1, lambda i,j: uniform(-1,+1))
    # the pivot is a
    piv = Matrix([[1], [0]])
    equ = polynomials(pt0,pt1,pt2,pt3,pt4,piv)
    for pol in equ:
        print(pol)

Then the output for random numbers as the parameters is

::

      c1^2 + s1^2 - 1;
      c2^2 + s2^2 - 1;
      c3^2 + s3^2 - 1;
      c4^2 + s4^2 - 1;
      x1*(-0.275586755195824*c1 + 0.325788266703467*s1 + 0.676839821431551) + x2*(0.325788266703467*c1 + 0.275586755195824*s1 - 0.0938422352018522) - 0.142416227311081;\n",
      x1*(0.902513020087508*c2 - 0.151712455719013*s2 + 0.676839821431551) + x2*(-0.151712455719013*c2 - 0.902513020087508*s2 - 0.0938422352018522) + 0.185313955832297;\n",
      x1*(-0.44237719048943*c3 - 0.542955104453471*s3 + 0.676839821431551) + x2*(-0.542955104453471*c3 + 0.44237719048943*s3 - 0.0938422352018522) + 0.0117896575671138;\n",
      x1*(-0.319438148253377*c4 - 0.397350378412077*s4 + 0.676839821431551) + x2*(-0.397350378412077*c4 + 0.319438148253377*s4 - 0.0938422352018522) - 0.103495227599703;\n",
      y1*(-1.27558675519582*c1 + 0.325788266703467*s1 + 1.67683982143155) + y2*(0.325788266703467*c1 + 1.27558675519582*s1 - 0.0938422352018522) - 0.543669293546807;\n",
      y1*(-0.0974869799124924*c2 - 0.151712455719013*s2 + 1.67683982143155) + y2*(-0.151712455719013*c2 + 0.0974869799124924*s2 - 0.0938422352018522) - 1.39403888568676;\n",
      y1*(-1.44237719048943*c3 - 0.542955104453471*s3 + 1.67683982143155) + y2*(-0.542955104453471*c3 + 1.44237719048943*s3 - 0.0938422352018522) - 0.222672973375007;\n",
      y1*(-1.31943814825338*c4 - 0.397350378412077*s4 + 1.67683982143155) + y2*(-0.397350378412077*c4 + 1.31943814825338*s4 - 0.0938422352018522) - 0.460896900777877;\n"

The solutions of the polynomial system define a mechanism
of which the coupler passes through the five points.

::

    sols = solve(equ)
    len(sols)

The number is `36` which is invariant for this problem.
Solving a general problem, for random precision points, 
shows that the number of solutions is 36.

a straight-line configuration
-----------------------------

Let us consider a special problem.
Observe the extraction of real solutions in the function below.

::

    def straight_line(verbose=True):
        """
        This function solves an instance where the five precision
        points lie on a line.  The coordinates are taken from Problem 7
        of the paper by A.P. Morgan and C.W. Wampler.
        Returns a list of solution dictionaries for the real solutions.
        """
        from phcpy.solutions import strsol2dict, is_real
        pt0 = Matrix([[ 0.50], [ 1.06]])
        pt1 = Matrix([[-0.83], [-0.27]])
        pt2 = Matrix([[-0.34], [ 0.22]])
        pt3 = Matrix([[-0.13], [ 0.43]])
        pt4 = Matrix([[ 0.22], [ 0.78]])
        piv = Matrix([[1], [0]])
        equ = polynomials(pt0,pt1,pt2,pt3,pt4,piv)
        if verbose:
            print('the polynomial system :')
            for pol in equ:
                print(pol)
        sols = solve(equ)
        if verbose:
            print('the solutions :')
            for (idx, sol) in enumerate(sols):
                print('Solution', idx+1, ':')
                print(sol)
            print('computed', len(sols), 'solutions')
        result = []
        for sol in sols:
            if is_real(sol, 1.0e-8):
                soldic = strsol2dict(sol)
                result.append(soldic)
        return result

Running the function

::

    sols = straight_line()

shows

::

    the polynomial system :
    c1^2 + s1^2 - 1;
    c2^2 + s2^2 - 1;
    c3^2 + s3^2 - 1;
    c4^2 + s4^2 - 1;
    x1*(-0.83*c1 - 0.27*s1 - 0.5) + x2*(-0.27*c1 + 0.83*s1 - 1.06) - 0.3059;
    x1*(-0.34*c2 + 0.22*s2 - 0.5) + x2*(0.22*c2 + 0.34*s2 - 1.06) - 0.6048;
    x1*(-0.13*c3 + 0.43*s3 - 0.5) + x2*(0.43*c3 + 0.13*s3 - 1.06) - 0.5859;
    x1*(0.22*c4 + 0.78*s4 - 0.5) + x2*(0.78*c4 - 0.22*s4 - 1.06) - 0.3584;
    y1*(-1.83*c1 - 0.27*s1 + 0.5) + y2*(-0.27*c1 + 1.83*s1 - 1.06) + 1.0241;
    y1*(-1.34*c2 + 0.22*s2 + 0.5) + y2*(0.22*c2 + 1.34*s2 - 1.06) + 0.2352;
    y1*(-1.13*c3 + 0.43*s3 + 0.5) + y2*(0.43*c3 + 1.13*s3 - 1.06) + 0.0440999999999999;
    y1*(-0.78*c4 + 0.78*s4 + 0.5) + y2*(0.78*c4 + 0.78*s4 - 1.06) - 0.0784;

and then continues with `the solutions :` which is skipped
as the output of the function gives the list of real solutions.

::

    for (idx, sol) in enumerate(sols):
        (x1v, x2v) = (sol['x1'].real, sol['x2'].real)
        (y1v, y2v) = (sol['y1'].real, sol['y2'].real)
        print('Solution', idx+1, ':')
        print('x = ', x1v, x2v)
        print('y = ', y1v, y2v)

The coordinates of the real solutions are shown below.

::

    Solution 1 :
    x =  -0.0877960434509403 -0.85138690751564
    y =  0.235837391307301 -1.41899202703639
    Solution 2 :
    x =  0.0193359267851516 -0.937757011012446
    y =  1.22226669109342 -1.08285087742709
    Solution 3 :
    x =  -0.595728628822183 -0.617010917712341
    y =  0.118171353650905 -1.82939267557673
    Solution 4 :
    x =  -0.158077261086826 -0.793782551346416
    y =  -0.548761782690284 0.278116829722178
    Solution 5 :
    x =  14.265306631912 -6.51576530896231
    y =  -0.621791031677556 -0.0713939584963069
    Solution 6 :
    x =  -1.79178664902321 1.04613207405924
    y =  -1.46486338398045 1.21676347168425
    Solution 7 :
    x =  0.130643755560844 -0.942516053801942
    y =  0.963729735050218 -1.01577587226827
    Solution 8 :
    x =  -0.358757861563373 -0.537230434093211
    y =  0.0870595124133798 1.5543474028655
    Solution 9 :
    x =  -11.0926159017278 0.450863935272926
    y =  -0.396207302280832 -1.04172821286545
    Solution 10 :
    x =  -0.154697709323186 -0.812626279169727
    y =  3.30145715645532 -2.31860323051595
    Solution 11 :
    x =  -0.0801573081756841 -0.855275240173407
    y =  -0.297321862562434 -2.18414388671793
    Solution 12 :
    x =  0.676178657404253 -0.613650952963839
    y =  0.356055523659319 0.310794500797803
    Solution 13 :
    x =  1.4739209688177 -1.71128474823024
    y =  -0.654679846479676 0.028907166911727
    Solution 14 :
    x =  -0.264640920049152 -0.69691152780256
    y =  0.370368746423895 -1.54221173415608
    Solution 15 :
    x =  -1.0856845753759 -0.352998488913482
    y =  0.319028475056347 0.687883260707162

a four-bar mechanism
--------------------

The code in the function below are applied to make the plots.

::

    def angle(csa, sna):
        """
        Given in csa and sna are the cosine and sine of an angle a,
        that is: csa = cos(a) and sna = sin(a).
        On return is the angle a, with the proper orientation.
        """
        from math import acos, pi
        agl = acos(csa)
        if sna >= 0:
            return agl
        else:
            dlt = pi - agl
            return pi + dlt

::
   
    def angles(soldic):
        """
        Given a solution dictionary, extracts the angles from
        the four cosines and sines of the angles.
        Returns None if the angles are not ordered increasingly.
        Otherwise, returns the sequence of ordered angles.
        """
        from math import acos, asin
        c1v, s1v = soldic['c1'].real, soldic['s1'].real
        c2v, s2v = soldic['c2'].real, soldic['s2'].real
        c3v, s3v = soldic['c3'].real, soldic['s3'].real
        c4v, s4v = soldic['c4'].real, soldic['s4'].real
        ag1 = angle(c1v, s1v)
        ag2 = angle(c2v, s2v)
        ag3 = angle(c3v, s3v)
        ag4 = angle(c4v, s4v)
        ordered = (ag1 > ag2) and (ag2 > ag3) and (ag3 > ag4)
        if ordered:
            print(ag1, ag2, ag3, ag4, 'ordered angles')
            return (ag1, ag2, ag3, ag4)
        return None

::

    def plotpoints(points):
        """
        Plots the precision points and the pivots.
        """
        xpt = [a for (a, b) in points]
        ypt = [b for (a, b) in points]
        plt.plot(xpt, ypt, 'ro')
        plt.text(xpt[0] - 0.01, ypt[0] + 0.08, \"0\")
        plt.text(xpt[1] - 0.01, ypt[1] + 0.08, \"1\")
        plt.text(xpt[2] - 0.01, ypt[2] + 0.08, \"2\")
        plt.text(xpt[3] - 0.01, ypt[3] + 0.08, \"3\")
        plt.text(xpt[4] - 0.01, ypt[4] + 0.08, \"4\")
        plt.plot([0, 1], [0, 0], 'w^') # pivots marked by white triangles
        plt.axis([-1.0, 1.5, -1.0, 1.5])

::

    def plotbar(fig, points, idx, x, y):
        """
        Plots a 4-bar with coordinates given in x and y,
        and the five precision points in the list points.
        The index idx is the position with respect to a point in points.
        """
        if idx < 0:
            fig.add_subplot(231, aspect='equal')
        if idx == 0:
            fig.add_subplot(232, aspect='equal')
        elif idx == 1:
            fig.add_subplot(233, aspect='equal')
        elif idx == 2:
            fig.add_subplot(234, aspect='equal')
        elif idx == 3:
            fig.add_subplot(235, aspect='equal')
        elif idx == 4:
            fig.add_subplot(236, aspect='equal')
        plotpoints(points)
        if idx >= 0:
            xpt = [a for (a, b) in points]
            ypt = [b for (a, b) in points]
            (xp0, xp1) = (x[0] + xpt[idx], x[1] + ypt[idx])
            (yp0, yp1) = (y[0] + xpt[idx], y[1] + ypt[idx])
            plt.plot([xp0, yp0], [xp1, yp1], 'go')
            plt.plot([xp0, yp0], [xp1, yp1], 'g')
            plt.text(xp0 - 0.04, xp1 - 0.22, \"x\")
            plt.text(yp0 - 0.04, yp1 - 0.22, \"y\")
            plt.plot([0, xp0], [0, xp1], 'g')
            plt.plot([yp0, 1], [yp1, 0], 'g')
            plt.plot([xp0, xpt[idx]], [xp1, ypt[idx]], 'b')
            plt.plot([yp0, xpt[idx]], [yp1, ypt[idx]], 'b')

::

    def rotate(x, y, a):
        """
        Applies a planar rotation defined by the angle a
        to the points x and y.
        """
        from sympy.matrices import Matrix
        from math import cos, sin
        rot = Matrix([[cos(a), -sin(a)], [sin(a), cos(a)]])
        xmt = Matrix([[x[0]], [x[1]]])
        ymt = Matrix([[y[0]], [y[1]]])
        rxm = rot*xmt
        rym = rot*ymt
        rox = (rxm[0], rxm[1])
        roy = (rym[0], rym[1])
        return (rox, roy)

::

    def show4bar():
        """
        Plots a 4-bar design, for the five precision points
        on a straight line, with coordinates taken from Problem 7
        of the Morgan-Wampler paper.
        """
        pt0 = ( 0.50,  1.06)
        pt1 = (-0.83, -0.27)
        pt2 = (-0.34,  0.22)
        pt3 = (-0.13,  0.43)
        pt4 = ( 0.22,  0.78)
        points = [pt0, pt1, pt2, pt3, pt4]
        ags = [1.44734213756, 0.928413708131, 0.751699211109, 0.387116282208]
        x =  (-0.0877960434509, -0.851386907516)
        y =  (0.235837391307, -1.41899202704)
        fig = plt.figure()
        plotbar(fig,points, -1, x, y)
        plotbar(fig,points, 0, x, y)
        rx1, ry1 = rotate(x, y, ags[0])
        plotbar(fig,points, 1, rx1, ry1)
        rx2, ry2 = rotate(x, y, ags[1])
        plotbar(fig,points, 2, rx2, ry2)
        rx3, ry3 = rotate(x, y, ags[2])
        plotbar(fig,points, 3, rx3, ry3)
        rx4, ry4 = rotate(x, y, ags[3])
        plotbar(fig,points, 4, rx4, ry4)
        fig.canvas.draw()
        plt.savefig('fourbarfig1')

The mechanism which passes through the precision points is shown in
:numref:`fourbarfig1` obtained as the output of

::

    show4bar()

.. _fourbarfig1:

.. figure:: ./fourbarfig1.png
   :align: center
    
   A mechanism passing through precision points.

::

    for sol in sols:
        agl = angles(sol)
        if agl != None:
            (x1v, x2v) = (sol['x1'].real, sol['x2'].real)
            (y1v, y2v) = (sol['y1'].real, sol['y2'].real)
            print('x = ', x1v, x2v)
            print('y = ', y1v, y2v)

The output is

::

    1.4473421375642717 0.9284137081314461 0.75169921110931 0.3871162822082786 ordered angles
    x =  -0.0877960434509403 -0.85138690751564
    y =  0.235837391307301 -1.41899202703639
    2.524711332238134 0.9038272905536054 0.7498546795650226 0.38277375732994035 ordered angles
    x =  -0.0801573081756841 -0.855275240173407
    y =  -0.297321862562434 -2.18414388671793
    5.771983513802544 3.9629563185486125 3.442223836627024 0.5242754656511442 ordered angles
    x =  0.676178657404253 -0.613650952963839
    y =  0.356055523659319 0.310794500797803

Observe that one of the lists of ordered angles is used in the `showbar()`.

the coupler curve
-----------------

The coupler curve is the curve drawn by the coupler point.

::
   
    def plotpoints2(points):
        """
        Plots the precision points and the pivots.
        """
        xpt = [a for (a, b) in points]
        ypt = [b for (a, b) in points]
        plt.plot(xpt, ypt, 'ro')
        plt.text(xpt[0] + 0.01, ypt[0] + 0.06, \"0\")
        plt.text(xpt[1] - 0.03, ypt[1] + 0.06, \"1\")
        plt.text(xpt[2] - 0.01, ypt[2] + 0.06, \"2\")
        plt.text(xpt[3] - 0.01, ypt[3] + 0.06, \"3\")
        plt.text(xpt[4] - 0.01, ypt[4] + 0.06, \"4\")
        plt.plot([0, 1], [0, 0], 'w^') # pivots marked by white triangles
        plt.axis([-1.2, 1.2, -1.0, 1.5])

::

    def plotbar2(fig, points, idx, x, y):
        """
        Plots a 4-bar with coordinates given in x and y,
        and the five precision points in the list points.
        The index idx is the position with respect to a point in points.
        """
        plotpoints2(points)
        xpt = [a for (a, b) in points]
        ypt = [b for (a, b) in points]
        (xp0, xp1) = (x[0] + xpt[0], x[1] + ypt[0])
        (yp0, yp1) = (y[0] + xpt[0], y[1] + ypt[0])
        if idx >= 0:
            (xp0, xp1) = (x[0] + xpt[idx], x[1] + ypt[idx])
            (yp0, yp1) = (y[0] + xpt[idx], y[1] + ypt[idx])
            plt.plot([xp0, yp0], [xp1, yp1], 'go')
            plt.plot([xp0, yp0], [xp1, yp1], 'g')
            plt.text(xp0 - 0.04, xp1 - 0.12, \"x\")
            plt.text(yp0 - 0.04, yp1 - 0.12, \"y\")
            plt.plot([0, xp0], [0, xp1], 'g')
            plt.plot([yp0, 1], [yp1, 0], 'g')
            plt.plot([xp0, xpt[idx]], [xp1, ypt[idx]], 'b')
            plt.plot([yp0, xpt[idx]], [yp1, ypt[idx]], 'b')

::

    def lenbar(pt0, x, y):
        """
        In pt0 are the coordinates of the first precision point
        and in x and y the coordinates of the solution design.
        Returns the length of the bar between x and y.
        """
        (xp0, xp1) = (x[0] + pt0[0], x[1] + pt0[1])
        (yp0, yp1) = (y[0] + pt0[0], y[1] + pt0[1])
        result = sqrt((xp0 - yp0)**2 + (xp1 - yp1)**2)
        return result

::

    def coupler(x, y, xr, yr):
        """
        In x and y are the coordinates of the solution design.
        In xr and yr are the distances to the coupler point.
        Computes the intersection between two circles, centered
        at x and y, with respective radii in xr and yr.
        """
        A = -2*x[0] + 2*y[0]
        B = -2*x[1] + 2*y[1]
        C = x[0]**2 + x[1]**2 - xr**2 - y[0]**2 - y[1]**2 + yr**2
        fail = True
        if A + 1.0 != 1.0: # eliminate z1
            (alpha, beta) = (-C/A, -B/A)
            a = beta**2 + 1
            b = 2*alpha*beta - 2*x[1] - 2*x[0]*beta
            c = alpha**2 + x[0]**2 + x[1]**2 - xr**2 - 2*x[0]*alpha
            if b**2 - 4*a*c >= 0:
                fail = False
                disc = sqrt(b**2 - 4*a*c)
                z2 = (-b + disc)/(2*a)
                z1 = alpha + beta*z2
        if fail:
            (alpha, beta) = (-C/B, -A/B)
            a = beta**2 + 1
            b = 2*alpha*beta - 2*y[1] - 2*y[0]*beta
            c = alpha**2 + y[0]**2 + y[1]**2 - yr**2 - 2*y[0]*alpha
            disc = sqrt(b**2 - 4*a*c)
            z1 = (-b + disc)/(2*a)
            z2 = alpha + beta*z1
            dxz = sqrt((x[0]-z1)**2 + (x[1]-z2)**2)
        return (z1, z2)

::

    def xcrank(pt0, x):
        """
        In pt0 are the coordinates of the first precision point
        and in x the coordinates of the solution design.
        This function computes the length of the crank
        and its initial angle with respect to the first point.
        """
        from math import atan
        (xp0, xp1) = (x[0] + pt0[0], x[1] + pt0[1])
        crklen = sqrt(xp0**2 + xp1**2)
        crkagl = atan(xp1/xp0)
        return (crklen, crkagl)

::

    def ycrank(pt0, y):
        """
        In pt0 are the coordinates of the first precision point
        and in y the coordinates of the solution design.
        This function computes the length of the crank
        and its initial angle with respect to the first point.
        """
        from math import cos, sin, acos, pi
        (yp0, yp1) = (y[0] + pt0[0], y[1] + pt0[1])
        crklen = sqrt((yp0 - 1)**2 + yp1**2)
        crkagl = acos((yp0-1)/crklen)
        if yp1 < 0:
            dlt = pi - crkagl
            crkagl = pi + dlt
        cx = 1 + crklen*cos(crkagl)
        cy = crklen*sin(crkagl)
        return (crklen, crkagl)

::

    def xpos(y1, y2, dxy, rad):
        """
        Given in y1 and y2 are the coordinates of the point y,
        in dxy is the distance between the points x and y,
        and rad is the distance between x and (1, 0).
        The coordinates of the point x are returned in a tuple.
        """
        A = -2*y1  # coefficient with y1
        B = -2*y2  # coefficient with y2
        C = y1**2 + y2**2 - dxy**2 + rad**2 # constant
        fail = True
        if abs(y2) < 1.0e-8:
            x1 = -C/A
            x2sqr = rad**2 - x1**2
            x2 = sqrt(x2sqr)
            fail = False
        else: # eliminate x2
            (alpha, beta) = (-C/B, -A/B)
            (a, b, c) = (1+beta**2, 2*alpha*beta, alpha**2 - rad**2)
            b4ac = b**2 - 4*a*c
            disc = sqrt(b4ac)
            x1m = (-b - disc)/(2*a)
            x2m = alpha + beta*x1m
            x1p = (-b + disc)/(2*a)
            x2p = alpha + beta*x1p
        return ((x1m, x2m), (x1p, x2p))

::

    def plotcrank(crk, agl, dxy, rad, xrd, yrd):
        """
        Plots several positions of the crank.  On input are:
        crk : length of the crank from the point y to (1, 0),
        agl : start angle,
        rad : length of the crank from (0, 0) to the point x,
        xrd : length from the point x to the coupler point,
        yrd : length from the point y to the coupler point.
        """
        from math import sin, cos, pi
        (xzm, yzm) = ([], [])
        (xzp, yzp) = ([], [])
        nbr = 205
        inc = (pi+0.11763)/nbr
        b = agl - 2.558 # 125
        for k in range(nbr):
            (y1, y2) = (1 + crk*cos(b), crk*sin(b))
            (xm, xp) = xpos(y1, y2, dxy, rad)
            (x1m, x2m) = xm
            (x1p, x2p) = xp
            (z1m, z2m) = coupler([x1m, x2m], [y1, y2], xrd, yrd)
            (z1p, z2p) = coupler([x1p, x2p], [y1, y2], xrd, yrd)
            xzm.append(z1m)
            yzm.append(z2m)
            xzp.append(z1p)
            yzp.append(z2p)
            if k < 0: # selective plot
                plt.plot([0, x1m], [0, x2m], 'g')
                plt.plot([x1m, y1], [x2m, y2], 'g')
                plt.plot([y1, 1], [y2, 0], 'g')
                dyp = sqrt((y1-1)**2 + y2**2)
                dyx = sqrt((x1m-y1)**2 + (x2m-y2)**2)
                print('dxy =', dxy, 'dyp =', dyp)
            if k < 0:
                print('y2 =', y2)\
                plt.plot([x1m, z1m], [x2m, z2m], 'b')
                plt.plot([y1, z1m], [y2, z2m], 'b')
                plt.plot([x1p, z1p], [x2p, z2p], 'b')
                plt.plot([y1, z1p], [y2, z2p], 'b')
            b = b + inc
        plt.plot(xzp[:1]+xzm[:102]+xzp[102:], \
                 yzp[:1]+yzm[:102]+yzp[102:], 'r')
        plt.plot(xzp[:102]+xzm[102:], yzp[:102]+yzm[102:], 'r')

::

    def plotcoupler():
        """
        Plots the coupler curve for a straight line 4-bar mechanism.
        """
        pt0 = ( 0.50,  1.06)
        pt1 = (-0.83, -0.27)
        pt2 = (-0.34,  0.22)
        pt3 = (-0.13,  0.43)
        pt4 = ( 0.22,  0.78)
        points = [pt0, pt1, pt2, pt3, pt4]
        ags = [1.44734213756, 0.928413708131, 0.751699211109, 0.387116282208]
        x = (-0.0877960434509, -0.851386907516)
        y = (0.235837391307, -1.41899202704)
        (xcrk, xagl) = xcrank(pt0, x)
        (ycrk, yagl) = ycrank(pt0, y)
        dxy = lenbar(pt0, x, y)
        fig = plt.figure()
        fig.add_subplot(111, aspect='equal')
        xrd = sqrt(x[0]**2 + x[1]**2) # distance from x to pt0
        yrd = sqrt(y[0]**2 + y[1]**2) # distance from y to pt0
        plotcrank(ycrk, yagl, dxy, xcrk, xrd, yrd)
        plotbar2(fig, points, 0, x, y)
        fig.canvas.draw()
        plt.savefig('fourbarfig2')

Running the function

:: 

    plotcoupler()

produces the plot in :numref:`fourbarfig2`.

.. _fourbarfig2:

.. figure:: ./fourbarfig2.png
   :align: center
    
   The coupler curve of a 4-bar mechanism.

Two Lines Meeting Four Given Lines
==================================

Given four lines in general position,
there are two lines which meet all four given lines.
With Pieri homotopies we can solve this Schubert problem.
For the verification of the intersection conditions, `numpy` is used.
The plots are made with `matplotlib`.

We use random numbers and for reproducible plots, fix the seed.

::

   from random import seed

From `numpy` we import the following.

::

   from numpy import zeros, array, concatenate, matrix
   from numpy.linalg import det, solve

The plots are in 3-space.

::

   import matplotlib.pyplot as plt
   from mpl_toolkits.mplot3d import Axes3D

From `phcpy` we import the following functions:

::

    "from phcpy.solutions import coordinates
    "from phcpy.schubert import random_complex_matrix
    "from phcpy.schubert import pieri_root_count, run_pieri_homotopies
    "from phcpy.schubert import real_osculating_planes
    "from phcpy.schubert import make_pieri_system
    "from phcpy.trackers import double_track as track

solving a general instance
--------------------------

A random instance of the four given lines will lead to two solution lines.
The formal root count run as

::

    (mdim, pdim, deg) = (2, 2, 0)
    pcnt = pieri_root_count(mdim, pdim, deg, False)
    pcnt

and outputs `2`.

To setup the problem, some auxiliary functions are first defined.

::

    def indices(name):
        """
        For the string name in the format xij
        return (i, j) as two integer indices.
        """
        return (int(name[1]), int(name[2]))

::

    def solution_plane(rows, cols, sol):
        """
        Returns a sympy matrix with as many rows
        as the value of rows and with as many columns
        as the value of columns, using the string
        represention of a solution in sol.
        """
        result = zeros((rows, cols), dtype=complex)
        for k in range(cols):
            result[k][k] = 1
        (vars, vals) = coordinates(sol)
        for (name, value) in zip(vars, vals):
            i, j = indices(name)
            result[i-1][j-1] = value
        return result

::

    def verify_determinants(inps, sols, verbose=True):
        """
        Verifies the intersection conditions with determinants,
        concatenating the planes in inps with those in the sols.
        Both inps and sols are lists of numpy arrays.
        Returns the sum of the absolute values of all determinants.
        If verbose, then for all solutions in sols, the computed
        determinants are printed to screen.
        """
        checksum = 0
        for sol in sols:
            if verbose:
                print('checking solution\\n', sol)
            for plane in inps:
                cat = concatenate([plane, sol], axis=-1)
                mat = matrix(cat)
                dcm = det(mat)
                if verbose:
                    print('the determinant :', dcm)
                checksum = checksum + abs(dcm)
        return checksum

::

    def solve_general(mdim, pdim, qdeg):
        """
        Solves a general instance of Pieri problem, computing the
        p-plane producing curves of degree qdeg which meet a number
        of general m-planes at general interpolation points,
        where p = pdim and m = mdim on input.
        For the problem of computing the two lines which meet
        four general lines, mdim = 2, pdim = 2, and qdeg = 0.
        Returns a tuple with four lists.
        The first two lists contain matrices with the input planes
        and the solution planes respectively.
        The third list is the list of polynomials solved
        and the last list is the solution list.
        """
        dim = mdim*pdim + qdeg*(mdim+pdim)
        ranplanes = [random_complex_matrix(mdim+pdim, mdim) for _ in range(0, dim)]
        (pols, sols) = run_pieri_homotopies(mdim, pdim, qdeg, ranplanes)
        inplanes = [array(plane) for plane in ranplanes]
        outplanes = [solution_plane(mdim+pdim, pdim, sol) for sol in sols]
        return (inplanes, outplanes, pols, sols)

::

    (inp, otp, pols, sols) = solve_general(mdim, pdim, deg)

The four input lines are represented as matrices.

::

    for plane in inp:
        print(plane)

shows 

::

    [[ 0.98771734-0.15625123j  0.52929265-0.84843933j]
     [ 0.0108879 -0.99994073j  0.43271012+0.90153311j]
     [ 0.670366  +0.74203061j  0.84995049-0.52686257j]
     [-0.99870177+0.05093886j  0.55311134-0.83310735j]]
    [[ 0.1176291 +0.9930576j   0.73982601-0.67279824j]
     [-0.4096813 -0.91222872j  0.98222659+0.18769903j]
     [ 0.49367521+0.86964635j -0.00101345-0.99999949j]
     [ 0.99603164-0.0889999j   0.37233497-0.92809841j]]
    [[-0.86632581+0.49947932j  0.99954174-0.03027052j]
     [ 0.26897023+0.96314849j  0.29943145+0.95411781j]
     [ 0.77919846-0.62677728j  0.52235751-0.85272659j]
     [ 0.4481898 +0.89393842j  0.97691942+0.21360816j]]
    [[ 0.40705515-0.91340358j -0.66900116+0.74326136j]
     [-0.11164153+0.99374854j -0.51718407-0.8558742j ]
     [-0.01384859+0.9999041j  -0.38779064+0.92174748j]
     [ 0.32407475-0.94603148j  0.87995025-0.47506584j]]

::

    print('The solution planes :')
    for plane in otp:
        print(plane)

has as output

::

    The solution planes :
    [[ 1.        +0.j          0.        +0.j        ]
     [-0.64379718+0.67758706j  1.        +0.j        ]
     [ 0.69735824-0.15805905j -1.46030164-0.68747669j]
     [ 0.        +0.j         -1.74595349+0.00175246j]]
    [[ 1.        +0.j          0.        +0.j        ]
     [ 1.4746012 +0.78327696j  1.        +0.j        ]
     [ 1.20071164-2.11957742j  0.91569812-1.31875637j]
     [ 0.        +0.j         -1.04202682+0.09584754j]]

To check the solutions, we use `numpy` as follows:

::

    check = verify_determinants(inp, otp)
    print('Sum of absolute values of determinants :', check)

The output of the check is

::

    checking solution
    [[ 1.        +0.j          0.        +0.j        ]
     [-0.64379718+0.67758706j  1.        +0.j        ]
     [ 0.69735824-0.15805905j -1.46030164-0.68747669j]
     [ 0.        +0.j         -1.74595349+0.00175246j]]
    the determinant : (2.9667224835639593e-15+1.3550262739027277e-15j)
    the determinant : (4.195866422887001e-15-1.4293281742484199e-15j)
    the determinant : (-1.8017495082844853e-15-1.5770416093056093e-15j)
    the determinant : (-2.0927676352675787e-16+1.091663409852285e-15j)
    checking solution
    [[ 1.        +0.j          0.        +0.j        ]
     [ 1.4746012 +0.78327696j  1.        +0.j        ]
     [ 1.20071164-2.11957742j  0.91569812-1.31875637j]
     [ 0.        +0.j         -1.04202682+0.09584754j]]
    the determinant : (1.0002339027616943e-14-3.132413944024583e-14j)
    the determinant : (2.8791053191246284e-14-3.6564204184655514e-15j)
    the determinant : (-3.605052372912635e-14+5.874582883240587e-15j)
    the determinant : (-2.6498852748806624e-14-2.7706915851697867e-15j)
    Sum of absolute values of determinants : 1.362741358344356e-13

Observe that all determines evaluate to numbers close to machine precision.

four real lines
---------------

We can generate inputs for which all solutions are real.

::

    def solve_real(mdim, pdim, start, sols):
        """
        Solves a real instance of Pieri problem, for input planes
        of dimension mdim osculating a rational normal curve.
        On return are the planes of dimension pdim.
        """
        oscplanes = real_osculating_planes(mdim, pdim, 0)
        target = make_pieri_system(mdim, pdim, 0, oscplanes, is_real=True)
        gamma, rtsols = track(target, start, sols)
        print('The solutions to the real problem :')
        for (idx, sol) in enumerate(rtsols):
            print('Solution', idx+1, ':')
            print(sol)
        inplanes = [array(plane) for plane in oscplanes]
        outplanes = [solution_plane(mdim+pdim, pdim, sol) for sol in rtsols]
        return (inplanes, outplanes, target, rtsols)

For visualization, the seed of the random number generators is set fixed."

::

    seed(400)"

The output of

::

    (oscp, otp2, pols2, sols2) = solve_real(mdim, pdim, pols, sols)

is

::

    The solutions to the real problem :
    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x21 : -2.84638025557899E-02   1.31371731030452E-46
     x32 : -1.19348750548289E-01  -2.62743462060903E-46
     x42 : -4.99706612461873E+00   2.38220738935219E-44
     x31 : -1.06771882518925E+00   3.15292154473084E-45
    == err :  5.410E-15 = rco :  5.611E-03 = res :  5.551E-16 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x21 : -5.52734869685360E-02   5.47382212626882E-48
     x32 : -1.19348750548290E-01   4.37905770101505E-47
     x42 : -2.57330433323918E+00   3.83167548838817E-47
     x31 : -6.60558824288729E-01   1.91583774419409E-47,
    == err :  6.174E-16 = rco :  1.324E-02 = res :  3.747E-16 =

::

    print('The input planes :')
    for plane in oscp:
        print(plane)

::

    The input planes :
    [[-0.63223829 -0.07958136]
     [ 0.24317589 -0.42625018]
     [ 0.44517428  0.75891681]
     [-0.58562795  0.48582185]]
    [[-0.63156273  0.07848797]
     [ 0.31098671 -0.49445305]
     [ 0.32529788  0.85734178]
     [-0.63134544  0.11966993]]
    [[-0.66765465 -0.21150281]
     [-0.41782225 -0.46153796]
     [ 0.14470336 -0.77816698]
     [ 0.59893469 -0.36973696]]
    [[-0.69033039  0.11246161]
     [-0.09104114 -0.32159814]
     [ 0.66631728 -0.28602371]
     [ 0.26678969  0.89561011]]

::

    print('The solution planes :')
    for plane in otp2:
        print(plane)

::

    The solution planes :
    [[ 1.        +0.00000000e+00j  0.        +0.00000000e+00j]
     [-0.0284638 +1.31371731e-46j  1.        +0.00000000e+00j]
     [-1.06771883+3.15292154e-45j -0.11934875-2.62743462e-46j]
     [ 0.        +0.00000000e+00j -4.99706612+2.38220739e-44j]]
    [[ 1.        +0.00000000e+00j  0.        +0.00000000e+00j]
     [-0.05527349+5.47382213e-48j  1.        +0.00000000e+00j]
     [-0.66055882+1.91583774e-47j -0.11934875+4.37905770e-47j]
     [ 0.        +0.00000000e+00j -2.57330433+3.83167549e-47j]]


Let us verify the real solution planes as well:

::

    check = verify_determinants(oscp, otp2)
    print('Sum of absolute values of determinants :', check)

Observe the output of the verification:

::

    checking solution
    [[ 1.        +0.00000000e+00j  0.        +0.00000000e+00j]
     [-0.0284638 +1.31371731e-46j  1.        +0.00000000e+00j]
     [-1.06771883+3.15292154e-45j -0.11934875-2.62743462e-46j]
     [ 0.        +0.00000000e+00j -4.99706612+2.38220739e-44j]]
    the determinant : (2.7334976213462325e-15-2.490244814186718e-45j)
    the determinant : (6.194410394095717e-15-2.4210378066256254e-45j)
    the determinant : (6.1256567148522274e-15-4.974841059851325e-47j)
    the determinant : (-1.7538510134158814e-15-1.6274706457865366e-45j)
    checking solution
    [[ 1.        +0.00000000e+00j  0.        +0.00000000e+00j]
     [-0.05527349+5.47382213e-48j  1.        +0.00000000e+00j]
     [-0.66055882+1.91583774e-47j -0.11934875+4.37905770e-47j]
     [ 0.        +0.00000000e+00j -2.57330433+3.83167549e-47j]]
    the determinant : (-6.163408511151722e-16-7.868415222942327e-49j)
    the determinant : (3.1253636658440115e-16-1.6674497687062525e-48j)
    the determinant : (-1.4639612348256832e-16-1.964057003033534e-47j)
    the determinant : (-1.3795665037633665e-15+8.091364203252659e-48j)
    Sum of absolute values of determinants : 1.926225558865557e-14

Observe the size of the values of the determinants.

visualization
-------------

The code in the functions below help visualizing the problem.

::

    def input_generators(plane):
        """
        Given in plane is a numpy matrix, with in its columns
        the coordinates of the points which span a line, in 4-space.
        The first coordinate must not be zero.
        Returns the affine representation of the line,
        after dividing each generator by its first coordinate.
        """
        pone = list(plane[:,0])
        ptwo = list(plane[:,1])
        aone = [x/pone[0] for x in pone]
        atwo = [x/ptwo[0] for x in ptwo]
        return (aone[1:], atwo[1:])

::

    def output_generators(plane):
        """
        Given in plane is a numpy matrix, with in its columns
        the coordinates of the points which span a line, in 4-space.
        The solution planes follow the localization pattern
        1, *, *, 0 for the first point and 0, 1, *, * for
        the second point, which means that the second point
        in standard projective coordinates lies at infinity.
        For the second generator, the sum of the points is taken.
        The imaginary part of each coordinate is omitted.
        """
        pone = list(plane[:,0])
        ptwo = list(plane[:,1])
        aone = [x.real for x in pone]
        atwo = [x.real + y.real for (x, y) in zip(pone, ptwo)]
        return (aone[1:], atwo[1:])

::

    def boxrange(inlines, outlines):
        """
        Returns a list of three lists with the [min, max]
        values of each coordinate of each generator in the lists
        inlines and outlines.
        The ranges are adjusted for the particular real case.
        """
        fst = inlines[0][0]
        result = {'xmin': fst[0], 'xmax': fst[0], \
                  'ymin': fst[1], 'ymax': fst[1], \
                  'zmin': fst[2], 'zmax': fst[2]} 
        pts = [x for (x, y) in inlines] + [y for (x, y) in inlines] \
            + [x for (x, y) in outlines] + [y for (x, y) in outlines]
        print('the points :\n', pts)
        for point in pts:
            result['xmin'] = min(result['xmin'], point[0])
            result['ymin'] = min(result['ymin'], point[1])
            result['zmin'] = min(result['zmin'], point[2])
            result['xmax'] = max(result['xmax'], point[0])
            result['ymax'] = max(result['ymax'], point[1])
            result['zmax'] = max(result['zmax'], point[2])
        return ((result['xmin']+3, result['xmax']-3), \
                (result['ymin']+8, result['ymax']-11), \
                (result['zmin']+3, result['zmax']-5))

::

    def inbox(point, lims):
        """
        Returns true if the coordinates of the point
        are in the box defined by the 3-tuple lims
        which contain the minima and maxima for the coordinates.
        """
        tol = 1.0e-8 # this is essential for roundoff
        (xlim, ylim, zlim) = lims
        if point[0] < xlim[0] - tol:
            return False
        elif point[0] > xlim[1] + tol:
            return False
        elif point[1] < ylim[0] - tol:
            return False
        elif point[1] > ylim[1] + tol:
            return False
        elif point[2] < zlim[0] - tol:
            return False
        elif point[2] > zlim[1] + tol:
            return False
        else:
            return True

::

    def equal(pt1, pt2):
        """
        Returns true if the all coordinates of pt1 and pt2
        match up to a tolerance of 1.0e-10.
        """
        tol = 1.0e-8
        if abs(pt1[0] - pt2[0]) > tol:
            return False
        elif abs(pt1[1] - pt2[1]) > tol:
            return False
        elif abs(pt1[2] - pt2[2]) > tol:
            return False
        return True

::

    def isin(points, pnt):
        """
        Returns true if pnt belongs to the list points.
        """
        if len(points) == 0:
            return False
        else:
            for point in points:
                if equal(point, pnt):
                    return True
            return False

::

    def plot_line(axs, line, lims, color):
        """
        Plots the line defined as a tuple of two points,
        using the axis object in axs.
        The 3-tuple lims contains three lists with limits [min, max]
        for the x, y, and z coordinates.
        """
        (fst, snd) = line
        axs.set_xlabel('x')
        axs.set_ylabel('y')
        axs.set_zlabel('z')
        axs.set_xlim(lims[0])
        axs.set_ylim(lims[1])
        axs.set_zlim(lims[2])
        dir = (fst[0] - snd[0], fst[1] - snd[1], fst[2] - snd[2])
        result = []
        for k in range(3):
            fac = (lims[k][1]-fst[k])/dir[k]
            pnt = (fst[0] + fac*dir[0], fst[1] + fac*dir[1], fst[2] + fac*dir[2])
            if inbox(pnt, lims):
                if not isin(result, pnt): result.append(pnt)
        for k in range(3):
            fac = (lims[k][0]-fst[k])/dir[k]
            pnt = (fst[0] + fac*dir[0], fst[1] + fac*dir[1], fst[2] + fac*dir[2])
            if inbox(pnt, lims):
                if not isin(result, pnt): result.append(pnt)
        (one, two) = (result[0], result[1])
        # axs.plot([fst[0], snd[0]], [fst[1], snd[1]], [fst[2], snd[2]], 'bo')
        # axs.plot([one[0], two[0]], [one[1], two[1]], [one[2], two[2]], 'ro')
        axs.plot([one[0], two[0]], [one[1], two[1]], [one[2], two[2]], color)
        plt.savefig('fourlinesfig1')

::

    def plot_lines(inlines, outlines, points, lims):
        """
        Generates coordinates of the points in a random line
        and then plots this line.  The intersection points are
        in the list points and limits for the bounding box in lims
        """
        fig = plt.figure()
        axs = fig.add_subplot(111, projection='3d')
        for line in inlines:
            plot_line(axs, line, lims, 'b')
        for line in outlines:
            plot_line(axs, line, lims, 'r')
        for point in points:
            axs.plot([point[0]], [point[1]], [point[2]], 'ro')
        axs.view_init(azim=5, elev=20)
        plt.show()
        plt.savefig('fourlinesfig2')

::

    def intersection_point(apl, bpl, check=True):
        """
        Given in apl the two points that define a line
        and in bpl the two points that define another line,
        returns the intersection point.
        If check, then additional tests are done
        and the outcome of the tests is written to screen.
        """
        (apt, bpt) = apl
        (cpt, dpt) = bpl
        mat = array([[apt[0], bpt[0], -cpt[0]], \
                     [apt[1], bpt[1], -cpt[1]], \
                     [apt[2], bpt[2], -cpt[2]]])
        rhs = array([[dpt[0]], [dpt[1]], [dpt[2]]])
        sol = solve(mat, rhs)
        cff = list(sol[:,0])
        csm = cff[0] + cff[1]
        result = ((cff[0]*apt[0] + cff[1]*bpt[0])/csm, \
                  (cff[0]*apt[1] + cff[1]*bpt[1])/csm, \
                  (cff[0]*apt[2] + cff[1]*bpt[2])/csm)
        if check:
            csm = cff[2] + 1.0
            verify = ((cff[2]*cpt[0] + dpt[0])/csm, \
                      (cff[2]*cpt[1] + dpt[1])/csm, \
                      (cff[2]*cpt[2] + dpt[2])/csm)
            print('the solution :\\n', result)
            print('the solution verified :\\n', verify)
            res = matrix(rhs) - matrix(mat)*matrix(sol)
            print('the residual :\n', res)
        return result

::

    def intersection_points(ipl, opl):
        """
        Returns the list of intersection points between
        the input planes in ipl and the output planes in opl.
        """
        result = []
        for inplane in ipl:
            for outplane in opl:
                result.append(intersection_point(inplane, outplane))
        return result

::

    def show_planes(ipl, opl):
        """
        Shows the input and the output planes.
        """
        (inlines, outlines) = ([], [])
        for plane in ipl:
            inlines.append(input_generators(plane))
        for plane in opl:
            outlines.append(output_generators(plane))
        print('The generators of the input lines :')
        for line in inlines:
            print(line)
        print('The generators of the output lines :')
        for line in outlines:
            print(line)
        brg = boxrange(inlines, outlines)
        print('the range:', brg)
        intpts = intersection_points(inlines, outlines)
        print('the intersection points :')
        for point in intpts:
            print(point)
        plot_lines(inlines, outlines, intpts, brg)
        plt.savefig('fourlinesfig3')

We end up with an interactive backend for the 3d plot.

::

    %matplotlib widget
    show_planes(oscp, otp2)

produces the following output:

::

    The generators of the input lines :
    ([-0.3846269613221122, -0.7041242012482366, 0.9262772651610497], [5.356155982058531, -9.53636379773747, -6.104719131981401])
    ([-0.4924082638753003, -0.5150682033346254, 0.9996559434380463], [-6.2997304815421655, 10.923225600528575, 1.5246914154709972])
    ([0.6258059455871108, -0.2167338369356443, -0.8970725931155779], [2.1821835337633937, 3.6792275561013374, 1.7481420529767318])
    ([0.13188052906200864, -0.9652150521086493, -0.3864666725234145], [-2.8596260453517486, -2.5433009310133032, 7.963696648872165])
    The generators of the output lines :
    ([-0.0284638025557899, -1.06771882518925, 0.0], [0.97153619744421, -1.187067575737539, -4.99706612461873])
    ([-0.055273486968536, -0.660558824288729, 0.0], [0.944726513031464, -0.779907574837019, -2.57330433323918])
    the points :
    [[-0.3846269613221122, -0.7041242012482366, 0.9262772651610497], [-0.4924082638753003, -0.5150682033346254, 0.9996559434380463], [0.6258059455871108, -0.2167338369356443, -0.8970725931155779], [0.13188052906200864, -0.9652150521086493, -0.3864666725234145], [5.356155982058531, -9.53636379773747, -6.104719131981401], [-6.2997304815421655, 10.923225600528575, 1.5246914154709972], [2.1821835337633937, 3.6792275561013374, 1.7481420529767318], [-2.8596260453517486, -2.5433009310133032, 7.963696648872165], [-0.0284638025557899, -1.06771882518925, 0.0], [-0.055273486968536, -0.660558824288729, 0.0], [0.97153619744421, -1.187067575737539, -4.99706612461873], [0.944726513031464, -0.779907574837019, -2.57330433323918]]
    the range: ((-3.2997304815421655, 2.3561559820585307), (-1.5363637977374704, -0.07677439947142517), (-3.104719131981401, 2.9636966488721654))
    the solution :
     (-0.15837537533646365, -1.052214041296111, 0.6491767195382462)
    the solution verified :
     (-0.15837537533646406, -1.0522140412961136, 0.6491767195382475)
    the residual :
     [[4.44089210e-16]
     [1.11022302e-15]
     [0.00000000e+00]]
    the solution :
     (-0.4430230234123302, -0.6142814015884848, 0.9977975623422988)
    the solution verified :
     (-0.44302302341232946, -0.6142814015884835, 0.997797562342297)
    the residual :
     [[ 0.00000000e+00]
     [-2.22044605e-16]
     [ 0.00000000e+00]]
    the solution :
     (-0.2236498742909531, -1.0444236114032255, 0.9753577070651858)
    the solution verified :
     (-0.22364987429095395, -1.0444236114032293, 0.9753577070651895)
    the residual :
     [[-1.11022302e-16]
     [-6.66133815e-16]
     [ 0.00000000e+00]]
    the solution :
     (-0.441973240857878, -0.6144066918247044, 0.9950961523459683)
    the solution verified :
     (-0.44197324085787826, -0.6144066918247048, 0.9950961523459688)
    the residual :
     [[1.11022302e-16]
     [2.22044605e-16]
     [0.00000000e+00]]
    the solution :
     (0.2715464337673154, -1.1035246720461052, -1.4991709889690488)
    the solution verified :
     (0.27154643376731663, -1.1035246720461096, -1.4991709889690552)
    the residual :
     [[1.11022302e-16]
     [2.22044605e-16]
     [0.00000000e+00]]
    the solution :
     (0.42557851238329614, -0.7179479096100174, -1.2373785335787928)
    the solution verified :
     (0.42557851238329597, -0.7179479096100173, -1.2373785335787926)
    the residual :
     [[-1.11022302e-16]
     [ 0.00000000e+00]
     [ 0.00000000e+00]]
    the solution :
     (-0.056164218290926694, -1.0644128151815966, 0.1384208091079073)
    the solution verified :
     (-0.05616421829092654, -1.0644128151815933, 0.1384208091079069)
    the residual :
     [[ 6.66133815e-16]
     [-2.44249065e-15]
     [-1.77635684e-15]]
    the solution :
     (0.5683194604437922, -0.7349838634131002, -1.6046944337535327)
    the solution verified :
     (0.5683194604438059, -0.7349838634131174, -1.6046944337535711)
    the residual :
     [[ 1.11022302e-16]
     [ 3.33066907e-16]
     [-4.44089210e-16]]
    the intersection points :
    (-0.15837537533646365, -1.052214041296111, 0.6491767195382462)
    (-0.4430230234123302, -0.6142814015884848, 0.9977975623422988)
    (-0.2236498742909531, -1.0444236114032255, 0.9753577070651858)
    (-0.441973240857878, -0.6144066918247044, 0.9950961523459683)
    (0.2715464337673154, -1.1035246720461052, -1.4991709889690488)
    (0.42557851238329614, -0.7179479096100174, -1.2373785335787928)
    (-0.056164218290926694, -1.0644128151815966, 0.1384208091079073)
    (0.5683194604437922, -0.7349838634131002, -1.6046944337535327)

The code produces :numref:`fourlinesfig3`.

.. _fourlinesfig3:

.. figure:: ./fourlinesfig3.png
   :align: center
    
   Two lines meeting four given lines.


The Circle Problem of Apollonius
================================

The circle problem of Apollonius has the following input/output specification:

Given three circles, find all circles that are tangent to the given circles.

the polynomial systems
----------------------

Without loss of generality, we take the first circle to be the unit circle,
centered at (0, 0) and with radius 1.  The origin of the second circle lies
on the first coordinate axis, so its center has coordinates (`c2x`, 0) and
radius `r2`.  The third circle has center (`c3x`, `c3y`) and radius `r3`.
So there are five parameters in this problem: `c2x`, `r2`, `c3x`, `c3y`,
and `r3`.
Values for the five parameters are defined by the first five equations.
The next three equations determine the center (`x`, `y`) and the radius `r`
of the circle which touches the three given circles.
The condition on the center of the touching circle is that its distance
to the center of the given circle is either the difference or the sum of
the radii of both circles.  So we arrive at eight polynomial systems.

The problem formulation is coded in the function `polynomials`.

::

    def polynomials(c2x, r2, c3x, c3y, r3):
        """
        On input are the five parameters of the circle problem of Apollonius:
        c2x : the x-coordinate of the center of the second circle,
        r2 : the radius of the second circle,
        c3x : the x-coordinate of the center of the third circle,
        c3y : the y-coordinate of the center of the third circle,
        r3 : the radius of the third circle.
        Returns a list of lists.  Each list contains a polynomial system.
        Solutions to each polynomial system define center (x, y) and radius r
        of a circle touching three given circles.
        """
        e1m = 'x^2 + y^2 - (r-1)^2;'
        e1p = 'x^2 + y^2 - (r+1)^2;'
        e2m = '(x-%.15f)^2 + y^2 - (r-%.15f)^2;' % (c2x, r2)
        e2p = '(x-%.15f)^2 + y^2 - (r+%.15f)^2;' % (c2x, r2)
        e3m = '(x-%.15f)^2 + (y-%.15f)^2 - (r-%.15f)^2;' % (c3x, c3y, r3)
        e3p = '(x-%.15f)^2 + (y-%.15f)^2 - (r+%.15f)^2;' % (c3x, c3y, r3)
        eqs0 = [e1m,e2m,e3m]
        eqs1 = [e1m,e2m,e3p]
        eqs2 = [e1m,e2p,e3m]
        eqs3 = [e1m,e2p,e3p]
        eqs4 = [e1p,e2m,e3m]
        eqs5 = [e1p,e2m,e3p]
        eqs6 = [e1p,e2p,e3m]
        eqs7 = [e1p,e2p,e3p]
        return [eqs0,eqs1,eqs2,eqs3,eqs4,eqs5,eqs6,eqs7]

As an example of a general problem, the center of the second circle 
is at `(2, 0)`, with radius `2/3`, and the third circle is centered
at `(1, 1)`, with a radius of `1/3`.

Let us look at the eight polynomial systems, 
computed as the output of the function `polynomials`.

::

    general_problem = polynomials(2, 2.0/3, 1, 1, 1.0/3)
    for pols in general_problem:
        print(pols)

The eight polynomial systems are shown below:

::

      ['x^2 + y^2 - (r-1)^2;', '(x-2.000000000000000)^2 + y^2 - (r-0.666666666666667)^2;', '(x-1.000000000000000)^2 + (y-1.000000000000000)^2 - (r-0.333333333333333)^2;']
      ['x^2 + y^2 - (r-1)^2;', '(x-2.000000000000000)^2 + y^2 - (r-0.666666666666667)^2;', '(x-1.000000000000000)^2 + (y-1.000000000000000)^2 - (r+0.333333333333333)^2;']
      ['x^2 + y^2 - (r-1)^2;', '(x-2.000000000000000)^2 + y^2 - (r+0.666666666666667)^2;', '(x-1.000000000000000)^2 + (y-1.000000000000000)^2 - (r-0.333333333333333)^2;']
      ['x^2 + y^2 - (r-1)^2;', '(x-2.000000000000000)^2 + y^2 - (r+0.666666666666667)^2;', '(x-1.000000000000000)^2 + (y-1.000000000000000)^2 - (r+0.333333333333333)^2;']
      ['x^2 + y^2 - (r+1)^2;', '(x-2.000000000000000)^2 + y^2 - (r-0.666666666666667)^2;', '(x-1.000000000000000)^2 + (y-1.000000000000000)^2 - (r-0.333333333333333)^2;']
      ['x^2 + y^2 - (r+1)^2;', '(x-2.000000000000000)^2 + y^2 - (r-0.666666666666667)^2;', '(x-1.000000000000000)^2 + (y-1.000000000000000)^2 - (r+0.333333333333333)^2;']
      ['x^2 + y^2 - (r+1)^2;', '(x-2.000000000000000)^2 + y^2 - (r+0.666666666666667)^2;', '(x-1.000000000000000)^2 + (y-1.000000000000000)^2 - (r-0.333333333333333)^2;']
      ['x^2 + y^2 - (r+1)^2;', '(x-2.000000000000000)^2 + y^2 - (r+0.666666666666667)^2;', '(x-1.000000000000000)^2 + (y-1.000000000000000)^2 - (r+0.333333333333333)^2;']

plotting circles
----------------

The package `matplotlib` has primitives to define circles.

::

    import matplotlib.pyplot as plt

The input to the three given circles of the general problem is codified 
in the list of tuples set below.

::

    crcdata = [((0, 0), 1), ((2, 0), 2.0/3), ((1, 1), 1.0/3)]

The input circles will be shown as blue disks.
Let us then render our general configuration.

::

    (xa, xb, ya, yb) = (-2, 4, -2, 3)

The code to make :numref:`apolloniusfig1` is below:

::

    fig = plt.figure()
    axs = fig.add_subplot(111, aspect='equal')
    for (center, radius) in crcdata:
        crc = plt.Circle(center, radius, edgecolor='blue', facecolor='blue')
        axs.add_patch(crc)
    plt.axis([xa, xb, ya, yb])
    fig.canvas.draw()

.. _apolloniusfig1:

.. figure:: ./apolloniusfig1.png
   :align: center
    
   Three input circles. 

solving polynomial systems
--------------------------

To solve the polynomial systems, we apply the blackbox solver."

::

    from phcpy.solver import solve"

and we need some functions to extract the real solutions.

::

    from phcpy.solutions import strsol2dict, is_real

The `solve4circles` calls the solver on the polynomial systems of the problem.

::

    def solve4circles(syst, verbose=True):
        """
        Given in syst is a list of polynomial systems.
        Returns a list of tuples.  Each tuple in the list of return
        consists of the coordinates of the center and the radius of
        a circle touching the three given circles.
        """
        (circle, eqscnt) = (0, 0)
        result = []
        for eqs in syst:
            eqscnt = eqscnt + 1
            if verbose:
                print('solving system', eqscnt, ':')
                for pol in eqs:
                    print(pol)
            sols = solve(eqs)
            if verbose:
                print('system', eqscnt, 'has', len(sols), 'solutions')
            for sol in sols:
                if is_real(sol, 1.0e-8):
                    soldic = strsol2dict(sol)
                    if soldic['r'].real > 0:
                        circle = circle + 1
                        ctr = (soldic['x'].real, soldic['y'].real)
                        rad = soldic['r'].real
                        result.append((ctr, rad))
                        if verbose:
                            print('solution circle', circle)
                            print('center =', ctr)
                            print('radius =', rad)
        return result

The function `solve4circles` puts the solutions of the polynomial 
system in the format of our problem.
Each solution is a circle, represented by a tuple of the coordinates 
of the center and the radius of the circle.

::

    sols = solve4circles(general_problem)

has as output    

::

    solving system 1 :
    x^2 + y^2 - (r-1)^2;
    (x-2.000000000000000)^2 + y^2 - (r-0.666666666666667)^2;
    (x-1.000000000000000)^2 + (y-1.000000000000000)^2 - (r-0.333333333333333)^2;
    system 1 has 2 solutions
    solution circle 1
    center = (0.792160611810177, -0.734629275680581)
    radius = 2.08036966247227
    solving system 2 :
    x^2 + y^2 - (r-1)^2;
    (x-2.000000000000000)^2 + y^2 - (r-0.666666666666667)^2;
    (x-1.000000000000000)^2 + (y-1.000000000000000)^2 - (r+0.333333333333333)^2;
    system 2 has 2 solutions
    solving system 3 :
    x^2 + y^2 - (r-1)^2;
    (x-2.000000000000000)^2 + y^2 - (r+0.666666666666667)^2;
    (x-1.000000000000000)^2 + (y-1.000000000000000)^2 - (r-0.333333333333333)^2;
    system 3 has 2 solutions
    solution circle 2
    center = (-0.200806137165905, 0.573494560766514)
    radius = 1.60763403126575
    solving system 4 :
    x^2 + y^2 - (r-1)^2;
    (x-2.000000000000000)^2 + y^2 - (r+0.666666666666667)^2;
    (x-1.000000000000000)^2 + (y-1.000000000000000)^2 - (r+0.333333333333333)^2;
    system 4 has 2 solutions
    solution circle 3
    center = (-0.0193166119185703, -0.389367744928919)
    radius = 1.38984660096895
    solving system 5 :
    x^2 + y^2 - (r+1)^2;
    (x-2.000000000000000)^2 + y^2 - (r-0.666666666666667)^2;
    (x-1.000000000000000)^2 + (y-1.000000000000000)^2 - (r-0.333333333333333)^2;
    system 5 has 2 solutions
    solution circle 4
    center = (5.35264994525194, 2.83381218937338)
    radius = 5.05651326763565
    solving system 6 :
    x^2 + y^2 - (r+1)^2;
    (x-2.000000000000000)^2 + y^2 - (r-0.666666666666667)^2;
    (x-1.000000000000000)^2 + (y-1.000000000000000)^2 - (r+0.333333333333333)^2;
    system 6 has 2 solutions
    solution circle 5
    center = (1.86747280383257, 0.159838772566819)
    radius = 0.874300697932419
    solving system 7 :
    x^2 + y^2 - (r+1)^2;
    (x-2.000000000000000)^2 + y^2 - (r+0.666666666666667)^2;
    (x-1.000000000000000)^2 + (y-1.000000000000000)^2 - (r-0.333333333333333)^2;
    system 7 has 2 solutions
    solution circle 6
    center = (1.43293571744453, 2.36388335544507)
    radius = 1.76428097133387
    solution circle 7
    center = (1.23373094922213, 0.96944997788827)
    radius = 0.56905236199947
    solving system 8 :
    x^2 + y^2 - (r+1)^2;
    (x-2.000000000000000)^2 + y^2 - (r+0.666666666666667)^2;
    (x-1.000000000000000)^2 + (y-1.000000000000000)^2 - (r+0.333333333333333)^2;
    system 8 has 2 solutions
    solution circle 8
    center = (1.1821983625488, 0.435483976535281)
    radius = 0.25985684195945

As a summary, let us print the solution circles.

::

    for (idx, circle) in enumerate(sols):
        print('Circle', idx+1, ':', circle)

::

    Circle 1 : ((0.792160611810177, -0.734629275680581), 2.08036966247227)
    Circle 2 : ((-0.200806137165905, 0.573494560766514), 1.60763403126575)
    Circle 3 : ((-0.0193166119185703, -0.389367744928919), 1.38984660096895)
    Circle 4 : ((5.35264994525194, 2.83381218937338), 5.05651326763565)
    Circle 5 : ((1.86747280383257, 0.159838772566819), 0.874300697932419)
    Circle 6 : ((1.43293571744453, 2.36388335544507), 1.76428097133387)
    Circle 7 : ((1.23373094922213, 0.96944997788827), 0.56905236199947)
    Circle 8 : ((1.1821983625488, 0.435483976535281), 0.25985684195945)

Observe that we have a constellation where all eight touching circles
have real coordinates as centers and a positive radius.

In :numref:`apolloniusfig2`
the given circles are plotted as blue disks,
while the eight solution circles are plotted in red,
done by the code below.

::

    fig = plt.figure()
    axs = fig.add_subplot(111, aspect='equal')
    for (center, radius) in crcdata:
        crc = plt.Circle(center, radius, edgecolor='blue', facecolor='blue')
        axs.add_patch(crc)
    for (center, radius) in sols:
        crc = plt.Circle(center, radius, edgecolor='red', facecolor='none')
        axs.add_patch(crc)
    plt.axis([xa, xb, ya, yb])
    fig.canvas.draw()

.. _apolloniusfig2:

.. figure:: ./apolloniusfig2.png
   :align: center
    
   Eight circles touching three given circles.

a special problem
-----------------

In a special configuration of three circles,
the three circles are mutually touching each other.

::

    from math import sqrt
    height = sqrt(3)

The output of

::

   special_problem = polynomials(2, 1, 1, height, 1)
   for pols in special_problem:
       print(pols)

is the following list of eight polynomial systems:

::

    ['x^2 + y^2 - (r-1)^2;', '(x-2.000000000000000)^2 + y^2 - (r-1.000000000000000)^2;', '(x-1.000000000000000)^2 + (y-1.732050807568877)^2 - (r-1.000000000000000)^2;']
    ['x^2 + y^2 - (r-1)^2;', '(x-2.000000000000000)^2 + y^2 - (r-1.000000000000000)^2;', '(x-1.000000000000000)^2 + (y-1.732050807568877)^2 - (r+1.000000000000000)^2;']
    ['x^2 + y^2 - (r-1)^2;', '(x-2.000000000000000)^2 + y^2 - (r+1.000000000000000)^2;', '(x-1.000000000000000)^2 + (y-1.732050807568877)^2 - (r-1.000000000000000)^2;']
    ['x^2 + y^2 - (r-1)^2;', '(x-2.000000000000000)^2 + y^2 - (r+1.000000000000000)^2;', '(x-1.000000000000000)^2 + (y-1.732050807568877)^2 - (r+1.000000000000000)^2;']
    ['x^2 + y^2 - (r+1)^2;', '(x-2.000000000000000)^2 + y^2 - (r-1.000000000000000)^2;', '(x-1.000000000000000)^2 + (y-1.732050807568877)^2 - (r-1.000000000000000)^2;']
    ['x^2 + y^2 - (r+1)^2;', '(x-2.000000000000000)^2 + y^2 - (r-1.000000000000000)^2;', '(x-1.000000000000000)^2 + (y-1.732050807568877)^2 - (r+1.000000000000000)^2;']
    ['x^2 + y^2 - (r+1)^2;', '(x-2.000000000000000)^2 + y^2 - (r+1.000000000000000)^2;', '(x-1.000000000000000)^2 + (y-1.732050807568877)^2 - (r-1.000000000000000)^2;']
    ['x^2 + y^2 - (r+1)^2;', '(x-2.000000000000000)^2 + y^2 - (r+1.000000000000000)^2;', '(x-1.000000000000000)^2 + (y-1.732050807568877)^2 - (r+1.000000000000000)^2;']

::

    specialinput = [((0, 0), 1), ((2, 0), 1), ((1, height), 1)]
    (xa, xb, ya, yb) = (-2, 4, -2, 4)

The code to show the special input is

::

    fig = plt.figure()
    axs = fig.add_subplot(111, aspect='equal')
    for (center, radius) in specialinput:
        crc = plt.Circle(center, radius, edgecolor='blue', facecolor='blue')
        axs.add_patch(crc)
    plt.axis([xa, xb, ya, yb])
    fig.canvas.draw()

which produces :numref:`apolloniusfig3`.

.. _apolloniusfig3:

.. figure:: ./apolloniusfig3.png
   :align: center
    
   Three touching input circles.

The output of

::

    specialsols = solve4circles(special_problem)

is

::

    solving system 1 :
    x^2 + y^2 - (r-1)^2;
    (x-2.000000000000000)^2 + y^2 - (r-1.000000000000000)^2;
    (x-1.000000000000000)^2 + (y-1.732050807568877)^2 - (r-1.000000000000000)^2;
    system 1 has 2 solutions
    solution circle 1
    center = (1.0, 0.577350269189626)
    radius = 2.15470053837925
    solving system 2 :
    x^2 + y^2 - (r-1)^2;
    (x-2.000000000000000)^2 + y^2 - (r-1.000000000000000)^2;
    (x-1.000000000000000)^2 + (y-1.732050807568877)^2 - (r+1.000000000000000)^2;
    system 2 has 1 solutions
    solving system 3 :
    x^2 + y^2 - (r-1)^2;
    (x-2.000000000000000)^2 + y^2 - (r+1.000000000000000)^2;
    (x-1.000000000000000)^2 + (y-1.732050807568877)^2 - (r-1.000000000000000)^2;
    system 3 has 1 solutions
    solving system 4 :
    x^2 + y^2 - (r-1)^2;
    (x-2.000000000000000)^2 + y^2 - (r+1.000000000000000)^2;
    (x-1.000000000000000)^2 + (y-1.732050807568877)^2 - (r+1.000000000000000)^2;
    system 4 has 1 solutions
    solution circle 2
    center = (2.89107059865923e-16, 1.15377761182971e-16)
    radius = 1.0
    solving system 5 :
    x^2 + y^2 - (r+1)^2;
    (x-2.000000000000000)^2 + y^2 - (r-1.000000000000000)^2;
    (x-1.000000000000000)^2 + (y-1.732050807568877)^2 - (r-1.000000000000000)^2;
    system 5 has 1 solutions
    solving system 6 :
    x^2 + y^2 - (r+1)^2;
    (x-2.000000000000000)^2 + y^2 - (r-1.000000000000000)^2;
    (x-1.000000000000000)^2 + (y-1.732050807568877)^2 - (r+1.000000000000000)^2;
    system 6 has 1 solutions
    solution circle 3
    center = (2.0, 7.69185074553436e-17)
    radius = 1.0
    solving system 7 :
    x^2 + y^2 - (r+1)^2;
    (x-2.000000000000000)^2 + y^2 - (r+1.000000000000000)^2;
    (x-1.000000000000000)^2 + (y-1.732050807568877)^2 - (r-1.000000000000000)^2;
    system 7 has 1 solutions
    solution circle 4
    center = (1.0, 1.73205080756888)
    radius = 0.999999999999999
    solving system 8 :
    x^2 + y^2 - (r+1)^2;
    (x-2.000000000000000)^2 + y^2 - (r+1.000000000000000)^2;
    (x-1.000000000000000)^2 + (y-1.732050807568877)^2 - (r+1.000000000000000)^2;
    system 8 has 2 solutions
    solution circle 5
    center = (1.0, 0.577350269189626)
    radius = 0.154700538379251

Let us look closer at the solutions :

::

    for (idx, circle) in enumerate(specialsols):
        print('Circle', idx+1, ':', circle)

::

    Circle 1 : ((1.0, 0.577350269189626), 2.15470053837925)
    Circle 2 : ((2.89107059865923e-16, 1.15377761182971e-16), 1.0)
    Circle 3 : ((2.0, 7.69185074553436e-17), 1.0)
    Circle 4 : ((1.0, 1.73205080756888), 0.999999999999999)
    Circle 5 : ((1.0, 0.577350269189626), 0.154700538379251)

We have five solutions?  Not eight?

The code for the next plot is in

::

    fig = plt.figure()
    axs = fig.add_subplot(111, aspect='equal')
    for (center, radius) in specialinput:
        crc = plt.Circle(center, radius, edgecolor='blue', facecolor='blue')
        axs.add_patch(crc)
    for (center, radius) in specialsols:
        crc = plt.Circle(center, radius, edgecolor='red', facecolor='none')
        axs.add_patch(crc)
    plt.axis([xa, xb, ya, yb])
    fig.canvas.draw()

The plot in :numref:`apolloniusfig4`
shows that the input circles are solutions as well.

.. _apolloniusfig4:

.. figure:: ./apolloniusfig4.png
   :align: center
    
   All circles touching three given touching circles.


a perturbed problem
-------------------

Consider a small perturbation of a special configuration of three circles,
where the three circles are mutually touching each other.

::

    perturbedinput = [((0, 0), 1), ((2.05, 0), 1), ((1.025, height+0.025), 1)]
    perturbed_problem = polynomials(2.05, 1, 1.025, height+0.025, 1)\n",
    perturbedsols = solve4circles(perturbed_problem)"

produces the following output :

::

    solving system 1 :
    x^2 + y^2 - (r-1)^2;
    (x-2.050000000000000)^2 + y^2 - (r-1.000000000000000)^2;
    (x-1.025000000000000)^2 + (y-1.757050807568877)^2 - (r-1.000000000000000)^2;
    system 1 has 2 solutions
    solution circle 1
    center = (1.025, 0.579551408418395)
    radius = 2.17749939915048
    solving system 2 :
    x^2 + y^2 - (r-1)^2;
    (x-2.050000000000000)^2 + y^2 - (r-1.000000000000000)^2;
    (x-1.025000000000000)^2 + (y-1.757050807568877)^2 - (r+1.000000000000000)^2;
    system 2 has 2 solutions
    solving system 3 :
    x^2 + y^2 - (r-1)^2;
    (x-2.050000000000000)^2 + y^2 - (r+1.000000000000000)^2;
    (x-1.025000000000000)^2 + (y-1.757050807568877)^2 - (r-1.000000000000000)^2;
    system 3 has 2 solutions
    solving system 4 :
    x^2 + y^2 - (r-1)^2;
    (x-2.050000000000000)^2 + y^2 - (r+1.000000000000000)^2;
    (x-1.025000000000000)^2 + (y-1.757050807568877)^2 - (r+1.000000000000000)^2;
    system 4 has 2 solutions
    solution circle 2
    center = (0.0248497799767383, -0.00390011791639834)
    radius = 1.02515397552384
    solution circle 3
    center = (-0.309008334843067, -0.198660887619915)
    radius = 1.36735854321414
    solving system 5 :
    x^2 + y^2 - (r+1)^2;
    (x-2.050000000000000)^2 + y^2 - (r-1.000000000000000)^2;
    (x-1.025000000000000)^2 + (y-1.757050807568877)^2 - (r-1.000000000000000)^2;
    system 5 has 2 solutions
    solving system 6 :
    x^2 + y^2 - (r+1)^2;
    (x-2.050000000000000)^2 + y^2 - (r-1.000000000000000)^2;
    (x-1.025000000000000)^2 + (y-1.757050807568877)^2 - (r+1.000000000000000)^2;
    system 6 has 2 solutions
    solution circle 4
    center = (2.35900833484306, -0.19866088761991)
    radius = 1.36735854321413
    solution circle 5
    center = (2.02515022002328, -0.00390011791640729)
    radius = 1.02515397552386
    solving system 7 :
    x^2 + y^2 - (r+1)^2;
    (x-2.050000000000000)^2 + y^2 - (r+1.000000000000000)^2;
    (x-1.025000000000000)^2 + (y-1.757050807568877)^2 - (r-1.000000000000000)^2;
    system 7 has 2 solutions
    solution circle 6
    center = (1.025, 1.73870496299037)
    radius = 1.01834584457851
    solution circle 7
    center = (1.025, 2.04075732867107)
    radius = 1.28370652110219
    solving system 8 :
    x^2 + y^2 - (r+1)^2;
    (x-2.050000000000000)^2 + y^2 - (r+1.000000000000000)^2;
    (x-1.025000000000000)^2 + (y-1.757050807568877)^2 - (r+1.000000000000000)^2;
    system 8 has 2 solutions
    solution circle 8
    center = (1.025, 0.579551408418395)
    radius = 0.177499399150482

Let us look at the solution circles :

::

    for (idx, circle) in enumerate(perturbedsols):
        print('circle', idx+1, ':', circle)

::

    circle 1 : ((1.025, 0.579551408418395), 2.17749939915048)
    circle 2 : ((0.0248497799767383, -0.00390011791639834), 1.02515397552384)
    circle 3 : ((-0.309008334843067, -0.198660887619915), 1.36735854321414)
    circle 4 : ((2.35900833484306, -0.19866088761991), 1.36735854321413)
    circle 5 : ((2.02515022002328, -0.00390011791640729), 1.02515397552386)
    circle 6 : ((1.025, 1.73870496299037), 1.01834584457851)
    circle 7 : ((1.025, 2.04075732867107), 1.28370652110219)
    circle 8 : ((1.025, 0.579551408418395), 0.177499399150482)

:numref:`apolloniusfig5` is made by the code below:

::

    fig = plt.figure()
    axs = fig.add_subplot(111, aspect='equal')
    for (center, radius) in perturbedinput:
        crc = plt.Circle(center, radius, edgecolor='blue', facecolor='blue')
        axs.add_patch(crc)
    for (center, radius) in perturbedsols:
        crc = plt.Circle(center, radius, edgecolor='red', facecolor='none')
        axs.add_patch(crc)
    plt.axis([xa, xb, ya, yb])
    fig.canvas.draw()

.. _apolloniusfig5:

.. figure:: ./apolloniusfig5.png
   :align: center
    
   The solution to the perturbed problem

The solution to the perturbed problem allows to account for the number five 
as the number of touching circles of the special problem: 
the original circles had to be counted twice, as their multiplicity equals two.  
And so we thus have :math:`3 \times 2 + 2 = 8`.

All Lines Tangent to Four Spheres
=================================

Consider all tangent lines to four mutually touching spheres.

The original formulation as polynomial system came from
Cassiano Durand, then at the CS department in Purdue.
The positioning of the centers of the spheres, each with radius
0.5 at the vertices of a tetrahedron came from Thorsten Theobald,
then at TU Muenich.  The centers of the four spheres are

.. math::

   c_1 = (0, 0, 0), \quad
   c_2 = (1, 0, 0), \quad
   c_3 = (1/2, \sqrt{3}/2, 0), \quad
   c_4 = (1/2, \sqrt{3}/6, \sqrt{6}/3).

Let :math:`t = (x_0, x_1, x_2)` be the tangent vector
and :math:`m = (x_3, x_4, x_5)` the moment vector.

The first equation is :math:`\|t\|=1`, the second :math:`m \cdot t = 0`,
the other equations are :math:`\|m - c_i \times t \|^2 - r^2 = 0`
where the radius :math:`r = 1/2`.

::

    from sympy import var, sqrt
    from sympy.vector import CoordSys3D, Vector
    import matplotlib.pyplot as plt
    import numpy as np

    from phcpy.solver import solve
    from phcpy.solutions import coordinates

centers and radii
-----------------

Choices of the centers and radii of four mutually tangent spheres are defined here.

::

    ctr1 = (0, 0, 0)
    ctr2 = (1, 0, 0)
    ctr3 = (0.5, sqrt(3.0)/2, 0)
    ctr4 = (0.5, sqrt(3.0)/6, sqrt(6.0)/3)
    radius = 0.5
    centers = [ctr1, ctr2, ctr3, ctr4]

The choices were made for the suitability of the plot.
Other choices can be found in the paper by Frank Sottile and Thorsten Theobald:
**Line problems in nonlinear computational geometry**.
In *Computational Geometry - Twenty Years Later*, pages 411-432,
edited by J.E. Goodman, J. Pach, and R. Pollack, AMS, 2008.

formulating the equations
-------------------------

We need some vector calculus, done with `sympy`.

::

    N = CoordSys3D('N')
    x0, x1, x2 = var('x0, x1, x2')
    vt = Vector.zero + x0*N.i + x1*N.j + x2*N.k
    normt = vt.dot(vt) - 1
    normt

which produces the first equation

::

    x0**2 + x1**2 + x2**2 - 1

The second equation is

::

    x0*x3 + x1*x4 + x2*x5

is computed by the code

::

    x3, x4, x5 = var('x3, x4, x5')
    vm = Vector.zero + x3*N.i + x4*N.j + x5*N.k
    momvt = vt.dot(vm)

The radii are `[0.5, 0.5, 0.5, 0.5]` defined by

::

    radii = [radius for _ in range(4)]

The polynomial system is constructed by

::

    eqs = [normt, momvt]
    for (ctr, rad) in zip(centers, radii):
        vc = Vector.zero + ctr[0]*N.i + ctr[1]*N.j + ctr[2]*N.k
        left = vm - vc.cross(vt)
        equ = left.dot(left) - rad**2
        eqs.append(equ)

To apply the blackbox solver, we have to convert the polynomials to strings.

::

    fourspheres = []
    print('the polynomial system :')
    for pol in eqs:
        print(pol)
        fourspheres.append(str(pol) + ';')

The output to the above code cell is

::

    the polynomial system :
    x0**2 + x1**2 + x2**2 - 1
    x0*x3 + x1*x4 + x2*x5
    x3**2 + x4**2 + x5**2 - 0.25
    x3**2 + (-x1 + x5)**2 + (x2 + x4)**2 - 0.25
    (-0.866025403784439*x2 + x3)**2 + (0.5*x2 + x4)**2 + (0.866025403784439*x0 - 0.5*x1 + x5)**2 - 0.25
    (-0.816496580927726*x0 + 0.5*x2 + x4)**2 + (0.288675134594813*x0 - 0.5*x1 + x5)**2 + (0.816496580927726*x1 - 0.288675134594813*x2 + x3)**2 - 0.25


So, we have six polynomial equations in six unknowns.

solving the problem
-------------------

Now we call the blackbox solver.

::

    sols = solve(fourspheres)

    for (idx, sol) in enumerate(sols):
        print('Solution', idx+1, ':')
        print(sol)

The solution list is shown below:

::

    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 4
    the solution for t :
     x0 :  1.82013100766029E-16   2.92227989168779E-16
     x1 : -8.16496580927726E-01  -2.50326444773076E-17
     x2 : -5.77350269189626E-01   3.54015053218724E-17
     x3 :  6.04879596033482E-17   3.06586029409515E-17
     x4 :  2.88675134594813E-01  -1.77007526609362E-17
     x5 : -4.08248290463863E-01  -1.25163222386536E-17,
    == err :  4.981E-16 = rco :  1.657E-17 = res :  3.821E-16 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 4
    the solution for t :
     x0 : -7.07106781186547E-01  -2.17839875856796E-16
     x1 : -4.08248290463863E-01   2.08243914343071E-16
     x2 :  5.77350269189626E-01  -1.19547586767375E-16
     x3 :  2.50000000000000E-01  -3.72860037233369E-31
     x4 : -4.33012701892219E-01   2.77333911991762E-31
     x5 : -1.99196604815539E-16   2.50510368981921E-16
    == err :  1.118E-15 = rco :  2.133E-17 = res :  4.441E-16 =
    Solution 3 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 4
    the solution for t :
     x0 :  7.07106781186547E-01   1.37982017054626E-16
     x1 : -4.08248290463863E-01   3.39894708038934E-18
     x2 :  5.77350269189626E-01  -1.66589349202482E-16
     x3 :  2.50000000000000E-01  -4.09837892165604E-31
     x4 : -1.44337567297407E-01   1.66589349202482E-16
     x5 : -4.08248290463863E-01  -5.88982292472640E-17
    == err :  1.023E-15 = rco :  4.667E-17 = res :  3.331E-16 =
    Solution 4 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 4
    the solution for t :
     x0 : -7.07106781186547E-01   8.51323534940577E-17
     x1 :  4.08248290463863E-01  -9.82010876877884E-18
     x2 : -5.77350269189626E-01  -1.11209278833498E-16
     x3 : -2.50000000000000E-01   6.61657084254124E-29
     x4 :  1.44337567297406E-01   1.11209278833581E-16
     x5 :  4.08248290463863E-01  -3.93184175971020E-17
    == err :  2.006E-14 = rco :  1.477E-17 = res :  5.551E-16 =
    Solution 5 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 4
    the solution for t :
     x0 :  7.07106781186548E-01   2.00971395228568E-16
     x1 :  4.08248290463863E-01  -2.45578336010813E-16
     x2 : -5.77350269189626E-01   7.24885788968468E-17
     x3 : -2.50000000000000E-01  -2.77333911991762E-31
     x4 :  4.33012701892219E-01  -3.82104500966428E-31
     x5 : -7.32462262249068E-17  -2.71206918859082E-16
    == err :  7.741E-16 = rco :  4.417E-17 = res :  3.886E-16 =
    Solution 6 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 4
    the solution for t :
     x0 : -3.22358540809185E-16   5.45289082605370E-16
     x1 :  8.16496580927726E-01  -1.64014136987415E-17
     x2 :  5.77350269189625E-01   2.31951016948522E-17
     x3 : -1.74561355056418E-16   2.00875473111053E-17
     x4 : -2.88675134594813E-01  -1.15975508474262E-17
     x5 :  4.08248290463863E-01  -8.20070684937081E-18
    == err :  5.471E-16 = rco :  1.448E-17 = res :  3.682E-16 =

Observe the `m : 4` which indicates the multiplicy four
of each solution.

the tangent lines
-----------------

The solutions contain the components of the tangent
and the moment vectors from which the tangent lines can be computed.

::

    def tangent_lines(solpts, verbose=True):
        """
        Given in solpts is the list of solution points,
        the tuples which respresent the tangent lines
        are returned in a list.
        Each tuple contains a point on the line
        and the tangent vector.
        """
        result = []
        for point in solpts:
            if verbose:
                print(point, end='')
            tan = Vector.zero + point[0]*N.i + point[1]*N.j + point[2]*N.k
            mom = Vector.zero + point[3]*N.i + point[4]*N.j + point[5]*N.k
            pnt = tan.cross(mom) # solves m = p x t
            pntcrd = (pnt.dot(N.i), pnt.dot(N.j), pnt.dot(N.k))
            tancrd = (tan.dot(N.i), tan.dot(N.j), tan.dot(N.k))
            if verbose:
                print(', appending :', pntcrd)
            result.append((pntcrd, tancrd))
        return result

The input to the `tangent_lines` function is computed below:

::

    crd = [coordinates(sol) for sol in sols]
    complexpoints = [values for (names, values) in crd]
    points = []
    for point in complexpoints:
        vals = []
        for values in point:
            vals.append(values.real)
        points.append(tuple(vals))

and then the tangents are computed as 

::

    tangents = tangent_lines(points)

plotting the lines
------------------

Let us first plot the four spheres...

::

    %matplotlib widget
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    R = float(radius)
    x1 = float(ctr1[0]) + R * np.outer(np.cos(u), np.sin(v))
    y1 = float(ctr1[1]) + R * np.outer(np.sin(u), np.sin(v))
    z1 = float(ctr1[2]) + R * np.outer(np.ones(np.size(u)), np.cos(v))
    x2 = float(ctr2[0]) + R * np.outer(np.cos(u), np.sin(v))
    y2 = float(ctr2[1]) + R * np.outer(np.sin(u), np.sin(v))
    z2 = float(ctr2[2]) + R * np.outer(np.ones(np.size(u)), np.cos(v))
    x3 = float(ctr3[0]) + R * np.outer(np.cos(u), np.sin(v))
    y3 = float(ctr3[1]) + R * np.outer(np.sin(u), np.sin(v))
    z3 = float(ctr3[2]) + R * np.outer(np.ones(np.size(u)), np.cos(v))
    x4 = float(ctr4[0]) + R * np.outer(np.cos(u), np.sin(v))
    y4 = float(ctr4[1]) + R * np.outer(np.sin(u), np.sin(v))
    z4 = float(ctr4[2]) + R * np.outer(np.ones(np.size(u)), np.cos(v))
    # Plot the surfaces
    sphere1 = ax.plot_surface(x1, y1, z1, alpha=0.8)
    sphere2 = ax.plot_surface(x2, y2, z2, alpha=0.8)
    sphere3 = ax.plot_surface(x3, y3, z3, alpha=0.8)
    sphere3 = ax.plot_surface(x4, y4, z4, alpha=0.8)
    # Set an equal aspect ratio
    ax.set_aspect('equal')
    plt.show()

The output of the code cell is in :numref:`fourspheresfig1`.

.. _fourspheresfig1:

.. figure:: ./fourspheresfig1.png
   :align: center
    
   Four touching spheres.

The second figure in :numref:`fourspheresfig2` shows the tangent lines.

::

    %matplotlib widget
    ax = plt.figure().add_subplot(projection='3d')
    # range of the tangent lines
    theta = np.linspace(-2.5, 2.5, 10)
    pnt1, tan1 = tangents[0]
    x1 = float(pnt1[0]) + theta*tan1[0]
    y1 = float(pnt1[1]) + theta*tan1[1]
    z1 = float(pnt1[2]) + theta*tan1[2]
    pnt2, tan2 = tangents[1]
    x2 = float(pnt2[0]) + theta*tan2[0]
    y2 = float(pnt2[1]) + theta*tan2[1]
    z2 = float(pnt2[2]) + theta*tan2[2]
    pnt3, tan3 = tangents[2]
    x3 = float(pnt3[0]) + theta*tan3[0]
    y3 = float(pnt3[1]) + theta*tan3[1]
    z3 = float(pnt3[2]) + theta*tan3[2]
    pnt4, tan4 = tangents[3]
    x4 = float(pnt4[0]) + theta*tan4[0]
    y4 = float(pnt4[1]) + theta*tan4[1]
    z4 = float(pnt4[2]) + theta*tan4[2]
    pnt5, tan5 = tangents[4]
    x5 = float(pnt5[0]) + theta*tan5[0]
    y5 = float(pnt5[1]) + theta*tan5[1]
    z5 = float(pnt5[2]) + theta*tan5[2]
    pnt6, tan6 = tangents[5]
    x6 = float(pnt6[0]) + theta*tan6[0]
    y6 = float(pnt6[1]) + theta*tan6[1]
    z6 = float(pnt6[2]) + theta*tan6[2]
    line1 = ax.plot(x1, y1, z1)
    line2 = ax.plot(x2, y2, z2)
    line3 = ax.plot(x3, y3, z3)
    line4 = ax.plot(x4, y4, z4)
    line5 = ax.plot(x5, y5, z5)
    line6 = ax.plot(x6, y6, z6)
    # Set an equal aspect ratio
    ax.set_aspect('equal')
    plt.show()

.. _fourspheresfig2:

.. figure:: ./fourspheresfig2.png
   :align: center
    
   The computed tangent lines.

And then we plot the spheres and the tangent lines:

::

    %matplotlib widget
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    R = float(radius)
    x1 = float(ctr1[0]) + R * np.outer(np.cos(u), np.sin(v))
    y1 = float(ctr1[1]) + R * np.outer(np.sin(u), np.sin(v))
    z1 = float(ctr1[2]) + R * np.outer(np.ones(np.size(u)), np.cos(v))
    x2 = float(ctr2[0]) + R * np.outer(np.cos(u), np.sin(v))
    y2 = float(ctr2[1]) + R * np.outer(np.sin(u), np.sin(v))
    z2 = float(ctr2[2]) + R * np.outer(np.ones(np.size(u)), np.cos(v))
    x3 = float(ctr3[0]) + R * np.outer(np.cos(u), np.sin(v))
    y3 = float(ctr3[1]) + R * np.outer(np.sin(u), np.sin(v))
    z3 = float(ctr3[2]) + R * np.outer(np.ones(np.size(u)), np.cos(v))
    x4 = float(ctr4[0]) + R * np.outer(np.cos(u), np.sin(v))
    y4 = float(ctr4[1]) + R * np.outer(np.sin(u), np.sin(v))
    z4 = float(ctr4[2]) + R * np.outer(np.ones(np.size(u)), np.cos(v))
    # Plot the surfaces
    sphere1 = ax.plot_surface(x1, y1, z1, alpha=0.8)
    sphere2 = ax.plot_surface(x2, y2, z2, alpha=0.8)
    sphere3 = ax.plot_surface(x3, y3, z3, alpha=0.8)
    sphere4 = ax.plot_surface(x4, y4, z4, alpha=0.8)
    # range of the tangent lines
    theta = np.linspace(-2.5, 2.5, 10)
    pnt1, tan1 = tangents[0]
    x1 = float(pnt1[0]) + theta*tan1[0]
    y1 = float(pnt1[1]) + theta*tan1[1]
    z1 = float(pnt1[2]) + theta*tan1[2
    pnt2, tan2 = tangents[1]
    x2 = float(pnt2[0]) + theta*tan2[0]
    y2 = float(pnt2[1]) + theta*tan2[1]
    z2 = float(pnt2[2]) + theta*tan2[2]
    pnt3, tan3 = tangents[2]
    x3 = float(pnt3[0]) + theta*tan3[0]
    y3 = float(pnt3[1]) + theta*tan3[1]
    z3 = float(pnt3[2]) + theta*tan3[2]
    pnt4, tan4 = tangents[3]
    x4 = float(pnt4[0]) + theta*tan4[0]
    y4 = float(pnt4[1]) + theta*tan4[1]
    z4 = float(pnt4[2]) + theta*tan4[2]
    pnt5, tan5 = tangents[4]
    x5 = float(pnt5[0]) + theta*tan5[0]
    y5 = float(pnt5[1]) + theta*tan5[1]
    z5 = float(pnt5[2]) + theta*tan5[2]
    pnt6, tan6 = tangents[5]
    x6 = float(pnt6[0]) + theta*tan6[0]
    y6 = float(pnt6[1]) + theta*tan6[1]
    z6 = float(pnt6[2]) + theta*tan6[2]
    line1 = ax.plot(x1, y1, z1)
    line2 = ax.plot(x2, y2, z2)
    line3 = ax.plot(x3, y3, z3)
    line4 = ax.plot(x4, y4, z4)
    line5 = ax.plot(x5, y5, z5)
    line6 = ax.plot(x6, y6, z6)
    # Set an equal aspect ratio
    ax.axes.set_xlim3d(-1.5, 1.5)
    ax.axes.set_ylim3d(-1.5, 1.5)
    ax.axes.set_zlim3d(-1.5, 1.5)
    ax.set_aspect('equal')
    ax.view_init(elev=30, azim=30, roll=0)
    plt.show()

which produces :numref:`fourspheresfig3`.

.. _fourspheresfig3:

.. figure:: ./fourspheresfig3.png
   :align: center
    
   All lines tangent to four spheres.

Tracking Paths Step by Step
===========================

One main benefit of phcpy is the interactivity and the inversion of control.
Instead of the path tracker deciding the pace of computing points on the path,
the user can ask for the next point on the path,
which is convenient to plot the points on the path.

plotting solutions paths
------------------------

In the example below, the real parts of the solution paths are plotted,
in :numref:`showpathsfig1`.

::

    from phcpy.dimension import set_seed, get_seed
    from phcpy.solutions import strsol2dict
    from phcpy.starters import total_degree_start_system
    from phcpy.trackers import initialize_double_tracker
    from phcpy.trackers import initialize_double_solution
    from phcpy.trackers import next_double_solution

The construction of the homotopy depends on the generation 
of random numbers.  To obtain consistently the same plots, 
the seed of the random number generator is fixed.

::

    set_seed(12871)
    print('the seed :', get_seed())

What is printed is `the seed : 12871`.

The system that will be solved is defined as follows:

::

    p = ['x^2 + y - 3;', 'x + 0.125*y^2 - 1.5;']

For this intersection of two quadrics, we construct 
a total degree start system and compute four start solutions:

::

    q, qsols = total_degree_start_system(p)
    for pol in q:
        print(pol)

Then here is the start system:

::

    x^2 - 1;
    y^2 - 1;
 
::

    print('number of start solutions :', len(qsols))

and `4` is printed.

::

    initialize_double_tracker(p, q, False)

    plt.ion()

The code


::

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
        initialize_double_solution(len(p),startsol)
        dictsol = strsol2dict(startsol)
        xpoints =  [dictsol['x']]
        ypoints =  [dictsol['y']]
        for k in range(300):
            ns = next_double_solution()
            dictsol = strsol2dict(ns)
            xpoints.append(dictsol['x'])
            ypoints.append(dictsol['y'])
            tval = dictsol['t'].real
            if(tval == 1.0):
                break
        print(ns)
        xre = [point.real for point in xpoints]
        yre = [point.real for point in ypoints]
        axs.set_xlim(min(xre)-0.3, max(xre)+0.3)
        axs.set_ylim(min(yre)-0.3, max(yre)+0.3)
        dots, = axs.plot(xre,yre,'r-')
        dots, = axs.plot(xre,yre,'ro')
        fig.canvas.draw()
    fig.canvas.draw()


prints the solutions at the end of the paths

::

    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  1.00000071115204E+00  -4.02767598012276E-06
     y :  1.99999857771171E+00   8.05535854116567E-06
    = err :  2.614E-06 = rco :  1.000E+00 = res :  9.361E-13 =
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  9.99994847664516E-01  -1.64638311049799E-06
     y :  2.00001030465627E+00   3.29275198635159E-06
    == err :  7.006E-06 = rco :  1.000E+00 = res :  1.187E-11 =
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  9.99999462998807E-01   3.01282709881038E-06
     y :  2.00000107400770E+00  -6.02565521719758E-06
    == err :  6.245E-06 = rco :  1.000E+00 = res :  7.729E-12 =
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x : -3.00000000000000E+00   1.39180902310149E-15
     y : -6.00000000000001E+00   5.55716009286787E-15
    == err :  4.185E-07 = rco :  1.000E+00 = res :  2.071E-14 =

and shows then the plot

.. _showpathsfig1:

.. figure:: ./showpathsfig1.png
   :align: center
    
   The real parts of four solution paths.

Typically, with an adaptive step size control, 
the points are closer to each other at the start 
and end of the paths, and where the paths turn.
The `trackers` module exports the original path trackers, 
which use *aposteriori step size control*.  
An aposteriori step size control algorithm 
determines the step size based on the performance of the corrector.

plotting paths and poles
------------------------

The *apriori step size control* determines the step size based
on the location of the nearest pole and the curvature.  
In addition to the (real parts of the paths),
the location of the nearest poles is plotted.

::

   import matplotlib.pyplot as plt

   from phcpy.dimension import set_seed, get_seed
   from phcpy.solutions import strsol2dict
   from phcpy.starters import total_degree_start_system
   from phcpy.curves import set_default_parameters, write_parameters
   from phcpy.curves import initialize_double_artificial_homotopy
   from phcpy.curves import set_double_solution, get_double_solution
   from phcpy.curves import double_predict_correct
   from phcpy.curves import double_t_value, double_closest_pole

The seed is fixed to obtain the same plots in each run.

::

    set_seed(12871)

The system that will be solved is defined as

::

    p = ['x^2 + y - 3;', 'x + 0.125*y^2 - 1.5;']

and a start system based on the total degree is constructed:

::

    x^2 - 1;
    y^2 - 1;

as the output of the code

::

    q, qsols = total_degree_start_system(p)
    for pol in q:
        print(pol)

The list `qsols` has the four solutions of 
the start system `q`.

Before lauching the path trackers,
the parameters must be set.
Default values are used 

::

   set_default_parameters()
   write_parameters()

and shown below:

::

   Values of the HOMOTOPY CONTINUATION PARAMETERS :
    1. gamma : (-0.8063005962200716-0.5915060004219376j)
    2. degree of numerator of Pade approximant    : 5
    3. degree of denominator of Pade approximant  : 1
    4. maximum step size                          : 0.1
    5. minimum step size                          : 1e-06
    6. multiplication factor for the pole radius  : 0.5
    7. multiplication factor for the curvature    : 0.005
    8. tolerance on the residual of the predictor : 0.001
    9. tolerance on the residual of the corrector : 1e-08
   10. tolerance on zero series coefficients      : 1e-12
   11. maximum number of corrector steps          : 4
   12. maximum steps on a path                    : 1000

Then the homotopy is constructed:

::

    initialize_double_artificial_homotopy(p, q, False)"

Then the code below

::

    plt.ion()
    fig1 = plt.figure()
    allpoles = []
    for k in range(len(qsols)):
        if(k == 0):
            axs = fig1.add_subplot(221)
        elif(k == 1):
            axs = fig1.add_subplot(222)
        elif(k == 2):
            axs = fig1.add_subplot(223)
        elif(k == 3):
            axs = fig1.add_subplot(224)
        startsol = qsols[k]
        set_double_solution(len(p), startsol)
        dictsol = strsol2dict(startsol)
        xpoints =  [dictsol['x']]
        ypoints =  [dictsol['y']]
        poles = []
        for k in range(100):
            ns = get_double_solution()
            dictsol = strsol2dict(ns)
            xpoints.append(dictsol['x'])
            ypoints.append(dictsol['y'])
            tval = dictsol['t'].real
            if(tval == 1.0):
                break
            double_predict_correct()
            pole = double_closest_pole()
            tval = double_t_value()
            locp = (tval+pole[0], pole[1])
            poles.append(locp)
        print(ns)
        xre = [point.real for point in xpoints]
        yre = [point.real for point in ypoints]
        axs.set_xlim(min(xre)-0.3, max(xre)+0.3)
        axs.set_ylim(min(yre)-0.3, max(yre)+0.3)
        dots, = axs.plot(xre,yre,'b-')
        dots, = axs.plot(xre,yre,'bo')
        fig1.canvas.draw()
        allpoles.append(poles)
    fig1.canvas.draw()

prints the solutions at the end of the path:

::

    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  1.00000458683462E+00   7.95321001653953E-06
     y :  1.99999082635867E+00  -1.59064683422953E-05
    == err :  1.458E-05 = rco :  3.155E-12 = res :  1.529E-12 =
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  9.99999434968695E-01   9.78667863138325E-07
     y :  2.00000113004483E+00  -1.95736651949408E-06
    == err :  1.658E-05 = rco :  6.664E-12 = res :  1.976E-12 =
    t :  1.00000000000000E+00   0.00000000000000E+00\
    m : 1
    the solution for t :
     x :  1.00000600449478E+00   2.32224203692682E-07
     y :  1.99998799104108E+00  -4.64448075210226E-07
    == err :  1.671E-05 = rco :  9.377E-13 = res :  3.579E-12 =
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x : -3.00000000000000E+00   0.00000000000000E+00
     y : -6.00000000000000E+00   0.00000000000000E+00
    == err :  5.551E-16 = rco :  1.965E-01 = res :  0.000E+00 =

and produces the plot in :numref:`showpolesfig1`.

.. _showpolesfig1:

.. figure:: ./showpolesfig1.png
   :align: center
    
   The real parts of four solution paths.

The poles are shown in :numref:`showpolesfig2`
plotted with the code in

::

    fig2 = plt.figure()
    for k in range(len(qsols)):
        if(k == 0):
            axs = fig2.add_subplot(221)
        elif(k == 1):
            axs = fig2.add_subplot(222)
        elif(k == 2):
            axs = fig2.add_subplot(223)
        elif(k == 3):
            axs = fig2.add_subplot(224)
        poles = allpoles[k]
        pl0 = [pole[0] for pole in poles]
        pl1 = [pole[1] for pole in poles]
        axs.set_xlim(-0.2, 1.2)
        axs.set_ylim(-0.5, 0.5)
        dots, = axs.plot(pl0,pl1,'r+')
        fig2.canvas.draw()
    fig2.canvas.draw()

.. _showpolesfig2:

.. figure:: ./showpolesfig2.png
   :align: center
    
   The nearest poles of four solution paths.

Observe that for this example, more poles are located closer 
to the middle and the end of the paths.

Adjacent 2-by-2 Minors
======================

This section documents the computation of the irreducible factors 
of a pure dimensional solution set, 
using the example of the 2-by-2 adjacent minors of a general 2-by-n matrix.
The 2-by-2 minors of a 2-by-n matrix of indeterminates
originates in algebraic statistics.

::

    from phcpy.families import adjacent_minors
    help(adjacent_minors)

shows

::

    adjacent_minors(rows, cols)
        Returns all adjacent 2-by-2 minors of a general
        matrix of dimensions rows by cols.
        This system originated in a paper on lattice walks and
        primary decomposition, written by P. Diaconis, D. Eisenbud,
        and B. Sturmfels, published by Birkhauser in 1998 in
        Mathematical Essays in Honor of Gian-Carlo Rota,
        edited by B. E. Sagan and R. P. Stanley,
        volume 161 of Progress in Mathematics, pages 173--193.
        See also the paper by S. Hosten and J. Shapiro on
        Primary decomposition of lattice basis ideals, published in 2000
        in the Journal of Symbolic Computation, volume 29, pages 625-639.

Let us look at the simplest nontrivial case:
the adjacent 2-by-2 minors of a 2-by-3 matrix.

::

    pols = adjacent_minors(2, 3)
    for pol in pols:
        print(pol)

shows

::

    x_1_1*x_2_2-x_2_1*x_1_2;
    x_1_2*x_2_3-x_2_2*x_1_3;

We have two polynomials in six variables.  
Therefore, we expect the solution set to be four dimensional.
The two polynomials are quadrics.  
So, the degree of the solution set is expected to be four.
The question is whether the four dimensional solution set 
of degree four is irreducible or not.

computing a witness set
-----------------------

A *witness set* of a positive dimensional solution set consists of

    1. the original polynomial system, 
       augmented with as many linear equations 
       as the dimension of the solution set; and

    2. generic points, as many as the degree of the solution set, 
       computed as solutions of the augmented polynomial system.
 
The *embedding* adds extra slack variables,
which are zero at the generic points.

The witness set data structure reduces the computation of
a positive dimensional solution set to computing the isolated
solutions of one embedded polynomial system.

::

    from phcpy.sets import double_embed
    epols = double_embed(6, 4, pols)
    for pol in epols:
        print(pol)

shows

::

    + x_1_1*x_2_2 - x_2_1*x_1_2 + (-9.99086911101846E-01-4.27240455127090E-02*i)*zz1 + (4.14818420957611E-01-9.09904213439104E-01*i)*zz2 + (-8.98599566975477E-01 + 4.38769664210604E-01*i)*zz3 + (-9.87637111749394E-01 + 1.56757569180295E-01*i)*zz4;
    - x_2_2*x_1_3 + x_1_2*x_2_3 + (-9.50915060423189E-01 + 3.09452012209265E-01*i)*zz1 + (8.94934777859038E-01-4.46196978226426E-01*i)*zz2 + (4.28909177508030E-01-9.03347617171477E-01*i)*zz3 + (9.99997394966914E-01 + 2.28255545098161E-03*i)*zz4;
    zz1;
    zz2;
    zz3;
    zz4;
    + (2.52642772862563E-01-1.64562207779935E-01*i)*x_1_1 + (2.38172268997908E-01-1.84886617118381E-01*i)*x_2_2 + (-1.91482090565462E-01-2.32902769201594E-01*i)*x_2_1 + (7.30954525688645E-02 + 2.92516915276440E-01*i)*x_1_2 + (-3.84959957593650E-02-2.99043724594892E-01*i)*x_2_3 + (-2.71276842664000E-01-1.31597741406691E-01*i)*x_1_3 + (-2.64408104251273E-04-3.01511228642393E-01*i)*zz1 + (2.72300014364985E-01 + 1.29467343704581E-01*i)*zz2 + (-2.99939451765406E-01-3.07476207820803E-02*i)*zz3 + (-2.95299363009813E-01 + 6.08882346195858E-02*i)*zz4 - 3.01511344577764E-01;
    + (-4.58249807046168E-01-5.40906757392526E-02*i)*x_1_1 + (4.14862726733922E-02 + 4.11506360074024E-01*i)*x_2_2 + (-2.32018034572266E-01 + 9.48639573945220E-02*i)*x_2_1 + (-2.80184386911322E-01-4.28818448828079E-01*i)*x_1_2 + (-1.20422832345040E-01-1.42772536416847E-01*i)*x_2_3 + (-1.75000947987830E-01-2.52077630538672E-02*i)*x_1_3 + (1.00853655248805E-01-1.59726738098879E-01*i)*zz1 + (1.47157089788836E-01 + 9.31710575065957E-02*i)*zz2 + (-3.18937022210321E-02-2.67242905524245E-01*i)*zz3 + (-1.70655540108348E-01 + 3.05332146451031E-02*i)*zz4+(-9.36346290234469E-02 + 2.17662695080044E-01*i);
    + (-8.19176308116668E-02 + 9.31779891695395E-02*i)*x_1_1 + (-2.23021019174348E-01 + 9.08742650844521E-02*i)*x_2_2 + (-2.55417845632647E-01 + 1.92893995983262E-01*i)*x_2_1 + (2.85202344419800E-02 + 4.06697473728353E-01*i)*x_1_2 + (-3.26980478298163E-01 + 5.39551406506719E-02*i)*x_2_3 + (-4.72487633926490E-03 + 3.03185406062837E-01*i)*x_1_3 + (2.49469923033290E-01 + 7.22359362814556E-02*i)*zz1 + (3.57054498446798E-01 + 1.84302643405729E-01*i)*zz2 + (2.44627852456931E-01-2.19409644826369E-02*i)*zz3 + (-6.45211959942027E-02-2.81391298463502E-01*i)*zz4+(1.78511770624327E-02-2.88585493131170E-01*i);
    + (-2.47934564140415E-01 + 3.78474051361417E-01*i)*x_1_1 + (2.57696009847760E-01-1.97353959339253E-01*i)*x_2_2 + (-1.32713498879035E-01-6.68083554253052E-02*i)*x_2_1 + (1.71277270038441E-01-1.69671423207279E-01*i)*x_1_2 + (2.45498609773496E-01-3.00565440004147E-01*i)*x_2_3 + (-1.33200278933953E-01 + 9.57961946142634E-02*i)*x_1_3 + (2.53856038005847E-01 + 2.28882119288700E-01*i)*zz1 + (3.08364764282092E-03-3.24002619621010E-01*i)*zz2 + (2.00774472927742E-01 + 2.60472818035180E-01*i)*zz3 + (-2.51374429242122E-01 + 8.28490424271253E-03*i)*zz4+(-7.60100938321149E-02-1.82187404809906E-01*i);

Now we compute the second part of the witness set, 
using the blackbox solver.

::

    from phcpy.solver import solve
    esols = solve(epols)
    for (idx, sol) in enumerate(esols):
        print('Solution', idx+1, ':')
        print(sol)

shows four generic points on the four dimensional solution set:

::

    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x_1_1 :  1.26098995657825E+00   1.23960253907080E-01
     x_2_2 :  6.08716792661783E-02  -3.72692628829057E-01
     x_2_1 : -1.49366498111755E+00  -1.72487935488432E+00
     x_1_2 :  1.17926528430272E-01   1.73403649424959E-01
     zz1 :  2.01661798954576E-33  -1.36621089687380E-31
     zz2 :  4.18934314423693E-32   1.19013821639201E-31
     zz3 :  9.09330358138888E-33  -3.97471598577921E-32
     zz4 :  9.41147475091108E-32  -5.25062461176711E-32
     x_1_3 : -5.43007038294243E-01   4.77552159454532E-01
     x_2_3 :  1.30126856409899E+00   4.91781165629461E-02
    == err :  4.606E-15 = rco :  1.659E-02 = res :  1.221E-15 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x_1_1 :  7.84232277025150E-01   2.84185732210647E-02
     x_2_2 :  1.09164662176987E-01  -6.94770373095709E-01
     x_2_1 :  3.30315020994409E-02  -9.45083305517748E-01
     x_1_2 :  5.76431528154133E-01   9.13299754350158E-02
     zz1 : -1.29868775020467E-32  -1.23245957744023E-32
     zz2 : -1.02927640426824E-32  -1.56443669914841E-32
     zz3 : -1.96403692130211E-32  -7.02597925370710E-33
     zz4 :  0.00000000000000E+00   0.00000000000000E+00
     x_1_3 : -6.32506552290242E-01   1.49836349424459E-01
     x_2_3 :  1.81539704759152E-01   7.61970172619267E-01
    == err :  1.716E-15 = rco :  2.536E-02 = res :  1.499E-15 =
    Solution 3 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x_1_1 : -1.07029652224429E+00  -2.84414365427269E-02
     x_2_2 :  4.26584686826893E-01   4.61346027545957E-01
     x_2_1 : -8.74524024788936E-02  -7.16513058512809E-01
     x_1_2 :  7.70137841287209E-01  -5.24903704207380E-01
     zz1 : -2.60172998762874E-31  -7.95558581219627E-32
     zz2 :  1.94241870643223E-31  -7.60476764044638E-32
     zz3 : -1.18161218057147E-31   6.47719980338975E-32
     zz4 :  3.21269292411328E-31   9.43964867510126E-32
     x_1_3 :  1.88771487997930E+00   2.52881775952061E+00
     x_2_3 : -1.49855102312210E+00   1.51018382382141E+00
    == err :  5.217E-15 = rco :  1.384E-02 = res :  1.749E-15 =
    Solution 4 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x_1_1 :  1.03149817185821E+00   1.48840565152077E-01
     x_2_2 : -1.97069785649029E-31  -4.19546322472859E-32
     x_2_1 : -1.96423581053286E+00  -2.24826401526688E+00
     x_1_2 : -6.05197132478501E-32   6.74819240986552E-32
     zz1 : -4.35083051092473E-32   8.73329103580651E-32
     zz2 :  6.61155389590362E-32   7.61859548953403E-32
     zz3 :  2.10093447239302E-32  -1.28492740537369E-32
     zz4 :  0.00000000000000E+00   0.00000000000000E+00
     x_1_3 : -1.40341064895811E-01   1.23712580175709E+00
     x_2_3 :  1.45874399778069E+00   6.42372628489355E-02\
    == err :  4.557E-15 = rco :  1.108E-02 = res :  1.332E-15 =

As expected, we find four solutions, 
equal to the degree of the solution set.

monodromy breakup
-----------------

The *numerical irreducible decomposition* 
of a pure dimensional solution set is a list of tuples, 
each tuple represents an irreducible component, with two elements

1. the list of labels to the generic points in the witness set; and

2. the certificate of the linear trace.
 
The *monodromy breakup* algorithm refines the witness set, 
partitioning the generic points in the witness set corresponding 
to the irreducible components.  The stop test in the monodromy looping
algorithm is provided by the linear trace, 
which serves as a certificate for the numerical computations.

The breakup algorithm runs in verbose mode:

::

    from phcpy.factor import double_monodromy_breakup
    deco = double_monodromy_breakup(epols, esols, 4, verbose=True)

In verbose mode, the progress of the algorithm is printed:

::

    ... running monodromy loops in double precision ...
    ... initializing the grid for the linear trace ...
    The diagnostics of the trace grid :
      largest error on the samples : 1.056863766729658e-14
      smallest distance between the samples : 1.1622396157009902
    ... starting loop 1 ...
    new permutation : [2, 1, 3, 4]
    number of factors : 4 -> 3
    the decomposition :
      factor 1 : ([1, 2], 0.2775442892800747)
      factor 2 : ([3], 0.2775442892800704)
      factor 3 : ([4], 9.992007221626409e-16)
    the permutation :  2 1 3 4 : 4 -> 3
    calculated sum at samples : -6.39839733359069E-02   1.09505498482705E-01
    value at the linear trace :  1.07936961323818E-02  -9.32611213290812E-02
    calculated sum at samples : -1.22710696579410E-01   2.89111792077992E-02
    value at the linear trace : -1.97488366047695E-01   2.31677799019584E-01
    calculated sum at samples :  2.68414486121697E-01  -5.27451814780890E-01
    value at the linear trace :  2.68414486121697E-01  -5.27451814780890E-01
    Certifying with linear trace test...
    calculated sum at samples : -6.39839733359069E-02   1.09505498482705E-01
    value at the linear trace :  1.07936961323818E-02  -9.32611213290812E-02
    The witness points 1 2 do not define a factor.
    The factorization cannot be certified.
    ... starting loop 2 ...
    new permutation : [1, 2, 3, 4]
    number of factors : 3 -> 3
    the decomposition :
      factor 1 : ([1, 2], 0.2775442892800747)
      factor 2 : ([3], 0.2775442892800704)
      factor 3 : ([4], 9.992007221626409e-16)
    ... starting loop 3 ...
    the permutation :  1 2 3 4 : 3 -> 3
    calculated sum at samples : -6.39839733359069E-02   1.09505498482705E-01
    value at the linear trace :  1.07936961323818E-02  -9.32611213290812E-02
    calculated sum at samples : -1.22710696579410E-01   2.89111792077992E-02
    value at the linear trace : -1.97488366047695E-01   2.31677799019584E-01
    calculated sum at samples :  2.68414486121697E-01  -5.27451814780890E-01
    value at the linear trace :  2.68414486121697E-01  -5.27451814780890E-01
    Certifying with linear trace test...
    calculated sum at samples : -6.39839733359069E-02   1.09505498482705E-01
    value at the linear trace :  1.07936961323818E-02  -9.32611213290812E-02
    The witness points 1 2 do not define a factor.
    The factorization cannot be certified.
    new permutation : [2, 3, 1, 4]
    number of factors : 3 -> 2
    the decomposition :
      factor 1 : ([1, 2, 3], 4.163336342344337e-15)
      factor 2 : ([4], 9.992007221626409e-16)
    the permutation :  2 3 1 4 : 3 -> 2
    calculated sum at samples : -1.86694669915317E-01   1.38416677690504E-01
    value at the linear trace : -1.86694669915313E-01   1.38416677690503E-01
    calculated sum at samples :  2.68414486121697E-01  -5.27451814780890E-01
    value at the linear trace :  2.68414486121697E-01  -5.27451814780890E-01
    Certifying with linear trace test...
    calculated sum at samples : -1.86694669915317E-01   1.38416677690504E-01
    value at the linear trace : -1.86694669915313E-01   1.38416677690503E-01
    The witness points 1 2 3 define a factor.
    calculated sum at samples :  2.68414486121697E-01  -5.27451814780890E-01
    value at the linear trace :  2.68414486121697E-01  -5.27451814780890E-01
    The witness points 4 define a factor.
    The factorization is certified.
    calculated sum at samples : -1.86694669915317E-01   1.38416677690504E-01
    value at the linear trace : -1.86694669915313E-01   1.38416677690503E-01
    calculated sum at samples :  2.68414486121697E-01  -5.27451814780890E-01
    value at the linear trace :  2.68414486121697E-01  -5.27451814780890E-01

As a summary, the contents of `deco` is written as

::

    from phcpy.factor import write_decomposition
    write_decomposition(deco)

which produces

::

    factor 1 : ([1, 2, 3], 4.163336342344337e-15)
    factor 2 : ([4], 9.992007221626409e-16)


There are two irreducible factors, one of degree three,
and another of degree one.  
The floating-point certificates are close to machine precision.

Design of a moving 7-bar mechanism
==================================

Laurent polynomial systems are systems that have negative exponents.
In this section, we consider a Laurent system with one irreducible
component of degree three and six isolated points.

A reference for the general case is the paper by Carlo Innocenti:
**Polynomial solution to the position analysis of the 7-line Assur kinematic
chain with one quaternary link**, in *Mech. Mach. Theory*, Vol. 30, No. 8,
pages 1295-1303, 1995.

The special case was introduced in the paper with title:
**Numerical decomposition of the solution sets of polynomial
systems into irreducible components**, *SIAM J. Numer. Anal.* 38(6):2022-2046,
2001, by Andrew Sommese, Jan Verschelde, and Charles Wampler.

This special sevenbar mechanism has 6 isolated solutions and a cubic curve."

SymPy is used to define the equations,
with complex arithmetic.

::

    from cmath import exp
    from sympy import var

From `phcpy`, the following functions are imported:

::

    from phcpy.solutions import coordinates, diagnostics, condition_tables
    from phcpy.solver import solve
    from phcpy.sets import double_laurent_membertest
    from phcpy.cascades import double_laurent_top_cascade
    from phcpy.cascades import double_laurent_cascade_filter
    from phcpy.factor import double_monodromy_breakup


a Laurent polynomial system
---------------------------

The code in this section defines the Laurent polynomial system
for a generic instance of the parameters.

::

    def symbolic_equations():
        """
        Returns the symbolic equations,
        with parameters a1, a2, a3, a4, a5, a6
        b0, b2, b3, b4, b5, and c0, with variables
        t1, t2, t3, t4, and t5. 
        """
        a0, a1, a2, a3, a4, a5, a6 = var('a0, a1, a2, a3, a4, a5, a6')
        b0, b2, b3, b4, b5, c0 = var('b0, b2, b3, b4, b5, c0')
        t1, t2, t3, t4, t5, t6 = var('t1, t2, t3, t4, t5, t6')
        eq1 = a1*t1 + a2*t2 - a3*t3 - a0
        eq2 = b2*t2 + a3*t3 - a4*t4 + a5*t5 - b0
        eq3 = a4*t4 + b5*t5 - a6*t6 - c0
        return [eq1, eq2, eq3]

Then the symbolic equations are computed via

::

    eqs = symbolic_equations()
    for equ in eqs:
        print(equ)

with output in       

::

    -a0 + a1*t1 + a2*t2 - a3*t3
    a3*t3 - a4*t4 + a5*t5 - b0 + b2*t2
    a4*t4 - a6*t6 + b5*t5 - c0

A generic instance of the problem is
defined in the following function:

::

    def generic_problem(eqs):
        """
        Given the symbolic equations in eqs,
        defines the equations for a generic problem,
        as a Laurent polynomial system.
        The system is returned as a list of string representations,
        suitable for input to the solve of phcpy.
        """
        i = complex(0, 1)
        subdict = {a0: 0.7 + 0.2*i, b0: 0.6, c0: 0.5 - 0.5*i, \
            a1: 0.7, a2: 0.8, b2: 0.6 + 0.5*i, a3: 0.4, a4: 0.6, \
            a5: 0.8, b5: 0.4 + 0.3*i, a6: 0.9}
        print(subdict)
        conjugates = {a0: 0.7 - 0.2*i, b0: 0.6, c0: 0.5 + 0.5*i, \
            a1: 0.7, a2: 0.8, b2: 0.6 - 0.5*i, a3: 0.4, a4: 0.6, \
            a5: 0.8, b5: 0.4 - 0.3*i, a6: 0.9}
        print(conjugates)
        result = []
        for equ in eqs:
            pol = equ.subs(subdict)
            result.append(str(pol) + ';')
        for equ in eqs:
            pol = str(equ.subs(conjugates))
            pol = pol.replace('t1', 't1**(-1)')
            pol = pol.replace('t2', 't2**(-1)')
            pol = pol.replace('t3', 't3**(-1)')
            pol = pol.replace('t4', 't4**(-1)')
            pol = pol.replace('t5', 't5**(-1)')
            pol = pol.replace('t6', 't6**(-1)')
            result.append(pol + ';')
        return result

Then the system is constructed symbolically via

::

    T1, T2, T3, T4, T5, T6 = var('T1, T2, T3, T4, T5, T6')
    generic = generic_problem(eqs)
    for equ in generic:
        print(equ)

with output

::

      0.7*t1 + 0.8*t2 - 0.4*t3 - 0.7 - 0.2*I;
      t2*(0.6 + 0.5*I) + 0.4*t3 - 0.6*t4 + 0.8*t5 - 0.6;
      0.6*t4 + t5*(0.4 + 0.3*I) - 0.9*t6 - 0.5 + 0.5*I;
      0.7*t1**(-1) + 0.8*t2**(-1) - 0.4*t3**(-1) - 0.7 + 0.2*I;
      t2**(-1)*(0.6 - 0.5*I) + 0.4*t3**(-1) - 0.6*t4**(-1) + 0.8*t5**(-1) - 0.6;
      0.6*t4**(-1) + t5**(-1)*(0.4 - 0.3*I) - 0.9*t6**(-1) - 0.5 - 0.5*I;

Observe the negative exponents.
Now, let us call the blackbox solver:

::

    sols = solve(generic)
    print('found', len(sols), 'solutions')

which prints `found 18 solutions`.

A condition table is a frequency table of the `err`,
`rco`, and `res` fields of the solutions.

::

    condition_tables(sols)

with output in

::

    ([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 2],
     [3, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11, 7])

The first and third row indicate that the forward and backward errors
of the solutions were all very small.  The second row indicates that
the estimates for the condition numbers were also in the good range.
So, 18 solutions are well conditioned.

a special problem
-----------------

As we have as many equations as variables, 
for general coefficients, we will have only isolated solutions.
For special parameters, the system has an irreducible cubic as a solution set.

::

    def special_parameters():
        """
        Returns a dictionary with special values for the parameters
        for the Assur7c in Roberts Cognate pattern.
        Before calling this function, the symbolic_equations()
        must have defined the variables for the parameters.
        """
        i = complex(0, 1)
        # start with the independent parameters
        result = {b0: 0.0, c0: 1.2, a2: 0.46, \
            b2: -0.11 + 0.49*i, a5: 0.41}
        theta4 = 0.6 + 0.8*i
        theta3 = exp(1.8*i)
        # add the derived parameters
        result[a3] = result[a5]
        beta = result[b2]/result[a2]
        result[a0] = result[c0]/beta
        result[b5] = result[a5]*beta
        result[a4] = abs(result[b2])
        result[a1] = abs(result[a0] + result[a3]*theta3 - result[a4]*theta4/beta)
        result[a6] = abs(result[a4]*theta4 - result[b5]*theta3-result[c0])
        return result

::

    def conjugates(dic):
        """
        Given on input a dictionary with variables as keys
        and complex numbers as values.
        Returns a dictionary with the same keys,
        but with values replaced by complex conjugates.
        """
        result = {}
        for key in list(dic.keys()):
            result[key] = dic[key].conjugate()
        return result

::

    def special_problem(eqs):
        """
        Given the symbolic equations in eqs,
        replaces the parameters with special values.
        """
        pars = special_parameters()
        conj = conjugates(pars)
        result = []
        for equ in eqs:
            pol = equ.subs(pars)
            result.append(str(pol) + ';')
        for equ in eqs:
            pol = str(equ.subs(conj))
            pol = pol.replace('t1', 't1**(-1)')
            pol = pol.replace('t2', 't2**(-1)')
            pol = pol.replace('t3', 't3**(-1)')
            pol = pol.replace('t4', 't4**(-1)')
            pol = pol.replace('t5', 't5**(-1)')
            pol = pol.replace('t6', 't6**(-1)')
            result.append(pol + ';')
        return result

Constructing the polynomials of the special problem

::

    special = special_problem(eqs)
    for equ in special:
        print(equ)

leads to

::

    0.710358341606049*t1 + 0.46*t2 - 0.41*t3 + 0.240761300555115 + 1.07248215701824*I;\n",
    t2*(-0.11 + 0.49*I) + 0.41*t3 - 0.502195181179589*t4 + 0.41*t5;\n",
    0.502195181179589*t4 + t5*(-0.0980434782608696 + 0.436739130434783*I) - 0.775518556663656*t6 - 1.2;\n",
    0.710358341606049*t1**(-1) + 0.46*t2**(-1) - 0.41*t3**(-1) + 0.240761300555115 - 1.07248215701824*I;\n",
    t2**(-1)*(-0.11 - 0.49*I) + 0.41*t3**(-1) - 0.502195181179589*t4**(-1) + 0.41*t5**(-1);\n",
    0.502195181179589*t4**(-1) + t5**(-1)*(-0.0980434782608696 - 0.436739130434783*I) - 0.775518556663656*t6**(-1) - 1.2;\n"


Running the `solve` of the `solver` module and printing
the number of solutions

::

    sols = solve(special)
    print('found', len(sols), 'solutions')"

shows `found 6 solutions`.
Let us look at all solutions, executing

::

   for (idx, sol) in enumerate(sols):
       print('Solution', idx+1, ':')
       print(sol)

which gives the output

::

    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -5.66058532229597E-01  -3.67759063106022E-01
     t2 :  5.87958784214731E-01  -2.31024744989669E-01
     t3 :  2.66141319725173E-01   1.71943916132782E+00
     t4 :  2.42136902539446E-01   1.56435563278512E+00
     t5 : -8.79136932379914E-02  -5.67977370543055E-01
     t6 : -1.05957839449340E+00   1.03531112257767E+00
    == err :  3.073E-15 = rco :  3.135E-02 = res :  9.853E-16 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  3.53717078322846E-01  -5.30879118400565E-01
     t2 : -5.99365365166588E-01  -1.57888836421937E+00
     t3 :  5.27607584716191E-01  -7.54168308853111E-02
     t4 :  5.86169848996842E-01  -8.37877878416834E-02
     t5 : -1.85739737768408E+00   2.65498503011424E-01
     t6 : -1.08247082572556E+00  -1.13383016827930E+00
    == err :  2.542E-15 = rco :  3.053E-02 = res :  4.510E-16 =
    Solution 3 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  8.69193907041793E-01  -1.30453665759482E+00
     t2 : -2.10146778358548E-01  -5.53582709999063E-01
     t3 :  1.85739737768408E+00  -2.65498503011422E-01
     t4 :  1.67183103323242E+00  -2.38973437749056E-01
     t5 : -5.27607584716191E-01   7.54168308853116E-02
     t6 : -4.40509781239126E-01  -4.61410384022367E-01
    == err :  8.299E-16 = rco :  2.577E-02 = res :  4.302E-16 =
    Solution 4 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -5.45421549881402E-01  -8.38161877518281E-01
     t2 : -5.49340801920855E-01  -8.35598398361888E-01
     t3 : -9.74098087752273E-01   2.26125884049935E-01
     t4 :  1.78912046655936E-01  -9.83865071827120E-01
     t5 :  4.72153176961024E-02  -9.98884734979395E-01
     t6 : -8.74935009483375E-01  -4.84240363022669E-01
    == err :  2.024E-15 = rco :  1.693E-01 = res :  6.661E-16 =
    Solution 5 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -1.24225818334613E+00  -8.07075027813303E-01
     t2 :  1.47332994921904E+00  -5.78910775622768E-01
     t3 :  8.79136932379908E-02   5.67977370543054E-01
     t4 :  9.66290808831226E-02   6.24284218493861E-01
     t5 : -2.66141319725174E-01  -1.71943916132781E+00
     t6 : -4.82817017275396E-01   4.71759173981634E-01
    == err :  2.046E-15 = rco :  3.752E-02 = res :  9.992E-16 =
    Solution 6 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  2.80836878017556E-01  -9.59755514673061E-01
     t2 : -9.99161738567148E-01   4.09367827689738E-02
     t3 : -4.72153176961015E-02   9.98884734979395E-01
     t4 :  9.35633636119240E-01  -3.52972660360954E-01
     t5 :  9.74098087752273E-01  -2.26125884049935E-01
     t6 : -9.37276397924231E-01   3.48587082225057E-01
    == err :  1.949E-15 = rco :  1.566E-01 = res :  7.910E-16 =

The output of

::

    condition_tables(sols)

is

::

    ([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 1],
     [2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6])

The `solve` of the solver module misses the component, but finds all isolated solutions.


a numerical irreducible decomposition
-------------------------------------

A numerical irreducible decomposition for this system augments 
the system with a linear equation and adds one slack variable.  
A cascade of homotopies find generic points 
on all positive dimensional components of the solution set.

The cascade is wrapped in the following function

::

    def embed_and_cascade(pols, topdim):
        """
        Computes and solves an embedding at top dimension topdim
        of the Laurent polynomials in pols, before running one
        step in the cascade homotopy.
        Returns the embedded system, the three generic points,
        and the filtered solutions at the end of the cascade.
        """
        (embpols, sols0, sols1) \
            = double_laurent_top_cascade(len(pols), topdim, pols, 1.0e-08)
        print('the top generic points :')
        for (idx, sol) in enumerate(sols0):
            print('Solution', idx+1, ':')
            print(sol)
        print('the nonsolutions :')
        for (idx, sol) in enumerate(sols1):
            print('Solution', idx+1, ':')
            print(sol)
        print('... running cascade step ...')
        (embdown, nsols1, sols2) = double_laurent_cascade_filter(1, embpols, \
            sols1, 1.0e-8)
        filtsols2 = []
        for (idx, sol) in enumerate(nsols1):
            err, rco, res = diagnostics(sol)
            if res < 1.0e-8 and rco > 1.0e-8:
                _, point = coordinates(sol)
                crdpt = []
                for pt in point:
                    crdpt.append(pt.real)
                    crdpt.append(pt.imag)
                onset = double_laurent_membertest(embpols, sols0, 1, crdpt)
                if not onset:
                    filtsols2.append(sol)
        print('... after running the cascade ...')
        for (idx, sol) in enumerate(filtsols2):
            print('Solution', idx+1, ':')
            print(sol)
        print('found %d isolated solutions' % len(filtsols2))
        return (embpols, sols0, filtsols2)

Running the code in

::

    (embpols, sols0, isosols) = embed_and_cascade(special, 1)

produces the following output:

::

    the top generic points :
    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -1.38546903964616E-01  -6.90082116709782E-01
     t2 :  1.73335780614296E+00  -1.87923014523194E+00
     t3 :  2.29192181084370E+00  -6.88217799479036E-01
     zz1 :  6.81370071005700E-17   1.65228741665612E-16
     t4 :  1.45392357364500E+00   2.10288883797136E+00
     t5 : -2.29192181084370E+00   6.88217799479037E-01
     t6 : -7.03671420729004E-01  -1.59719770241003E-02
    == err :  1.018E-14 = rco :  7.104E-03 = res :  1.619E-15 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  9.95322784254993E-01   2.68098302301308E-02
     t2 : -2.29577075969873E-01  -2.62541647600257E-01
     t3 :  2.05412606828065E+00   2.36770142844668E+00
     zz1 :  7.58249556746807E-16  -7.20127010307089E-16
     t4 :  3.06452334567059E-01  -1.66495396855090E-01
     t5 : -2.05412606828065E+00  -2.36770142844668E+00
     t6 :  2.44172639794733E-01  -9.65280235905487E-01
    == err :  3.577E-15 = rco :  2.031E-02 = res :  1.998E-15 =
    Solution 3 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -2.82166573126121E-01   3.94770824953965E-01
     t2 : -1.34882663123749E-01  -2.78088349356538E+00
     t3 : -5.29856181081608E-02   1.79767069529071E-01
     zz1 :  2.43319655939393E-15   2.43565365274284E-15
     t4 :  2.74289769478703E+00   4.77512904043170E-01
     t5 :  5.29856181081606E-02  -1.79767069529072E-01
     t6 :  3.23379011328058E-01   3.61784458285676E-01
    == err :  6.752E-15 = rco :  2.149E-03 = res :  2.746E-15 =
    the nonsolutions :
    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  2.06608664978716E-01  -5.16879012501813E-01
     t2 : -2.63044555866456E+00  -3.32754563190078E+00
     t3 :  2.63517362859822E-02  -3.68638950244269E-01
     zz1 :  2.85704892320864E-01  -1.03310141444983E+00
     t4 :  3.69849125743338E+00  -3.24986222144375E-01
     t5 :  2.43555779349728E+00   2.20779169965979E+00
     t6 : -9.87691023055072E-01   2.23467268352943E+00
    == err :  3.403E-15 = rco :  2.764E-02 = res :  1.561E-15 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  3.09054535971898E-01   3.45124987046649E-02
     t2 : -2.62571332002632E-01   3.96698883381851E-02
     t3 :  4.32699725536785E+00   5.00970765602759E+00
     zz1 :  2.88323234772462E-01  -1.68997834567399E+00
     t4 : -1.07803480589631E+00  -1.04944159716503E+00
     t5 : -1.51252934187501E+00  -6.41341577084612E+00
     t6 :  1.32322533488782E+00   1.47762805502765E+00
    == err :  6.549E-15 = rco :  1.025E-02 = res :  3.109E-15 =
    Solution 3 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -4.38810628559853E-02  -3.73466290987713E-01
     t2 :  1.41948412243664E+00   2.64381540059301E-01
     t3 : -1.51714321706050E-01   2.82809534841746E+00
     zz1 :  5.80561465793845E-01   7.55877838451024E-01
     t4 :  6.84618807253357E-01  -1.45281628643651E-01
     t5 :  2.97376151600431E-01  -2.76806627220378E+00
     t6 : -3.91051079512042E-01  -5.02441609168954E-01
    == err :  3.533E-15 = rco :  1.180E-02 = res :  1.086E-15 =
    Solution 4 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -5.56053317223021E-01  -1.10476501709143E+00
     t2 :  1.08198361813419E+00  -6.10107625522898E-01
     t3 :  1.63797645826712E+00   2.65500692732987E-01
     zz1 : -3.75750870070081E-02  -3.41463274433182E-01
     t4 :  1.13327119828200E+00   1.00882972129605E+00
     t5 :  8.79031826440765E-02  -8.01012123191631E-01
     t6 : -2.97523163734206E-01   1.24044339222691E+00
    == err :  1.312E-15 = rco :  2.758E-02 = res :  1.305E-15 =
    Solution 5 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  1.42858931075131E-01   4.57096616580137E-01
     t2 : -1.11064648867271E-01  -1.72825737676934E-01
     t3 :  2.13308888441550E+00   1.81873215330586E+00
     zz1 : -7.57293747789502E-01  -3.06703567333995E-01
     t4 :  2.49672678368341E-01  -1.36328639570869E-01
     t5 : -1.84534929207840E+00  -3.88014927109987E+00
     t6 :  2.03214759107748E+00  -3.03532919113466E-01
    == err :  7.406E-15 = rco :  6.149E-03 = res :  1.971E-15 =
    Solution 6 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -4.75108750576191E-02  -1.75815265025158E-02
     t2 : -3.39314399121161E+00   1.41556569134256E+01
     t3 : -8.54021034392022E+00  -1.38811633117136E+01
     zz1 : -1.13042166796985E+01   7.26160237898825E+00
     t4 :  4.31260891858428E-03  -2.20764606490238E-02
     t5 :  1.89710631587227E-02  -1.71967991631555E-02
     t6 :  1.24227998831989E+01  -1.02611559314054E+01
    == err :  1.583E-14 = rco :  1.029E-05 = res :  1.698E-14 =
    Solution 7 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  9.24280274976326E-01  -1.03328843987402E+00
     t2 : -5.37485466298478E-01  -6.93996812904627E-01
     t3 :  2.16870928120656E+00   3.05270120118574E-01
     zz1 :  1.72164678001449E-03  -2.61490044232207E-01
     t4 :  1.12688433550202E+00   6.04135856872397E-01
     t5 : -1.14709448190074E+00   7.21699657198301E-01
     t6 : -1.06009576375345E+00  -9.36416226240182E-03
    == err :  1.235E-15 = rco :  1.803E-02 = res :  1.473E-15 =
    Solution 8 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  2.56141977779414E-01  -1.26192262870549E-01
     t2 :  2.72348811313250E-02   2.27695964490337E+00
     t3 :  8.94179859360769E-02   1.64582692269313E-02
     zz1 : -1.69636690132142E+00   1.17291673807606E+00
     t4 :  1.44816411386502E-01  -5.21926135720004E-02
     t5 : -1.06047458297333E+00  -2.70699717824994E+00
     t6 :  2.29314581678167E+00  -1.93548691562723E+00
    == err :  3.630E-15 = rco :  1.266E-03 = res :  2.873E-15 =
    Solution 9 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -5.55849610226793E-01  -6.71303249412053E-01
     t2 :  2.22897995402302E+00  -6.36069101208733E-01
     t3 : -6.01231159130768E-03   3.20619531479791E-01
     zz1 :  1.91319173054801E-01   8.69593783238856E-01
     t4 :  1.46082467405615E+00   1.30115181744724E+00
     t5 : -2.81326700289881E-01  -5.36286320208369E-01
     t6 : -5.80370651666996E-01  -3.51674733202245E-01
    == err :  1.098E-15 = rco :  2.893E-02 = res :  1.563E-15 =
    Solution 10 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  1.96420789232812E-01  -7.25547806896222E-01
     t2 :  4.99261380871303E-01   6.57994832350889E-02
     t3 :  1.79119494447843E-01   8.39067157393641E-02
     zz1 : -2.92932332587968E-01   7.12592654282584E-01
     t4 :  1.31985073905489E+00   5.86875406583772E-01
     t5 : -2.16641539477234E-01  -1.59643525554850E-01
     t6 : -2.56051374441468E-01  -6.62533448829890E-01
    == err :  4.769E-16 = rco :  7.725E-03 = res :  1.402E-15 =
    Solution 11 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -2.60846935289829E-01  -1.87501112823166E+00
     t2 :  8.91205034115317E-01   2.18216656676799E+00
     t3 :  1.35354427992550E-01   2.14130902473316E-01
     zz1 : -4.38444725918741E-01   6.37863637010385E-01
     t4 :  4.58998133742581E-01   3.40221617387902E-01
     t5 :  1.48618652575694E+00  -8.83671361618052E-01
     t6 : -4.27731643521526E-01   3.12634407365499E-01
    == err :  4.740E-15 = rco :  1.451E-02 = res :  2.900E-15 =
    Solution 12 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  2.58526368145855E-02  -5.89746935862251E-01
     t2 :  5.89197128653301E-01  -8.56471323438483E-01
     t3 :  4.99885121099350E-01  -8.51253650022305E-02
     zz1 : -1.40247285577845E-01   4.15696687647816E-01
     t4 :  1.66773190687860E+00   8.55179987286000E-01
     t5 : -3.91334700102683E-01   1.44828575079480E-01
     t6 : -3.52629174262130E-01  -2.31227315105538E-01
    == err :  3.496E-15 = rco :  1.796E-02 = res :  1.877E-15 =
    Solution 13 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  1.31675077902931E+00  -2.33335629490662E+00
     t2 :  2.60288716555775E-01  -6.40025748193045E-02
     t3 :  3.47791530603780E+00   3.39861528263726E+00
     zz1 :  1.78921624866824E+00  -9.20506290946835E-01
     t4 : -3.11182618390869E-02   1.27100852726915E-01
     t5 : -1.76684760979132E-01   1.89283790316579E-02
     t6 : -3.78392322972795E+00   1.30979364462129E+00
    == err :  5.729E-15 = rco :  1.880E-03 = res :  3.067E-15 =
    Solution 14 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  1.09751178036435E+00  -1.40574539820329E+00
     t2 : -6.26603734750229E-01   2.77315044013342E+00
     t3 :  1.16737341994583E-01   1.23560628499014E-01
     zz1 : -9.17943053708135E-01   1.14574950444327E+00
     t4 :  2.24794207999588E-01  -2.20864560734679E-01
     t5 :  7.19265347716089E-03  -2.96930745959293E-01
     t6 : -1.46858876699108E-01  -1.65019111957513E+00
    == err :  9.467E-16 = rco :  1.061E-02 = res :  1.874E-15 =
    Solution 15 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -1.30458602699707E+00  -2.33541493076135E+00
     t2 : -1.95199047692026E-01  -4.46316027103696E-01
     t3 : -2.63391619655473E-01  -2.31639399176850E-01
     zz1 :  3.72490200270314E-01  -8.90356341436546E-01
     t4 :  2.41781694569341E+00   1.09012175794204E-01
     t5 :  4.97584553326169E+00   7.63254192555898E-01
     t6 : -1.44788918406800E+00   3.95223854445634E+00
    == err :  1.808E-15 = rco :  2.602E-02 = res :  1.971E-15 =
    Solution 16 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  5.06717666964459E-01  -7.97241780193331E-01
     t2 : -4.84953123213085E-01  -4.42303335072244E-01
     t3 :  9.44820387296986E-01  -1.44512677763047E-02
     zz1 : -2.86872573899956E-01   1.14217012524121E-01
     t4 :  1.47330976885270E+00   3.99769298368997E-01
     t5 : -2.57072989456668E-01   3.67302107948625E-01
     t6 : -4.07703118258568E-01  -1.02532173760480E-01
    == err :  1.760E-15 = rco :  2.381E-02 = res :  1.249E-15 =
    Solution 17 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -4.44607365409264E-02  -2.59406017880835E-01
     t2 :  3.78176637095005E-01   4.39268560397966E-01
     t3 : -1.79503189927382E+00   3.20305158897225E+00
     zz1 :  6.51015867427647E-01   9.37169953082408E-01
     t4 :  1.30934173849756E+00  -1.12550388002783E+00
     t5 :  2.25665017724524E+00  -2.76677636166897E+00
     t6 : -3.40266607063703E-01  -2.61585700440509E-01
    == err :  7.075E-15 = rco :  6.393E-03 = res :  2.526E-15 =
    Solution 18 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -2.53488386324696E-01  -2.83321822997287E-01
     t2 :  4.84338946463415E-01  -1.67152891112537E+00
     t3 : -4.67320096443835E-01  -4.71307905321838E-02
     zz1 :  7.80393042398170E-02   4.84167415461090E-01
     t4 :  2.41580881604416E+00   6.84660143775337E-01
     t5 :  4.73818299916806E-01   3.62365346793219E-01
     t6 : -3.86544038004814E-01   4.76118474348958E-02
    == err :  3.492E-15 = rco :  2.148E-02 = res :  1.291E-15 =
    Solution 19 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -1.33188225629716E-01   2.10966347401895E-02
     t2 :  7.81125047936108E-02  -6.85635891220529E-02
     t3 : -1.17563370850181E+00  -1.23520999866901E+00
     zz1 : -1.16759331681420E+00   1.23236695756379E+00
     t4 :  2.02453946582577E+00   5.47406959303102E-01
     t5 : -7.18144929758816E-02  -1.30397505946681E-01
     t6 :  1.24906033996779E+00  -1.34990480788646E+00
    == err :  2.561E-15 = rco :  1.448E-03 = res :  5.407E-15 =
    Solution 20 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  7.93825359960543E-01   3.42942332165835E-01
     t2 :  1.61651669536702E-01  -2.20420239291763E+00
     t3 :  4.64285871839968E+00   3.38285150308267E+00
     zz1 :  5.85868505889160E-01  -1.37232304716824E+00
     t4 :  4.40119366407873E-02   3.19449005533791E-01
     t5 : -3.57045518779285E+00  -3.30999449881425E+00
     t6 :  1.53654644664973E-01   4.28064658906619E-01
    == err :  3.056E-15 = rco :  1.438E-02 = res :  3.844E-15 =
    Solution 21 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -7.60489144456182E-02  -1.04300898291791E+00
     t2 :  1.30846500932422E+00   6.09095986219939E-02
     t3 :  2.60207018922791E+00   2.05366973304254E+00
     zz1 :  3.31325254376896E-01  -4.47605721252978E-01
     t4 :  5.02974219619557E-02   5.96655221299233E-01
     t5 : -8.46448240981595E-01  -2.38906321758213E+00
     t6 : -4.52524421602189E-01   8.14559254657865E-01
    == err :  1.498E-15 = rco :  2.468E-02 = res :  2.047E-15 =
    Solution 22 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  1.17748055149989E+00  -2.43215785209726E+00
     t2 : -1.11238676756403E+00   1.89832325651549E+00
     t3 :  3.71807527705904E+00   1.18513377244720E+00
     zz1 : -1.36994965965414E-01  -9.86163388302379E-01
     t4 : -5.27247356627358E-01  -7.14664957510453E-02
     t5 : -1.69492044434729E-01  -4.08782307130586E-01
     t6 : -1.38104781392739E+00   1.16797791507032E+00
    == err :  2.038E-15 = rco :  4.306E-02 = res :  1.249E-15 =
    Solution 23 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  3.12364985327808E+00  -2.05920731130273E+00
     t2 : -7.83955013846558E+00   3.44568114237919E+00
     t3 :  2.41311468318769E-01  -1.43824031559247E-02
     zz1 : -1.59787001054699E+00  -6.62918168050137E-01
     t4 :  1.10281361738592E+00  -3.47888278489006E+00
     t5 :  3.62179107231530E+00   1.85699344604169E+00
     t6 : -2.26905454742053E-01   2.75923139264321E-01
    == err :  1.042E-14 = rco :  4.252E-03 = res :  2.163E-15 =
    Solution 24 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  3.08355119326816E+00  -8.60649488234224E+00
     t2 : -5.82366068292758E-01   6.35509246767222E+00
     t3 : -6.51342246123743E-02   9.06094708166825E-02
     zz1 :  2.84990305286854E+00   1.14816597303183E+00
     t4 :  4.34422945742045E-02   1.52135217854099E-01
     t5 :  6.74998864667608E+00   9.94704932815523E+00
     t6 : -1.17348229611637E+01   1.39527617933237E+00
    == err :  3.194E-14 = rco :  5.571E-04 = res :  6.391E-15 =
    Solution 25 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  3.99353498759161E-01  -1.05992514732995E+00
     t2 : -6.97576291634198E-01   2.90503809343657E-01
     t3 :  7.12853582076897E-02   4.30854721150634E-01
     zz1 : -1.84001525983760E-01   2.70200187419495E-01
     t4 :  1.46156773352004E+00  -1.50613908910889E-01
     t5 :  1.12282954533061E+00   4.33717117615447E-02
     t6 : -5.52341835801268E-01   1.66701841647516E-01
    == err :  9.415E-16 = rco :  4.214E-02 = res :  7.841E-16 =
    Solution 26 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  1.68976453558350E-01   4.67689996418420E-01
     t2 : -1.11289042160100E+00  -2.36087487730329E+00
     t3 :  7.14519021473384E-02   1.63214961074351E-01
     zz1 : -3.02872131405798E-01  -6.49699828630491E-02
     t4 :  2.59822631620836E+00   2.73920748230829E-02
     t5 : -5.72162586858151E-02  -1.87011959227602E-01
     t6 :  6.42736202016055E-01   6.82645862607672E-02
    == err :  2.532E-15 = rco :  2.192E-03 = res :  2.207E-15 =
    Solution 27 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  5.98088716113923E-01  -1.09928306831808E+00
     t2 : -1.32166665802307E+00   9.68674190067784E-02
     t3 :  1.25565072505728E+00   7.88279907441350E-01
     zz1 : -1.94309421837523E-01  -4.14018780046576E-01
     t4 :  1.15900115504862E+00  -3.40883097941626E-01
     t5 :  7.68268330674074E-01  -3.30650582363870E-01
     t6 : -4.24188159828689E-01   7.70802547523126E-01
    == err :  1.539E-15 = rco :  3.512E-02 = res :  7.355E-16 =
    Solution 28 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
    " t1 :  1.74056293633403E+00   1.02810794797457E+00
     t2 : -4.53362942425264E+00  -3.21821306993922E+00
     t3 :  9.07980903275410E-01   2.06897987372685E+00
     zz1 :  9.08970112275608E-02  -1.10894123925719E+00
     t4 :  2.76795451024622E+00  -1.85436435929005E+00
     t5 :  8.31726772379002E-02  -3.06761087764399E-01
     t6 :  3.80053540200631E-01   3.19279779559526E-01
    == err :  5.467E-15 = rco :  1.088E-02 = res :  4.621E-15 =
    Solution 29 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  9.98996829187838E-02   5.45879605814945E+00
     t2 :  4.77713088756544E-02  -1.65064542093410E-01
     t3 :  7.03293595992006E-02  -5.75960210248865E-03
     zz1 : -4.34989172011554E+00   2.22547996626009E+00
     t4 :  7.42601541252893E-02  -4.06346066756681E-02
     t5 : -8.26936045369789E+00  -8.88179809231708E+00
     t6 :  9.96592014933793E+00  -6.77636361154879E+00
    == err :  2.141E-14 = rco :  3.614E-04 = res :  4.914E-15 =
    Solution 30 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  1.26911098400795E+00  -1.26566539297069E+00
     t2 : -8.06762092689544E-01   6.05358763159280E-01
     t3 :  4.36774346983282E+00   1.58729813776646E+00
     zz1 : -2.24457954433668E-01  -1.01428467390231E+00
     t4 : -5.06895216434008E-01   2.00062592651442E-01
     t5 : -2.24941445482227E+00  -1.41444572972009E+00
     t6 : -4.23738271342406E-01   3.28737704966059E-01
    == err :  2.782E-15 = rco :  2.863E-02 = res :  1.360E-15 =
    ... running cascade step ...
    Solution at position 14 is not appended.
    Solution at position 22 is not appended.
    ... after running the cascade ...
    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -5.45421549881402E-01  -8.38161877518281E-01
     t2 : -5.49340801920855E-01  -8.35598398361888E-01
     t3 : -9.74098087752273E-01   2.26125884049935E-01
     t6 : -8.74935009483375E-01  -4.84240363022669E-01
     t4 :  1.78912046655936E-01  -9.83865071827120E-01
     t5 :  4.72153176961023E-02  -9.98884734979395E-01
    == err :  3.858E-16 = rco :  3.659E-02 = res :  1.388E-16 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -5.66058532229596E-01  -3.67759063106022E-01
     t2 :  5.87958784214732E-01  -2.31024744989669E-01
     t3 :  2.66141319725174E-01   1.71943916132782E+00
     t6 : -1.05957839449340E+00   1.03531112257767E+00
     t4 :  2.42136902539447E-01   1.56435563278512E+00
     t5 : -8.79136932379912E-02  -5.67977370543055E-01
    == err :  1.140E-15 = rco :  2.128E-02 = res :  4.441E-16 =
    Solution 3 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  2.80836878017556E-01  -9.59755514673061E-01
     t2 : -9.99161738567148E-01   4.09367827689737E-02
     t3 : -4.72153176961024E-02   9.98884734979395E-01
     t6 : -9.37276397924231E-01   3.48587082225057E-01
     t4 :  9.35633636119240E-01  -3.52972660360955E-01
     t5 :  9.74098087752273E-01  -2.26125884049935E-01
    == err :  2.120E-13 = rco :  2.187E-02 = res :  2.220E-16 =
    Solution 4 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -1.24225818334613E+00  -8.07075027813304E-01
     t2 :  1.47332994921904E+00  -5.78910775622768E-01
     t3 :  8.79136932379910E-02   5.67977370543054E-01
     t6 : -4.82817017275396E-01   4.71759173981634E-01
     t4 :  9.66290808831228E-02   6.24284218493860E-01
     t5 : -2.66141319725174E-01  -1.71943916132781E+00
    == err :  2.043E-13 = rco :  1.402E-02 = res :  2.637E-16 =
    Solution 5 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  8.69193907041793E-01  -1.30453665759482E+00
     t2 : -2.10146778358548E-01  -5.53582709999063E-01
     t3 :  1.85739737768408E+00  -2.65498503011422E-01
     t6 : -4.40509781239126E-01  -4.61410384022367E-01
     t4 :  1.67183103323242E+00  -2.38973437749056E-01
     t5 : -5.27607584716191E-01   7.54168308853117E-02
    == err :  6.392E-16 = rco :  1.136E-02 = res :  2.220E-16 =
    Solution 6 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  3.53717078322846E-01  -5.30879118400565E-01
     t2 : -5.99365365166589E-01  -1.57888836421937E+00
     t3 :  5.27607584716190E-01  -7.54168308853108E-02
     t6 : -1.08247082572556E+00  -1.13383016827931E+00
     t4 :  5.86169848996842E-01  -8.37877878416832E-02
     t5 : -1.85739737768408E+00   2.65498503011425E-01
    == err :  5.809E-16 = rco :  1.811E-02 = res :  6.661E-16 =
    found 6 isolated solutions

The polynomials in the embedded system are printed by

::

    for pol in embpols:
        print(pol)

with output in

::

     + 7.10358341606049E-01*t1 + 4.60000000000000E-01*t2 - 4.10000000000000E-01*t3 + (-3.99034499614690E-01 + 9.16935912764493E-01*i)*zz1+(2.40761300555115E-01 + 1.07248215701824E+00*i);
     + (-1.10000000000000E-01 + 4.90000000000000E-01*i)*t2 + 4.10000000000000E-01*t3 - 5.02195181179589E-01*t4 + 4.10000000000000E-01*t5 + (-2.71603706128330E-01-9.62409178477302E-01*i)*zz1;
     + 5.02195181179589E-01*t4 + (-9.80434782608696E-02 + 4.36739130434783E-01*i)*t5 - 7.75518556663656E-01*t6 + (-9.98029449494905E-01 + 6.27472544490738E-02*i)*zz1 - 1.20000000000000E+00;
     + (-9.58594885352021E-01-2.84773323499490E-01*i)*zz1+(2.40761300555115E-01-1.07248215701824E+00*i) - 4.10000000000000E-01*t3^-1 + 4.60000000000000E-01*t2^-1 + 7.10358341606049E-01*t1^-1;
     + (9.19472457329587E-01-3.93154422857342E-01*i)*zz1 + 4.10000000000000E-01*t5^-1 - 5.02195181179589E-01*t4^-1 + 4.10000000000000E-01*t3^-1 + (-1.10000000000000E-01-4.90000000000000E-01*i)*t2^-1;
     + (4.61876615285503E-01 + 8.86944187788841E-01*i)*zz1 - 1.20000000000000E+00 - 7.75518556663656E-01*t6^-1 + (-9.80434782608696E-02-4.36739130434783E-01*i)*t5^-1 + 5.02195181179589E-01*t4^-1;
     + (9.04264731472024E-01-4.26972241973442E-01*i)*t1 + (9.93763178327739E-01-1.11511189572842E-01*i)*t2 + (8.26568833449879E-01 + 5.62835645254728E-01*i)*t3 + (8.13748932516514E-01 + 5.81216547276687E-01*i)*t4 + (9.60611770393053E-01 + 2.77893912459997E-01*i)*t5 + (-5.80987954131586E-02 + 9.98310838352234E-01*i)*t6 + (1.86846899315945E-02 + 9.99825425943029E-01*i)*zz1+(-9.99671453748700E-01 + 2.56317100475435E-02*i);

Observe that we have seven equation in seven variables,
where the last variable is the slack variable `zz1`.

The code in

::

    print('the isolated solutions :')
    for (idx, sol) in enumerate(isosols):
        print('Solution', idx+1, ':')
        print(sol)

produces the following output:

::

    the isolated solutions :
    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -5.45421549881402E-01  -8.38161877518281E-01
     t2 : -5.49340801920855E-01  -8.35598398361888E-01
     t3 : -9.74098087752273E-01   2.26125884049935E-01
     t6 : -8.74935009483375E-01  -4.84240363022669E-01
     t4 :  1.78912046655936E-01  -9.83865071827120E-01
     t5 :  4.72153176961023E-02  -9.98884734979395E-01
    == err :  3.858E-16 = rco :  3.659E-02 = res :  1.388E-16 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -5.66058532229596E-01  -3.67759063106022E-01
     t2 :  5.87958784214732E-01  -2.31024744989669E-01
     t3 :  2.66141319725174E-01   1.71943916132782E+00
     t6 : -1.05957839449340E+00   1.03531112257767E+00
     t4 :  2.42136902539447E-01   1.56435563278512E+00
     t5 : -8.79136932379912E-02  -5.67977370543055E-01
    == err :  1.140E-15 = rco :  2.128E-02 = res :  4.441E-16 =
    Solution 3 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  2.80836878017556E-01  -9.59755514673061E-01
     t2 : -9.99161738567148E-01   4.09367827689737E-02
     t3 : -4.72153176961024E-02   9.98884734979395E-01
     t6 : -9.37276397924231E-01   3.48587082225057E-01
     t4 :  9.35633636119240E-01  -3.52972660360955E-01
     t5 :  9.74098087752273E-01  -2.26125884049935E-01
    == err :  2.120E-13 = rco :  2.187E-02 = res :  2.220E-16 =
    Solution 4 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 : -1.24225818334613E+00  -8.07075027813304E-01
     t2 :  1.47332994921904E+00  -5.78910775622768E-01
     t3 :  8.79136932379910E-02   5.67977370543054E-01
     t6 : -4.82817017275396E-01   4.71759173981634E-01
     t4 :  9.66290808831228E-02   6.24284218493860E-01
     t5 : -2.66141319725174E-01  -1.71943916132781E+00
    == err :  2.043E-13 = rco :  1.402E-02 = res :  2.637E-16 =
    Solution 5 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  8.69193907041793E-01  -1.30453665759482E+00
     t2 : -2.10146778358548E-01  -5.53582709999063E-01
     t3 :  1.85739737768408E+00  -2.65498503011422E-01
     t6 : -4.40509781239126E-01  -4.61410384022367E-01
     t4 :  1.67183103323242E+00  -2.38973437749056E-01
     t5 : -5.27607584716191E-01   7.54168308853117E-02
    == err :  6.392E-16 = rco :  1.136E-02 = res :  2.220E-16 =
    Solution 6 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     t1 :  3.53717078322846E-01  -5.30879118400565E-01
     t2 : -5.99365365166589E-01  -1.57888836421937E+00
     t3 :  5.27607584716190E-01  -7.54168308853108E-02
     t6 : -1.08247082572556E+00  -1.13383016827931E+00
     t4 :  5.86169848996842E-01  -8.37877878416832E-02
     t5 : -1.85739737768408E+00   2.65498503011425E-01
    == err :  5.809E-16 = rco :  1.811E-02 = res :  6.661E-16 =

::

    fac = double_monodromy_breakup(embpols, sols0, 1, islaurent=True, verbose=True)

produces the following output:

::

    ... running monodromy loops in double precision ...
    ... initializing the grid for the linear trace ...
    The diagnostics of the trace grid :
      largest error on the samples : 1.1914895647759973e-14
      smallest distance between the samples : 2.3530539196335214
    ... starting loop 1 ...
    new permutation : [3, 2, 1]
    number of factors : 3 -> 2
    the decomposition :
      factor 1 : ([1, 3], 0.33179448122518296)
      factor 2 : ([2], 0.3317944812254451)
    ... starting loop 2 ...
    new permutation : [2, 3, 1]
    number of factors : 2 -> 1
    the decomposition :
      factor 1 : ([1, 2, 3], 4.11226608321158e-13)

The grouping of the generic sets in one set 
shows that the solution curve is an irreducible cubic.
