Tangent Lines via Witness Set Intersection
==========================================

Let us compute the tangents lines to a circle
via a witness set intersection.
With matplotlib a plot is made of the tangent lines.

::

    from math import cos, sin, pi
    from random import uniform
    import matplotlib.pyplot as plt
    from phcpy.solver import solve
    from phcpy.solutions import make_solution, coordinates, strsol2dict
    from phcpy.sets import double_embed
    from phcpy.diagonal import double_diagonal_solve

a witness set of the circle
---------------------------

Consider the following system of polynomial equations:

.. math::

   \left\{
      \begin{array}{rcl}
         (x-a)^2 + (y-b)^2 & = & r^2 \\
                         y & = & s x
      \end{array}
   \right.

which represents a circle centered 
at :math:`(a, b)` with radius :math:`r`, 
intersected by a line with slope :math:`s`.
The code below plots an example, 
a circle centered at :math:`(3, 2)`, with radius :math:`1`, 
intersected with the line :math:`y = x`.

:numref:`touchcirclefig1` is made by the code below

::

    fig1 = plt.figure()
    axs = fig1.gca()
    center = (3, 2)
    radius = 1
    circle = plt.Circle(center, radius, edgecolor='blue', facecolor='none')
    axs.add_artist(circle)
    plt.plot([0, 4], [0, 4], 'r')
    plt.plot([2, 3], [2, 3], 'go')
    plt.axis([0, 5, 0, 4])
    plt.show()

.. _touchcirclefig1:

.. figure:: ./touchcirclefig1.png
   :align: center
    
   A circle intersected with a line.

The system above has two equations in three variables 
:math:`(x, y, s)` and thus defines a curve.  
What is the degree of this curve?

::

    def polynomials(a, b, r):
        """
        Returns string representations of two polynomials:
        1) a circle with radius r centered at (a, b);
        2) a line through the origin with slope s.
        """
        crc = '(x - %.15e)^2 + (y - %.15e)^2 - %.15e;' % (a, b, r**2)
        lin = 'y - s*x;'
        return [crc, lin]

::

    def make_witness_set(pols, verbose=True):
        """
        We have two equations in three variables in pols
        and therefore we expect a one dimensional set.
        """
        embpols = double_embed(3, 1, pols)
        if verbose:
            print('the embedded system :')
            for pol in embpols:
                print(pol)
        embsols = solve(embpols)
        if verbose:
            print('the witness points :')
            for sol in embsols:
                print(sol)
        return (embpols, embsols)

::

    def special_solutions(pols, slope):
        """
        Given in pols the polynomials for the line intersecting a circle
        for a general slope s, solves the problem for a special numerical
        value for s, given in the slope.
        The value of the slope is added as last coordinate to each solution.
        """
        special = []
        for pol in pols:
            rpl = pol.replace('s', '(' + str(slope) + ')')
            special.append(rpl)
        sols = solve(special)
        result = []
        for sol in sols:
            (vars, vals) = coordinates(sol)
            vars.append('s')
            vals.append(slope)
            extsol = make_solution(vars, vals)
            result.append(extsol)
        return result

Consider a circle with radius :math:`1`, centered at :math:`(3, 2)`, 
intersected with the line :math:`y = x`.

::

    def circle_line_set():
        """
        Generates the system and its witness set.
        Returns the witness set for a fixed circle
        intersect with a one parameter family of lines
        as a tuple of polynomials and solutions.
        """
        syst = polynomials(3, 2, 1)
        for pol in syst:
            print(pol)
        spsols = special_solutions(syst, 1)
        for (idx, sol) in enumerate(spsols):
            print('Solution', idx+1, ':')
            print(sol)
        (embsyst, embsols) = make_witness_set(syst, False)
        print('the polynomials in the witness set:')
        for pol in embsyst:
            print(pol)
        print('the solutions :')
        for (idx, sol) in enumerate(embsols):
            print('Solution', idx+1, ':')
            print(sol)
        print('degree of the set :', len(embsols))
        return (embsyst, embsols)

The output of

::

    (witpols, witsols) = circle_line_set()

is

::

    (x - 3.000000000000000e+00)^2 + (y - 2.000000000000000e+00)^2 - 1.000000000000000e+00;
    y - s*x;
    Solution 1 :
    t : 0.000000000000000E+00 0.000000000000000E+00
    m : 1
    the solution for t :
     x : 3.000000000000000E+00  0.000000000000000E+00
     y : 3.000000000000000E+00  0.000000000000000E+00
     s : 1.000000000000000E+00  0.0
    == err : 0.000E+00 = rco : 1.000E+00 = res : 0.000E+00 =
    Solution 2 :
    t : 0.000000000000000E+00 0.000000000000000E+00
    m : 1
    the solution for t :
     x : 2.000000000000000E+00  0.000000000000000E+00
     y : 2.000000000000000E+00  0.000000000000000E+00
     s : 1.000000000000000E+00  0.0
    == err : 0.000E+00 = rco : 1.000E+00 = res : 0.000E+00 =
    the polynomials in the witness set:
     + x^2 + y^2 - 6*x - 4*y + (-2.25888608394754E-01 + 9.74153138165392E-01*i)*zz1 + 12;
     - x*s + y + (4.03414981775719E-02 + 9.99185950424038E-01*i)*zz1;
    zz1;
     + (-1.19285340138354E-01 + 9.92860014114818E-01*i)*x + (-8.44394833617596E-01-5.35721350106483E-01*i)*y + (2.84736834956935E-01-9.58605724382401E-01*i)*s + (-8.64325704104395E-01-5.02932477798403E-01*i)*zz1+(4.42123925817800E-01-8.96953975530214E-01*i);
    the solutions :
    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  4.71992186119667E+00  -2.06479191154844E+00
     y :  4.21048963164569E+00   1.60655842789469E+00
     zz1 :  1.47124630003083E-32  -4.65633532178132E-32
     s :  6.23787940798390E-01   6.13262424188266E-01
    == err :  3.253E-15 = rco :  4.216E-02 = res :  1.554E-15 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  2.22786304688919E+00  -3.25328777253638E-01
     y :  1.21728610339575E+00   3.20932555200177E-01
     zz1 :  8.77875924683185E-33  -8.43345264643039E-33
     s :  5.14387214185304E-01   2.19168552262572E-01
    == err :  6.116E-16 = rco :  1.366E-01 = res :  2.220E-16 =
    Solution 3 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x : -3.19419508001181E-01  -6.36087292918970E-01
     y :  1.33420600110974E+00   3.17131210618635E+00
     zz1 :  0.00000000000000E+00   0.00000000000000E+00
     s : -4.82279862042460E+00  -3.24310772611730E-01
    == err :  3.542E-16 = rco :  7.767E-02 = res :  8.327E-16 =
    Solution 4 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x : -1.02740100906214E+00  -1.33614560499347E+00
     y :  3.37463784588945E+00  -3.91462680435862E+00
     zz1 :  0.00000000000000E+00   0.00000000000000E+00
     s :  6.20734137916362E-01   3.00295170717090E+00
    == err :  5.028E-16 = rco :  5.750E-02 = res :  1.110E-15 =\
    degree of the set : 4

Thus, the degree of the algebraic curve defined by the system equals four.

intersecting with the Jacobian
------------------------------

For any random slope :math:`s`, the line :math:`y = s x` will intersect 
the circle in exactly two (complex) points, 
as solutions of the polynomial system.  
A tangent line to the circle intersects the circle in one double solution, 
that is: a solution that has to be counted twice.  
We look for values of the slope 
for which the intersection points are singular. 
Singular points satisfy the original equations and the equations 
defined by all partial derivatives of the system, called the Jacobian.

First, a witness set for the Jacobian is constructed.

Let :math:`f_1` and :math:`f_2` denote the two polynomials 
in :math:`(x, y, s)`, then the Jacobian is

.. math::

   J =
   \left[
          \begin{array}{ccc}
             \frac{\partial f_1}{\partial x}
             & \frac{\partial f_1}{\partial y}
             & \frac{\partial f_1}{\partial s} \\
             \frac{\partial f_2}{\partial x}
             & \frac{\partial f_2}{\partial y}
             & \frac{\partial f_2}{\partial s}
          \end{array}
   \right]

which for the polynomials in our problem becomes

.. math::

   J =
   \left[
      \begin{array}{ccc}
         2 (x - a) & 2(y-b) & 0 \\
                -s &    1   & -x
      \end{array}
   \right]

where the :math:`-` on the second row appears 
from the equation :math:`-sx + y = 0`.

As we look for a singular solution,
the rank of the Jacobian :math:`J` must be one, not two.
It suffices here to consider the first 2-by-2 minor and 
requiring that the first two columns are linearly dependent
can be expressed by

.. math::

   \left\{
      \begin{array}{rlcrlcl}
         2(x-a) & L_1 & + & 2(y-b) & L_2 & = & 0 \\
             -s & L_1 & + &        & L_2 & = & 0 \\
      c_0 + c_1 & L_1 & + &   c_2  & L_2 & = & 0
      \end{array}
   \right.

where :math:`L_1` and :math:`L_2` are two new variables 
and :math:`c_0`, :math:`c_1`, and :math:`c_2` are random constants.

::

    def random_complex():\n",
        """
        Returns a random complex number on the unit circle.
        """
        theta = uniform(0, 2*pi)\n",
        return complex(cos(theta), sin(theta))"

::

    def random_hyperplane(vars):
        """
        Returns a linear equation in the variables
        in the list vars, with random complex coefficients.
        """
        cf0 = str(random_complex())
        tf0 = cf0.replace('j', '*i')
        result = tf0
        for var in vars:
            cff = str(random_complex())
            tcf = cff.replace('j', '*i')
            result = result + '+' + tcf + '*' + var
        return result + ';'

::

    def jacobian(a, b):
        """
        Returns the equations which define the points
        where the Jacobian matrix is singular,
        as a random linear combination of the columns.
        Random complex coefficients are generated to
        scale the multiplier variables.
        """
        eq1 = '2*(x-%.15e)*L1 + 2*(y-%.15e)*L2;' % (a, b)
        eq2 = '-s*L1 + L2;'
        eq3 = random_hyperplane(['L1', 'L2'])
        return [eq1, eq2, eq3]

Then here is an example for the center :math:`(3, 2)`.

::

    jacpols = jacobian(3, 2)
    for pol in jacpols:
        print(pol)

which produces the polynomials

::

    2*(x-3.000000000000000e+00)*L1 + 2*(y-2.000000000000000e+00)*L2;
    -s*L1 + L2;
    (0.7467754449282183-0.6650762624333104*i)+(-0.05113246140388107-0.9986918801065625*i)*L1+(0.8583591351915719-0.5130493105279226*i)*L2;

Let us now make a witness set, using the auxiliary function below.

::

    def witset(pols, verbose=True):
        """
        We have three equations in pols in five variables:
        x, y, s, L1, and L2.
        Therefore we expect a two dimensional set.
        """
        embpols = double_embed(5, 2, pols)
        if verbose:
            print('the embedded system :')
            for pol in embpols:
                print(pol)
        embsols = solve(embpols)
        if verbose:
            print('the witness points :')
            for sol in embsols:
                print(sol)
        return (embpols, embsols)

As we have three equations in five variables 
:math:`(x, y, s, L_1, L_3)`, the solution set is two dimensional.
The construction of this witness set 
will compute the degree of this two dimensional solution set

::

    def singular_locus_set():
        """
        Generates a witness set for the singular locus
        of the algebraic set of a fixed circle intersected
        with a one parameter family of lines.
        """
        syst = jacobian(3, 2)
        for pol in syst:
            print(pol)
        (embsyst, embsols) = witset(syst, False)
        print('the polynomials in the witness set :')
        for pol in embsyst:
            print(pol)
        print('the solutions :')
        for sol in embsols:
            print(sol)
        print('degree of the singular locus set :', len(embsols))
        return (embsyst, embsols)

The output of

::

    witset2 = singular_locus_set()

is

::

    2*(x-3.000000000000000e+00)*L1 + 2*(y-2.000000000000000e+00)*L2;
    -s*L1 + L2;
    (-0.6742208883330486-0.7385297514219686*i)+(-0.9307410468108358-0.36567896272751255*i)*L1+(0.4913219172246197+0.8709780557825346*i)*L2;
    the polynomials in the witness set :
     + 2*x*L1 + 2*y*L2 - 6*L1 - 4*L2 + (4.93184439043922E-01 + 8.69924772083731E-01*i)*zz1 + (4.77962297913147E-01 + 8.78380351427321E-01*i)*zz2;
     - L1*s + L2 + (-3.50715092194063E-03-9.99993849927294E-01*i)*zz1 + (6.78305496062243E-01-7.34780003818663E-01*i)*zz2;
     + (-9.30741046810836E-01-3.65678962727513E-01*i)*L1 + (4.91321917224620E-01 + 8.70978055782535E-01*i)*L2 + (-1.73576955750725E-01 + 9.84820308702207E-01*i)*zz1 + (-8.69131200089510E-01 + 4.94581597950195E-01*i)*zz2+(-6.74220888333049E-01-7.38529751421969E-01*i);
    zz1;
    zz2;
     + (-1.94485007680957E-01-2.95255113058755E-01*i)*x + (-9.05619807296798E-02-3.41757995731361E-01*i)*L1 + (2.35555493200339E-01-2.63654337387317E-01*i)*y + (2.27470946523190E-01 + 2.70660245488406E-01*i)*L2 + (3.53287256092446E-01 + 1.37154906099080E-02*i)*s + (-2.50508068634215E-01-2.49490896726024E-01*i)*zz1 + (3.37269024813124E-01 + 1.06064154649929E-01*i)*zz2+(-3.53553390593274E-01-1.38777878078145E-17*i);
     + (-4.40566415427727E-01 + 1.40592718839002E-02*i)*x + (-2.45676058279196E-01-2.88340078644854E-01*i)*L1 + (2.53550532539115E-01 + 6.86270383050135E-02*i)*y + (-2.85183161836763E-01-2.03395928752962E-01*i)*L2 + (2.39045196476971E-01 + 1.56472612484202E-01*i)*s + (-9.34949306113738E-02 + 4.61532488129392E-01*i)*zz1 + (-1.98366575492475E-01 + 1.97233362690695E-01*i)*zz2+(2.06250450361319E-01-2.15268649297659E-01*i);
    the solutions :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  1.89091316713556E+00  -2.07966565089632E+00
     L1 : -1.55981268759518E-01   1.28195421958451E+00
     y :  9.29382853428914E-01   1.21952076881767E+00
     L2 :  1.66239264497059E+00   8.68576316936316E-01
     zz1 :  1.23844424223158E-33   2.14768913973412E-32
     zz2 : -3.68597450920629E-33  -1.56050619802567E-32
     s :  5.12174926048101E-01  -1.35908311946356E+00
    == err :  1.287E-15 = rco :  3.154E-02 = res :  1.332E-15 =
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  1.30693118366194E+00  -1.84306554595479E-01
     L1 : -7.81141131311130E-01   2.88897684810126E-01
     y :  9.23707475868645E-01   1.58931255339032E+00
     L2 :  5.50791044105935E-01   4.92640130916949E-01
     zz1 : -5.01359801141478E-33  -1.18066233289204E-32
     zz2 :  0.00000000000000E+00   0.00000000000000E+00
     s : -4.15087884109074E-01  -7.84183593815662E-01
    == err :  6.146E-16 = rco :  3.225E-02 = res :  8.882E-16 =
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  2.56295575999062E-01   5.00768723718571E-01
     L1 : -1.54744896404770E-01  -3.07935271191708E-01
     y :  1.82176565498780E+00  -1.26116846711605E+00
     L2 :  6.60149130688637E-01  -3.65627094004125E-01
     zz1 : -8.24676473579955E-32  -3.81646685718875E-32
     zz2 :  0.00000000000000E+00   0.00000000000000E+00
     s :  8.78568617765025E-02   2.18794206021058E+00
    == err :  8.086E-16 = rco :  2.783E-02 = res :  5.135E-16 =
    degree of the singular locus set : 3

Indeed, the singular locus has degree three.

extending with slack variables
------------------------------

The first witness set of the circle intersected 
with a line is a cubic curve.  
In order to intersect this set with the singular locus, 
we have to add :math:`L_1` and :math:`L_2`, 
increasing the dimension by two, using one as the values 
of the witness points for :math:`L_1` and :math:`L_2`.

The function below adds to the solutions the values 
for the :math:`L_1` and :math:`L_2`, 
and adds two additional slack variables ``zz2`` and ``zz3``.

::

    def extend_solutions(sols):
        """
        To each solution in sols, adds L1 and L2 with values 1,
        and zz2 and zz3 with values zero.
        """
        result = []
        for sol in sols:
            (vars, vals) = coordinates(sol)
            vars = vars + ['L1', 'L2', 'zz2', 'zz3']
            vals = vals + [1, 1, 0, 0]
            extsol = make_solution(vars, vals)
            result.append(extsol)
        return result

::

    def extend(pols, sols, verbose=True):
        """
        Extends the witness set with two free variables
        L1 and L2, addition two linear equations,
        and two slack variables zz2 and zz3.
        """
        vars = ['zz2', 'zz3']
        eq1 = 'zz2;'
        eq2 = 'zz3;'
        eq3 = 'L1 - 1;'
        eq4 = 'L2 - 1;'
        extpols = pols[:-1] + [eq1, eq2, eq3, eq4, pols[-1]]
        extsols = extend_solutions(sols)
        if verbose:
            print('the extended polynomials :')
            for pol in extpols:
                print(pol)
            print('the extended solutions :')
            for sol in extsols:
                print(sol)
        return (extpols, extsols)

The output of 

::

    witset1 = extend(witpols, witsols)

is

::

    the extended polynomials :
     + x^2 + y^2 - 6*x - 4*y + (-2.25888608394754E-01 + 9.74153138165392E-01*i)*zz1 + 12;
     - x*s + y + (4.03414981775719E-02 + 9.99185950424038E-01*i)*zz1;
    zz1;
    zz2;
    zz3;
    L1 - 1;
    L2 - 1;
     + (-1.19285340138354E-01 + 9.92860014114818E-01*i)*x + (-8.44394833617596E-01-5.35721350106483E-01*i)*y + (2.84736834956935E-01-9.58605724382401E-01*i)*s + (-8.64325704104395E-01-5.02932477798403E-01*i)*zz1+(4.42123925817800E-01-8.96953975530214E-01*i);
    the extended solutions :
    t : 0.000000000000000E+00 0.000000000000000E+00
    m : 1
    the solution for t :
     x : 4.719921861196670E+00  -2.064791911548440E+00
     y : 4.210489631645690E+00  1.606558427894690E+00
     zz1 : 1.471246300030830E-32  -4.656335321781320E-32
     s : 6.237879407983900E-01  6.132624241882660E-01
     L1 : 1.000000000000000E+00  0.0
     L2 : 1.000000000000000E+00  0.0
     zz2 : 0.000000000000000E+00  0.0
     zz3 : 0.000000000000000E+00  0.0
    == err : 0.000E+00 = rco : 1.000E+00 = res : 0.000E+00 =
    t : 0.000000000000000E+00 0.000000000000000E+00
    m : 1
    the solution for t :
     x : 2.227863046889190E+00  -3.253287772536380E-01
     y : 1.217286103395750E+00  3.209325552001770E-01
     zz1 : 8.778759246831849E-33  -8.433452646430390E-33
     s : 5.143872141853040E-01  2.191685522625720E-01
     L1 : 1.000000000000000E+00  0.0
     L2 : 1.000000000000000E+00  0.0
     zz2 : 0.000000000000000E+00  0.0
     zz3 : 0.000000000000000E+00  0.0
    == err : 0.000E+00 = rco : 1.000E+00 = res : 0.000E+00 =
    t : 0.000000000000000E+00 0.000000000000000E+00
    m : 1
    the solution for t :
     x : -3.194195080011810E-01  -6.360872929189700E-01
     y : 1.334206001109740E+00  3.171312106186350E+00
     zz1 : 0.000000000000000E+00  0.000000000000000E+00
     s : -4.822798620424600E+00  -3.243107726117300E-01
     L1 : 1.000000000000000E+00  0.0
     L2 : 1.000000000000000E+00  0.0
     zz2 : 0.000000000000000E+00  0.0
     zz3 : 0.000000000000000E+00  0.0
    == err : 0.000E+00 = rco : 1.000E+00 = res : 0.000E+00 =
    t : 0.000000000000000E+00 0.000000000000000E+00
    m : 1
    the solution for t :
     x : -1.027401009062140E+00  -1.336145604993470E+00
     y : 3.374637845889450E+00  -3.914626804358620E+00
     zz1 : 0.000000000000000E+00  0.000000000000000E+00
     s : 6.207341379163620E-01  3.002951707170900E+00
     L1 : 1.000000000000000E+00  0.0
     L2 : 1.000000000000000E+00  0.0
     zz2 : 0.000000000000000E+00  0.0
     zz3 : 0.000000000000000E+00  0.0
    == err : 0.000E+00 = rco : 1.000E+00 = res : 0.000E+00 =

intersecting two witness sets
-----------------------------

When intersecting two witness sets, the symbols must line up.  
The function below ensures that the order of symbols is 
:math:`(x, y, s, L1, L2)` by adding a zero polynomial to a given polynomial.

::

    def insert_symbols(pol):
        """
        To the string pol, adds the sequence of symbols.
        """
        q = pol.lstrip()
        if q[0] == '+' or q[0] == '-':
            smb = 'x - x + y - y + s + L1 - L1 + L2 - L2 - s '
        else:
            smb = 'x - x + y - y + s + L1 - L1 + L2 - L2 - s + '
        return smb + pol

::

    def intersect(dim, w1d, w2d, ws1, ws2, verbose=True):
        """
        Applies the diagonal homotopy to intersect two witness sets
        w1 and w2 of dimensions w1d and w2d in a space of dimension dim.
        """
        w1eqs, w1sols = ws1
        w2eqs, w2sols = ws2
        nw1eq0 = insert_symbols(w1eqs[0])
        nw2eq0 = insert_symbols(w2eqs[0])
        nw1eqs = [nw1eq0] + w1eqs[1:]
        nw2eqs = [nw2eq0] + w2eqs[1:]
        if verbose:
            print('number of equations in first witness set :', len(w1eqs))
            print('number of equations in second witness set :', len(w2eqs))
        result = double_diagonal_solve(dim, w1d, nw1eqs, w1sols, w2d, nw2eqs, w2sols, \
            vrblvl=int(verbose))
        (eqs, sols) = result
        if verbose:
            print('the equations :')
            for pol in eqs:
                print(pol)
        return result

In a five dimensional space, 
we intersect a three dimensional set with a two dimensional one.

::

    (eqs, sols) = intersect(5, 3, 2, witset1, witset2, verbose=False)
    print('the solutions :')
    for (idx, sol) in enumerate(sols):
        print('Solution', idx+1, ':')
        print(sol)

The output of the above code is
  
::

    the solutions :
    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  3.30216947925196E+00   1.91366021429925E-16
     y :  1.04674578112206E+00   6.06605980827076E-17
     s :  3.16987298107781E-01  -4.27079988595092E-32
     L1 : -3.06590066624447E-01  -2.66629188670126E-01
     L2 : -9.67199848241872E-01  -8.41135245045269E-01
    == err :  2.341E-15 = rco :  1.884E-02 = res :  1.332E-15 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  2.23629205920958E+00   2.45132855242376E-17
     y :  2.64556191118564E+00   2.89995281402837E-17
     s :  1.18301270189222E+00   1.27877769258935E-33
     L1 :  5.35439409455237E-01  -1.48149502880241E+00
     L2 :  4.52606644543044E-01  -1.25230695024049E+00
    == err :  6.349E-15 = rco :  2.003E-02 = res :  2.220E-16 =


Observe that we computed two real solutions, 
with the values of the intersection points 
and the slopes of the tangent lines.
The two real solutions are regular solutions.

plotting the tangent lines
--------------------------

From the solution we extract the (real) coordinates 
of the intersection points and the slopes of the tangent lines,
as codified in the following function.

::

    def coordinates_and_slopes(sol):
        """
        Given a solution, return the 3-tuple with the x and y coordinates
        of the tangent point and the slope of the tangent line.
        The real parts of the coordinates are selected.
        """
        sdc = strsol2dict(sol)
        return (sdc['x'].real, sdc['y'].real, sdc['s'].real)

We apply the function to the first solution:

::

    sol1 = sols[0]
    coordinates_and_slopes(sol1)

with the output

::

    (3.30216947925196, 1.04674578112206, 0.316987298107781)

and we do this also for the second solution

::

    sol2 = sols[1]
    coordinates_and_slopes(sol2)

which gives the output

::

    (2.23629205920958, 2.64556191118564, 1.18301270189222)"

The code below
plots the circle with radius 1, centered at (3, 2),
along with its two lines tangent through (0, 0),
with the plot shown in :numref:`touchcirclefig2`.
The two solutions ``sol1`` and ``sol2`` define 3-tuples of
x and y coordinates and a slope.

::

    fig2 = plt.figure()
    (xp1, yp1, s1) = coordinates_and_slopes(sol1)
    (xp2, yp2, s2) = coordinates_and_slopes(sol2)
    axs = fig2.gca()
    center = (3, 2)
    radius = 1
    circle = plt.Circle(center, radius, edgecolor='blue', facecolor='none')
    axs.add_artist(circle)
    y1 = 5*s1 # first tangent line
    y2 = 5*s2 # second tangent line
    plt.plot([0, 5], [0, y1], 'r')
    plt.plot([0, 5], [0, y2], 'r')
    plt.plot([xp1, xp2], [yp1, yp2], 'go')
    plt.axis([0, 5, 0, 4])
    plt.show()

.. _touchcirclefig2:

.. figure:: ./touchcirclefig2.png
   :align: center
    
   Two lines tangent to a circle.
