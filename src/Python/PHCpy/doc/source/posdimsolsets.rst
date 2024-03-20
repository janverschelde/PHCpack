Positive Dimensional Solution Sets
==================================

A more complete solver returns a numerical irreducible decomposition,
with not only the isolated solutions, but the dimension and the degrees 
of all positive dimensional solution sets, for all dimensions.

witness sets
------------

A witness set of a pure dimensional solution set consists of

1. the original polynomial system augmented with as many random hyperplanes
   as the dimension of the solution set; and

2. as many solutions to the augmented system
   as the degree of the solution set.

Because of the random coefficients in the hyperplanes,
the solutions are generic points.

Via an embedding of the augmented system, the computation of generic points
is reduced to the computation of isolated solutions, which can be handled
well by the blackbox solver.

::

    from phcpy.sets import double_embed
    from phcpy.solver import solve

Let us make a witness set of the twisted cubic.

::

    twisted = ['x^2 - y;', 'x^3 - z;']

The twisted cubic is the space curve with parametric 
form :math:`(t, t^2, t^3)`.

The dimension is ``1`` and we are in dimension ``3``.
The embedded system is constructed as:

::

    embtwist = double_embed(3, 1, twisted)
    for pol in embtwist:
        print(pol)

and the embedded polynomials are

::

    x^2 - y + (5.62891189304868E-01-8.26531009099448E-01*i)*zz1;
    x^3 - z + (-9.84048066692442E-01-1.77902789294791E-01*i)*zz1;
    zz1;
     + (8.39559929055672E-01-5.43267084889224E-01*i)*x + (-1.14034198908312E-01-9.93476824832537E-01*i)*y + (-9.45117489397468E-01-3.26730670790218E-01*i)*z + (4.57472148097901E-01-8.89223950259265E-01*i)*zz1+(-9.21723254199187E-01-3.87848221174806E-01*i);

The symbol ``zz1`` is used for the slack variable.

::

    sols = solve(embtwist)
    for (idx, sol) in enumerate(sols):
        print('Solution', idx+1, ':')
        print(sol)

and the generic points on the twisted cubic are

::

    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  6.31159601615322E-01  -1.19854989224432E+00
     y : -1.03815940148766E+00  -1.51295254501003E+00
     zz1 : -4.04959151575673E-33   1.16188546226650E-32
     z : -2.46859338404870E+00   2.89371313214052E-01
    == err :  6.130E-16 = rco :  2.728E-02 = res :  3.331E-16 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x : -1.34539355262783E+00  -1.68943360753041E-01
     y :  1.78154195231000E+00   4.54590616632837E-01
     zz1 :  0.00000000000000E+00   0.00000000000000E+00
     z : -2.32007498983311E+00  -9.12582969448712E-01
    == err :  7.943E-16 = rco :  3.560E-02 = res :  1.166E-15 =
    Solution 3 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  2.81858885842759E-01   4.65799400839404E-01
     y : -1.37524650293826E-01   2.62579400293639E-01
     zz1 :  5.64466889738308E-34  -9.64221135309633E-34
     z : -1.61071872037280E-01   9.95143750451212E-03
    == err :  5.737E-17 = rco :  9.995E-02 = res :  4.163E-17 =

As the degree of the twisted cubic is three,
we have three generic points.
Observe the tiny values for the slack variable ``zz1``,
as the generic points must satisfy the added random
hyperplane of the augmented system.

homotopy membership test
------------------------

A witness set can be used to decide whether any point belongs 
to the algebraic set represented by the witness set.  
The homotopy membership test is illustrated via the solution set 
of the cyclic 4-roots system.

::

    from phcpy.families import cyclic
    c4 = cyclic(4)
    for pol in c4:
        print(pol)

Here are the four polynomials:

::

    x0 + x1 + x2 + x3;
    x0*x1 + x1*x2 + x2*x3 + x3*x0;
    x0*x1*x2 + x1*x2*x3 + x2*x3*x0 + x3*x0*x1;
    x0*x1*x2*x3 - 1;

Although we have four equations and four unknowns
(we expect thus only isolated solutions), we know that
for this particular system,
the solution set is pure dimensional, of dimension 1.

To construct a witness set, we import the following functions:

::

    from phcpy.sets import double_embed
    from phcpy.solver import solve
    from phcpy.solutions import filter_zero_coordinates as filter

and then execute

::
  
    c4e1 = double_embed(4, 1, c4)
    sols = solve(c4e1)
    genpts = filter(sols, 'zz1', 1.0e-8, 'select')
    print('generic points on the cyclic 4-roots set :')
    for (idx, sol) in enumerate(genpts):
        print('Solution', idx+1, ':')
        print(sol)

to see the generic points on the solution curve
of the cyclic 4-roots:

::

    generic points on the cyclic 4-roots set :
    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x0 :  9.65599349076935E-01  -1.16989010460731E+00
     x1 : -4.19638798339057E-01  -5.08421301397389E-01
     x2 : -9.65599349076935E-01   1.16989010460732E+00
     x3 :  4.19638798339057E-01   5.08421301397389E-01
     zz1 :  1.76873803944398E-16  -1.05541954650188E-16
    == err :  1.859E-15 = rco :  4.629E-02 = res :  9.649E-16 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x0 : -4.72839263499989E-01  -1.41379607496008E+00
     x1 : -2.12761000919484E-01   6.36158397206689E-01
     x2 :  4.72839263499989E-01   1.41379607496008E+00
     x3 :  2.12761000919484E-01  -6.36158397206689E-01
     zz1 : -8.45579970922059E-17   3.34398025784039E-17
    == err :  1.268E-15 = rco :  5.880E-02 = res :  5.892E-16 =
    Solution 3 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x0 : -7.78676715642733E-01   2.78291443186817E-01
     x1 :  1.13877660574983E+00   4.06987622353548E-01
     x2 :  7.78676715642734E-01  -2.78291443186816E-01
     x3 : -1.13877660574983E+00  -4.06987622353548E-01
     zz1 : -8.30236390890514E-17  -4.21955685140459E-17
    == err :  1.296E-15 = rco :  1.051E-01 = res :  1.496E-15 =
    Solution 4 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x0 :  6.59761896934191E-01   5.22197413539580E-01
     x1 :  9.31898808330260E-01  -7.37592076250530E-01
     x2 : -6.59761896934191E-01  -5.22197413539580E-01
     x3 : -9.31898808330260E-01   7.37592076250530E-01
     zz1 : -6.83231299973127E-17   7.95281480444320E-17
    == err :  1.022E-15 = rco :  6.324E-02 = res :  9.147E-16 =

For the membership test in double precision, we use the function:

::

    from phcpy.sets import double_membertest

Consider two test points ``pt0`` and ``pt1``.

The first test

::

    pt0 = [1, 0, -1, 0, 1, 0, -1, 0]
    ismbr = double_membertest(c4e1, sols, 1, pt0)
    print('Is', pt0, 'a member?', ismbr)

gives as output

::

    Is [1, 0, -1, 0, 1, 0, -1, 0] a member? False

and the second test

::

    pt1 = [1, 0, 1, 0, -1, 0, -1, 0]
    ismbr = double_membertest(c4e1, sols, 1, pt1)
    print('Is', pt1, 'a member?', ismbr)

yields

::

     Is [1, 0, 1, 0, -1, 0, -1, 0] a member? True

monodromy breakup
-----------------

The factorization into irreducible components is illustrated 
on a cubic curve.

::

    cubic = '(x+1)*(x^2 + y^2 + 1);'

The input to the factorization function is a witness set.

::

    from phcpy.sets import double_hypersurface_set
    from phcpy.factor import double_monodromy_breakup, write_factorization

The construction of the witness set happens via

::

    (wit, pts) = double_hypersurface_set(2, cubic)
    for pol in wit:
        print(pol)\n",
    print('number of witness points :', len(pts))

and the output is

::

    x^3 + x*y^2 + x^2 + y^2 + x + (5.56101869358167E-01-8.31114138308544E-01*i)*zz1 + 1;
    zz1;
     + (9.85343874390340E-01 + 1.70579744405467E-01*i)*x + y + zz1+(-1.26905195457699E+00 + 1.50366546205483E+00*i);
    number of witness points : 3

To see the witness points, execute

::

    for (idx, sol) in enumerate(pts):
        print('Solution', idx+1, ':')
        print(sol)

which then prints

::

    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x : -4.30394583533542E-02  -1.45862430717631E+00
     y :  1.06264885972081E+00  -5.90772761365393E-02
     zz1 :  0.00000000000000E+00   0.00000000000000E+00
    == err :  0.000E+00 = rco :  1.000E+00 = res :  0.000E+00 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x :  1.33096744731273E+00  -6.74068593392326E-02
     y : -5.39069114828172E-02  -1.66428261308763E+00
     zz1 :  0.00000000000000E+00   0.00000000000000E+00
    == err :  0.000E+00 = rco :  1.000E+00 = res :  0.000E+00 =
    Solution 3 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x : -1.00000000000000E+00   1.11022302462516E-16
     y :  2.25439582896733E+00  -1.33308571764937E+00
     zz1 :  0.00000000000000E+00   0.00000000000000E+00
    == err :  0.000E+00 = rco :  1.000E+00 = res :  0.000E+00 =

The factorization is then computed via

::

    deco = double_monodromy_breakup(wit, pts, dim=1)

To see the grouping of the generic points according 
to the irreducible factors, we do:

::

    write_factorization(deco)

which prints

::

      factor 1 : ([1, 2], 3.683680191092806e-15)
      factor 2 : ([3], 9.742207041085749e-15)

cascade of homotopies
---------------------

A cascade of homotopies computes generic points on all positive
components of the solution set, for all dimensions.

::

    pol1 = '(x^2 + y^2 + z^2 - 1)*(y - x^2)*(x - 0.5);'
    pol2 = '(x^2 + y^2 + z^2 - 1)*(z - x^3)*(y - 0.5);'
    pol3 = '(x^2 + y^2 + z^2 - 1)*(z - x*y)*(z - 0.5);'
    pols = [pol1, pol2, pol3]"

The solution set of ``pols`` contains the sphere, the twisted cubic,
some lines, and an isolated point.

To run a cascade of homotopies, import the following functions:

::

    from phcpy.cascades import double_top_cascade, double_cascade_filter

We start at the top dimension of the solution set:

::

    (embpols, sols0, sols1) = double_top_cascade(3, 2, pols)
    print('at dimension 2, degree :', len(sols0))"

to find two witness points on the sphere:

::

    at dimension 2, degree : 2

and then continue to the one dimensional solution sets:

::

    (wp1, ws0, ws1) = double_cascade_filter(2, embpols, sols1, tol=1.0e-8)
    print('at dimension 1, candidate generic points :', len(ws0))

to obtain:

::

    at dimension 1, candidate generic points : 9

and then continue to the isolated solutions:

::

    (wp0, ws0, ws1) = double_cascade_filter(1, wp1, ws1, tol=1.0e-8)
    print('candidate isolated points :', len(ws0))

to find

::

    candidate isolated points : 24

The output of the cascade needs further processing.
The ``solve`` in the next section does all steps.

numerical irreducible decomposition
-----------------------------------

The computation of a numerical irreducible decomposition starts 
by solving the top dimensional system in an embedding and then 
applies a cascade of homotopies to compute candidate generic points on
each positive dimensional solution set, ending at the isolated solutions.
After each step in the cascade, the candidate generic points are filtered
as some may lie on higher dimensional sets, and the generic points are 
grouped according to their irreducible factors.

::

    from phcpy.decomposition import solve, write_decomposition

The second blackbox solver is illustrated on the following example:

::

    pol0 = '(x1-1)*(x1-2)*(x1-3)*(x1-4);'
    pol1 = '(x1-1)*(x2-1)*(x2-2)*(x2-3);'
    pol2 = '(x1-1)*(x1-2)*(x3-1)*(x3-2);'
    pol3 = '(x1-1)*(x2-1)*(x3-1)*(x4-1);'
    pols = [pol0, pol1, pol2, pol3]

For this small example, we ask the ``solve`` to be silent:

::

    deco = solve(pols, verbose=False)

and then write the decomposition in ``deco`` as

::

    write_decomposition(deco)

The output is rather extensive ...

::

    set of dimension 0 has degree 4
    the polynomials :
     + x1^4 - 10*x1^3 + 35*x1^2 - 50*x1 + 24;
    x1*x2^3 - 6*x1*x2^2 - x2^3 + 11*x1*x2 + 6*x2^2 - 6*x1 - 11*x2 + 6;
    x1^2*x3^2 - 3*x1^2*x3 - 3*x1*x3^2 + 2*x1^2 + 9*x1*x3 + 2*x3^2 - 6*x1 - 6*x3 + 4;
    x1*x2*x3*x4 - x1*x2*x3 - x1*x2*x4 - x1*x3*x4 - x2*x3*x4
    + x1*x2 + x1*x3 + x1*x4 + x2*x3 + x2*x4 + x3*x4 - x1 - x2 - x3 - x4 + 1;
    the generic points :
    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x1 :  4.00000000000000E+00   0.00000000000000E+00
     x2 :  3.00000000000000E+00   0.00000000000000E+00
     x3 :  2.00000000000000E+00   0.00000000000000E+00
     x4 :  1.00000000000000E+00   0.00000000000000E+00
    == err :  4.956E-26 = rco :  1.359E-01 = res :  0.000E+00 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x1 :  4.00000000000000E+00   3.85357077689325E-45
     x2 :  2.00000000000000E+00  -2.94004118110142E-50
     x3 :  2.00000000000000E+00   1.64214663788065E-46
     x4 :  1.00000000000000E+00   1.77899219103737E-47
    == err :  7.105E-15 = rco :  1.191E-01 = res :  2.416E-44 =
    Solution 3 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x1 :  3.00000000000000E+00   4.27642353614751E-50
     x2 :  3.00000000000000E+00   0.00000000000000E+00
     x3 :  2.00000000000000E+00   3.20731765211063E-50
     x4 :  1.00000000000000E+00   0.00000000000000E+00
    == err :  1.204E-34 = rco :  1.043E-01 = res :  1.497E-49 =
    Solution 4 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x1 :  2.99999999999999E+00   3.50324616081204E-45
     x2 :  2.00000000000000E+00   2.83543986140725E-45
     x3 :  2.00000000000000E+00  -7.18165462966469E-45
     x4 :  1.00000000000000E+00   6.30584308946168E-45
    == err :  9.770E-15 = rco :  9.066E-02 = res :  1.066E-14 =
    set of dimension 1 has degree 12
    the polynomials :
     + x1^4 - 10*x1^3 + 35*x1^2 - 50*x1
    + (9.98263256449285E-01-5.89107021114901E-02*i)*zz1 + 24;
    x1*x2^3 - 6*x1*x2^2 - x2^3 + 11*x1*x2 + 6*x2^2 - 6*x1 - 11*x2
    + (8.25826300922299E-02-9.96584220829855E-01*i)*zz1 + 6;
    x1^2*x3^2 - 3*x1^2*x3 - 3*x1*x3^2 + 2*x1^2 + 9*x1*x3 + 2*x3^2 - 6*x1
    - 6*x3 + (1.55033404443160E-01-9.87909228374127E-01*i)*zz1 + 4;
    x1*x2*x3*x4 - x1*x2*x3 - x1*x2*x4 - x1*x3*x4 - x2*x3*x4
    + x1*x2 + x1*x3 + x1*x4 + x2*x3 + x2*x4 + x3*x4 - x1 - x2 - x3 - x4
    + (-8.97158012682676E-01-4.41709746642829E-01*i)*zz1 + 1;
    + (-3.53404925407865E-01-1.02449352102111E-02*i)*x1
    + (-3.23378971672498E-01-1.42919700111768E-01*i)*x2
    + (2.20216662912369E-01-2.76594687902244E-01*i)*x3
    + (8.35397253138817E-02-3.43542012415485E-01*i)*x4
    + (-2.67791727313329E-01-2.30841050904175E-01*i)*zz1 - 3.53553390593274E-01;
    the generic points :
    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x1 :  2.00000000000000E+00   1.66274442058684E-16
     x2 :  2.00000000000000E+00   1.93730374623830E-15
     x3 :  1.42230996999871E+00   4.73749193868512E+00
     x4 :  1.00000000000000E+00   3.09283683575416E-16
     zz1 : -1.98123724512471E-15  -4.50046576766827E-16
    == err :  2.605E-14 = rco :  1.263E-02 = res :  6.855E-14 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x1 :  4.00000000000001E+00   1.93485385394018E-17
     x2 :  1.00000000000000E+00  -5.10382551856806E-18
     x3 :  2.00000000000000E+00  -3.61201796027518E-18
     x4 : -9.22964074590540E-01   5.02769047428598E+00
     zz1 : -4.05629881527961E-17  -1.18686954150410E-16
    == err :  2.775E-14 = rco :  1.246E-02 = res :  7.817E-14 =
    Solution 3 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x1 :  2.99999999999999E+00   3.50163623105219E-16
     x2 :  3.00000000000000E+00   8.96756225629632E-19
     x3 :  1.00000000000000E+00   2.40219551951397E-17
     x4 : -5.76987365375915E-01   6.43848410091975E+00
     zz1 :  6.20367476284742E-17   7.05206637649815E-16
    == err :  3.038E-14 = rco :  7.888E-03 = res :  3.206E-14 =
    Solution 4 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x1 :  3.00000000000000E+00   1.29780826664201E-16
     x2 :  1.00000000000000E+00   9.36913548588266E-18
     x3 :  2.00000000000000E+00   8.93428893566725E-18
     x4 : -1.13099435246225E+00   4.04956808752212E+00
     zz1 :  5.94418462635302E-17   2.63521082767564E-16
    == err :  2.920E-15 = rco :  2.162E-02 = res :  2.671E-16 =
    Solution 5 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :\n",
     x1 :  2.00000000000000E+00  -3.65568752932322E-16
     x2 :  3.00000000000000E+00   2.82574614718779E-15
     x3 :  1.67577083520075E+00   5.70483758002075E+00
     x4 :  1.00000000000000E+00  -3.72113707233408E-16
     zz1 :  5.75972032456864E-15   1.07230899989073E-15
    == err :  2.494E-14 = rco :  8.760E-03 = res :  2.468E-14 =
    Solution 6 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x1 :  2.00000000000000E+00   1.76757984029786E-16
     x2 :  3.00000000000000E+00  -1.31917536806969E-17
     x3 :  1.00000000000000E+00  -4.99155351560570E-17
     x4 : -7.85017643247620E-01   5.46036171415591E+00
     zz1 : -5.60935500844177E-17  -3.57441262286022E-16
    == err :  2.215E-14 = rco :  1.929E-02 = res :  1.493E-14 =
    Solution 7 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x1 :  4.00000000000000E+00  -2.92505382506395E-18
     x2 :  2.00000000000000E+00  -1.74547448231340E-18
     x3 :  1.00000000000000E+00  -6.45889285731916E-19
     x4 : -1.92285640108937E-01   6.43233660615961E+00
     zz1 :  6.74420113461971E-18   1.79788532318047E-17
    == err :  1.018E-14 = rco :  1.020E-02 = res :  6.776E-14 =
    Solution 8 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x1 :  2.99999999999999E+00   1.21726975202562E-16
     x2 :  1.00000000000000E+00  -5.19731533632551E-18
     x3 :  1.00000000000000E+00   1.92249486589598E-17
     x4 : -2.23644470585381E-01   4.46994433787176E+00
     zz1 : -6.54642209058018E-19   2.43838870558812E-16
    == err :  2.701E-14 = rco :  3.711E-02 = res :  1.071E-14 =
    Solution 9 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x1 :  3.00000000000001E+00   7.59453204499381E-16
     x2 :  2.00000000000000E+00  -3.13045199988346E-17
     x3 :  1.00000000000000E+00   2.50434205156814E-17
     x4 : -4.00315917980645E-01   5.45421421939577E+00
     zz1 :  1.89836152203951E-16   1.53275178679189E-15
    == err :  6.146E-14 = rco :  1.385E-02 = res :  3.244E-14 =
    Solution 10 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x1 :  4.00000000000000E+00  -1.25836606607375E-17
     x2 :  3.00000000000000E+00   1.36281411894719E-17
     x3 :  1.00000000000000E+00  -1.25232978409423E-17
     x4 : -3.68957087504202E-01   7.41660648768361E+00
     zz1 :  8.87505203002171E-17   8.08707712184250E-17
    == err :  1.058E-14 = rco :  5.769E-03 = res :  6.754E-14 =
    Solution 11 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x1 :  4.00000000000000E+00  -4.37056700414583E-17
     x2 :  1.00000000000000E+00   2.28708329195796E-17
     x3 :  1.00000000000000E+00  -1.93529333821221E-17
     x4 : -1.56141927136710E-02   5.44806672463562E+00
     zz1 :  1.60246975171774E-16   2.72146931495479E-16
    == err :  8.019E-16 = rco :  5.425E-02 = res :  2.567E-16 =
    Solution 12 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x1 :  2.00000000000000E+00   3.54994546241256E-16
     x2 :  2.00000000000000E+00   5.29877124103883E-17
     x3 :  1.00000000000000E+00   8.46409944528661E-17
     x4 : -6.08346195852355E-01   4.47609183263191E+00
     zz1 : -1.12656321968022E-16  -7.17872515969298E-16
    == err :  6.181E-15 = rco :  3.233E-02 = res :  1.200E-14 =
    set of dimension 2 has degree 1
    the polynomials :
    x1^4 - 10*x1^3 + 35*x1^2 - 50*x1
    + (9.98263256449285E-01-5.89107021114901E-02*i)*zz1
    + (-4.67731979750726E-01 + 8.83870349722439E-01*i)*zz2 + 24;
    x1*x2^3 - 6*x1*x2^2 - x2^3 + 11*x1*x2 + 6*x2^2 - 6*x1 - 11*x2
    + (8.25826300922299E-02-9.96584220829855E-01*i)*zz1
    + (-8.21936276409882E-01 + 5.69579456723519E-01*i)*zz2 + 6;
    + x1^2*x3^2 - 3*x1^2*x3 - 3*x1*x3^2 + 2*x1^2 + 9*x1*x3 + 2*x3^2 - 6*x1 - 6*x3
    + (1.55033404443160E-01-9.87909228374127E-01*i)*zz1
    + (-6.80008285278479E-01-7.33204427122902E-01*i)*zz2 + 4;
    x1*x2*x3*x4 - x1*x2*x3 - x1*x2*x4 - x1*x3*x4 - x2*x3*x4
    + x1*x2 + x1*x3 + x1*x4 + x2*x3 + x2*x4 + x3*x4
    - x1 - x2 - x3 - x4 + (-8.97158012682676E-01-4.41709746642829E-01*i)*zz1
    + (-2.11981526746876E-01 + 9.77273673193985E-01*i)*zz2 + 1;
    + (-3.53404925407865E-01-1.02449352102111E-02*i)*x1
    + (-3.23378971672498E-01-1.42919700111768E-01*i)*x2
    + (2.20216662912369E-01-2.76594687902244E-01*i)*x3
    + (8.35397253138817E-02-3.43542012415485E-01*i)*x4
    + (-2.67791727313329E-01-2.30841050904175E-01*i)*zz1
    + (3.28014183926858E-01-1.31934435015266E-01*i)*zz2 - 3.53553390593274E-01;
    + (-4.95361781887585E-01 + 5.82835865995319E-03*i)*x1
    + (1.04868230499720E-01 + 2.82848378492330E-01*i)*x2i
    + (-1.92957771968200E-01 + 1.72592823677708E-01*i)*x3
    + (1.50836540788987E-01-4.65440876777183E-01*i)*x4
    + (2.59311098918677E-01-3.82322947183237E-02*i)*zz1
    + (-2.61338053548032E-01-1.91269889240060E-01*i)*zz2
    + (2.54532009330191E-01 + 1.49441343387202E-02*i);
    the generic points :
    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x1 :  2.00000000000000E+00  -1.02870981636522E-15
     x2 :  1.00000000000000E+00  -6.89875692784332E-17
     x3 :  5.41257790254410E-01   1.83800265575193E+00
     x4 :  7.57217605925192E-01   2.01695478853197E+00
     zz1 : -1.53177350208010E-16   1.50858780682761E-15
     zz2 :  1.11064312344268E-15   9.39077236860072E-16
    == err :  4.130E-14 = rco :  4.252E-02 = res :  6.734E-15 =
    set of dimension 3 has degree 1
    the polynomials :
      x1^4 - 10*x1^3 + 35*x1^2 - 50*x1
    + (9.98263256449285E-01-5.89107021114901E-02*i)*zz1
    + (-4.67731979750726E-01 + 8.83870349722439E-01*i)*zz2
    + (-3.46149707492286E-01-9.38179289903057E-01*i)*zz3 + 24;
     + x1*x2^3 - 6*x1*x2^2 - x2^3 + 11*x1*x2 + 6*x2^2 - 6*x1 - 11*x2
    + (8.25826300922299E-02-9.96584220829855E-01*i)*zz1
    + (-8.21936276409882E-01 + 5.69579456723519E-01*i)*zz2
    + (-9.94974897178992E-01-1.00124692177572E-01*i)*zz3 + 6;
     + x1^2*x3^2 - 3*x1^2*x3 - 3*x1*x3^2 + 2*x1^2 + 9*x1*x3 + 2*x3^2
    - 6*x1 - 6*x3 + (1.55033404443160E-01-9.87909228374127E-01*i)*zz1
    + (-6.80008285278479E-01-7.33204427122902E-01*i)*zz2
    + (9.94504055917016E-01 + 1.04698055209280E-01*i)*zz3 + 4;
      x1*x2*x3*x4 - x1*x2*x3 - x1*x2*x4 - x1*x3*x4 - x2*x3*x4
    + x1*x2 + x1*x3 + x1*x4 + x2*x3 + x2*x4 + x3*x4 - x1 - x2 - x3 - x4
    + (-8.97158012682676E-01-4.41709746642829E-01*i)*zz1
    + (-2.11981526746876E-01 + 9.77273673193985E-01*i)*zz2
    + (7.37616926767541E-01 + 6.75219423110746E-01*i)*zz3 + 1;
     + (-3.53404925407865E-01-1.02449352102111E-02*i)*x1
    + (-3.23378971672498E-01-1.42919700111768E-01*i)*x2
    + (2.20216662912369E-01-2.76594687902244E-01*i)*x3
    + (8.35397253138817E-02-3.43542012415485E-01*i)*x4
    + (-2.67791727313329E-01-2.30841050904175E-01*i)*zz1
    + (3.28014183926858E-01-1.31934435015266E-01*i)*zz2
    + (-1.52787277524516E-01 + 3.18835455723868E-01*i)*zz3 - 3.53553390593274E-01;
     + (-4.95361781887585E-01 + 5.82835865995319E-03*i)*x1
    + (1.04868230499720E-01 + 2.82848378492330E-01*i)*x2
    + (-1.92957771968200E-01 + 1.72592823677708E-01*i)*x3
    + (1.50836540788987E-01-4.65440876777183E-01*i)*x4
    + (2.59311098918677E-01-3.82322947183237E-02*i)*zz1
    + (-2.61338053548032E-01-1.91269889240060E-01*i)*zz2
    + (-3.36536283081467E-01-7.29526150129006E-02*i)*zz3
    + (2.54532009330191E-01 + 1.49441343387202E-02*i);
     + (1.15184561011707E-01 + 1.13377932381136E-01*i)*x1
    + (-5.17726679477236E-01-8.32947479684604E-02*i)*x2
    + (-3.72758479505006E-01-2.38502017467544E-02*i)*x3
    + (-9.43470497132125E-02-1.99900948005935E-01*i)*x4
    + (-2.38409509176001E-01-2.75632135792588E-01*i)*zz1
    + (-1.20802218865972E-01 + 2.04083624753821E-01*i)*zz2
    + (1.70161843122080E-02-4.70880850401632E-01*i)*zz3
    + (8.75322932557317E-02-3.02958515664302E-01*i);
      the generic points :
      Solution 1 :
      t :  1.00000000000000E+00   0.00000000000000E+00
      m : 1
      the solution for t :
       x1 :  1.00000000000000E+00   8.42939386835822E-17
       x2 : -4.63024527132596E-01  -8.07963582510971E-01
       x3 :  1.38259558372538E+00   3.26937976582811E-01
       x4 :  2.04312757081212E-01   7.58951450601146E-01
       zz1 : -7.86154493035376E-16   7.63862347412176E-16
       zz2 :  2.83697554118352E-16   9.21340974182504E-16
       zz3 :  1.05678537850043E-16   6.86142378280738E-17
      == err :  6.527E-15 = rco :  1.666E-02 = res :  3.737E-15 =

diagonal homotopies
-------------------

An alternative to solving polynomial systems from the top to the bottom,
is to start intersection the equations one after the other.
Consider the intersection of a sphere with a cylinder.

::

    sphere = 'X^2 + Y^2 + Z^2 - 1;'
    cylinder = 'X^2 + 1.0e-14*Y^2 + (Z - 0.5)^2 - 1;'

Observe the tiny coefficient of ``Y^2`` which is a trick to align 
the symbols of the two equations.
The upper case letters of the variables are needed for the ``verify``
not to be confused by ``zz1``.

First, witness sets are constructed for the two hypersurfaces.

::

    from phcpy.sets import double_hypersurface_set

In each instance we check the number of generic points computed, 
which should be ``2`` in both.

::

    (spheqs, sphpts) = double_hypersurface_set(3, sphere)
    len(sphpts)

which indeed shows ``2`` and then also

::

    (cyleqs, cylpts) = double_hypersurface_set(3, cylinder)
    len(cylpts)

shows ``2`` as the number of generic points.

::

    from phcpy.diagonal import double_diagonal_solve
    quaeqs, quapts = double_diagonal_solve(3, 2, spheqs, sphpts, 2, cyleqs, cylpts)

The polynomials in the computed witness set of the intersection are ..."

::

    for pol in quaeqs:
        print(pol)

shown below:

::

     + X^2 + Y^2 + Z^2 + (1.66568720348346E-01 + 9.86029848129109E-01*i)*zz1 - 1;
     + X^2 + 1.00000000000000E-14*Y^2 + Z^2 - Z + (7.86376622880735E-01 + 6.17747365018006E-01*i)*zz1 - 7.50000000000000E-01;
    zz1;
     + (1.58876905105528E-01-7.26285688514595E-01*i)*X + (-8.16467339319674E-01 + 1.71502319167546E+00*i)*Y + (-4.76416867298768E-02 + 4.42262587408809E-01*i)*Z + (5.34984994186327E-01 + 8.44861560254374E-01*i)*zz1+(-8.46146634046559E-01 + 5.32950160607611E-01*i);

and the generic points in the computed witness set are printed with

::

    for (idx, sol) in enumerate(quapts):
        print('Solution', idx+1, ':')
        print(sol)

and the output is

::

    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     X :  9.91243806839434E-01  -1.45160013935997E-02
     Y : -1.36376604565112E-01  -3.29060036857026E-01
     Z :  3.39681929583638E-01  -8.97521810492629E-02
     zz1 :  0.00000000000000E+00   0.00000000000000E+00
    == err :  3.374E-16 = rco :  4.134E-02 = res :  1.769E-16 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     X : -7.67613124615721E-01   1.54768674480021E-01
     Y : -6.69973298698680E-01  -1.30010400626017E-01
     Z : -1.81961516698249E-01  -1.74206993945098E-01
     zz1 :  0.00000000000000E+00   0.00000000000000E+00
    == err :  3.448E-16 = rco :  5.434E-02 = res :  3.886E-16 =
    Solution 3 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     X :  4.60320208343188E+00   4.38623359269609E+00
     Y :  9.59637373859308E-01   2.36891289291490E+00
     Z :  4.94084440491082E+00  -4.54659469491659E+00
     zz1 :  0.00000000000000E+00   0.00000000000000E+00
    == err :  3.427E-14 = rco :  6.196E-04 = res :  2.465E-14 =
    Solution 4 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     X :  6.32992398543307E+00   2.24044198839475E+00
     Y :  2.07264148009942E+00  -1.51003268649074E+00
     Z : -1.76564399075824E+00   6.25951276465330E+00
     zz1 :  0.00000000000000E+00   0.00000000000000E+00
    == err :  4.575E-14 = rco :  1.103E-03 = res :  2.687E-14 =

Four generic points are obtained in the computed witness set,
as the intersection of a sphere with a cylinder is a quartic.
