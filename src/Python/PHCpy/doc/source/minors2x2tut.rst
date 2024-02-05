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

As a summary, the contents of ``deco`` is written as

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
