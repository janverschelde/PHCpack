Homotopies for Problems in Enumerative Geometry
===============================================

Problems in enumerative geometry have a particular structure, 
well suited for polynomial homotopies.  
Based on the Pieri root counts and the Littlewood-Richardson rule, 
there are homotopies which are *generical optimal* in the sense that 
the number of solution paths matches the number of solutions, 
for generic problems.

Pieri homotopies
----------------

A classical problem in enumerative goes as follows.
Given four lines in 3-space, how many lines do meet those four given lines
in a point?  Although seemingly a linear problem, 
it turns out that the answer is two.  
For four given lines in general position, 
there are two lines meeting the four given lines.

The code below applies Pieri homotopies to compute all lines meeting
four random lines in 3-space.

::

    from phcpy.schubert import pieri_root_count

A line is represented by a matrix with 2 columns.  
Therefore, both ``m``. the dimension of the input, and ``p``,
the dimension of the output, are both equal to ``2``.

::

    pieri_root_count(2, 2)

returns ``2``.  Let us compute those lines.

::

    from phcpy.schubert import random_complex_matrix, run_pieri_homotopies

The setup for the Pieri problems are formulated for general
``m``, ``p`` which requires ``m*p`` inputs for the problem
to have finitely many isolated solutions.

::

    (m, p) = (2, 2)
    dim = m*p
    L = [random_complex_matrix(m+p, m) for _ in range(dim)]

With the setup done, the homotopies are defined and all paths
are tracked in the following code cell.

::

    (fsys, fsols) = run_pieri_homotopies(m, p, 0, L)

The code to print the solutions

::

    for (idx, sol) in enumerate(fsols):
        print('Solution', idx+1, ':')
        print(sol)

shows

::

    Solution 1 :
    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x21 : -5.90757451491546E-01  -2.73796238095519E-02
     x31 : -9.20374122953531E-01   2.18442665351047E-01
     x32 :  8.14344855671798E-01  -5.32153788006048E-01
     x42 : -5.75636214663051E-01  -2.78447850739093E-01
    == err :  2.304E-15 = rco :  3.053E-02 = res :  1.388E-15 =
    Solution 2 :
    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x21 :  6.13353565245472E-01   5.36374877135208E-01
     x31 : -1.40790300814705E-01   3.49424087655845E-02
     x32 : -5.00597509424303E-01  -1.02203061067751E+00
     x42 : -1.74780414053227E+00   6.25395371456727E-01
    == err :  2.558E-15 = rco :  6.714E-02 = res :  2.331E-15 =

The names of the variables use the indexing to denote the position 
in the matrix of the generators of the line.  
The formal root count is summarized in the poset of localization patterns.

::

    from phcpy.schubert import pieri_localization_poset
    poset22 = pieri_localization_poset(2, 2)
    print(poset22)

and the poset for this problem is

::

    n = 0 : ([2 3],[2 3],1)([1 4],[1 4],1)
    n = 1 : 
    n = 2 : ([1 3],[2 4],2)
    n = 3 :
    n = 4 : ([1 2],[3 4],2)

If the degree ``q`` is nonzero, the problem of p-planes 
meeting m-planes is extended to computing curves of degree ``q``
that produce p-planes that meet m-planes at interpolation points.
For example, consider line producing interpolating curves.

::

    poset221 = pieri_localization_poset(2, 2, 1)
    print(poset221)

::

    n = 0 : ([3 4],[3 4],1)([2 5],[2 5],1)
    n = 1 :
    n = 2 : ([2 4],[3 5],2)
    n = 3 :
    n = 4 : ([2 3],[3 6],2)([2 3],[4 5],2)([1 4],[3 6],2)([1 4],[4 5],2)
    n = 5 :
    n = 6 : ([1 3],[4 6],8)
    n = 7 :
    n = 8 : ([1 2],[4 7],8)

There are ``8`` solutions as can be computed by the following code:

::

    (m, p, q) = (2, 2, 1)
    dim = m*p + q*(m+p)
    roco = pieri_root_count(m, p, q)

In addition to m-planes, interpolation points must be defined.

::

    L = [random_complex_matrix(m+p, m) for _ in range(dim)]
    points = random_complex_matrix(dim, 1)

Then all paths defined by the Pieri homotopies to compute 
line producing curves are tracked by the statement.

::

    (fsys, fsols) = run_pieri_homotopies(m, p, q, L, 0, points)

Then ``len(fsols)`` returns ``8``, maching the root count.

::

    for (idx, sol) in enumerate(fsols):
        print('Solution', idx+1, ':')
        print(sol)

and the solutions are

::

    Solution 1 :
    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x11s0 : -2.14906949033699E-01  -6.39278480067621E-01
     x21s0 : -2.15094936097349E-01   7.13137735748634E-01
     x31s0 :  4.89419546169109E-01   7.96444688560527E-02
     x41s0 :  3.20065497988893E-01   5.85095426596230E-01
     x22s0 :  2.31759428669014E-01  -2.62654162829177E-01
     x32s0 : -9.91726822253624E-01  -6.97520762259368E-01
     x42s0 : -9.55529418833548E-01   7.00387732447402E-01
     x12s1 :  9.33661460517070E-01  -6.12105981775290E-01
    == err :  0.000E+00 = rco :  0.000E+00 = res :  0.000E+00 =
    Solution 2 :
    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x11s0 :  1.10520066466899E+00  -3.50325368846109E-01
     x21s0 :  5.18574648064343E-01   1.32725970575267E-01
     x31s0 :  2.62487605641133E-01   3.74101095634462E-01
     x41s0 :  5.47831843633738E-01   7.39727656004381E-01
     x22s0 : -3.83564646759339E-01  -3.63155591886385E-01
     x32s0 : -4.18752363921196E-01  -1.00720592449952E-01
     x42s0 : -1.73131137487804E+00   2.32539054040540E-01
     x12s1 :  2.63582010113658E-01  -5.40149341401405E-01
    == err :  0.000E+00 = rco :  0.000E+00 = res :  0.000E+00 =
    Solution 3 :
    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x11s0 : -1.09406342297013E+00   1.40229579344629E+00
     x21s0 : -9.14923921567758E-01  -5.69034530991378E-02
     x31s0 : -2.57253733233100E-01   1.11315272481983E+00
     x41s0 : -1.23385040712015E+00  -3.62148913829782E+00
     x22s0 : -1.48508415565715E+00  -1.46881674145667E+00
     x32s0 : -2.18057639295167E+00  -1.82823489900416E+00
     x42s0 :  3.49088194334639E+00   2.33832870569806E+00
     x12s1 : -7.99089515647370E-01   2.01242656315116E+00
    == err :  0.000E+00 = rco :  0.000E+00 = res :  0.000E+00 =
    Solution 4 :
    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x11s0 :  6.20220385308888E-01   7.14013088945019E-01
     x21s0 :  1.96622054165017E-01   6.37755854895062E-01
     x31s0 :  3.84896917178092E-01  -2.74979883047354E-01
     x41s0 :  2.85390395538486E-01   2.21952044787881E-01
     x22s0 : -3.21807655249411E-01  -5.36746050197169E-01
     x32s0 :  1.15438122582161E-01  -9.43108020255466E-01
     x42s0 : -4.47775301334196E-01  -1.25556346180423E-01
     x12s1 :  7.52973406987737E-01  -4.13342381841601E-01,
    == err :  0.000E+00 = rco :  0.000E+00 = res :  0.000E+00 =
    Solution 5 :
    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x11s0 : -6.40605568129237E-01   5.53687431204576E-01
     x21s0 : -3.17434054891094E-01  -1.17563450907422E+00
     x31s0 : -1.27461906338159E+00   5.02576168841663E-01
     x41s0 : -1.10280584587259E+00   4.85210021219730E-01
     x22s0 :  8.42063423093745E-01   8.61583481543135E-01
     x32s0 :  6.33194409804748E-01   6.14896771183676E-01
     x42s0 :  3.64149459917250E-01   2.12619523294958E-01
     x12s1 : -2.47538557540740E-01   4.66198793734136E-02,
    == err :  0.000E+00 = rco :  0.000E+00 = res :  0.000E+00 =
    Solution 6 :
    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x11s0 :  6.09038672836626E-01  -8.98742630676646E-01
     x21s0 :  7.42293129431968E-02   9.26153607979929E-01
     x31s0 :  8.37592049308958E-01  -4.43304301590072E-01
     x41s0 :  2.19220600974304E-02   1.14241100871086E+00,
     x22s0 :  9.44017626164027E-01  -2.71712341229773E-01
     x32s0 : -5.64425652511392E-01  -2.36806482566213E-01
     x42s0 : -7.39375051122763E-01  -4.58006540806617E-01
     x12s1 :  5.62059852022753E-01   2.93332939159073E-02
    == err :  0.000E+00 = rco :  0.000E+00 = res :  0.000E+00 =
    Solution 7 :
    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x11s0 : -7.07534912238552E-01   4.92999096377366E-02
     x21s0 :  4.79149192690524E-02   6.23651158687302E-01
     x31s0 :  4.35563235063921E-01   2.57626463594513E-01
     x41s0 :  2.93479058252621E-01   2.14109168395354E-01
     x22s0 : -1.20085603380698E-01  -3.52342158354991E-01
     x32s0 : -1.14999164093347E+00  -3.86601319323260E-01
     x42s0 : -1.50274237644790E-01   5.36877541048399E-01
     x12s1 :  2.25543133569722E-01  -1.00744217331902E+00
    == err :  0.000E+00 = rco :  0.000E+00 = res :  0.000E+00 =
    Solution 8 :
    t :  0.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x11s0 : -5.04261904420292E-02   1.59417581791122E+00
     x21s0 : -6.78725206430762E-01  -2.02091391977252E-01
     x31s0 :  4.63632342073102E-01   6.74339361230072E-01
     x41s0 :  1.10851072800672E+00   1.07992063150119E+00
     x22s0 :  9.14318142078649E-01   2.67548922924622E-01
     x32s0 :  2.66219598193448E-01   1.24867156045807E+00
     x42s0 : -2.64481349873839E+00  -2.05243627160146E-01
     x12s1 :  3.30575444135274E-01  -1.42397007492760E+00
    == err :  0.000E+00 = rco :  0.000E+00 = res :  0.000E+00 =

The index following the ``s`` in the variable name represents
the degree of the parameter ``s`` 
in the curve that produces lines in 3-space.

Littlewood-Richardson homotopies
--------------------------------

A Schubert condition is represented by a sequence of brackets.
Each bracket represents conditions on the dimensions 
of the intersections with the given inputs.

With Littlewood-Richardson rule, we count the number of solutions,
resolving the Schubert condition.

::

    from phcpy.schubert import resolve_schubert_conditions

The intersection conditions are defined below.

::

   brackets = [[2, 4, 6], [2, 4, 6], [2, 4, 6]]

We are looking for 3-planes $X$ in 6-planes that meet flags as follows:

.. math::

   1. \mbox{dim}(X \cap \langle f_1, f_2 \rangle) = 1.

   2. \mbox{dim}(X \cap \langle f_1, f_2, f_3, f_4 \rangle) = 2.

   3. \mbox{dim}(X \cap \langle f_1, f_2, f_3, f_4, f_5, f_6 \rangle) = 3.

For these conditions, there are finitely many solutions :math:`X`.
The number of solution is computed as follows.

::

    roco = resolve_schubert_conditions(6, 3, brackets)

Littlewood-Richardson homotopies track exactly as many paths
as the value of ``roco``, which is ``2`` for this problem.

::

    from phcpy.schubert import double_littlewood_richardson_homotopies as lrh

Tracking all paths of 3-planes in 6-space defined by
Littlewood-Richardson homotopies is done by the execution
of the statement.

::

    (count, flags, sys, sols) = lrh(6, 3, brackets, verbose=False)

The value of ``count`` is ``2`` and the solutions are

::

    Solution 1 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x11 : -1.71828539203956E+00   3.70396971521702E-01
     x32 : -9.38154978067327E-01   4.39465496011351E-01
     x53 : -4.43650959809938E-01   9.55468566341054E-02
    == err :  0.000E+00 = rco :  1.000E+00 = res :  4.785E-16 =
    Solution 2 :
    t :  1.00000000000000E+00   0.00000000000000E+00
    m : 1
    the solution for t :
     x11 : -6.49381975027210E-01  -4.99975537206415E-01
     x32 : -1.40857387994158E+00   4.74177393405449E-01
     x53 : -7.94711695711224E-01  -1.11583537216770E-01
    == err :  0.000E+00 = rco :  1.000E+00 = res :  4.785E-16 =

and once again, the indices of the variable names indicate
the position of the numbers in the 3-planes.
