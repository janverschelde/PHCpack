numerical Schubert calculus
===========================

The module schubert.py exports the hypersurface and quantum
Pieri homotopies to solve the following Schubert problem:
Given a sequence of generic m-planes and a corresponding
sequence of interpolation points, compute all maps of
degree q that meet the given m-planes nontrivially at
the interpolation points.

Pieri homotopies
----------------

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

   >>> from phcpy.schubert import pieri_root_count
   >>> (m,p,q) = (2,2,1)
   >>> n = m*p + q*(m+p)
   >>> r = pieri_root_count(m,p,q)
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
   >>> from phcpy.schubert import random_complex_matrix
   >>> L = [random_complex_matrix(m+p,m) for k in range(n)]
   >>> points = random_complex_matrix(n,1)
   >>> from phcpy.schubert run_pieri_homotopies
   >>> (f,fsols) = run_pieri_homotopies(m,p,q,L,points)

The function **test()** of the module ``schubert``
runs an interactive session to solve instances that
are fully real (in case q = 0).

Littlewood-Richardson homotopies 
--------------------------------

With the Littlewood-Richardson homotopies 
we can solve general Schubert problems.
The input to a Schubert problem is a sequence of n-by-n matrices and a
corresponding list of intersection conditions, represented by brackets.
For example, the bracket [2, 4, 6] imposes on a 3-plane in 6-space that it
meets nontrivially the space spanned by the first two columns of the 
corresponding matrix in a line and that it meets 
the space spanned by the first four columns of
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
