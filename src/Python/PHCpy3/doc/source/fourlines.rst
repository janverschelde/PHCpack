lines meeting four given lines
==============================

Consider as given four lines, our problem is to compute all
lines which meet the four given lines in a point.

A line in projective 3-space is represented by two points,
stored in the columns of a 4-by-2 matrix.
So the space we work in is the complex 4-space.
For this problem we have a formal root count,
named after Pieri.

::

   from phcpy.schubert import pieri_root_count
   rc = pieri_root_count(2, 2, 0, verbose=False)

In 4-space, the dimension of the input planes equals two
and also the dimension of the output planes is two.
The value returned in ``rc`` is two for this problem.

a general configuration
-----------------------

In a general configuration, random number generators are applied
to determine the points which span the input lines.
The solving of a general configuration is encapsulated in the
function ``solve_general``.

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
       from numpy import array
       from phcpy.schubert import random_complex_matrix
       from phcpy.schubert import run_pieri_homotopies
       dim = mdim*pdim + qdeg*(mdim+pdim)
       ranplanes = [random_complex_matrix(mdim+pdim, mdim) \
           for _ in range(0, dim)]
       (pols, sols) = run_pieri_homotopies(mdim, pdim, qdeg, ranplanes, \
           verbose=False)
       inplanes = [array(plane) for plane in ranplanes]
       outplanes = [solution_plane(mdim+pdim, pdim, sol) for sol in sols]
       return (inplanes, outplanes, pols, sols)

The solutions returned by ``run_pieri_homotopies`` are converted into
numpy matrices, as defined by the function ``solution_plane``.

::

   def solution_plane(rows, cols, sol):
       """
       Returns a sympy array with as many rows
       as the value of rows and with as many columns
       as the value of columns, using the string
       represention of a solution in sol.
       """
       from numpy import zeros
       from phcpy.solutions import coordinates
       result = zeros((rows, cols), dtype=complex)
       for k in range(cols):
           result[k][k] = 1
       (vars, vals) = coordinates(sol)
       for (name, value) in zip(vars, vals):
           i, j = (int(name[1]), int(name[2]))
           result[i-1][j-1] = value
       return result

For the verification of the intersection conditions,
the matrices of the input planes are concatenated to the solution planes
and the determinant of the concatenated matrix is computed.

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
       from numpy import matrix
       from numpy.linalg import det
       checksum = 0
       for sol in sols:
           if verbose:
               print('checking solution\n', sol)
           for plane in inps:
               cat = concatenate([plane, sol], axis=-1)
               mat = matrix(cat)
               dcm = det(mat)
               if verbose:
                   print('the determinant :', dcm)
               checksum = checksum + abs(dcm)
       return checksum

Then the ``main()`` function contains the following code.

::

   (inp, otp, pols, sols) = solve_general(mdim, pdim, deg)
   print('The input planes :')
   for plane in inp:
       print(plane)
   print('The solution planes :')
   for plane in otp:
       print(plane)
   check = verify_determinants(inp, otp)
   print('Sum of absolute values of determinants :', check)

The polynomial system in ``pols`` with corresponding solutions
in ``sols`` can be used as start system to solve specific problems,
as will be done in the next section.

a real configuration
--------------------

The solution of a real instance takes on input the system
and corresponding solutions of a general instance.

::

   def solve_real(mdim, pdim, start, sols):
       """
       Solves a real instance of Pieri problem, for input planes
       of dimension mdim osculating a rational normal curve.
       On return are the planes of dimension pdim.
       """
       from phcpy.schubert import real_osculating_planes
       from phcpy.schubert import make_pieri_system
       from phcpy.trackers import track
       oscplanes = real_osculating_planes(mdim, pdim, 0)
       target = make_pieri_system(mdim, pdim, 0, oscplanes, False)
       rtsols = track(target, start, sols)
       inplanes = [array(plane) for plane in oscplanes]
       outplanes = [solution_plane(mdim+pdim, pdim, sol) for sol in rtsols]
       return (inplanes, outplanes, target, rtsols)

The code for the ``main()`` is similar as when calling
``solve_general()``, as shown above at the end of the previous section.
