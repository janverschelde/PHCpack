"""
Given four lines in general position,
there are two lines which meet all four given lines.
With Pieri homotopies we can solve this Schubert problem.
For the verification of the intersection conditions, numpy is used.
"""
from numpy import zeros, array, concatenate, matrix

def indices(name):
    """
    For the string name in the format xij
    return (i, j) as two integer indices.
    """
    return (int(name[1]), int(name[2]))

def solution_plane(rows, cols, sol):
    """
    Returns a sympy matrix with as many rows
    as the value of rows and with as many columns
    as the value of columns, using the string
    represention of a solution in sol.
    """
    from phcpy.solutions import coordinates
    result = zeros((rows, cols), dtype=complex)
    for k in range(cols):
        result[k][k] = 1
    (vars, vals) = coordinates(sol)
    for (name, value) in zip(vars, vals):
        i, j = indices(name)
        result[i-1][j-1] = value
    return result

def verify_determinants(inps, sols, verbose=True):
    """
    Verifies the intersection conditions with determinants,
    concatenating the planes in inps with those in the sols.
    Both inps and sols are lists of numpy arrays.
    Returns the sum of the absolute values of all determinants.
    If verbose, then for all solutions in sols, the computed
    determinants are printed to screen.
    """
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
    from phcpy.schubert import random_complex_matrix
    from phcpy.schubert import run_pieri_homotopies
    dim = mdim*pdim + qdeg*(mdim+pdim)
    ranplanes = [random_complex_matrix(mdim+pdim, mdim) for _ in range(0, dim)]
    (pols, sols) = run_pieri_homotopies(mdim, pdim, qdeg, ranplanes, \
        verbose=False)
    inplanes = [array(plane) for plane in ranplanes]
    outplanes = [solution_plane(mdim+pdim, pdim, sol) for sol in sols]
    return (inplanes, outplanes, pols, sols)

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

def main():
    """
    We start with the formalism of the root count,
    solve a general configuration and then a special problem.
    """
    from phcpy.schubert import pieri_root_count
    (mdim, pdim, deg) = (2, 2, 0)
    pcnt = pieri_root_count(mdim, pdim, deg, False)
    print('The Pieri root count :', pcnt)
    print('Solving a general case ...')
    (inp, otp, pols, sols) = solve_general(mdim, pdim, deg)
    print('The input planes :')
    for plane in inp:
        print(plane)
    print('The solution planes :')
    for plane in otp:
        print(plane)
    check = verify_determinants(inp, otp)
    print('Sum of absolute values of determinants :', check)
    input('Hit enter to continue.')
    (oscp, otp2, pols2, sols2) = solve_real(mdim, pdim, pols, sols)
    print('The input planes :')
    for plane in oscp:
        print(plane)
    print('The solution planes :')
    for plane in otp2:
        print(plane)
    check = verify_determinants(oscp, otp2)
    print('Sum of absolute values of determinants :', check)

if __name__ == "__main__":
    main()
