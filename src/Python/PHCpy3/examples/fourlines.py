"""
Given four lines in general position,
there are two lines which meet all four given lines.
With Pieri homotopies we can solve this Schubert problem.
"""

def solve_general(mdim, pdim, qdeg):
    """
    Solves a general instance of Pieri problem, computing the
    p-plane producing curves of degree qdeg which meet a number
    of general m-planes at general interpolation points,
    where p = pdim and m = mdim on input.
    For the problem of computing the two lines which meet
    four general lines, mdim = 2, pdim = 2, and qdeg = 0.
    Returns the number of solutions computed.
    """
    from phcpy.schubert import random_complex_matrix
    from phcpy.schubert import run_pieri_homotopies
    dim = mdim*pdim + qdeg*(mdim+pdim)
    planes = [random_complex_matrix(mdim+pdim, mdim) for _ in range(0, dim)]
    (pols, sols) = run_pieri_homotopies(mdim, pdim, qdeg, planes)
    return len(sols)

def main():
    """
    We start with the formalism of the root count,
    solve a general configuration and then a special problem.
    """
    from phcpy.schubert import pieri_root_count
    (mdim, pdim, deg) = (2, 2, 0)
    pcnt = pieri_root_count(mdim, pdim, deg, False)
    print('The Pieri root count :', pcnt)
    nbr = solve_general(mdim, pdim, deg)
    print('The number of solutions :', nbr)

if __name__ == "__main__":
    main()
