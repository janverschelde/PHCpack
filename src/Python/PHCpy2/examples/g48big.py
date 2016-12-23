"""
The script tests the Littlewood-Richardson homotopies on some
relatively big problem in G(4,8).  The ten problems were formulated
by Frank Sottile at an AIM squares meeting in the Fall of 2015.
(0) [3 6 7 8]^7*[4 6 7 8]^2 = 231
(1) [2 5 6 8]*[3 6 7 8]*[4 5 7 8]*[4 6 7 8]^7 = 294
(2) [2 4 7 8]*[3 6 7 8]^2*[4 6 7 8]^7 = 356
(3) [2 5 7 8]^2*[4 6 7 8]^8 = 398
(4) [3 5 7 8]*[4 5 7 8]^2*[3 6 7 8]^2*[4 6 7 8]^5 = 526
(5) [4 5 7 8]^3*[3 6 7 8]^2*[4 6 7 8]^6 = 734
(6) [3 5 7 8]*[3 6 7 8]^3*[4 6 7 8]^7 = 1104
(7) [3 6 7 8]^4*[4 6 7 8]^8 = 1600
(8) [3 5 7 8]*[3 6 7 8]^2*[4 6 7 8]^9 = 2166
(9) [3 6 7 8]^3*[4 6 7 8]^10 = 3102
"""

def ten_problems():
    """
    Returns a list of ten Schubert problems in G(4, 8) as a list
    of tuples.  The first element in each tuple contains the brackets
    and the second element the corresponding root count.
    """
    result = []
    b0 = [[3, 6, 7, 8] for _ in range(7)] \
       + [[4, 6, 7, 8] for _ in range(2)]
    r0 = 231
    result.append((b0, r0))
    b1 = [[2, 5, 6, 8], [3, 6, 7, 8], [4, 5, 7, 8]] \
       + [[4, 6, 7, 8] for _ in range(7)]
    r1 = 294
    result.append((b1, r1))
    b2 = [[2, 4, 7, 8]] + [[3, 6, 7, 8] for _ in range(2)] \
       + [[4, 6, 7, 8] for _ in range(7)]
    r2 = 356
    result.append((b2, r2))
    b3 = [[2, 5, 7, 8] for _ in range(2)] \
       + [[4, 6, 7, 8] for _ in range(8)]
    r3 = 398
    result.append((b3, r3))
    b4 = [[3, 5, 7, 8]] \
       + [[4, 5, 7, 8] for _ in range(2)] \
       + [[3, 6, 7, 8] for _ in range(2)] \
       + [[4, 6, 7, 8] for _ in range(5)]
    r4 = 526
    result.append((b4, r4))
    b5 = [[4, 5, 7, 8] for _ in range(3)] \
       + [[3, 6, 7, 8] for _ in range(2)] \
       + [[4, 6, 7, 8] for _ in range(6)]
    r5 = 734
    result.append((b5, r5))
    b6 = [[3, 5, 7, 8]] \
       + [[3, 6, 7, 8] for _ in range(3)] \
       + [[4, 6, 7, 8] for _ in range(7)]
    r6 = 1104
    result.append((b6, r6))
    b7 = [[3, 6, 7, 8] for _ in range(4)] \
       + [[4, 6, 7, 8] for _ in range(8)]
    r7 = 1600
    result.append((b7, r7))
    b8 = [[3, 5, 7, 8]] \
       + [[3, 6, 7, 8] for _ in range(2)] \
       + [[4, 6, 7, 8] for _ in range(9)]
    r8 = 2166
    result.append((b8, r8))
    b9 = [[3, 6, 7, 8] for _ in range(3)] \
       + [[4, 6, 7, 8] for _ in range(10)]
    r9 = 3102
    result.append((b9, r9))
    return result

def count_roots(problems):
    """
    Runs the root count on the ten problems.
    Checks whether the computed root corresponds with the exact one.
    An exception will be raised if an error occurs in the check.
    """
    from phcpy.schubert import resolve_schubert_conditions as root_count
    (k, n) = (4, 8) # in G(4,8)
    for (brackets, count) in problems:
        c = root_count(n, k, brackets)
        print brackets, '=', c
        assert(count == c)

def solve(problems):
    """
    Runs the Littlewood-Richardson homotopies on the problems.
    """
    from phcpy.schubert import littlewood_richardson_homotopies as lrhom
    (k, n) = (4, 8) # in G(4,8)
    for (brackets, count) in problems:
        print brackets, '=', count
        (rc, flags, pols, sols) = lrhom(n, k, brackets)
        assert(len(sols) == count)

def main():
    """
    Generates ten problems and counts the roots.
    """
    prbs = ten_problems()
    count_roots(prbs)
    solve(prbs)

main()
