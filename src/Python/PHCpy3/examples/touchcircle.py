"""
This script computes the tangents lines to a circle
via a witness set intersection.
"""
def polynomials(a, b, r):
    """
    Returns string representations of two polynomials:
    1) a circle with radius r centered at (a, b);
    2) a line through the origin with slope s.
    """
    crc = '(x - %.15e)^2 + (y - %.15e)^2 - %.15e;' % (a, b, r**2)
    lin = 'y - s*x;'
    return [crc, lin]

def special_solutions(pols, slope):
    """
    Given in pols the polynomials for the line intersecting a circle
    for a general slope s, solves the problem for a special numerical
    value for s, given in the slope.
    The value of the slope is added as last coordinate to each solution.
    """
    from phcpy.solver import solve
    from phcpy.solutions import make_solution, coordinates
    special = []
    for pol in pols:
        rpl = pol.replace('s', '(' + str(slope) + ')')
        special.append(rpl)
    sols = solve(special, silent=True)
    result = []
    for sol in sols:
        (vars, vals) = coordinates(sol)
        vars.append('s')
        vals.append(slope)
        extsol = make_solution(vars, vals)
        result.append(extsol)
    return result

def make_witness_set(pols, verbose=True):
    """
    We have two equations in three variables in pols
    and therefore we expect a one dimensional set.
    """
    from phcpy.sets import embed
    from phcpy.solver import solve
    embpols = embed(3, 1, pols)
    if verbose:
        print('the embedded system :')
        for pol in embpols:
            print(pol)
    embsols = solve(embpols, silent=not verbose)
    if verbose:
        print('the witness points :')
        for sol in embsols:
            print(sol)
    return (embpols, embsols)

def membership_test(witsys, witsols, point, verbose=True):
    """
    Given a witness sets in the tuple witset,
    runs a membership test on a solution.
    """
    from phcpy.solutions import make_solution
    from phcpy.sets import is_member
    print('testing the point\n', point)
    ismb = is_member(witsys, witsols, 1, point, verbose=False)
    if ismb:
        print('the point is a member')
    else:
        print('the point is NOT a member')
    return ismb

def main():
    """
    Generates the system and its witness set.
    As a sanity check, a membership test is done.
    """
    syst = polynomials(3, 2, 1)
    for pol in syst:
        print(pol)
    spsols = special_solutions(syst, 1)
    for sol in spsols:
        print(sol)
    input('hit enter to continue')
    (embsyst, embsols) = make_witness_set(syst, False)
    print('the polynomials in the witness set:')
    for pol in embsyst:
        print(pol)
    print('the solutions :')
    for sol in embsols:
        print(sol)
    print('degree of the set :', len(embsols))
    for sol in spsols:
        membership_test(embsyst, embsols, sol)

if __name__ == "__main__":
    main()
