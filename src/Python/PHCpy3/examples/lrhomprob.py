"""
Problems for the Littlewood-Richardson homotopies.
"""
def verify(pols, sols):
    """
    Verifies whether the solutions in sols
    satisfy the polynomials of the system in pols.
    """
    from phcpy.solutions import strsol2dict, evaluate
    dictsols = [strsol2dict(sol) for sol in sols]
    checksum = 0
    for sol in dictsols:
        sumeval = sum(evaluate(pols, sol))
        print(sumeval)
        checksum = checksum + sumeval
    print('the total check sum :', checksum)

def test_basic(prc='d'):
    """
    Performs a basic test on the Littlewood-Richardson homotopies.
    """
    from phcpy.schubert import littlewood_richardson_homotopies as lrhom
    brk = [[2, 4, 6], [2, 4, 6], [2, 4, 6]]
    print('basic test with brackets', brk)
    (roco, flags, fsys, sols) = lrhom(6, 3, brk, precision=prc)
    print('the root count :', roco)
    print('the flags :', flags)
    print('the solutions :')
    for sol in sols:
        print(sol)
    verify(fsys, sols)
    print('number of solutions :', len(sols))

def test_p1(prc='d'):
    """
    First big problem in G(4,8).
    """
    from phcpy.schubert import littlewood_richardson_homotopies as lrhom
    brk = [[3, 6, 7, 8], [3, 6, 7, 8], [3, 6, 7, 8], [3, 6, 7, 8], \
        [3, 6, 7, 8], [3, 6, 7, 8], [3, 6, 7, 8], [4, 6, 7, 8], \
        [4, 6, 7, 8]]
    print('basic test with brackets', brk)
    (roco, flags, fsys, sols) = lrhom(8, 4, brk, precision=prc)
    print('the root count :', roco)
    print('the flags :', flags)
    print('the solutions :')
    for sol in sols:
        print(sol)
    verify(fsys, sols)
    print('number of solutions :', len(sols))

def main():
    """
    Tests the Littlewood-Richardson homotopies.
    """
    print("\nTesting the Littlewood-Richardson homotopies ...")
    test_basic()
    test_p1()

if __name__ == "__main__":
    main()
