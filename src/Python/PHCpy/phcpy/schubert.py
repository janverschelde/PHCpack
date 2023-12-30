"""
Numerical Schubert calculus defines homotopies for enumerative geometry.
"""
from ctypes import c_int, c_int32, c_double, pointer, create_string_buffer
from random import uniform
from phcpy.version import get_phcfun
from phcpy.version import int4a2nbr, nbr2int4a, int4a2str, str2int4a

def pieri_root_count(mdim, pdim, qdeg, vrblvl=0):
    r"""
    Computes the number of *pdim*-plane producing maps of degree *qdeg*
    that meet *mdim*-planes at mdim*pdim + qdeg*(mdim+pdim) points.
    """
    if vrblvl > 0:
        print('in pieri_root_count, m :', mdim, end='')
        print(', p :', pdim, ', q :', qdeg)
    phc = get_phcfun()
    apars = int4a2nbr([mdim, pdim, qdeg], (vrblvl > 0))
    roco = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    if vrblvl > 0:
        print('-> pieri_root_count calls phc ...')
        print('apars =', nbr2int4a(apars))
    retval = phc(223, apars, roco, ccc, vrb)
    if vrblvl > 0:
        print('the return value of pieri_root_count :', retval)
        print('the Pieri root count :', roco[0])
    return roco[0]

def pieri_localization_poset(mdim, pdim, qdeg, size=10240, vrblvl=0):
    r"""
    Returns the string representation of the localization poset
    used to compute the Pieri root count, for the number of *pdim*-plane
    producing maps of degree *qdeg* that meet *mdim*-planes at
    mdim*pdim + qdeg*(mdim+pdim) points.
    The size of the string representation is given on input in size.
    """
    if vrblvl > 0:
        print('in pieri_root_count, m :', mdim, end='')
        print(', p :', pdim, ', q :', qdeg)
    phc = get_phcfun()
    apars = int4a2nbr([mdim, pdim, qdeg], (vrblvl > 0))
    poset = create_string_buffer(b"", 4*size)
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    if vrblvl > 0:
        print('-> pieri_localization_poset calls phc ...')
        print('apars =', nbr2int4a(apars))
    retval = phc(224, apars, poset, ccc, vrb)
    if vrblvl > 0:
        print('the return value of pieri_localization_poset :', retval)
    result = int4a2str(poset, (vrblvl > 0))
    return result

def resolve_schubert_conditions(ndim, kdim, brackets, vrblvl=0):
    r"""
    In n-dimensional space we consider k-dimensional planes,
    subject to intersection conditions represented by brackets.
    The *brackets* is a list of brackets.  A bracket is a list
    of as many natural numbers (in the range 1..*ndim*) as *kdim*.
    On return is the formal root count, which is sharp for general flags.
    and the coordinates of the flags, stored row wise in a list
    of real and imaginary parts.
    """
    if vrblvl > 0:
        print('in resolve_schubert_conditions, n :', ndim, end='')
        print(', k :', kdim, end='')
        print(', #brackets :', len(brackets))
    size = kdim*len(brackets)
    cds = (c_int32 * size)()
    if vrblvl > 0:
        print('conditions :', cds[:])
    idx = 0
    for bracket in brackets:
        for nbr in bracket:
            cds[idx] = c_int32(nbr)
            idx = idx + 1
    if vrblvl > 0:
        print('conditions :', cds[:])
    phc = get_phcfun()
    apars = int4a2nbr([ndim, kdim, len(brackets)], (vrblvl > 0))
    brk = pointer(cds)
    roco = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    if vrblvl > 0:
        print('-> resolve_schubert_conditions calls phc ...')
        print('apars =', nbr2int4a(apars))
    retval = phc(228, apars, brk, roco, vrb)
    if vrblvl > 0:
        print('the return value of resolve_schubert_conditions :', retval)
        print('the root count :', int(roco[0]))
    return int(roco[0])

def test_pieri(vrblvl=0):
    """
    Tests the Pieri root count.
    """
    roco = pieri_root_count(2, 2, 1, vrblvl)
    poset = pieri_localization_poset(2, 2, 1, vrblvl=vrblvl)
    if vrblvl > 0:
        print('The pieri localization poset :')
        print(poset)
    if roco == 8:
        if vrblvl > 0:
            print('Pieri root count', roco, 'is correct.')
        return 0
    if vrblvl > 0:
        print('Pieri root count', roco, 'is wrong.')
    return 1

def test_littlewood_richardson_rule(vrblvl=0):
    """
    Tests the Littlewood-Richardson rule.
    """
    brk = [[2, 4, 6], [2, 4, 6], [2, 4, 6]]
    roco = resolve_schubert_conditions(6, 3, brk, vrblvl)
    if roco == 2:
        if vrblvl > 0:
            print('Littlewood-Richardson root count :', roco, 'is correct.')
        return 0
    if vrblvl > 0:
        print('Littlewood-Richardson root count :', roco, 'is wrong.')
    return 1

def main():
    """
    Runs some tests.
    """
    lvl = 10
    fail = test_pieri(lvl)
    fail = fail + test_littlewood_richardson_rule(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=="__main__":
    main()
