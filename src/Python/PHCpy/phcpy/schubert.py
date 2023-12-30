"""
Numerical Schubert calculus defines homotopies for enumerative geometry.
"""
from ctypes import c_int, c_double, pointer, create_string_buffer
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

def main():
    """
    Runs some tests.
    """
    lvl = 10
    fail = test_pieri(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=="__main__":
    main()
