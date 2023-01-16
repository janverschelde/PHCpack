"""
Exports the definition of polynomial systems.
"""

def setDoublePolynomial(idx, nvr, pol, vrblvl=0):
    """
    Sets the polynomial by the string pol, in double precision,
    with a number of variables no more than nvr, at index idx.
    The index starts at one, not at zero.
    This function does not set the dimension,
    which must be set to a value at least idx.
    """
    from ctypes import c_int, c_double, pointer
    from version import getPHCmod, int4a2nbr, str2int4a
    if(vrblvl > 0):
        print("-> setDoublePolynomial, pol = ", pol)
    modPHCpack = getPHCmod()
    f = modPHCpack._ada_use_c2phc
    vrb = (vrblvl > 0)
    a = int4a2nbr([len(pol), nvr, idx], vrb)
    b = str2int4a(pol)
    c = pointer(c_double(0.0))
    v = c_int(vrblvl)
    r = f(76, a, b, c, v)
    if(vrblvl > 0):
        print('-> setDoublePolynomial return value :', r)
    return r

def setDoubleSystem(nvr, pols, vrblvl=0):
    """
    Sets the system defines by the strings in pols, in double precision,
    with a number of variables no more than nvr.
    The dimension of the system is set to len(pols).
    Returns the sum of the return values of setDoublePolynomial.
    """
    from dimension import setDoubleDimension
    if(vrblvl > 0):
        print('-> setDoubleSystem, nvr =', nvr)
    dim = len(pols)
    setDoubleDimension(dim, vrblvl)
    result = 0
    for k in range(dim):
        if(vrblvl > 0):
            print('-> setDoubleSystem, pols[%d] = %s' % (k, pols[k]))
        r = setDoublePolynomial(k+1, nvr, pols[k], vrblvl)
        if(vrblvl > 0):
            print('-> setDoublePolynomial returns', r)
        result = result + r
    return result

def getDoublePolynomial(idx, vrblvl=0):
    """
    Returns the string representation of the polynomial
    in double precision, at index idx.
    """
    from ctypes import c_int, c_double, pointer, create_string_buffer
    from version import getPHCmod, int4a2str
    if(vrblvl > 0):
        print("-> getDoublePolynomial idx = ", idx)
    modPHCpack = getPHCmod()
    f = modPHCpack._ada_use_c2phc
    a = pointer(c_int(idx))
    b = pointer(c_int(0))
    c = pointer(c_double(0.0))
    v = c_int(vrblvl)
    r = f(600, a, b, c, v)
    if(vrblvl > 0):
        print('-> the return value of getting the size :', r)
        print('-> size of the polynomial :', b[0])
    szd = b[0]
    poldata = create_string_buffer(szd)
    r = f(67, a, poldata, c, v)
    if(vrblvl > 0):
        print('-> the return value of getting the polynomial :', r)
    result = int4a2str(poldata, True)
    return result

def getDoubleSystem(vrblvl=0):
    """
    Returns the string representation of the polynomials
    in double precision.
    """
    from dimension import getDoubleDimension
    dim = getDoubleDimension(vrblvl)
    if(vrblvl > 0):
        print('-> getDoubleSystem, dim =', dim)
    result = []
    for k in range(1, dim+1):
        pol = getDoublePolynomial(k, vrblvl)
        print('-> getDoubleSystem, pol =', pol)
        result.append(pol)
    return result

def showPolynomial():
    """
    Tests the setting and getting of a polynomial.
    """
    from dimension import setDoubleDimension, getDoubleDimension
    lvl = 10
    setDoubleDimension(2, lvl)
    dim = getDoubleDimension(lvl)
    print('the dimension :', dim)
    org = "x*y - 1;"
    idx = 1
    setDoublePolynomial(idx, dim, org, lvl) 
    pol = getDoublePolynomial(idx, lvl)
    print('the retrieved polynomial :', pol)

def showSystem():
    """
    Tests the setting and getting of a system.
    """
    lvl = 10
    polynomials = ["x^3 + 2*x*y - 1;", "x + y - 1;"]
    r = setDoubleSystem(2, polynomials, lvl)
    pols = getDoubleSystem(lvl)
    print('the retrieved polynomials :')
    for pol in pols:
        print(pol)

if __name__=="__main__":
    # showPolynomial()
    showSystem()
