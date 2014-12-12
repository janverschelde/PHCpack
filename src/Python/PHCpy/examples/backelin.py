"""
Backelin's Lemma states that the cyclic n-roots problem has an m-dimensional
solution set for n = L*m^2, where L is no multiple of k^2, m >= 2.
This script provides an exact representation of the solution set,
following the tropical formulation with m-1 free parameters.
"""

def constants(m, L):
    """
    Returns the tuple of constants u and gamma
    for n = m*L**2.
    """
    from cmath import exp, pi
    alpha = m*(m*L - 1)
    beta = alpha % 2
    i = complex(0, 1)
    u = exp(i*2*pi/(m*L))
    gamma = exp(i*pi*beta/(m*L))
    return (u, gamma)

def component(m, L, u, gamma):
    """
    Returns a list of strings, where the k-th string
    stores the representation of the k-th coordinate.
    The m-1 parameters are t0, t1, .., tm-2.
    """
    result = []
    for ell in range(L):
        for k in range(m):
            if(ell % 2 == 1):
                s = '-' + str(u) + '**' + str(k)
            else:
                s = '+' + str(u) + '**' + str(k)
            for i in range(m-1):
                s = s + '*t' + str(i)
                result.append(s)
            if(ell % 2 == 1):
                s = '-' + str(gamma) + '*' + str(u) + '**' + str(k)
            else:
                s = '+' + str(gamma) + '*' + str(u) + '**' + str(k)
            for i in range(m-1):
                e = -m + 1 + i
                s = s + '*t' + str(i) + '**(' + str(e) + ')'
            result.append(s)
    return result

def sample(m, roots):
    """
    Generates m-1 random complex numbers as parameters
    and evaluates the representation for the roots.
    """
    from random import uniform
    from cmath import exp, pi
    i = complex(0, 1)
    d = globals()
    for k in range(m-1):
        var = 't' + str(k)
        u = uniform(0, 2*pi)
        d[var] = exp(i*u)
    result = []
    for r in roots:
        result.append(eval(r))
    return result

def evaluate(point):
    """
    Evaluates the point in the cyclic n-roots problem,
    where n = len(point).
    """
    from phcpy.families import cyclic
    n = len(point)
    f = cyclic(n)
    d = globals()
    for k in range(n):
        var = 'x' + str(k)
        d[var] = point[k]
    for k in range(n):
        var = 'x' + str(k)
        print var, '=', d[var]
    result = []
    for r in f:
        result.append(eval(r[0:-1]))
    return result

def main():
    """
    Prompts the user for the parameters m and L.
    """
    m = input("Give the m in n = m**2*L : ")
    L = input("Give the L in n = m**2*L : ")
    n = m**2*L
    print 'The dimension n = %d**2*%d = %d.' % (m, L, n)
    u, gamma = constants(m, L)
    print 'u =', u
    print 'gamma =', gamma
    coords = component(m, L, u, gamma)
    print 'coordinates of a cyclic %d-roots component :' % n
    for c in coords:
        print c
    point = sample(m, coords)
    print 'a random point:'
    for p in point:
        print p
    residual = evaluate(point)
    print 'the residual :'
    for r in residual:
        print r
    print 'the sum :', sum([abs(r) for r in residual])

main()
