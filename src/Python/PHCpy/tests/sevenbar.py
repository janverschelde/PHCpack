"""
Setup of a moving sevenbar example with sympy.
There is one irreducible component of degree three
and six isolated points.
"""
from sympy import var
from cmath import exp
from phcpy.solutions import coordinates, diagnostics
from phcpy.sets import double_laurent_membertest
from phcpy.cascades import double_laurent_top_cascade
from phcpy.cascades import double_laurent_cascade_filter
from phcpy.factor import double_monodromy_breakup

def symbolic_equations():
    """
    Returns the symbolic equations,
    with parameters a1, a2, a3, a4, a5, a6
    b0, b2, b3, b4, b5, and c0, with variables
    t1, t2, t3, t4, and t5. 
    """
    a0, a1, a2, a3, a4, a5, a6 = var('a0, a1, a2, a3, a4, a5, a6')
    b0, b2, b3, b4, b5, c0 = var('b0, b2, b3, b4, b5, c0')
    t1, t2, t3, t4, t5, t6 = var('t1, t2, t3, t4, t5, t6')
    eq1 = a1*t1 + a2*t2 - a3*t3 - a0
    eq2 = b2*t2 + a3*t3 - a4*t4 + a5*t5 - b0
    eq3 = a4*t4 + b5*t5 - a6*t6 - c0
    return [eq1, eq2, eq3]

def generic_problem(eqs):
    """
    Given the symbolic equations in eqs,
    defines the equations for a generic problem,
    as a Laurent polynomial system.
    The system is returned as a list of string representations,
    suitable for input to the solve of phcpy.
    """
    i = complex(0, 1)
    subdict = {a0: 0.7 + 0.2*i, b0: 0.6, c0: 0.5 - 0.5*i, \
        a1: 0.7, a2: 0.8, b2: 0.6 + 0.5*i, a3: 0.4, a4: 0.6, \
        a5: 0.8, b5: 0.4 + 0.3*i, a6: 0.9}
    print(subdict)
    conjugates = {a0: 0.7 - 0.2*i, b0: 0.6, c0: 0.5 + 0.5*i, \
        a1: 0.7, a2: 0.8, b2: 0.6 - 0.5*i, a3: 0.4, a4: 0.6, \
        a5: 0.8, b5: 0.4 - 0.3*i, a6: 0.9}
    print(conjugates)
    result = []
    for equ in eqs:
        pol = equ.subs(subdict)
        result.append(str(pol) + ';')
    for equ in eqs:
        pol = str(equ.subs(conjugates))
        pol = pol.replace('t1', 't1**(-1)')
        pol = pol.replace('t2', 't2**(-1)')
        pol = pol.replace('t3', 't3**(-1)')
        pol = pol.replace('t4', 't4**(-1)')
        pol = pol.replace('t5', 't5**(-1)')
        pol = pol.replace('t6', 't6**(-1)')
        result.append(pol + ';')
    return result

def special_parameters():
    """
    Returns a dictionary with special values for the parameters
    for the Assur7c in Roberts Cognate pattern.
    Before calling this function, the symbolic_equations()
    must have defined the variables for the parameters.
    """
    i = complex(0, 1)
    # start with the independent parameters
    result = {b0: 0.0, c0: 1.2, a2: 0.46, \
        b2: -0.11 + 0.49*i, a5: 0.41}
    theta4 = 0.6 + 0.8*i
    theta3 = exp(1.8*i)
    # add the derived parameters
    result[a3] = result[a5]
    beta = result[b2]/result[a2]
    result[a0] = result[c0]/beta
    result[b5] = result[a5]*beta
    result[a4] = abs(result[b2])
    result[a1] = abs(result[a0] + result[a3]*theta3 - result[a4]*theta4/beta)
    result[a6] = abs(result[a4]*theta4 - result[b5]*theta3-result[c0])
    return result

def conjugates(dic):
    """
    Given on input a dictionary with variables as keys
    and complex numbers as values.
    Returns a dictionary with the same keys,
    but with values replaced by complex conjugates.
    """
    result = {}
    for key in list(dic.keys()):
        result[key] = dic[key].conjugate()
    return result
  
def special_problem(eqs):
    """
    Given the symbolic equations in eqs,
    replaces the parameters with special values.
    """
    pars = special_parameters()
    conj = conjugates(pars)
    result = []
    for equ in eqs:
        pol = equ.subs(pars)
        result.append(str(pol) + ';')
    for equ in eqs:
        pol = str(equ.subs(conj))
        pol = pol.replace('t1', 't1**(-1)')
        pol = pol.replace('t2', 't2**(-1)')
        pol = pol.replace('t3', 't3**(-1)')
        pol = pol.replace('t4', 't4**(-1)')
        pol = pol.replace('t5', 't5**(-1)')
        pol = pol.replace('t6', 't6**(-1)')
        result.append(pol + ';')
    return result

def embed_and_cascade(pols, topdim):
    """
    Computes and solves an embedding at top dimension topdim
    of the Laurent polynomials in pols, before running one
    step in the cascade homotopy.
    Returns the embedded system, the three generic points,
    and the filtered solutions at the end of the cascade.
    """
    (embpols, sols0, sols1) \
        = double_laurent_top_cascade(len(pols), topdim, pols, 1.0e-08)
    print('the top generic points :')
    for (idx, sol) in enumerate(sols0):
        print('Solution', idx+1, ':')
        print(sol)
    print('the nonsolutions :')
    for (idx, sol) in enumerate(sols1):
        print('Solution', idx+1, ':')
        print(sol)
    input('hit enter to continue')
    print('... running cascade step ...')
    (embdown, nsols1, sols2) = double_laurent_cascade_filter(1, embpols, \
        sols1, 1.0e-8)
    filtsols2 = []
    for (idx, sol) in enumerate(nsols1):
        err, rco, res = diagnostics(sol)
        if res < 1.0e-8 and rco > 1.0e-8:
            _, point = coordinates(sol)
            crdpt = [] 
            for pt in point:
                crdpt.append(pt.real)
                crdpt.append(pt.imag)
            onset = double_laurent_membertest(embpols, sols0, 1, crdpt) 
            if not onset:
                filtsols2.append(sol)
    print('... after running the cascade ...')
    for (idx, sol) in enumerate(filtsols2):
        print('Solution', idx+1, ':')
        print(sol)
    print('found %d isolated solutions' % len(filtsols2))
    return (embpols, sols0, filtsols2)

def monodromy_factor(witpols, witpnts):
    """
    Given the polynomials and points in the witness set,
    in witpols and witpnts respectively, applies monodromy
    to factor the witness set.
    """
    print('the enbedded Laurent polynomial system :')
    for pol in witpols:
        print(pol)
    input('hit enter to call the factor method')
    fac = double_monodromy_breakup(witpols, witpnts, 1, \
        islaurent=True, verbose=True)
    print('the factorization :', fac)

def main():
    """
    Defines the equations for the sevenbar problem.
    """
    eqs = symbolic_equations()
    for equ in eqs:
        print(equ)
    T1, T2, T3, T4, T5, T6 = var('T1, T2, T3, T4, T5, T6')
    generic = generic_problem(eqs)
    for equ in generic:
        print(equ)
    from phcpy.solver import solve
    sols = solve(generic)
    print('found', len(sols), 'solutions')
    input('hit enter to continue')
    for sol in sols:
        print(sol)
    special = special_problem(eqs)
    for equ in special:
        print(equ)
    sols = solve(special)
    print('found', len(sols), 'solutions')
    input('hit enter to continue')
    for sol in sols:
        print(sol)
    input('hit enter to continue')
    (embpols, sols0, isosols) = embed_and_cascade(special, 1)
    monodromy_factor(embpols, sols0)
    print('the isolated solutions :')
    for (idx, sol) in enumerate(isosols):
        print('Solution', idx+1, ':')
        print(sol)
   
main()
