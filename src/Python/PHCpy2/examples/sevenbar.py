"""
Setup of a moving sevenbar example with sympy.
"""
from sympy import var

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
    subdict = {a0:0.7+.2*i, b0:0.6, c0:0.5-0.5*i, \
        a1:0.7, a2:0.8, b2:0.6+0.5*i, a3:0.4, a4:0.6, \
        a5:0.8, b5:0.4+0.3*i, a6:0.9}
    print subdict
    conjugates = {a0:0.7-.2*i, b0:0.6, c0:0.5+0.5*i, \
        a1:0.7, a2:0.8, b2:0.6-0.5*i, a3:0.4, a4:0.6, \
        a5:0.8, b5:0.4-0.3*i, a6:0.9}
    print conjugates
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

def main():
    """
    Defines the equations for the sevenbar problem.
    """
    eqs = symbolic_equations()
    for equ in eqs:
        print equ
    T1, T2, T3, T4, T5, T6 = var('T1, T2, T3, T4, T5, T6')
    generic = generic_problem(eqs)
    for equ in generic:
        print equ
    from phcpy.solver import solve
    sols = solve(generic)
    for sol in sols:
        print sol
   
main()
