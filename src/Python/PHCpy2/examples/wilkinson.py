"""
The perfidious polynomial of Wilkinson of degree d is the product of the
first d factor x - k, for k ranging from 1 to d.
"""
def wilkpol(deg):
    """
    Returns the string representation of the Wilkinson polynomial
    of degree deg > 0.
    """
    result = "(x - 1)"
    for k in range(2, deg+1):
        factor = '*(x - %d)' % k
        result = result + factor
    result = result + ';'
    return result

def solve(pol):
    """
    Applies the blackbox solver of phcpy to the polynomial pol.
    """
    from phcpy.solver import solve
    try:
        sols = solve([pol], verbose=False)
        print sols
        return sols
    except:
        print 'an exception happened in the blackbox solver'
        return []

def track(pol):
    """
    Applies the series-Pade tracker to solve the polynomial,
    either in a step-by-step manner or without this interaction.
    """
    from phcpy.solver import total_degree_start_system as tdss
    from phcpy.curves import tune_homotopy_continuation_parameters
    from phcpy.curves import standard_next_track, standard_track
    (startpol, startsols) = tdss([pol])
    print 'the start polynomial :'
    print startpol
    print 'the start solutions :'
    for sol in startsols:
        print sol
    tune_homotopy_continuation_parameters()
    ans = raw_input("Interactive, step-by-step track ? (y/n) ");
    if(ans == 'y'):
        sols = standard_next_track([pol], startpol, startsols, True)
    else:
        sols = standard_track([pol], startpol, startsols, \
                              filename="/tmp/outoftrack", verbose=True)
    print 'the computed solutions :'
    for sol in sols:
        print sol

def main():
    """
    Prompts the user for a degree and then computes the roots
    of the Wilkinson polynomial.
    """
    deg = int(raw_input("Give the degree : "))
    pol = wilkpol(deg)
    print pol
    ans = raw_input("Apply the blackbox solver ? (y/n) ")
    if(ans == 'y'):
        sols = solve(pol)
        for sol in sols:
            print sol
    track(pol)

main()
