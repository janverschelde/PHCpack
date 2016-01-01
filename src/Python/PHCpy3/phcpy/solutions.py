"""
module to convert list of PHCpack solution strings into dictionaries
"""

def diagnostics(sol):
    """
    Extracts the diagnostics (err,rco,res)
    from the PHCpack string solution in sol and
    returns a triplet of three floats.
    """
    from ast import literal_eval
    banner = sol.split("==")
    data = banner[1]
    diag = data.split("=")
    str_err = diag[0].split(':')
    str_rco = diag[1].split(':')
    str_res = diag[2].split(':')
    val_err = literal_eval(str_err[1].lstrip())
    val_rco = literal_eval(str_rco[1].lstrip())
    val_res = literal_eval(str_res[1].lstrip())
    # print 'err =', val_err, 'rco =', val_rco, 'res =', val_res
    return (val_err, val_rco, val_res)

def str2complex(scn):
    """
    The string scn contains a complex number,
    the real and imaginary part separated by spaces.
    On return is the Python complex number.
    """
    from ast import literal_eval
    stripped = scn.strip()
    realimag = stripped.split(' ')
    realpart = realimag[0].replace('E', 'e')
    imagpart = realimag[len(realimag)-1].replace('E', 'e')
    return complex(literal_eval(realpart), literal_eval(imagpart))

def coordinates(sol):
    """
    Returns the coordinates of the solution
    in the PHCpack solution string sol.
    """
    banner = sol.split("==")
    data = banner[0]
    nums = data.split("the solution for t :")
    firstnums = nums[1]
    lines = firstnums.split('\n')
    vard = []
    vals = []
    for line in lines:
        if line != '':
            xval = line.split(':')
            vard.append(xval[0].strip())
            vals.append(str2complex(xval[1]))
    return (vard, vals)

def endmultiplicity(sol):
    """
    Returns the value of t at the end
    and the multiplicity as (t,m)
    for the PHCpack solution string sol.
    """
    from ast import literal_eval
    data = sol.split("the solution for t :")
    tstr = data[0]
    line = tstr.split('\n')
    tstr = line[0].split(':')
    tval = tstr[1]
    mstr = line[1].split(':')
    mval = literal_eval(mstr[1].lstrip())
    return (tval, mval)

def strsol2dict(sol):
    """
    Converts the solution in the string sol
    into a dictionary format.
    """
    result = {}
    (tval, mult) = endmultiplicity(sol)
    result['t'] = tval
    result['m'] = mult
    (err, rco, res) = diagnostics(sol)
    result['err'] = err
    result['rco'] = rco
    result['res'] = res
    (var, val) = coordinates(sol)
    for i in range(len(var)):
        result[var[i]] = val[i]
    return result

def variables(dsol):
    """
    Given the dictionary format of a solution,
    returns the list of variables.
    """
    solkeys = list(dsol.keys())
    solkeys.remove('t')
    solkeys.remove('m')
    solkeys.remove('err')
    solkeys.remove('rco')
    solkeys.remove('res')
    return solkeys

def evaluate_polynomial(pol, dsol):
    """
    Evaluates the polynomial pol at the solution
    dictionary dsol by string substitution.
    """
    from ast import literal_eval
    varsd = variables(dsol)
    rpol = pol
    rpol = rpol.replace('i', 'j')
    rpol = rpol.replace('E', 'e')
    rpol = rpol.replace('^', '**')
    j = complex(0, 1) # j is used in literal_eval(result)
    for varname in varsd:
        xre = '%+.17f' % dsol[varname].real
        xim = '%+.17f' % dsol[varname].imag
        value = '(' + xre + xim + '*j)'
        rpol = rpol.replace(varname, value)
    result = rpol[:-1]
    return eval(result)

def evaluate(pols, dsol):
    """
    Evaluates a list of polynomials given as string in pols
    at the solution in dictionary format in dsol.
    """
    result = []
    for pol in pols:
        result.append(evaluate_polynomial(pol, dsol))
    return result

def make_solution(sol, vals):
    """
    Makes the string representation in PHCpack format
    with in sol a list of strings for the variables names
    and in vals a list of complex values for the coordinates.
    For example: s = make_solution(['x','y'],[1,2])
    returns the string s to represent the solution with
    coordinates 1 and 2 for the variables x and y.
    Applying the function coordinates on the result of
    make_solution returns the tuple of arguments given
    on input to make_solution.
    """
    result = 't : 0.0 0.0\nm : 1\n'
    result = result + 'the solution for t :\n'
    for k in range(len(sol)):
        result = result + ' ' + sol[k] + ' : '
        if isinstance(vals[k], complex):
            c_re = '%.15E' % vals[k].real
            c_im = '%.15E' % vals[k].imag
            result = result + c_re + '  ' + c_im + '\n'
        elif isinstance(vals[k], float):
            flt = '%.15E' % vals[k]
            result = result + flt + '  ' + '0.0' + '\n'
        elif isinstance(vals[k], int):
            i = '%.15E' % vals[k]
            result = result + i + '  ' + '0.0' + '\n'
    result = result + '== err : 0.0 = rco : 1.0 = res : 0.0 =='
    return result

def filter_regular(sols, tol, oper):
    """
    Filters solutions in sols for the estimate of
    the inverse of the condition number.
    The input parameters are
    (1) sols is a list of solution strings in PHCpack format,
    (2) tol is the tolerance on the value for the estimate rco
    for the inverse of the condition number to decide whether
    a solution is singular (if rco < tol) or not.
    (3) oper is either 'select' or 'remove'
    if oper == 'select' then solutions with value rco > tol
    are selected and in the list on return,
    if oper == 'remove' then solutions with value rco <= tol
    are in the list on return.
    """
    result = []
    for sol in sols:
        rco = diagnostics(sol)[1]
        if(oper == 'select'):
            if(rco > tol):
                result.append(sol)
        if(oper == 'remove'):
            if(rco <= tol):
                result.append(sol)
    return result

def filter_zero_coordinates(sols, varname, tol, oper):
    """
    Filters the solutions in sols for variables
    that have a value less than the tolerance.
    The input parameters are
    (1) sols is a list of solution strings in PHCpack format,
    (2) varname is a string with the name of the variable,
    (3) tol is the tolerance to decide whether a complex
    number equals zero or not, and
    (4) oper is either 'select' or 'remove'
    if oper == 'select' then solutions with value for the
    variable v that is less than tol are selected and
    in the list on return,
    if oper == 'remove' then solutions with value for the
    variable v that is less than tol are removed and
    not in the list on return.
    """
    result = []
    for sol in sols:
        dsol = strsol2dict(sol)
        if(oper == 'select'):
            if(abs(dsol[varname]) < tol):
                result.append(sol)
        if(oper == 'remove'):
            if(abs(dsol[varname]) >= tol):
                result.append(sol)
    return result

def test():
    """
    Generates a random trinomial system,
    solves it, converts the solutions,
    and then sums the multiplicities.
    """
    from phcpy import solver
    pols = solver.random_trinomials()
    sols = solver.solve(pols)
    dsols = [strsol2dict(sol) for sol in sols]
    mult = 0
    s0d = strsol2dict(sols[0])
    print('variables :', variables(s0d))
    print(evaluate(pols, s0d))
    for sol in dsols:
        mult = mult + sol['m']
    print('sum of multiplicities :', mult)
    print('sum of values at the solutions :')
    for sol in dsols:
        print(sum(evaluate(pols, sol)))

if __name__ == "__main__":
    test()
