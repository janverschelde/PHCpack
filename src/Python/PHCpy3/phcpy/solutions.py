"""
The module solutions exports functions to convert a list of
PHCpack solution strings into Python dictionaries.
The module exports the definition of the class Solution,
as an object-oriented representation of a solution.
"""

def diagnostics(sol):
    r"""
    Extracts the diagnostics (err, rco, res)
    from the PHCpack string solution in *sol* and
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

def map_double(freqtab, nbr):
    """
    On input in freqtab is a list of integer numbers and nbr is a double.
    The list freqtab represents a frequency table of magnitudes.
    The magnitude of the double nbr is mapped into the frequency table.
    The counter in freqtab that will be updated is at position
    floor(-log10(nbr)) 
    """
    if nbr > 1.0:
       freqtab[0] = freqtab[0] + 1
    else:
       from math import log10, floor
       tol = 10.0**(-len(freqtab)+1)
       if nbr < tol:
          freqtab[len(freqtab)-1] = freqtab[len(freqtab)-1] + 1
       else:
          idx = floor(-log10(nbr))
          if idx < 0:
             freqtab[0] = freqtab[0] + 1
          elif idx >= len(freqtab):
             freqtab[len(freqtab)-1] = freqtab[len(freqtab)-1] + 1
          else:
             freqtab[idx] = freqtab[idx] + 1

def condition_tables(sols):
    """
    The input in sols is a list of PHCpack string solutions.
    A condition table is triplet of three frequency tables,
    computed from the diagnostics (err, rco, res) of each solution.
    The i-th entry in each frequency table counts the number of
    doubles x which floor(-log10(x)) mapped to the index i.
    Small numbers are mapped to the right of the table,
    large numbers are mapped to the left of the table.
    """
    errtab = [0 for _ in range(16)]
    rcotab = [0 for _ in range(16)]
    restab = [0 for _ in range(16)]
    for sol in sols:
        (err, rco, res) = diagnostics(sol)
        map_double(errtab, err)
        map_double(rcotab, rco)
        map_double(restab, res)
    return (errtab, rcotab, restab);

def str2complex(scn):
    r"""
    The string *scn* contains a complex number,
    the real and imaginary part separated by spaces.
    On return is the Python complex number.
    """
    from ast import literal_eval
    stripped = scn.strip()
    realimag = stripped.split(' ')
    realpart = realimag[0].replace('E', 'e')
    imagpart = realimag[len(realimag)-1].replace('E', 'e')
    return complex(literal_eval(realpart), literal_eval(imagpart))

def string_complex(scn):
    r"""
    The string *scn* contains a complex number,
    the real and imaginary part separated by spaces.
    On return is the string representation of a complex number,
    in Python format.  The use of this function is for when
    the coordinates are calculated in higher precision.
    """
    stripped = scn.strip()
    realimag = stripped.split(' ')
    realpart = realimag[0].replace('E', 'e')
    imagpart = realimag[len(realimag)-1].replace('E', 'e')
    if imagpart[0] == '-':
        result = '(' + realpart + imagpart + 'j)'
    else:
        result = '(' + realpart + '+' + imagpart + 'j)'
    return result

def coordinates(sol):
    r"""
    Returns the coordinates of the solution
    in the PHCpack solution string *sol*,
    as a tuple of two lists: (names, values).
    The list names contains the strings of the variable names.
    The list values contains the complex values for the
    coordinates of the solution.  The entries in the list
    names correspond to the entries in the list values.
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

def string_coordinates(sol):
    r"""
    Returns the coordinates of the solution
    in the PHCpack solution string *sol*,
    as a tuple of two lists: (names, values).
    For each name in names there is a value in values.
    The list names contains the strings of the variable names.
    The list values contains the values of the coordinates,
    represented as strings.  This function is useful for
    when the coordinates are computed in higher precision.
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
            vals.append(string_complex(xval[1]))
    return (vard, vals)

def endmultiplicity(sol):
    r"""
    Returns the value of t at the end
    and the multiplicity as (t,m)
    for the PHCpack solution string *sol*.
    """
    from ast import literal_eval
    data = sol.split("the solution for t :")
    tstr = data[0]
    line = tstr.split('\n')
    tstr = line[0].split(':')
    tval = str2complex(tstr[1])
    mstr = line[1].split(':')
    mval = literal_eval(mstr[1].lstrip())
    return (tval, mval)

def strsol2dict(sol, precision='d'):
    r"""
    Converts the solution in the string *sol*
    into a dictionary format.
    By default, the precision of the coordinates is assumed
    to be double float ('d' on input).
    If the precision is not 'd', then the coordinates of the solution
    are returned as Python complex number string representations.
    """
    result = {}
    (tval, mult) = endmultiplicity(sol)
    result['t'] = tval
    result['m'] = mult
    (err, rco, res) = diagnostics(sol)
    result['err'] = err
    result['rco'] = rco
    result['res'] = res
    if precision == 'd':
        (var, val) = coordinates(sol)
    else:
        (var, val) = string_coordinates(sol)
    for idx, name in enumerate(var):
        result[name] = val[idx]
    return result

def formdictlist(sols, precision='d'):
    r"""
    Given in *sols* is a list of strings.
    Each string in *sols* represents a solution,
    in the PHCpack format.
    On return is the list of dictionaries.
    Each dictionary in the list of return
    stores each solution of the list *sols*
    in the dictionary format.
    By default, the precision of the coordinates is assumed
    to be double float ('d' on input).
    If the precision is not 'd', then the coordinates of the solution
    are returned as Python complex number string representations.
    """
    return [strsol2dict(sol, precision) for sol in sols]

def variables(dsol):
    r"""
    Given the dictionary format of a solution in *dsol*,
    returns the list of variables.
    """
    solkeys = list(dsol.keys())
    solkeys.remove('t')
    solkeys.remove('m')
    solkeys.remove('err')
    solkeys.remove('rco')
    solkeys.remove('res')
    return solkeys

def numerals(dsol):
    r"""
    Given the dictionary format of a solution *dsol*,
    returns the list of numeric values of the variables in the solution.
    """
    names = variables(dsol)
    return [dsol[name] for name in names]

def evaluate_polynomial(pol, dsol):
    r"""
    Evaluates the polynomial *pol* at the solution
    dictionary *dsol* by string substitution.
    """
    j = complex(0, 1) # needed for eval(result) in test()
    varsd = variables(dsol)
    rpol = pol
    rpol = rpol.replace('i', 'j')
    rpol = rpol.replace('E', 'e')
    rpol = rpol.replace('^', '**')
    for varname in varsd:
        xre = '%+.17f' % dsol[varname].real
        xim = '%+.17f' % dsol[varname].imag
        value = '(' + xre + xim + 'j)'
        rpol = rpol.replace(varname, value)
    result = rpol[:-1]
    return eval(result)

def evaluate(pols, dsol):
    r"""
    Evaluates a list of polynomials given as string in *pols*
    at the solution in dictionary format in *dsol*.
    """
    result = []
    for pol in pols:
        result.append(evaluate_polynomial(pol, dsol))
    return result

def make_solution(names, values, \
    err=0.0, rco=1.0, res=0.0, tval=0, multiplicity=1):
    r"""
    Makes the string representation in PHCpack format
    with in *names* a list of strings for the variables names
    and in *values* a list of (complex) values for the coordinates.
    For example:

    s = make_solution(['x','y'],[(1+2j), 3])

    returns the string s to represent the solution with
    coordinates (1+2j) and 3 for the variables x and y.
    The imaginary unit is the Python j instead of i.
    Other arguments for this function are

    1. *err* is the magnitude of an update, or the forward error,

    2. *rco* is the estimate for the inverse condition number,

    3. *res* is the value of the residual, or backward error,

    4. *tval* is the value for the continuation parameter t,

    5. *multiplicity* is the multiplicity of the solution.

    For those above arguments, default values are provided.
    Applying the function coordinates on the result of
    make_solution returns the tuple of arguments given
    on input to **make_solution()**.
    """
    if isinstance(tval, complex):
        (tre, tim) = (tval.real, tval.imag)
    elif isinstance(tval, float):
        (tre, tim) = (tval, 0.0)
    elif isinstance(tval, int):
        (tre, tim) = (tval, 0.0)
    else:
        print('wrong type for the value of the continuation parameter t')
        return ""
    result = 't : %.15E  %.15E\n' % (tre, tim)
    mstr = 'm : %d\n' % multiplicity
    result = result + mstr
    result = result + 'the solution for t :\n'
    for idx, name in enumerate(names):
        result = result + ' ' + name + ' : '
        if isinstance(values[idx], complex):
            c_re = '%.15E' % values[idx].real
            c_im = '%.15E' % values[idx].imag
            result = result + c_re + '  ' + c_im + '\n'
        elif isinstance(values[idx], float):
            flt = '%.15E' % values[idx]
            result = result + flt + '  ' + '0.0' + '\n'
        elif isinstance(values[idx], int):
            i = '%.15E' % values[idx]
            result = result + i + '  ' + '0.0' + '\n'
        else:
            print('wrong type for coordinate value')
            return result
    lastline = '== err :  %.3E = rco :  %.3E = res :  %.3E =' % (err, rco, res)
    result = result + lastline
    return result

def is_real(sol, tol):
    r"""
    Returns True if the solution in *sol* is real with respect
    to the given tolerance *tol*: if the absolute value of the imaginary
    part of all coordinates are less than *tol*.
    """
    (_, vals) = coordinates(sol)
    for value in vals:
        if abs(value.imag) > tol:
            return False
    return True

def filter_real(sols, tol, oper):
    r"""
    Filters the real solutions in *sols*.
    The input parameters are

    1. *sols* is a list of solution strings in PHCpack format,

    2. *tol* is the tolerance on the absolute value of the
       imaginary parts of the coordinates of the solution.

    3. *oper* is either 'select' or 'remove'

       if *oper* == 'select' then solutions that are considered real
       are selected and in the list on return,

       if *oper* == 'remove' then solutions that are considered real
       are in the list on return.
    """
    result = []
    for sol in sols:
        isreal = is_real(sol, tol)
        if oper == 'select':
            if isreal:
                result.append(sol)
        if oper == 'remove':
            if not isreal:
                result.append(sol)
    return result

def filter_regular(sols, tol, oper):
    r"""
    Filters solutions in *sols* for the estimate of
    the inverse of the condition number.
    The input parameters are

    1. *sols* is a list of solution strings in PHCpack format,

    2. *tol* is the tolerance on the value for the estimate rco
       for the inverse of the condition number to decide whether
       a solution is singular (if rco < *tol*) or not.

    3. *oper* is either 'select' or 'remove'

       if *oper* == 'select' then solutions with value rco > *tol*
       are selected and in the list on return,

       if *oper* == 'remove' then solutions with value rco <= *tol*
       are in the list on return.
    """
    result = []
    for sol in sols:
        rco = diagnostics(sol)[1]
        if oper == 'select':
            if rco > tol:
                result.append(sol)
        if oper == 'remove':
            if rco <= tol:
                result.append(sol)
    return result

def filter_zero_coordinates(sols, varname, tol, oper):
    r"""
    Filters the solutions in *sols* for variables
    that have a value less than the tolerance.
    The input parameters are

    1. *sols* is a list of solution strings in PHCpack format,

    2. *varname* is a string with the name of the variable,

    3. *tol* is the tolerance to decide whether a complex
       number equals zero or not, and

    4. *oper* is either 'select' or 'remove'

       if *oper* == 'select' then solutions with value for the
       variable v that is less than *tol* are selected and
       in the list on return,

       if *oper* == 'remove' then solutions with value for the
       variable v that is less than *tol* are removed and
       not in the list on return.
    """
    result = []
    for sol in sols:
        dsol = strsol2dict(sol)
        if oper == 'select':
            if abs(dsol[varname]) < tol:
                result.append(sol)
        if oper == 'remove':
            if abs(dsol[varname]) >= tol:
                result.append(sol)
    return result

def is_vanishing(sol, tol):
    r"""
    Given in *sol* is a solution string and
    *tol* is the tolerance on the residual.
    Returns True if the residual of *sol*
    is less than or equal to *tol*.
    Returns False otherwise.
    """
    dgn = diagnostics(sol)
    return dgn[2] <= tol

def filter_vanishing(sols, tol):
    r"""
    Returns the list of solutions in *sols* that have a residual
    less than or equal to the given tolerance in *tol*.
    """
    result = []
    for sol in sols:
        if is_vanishing(sol, tol):
            result.append(sol)
    return result

def test_functions():
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

class Solution(object):
    """
    Wraps the functions on solution strings.
    """
    def __init__(self, sol):
        """
        A solution is constructed from the
        string representation in PHCpack format,
        or otherwise from a dictionary with keys
        the names of the variables and values
        the complex numbers for the corresponding
        coordinates of the solution.
        """
        if isinstance(sol, str):
            self.dict = strsol2dict(sol)
        elif isinstance(sol, dict):
            self.dict = sol
            if 'err' not in self.dict:
                self.dict['err'] = 0.0
            if 'rco' not in self.dict:
                self.dict['rco'] = 1.0
            if 'res' not in self.dict:
                self.dict['res'] = 0.0
            if 'm' not in self.dict:
                self.dict['m'] = 1
            if 't' not in self.dict:
                self.dict['t'] = complex(0.0)
        else:
            print('Wrong argument type, provide string or dict.')

    def __str__(self):
        """
        Returns the string representation of a solution.
        """
        result = make_solution(variables(self.dict), numerals(self.dict), \
            err=self.dict['err'], rco=self.dict['rco'], res=self.dict['res'], \
            tval=self.dict['t'], multiplicity=self.dict['m'])
        return result

    def __repr__(self):
        """
        Defines the representation as the string representation.
        """
        return str(self)

    def coordinates(self):
        """
        Returns the values of the coordinates of the solution,
        as a tuple of variable names and corresponding values.
        """
        return (variables(self.dict), numerals(self.dict))

    def dictionary(self):
        """
        Returns the dictionary format of the solution.
        """
        return self.dict

    def variables(self):
        """
        Returns the variable names of the coordinates.
        """
        return variables(self.dict)

    def numerals(self):
        """
        Returns the numerical values of the coordinates.
        """
        return numerals(self.dict)

    def diagnostics(self):
        """
        Returns the numerical diagnostics.
        """
        return (self.dict['err'], self.dict['rco'], self.dict['res'])

    def multiplicity(self):
        """
        Returns the multiplicity.
        """
        return self.dict['m']

    def timevalue(self):
        """
        Returns the value of the continuation parameter.
        """
        return self.dict['t']

def test_class():
    """
    Tests the methods in the class Solution.
    """
    from phcpy import solver
    pols = solver.random_trinomials()
    sols = solver.solve(pols)
    print('sols[0] :')
    print(sols[0])
    s = Solution(sols[0])
    print('the first solution :')
    print(s)
    print('its coordinates :\n', s.coordinates())
    print('its variable names :\n ', s.variables())
    print('its numerical values :\n', s.numerals())
    print('its diagnostics :\n', s.diagnostics())
    print('its continuation parameter :', s.timevalue())
    print('its multiplicity :', s.multiplicity())
    print('the dictionary :\n', s.dictionary())
    mysol = Solution({'x': complex(1,2), 'y': complex(-7,0)})
    print('my solution :')
    print(mysol)
    (errtab, rcotab, restab) = condition_tables(sols)
    print('The frequency tables for err, rco, and res:')
    print(errtab)
    print(rcotab)
    print(restab)

if __name__ == "__main__":
    # test_functions()
    test_class()
