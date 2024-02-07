"""
Exports functions on solutions.
"""
from math import log10, floor
from ast import literal_eval
from ctypes import c_int32, c_double, pointer, create_string_buffer
from phcpy.version import get_phcfun, int4a2nbr, int4a2str, str2int4a

def diagnostics(sol, vrblvl=0):
    r"""
    Extracts the diagnostics (err, rco, res)
    from the PHCpack string solution in *sol* and
    returns a triplet of three floats.
    """
    if vrblvl > 0:
        print('in diagnostics, sol :')
        print(sol)
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

def map_double(freqtab, nbr, vrblvl=0):
    """
    On input in freqtab is a list of integer numbers and nbr is a double.
    The list freqtab represents a frequency table of magnitudes.
    The magnitude of the double nbr is mapped into the frequency table.
    The counter in freqtab that will be updated is at position
    floor(-log10(nbr)) 
    """
    if vrblvl > 0:
        print('in map_double, nbr :', nbr)
        print('freqtab :', freqtab)
    if nbr > 1.0:
        freqtab[0] = freqtab[0] + 1
    else:
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

def condition_tables(sols, vrblvl=0):
    """
    The input in sols is a list of PHCpack string solutions.
    A condition table is triplet of three frequency tables,
    computed from the diagnostics (err, rco, res) of each solution.
    The i-th entry in each frequency table counts the number of
    doubles x which floor(-log10(x)) mapped to the index i.
    Small numbers are mapped to the right of the table,
    large numbers are mapped to the left of the table.
    """
    if vrblvl > 0:
        print('in condition_tables, sols :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    errtab = [0 for _ in range(16)]
    rcotab = [0 for _ in range(16)]
    restab = [0 for _ in range(16)]
    for sol in sols:
        (err, rco, res) = diagnostics(sol)
        map_double(errtab, err)
        map_double(rcotab, rco)
        map_double(restab, res)
    return (errtab, rcotab, restab)

def str2complex(scn, vrblvl=0):
    r"""
    The string *scn* contains a complex number,
    the real and imaginary part separated by spaces.
    On return is the Python complex number.
    """
    if vrblvl > 0:
        print('in str2complex, scn :', scn)
    stripped = scn.strip()
    realimag = stripped.split(' ')
    realpart = realimag[0].replace('E', 'e')
    imagpart = realimag[len(realimag)-1].replace('E', 'e')
    return complex(literal_eval(realpart), literal_eval(imagpart))

def string_complex(scn, vrblvl=0):
    r"""
    The string *scn* contains a complex number,
    the real and imaginary part separated by spaces.
    On return is the string representation of a complex number,
    in Python format.  The use of this function is for when
    the coordinates are calculated in higher precision.
    """
    if vrblvl > 0:
        print('in string_complex, scn :', scn)
    stripped = scn.strip()
    realimag = stripped.split(' ')
    realpart = realimag[0].replace('E', 'e')
    imagpart = realimag[len(realimag)-1].replace('E', 'e')
    if imagpart[0] == '-':
        result = '(' + realpart + imagpart + 'j)'
    else:
        result = '(' + realpart + '+' + imagpart + 'j)'
    return result

def coordinates(sol, vrblvl=0):
    r"""
    Returns the coordinates of the solution
    in the PHCpack solution string *sol*,
    as a tuple of two lists: (names, values).
    The list names contains the strings of the variable names.
    The list values contains the complex values for the
    coordinates of the solution.  The entries in the list
    names correspond to the entries in the list values.
    """
    if vrblvl > 0:
        print('in coordinates, sol :')
        print(sol)
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

def string_coordinates(sol, vrblvl=0):
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
    if vrblvl > 0:
        print('in string_coordinates, sol:')
        print(sol)
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

def endmultiplicity(sol, vrblvl=0):
    r"""
    Returns the value of t at the end
    and the multiplicity as (t,m)
    for the PHCpack solution string *sol*.
    """
    if vrblvl > 0:
        print('in endmultiplicity, sol :')
        print(sol)
    data = sol.split("the solution for t :")
    tstr = data[0]
    line = tstr.split('\n')
    tstr = line[0].split(':')
    tval = str2complex(tstr[1])
    mstr = line[1].split(':')
    mval = literal_eval(mstr[1].lstrip())
    return (tval, mval)

def strsol2dict(sol, precision='d', vrblvl=0):
    r"""
    Converts the solution in the string *sol*
    into a dictionary format.
    By default, the precision of the coordinates is assumed
    to be double float ('d' on input).
    If the precision is not 'd', then the coordinates of the solution
    are returned as Python complex number string representations.
    """
    if vrblvl > 0:
        print('in strsol2dict, precision :', precision, 'sol :')
        print(sol)
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

def formdictlist(sols, precision='d', vrblvl=0):
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
    if vrblvl > 0:
        print('in formdictlist, precision :', precision)
        print('sols :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    return [strsol2dict(sol, precision, vrblvl) for sol in sols]

def variables(dsol, vrblvl=0):
    r"""
    Given the dictionary format of a solution in *dsol*,
    returns the list of variables.
    """
    if vrblvl > 0:
        print('in veriables, dsol :')
        print(dsol)
    solkeys = list(dsol.keys())
    solkeys.remove('t')
    solkeys.remove('m')
    solkeys.remove('err')
    solkeys.remove('rco')
    solkeys.remove('res')
    return solkeys

def numerals(dsol, vrblvl=0):
    r"""
    Given the dictionary format of a solution *dsol*,
    returns the list of numeric values of the variables in the solution.
    """
    if vrblvl > 0:
        print('in numerals, dsol :')
        print(dsol)
    names = variables(dsol, vrblvl)
    return [dsol[name] for name in names]

def evaluate_polynomial(pol, dsol, vrblvl=0):
    r"""
    Evaluates the polynomial *pol* at the solution
    dictionary *dsol* by string substitution.
    """
    if vrblvl > 0:
        print('in evaluate_polynomial, pol :')
        print(pol)
        print('dsol :')
        print(dsol)
    j = complex(0, 1) # needed for eval(result) in test()
    varsd = variables(dsol, vrblvl)
    rpol = pol
    rpol = rpol.replace('i', 'j')
    rpol = rpol.replace('E', 'e')
    rpol = rpol.replace('^', '**')
    for varname in varsd:
        if isinstance(dsol[varname], complex):
            xre = f"{dsol[varname].real:+.17f}"
            xim = f"{dsol[varname].imag:+.17f}"
            value = '(' + xre + xim + 'j)'
        else:
            value = dsol[varname]
        rpol = rpol.replace(varname, value)
    result = rpol[:-1]
    return eval(result)

def evaluate(pols, dsol, vrblvl=0):
    r"""
    Evaluates a list of polynomials given as string in *pols*
    at the solution in dictionary format in *dsol*.
    """
    if vrblvl > 0:
        print('in evaluate, pols :')
        for pol in pols:
            print(pol)
        print('dsol :')
        print(dsol)
    result = []
    for pol in pols:
        result.append(evaluate_polynomial(pol, dsol, vrblvl))
    return result

def verify(pols, sols, vrblvl=0):
    r"""
    Verifies whether the solutions in *sols*
    satisfy the polynomials of the system in *pols*.
    Returns the sum of the absolute values of the 
    residuals in all polynomials.
    """
    if vrblvl > 0:
        print('in verify, pols :')
        for pol in pols:
            print(pol)
        print('sols :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    dictsols = [strsol2dict(sol, vrblvl=vrblvl) for sol in sols]
    checksum = 0
    for (idx, sol) in enumerate(dictsols):
        eva = evaluate(pols, sol, vrblvl-1)
        sumeval = sum(abs(nbr) for nbr in eva)
        if vrblvl > 0:
            print('sum at solution', idx, ':', sumeval)
        checksum = checksum + sumeval
    if vrblvl > 0:
        print('the total check sum :', checksum)
    return checksum.real

def make_solution(names, values, \
    err=0.0, rco=1.0, res=0.0, tval=0, multiplicity=1, vrblvl=0):
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
    if vrblvl > 0:
        print('in make_solution, names :', names)
        print('values :', values)
    if isinstance(tval, complex):
        (tre, tim) = (tval.real, tval.imag)
    elif isinstance(tval, float):
        (tre, tim) = (tval, 0.0)
    elif isinstance(tval, int):
        (tre, tim) = (tval, 0.0)
    else:
        print('wrong type for the value of the continuation parameter t')
        return ""
    result = f"t : {tre:.15E} {tim:.15E}\n"
    mstr = f"m : {multiplicity}\n"
    result = result + mstr
    result = result + 'the solution for t :\n'
    for idx, name in enumerate(names):
        result = result + ' ' + name + ' : '
        if isinstance(values[idx], complex):
            c_re = f"{values[idx].real:.15E}"
            c_im = f"{values[idx].imag:.15E}"
            result = result + c_re + '  ' + c_im + '\n'
        elif isinstance(values[idx], float):
            flt = f"{values[idx]:.15E}"
            result = result + flt + '  ' + '0.0' + '\n'
        elif isinstance(values[idx], int):
            i = f"{values[idx]:.15E}"
            result = result + i + '  ' + '0.0' + '\n'
        else:
            print('wrong type for coordinate value')
            return result
    lastline = f"== err : {err:.3E} = rco : {rco:.3E} = res : {res:.3E} ="
    result = result + lastline
    return result

def is_real(sol, tol, vrblvl=0):
    r"""
    Returns True if the solution in *sol* is real with respect
    to the given tolerance *tol*: if the absolute value of the imaginary
    part of all coordinates are less than *tol*.
    """
    if vrblvl > 0:
        print('in is_real, tol :', tol)
        print('sol :')
        print(sol)
    (_, vals) = coordinates(sol, vrblvl)
    for value in vals:
        if abs(value.imag) > tol:
            return False
    return True

def filter_real(sols, tol, oper, vrblvl=0):
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
    if vrblvl > 0:
        print('in filter_real, tol :', tol, ', oper :', oper)
        print('sols :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    result = []
    for sol in sols:
        isreal = is_real(sol, tol, vrblvl)
        if oper == 'select':
            if isreal:
                result.append(sol)
        if oper == 'remove':
            if not isreal:
                result.append(sol)
    return result

def filter_regular(sols, tol, oper, vrblvl=0):
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
    if vrblvl > 0:
        print('in filter_regular, tol :', tol, ', oper :', oper)
        print('sols :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    result = []
    for sol in sols:
        rco = diagnostics(sol, vrblvl)[1]
        if oper == 'select':
            if rco > tol:
                result.append(sol)
        if oper == 'remove':
            if rco <= tol:
                result.append(sol)
    return result

def filter_zero_coordinates(sols, varname, tol, oper, vrblvl=0):
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
    if vrblvl > 0:
        print('in filter_zero_coordinates, tol :', tol)
        print(', varname :', varname, end='')
        print(', oper :', oper)
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    result = []
    for sol in sols:
        dsol = strsol2dict(sol, vrblvl=vrblvl)
        if oper == 'select':
            if abs(dsol[varname]) < tol:
                result.append(sol)
        if oper == 'remove':
            if abs(dsol[varname]) >= tol:
                result.append(sol)
    return result

def is_vanishing(sol, tol, vrblvl=0):
    r"""
    Given in *sol* is a solution string and
    *tol* is the tolerance on the residual.
    Returns True if the residual of *sol*
    is less than or equal to *tol*.
    Returns False otherwise.
    """
    if vrblvl > 0:
        print('in is_vanishing, tol :', tol)
        print('sol :')
        print(sol)
    dgn = diagnostics(sol, vrblvl)
    return dgn[2] <= tol

def filter_vanishing(sols, tol, vrblvl=0):
    r"""
    Returns the list of solutions in *sols* that have a residual
    less than or equal to the given tolerance in *tol*.
    """
    if vrblvl > 0:
        print('in filter_vanishing, tol :', tol)
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    result = []
    for sol in sols:
        if is_vanishing(sol, tol, vrblvl):
            result.append(sol)
    return result

def append_double_solution_string(nvr, sol, vrblvl=0):
    """
    Appends the string in sol to the solutions in double precision,
    where the number of variables equals nvr.
    """
    if vrblvl > 0:
        print('-> append_double_solution_string, nvr =', nvr)
        print('The solution string :')
        print(sol)
    phc = get_phcfun(vrblvl-1)
    apars = int4a2nbr([nvr, len(sol)], vrblvl=vrblvl-1)
    bsol = str2int4a(sol, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> append_double_solution_string calls phc', end='')
    retval = phc(208, apars, bsol, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def append_double_double_solution_string(nvr, sol, vrblvl=0):
    """
    Appends the string in sol to the solutions in double double precision,
    where the number of variables equals nvr.
    """
    if vrblvl > 0:
        print('-> append_double_double_solution_string, nvr =', nvr)
        print('The solution string :')
        print(sol)
    phc = get_phcfun(vrblvl-1)
    apars = int4a2nbr([nvr, len(sol)], vrblvl=vrblvl-1)
    bsol = str2int4a(sol, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> append_double_double_solution_string calls phc', end='')
    retval = phc(378, apars, bsol, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def append_quad_double_solution_string(nvr, sol, vrblvl=0):
    """
    Appends the string in sol to the solutions in quad double precision,
    where the number of variables equals nvr.
    """
    if vrblvl > 0:
        print('-> append_quad_double_solution_string, nvr =', nvr)
        print('The solution string :')
        print(sol)
    phc = get_phcfun(vrblvl-1)
    apars = int4a2nbr([nvr, len(sol)], vrblvl=vrblvl-1)
    bsol = str2int4a(sol, vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> append_quad_double_solution_string calls phc', end='')
    retval = phc(428, apars, bsol, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_double_solutions(vrblvl=0):
    """
    Clears the solutions defined in double precision.
    """
    if vrblvl > 0:
        print('in clear_double_solutions ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_double_solutions calls phc', end='')
    retval = phc(37, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_double_double_solutions(vrblvl=0):
    """
    Clears the solutions defined in double double precision.
    """
    if vrblvl > 0:
        print('in clear_double_double_solutions ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_double_double_solutions calls phc', end='')
    retval = phc(347, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_quad_double_solutions(vrblvl=0):
    """
    Clears the solutions defined in quad double precision.
    """
    if vrblvl > 0:
        print('in clear_quad_double_solutions ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> clear_quad_double_solutions calls phc', end='')
    retval = phc(397, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_solutions(nvr, sols, vrblvl=0):
    """
    Sets the solutions in double precision, with the strings in sols,
    where the number of variables equals nvr.
    """
    if vrblvl > 0:
        print('-> set_double_solutions, nvr :', nvr)
    clear_double_solutions(vrblvl)
    fail = 0
    for (ind, sol) in enumerate(sols):
        fail = append_double_solution_string(nvr, sol, vrblvl)
        if fail != 0:
            print('Solution at position', ind, 'is not appended.')
    return fail

def set_double_double_solutions(nvr, sols, vrblvl=0):
    """
    Sets the solutions in double double precision, with the strings in sols,
    where the number of variables equals nvr.
    """
    if vrblvl > 0:
        print('in set_double_double_solutions, nvr :', nvr)
    clear_double_double_solutions(vrblvl)
    fail = 0
    for (ind, sol) in enumerate(sols):
        fail = append_double_double_solution_string(nvr, sol, vrblvl)
        if fail != 0:
            print('Solution at position', ind, 'is not appended.')
    return fail

def set_quad_double_solutions(nvr, sols, vrblvl=0):
    """
    Sets the solutions in quad double precision, with the strings in sols,
    where the number of variables equals nvr.
    """
    if vrblvl > 0:
        print('in set_quad_double_solutions, nvr :', nvr)
    clear_quad_double_solutions(vrblvl)
    fail = 0
    for (ind, sol) in enumerate(sols):
        fail = append_quad_double_solution_string(nvr, sol, vrblvl)
        if fail != 0:
            print('Solution at position', ind, 'is not appended.')
    return fail

def number_double_solutions(vrblvl=0):
    """
    Returns the number of solutions in double precision.
    The vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print('in number_double_solutions ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> number_double_solutions calls phc', end='')
    retval = phc(32, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return bbb[0]

def number_double_double_solutions(vrblvl=0):
    """
    Returns the number of solutions in double double precision.
    The vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print('in number_double_double_solutions ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> number_double_double_solutions calls phc', end='')
    retval = phc(342, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return bbb[0]

def number_quad_double_solutions(vrblvl=0):
    """
    Returns the number of solutions in quad double precision.
    The vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print('in number_quad_double_solutions ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> number_quad_double_solutions calls phc', end='')
    retval = phc(392, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return bbb[0]

def get_next_double_solution(idx, vrblvl=0):
    """
    Returns the string representation of the next solution
    in double precision, at the index idx.
    The vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print('in get_next_double_solution, idx :', idx)
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(idx)) # at the given index
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_next_double_solution calls phc', end='')
    retval = phc(525, aaa, bbb, ccc, vrb)
    size = bbb[0]
    if vrblvl > 0:
        print(', return value :', retval)
        print('-> get_next_double_solution, size :', size)
    soldata = create_string_buffer(b"", 4*size)
    if vrblvl > 0:
        print('-> get_next_double_solution calls phc', end='')
    retval = phc(533, bbb, soldata, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    result = int4a2str(soldata, vrblvl=vrblvl-1)
    return result

def get_next_double_double_solution(idx, vrblvl=0):
    """
    Returns the string representation of the next solution
    in double double precision, at the index idx.
    The vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print('in get_next_double_double_solution, idx :', idx)
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(idx)) # at the given index
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_next_double_double_solution calls phc', end='')
    retval = phc(526, aaa, bbb, ccc, vrb)
    size = bbb[0]
    if vrblvl > 0:
        print(', return value :', retval)
        print('-> get_next_double_double_solution, size :', size)
    soldata = create_string_buffer(b"", 4*size)
    if vrblvl > 0:
        print('-> get_next_double_double_solution calls phc', end='')
    retval = phc(534, bbb, soldata, ccc, vrb)
    if vrblvl > 0:
        print(", return value :", retval)
    result = int4a2str(soldata, vrblvl=vrblvl-1)
    return result

def get_next_quad_double_solution(idx, vrblvl=0):
    """
    Returns the string representation of the next solution
    in quad double precision, at the index idx.
    The vrblvl is the verbose level.
    """
    if vrblvl > 0:
        print('in get_next_quad_double_solution, idx :', idx)
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(idx)) # at the given index
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> get_next_quad_double_solution calls phc', end='')
    retval = phc(527, aaa, bbb, ccc, vrb)
    size = bbb[0]
    if vrblvl > 0:
        print(', return value :', retval)
        print('-> get_next_quad_double_solution, size :', size)
    soldata = create_string_buffer(b"", 4*size)
    if vrblvl > 0:
        print('-> get_next_quad_double_solution calls phc', end='')
    retval = phc(535, bbb, soldata, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    result = int4a2str(soldata, vrblvl=vrblvl-1)
    return result

def move_double_solution_cursor(idx, vrblvl=0):
    """
    Moves the cursor to the next solution, following the index,
    in double precision.
    """
    if vrblvl > 0:
        print('in move_double_solution_cursor, idx :', idx)
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(idx)) # at the given index
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> move_double_solution_cursor calls phc', end='')
    retval = phc(454, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return aaa[0]

def move_double_double_solution_cursor(idx, vrblvl=0):
    """
    Moves the cursor to the next solution, following the index,
    in double double precision.
    """
    if vrblvl > 0:
        print('in move_double_double_solution_cursor, idx :', idx)
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(idx)) # at the given index
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> move_double_double_solution_cursor calls phc', end='')
    retval = phc(455, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return aaa[0]

def move_quad_double_solution_cursor(idx, vrblvl=0):
    """
    Moves the cursor to the next solution, following the index,
    in quad double precision.
    """
    if vrblvl > 0:
        print('in move_quad_double_solution_cursor, idx :', idx)
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(idx)) # at the given index
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> move_quad_double_solution_cursor calls phc', end='')
    retval = phc(456, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return aaa[0]

def get_double_solutions(vrblvl=0):
    """
    Returns the solution strings in double precision.
    """
    if vrblvl > 0:
        print('in get_double_solutions ...')
    nbrsols = number_double_solutions(vrblvl)
    if vrblvl > 0:
        print('number of solutions retrieved :', nbrsols)
    result = []
    if nbrsols > 0:
        sol = get_next_double_solution(1, vrblvl)
        if vrblvl > 0:
            print('the first solution :\n', sol)
        result.append(sol)
        idx = 1
        for _ in range(1, nbrsols):
            idx = move_double_solution_cursor(idx, vrblvl)
            if vrblvl > 0:
                print('the next index :', idx)
            sol = get_next_double_solution(idx, vrblvl)
            if vrblvl > 0:
                print('the solution at index', idx, ':\n', sol)
            result.append(sol)
    return result

def get_double_double_solutions(vrblvl=0):
    """
    Returns the solution strings in double double precision.
    """
    if vrblvl > 0:
        print('in get_double_double_solutions ...')
    nbrsols = number_double_double_solutions(vrblvl)
    if vrblvl > 0:
        print('number of solutions retrieved :', nbrsols)
    result = []
    if nbrsols > 0:
        sol = get_next_double_double_solution(1, vrblvl)
        if vrblvl > 0:
            print('the first solution :\n', sol)
        result.append(sol)
        idx = 1
        for _ in range(1, nbrsols):
            idx = move_double_double_solution_cursor(idx, vrblvl)
            if vrblvl > 0:
                print('the next index :', idx)
            sol = get_next_double_double_solution(idx, vrblvl)
            if vrblvl > 0:
                print('the solution at index', idx, ':\n', sol)
            result.append(sol)
    return result

def get_quad_double_solutions(vrblvl=0):
    """
    Returns the solution strings in quad double precision.
    """
    if vrblvl > 0:
        print('in get_quad_double_solutions ...')
    nbrsols = number_quad_double_solutions(vrblvl)
    if vrblvl > 0:
        print('number of solutions retrieved :', nbrsols)
    result = []
    if nbrsols > 0:
        sol = get_next_quad_double_solution(1, vrblvl)
        if vrblvl > 0:
            print('the first solution :\n', sol)
        result.append(sol)
        idx = 1
        for _ in range(1, nbrsols):
            idx = move_quad_double_solution_cursor(idx, vrblvl)
            if vrblvl > 0:
                print('the next index :', idx)
            sol = get_next_quad_double_solution(idx, vrblvl)
            if vrblvl > 0:
                print('the solution at index', idx, ':\n', sol)
            result.append(sol)
    return result

def write_double_solutions(vrblvl=0):
    """
    Writes the solutions stored in double precision.
    """
    if vrblvl > 0:
        print('in write_double_solutions, vrblvl :', vrblvl)
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> write_double_solutions calls phc', end='')
    retval = phc(31, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def write_double_double_solutions(vrblvl=0):
    """
    Writes the solutions stored in double double precision.
    """
    if vrblvl > 0:
        print('in write_double_double_solutions, vrblvl :', vrblvl)
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> write_double_double_solutions calls phc', end='')
    retval = phc(341, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def write_quad_double_solutions(vrblvl=0):
    """
    Writes the solutions stored in quad double precision.
    """
    if vrblvl > 0:
        print('in write_quad_double_solutions, vrblvl :', vrblvl)
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> write_quad_double_solutions calls phc', end='')
    retval = phc(391, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

class DoubleSolution():
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

def test_double_functions(vrblvl=0):
    """
    Generates a random trinomial system,
    solves it, converts the solutions,
    and then sums the multiplicities.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in test_double_functions ...')
    pols = ['(x - 1)*(y - 1);', '(x + 1)*(y + 1);']
    names = ['x', 'y']
    sol1 = make_solution(names, [1, -1], vrblvl=vrblvl)
    sol2 = make_solution(names, [-1, 1], vrblvl=vrblvl)
    sols = [sol1, sol2]
    dsols = [strsol2dict(sol, vrblvl=vrblvl) for sol in sols]
    mult = 0
    s0d = strsol2dict(sols[0], vrblvl=vrblvl)
    if vrblvl > 0:
        print('variables :', variables(s0d))
        print(evaluate(pols, s0d, vrblvl))
    for sol in dsols:
        mult = mult + sol['m']
    if vrblvl > 0:
        print('sum of multiplicities :', mult)
        print('sum of values at the solutions :')
    errsum = 0
    for sol in dsols:
        if vrblvl > 0:
            print(sum(evaluate(pols, sol, vrblvl)))
        errsum = errsum + abs(sum(evaluate(pols, sol, vrblvl)))
    if vrblvl > 0:
        print('error sum :', errsum)
    return errsum > 1.0e-14

def test_double_solution_class(vrblvl=0):
    """
    Tests the methods in the class DoubleSolution.
    The verbose level is given by vrblvl.
    """
    mysol = DoubleSolution({'x': complex(1,2), 'y': complex(-7,0)})
    if vrblvl > 0:
        print('my solution :')
        print(mysol)
    names = ['x', 'y']
    sol1 = make_solution(names, [1, 1], err=0.01, rco= 1.2, res=0.01)
    sol2 = make_solution(names, [1, -1], err=0.1, rco= 1.1, res=0.03)
    sols = [sol1, sol2]
    sol = DoubleSolution(sols[0])
    if vrblvl > 0:
        print('the first solution :')
        print(sol)
        print('its coordinates :\n', sol.coordinates())
        print('its variable names :\n ', sol.variables())
        print('its numerical values :\n', sol.numerals())
        print('its diagnostics :\n', sol.diagnostics())
        print('its continuation parameter :', sol.timevalue())
        print('its multiplicity :', sol.multiplicity())
        print('the dictionary :\n', sol.dictionary())
    (errtab, rcotab, restab) = condition_tables(sols)
    if vrblvl > 0:
        print('The frequency tables for err, rco, and res:')
        print(errtab)
        print(rcotab)
        print(restab)
    return 0

def main():
    """
    Runs some tests on solutions.
    """
    lvl = 1
    fail = test_double_functions(lvl)
    fail = fail + test_double_solution_class(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()
