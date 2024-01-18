"""
Exports functions to compute a numerical irreducible decomposition.
"""
from ctypes import c_int32, c_double, pointer
from phcpy.version import get_phcfun
from phcpy.polynomials import number_of_symbols
from phcpy.polynomials import set_double_system, get_double_system
from phcpy.polynomials import set_double_laurent_system
from phcpy.polynomials import get_double_laurent_system
from phcpy.polynomials import set_double_double_system
from phcpy.polynomials import get_double_double_system
from phcpy.polynomials import set_double_double_laurent_system
from phcpy.polynomials import get_double_double_laurent_system
from phcpy.polynomials import set_quad_double_system, get_quad_double_system
from phcpy.polynomials import set_quad_double_laurent_system
from phcpy.polynomials import get_quad_double_laurent_system
from phcpy.solutions import get_double_solutions
from phcpy.solutions import get_double_double_solutions
from phcpy.solutions import get_quad_double_solutions

def copy_double_witset(dim, vrblvl=0):
    """
    Copies the witness set at dimension dim to the polynomials
    and the solutions in double precision.
    """
    if vrblvl > 0:
        print('in copy_double_witset, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = pointer(c_int32(vrblvl-1))
    if vrblvl > 0:
        print('-> copy_double_witset calls phc', end='')
    retval = phc(851, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_witset(dim, vrblvl=0):
    """
    Copies the witness set at dimension dim to the polynomials
    and the solutions in double double precision.
    """
    if vrblvl > 0:
        print('in copy_double_double_witset, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = pointer(c_int32(vrblvl-1))
    if vrblvl > 0:
        print('-> copy_double_double_witset calls phc', end='')
    retval = phc(853, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_witset(dim, vrblvl=0):
    """
    Copies the witness set at dimension dim to the polynomials
    and the solutions in quad double precision.
    """
    if vrblvl > 0:
        print('in copy_quad_double_witset, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = pointer(c_int32(vrblvl-1))
    if vrblvl > 0:
        print('-> copy_quad_double_witset calls phc', end='')
    retval = phc(855, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_laurent_witset(dim, vrblvl=0):
    """
    Copies the witness set at dimension dim to the Laurent polynomials
    and the solutions in double precision.
    """
    if vrblvl > 0:
        print('in copy_double_laurent_witset, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = pointer(c_int32(vrblvl-1))
    if vrblvl > 0:
        print('-> copy_double_laurent_witset calls phc', end='')
    retval = phc(852, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_laurent_witset(dim, vrblvl=0):
    """
    Copies the witness set at dimension dim to the Laurent polynomials
    and the solutions in double double precision.
    """
    if vrblvl > 0:
        print('in copy_double_double_laurent_witset, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = pointer(c_int32(vrblvl-1))
    if vrblvl > 0:
        print('-> copy_double_double_laurent_witset calls phc', end='')
    retval = phc(854, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_laurent_witset(dim, vrblvl=0):
    """
    Copies the witness set at dimension dim to the Laurent polynomials
    and the solutions in quad double precision.
    """
    if vrblvl > 0:
        print('in copy_quad_double_laurent_witset, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = pointer(c_int32(vrblvl-1))
    if vrblvl > 0:
        print('-> copy_quad_double_laurent_witset calls phc', end='')
    retval = phc(856, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_solve(pols, topdim=-1, \
    filtsols=True, factor=True, tasks=0, verbose=True, vrblvl=0):
    """
    Runs the cascades of homotopies on the polynomial system in pols
    in double precision.  The default top dimension topdim
    is the number of variables in pols minus one.
    """
    if vrblvl > 0:
        print('in double_solve, topdim :', topdim, end='')
        print(', filtsols :', filtsols, ', factor :', factor)
        print('tasks :', tasks, ', verbose :', verbose)
        print('the polynomials on input :')
        for pol in pols:
            print(pol)
    nvr = number_of_symbols(pols, vrblvl-1)
    if(topdim == -1):
        topdim = nvr - 1
    set_double_system(nvr, pols,vrblvl-1) 
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 4)()
    apars[0] = c_int32(tasks)
    apars[1] = c_int32(topdim)
    apars[2] = c_int32(filtsols)
    apars[3] = c_int32(factor)
    pars = pointer(apars)
    bvrb = pointer(c_int32(verbose))
    ccc = pointer(c_double(0.0))
    vrb = pointer(c_int32(vrblvl-1))
    if vrblvl > 0:
        print('-> double_solve calls phc', end='')
    retval = phc(845, pars, bvrb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    witsols = []
    for soldim in range(topdim+1):
        copy_double_witset(soldim)
        witpol = get_double_system(vrblvl-1)
        witpts = get_double_solutions(vrblvl-1)
        witset = (witpol, witpts)
        witsols.append(witset)
    return witsols

def double_laurent_solve(pols, topdim=-1, \
    filtsols=True, factor=True, tasks=0, verbose=True, vrblvl=0):
    """
    Runs the cascades of homotopies on the Laurent system in pols
    in double precision.  The default top dimension topdim
    is the number of variables in pols minus one.
    """
    if vrblvl > 0:
        print('in double_laurent_solve, topdim :', topdim, end='')
        print(', filtsols :', filtsols, ', factor :', factor)
        print('tasks :', tasks, ', verbose :', verbose)
        print('the polynomials on input :')
        for pol in pols:
            print(pol)
    nvr = number_of_symbols(pols, vrblvl-1)
    if(topdim == -1):
        topdim = nvr - 1
    set_double_laurent_system(nvr, pols,vrblvl-1) 
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 4)()
    apars[0] = c_int32(tasks)
    apars[1] = c_int32(topdim)
    apars[2] = c_int32(filtsols)
    apars[3] = c_int32(factor)
    pars = pointer(apars)
    bvrb = pointer(c_int32(verbose))
    ccc = pointer(c_double(0.0))
    vrb = pointer(c_int32(vrblvl-1))
    if vrblvl > 0:
        print('-> double_laurent_solve calls phc', end='')
    retval = phc(846, pars, bvrb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    witsols = []
    for soldim in range(topdim+1):
        copy_double_laurent_witset(soldim)
        witpol = get_double_laurent_system(vrblvl-1)
        witpts = get_double_solutions(vrblvl-1)
        witset = (witpol, witpts)
        witsols.append(witset)
    return witsols

def double_double_solve(pols, topdim=-1, \
    filtsols=True, factor=True, tasks=0, verbose=True, vrblvl=0):
    """
    Runs the cascades of homotopies on the polynomial system in pols
    in double double precision.  The default top dimension topdim
    is the number of variables in pols minus one.
    """
    if vrblvl > 0:
        print('in double_double_solve, topdim :', topdim, end='')
        print(', filtsols :', filtsols, ', factor :', factor)
        print('tasks :', tasks, ', verbose :', verbose)
        print('the polynomials on input :')
        for pol in pols:
            print(pol)
    nvr = number_of_symbols(pols, vrblvl-1)
    if(topdim == -1):
        topdim = nvr - 1
    set_double_double_system(nvr, pols,vrblvl-1) 
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 4)()
    apars[0] = c_int32(tasks)
    apars[1] = c_int32(topdim)
    apars[2] = c_int32(filtsols)
    apars[3] = c_int32(factor)
    pars = pointer(apars)
    bvrb = pointer(c_int32(verbose))
    ccc = pointer(c_double(0.0))
    vrb = pointer(c_int32(vrblvl-1))
    if vrblvl > 0:
        print('-> double_double_solve calls phc', end='')
    retval = phc(847, pars, bvrb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    witsols = []
    for soldim in range(topdim+1):
        copy_double_double_witset(soldim)
        witpol = get_double_double_system(vrblvl-1)
        witpts = get_double_double_solutions(vrblvl-1)
        witset = (witpol, witpts)
        witsols.append(witset)
    return witsols

def double_double_laurent_solve(pols, topdim=-1, \
    filtsols=True, factor=True, tasks=0, verbose=True, vrblvl=0):
    """
    Runs the cascades of homotopies on the Laurent system in pols
    in double double precision.  The default top dimension topdim
    is the number of variables in pols minus one.
    """
    if vrblvl > 0:
        print('in double_double_laurent_solve, topdim :', topdim, end='')
        print(', filtsols :', filtsols, ', factor :', factor)
        print('tasks :', tasks, ', verbose :', verbose)
        print('the polynomials on input :')
        for pol in pols:
            print(pol)
    nvr = number_of_symbols(pols, vrblvl-1)
    if(topdim == -1):
        topdim = nvr - 1
    set_double_double_laurent_system(nvr, pols,vrblvl-1) 
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 4)()
    apars[0] = c_int32(tasks)
    apars[1] = c_int32(topdim)
    apars[2] = c_int32(filtsols)
    apars[3] = c_int32(factor)
    pars = pointer(apars)
    bvrb = pointer(c_int32(verbose))
    ccc = pointer(c_double(0.0))
    vrb = pointer(c_int32(vrblvl-1))
    if vrblvl > 0:
        print('-> double_double_laurent_solve calls phc', end='')
    retval = phc(848, pars, bvrb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    witsols = []
    for soldim in range(topdim+1):
        copy_double_double_laurent_witset(soldim)
        witpol = get_double_double_laurent_system(vrblvl-1)
        witpts = get_double_double_solutions(vrblvl-1)
        witset = (witpol, witpts)
        witsols.append(witset)
    return witsols

def quad_double_solve(pols, topdim=-1, \
    filtsols=True, factor=True, tasks=0, verbose=True, vrblvl=0):
    """
    Runs the cascades of homotopies on the polynomial system in pols
    in quad double precision.  The default top dimension topdim
    is the number of variables in pols minus one.
    """
    if vrblvl > 0:
        print('in quad_double_solve, topdim :', topdim, end='')
        print(', filtsols :', filtsols, ', factor :', factor)
        print('tasks :', tasks, ', verbose :', verbose)
        print('the polynomials on input :')
        for pol in pols:
            print(pol)
    nvr = number_of_symbols(pols, vrblvl-1)
    if(topdim == -1):
        topdim = nvr - 1
    set_quad_double_system(nvr, pols,vrblvl-1) 
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 4)()
    apars[0] = c_int32(tasks)
    apars[1] = c_int32(topdim)
    apars[2] = c_int32(filtsols)
    apars[3] = c_int32(factor)
    pars = pointer(apars)
    bvrb = pointer(c_int32(verbose))
    ccc = pointer(c_double(0.0))
    vrb = pointer(c_int32(vrblvl-1))
    if vrblvl > 0:
        print('-> quad_double_solve calls phc', end='')
    retval = phc(849, pars, bvrb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    witsols = []
    for soldim in range(topdim+1):
        copy_quad_double_witset(soldim)
        witpol = get_quad_double_system(vrblvl-1)
        witpts = get_quad_double_solutions(vrblvl-1)
        witset = (witpol, witpts)
        witsols.append(witset)
    return witsols

def quad_double_laurent_solve(pols, topdim=-1, \
    filtsols=True, factor=True, tasks=0, verbose=True, vrblvl=0):
    """
    Runs the cascades of homotopies on the Laurent system in pols
    in quad double precision.  The default top dimension topdim
    is the number of variables in pols minus one.
    """
    if vrblvl > 0:
        print('in quad_double_laurent_solve, topdim :', topdim, end='')
        print(', filtsols :', filtsols, ', factor :', factor)
        print('tasks :', tasks, ', verbose :', verbose)
        print('the polynomials on input :')
        for pol in pols:
            print(pol)
    nvr = number_of_symbols(pols, vrblvl-1)
    if(topdim == -1):
        topdim = nvr - 1
    set_quad_double_laurent_system(nvr, pols,vrblvl-1) 
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 4)()
    apars[0] = c_int32(tasks)
    apars[1] = c_int32(topdim)
    apars[2] = c_int32(filtsols)
    apars[3] = c_int32(factor)
    pars = pointer(apars)
    bvrb = pointer(c_int32(verbose))
    ccc = pointer(c_double(0.0))
    vrb = pointer(c_int32(vrblvl-1))
    if vrblvl > 0:
        print('-> quad_double_laurent_solve calls phc', end='')
    retval = phc(850, pars, bvrb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    witsols = []
    for soldim in range(topdim+1):
        copy_quad_double_laurent_witset(soldim)
        witpol = get_quad_double_laurent_system(vrblvl-1)
        witpts = get_quad_double_solutions(vrblvl-1)
        witset = (witpol, witpts)
        witsols.append(witset)
    return witsols

def test_double_solve(vrblvl=0):
    """
    Runs a test on solving in double precision.
    """
    if vrblvl > 0:
        print('in test_double_solve ...')
    pols = ['(x1-1)*(x1-2)*(x1-3)*(x1-4);', \
            '(x1-1)*(x2-1)*(x2-2)*(x2-3);', \
            '(x1-1)*(x1-2)*(x3-1)*(x3-2);', \
            '(x1-1)*(x2-1)*(x3-1)*(x4-1);']
    sols = double_solve(pols, vrblvl=vrblvl)
    fail = 0
    degs = [4, 12, 1, 1] # degrees of the components
    for (dim, witset) in enumerate(sols):
        deg = len(witset[1])
        fail = fail + int(deg != degs[dim])
        print('degree of solution set at dimension', dim, ':', deg)
    return fail

def test_double_double_solve(vrblvl=0):
    """
    Runs a test on solving in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_solve ...')
    pols = ['(x - 1)*(y-x^2);', \
            '(x - 1)*(z-x^3);', \
            '(x^2 - 1)*(y-x^2);' ]
    sols = double_double_solve(pols, tasks=2, vrblvl=vrblvl)
    fail = 0
    degs = [0, 3, 1] # degrees of the components
    for (dim, witset) in enumerate(sols):
        deg = len(witset[1])
        fail = fail + int(deg != degs[dim])
        print('degree of solution set at dimension', dim, ':', deg)
    return fail

def test_quad_double_solve(vrblvl=0):
    """
    Runs a test on solving in quad double precision.
    """
    if vrblvl > 0:
        print('in test_double_solve ...')
    pols = ['(x - 1)*(y-x^2);', \
            '(x - 1)*(z-x^3);', \
            '(x^2 - 1)*(y-x^2);' ]
    sols = quad_double_solve(pols, tasks=2, vrblvl=vrblvl)
    fail = 0
    degs = [0, 3, 1] # degrees of the components
    for (dim, witset) in enumerate(sols):
        deg = len(witset[1])
        fail = fail + int(deg != degs[dim])
        print('degree of solution set at dimension', dim, ':', deg)
    return fail

def test_double_laurent_solve(vrblvl=0):
    """
    Runs a test on solving a Laurent system in double precision.
    """
    if vrblvl > 0:
        print('in test_double_laurent_solve ...')
    pols = ['(x^(-1) - 1)*(y-x^2);', \
            '(x^(-1) - 1)*(z-x^3);', \
            '(x^(-2) - 1)*(y-x^2);' ]
    sols = double_laurent_solve(pols, vrblvl=vrblvl)
    fail = 0
    degs = [0, 3, 1] # degrees of the components
    for (dim, witset) in enumerate(sols):
        deg = len(witset[1])
        fail = fail + int(deg != degs[dim])
        print('degree of solution set at dimension', dim, ':', deg)
    return fail

def test_double_double_laurent_solve(vrblvl=0):
    """
    Runs a test on solving a Laurent system in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_laurent_solve ...')
    pols = ['(x^(-1) - 1)*(y-x^2);', \
            '(x^(-1) - 1)*(z-x^3);', \
            '(x^(-2) - 1)*(y-x^2);' ]
    sols = double_double_laurent_solve(pols, vrblvl=vrblvl)
    fail = 0
    degs = [0, 3, 1] # degrees of the components
    for (dim, witset) in enumerate(sols):
        deg = len(witset[1])
        fail = fail + int(deg != degs[dim])
        print('degree of solution set at dimension', dim, ':', deg)
    return fail

def test_quad_double_laurent_solve(vrblvl=0):
    """
    Runs a test on solving a Laurent system in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_laurent_solve ...')
    pols = ['(x^(-1) - 1)*(y-x^2);', \
            '(x^(-1) - 1)*(z-x^3);', \
            '(x^(-2) - 1)*(y-x^2);' ]
    sols = quad_double_laurent_solve(pols, vrblvl=vrblvl)
    fail = 0
    degs = [0, 3, 1] # degrees of the components
    for (dim, witset) in enumerate(sols):
        deg = len(witset[1])
        fail = fail + int(deg != degs[dim])
        print('degree of solution set at dimension', dim, ':', deg)
    return fail

def main():
    """
    Runs some tests.
    """
    lvl = 0
    fail = test_double_solve(lvl)
    fail = fail + test_double_laurent_solve(lvl)
    fail = fail + test_double_double_solve(lvl)
    fail = fail + test_double_double_laurent_solve(lvl)
    fail = fail + test_quad_double_solve(lvl)
    fail = fail + test_quad_double_laurent_solve(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()