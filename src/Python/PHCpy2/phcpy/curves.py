"""
The module curves exports functions to approximate algebraic space curves
with rational expressions, in particular Pade approximants.
"""

def get_homotopy_continuation_parameter(idx):
    """
    Returns the current value of the homotopy continuation parameter
    with index idx, where idx is an integer in range(1, 13).
    """
    from phcpy.phcpy2c2 import py2c_padcon_get_homotopy_continuation_parameter
    return py2c_padcon_get_homotopy_continuation_parameter(idx)

def set_homotopy_continuation_parameter(idx, val):
    """
    Sets the value of the homotopy continuation parameter with index idx,
    where idx is in an integer in range(2, 13), to the value val.
    """
    from phcpy.phcpy2c2 import py2c_padcon_set_homotopy_continuation_parameter
    return py2c_padcon_set_homotopy_continuation_parameter(idx, val)

def set_homotopy_continuation_gamma(regamma=0, imgamma=0):
    """
    Sets the value of the homotopy continuation gamma constant 
    to the complex number with real part in regamma
    and the imaginary part in imgamma.
    If both regamma and imgamma are zero,
    then the user is prompted to provide values for regamma and imgamma.
    """
    from phcpy.phcpy2c2 import py2c_padcon_set_homotopy_continuation_gamma
    if((regamma == 0) and (imgamma == 0)):
        regm = float(raw_input('-> give the real part of gamma : '))
        imgm = float(raw_input('-> give the imaginary part of gamma : '))
        return py2c_padcon_set_homotopy_continuation_gamma(regm, imgm)
    else:
        return py2c_padcon_set_homotopy_continuation_gamma(regamma, imgamma)

def write_homotopy_continuation_parameters():
    """
    Writes the values of the homotopy continuation parameters.
    """
    pars = [get_homotopy_continuation_parameter(k) for k in range(1, 13)]
    regamma, imgamma = pars[0]
    print "Values of the HOMOTOPY CONTINUATION PARAMETERS :"
    print " 1. gamma :", regamma + imgamma*complex(0,1)
    print " 2. degree of numerator of Pade approximant    :", pars[1]
    print " 3. degree of denominator of Pade approximant  :", pars[2]
    print " 4. maximum step size                          :", pars[3]
    print " 5. minimum step size                          :", pars[4]
    print " 6. multiplication factor of the series step   :", pars[5]
    print " 7. multiplication factor of the pole radius   :", pars[6]
    print " 8. tolerance on the residual of the predictor :", pars[7]
    print " 9. tolerance on the residual of the corrector :", pars[8]
    print "10. tolerance on zero series coefficients      :", pars[9]
    print "11. maximum number of corrector steps          :", pars[10]
    print "12. maximum steps on a path                    :", pars[11]

def tune_homotopy_continuation_parameters():
    """
    The interactive loop prompts the user to tune the parameters.
    """
    idx = 1
    while(idx > 0):
        write_homotopy_continuation_parameters()
        idxprompt = 'To change a value, give an index (0 to exit) : '
        idx = int(raw_input(idxprompt))
        if(idx > 0 and idx < 13):
            if(idx == 1):
                set_homotopy_continuation_gamma()
            else:
                valprompt = '-> give a value for parameter %d : ' % idx
                val = float(raw_input(valprompt))
                set_homotopy_continuation_parameter(idx, val);

def standard_track(target, start, sols, filename="", verbose=False):
    """
    Wraps the tracker for Pade continuation in standard double precision.
    On input are a target system, a start system with solutions,
    optionally: a string *filename* and the *verbose* flag.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in *sols*.
    The *sols* is a list of strings representing start solutions.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy.phcpy2c2 import py2c_copy_standard_container_to_target_system
    from phcpy.phcpy2c2 import py2c_copy_standard_container_to_start_system
    from phcpy.phcpy2c2 import py2c_copy_standard_container_to_start_solutions
    from phcpy.phcpy2c2 import py2c_create_standard_homotopy
    from phcpy.phcpy2c2 import py2c_create_standard_homotopy_with_gamma
    from phcpy.phcpy2c2 import py2c_solve_by_standard_homotopy_continuation
    from phcpy.phcpy2c2 import py2c_clear_standard_homotopy
    from phcpy.phcpy2c2 import py2c_clear_standard_operations_data
    from phcpy.interface import store_standard_system
    from phcpy.interface import store_standard_solutions
    from phcpy.interface import load_standard_solutions
    from phcpy.solver import number_of_symbols
    from phcpy.phcpy2c2 import py2c_padcon_standard_track
    dim = number_of_symbols(start)
    store_standard_system(target, nbvar=dim)
    py2c_copy_standard_container_to_target_system()
    store_standard_system(start, nbvar=dim)
    py2c_copy_standard_container_to_start_system()
    (regamma, imgamma) = get_homotopy_continuation_parameter(1)
    py2c_create_standard_homotopy_with_gamma(regamma, imgamma)
    store_standard_solutions(dim, sols)
    py2c_copy_standard_container_to_start_solutions()
    nbc = len(filename)
    fail = py2c_padcon_standard_track(nbc,filename,int(verbose))
    # py2c_clear_standard_homotopy()
    # py2c_clear_standard_operations_data()
    return load_standard_solutions()

def dobldobl_track(target, start, sols, filename="", verbose=False):
    """
    Wraps the tracker for Pade continuation in double double precision.
    On input are a target system, a start system with solutions,
    optionally: a string *filename* and the *verbose* flag.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in *sols*.
    The *sols* is a list of strings representing start solutions.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy.phcpy2c2 import py2c_copy_dobldobl_container_to_target_system
    from phcpy.phcpy2c2 import py2c_copy_dobldobl_container_to_start_system
    from phcpy.phcpy2c2 import py2c_copy_dobldobl_container_to_start_solutions
    from phcpy.phcpy2c2 import py2c_create_dobldobl_homotopy
    from phcpy.phcpy2c2 import py2c_create_dobldobl_homotopy_with_gamma
    from phcpy.phcpy2c2 import py2c_solve_by_dobldobl_homotopy_continuation
    from phcpy.phcpy2c2 import py2c_clear_dobldobl_homotopy
    from phcpy.phcpy2c2 import py2c_clear_dobldobl_operations_data
    from phcpy.interface import store_dobldobl_system
    from phcpy.interface import store_dobldobl_solutions
    from phcpy.interface import load_dobldobl_solutions
    from phcpy.solver import number_of_symbols
    from phcpy.phcpy2c2 import py2c_padcon_dobldobl_track
    dim = number_of_symbols(start)
    store_dobldobl_system(target, nbvar=dim)
    py2c_copy_dobldobl_container_to_target_system()
    store_dobldobl_system(start, nbvar=dim)
    py2c_copy_dobldobl_container_to_start_system()
    (regamma, imgamma) = get_homotopy_continuation_parameter(1)
    py2c_create_dobldobl_homotopy_with_gamma(regamma, imgamma)
    store_dobldobl_solutions(dim, sols)
    py2c_copy_dobldobl_container_to_start_solutions()
    nbc = len(filename)
    fail = py2c_padcon_dobldobl_track(nbc,filename,int(verbose))
    # py2c_clear_dobldobl_homotopy()
    # py2c_clear_dobldobl_operations_data()
    return load_dobldobl_solutions()

def quaddobl_track(target, start, sols, filename="", verbose=False):
    """
    Wraps the tracker for Pade continuation in quad double precision.
    On input are a target system, a start system with solutions,
    optionally: a string *filename* and the *verbose* flag.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in *sols*.
    The *sols* is a list of strings representing start solutions.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy.phcpy2c2 import py2c_copy_quaddobl_container_to_target_system
    from phcpy.phcpy2c2 import py2c_copy_quaddobl_container_to_start_system
    from phcpy.phcpy2c2 import py2c_copy_quaddobl_container_to_start_solutions
    from phcpy.phcpy2c2 import py2c_create_quaddobl_homotopy
    from phcpy.phcpy2c2 import py2c_create_quaddobl_homotopy_with_gamma
    from phcpy.phcpy2c2 import py2c_solve_by_quaddobl_homotopy_continuation
    from phcpy.phcpy2c2 import py2c_clear_quaddobl_homotopy
    from phcpy.phcpy2c2 import py2c_clear_quaddobl_operations_data
    from phcpy.interface import store_quaddobl_system
    from phcpy.interface import store_quaddobl_solutions
    from phcpy.interface import load_quaddobl_solutions
    from phcpy.solver import number_of_symbols
    from phcpy.phcpy2c2 import py2c_padcon_quaddobl_track
    dim = number_of_symbols(start)
    store_quaddobl_system(target, nbvar=dim)
    py2c_copy_quaddobl_container_to_target_system()
    store_quaddobl_system(start, nbvar=dim)
    py2c_copy_quaddobl_container_to_start_system()
    (regamma, imgamma) = get_homotopy_continuation_parameter(1)
    py2c_create_quaddobl_homotopy_with_gamma(regamma, imgamma)
    store_quaddobl_solutions(dim, sols)
    py2c_copy_quaddobl_container_to_start_solutions()
    nbc = len(filename)
    fail = py2c_padcon_quaddobl_track(nbc,filename,int(verbose))
    # py2c_clear_quaddobl_homotopy()
    # py2c_clear_quaddobl_operations_data()
    return load_quaddobl_solutions()

def test(precision='d'):
    """
    Tunes the parameters and runs a simple test on the trackers.
    The precision is either double ('d'), double double ('dd'),
    or quad double ('qd').
    """
    pars = [get_homotopy_continuation_parameter(k) for k in range(1, 13)]
    print pars
    tune_homotopy_continuation_parameters()
    from phcpy.families import katsura
    k3 = katsura(3)
    from phcpy.solver import total_degree_start_system as tdss
    (k3q, k3qsols) = tdss(k3)
    print 'tracking', len(k3qsols), 'paths ...'
    if(precision == 'd'):
        k3sols = standard_track(k3, k3q, k3qsols,"",True)
    elif(precision == 'dd'):
        k3sols = dobldobl_track(k3, k3q, k3qsols,"",True)
    elif(precision == 'qd'):
        k3sols = quaddobl_track(k3, k3q, k3qsols,"",True)
    else:
        print 'wrong precision'
    for sol in k3sols:
        print sol

if __name__ == "__main__":
    test('dd')
    # test('qd')
