"""
The module curves exports functions to approximate algebraic space curves
with rational expressions, in particular Pade approximants.
"""

def get_homotopy_continuation_parameter(idx):
    """
    Returns the current value of the homotopy continuation parameter
    with index idx, where idx is an integer in range(1, 13).
    """
    from phcpy.phcpy2c3 import py2c_padcon_get_homotopy_continuation_parameter
    return py2c_padcon_get_homotopy_continuation_parameter(idx)

def get_gamma_constant():
    """
    Returns the current value of the gamma constant in the homotopy.
    A random value for gamma will guarantee the absence of singular
    solutions along a path, as unlucky choices belong to an algebraic set.
    A tuple of two floats is returned, respectively with the real
    and imaginary parts of the complex value for the gamma constant.
    """
    return get_homotopy_continuation_parameter(1)

def get_degree_of_numerator():
    """
    Returns the current value of the degree of the numerator of the
    Pade approximant, evaluated to predict the next solution on a path.
    """
    return get_homotopy_continuation_parameter(2)

def get_degree_of_denominator():
    """
    Returns the current value of the degree of the denominator of the
    Pade approximant, evaluated to predict the next solution on a path.
    """
    return get_homotopy_continuation_parameter(3)

def get_maximum_step_size():
    """
    Returns the current value of the maximum step size.
    The step size is the increment to the continuation parameter.
    """
    return get_homotopy_continuation_parameter(4)

def get_minimum_step_size():
    """
    Returns the current value of the minimum step size.
    The path tracking will stop if the step size is larger than the
    minimum step size and if the predictor residual is larger than
    the value for alpha, the tolerance on the predictor residual.
    """
    return get_homotopy_continuation_parameter(5)

def get_series_beta_factor():
    """
    Returns the current multiplication factor of the step size set
    by the power series approximation for the algebraic curve.
    """
    return get_homotopy_continuation_parameter(6)

def get_pole_radius_beta_factor():
    """
    Returns the current multiplication factor of the smallest pole radius.
    The smallest radius of the poles of the Pade approximant gives
    an upper bound on a safe step size.  The step size is set by
    multiplication of the smallest pole radius with the beta factor. 
    """
    return get_homotopy_continuation_parameter(7)

def get_predictor_residual_alpha():
    """
    Returns the current tolerance on the residual of the predictor.
    This alpha parameter controls the accuracy of the tracking.
    As long as the residual of the evaluated predicted solution
    is larger than alpha, the step size is cut in half.
    """
    return get_homotopy_continuation_parameter(8)

def get_corrector_residual_tolerance():
    """
    Returns the current tolerance on the corrector residual.
    The corrector stops if the residual of the current approximation
    drops below this tolerance.
    """
    return get_homotopy_continuation_parameter(9)

def get_zero_series_coefficient_tolerance():
    """
    Returns the current tolerance on the series coefficient to be zero.
    A coefficient in a power series will be considered as zero if
    its absolute value drops below this tolerance.
    """
    return get_homotopy_continuation_parameter(10)

def get_maximum_corrector_steps():
    """
    Returns the current value of the maximum number of corrector steps
    executed after the predictor stage.
    """
    return get_homotopy_continuation_parameter(11)

def get_maximum_steps_on_path():
    """
    Returns the current value of the maximum number of steps on a path.
    The path trackers abandons the tracking of a path once the number
    of steps reaches this maximum number.
    """
    return get_homotopy_continuation_parameter(12)

def set_homotopy_continuation_parameter(idx, val):
    """
    Sets the value of the homotopy continuation parameter with index idx,
    where idx is in an integer in range(2, 13), to the value val.
    """
    from phcpy.phcpy2c3 import py2c_padcon_set_homotopy_continuation_parameter
    return py2c_padcon_set_homotopy_continuation_parameter(idx, val)

def set_homotopy_continuation_gamma(regamma=0, imgamma=0):
    """
    Sets the value of the homotopy continuation gamma constant 
    to the complex number with real part in regamma
    and the imaginary part in imgamma.
    If both regamma and imgamma are zero,
    then the user is prompted to provide values for regamma and imgamma.
    """
    from phcpy.phcpy2c3 import py2c_padcon_set_homotopy_continuation_gamma
    if((regamma == 0) and (imgamma == 0)):
        regm = float(input('-> give the real part of gamma : '))
        imgm = float(input('-> give the imaginary part of gamma : '))
        return py2c_padcon_set_homotopy_continuation_gamma(regm, imgm)
    else:
        return py2c_padcon_set_homotopy_continuation_gamma(regamma, imgamma)

def write_homotopy_continuation_parameters():
    """
    Writes the values of the homotopy continuation parameters.
    """
    pars = [get_homotopy_continuation_parameter(k) for k in range(1, 13)]
    regamma, imgamma = pars[0]
    print("Values of the HOMOTOPY CONTINUATION PARAMETERS :")
    print(" 1. gamma :", regamma + imgamma*complex(0,1))
    print(" 2. degree of numerator of Pade approximant    :", pars[1])
    print(" 3. degree of denominator of Pade approximant  :", pars[2])
    print(" 4. maximum step size                          :", pars[3])
    print(" 5. minimum step size                          :", pars[4])
    print(" 6. multiplication factor of the series step   :", pars[5])
    print(" 7. multiplication factor of the pole radius   :", pars[6])
    print(" 8. tolerance on the residual of the predictor :", pars[7])
    print(" 9. tolerance on the residual of the corrector :", pars[8])
    print("10. tolerance on zero series coefficients      :", pars[9])
    print("11. maximum number of corrector steps          :", pars[10])
    print("12. maximum steps on a path                    :", pars[11])

def tune_homotopy_continuation_parameters():
    """
    The interactive loop prompts the user to tune the parameters.
    """
    idx = 1
    while(idx > 0):
        write_homotopy_continuation_parameters()
        idxprompt = 'To change a value, give an index (0 to exit) : '
        idx = int(input(idxprompt))
        if(idx > 0 and idx < 13):
            if(idx == 1):
                set_homotopy_continuation_gamma()
            else:
                valprompt = '-> give a value for parameter %d : ' % idx
                val = float(input(valprompt))
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
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_start_system
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_start_solutions
    from phcpy.phcpy2c3 import py2c_create_standard_homotopy_with_gamma
    from phcpy.phcpy2c3 import py2c_clear_standard_operations_data
    from phcpy.interface import store_standard_system
    from phcpy.interface import store_standard_solutions
    from phcpy.interface import load_standard_solutions
    from phcpy.solver import number_of_symbols
    from phcpy.phcpy2c3 import py2c_padcon_standard_track
    dim = number_of_symbols(start)
    store_standard_system(target, nbvar=dim)
    py2c_copy_standard_container_to_target_system()
    store_standard_system(start, nbvar=dim)
    py2c_copy_standard_container_to_start_system()
    store_standard_solutions(dim, sols)
    py2c_copy_standard_container_to_start_solutions()
    (regamma, imgamma) = get_homotopy_continuation_parameter(1)
    py2c_create_standard_homotopy_with_gamma(regamma, imgamma)
    nbc = len(filename)
    fail = py2c_padcon_standard_track(nbc,filename,int(verbose))
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
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_start_system
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_start_solutions
    from phcpy.phcpy2c3 import py2c_create_dobldobl_homotopy_with_gamma
    from phcpy.phcpy2c3 import py2c_clear_dobldobl_operations_data
    from phcpy.interface import store_dobldobl_system
    from phcpy.interface import store_dobldobl_solutions
    from phcpy.interface import load_dobldobl_solutions
    from phcpy.solver import number_of_symbols
    from phcpy.phcpy2c3 import py2c_padcon_dobldobl_track
    dim = number_of_symbols(start)
    store_dobldobl_system(target, nbvar=dim)
    py2c_copy_dobldobl_container_to_target_system()
    store_dobldobl_system(start, nbvar=dim)
    py2c_copy_dobldobl_container_to_start_system()
    store_dobldobl_solutions(dim, sols)
    py2c_copy_dobldobl_container_to_start_solutions()
    (regamma, imgamma) = get_homotopy_continuation_parameter(1)
    py2c_create_dobldobl_homotopy_with_gamma(regamma, imgamma)
    nbc = len(filename)
    fail = py2c_padcon_dobldobl_track(nbc,filename,int(verbose))
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
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_start_system
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_start_solutions
    from phcpy.phcpy2c3 import py2c_create_quaddobl_homotopy_with_gamma
    from phcpy.phcpy2c3 import py2c_clear_quaddobl_operations_data
    from phcpy.interface import store_quaddobl_system
    from phcpy.interface import store_quaddobl_solutions
    from phcpy.interface import load_quaddobl_solutions
    from phcpy.solver import number_of_symbols
    from phcpy.phcpy2c3 import py2c_padcon_quaddobl_track
    dim = number_of_symbols(start)
    store_quaddobl_system(target, nbvar=dim)
    py2c_copy_quaddobl_container_to_target_system()
    store_quaddobl_system(start, nbvar=dim)
    py2c_copy_quaddobl_container_to_start_system()
    store_quaddobl_solutions(dim, sols)
    py2c_copy_quaddobl_container_to_start_solutions()
    (regamma, imgamma) = get_homotopy_continuation_parameter(1)
    py2c_create_quaddobl_homotopy_with_gamma(regamma, imgamma)
    nbc = len(filename)
    fail = py2c_padcon_quaddobl_track(nbc,filename,int(verbose))
    # py2c_clear_quaddobl_operations_data()
    return load_quaddobl_solutions()

def test(precision='d'):
    """
    Tunes the parameters and runs a simple test on the trackers.
    The precision is either double ('d'), double double ('dd'),
    or quad double ('qd').
    """
    pars = [get_homotopy_continuation_parameter(k) for k in range(1, 13)]
    print(pars)
    print('gamma constant :', get_gamma_constant())
    print('degree of numerator :', get_degree_of_numerator())
    print('degree of denominator :', get_degree_of_denominator())
    print('maximum step size :', get_maximum_step_size())
    print('minimum step size :', get_minimum_step_size())
    print('series beta factor :', get_series_beta_factor())
    print('pole radius beta factor :', get_pole_radius_beta_factor())
    print('predictor residual alpha :', get_predictor_residual_alpha())
    print('corrector residual tolerance :', get_corrector_residual_tolerance())
    print('coefficient tolerance :',  get_zero_series_coefficient_tolerance())
    print('maximum #corrector steps :', get_maximum_corrector_steps())
    print('maximum #steps on path :',  get_maximum_steps_on_path())
    tune_homotopy_continuation_parameters()
    from phcpy.families import katsura
    k3 = katsura(3)
    from phcpy.solver import total_degree_start_system as tdss
    (k3q, k3qsols) = tdss(k3)
    print('tracking', len(k3qsols), 'paths ...')
    if(precision == 'd'):
        k3sols = standard_track(k3, k3q, k3qsols,"",True)
    elif(precision == 'dd'):
        k3sols = dobldobl_track(k3, k3q, k3qsols,"",True)
    elif(precision == 'qd'):
        k3sols = quaddobl_track(k3, k3q, k3qsols,"",True)
    else:
        print('wrong precision')
    for sol in k3sols:
        print(sol)

if __name__ == "__main__":
    test('qd')
