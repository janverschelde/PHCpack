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

def get_curvature_beta_factor():
    """
    Returns the current multiplication factor of the curvature bound.
    This curvature bound gives an upper bound on a safe step size.
    The step size is set by multiplication of the curvature bound
    with the beta factor. 
    """
    return get_homotopy_continuation_parameter(8)

def get_predictor_residual_alpha():
    """
    Returns the current tolerance on the residual of the predictor.
    This alpha parameter controls the accuracy of the tracking.
    As long as the residual of the evaluated predicted solution
    is larger than alpha, the step size is cut in half.
    """
    return get_homotopy_continuation_parameter(9)

def get_corrector_residual_tolerance():
    """
    Returns the current tolerance on the corrector residual.
    The corrector stops if the residual of the current approximation
    drops below this tolerance.
    """
    return get_homotopy_continuation_parameter(10)

def get_zero_series_coefficient_tolerance():
    """
    Returns the current tolerance on the series coefficient to be zero.
    A coefficient in a power series will be considered as zero if
    its absolute value drops below this tolerance.
    """
    return get_homotopy_continuation_parameter(11)

def get_maximum_corrector_steps():
    """
    Returns the current value of the maximum number of corrector steps
    executed after the predictor stage.
    """
    return get_homotopy_continuation_parameter(12)

def get_maximum_steps_on_path():
    """
    Returns the current value of the maximum number of steps on a path.
    The path trackers abandons the tracking of a path once the number
    of steps reaches this maximum number.
    """
    return get_homotopy_continuation_parameter(13)

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
    print(" 6. multiplication factor for the pole radius  :", pars[5])
    print(" 7. multiplication factor for the curvature    :", pars[6])
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

def standard_track\
    (target, start, sols, filename="", verbose=False, mhom=0, partition=None):
    r"""
    Wraps the tracker for Pade continuation in standard double precision.
    On input are a target system, a start system with solutions,
    optionally: a string *filename* and the *verbose* flag.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in *sols*.
    The *sols* is a list of strings representing start solutions.
    By default *mhom* is zero and tracking happens in the original coordinates,
    if *mhom* equals one, then 1-homogeneous coordinates are used, and
    if *mhom* is two or higher, then multi-homogenization applies and
    *partition* contains the index representation of the partition of
    the set of variables.  This index representation is a list of as many
    indices as the number of variables, defining which set of the partition
    each variables belongs to.
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
    fail = py2c_padcon_standard_track\
               (nbc,filename,0,int(verbose),mhom,dim,str(partition))
    # py2c_clear_standard_operations_data()
    return load_standard_solutions()

def dobldobl_track\
    (target, start, sols, filename="", verbose=False, mhom=0, partition=None):
    r"""
    Wraps the tracker for Pade continuation in double double precision.
    On input are a target system, a start system with solutions,
    optionally: a string *filename* and the *verbose* flag.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in *sols*.
    The *sols* is a list of strings representing start solutions.
    By default *mhom* is zero and tracking happens in the original coordinates,
    if *mhom* equals one, then 1-homogeneous coordinates are used, and
    if *mhom* is two or higher, then multi-homogenization applies and
    *partition* contains the index representation of the partition of
    the set of variables.  This index representation is a list of as many
    indices as the number of variables, defining which set of the partition
    each variables belongs to.
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
    fail = py2c_padcon_dobldobl_track\
               (nbc,filename,0,int(verbose),mhom,dim,str(partition))
    # py2c_clear_dobldobl_operations_data()
    return load_dobldobl_solutions()

def quaddobl_track\
    (target, start, sols, filename="", verbose=False, mhom=0, partition=None):
    r"""
    Wraps the tracker for Pade continuation in quad double precision.
    On input are a target system, a start system with solutions,
    optionally: a string *filename* and the *verbose* flag.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in *sols*.
    The *sols* is a list of strings representing start solutions.
    By default *mhom* is zero and tracking happens in the original coordinates,
    if *mhom* equals one, then 1-homogeneous coordinates are used, and
    if *mhom* is two or higher, then multi-homogenization applies and
    *partition* contains the index representation of the partition of
    the set of variables.  This index representation is a list of as many
    indices as the number of variables, defining which set of the partition
    each variables belongs to.
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
    fail = py2c_padcon_quaddobl_track\
               (nbc,filename,0,int(verbose),mhom,dim,str(partition))
    # py2c_clear_quaddobl_operations_data()
    return load_quaddobl_solutions()

def standard_set_homotopy(target, start, verbose=False, homogeneous=False):
    """
    Initializes the homotopy with the target and start system for a
    step-by-step run of the series-Pade tracker, in double precision.
    If verbose, then extra output is written.
    If homogeneous, then path tracking happens in projective space,
    otherwise the original affine coordinates are used.
    Returns the failure code of the homotopy initializer.
    """
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_start_system
    from phcpy.interface import store_standard_system
    from phcpy.solver import number_of_symbols
    dim = number_of_symbols(start)
    store_standard_system(target, nbvar=dim)
    py2c_copy_standard_container_to_target_system()
    store_standard_system(start, nbvar=dim)
    py2c_copy_standard_container_to_start_system()
    from phcpy.phcpy2c3 import py2c_padcon_standard_initialize_homotopy
    (vrb, hmg) = (int(verbose), int(homogeneous))
    return py2c_padcon_standard_initialize_homotopy(vrb, hmg)

def dobldobl_set_homotopy(target, start, verbose=False, homogeneous=False):
    """
    Initializes the homotopy with the target and start system for a
    step-by-step run of the series-Pade tracker, in double double precision.
    If verbose, then extra output is written.
    If homogeneous, then path tracking happens in projective space,
    otherwise the original affine coordinates are used.
    Returns the failure code of the homotopy initializer.
    """
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_start_system
    from phcpy.interface import store_dobldobl_system
    from phcpy.solver import number_of_symbols
    dim = number_of_symbols(start)
    store_dobldobl_system(target, nbvar=dim)
    py2c_copy_dobldobl_container_to_target_system()
    store_dobldobl_system(start, nbvar=dim)
    py2c_copy_dobldobl_container_to_start_system()
    from phcpy.phcpy2c3 import py2c_padcon_dobldobl_initialize_homotopy
    (vrb, hmg) = (int(verbose), int(homogeneous))
    return py2c_padcon_dobldobl_initialize_homotopy(vrb, hmg)

def quaddobl_set_homotopy(target, start, verbose=False, homogeneous=False):
    """
    Initializes the homotopy with the target and start system for a
    step-by-step run of the series-Pade tracker, in quad double precision.
    If verbose, then extra output is written.
    If homogeneous, then path tracking happens in projective space,
    otherwise the original affine coordinates are used.
    Returns the failure code of the homotopy initializer.
    """
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_target_system
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_start_system
    from phcpy.interface import store_quaddobl_system
    from phcpy.solver import number_of_symbols
    dim = number_of_symbols(start)
    store_quaddobl_system(target, nbvar=dim)
    py2c_copy_quaddobl_container_to_target_system()
    store_quaddobl_system(start, nbvar=dim)
    py2c_copy_quaddobl_container_to_start_system()
    from phcpy.phcpy2c3 import py2c_padcon_quaddobl_initialize_homotopy
    (vrb, hmg) = (int(verbose), int(homogeneous))
    return py2c_padcon_quaddobl_initialize_homotopy(vrb, hmg)

def standard_set_parameter_homotopy(hom, idx, verbose=False):
    """
    Initializes the homotopy with the polynomials in hom for a
    step-by-step run of the series-Pade tracker, in double precision.
    The value idx gives the index of the continuation parameter
    and is the index of one of the variables in the homotopy hom.
    If verbose, then extra output is written.
    Returns the failure code of the homotopy initializer.
    """
    from phcpy.phcpy2c3 import py2c_copy_standard_container_to_target_system
    from phcpy.interface import store_standard_system
    from phcpy.solver import number_of_symbols
    dim = number_of_symbols(hom)
    store_standard_system(hom, nbvar=dim)
    py2c_copy_standard_container_to_target_system()
    from phcpy.phcpy2c3  \
        import py2c_padcon_standard_initialize_parameter_homotopy
    vrb = int(verbose)
    return py2c_padcon_standard_initialize_parameter_homotopy(idx, vrb)

def dobldobl_set_parameter_homotopy(hom, idx, verbose=False):
    """
    Initializes the homotopy with the polynomials in hom for a
    step-by-step run of the series-Pade tracker, in double double precision.
    The value idx gives the index of the continuation parameter
    and is the index of one of the variables in the homotopy hom.
    If verbose, then extra output is written.
    Returns the failure code of the homotopy initializer.
    """
    from phcpy.phcpy2c3 import py2c_copy_dobldobl_container_to_target_system
    from phcpy.interface import store_dobldobl_system
    from phcpy.solver import number_of_symbols
    dim = number_of_symbols(hom)
    store_dobldobl_system(hom, nbvar=dim)
    py2c_copy_dobldobl_container_to_target_system()
    from phcpy.phcpy2c3  \
        import py2c_padcon_dobldobl_initialize_parameter_homotopy
    vrb = int(verbose)
    return py2c_padcon_dobldobl_initialize_parameter_homotopy(idx, vrb)

def quaddobl_set_parameter_homotopy(hom, idx, verbose=False):
    """
    Initializes the homotopy with the polynomials in hom for a
    step-by-step run of the series-Pade tracker, in quad double precision.
    The value idx gives the index of the continuation parameter
    and is the index of one of the variables in the homotopy hom.
    If verbose, then extra output is written.
    Returns the failure code of the homotopy initializer.
    """
    from phcpy.phcpy2c3 import py2c_copy_quaddobl_container_to_target_system
    from phcpy.interface import store_quaddobl_system
    from phcpy.solver import number_of_symbols
    dim = number_of_symbols(hom)
    store_quaddobl_system(hom, nbvar=dim)
    py2c_copy_quaddobl_container_to_target_system()
    from phcpy.phcpy2c3  \
        import py2c_padcon_quaddobl_initialize_parameter_homotopy
    vrb = int(verbose)
    return py2c_padcon_quaddobl_initialize_parameter_homotopy(idx, vrb)

def standard_set_solution(nvar, sol, verbose=False):
    r"""
    Sets the start solution in *sol* for the step-by-step run of
    the series-Pade tracker, in standard double precision.
    If verbose, then extra output is written.
    The number of variables is in *nvar*.
    """
    from phcpy.phcpy2c3 import py2c_padcon_initialize_standard_solution
    from phcpy.interface import store_standard_solutions
    store_standard_solutions(nvar, [sol])
    return py2c_padcon_initialize_standard_solution(1, int(verbose))

def dobldobl_set_solution(nvar, sol, verbose=False):
    r"""
    Sets the start solution in *sol* for the step-by-step run of
    the series-Pade tracker, in double double precision.
    If verbose, then extra output is written.
    The number of variables is in *nvar*.
    """
    from phcpy.phcpy2c3 import py2c_padcon_initialize_dobldobl_solution
    from phcpy.interface import store_dobldobl_solutions
    store_dobldobl_solutions(nvar, [sol])
    return py2c_padcon_initialize_dobldobl_solution(1, int(verbose))

def quaddobl_set_solution(nvar, sol, verbose=False):
    r"""
    Sets the start solution in *sol* for the step-by-step run of
    the series-Pade tracker, in quad double precision.
    If verbose, then extra output is written.
    The number of variables is in *nvar*.
    """
    from phcpy.phcpy2c3 import py2c_padcon_initialize_quaddobl_solution
    from phcpy.interface import store_quaddobl_solutions
    store_quaddobl_solutions(nvar, [sol])
    return py2c_padcon_initialize_quaddobl_solution(1, int(verbose))

def standard_get_solution(verbose=False):
    """
    Returns the current solution on the path, in double precision,
    which starts at the solution set with standard_set_solution().
    If verbose, then extra output is written.
    """
    from phcpy.phcpy2c3 import py2c_padcon_get_standard_solution
    from phcpy.phcpy2c3 import py2c_solcon_length_standard_solution_string
    from phcpy.phcpy2c3 import py2c_solcon_write_standard_solution_string
    py2c_padcon_get_standard_solution(1, int(verbose))
    lns = py2c_solcon_length_standard_solution_string(1)
    return py2c_solcon_write_standard_solution_string(1, lns)

def dobldobl_get_solution(verbose=False):
    """
    Returns the current solution on the path, in double double precision,
    which starts at the solution set with dobldobl_set_solution().
    If verbose, then extra output is written.
    """
    from phcpy.phcpy2c3 import py2c_padcon_get_dobldobl_solution
    from phcpy.phcpy2c3 import py2c_solcon_length_dobldobl_solution_string
    from phcpy.phcpy2c3 import py2c_solcon_write_dobldobl_solution_string
    py2c_padcon_get_dobldobl_solution(1, int(verbose))
    lns = py2c_solcon_length_dobldobl_solution_string(1)
    return py2c_solcon_write_dobldobl_solution_string(1, lns)

def quaddobl_get_solution(verbose=False):
    """
    Returns the current solution on the path, in quad double precision,
    which starts at the solution set with quaddobl_set_solution().
    If verbose, then extra output is written.
    """
    from phcpy.phcpy2c3 import py2c_padcon_get_quaddobl_solution
    from phcpy.phcpy2c3 import py2c_solcon_length_quaddobl_solution_string
    from phcpy.phcpy2c3 import py2c_solcon_write_quaddobl_solution_string
    py2c_padcon_get_quaddobl_solution(1, int(verbose))
    lns = py2c_solcon_length_quaddobl_solution_string(1)
    return py2c_solcon_write_quaddobl_solution_string(1, lns)

def standard_predict_correct(verbose=False):
    """
    Performs one predictor and one corrector step on the set homotopy
    and the set solution, in standard double precision.
    If verbose, then extra output is written.
    """
    from phcpy.phcpy2c3 import py2c_padcon_standard_predict_correct
    py2c_padcon_standard_predict_correct(int(verbose))

def dobldobl_predict_correct(verbose=False):
    """
    Performs one predictor and one corrector step on the set homotopy
    and the set solution, in double double precision.
    If verbose, then extra output is written.
    """
    from phcpy.phcpy2c3 import py2c_padcon_dobldobl_predict_correct
    py2c_padcon_dobldobl_predict_correct(int(verbose))

def quaddobl_predict_correct(verbose=False):
    """
    Performs one predictor and one corrector step on the set homotopy
    and the set solution, in quad double precision.
    If verbose, then extra output is written.
    """
    from phcpy.phcpy2c3 import py2c_padcon_quaddobl_predict_correct
    py2c_padcon_quaddobl_predict_correct(int(verbose))

def standard_t_value():
    """
    Returns the current t value in the tracker in double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_standard_t_value
    return py2c_padcon_standard_t_value()

def dobldobl_t_value():
    """
    Returns the current t value in the tracker in double double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_dobldobl_t_value
    return py2c_padcon_dobldobl_t_value()

def quaddobl_t_value():
    """
    Returns the current t value in the tracker in quad double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_quaddobl_t_value
    return py2c_padcon_quaddobl_t_value()

def standard_step_size():
    """
    Returns the current step size in the tracker in double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_standard_step_size
    return py2c_padcon_standard_step_size()

def dobldobl_step_size():
    """
    Returns the current step size in the tracker in double double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_dobldobl_step_size
    return py2c_padcon_dobldobl_step_size()

def quaddobl_step_size():
    """
    Returns the current step size in the tracker in quad double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_quaddobl_step_size
    return py2c_padcon_quaddobl_step_size()

def standard_series_step():
    """
    Returns the current series step in the tracker in double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_standard_series_step
    return py2c_padcon_standard_series_step()

def dobldobl_series_step():
    """
    Returns the current series step in the tracker in double double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_dobldobl_series_step
    return py2c_padcon_dobldobl_series_step()

def quaddobl_series_step():
    """
    Returns the current series step in the tracker in quad double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_quaddobl_series_step
    return py2c_padcon_quaddobl_series_step()

def standard_pole_step():
    """
    Returns the current pole step in the tracker in double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_standard_pole_step
    return py2c_padcon_standard_pole_step()

def dobldobl_pole_step():
    """
    Returns the current pole step in the tracker in double double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_dobldobl_pole_step
    return py2c_padcon_dobldobl_pole_step()

def quaddobl_pole_step():
    """
    Returns the current pole step in the tracker in quad double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_quaddobl_pole_step
    return py2c_padcon_quaddobl_pole_step()

def standard_estimated_distance():
    """
    Returns the current estimated distance to the closest solution computed
    by the tracker in double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_standard_estimated_distance
    return py2c_padcon_standard_estimated_distance()

def dobldobl_estimated_distance():
    """
    Returns the current estimated distance to the closest solution computed
    by the tracker in double double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_dobldobl_estimated_distance
    return py2c_padcon_dobldobl_estimated_distance()

def quaddobl_estimated_distance():
    """
    Returns the current estimated distance to the closest solution computed
    by the tracker in quad double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_quaddobl_estimated_distance
    return py2c_padcon_quaddobl_estimated_distance()

def standard_hessian_step():
    """
    Returns the current Hessian step in the tracker in double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_standard_hessian_step
    return py2c_padcon_standard_hessian_step()

def dobldobl_hessian_step():
    """
    Returns the current Hessian step in the tracker in double double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_dobldobl_hessian_step
    return py2c_padcon_dobldobl_hessian_step()

def quaddobl_hessian_step():
    """
    Returns the current Hessian step in the tracker in quad double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_quaddobl_hessian_step
    return py2c_padcon_quaddobl_hessian_step()

def standard_pole_radius():
    """
    Returns the smallest pole radius,
    used in the predictor in standard double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_standard_pole_radius
    return py2c_padcon_standard_pole_radius()

def dobldobl_pole_radius():
    """
    Returns the smallest pole radius,
    used in the predictor in double double precision.
    The double on return is the high part of a double double.
    """
    from phcpy.phcpy2c3 import py2c_padcon_dobldobl_pole_radius
    return py2c_padcon_dobldobl_pole_radius()

def quaddobl_pole_radius():
    """
    Returns the smallest pole radius,
    used in the predictor in quad double precision.
    The double on return is the highest part of a quad double.
    """
    from phcpy.phcpy2c3 import py2c_padcon_quaddobl_pole_radius
    return py2c_padcon_quaddobl_pole_radius()

def standard_closest_pole():
    """
    Returns a tuple with the real and imaginary part of the closest pole
    used in the predictor in standard double precision.
    The result is only meaningful is the real part is positive.
    """
    from phcpy.phcpy2c3 import py2c_padcon_standard_closest_pole
    return py2c_padcon_standard_closest_pole()

def dobldobl_closest_pole():
    """
    Returns a tuple with the real and imaginary part of the closest pole
    used in the predictor in double double precision.
    The result is only meaningful is the real part is positive.
    The double coefficients are the high parts of the double doubles.
    """
    from phcpy.phcpy2c3 import py2c_padcon_dobldobl_closest_pole
    return py2c_padcon_dobldobl_closest_pole()

def quaddobl_closest_pole():
    """
    Returns a tuple with the real and imaginary part of the closest pole
    used in the predictor in quad double precision.
    The result is only meaningful is the real part is positive.
    The double coefficients are the highest parts of the quad doubles.
    """
    from phcpy.phcpy2c3 import py2c_padcon_quaddobl_closest_pole
    return py2c_padcon_quaddobl_closest_pole()

def standard_series_coefficients(dim):
    """
    Returns a list of lists with the coefficients of the series
    computed by the predictor in standard double precision.
    On entry in dim is the number of variables.
    """
    from phcpy.phcpy2c3 import py2c_padcon_standard_series_coefficient
    deg = get_degree_of_numerator() + get_degree_of_denominator()
    result = []
    for row in range(1, dim+1):
        cfs = []
        for col in range(deg+1):
            (rcf, icf) = py2c_padcon_standard_series_coefficient(row, col, 0)
            cfs.append(complex(rcf, icf))
        result.append(cfs)
    return result

def dobldobl_series_coefficients(dim):
    """
    Returns a list of lists with the coefficients of the series
    computed by the predictor in double double precision.
    The double coefficients are the high parts of the double doubles.
    On entry in dim is the number of variables.
    """
    from phcpy.phcpy2c3 import py2c_padcon_dobldobl_series_coefficient
    deg = get_degree_of_numerator() + get_degree_of_denominator()
    result = []
    for row in range(1, dim+1):
        cfs = []
        for col in range(deg+1):
            (rcf, icf) = py2c_padcon_dobldobl_series_coefficient(row, col, 0)
            cfs.append(complex(rcf, icf))
        result.append(cfs)
    return result

def quaddobl_series_coefficients(dim):
    """
    Returns a list of lists with the coefficients of the series
    computed by the predictor in quad double precision.
    The double coefficients are the highest parts of the quad doubles.
    On entry in dim is the number of variables.
    """
    from phcpy.phcpy2c3 import py2c_padcon_quaddobl_series_coefficient
    deg = get_degree_of_numerator() + get_degree_of_denominator()
    result = []
    for row in range(1, dim+1):
        cfs = []
        for col in range(deg+1):
            (rcf, icf) = py2c_padcon_quaddobl_series_coefficient(row, col, 0)
            cfs.append(complex(rcf, icf))
        result.append(cfs)
    return result

def standard_pade_coefficients(idx):
    """
    Returns a tuple of lists with the coefficients of the Pade approximants
    computed by the predictor in standard double precision.
    The first list in the tuple holds the coefficients of the numerator,
    the second list in the tuple holds the denominator coefficients.
    On entry in idx is the index of a variable.
    """
    from phcpy.phcpy2c3 import py2c_padcon_standard_numerator_coefficient
    from phcpy.phcpy2c3 import py2c_padcon_standard_denominator_coefficient
    numdeg = get_degree_of_numerator()
    dendeg = get_degree_of_denominator()
    numcfs = []
    for col in range(numdeg+1):
        (rcf, icf) = py2c_padcon_standard_numerator_coefficient(idx, col, 0)
        numcfs.append(complex(rcf, icf))
    dencfs = []
    for col in range(dendeg+1):
        (rcf, icf) = py2c_padcon_standard_denominator_coefficient(idx, col, 0)
        dencfs.append(complex(rcf, icf))
    return (numcfs, dencfs)

def dobldobl_pade_coefficients(idx):
    """
    Returns a tuple of lists with the coefficients of the Pade approximants
    computed by the predictor in double double precision.
    The double coefficients are the high parts of the double doubles.
    The first list in the tuple holds the coefficients of the numerator,
    the second list in the tuple holds the denominator coefficients.
    On entry in idx is the index of a variable.
    """
    from phcpy.phcpy2c3 import py2c_padcon_dobldobl_numerator_coefficient
    from phcpy.phcpy2c3 import py2c_padcon_dobldobl_denominator_coefficient
    numdeg = get_degree_of_numerator()
    dendeg = get_degree_of_denominator()
    numcfs = []
    for col in range(numdeg+1):
        (rcf, icf) = py2c_padcon_dobldobl_numerator_coefficient(idx, col, 0)
        numcfs.append(complex(rcf, icf))
    dencfs = []
    for col in range(dendeg+1):
        (rcf, icf) = py2c_padcon_dobldobl_denominator_coefficient(idx, col, 0)
        dencfs.append(complex(rcf, icf))
    return (numcfs, dencfs)

def quaddobl_pade_coefficients(idx):
    """
    Returns a tuple of lists with the coefficients of the Pade approximants
    computed by the predictor in quad double precision.
    The double coefficients are the highest parts of the quad doubles.
    The first list in the tuple holds the coefficients of the numerator,
    the second list in the tuple holds the denominator coefficients.
    On entry in idx is the index of a variable.
    """
    from phcpy.phcpy2c3 import py2c_padcon_quaddobl_numerator_coefficient
    from phcpy.phcpy2c3 import py2c_padcon_quaddobl_denominator_coefficient
    numdeg = get_degree_of_numerator()
    dendeg = get_degree_of_denominator()
    numcfs = []
    for col in range(numdeg+1):
        (rcf, icf) = py2c_padcon_quaddobl_numerator_coefficient(idx, col, 0)
        numcfs.append(complex(rcf, icf))
    dencfs = []
    for col in range(dendeg+1):
        (rcf, icf) = py2c_padcon_quaddobl_denominator_coefficient(idx, col, 0)
        dencfs.append(complex(rcf, icf))
    return (numcfs, dencfs)

def standard_pade_vector(dim):
    """
    Returns the list of all coefficients over all dim variables,
    computed by the predictor in standard double precision.
    """
    result = []
    for i in range(1, dim+1):
        result.append(standard_pade_coefficients(i))
    return result

def dobldobl_pade_vector(dim):
    """
    Returns the list of all coefficients over all dim variables,
    computed by the predictor in double double precision.
    The double coefficients are the high parts of the double doubles.
    """
    result = []
    for i in range(1, dim+1):
        result.append(dobldobl_pade_coefficients(i))
    return result

def quaddobl_pade_vector(dim):
    """
    Returns the list of all coefficients over all dim variables,
    computed by the predictor in quad double precision.
    The double coefficients are the highest parts of the quad doubles.
    """
    result = []
    for i in range(1, dim+1):
        result.append(quaddobl_pade_coefficients(i))
    return result

def standard_poles(dim):
    """
    Returns a list of lists of all poles of the vector of length dim,
    computed by the predictor in standard double precision.
    """
    from phcpy.phcpy2c3 import py2c_padcon_standard_pole
    dendeg = get_degree_of_denominator()
    result = []
    for i in range(1, dim+1):
        poles = []
        for j in range(1, dendeg+1):
            (repole, impole) = py2c_padcon_standard_pole(i, j, 0)
            poles.append(complex(repole, impole))
        result.append(poles)
    return result

def dobldobl_poles(dim):
    """
    Returns a list of lists of all poles of the vector of length dim,
    computed by the predictor in double double precision.
    The doubles in the poles are the high parts of the double doubles.
    """
    from phcpy.phcpy2c3 import py2c_padcon_dobldobl_pole
    dendeg = get_degree_of_denominator()
    result = []
    for i in range(1, dim+1):
        poles = []
        for j in range(1, dendeg+1):
            (repole, impole) = py2c_padcon_dobldobl_pole(i, j, 0)
            poles.append(complex(repole, impole))
        result.append(poles)
    return result

def quaddobl_poles(dim):
    """
    Returns a list of lists of all poles of the vector of length dim,
    computed by the predictor in quad double precision.
    The doubles in the poles are the highest parts of the quad doubles.
    """
    from phcpy.phcpy2c3 import py2c_padcon_quaddobl_pole
    dendeg = get_degree_of_denominator()
    result = []
    for i in range(1, dim+1):
        poles = []
        for j in range(1, dendeg+1):
            (repole, impole) = py2c_padcon_quaddobl_pole(i, j, 0)
            poles.append(complex(repole, impole))
        result.append(poles)
    return result

def symbolic_pade_approximant(cff):
    """
    Given in cff are the coefficients of numerator and denominator
    of a Pade approximant, given as a tuple of two lists.
    On return is the string representation of the Pade approximant,
    using 't' as the variable.
    """
    (nq, dq) = cff
    num = ['+' + str(nq[k]) + ('*t**%d' % k) for k in range(len(nq))]
    numerator = ''.join(num)
    den = ['+' + str(nq[k]) + ('*t**%d' % k) for k in range(len(dq))]
    denominator = ''.join(den)
    result = '(' + numerator + ')/(' + denominator + ')'
    return result

def symbolic_pade_vector(cff):
    """
    Given in cff are the coefficients of numerator and denominator
    of a Pade vector, given as a list of tuples of two lists each.
    On return is the string representation of the Pade approximant,
    using 't' as the variable.
    """
    result = []
    for coefficients in cff:
        result.append(symbolic_pade_approximant(coefficients))
    return result

def standard_next_track(target, start, sols, verbose=False, homogeneous=False):
    r"""
    Runs the series-Pade tracker step by step in double precision,
    for an artificial-parameter homotopy.
    On input are a target system and a start system with solutions.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in *sols*.
    The *sols* is a list of strings representing start solutions.
    The function is interactive, prompting the user each time
    before performing the next predictor-corrector step.
    If *verbose*, then extra output is written.
    If *homogeneous*, then path tracking happens in projective space,
    otherwise the original affine coordinates are used.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy.solver import number_of_symbols
    result = []
    dim = number_of_symbols(start)
    standard_set_homotopy(target, start, verbose, homogeneous)
    (idx, tval) = (0, 0.0)
    fmt = 'pole step : %.3e, estimated distance : %.3e, Hessian step : %.3e'
    for sol in sols:
        idx = idx + 1
        standard_set_solution(dim, sol, verbose)
        while(True):
            answer = input('next predictor-corrector step ? (y/n) ')
            if(answer != 'y'):
                result.append(sol)
                break
            else:
                standard_predict_correct(verbose)
                sol = standard_get_solution(verbose)
                print(sol)
                polestep = standard_pole_step()
                estidist = standard_estimated_distance()
                curvstep = standard_hessian_step()
                print(fmt % (polestep, estidist, curvstep))
                previoustval = tval
                (tval, step) = (standard_t_value(), standard_step_size())
                frp = standard_pole_radius()
                print("t : %.3e, step : %.3e, frp : %.3e" % (tval, step, frp))
                cfp = standard_closest_pole()
                print('For the previous t value', previoustval, ':')
                print('1) closest pole : ', cfp)
                print('2) the series:', standard_series_coefficients(dim))
                print('3) Pade vector:', standard_pade_vector(dim))
                print('4) poles:', standard_poles(dim))
    return result

def dobldobl_next_track(target, start, sols, verbose=False, homogeneous=False):
    r"""
    Runs the series-Pade tracker step by step in double double precision,
    for an artificial-parameter homotopy.
    On input are a target system and a start system with solutions.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in *sols*.
    The *sols* is a list of strings representing start solutions.
    The function is interactive, prompting the user each time
    before performing the next predictor-corrector step.
    If *verbose*, then extra output is written.
    If *homogeneous*, then path tracking happens in projective space,
    otherwise the original affine coordinates are used.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy.solver import number_of_symbols
    result = []
    dim = number_of_symbols(start)
    dobldobl_set_homotopy(target, start, verbose, homogeneous)
    (idx, tval) = (0, 0.0)
    fmt = 'pole step : %.3e, estimated distance : %.3e, Hessian step : %.3e'
    for sol in sols:
        idx = idx + 1
        dobldobl_set_solution(dim, sol, verbose)
        while(True):
            answer = input('next predictor-corrector step ? (y/n) ')
            if(answer != 'y'):
                result.append(sol)
                break
            else:
                dobldobl_predict_correct(verbose)
                sol = dobldobl_get_solution(verbose)
                print(sol)
                polestep = dobldobl_pole_step()
                estidist = dobldobl_estimated_distance()
                curvstep = dobldobl_hessian_step()
                print(fmt % (polestep, estidist, curvstep))
                previoustval = tval
                (tval, step) = (dobldobl_t_value(), dobldobl_step_size())
                frp = dobldobl_pole_radius()
                print("t : %.3e, step : %.3e, frp : %.3e" % (tval, step, frp))
                cfp = dobldobl_closest_pole()
                print('For the previous t value', previoustval, ':')
                print('1) closest pole : ', cfp)
                print('2) the series:', dobldobl_series_coefficients(dim))
                print('3) Pade vector:', dobldobl_pade_vector(dim))
                print('4) poles:', dobldobl_poles(dim))
    return result

def quaddobl_next_track(target, start, sols, verbose=False, homogeneous=False):
    r"""
    Runs the series-Pade tracker step by step in quad double precision,
    for an artificial-parameter homotopy.
    On input are a target system and a start system with solutions.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in *sols*.
    The *sols* is a list of strings representing start solutions.
    The function is interactive, prompting the user each time
    before performing the next predictor-corrector step.
    If *verbose*, then extra output is written.
    If *homogeneous*, then path tracking happens in projective space,
    otherwise the original affine coordinates are used.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy.solver import number_of_symbols
    result = []
    dim = number_of_symbols(start)
    quaddobl_set_homotopy(target, start, verbose, homogeneous)
    (idx, tval) = (0, 0.0)
    fmt = 'pole step : %.3e, estimated distance : %.3e, Hessian step : %.3e'
    for sol in sols:
        idx = idx + 1
        quaddobl_set_solution(dim, sol, verbose)
        while(True):
            answer = input('next predictor-corrector step ? (y/n) ')
            if(answer != 'y'):
                result.append(sol)
                break
            else:
                quaddobl_predict_correct(verbose)
                sol = quaddobl_get_solution(verbose)
                print(sol)
                polestep = quaddobl_pole_step()
                estidist = quaddobl_estimated_distance()
                curvstep = quaddobl_hessian_step()
                print(fmt % (polestep, estidist, curvstep))
                previoustval = tval
                (tval, step) = (quaddobl_t_value(), quaddobl_step_size())
                frp = quaddobl_pole_radius()
                print("t : %.3e, step : %.3e, frp : %.3e" % (tval, step, frp))
                cfp = quaddobl_closest_pole()
                print('For the previous t value', previoustval, ':')
                print('1) closest pole : ', cfp)
                print('2) the series:', quaddobl_series_coefficients(dim))
                print('3) Pade vector:', quaddobl_pade_vector(dim))
                print('4) poles:', quaddobl_poles(dim))
    return result

def standard_next_loop(hom, idx, sols, verbose=False):
    """
    Runs the series-Pade tracker step by step in double precision.
    On input is a natural parameter homotopy with solutions.
    The *hom* is a list of strings representing the polynomials
    of the natural parameter homotopy.
    The *idx* is the index of the variable in hom which is
    the continuation parameter.
    The *sols* is a list of strings representing start solutions.
    The start solutions do *not* contain the value of the continuation
    parameter, which is assumed to be equal to zero.
    The function is interactive, prompting the user each time
    before performing the next predictor-corrector step.
    If verbose, then extra output is written.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy.solver import number_of_symbols
    result = []
    dim = number_of_symbols(hom) - 1
    standard_set_parameter_homotopy(hom, idx,verbose)
    (idx, tval) = (0, 0.0)
    fmt = 'pole step : %.3e, estimated distance : %.3e, Hessian step : %.3e'
    for sol in sols:
        idx = idx + 1
        print('tracking solution path', idx, '...')
        standard_set_solution(dim, sol, verbose)
        while(True):
            answer = input('next predictor-corrector step ? (y/n) ')
            if(answer != 'y'):
                result.append(sol)
                break
            else:
                standard_predict_correct(verbose)
                sol = standard_get_solution(verbose)
                print(sol)
                polestep = standard_pole_step()
                estidist = standard_estimated_distance()
                curvstep = standard_hessian_step()
                print(fmt % (polestep, estidist, curvstep))
                previoustval = tval
                (tval, step) = (standard_t_value(), standard_step_size())
                frp = standard_pole_radius()
                print('t : %.3e, step : %.3e, frp : %.3e' % (tval, step, frp))
                cfp = standard_closest_pole()
                print('For the previous t value', previoustval, ':')
                print('1) closest pole : ', cfp)
                print('2) the series:', standard_series_coefficients(dim))
                print('3) Pade vector:', standard_pade_vector(dim))
                print('4) poles:', standard_poles(dim))
    return result

def dobldobl_next_loop(hom, idx, sols, verbose=False):
    """
    Runs the series-Pade tracker step by step in double double precision.
    On input is a natural parameter homotopy with solutions.
    The *hom* is a list of strings representing the polynomials
    of the natural parameter homotopy.
    The *idx* is the index of the variable in hom which is
    the continuation parameter.
    The *sols* is a list of strings representing start solutions.
    The start solutions do *not* contain the value of the continuation
    parameter, which is assumed to be equal to zero.
    The function is interactive, prompting the user each time
    before performing the next predictor-corrector step.
    If verbose, then extra output is written.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy.solver import number_of_symbols
    result = []
    dim = number_of_symbols(hom) - 1
    dobldobl_set_parameter_homotopy(hom, idx,verbose)
    (idx, tval) = (0, 0.0)
    fmt = 'pole step : %.3e, estimated distance : %.3e, Hessian step : %.3e'
    for sol in sols:
        idx = idx + 1
        print('tracking solution path', idx, '...')
        dobldobl_set_solution(dim, sol, verbose)
        while(True):
            answer = input('next predictor-corrector step ? (y/n) ')
            if(answer != 'y'):
                result.append(sol)
                break
            else:
                dobldobl_predict_correct(verbose)
                sol = dobldobl_get_solution(verbose)
                print(sol)
                polestep = dobldobl_pole_step()
                estidist = dobldobl_estimated_distance()
                curvstep = dobldobl_hessian_step()
                print(fmt % (polestep, estidist, curvstep))
                previoustval = tval
                (tval, step) = (dobldobl_t_value(), dobldobl_step_size())
                frp = dobldobl_pole_radius()
                print('t : %.3e, step : %.3e, frp : %.3e' % (tval, step, frp))
                cfp = dobldobl_closest_pole()
                print('For the previous t value', previoustval, ':')
                print('1) closest pole : ', cfp)
                print('2) the series:', dobldobl_series_coefficients(dim))
                print('3) Pade vector:', dobldobl_pade_vector(dim))
                print('4) poles:', dobldobl_poles(dim))
    return result

def quaddobl_next_loop(hom, idx, sols, verbose=False):
    """
    Runs the series-Pade tracker step by step in quad double precision.
    On input is a natural parameter homotopy with solutions.
    The *hom* is a list of strings representing the polynomials
    of the natural parameter homotopy.
    The *idx* is the index of the variable in hom which is
    the continuation parameter.
    The *sols* is a list of strings representing start solutions.
    The start solutions do *not* contain the value of the continuation
    parameter, which is assumed to be equal to zero.
    The function is interactive, prompting the user each time
    before performing the next predictor-corrector step.
    If verbose, then extra output is written.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    from phcpy.solver import number_of_symbols
    result = []
    dim = number_of_symbols(hom) - 1
    quaddobl_set_parameter_homotopy(hom, idx,verbose)
    (idx, tval) = (0, 0.0)
    fmt = 'pole step : %.3e, estimated distance : %.3e, Hessian step : %.3e'
    for sol in sols:
        idx = idx + 1
        print('tracking solution path', idx, '...')
        quaddobl_set_solution(dim, sol, verbose)
        while(True):
            answer = input('next predictor-corrector step ? (y/n) ')
            if(answer != 'y'):
                result.append(sol)
                break
            else:
                quaddobl_predict_correct(verbose)
                sol = quaddobl_get_solution(verbose)
                print(sol)
                polestep = quaddobl_pole_step()
                estidist = quaddobl_estimated_distance()
                curvstep = quaddobl_hessian_step()
                print(fmt % (polestep, estidist, curvstep))
                previoustval = tval
                (tval, step) = (quaddobl_t_value(), quaddobl_step_size())
                frp = quaddobl_pole_radius()
                print('t : %.3e, step : %.3e, frp : %.3e' % (tval, step, frp))
                cfp = quaddobl_closest_pole()
                print('For the previous t value', previoustval, ':')
                print('1) closest pole : ', cfp)
                print('2) the series:', quaddobl_series_coefficients(dim))
                print('3) Pade vector:', quaddobl_pade_vector(dim))
                print('4) poles:', quaddobl_poles(dim))
    return result

def test_hyperbola(precision='d'):
    """
    Tests the step-by-step Pade tracker on a hyperbola,
    by default in double precision (when precision equals 'd').
    Other supported precisions are double double ('dd')
    and quad double ('qd').
    """
    from math import sqrt
    from phcpy.curves import tune_homotopy_continuation_parameters
    from phcpy.solutions import make_solution
    tune_homotopy_continuation_parameters()
    par = 0.1
    xtp = ['x^2 - (t - 0.5)^2 - 0.01;']
    # print 'initializing the homotopy :', xtp[0]
    # standard_set_parameter_homotopy(hom=xtp, idx=2, verbose=True)
    solx = sqrt(4*par**2+1)/2
    print('\nvalue of the first start solution :', solx)
    sol1 = make_solution(['x'], [solx])
    print('the first start solution :\n', sol1)
    sol2 = make_solution(['x'], [-solx])
    print('the second start solution :\n', sol2)
    # standard_set_solution(nvar=1, sol=sol1, verbose=True)
    # standard_predict_correct(verbose=True)
    ans = input('Verbose ? (y/n) ')
    vrb = (ans == 'y')
    if precision == 'd':
        print('tracking in double precision ...')
        standard_next_loop(xtp, 2, [sol1, sol2], vrb)
    elif precision == 'dd':
        print('tracking in double double precision ...')
        dobldobl_next_loop(xtp, 2, [sol1, sol2], vrb)
    elif precision == 'qd':
        print('tracking in quad double precision ...')
        quaddobl_next_loop(xtp, 2, [sol1, sol2], vrb)
    else:
        print('wrong value for the precision')

def test_simple_track(precision='d'):
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
        k3sols = standard_track(k3, k3q, k3qsols, "/tmp/out", True, 1)
    elif(precision == 'dd'):
        k3sols = dobldobl_track(k3, k3q, k3qsols, "/tmp/out", True, 1)
    elif(precision == 'qd'):
        k3sols = quaddobl_track(k3, k3q, k3qsols, "/tmp/out", True, 1)
    else:
        print('wrong precision')
    for sol in k3sols:
        print(sol)

def test_next_track(precision='d'):
    """
    Tunes the parameters and runs the step-by-step tracker.
    The precision is either double ('d'), double double ('dd'),
    or quad double ('qd').
    """
    tune_homotopy_continuation_parameters()
    from phcpy.families import katsura
    k3 = katsura(3)
    from phcpy.solver import total_degree_start_system as tdss
    (k3q, k3qsols) = tdss(k3)
    ans = input('homogeneous coordinates ? (y/n) ')
    homo = (ans == 'y')
    print('tracking', len(k3qsols), 'paths ...')
    if(precision == 'd'):
        k3sols = standard_next_track(k3, k3q, k3qsols, True, homo)
    elif(precision == 'dd'):
        k3sols = dobldobl_next_track(k3, k3q, k3qsols, True, homo)
    elif(precision == 'qd'):
        k3sols = quaddobl_next_track(k3, k3q, k3qsols, True, homo)
    else:
        print('wrong precision')
    for sol in k3sols:
        print(sol)

def test_twohom(precision='d'):
    """
    Test on a 2-homogeneous homotopy.
    """
    from phcpy.solver import linear_product_root_count as lprc
    from phcpy.solver import random_linear_product_system as rlps
    tws = ['x^2 + x - 1;', 'x*y + x - y + 3;']
    roco = lprc(tws)
    twsq, twsqsols = rlps(tws)
    print('the start system :')
    for pol in twsq: print(pol)
    print('the start solutions :')
    for sol in twsqsols: print(sol)
    if(precision == 'd'):
        sols = standard_track(tws, twsq, twsqsols, "/tmp/out", True, 2, [1, 2])
    elif(precision == 'dd'):
        sols = dobldobl_track(tws, twsq, twsqsols, "/tmp/out", True, 2, [1, 2])
    elif(precision == 'qd'):
        sols = quaddobl_track(tws, twsq, twsqsols, "/tmp/out", True, 2, [1, 2])
    else:
        print('wrong value for the precision')
    for sol in sols:
        print(sol)

if __name__ == "__main__":
    test_hyperbola()
    # test_next_track()
    # test_simple_track('qd')
