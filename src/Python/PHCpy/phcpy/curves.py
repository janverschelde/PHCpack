"""
The module curves exports functions to approximate algebraic space curves
with rational expressions, also known as Pade approximants,
for use in a path tracker with apriori step size control.
"""
from ctypes import c_int32, c_double, pointer
from math import sqrt
from phcpy.version import get_phcfun, str2int4a
from phcpy.polynomials import number_of_symbols
from phcpy.solutions import make_solution
from phcpy.solutions import set_double_solutions, verify
from phcpy.solutions import get_next_double_solution
from phcpy.solutions import get_next_double_double_solution
from phcpy.solutions import get_next_quad_double_solution
from phcpy.solutions import get_double_solutions, clear_double_solutions
from phcpy.solutions import set_double_double_solutions
from phcpy.solutions import get_double_double_solutions
from phcpy.solutions import clear_double_double_solutions
from phcpy.solutions import set_quad_double_solutions
from phcpy.solutions import get_quad_double_solutions
from phcpy.solutions import clear_quad_double_solutions
from phcpy.starters import total_degree_start_system
from phcpy.homotopies import set_double_target_system
from phcpy.homotopies import set_double_start_system
from phcpy.homotopies import set_double_start_solutions
from phcpy.homotopies import set_double_double_target_system
from phcpy.homotopies import set_double_double_start_system
from phcpy.homotopies import set_double_double_start_solutions
from phcpy.homotopies import set_quad_double_target_system
from phcpy.homotopies import set_quad_double_start_system
from phcpy.homotopies import set_quad_double_start_solutions
from phcpy.homotopies import set_double_homotopy
from phcpy.homotopies import set_double_double_homotopy
from phcpy.homotopies import set_quad_double_homotopy

def set_default_parameters(vrblvl=0):
    """
    Sets the default values of the parameters.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_default_parameters ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_default_parameters calls phc', end='')
    retval = phc(735, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_parameter_value(idx, value, vrblvl=0):
    """
    Sets the parameter with index idx to the given value.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_parameter_value, idx :', idx, end='')
        print(', value :', value)
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(idx))
    bval = pointer(c_int32(0))
    cpar = (c_double * 2)()
    if idx == 1:
        cpar[0] = c_double(value.real)
        cpar[1] = c_double(value.imag)
    else:
        if idx in (2, 3, 11, 12):
            bval = pointer(c_int32(value))
        else:
            bval = pointer(c_int32(0))
            cpar[0] = c_double(value)
    cval = pointer(cpar)
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_parameter_value calls phc', end='')
    retval = phc(738, aidx, bval, cval, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def get_parameter_value(idx, vrblvl=0):
    """
    Returns the value of the parameter with index idx, 
    where idx is an integer in 1, 2, .., 12.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_parameter_value, idx :', idx)
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(idx))
    bval = pointer(c_int32(0))
    cpar = (c_double * 2)()
    cval = pointer(cpar)
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> get_parameter_value calls phc', end='')
    retval = phc(737, aidx, bval, cval, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    if idx == 1:
        vals = cval[:2]
        return complex(vals[0][0], vals[0][1])
    if idx in (2, 3, 11, 12):
        return bval[0]
    vals = cval[:2]
    return vals[0][0]

def set_gamma_constant(gamma, vrblvl=0):
    """
    Sets the gamma constant to the value of the complex number
    given by gamma.  The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_gamma_constant, gamma :', gamma)
    return set_parameter_value(1, gamma, vrblvl)

def get_gamma_constant(vrblvl=0):
    """
    Returns the current value of the gamma constant in the homotopy.
    A random value for gamma will guarantee the absence of singular
    solutions along a path, as unlucky choices belong to an algebraic set.
    A tuple of two floats is returned, respectively with the real
    and imaginary parts of the complex value for the gamma constant.
    """
    if vrblvl > 0:
        print('in get_gamma_constant ...')
    return get_parameter_value(1, vrblvl)

def set_degree_of_numerator(deg, vrblvl=0):
    """
    Set the value of the degree of the numerator of the Pade approximant
    to deg.  The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_degree_of_numerator, deg :', deg)
    return set_parameter_value(2, deg, vrblvl)

def get_degree_of_numerator(vrblvl=0):
    """
    Returns the current value of the degree of the numerator of the
    Pade approximant, evaluated to predict the next solution on a path.
    """
    if vrblvl > 0:
        print('in get_degree_of_numerator ...')
    return get_parameter_value(2, vrblvl)

def set_degree_of_denominator(deg, vrblvl=0):
    """
    Set the value of the degree of the denominator of the Pade approximant
    to deg.  The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_degree_of_denominator, deg :', deg)
    return set_parameter_value(3, deg, vrblvl)

def get_degree_of_denominator(vrblvl=0):
    """
    Returns the current value of the degree of the denominator of the
    Pade approximant, evaluated to predict the next solution on a path.
    """
    if vrblvl > 0:
        print('in get_degree_of_denominator ...')
    return get_parameter_value(3, vrblvl)

def set_maximum_step_size(maxstep, vrblvl=0):
    """
    Sets the maximum value of the step size to maxstep.
    """
    if vrblvl > 0:
        print('in set_maximum_step_size, maxstep :', maxstep)
    return set_parameter_value(4, maxstep, vrblvl)

def get_maximum_step_size(vrblvl=0):
    """
    Returns the current value of the maximum step size.
    The step size is the increment to the continuation parameter.
    """
    if vrblvl > 0:
        print('in get_maximum_step_size ...')
    return get_parameter_value(4, vrblvl)

def set_minimum_step_size(minstep, vrblvl=0):
    """
    Sets the minimum value of the step size to minstep.
    """
    if vrblvl > 0:
        print('in set_minimum_step_size, minstep :', minstep)
    return set_parameter_value(5, minstep, vrblvl)

def get_minimum_step_size(vrblvl=0):
    """
    Returns the current value of the minimum step size.
    The path tracking will stop if the step size is larger than the
    minimum step size and if the predictor residual is larger than
    the value for alpha, the tolerance on the predictor residual.
    """
    if vrblvl > 0:
        print('in get_minimum_step_size ...')
    return get_parameter_value(5, vrblvl)

def set_pole_radius_beta_factor(beta, vrblvl=0):
    """
    Sets the pole radius beta factor to the value of beta.
    """
    if vrblvl > 0:
        print('in set_pole_radius_beta_factor, beta :', beta)
    return set_parameter_value(6, beta, vrblvl)

def get_pole_radius_beta_factor(vrblvl=0):
    """
    Returns the current multiplication factor of the smallest pole radius.
    The smallest radius of the poles of the Pade approximant gives
    an upper bound on a safe step size.  The step size is set by
    multiplication of the smallest pole radius with the beta factor. 
    """
    if vrblvl > 0:
        print('in get_pole_radius_beta_factor ...')
    return get_parameter_value(6, vrblvl)

def set_curvature_beta_factor(beta, vrblvl=0):
    """
    Sets the curvature beta factor to the value of beta.
    """
    if vrblvl > 0:
        print('in set_curvature_beta_factor, beta :', beta)
    return set_parameter_value(7, beta, vrblvl)

def get_curvature_beta_factor(vrblvl=0):
    """
    Returns the current multiplication factor of the curvature bound.
    This curvature bound gives an upper bound on a safe step size.
    The step size is set by multiplication of the curvature bound
    with the beta factor. 
    """
    if vrblvl > 0:
        print('in get_curvature_beta_factor ...')
    return get_parameter_value(7, vrblvl)

def set_predictor_residual_alpha(alpha, vrblvl=0):
    """
    Sets the tolerance on the residual of the predictor to alpha.
    """
    if vrblvl > 0:
        print('in set_predictor_residual_alpha, alpha :', alpha)
    return set_parameter_value(8, alpha, vrblvl)

def get_predictor_residual_alpha(vrblvl=0):
    """
    Returns the current tolerance on the residual of the predictor.
    This alpha parameter controls the accuracy of the tracking.
    As long as the residual of the evaluated predicted solution
    is larger than alpha, the step size is cut in half.
    """
    if vrblvl > 0:
        print('in get_predictor_residual_alpha ...')
    return get_parameter_value(8, vrblvl)

def set_corrector_residual_tolerance(tol, vrblvl=0):
    """
    Sets the tolerance on the corrector residual to tol.
    """
    if vrblvl > 0:
        print('in set_corrector_residual_tolerance, tol :', tol)
    return set_parameter_value(9, tol, vrblvl)

def get_corrector_residual_tolerance(vrblvl=0):
    """
    Returns the current tolerance on the corrector residual.
    The corrector stops if the residual of the current approximation
    drops below this tolerance.
    """
    if vrblvl > 0:
        print('in get_corrector_residual_tolerance ...')
    return get_parameter_value(9, vrblvl)

def set_zero_series_coefficient_tolerance(tol, vrblvl=0):
    """
    Sets the tolerance on the zero series coefficient to tol.
    """
    if vrblvl > 0:
        print('in set_zero_series_coefficient_tolerance, tol :', tol)
    return set_parameter_value(10, tol, vrblvl)

def get_zero_series_coefficient_tolerance(vrblvl=0):
    """
    Returns the current tolerance on the series coefficient to be zero.
    A coefficient in a power series will be considered as zero if
    its absolute value drops below this tolerance.
    """
    if vrblvl > 0:
        print('in get_zero_series_coefficient_tolerance ...')
    return get_parameter_value(10, vrblvl)

def set_maximum_corrector_steps(maxsteps, vrblvl=0):
    """
    Sets the maximum corrector steps to maxsteps.
    """
    if vrblvl > 0:
        print('in set_maximum_corrector_steps, maxsteps :', maxsteps)
    return set_parameter_value(11, maxsteps, vrblvl)

def get_maximum_corrector_steps(vrblvl=0):
    """
    Returns the current value of the maximum number of corrector steps
    executed after the predictor stage.
    """
    if vrblvl > 0:
        print('in get_maximum_corrector_steps ...')
    return get_parameter_value(11, vrblvl)

def set_maximum_steps_on_path(maxsteps, vrblvl=0):
    """
    Sets the maximum steps on the path to maxsteps.
    """
    if vrblvl > 0:
        print('in set_maximum_steps_on_path, maxsteps :', maxsteps)
    return set_parameter_value(12, maxsteps, vrblvl)

def get_maximum_steps_on_path(vrblvl=0):
    """
    Returns the current value of the maximum number of steps on a path.
    The path trackers abandons the tracking of a path once the number
    of steps reaches this maximum number.
    """
    if vrblvl > 0:
        print('in get_maximum_steps_on_path ...')
    return get_parameter_value(12, vrblvl)

def write_parameters(vrblvl=0):
    """
    Writes the values of the homotopy continuation parameters.
    """
    pars = [get_parameter_value(k, vrblvl) for k in range(1, 13)]
    print("Values of the HOMOTOPY CONTINUATION PARAMETERS :")
    print(" 1. gamma :", pars[0])
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

def reset_parameters(precision=0, vrblvl=0):
    """
    Resets the homotopy continuation parameters for the step-by-step
    path trackers, using the value of the precision, 0 for double,
    1 for double double, or 2 for quad double.
    """
    if vrblvl > 0:
        print('in reset_parameters, with precision :', precision)
    phc = get_phcfun(vrblvl-1)
    aprc = pointer(c_int32(precision))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> reset_parameters calls phc', end='')
    retval = phc(740, aprc, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_track(target, start, startsols, filename="", \
    mhom=0, partition=None, vrblvl=0):
    r"""
    Does path tracking in double precision.
    On input are a target system, a start system with solutions.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in sols.
    The *startsols* is a list of strings representing start solutions.
    The first optional argument is the *filename*
    for writing extra diagnostics during the tracking.
    By default *mhom* is zero and all happens in the original coordinates,
    if *mhom* equals one, then 1-homogeneous coordinates are used, and
    if *mhom* is two or higher, then multi-homogenization applies and
    *partition* contains the index representation of the partition of
    the set of variables.  This index representation is a list of as many
    indices as the number of variables, defining which set of the partition
    each variables belongs to.
    On return is a tuple, with first the gamma used in the homotopy
    and then second, the string representations of the solutions
    computed at the end of the paths.
    Note: for mhom > 0 to work, the target, start system and solution
    must be provided in homogeneneous coordinates.
    """
    if vrblvl > 0:
        print('in double_track, with target :')
        for pol in target:
            print(pol)
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    set_double_target_system(target, vrblvl)
    set_double_start_system(start, vrblvl)
    set_double_start_solutions(len(target), startsols, vrblvl)
    usedgamma = get_gamma_constant(vrblvl)
    usedgamma = set_double_homotopy(usedgamma, pwt=1, vrblvl=vrblvl)
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 6)()
    apars[0] = c_int32(0)
    apars[1] = c_int32(len(filename))
    apars[2] = c_int32(vrblvl)
    apars[3] = c_int32(mhom)
    apars[4] = c_int32(0)
    if mhom < 2:
        apars[5] = c_int32(len(target))
    else:
        apars[5] = c_int32(len(partition))
    pars = pointer(apars)
    bname = str2int4a(filename, (vrblvl > 0))
    if mhom < 2:
        ccc = pointer(c_double(0.0))
    else:
        cpart = (c_double * len(partition))()
        for (idx, nbr) in enumerate(partition):
            cpart[idx] = c_double(nbr)
        ccc = pointer(cpart)
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_track calls phc', end='')
    retval = phc(739, pars, bname, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    sols = get_double_solutions(vrblvl)
    clear_double_solutions(vrblvl)
    return (usedgamma, sols)

def double_double_track(target, start, startsols, filename="", \
    mhom=0, partition=None, vrblvl=0):
    r"""
    Does path tracking in double double precision.
    On input are a target system, a start system with solutions.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in sols.
    The *startsols* is a list of strings representing start solutions.
    The first optional argument is the *filename*
    for writing extra diagnostics during the tracking.
    By default *mhom* is zero and all happens in the original coordinates,
    if *mhom* equals one, then 1-homogeneous coordinates are used, and
    if *mhom* is two or higher, then multi-homogenization applies and
    *partition* contains the index representation of the partition of
    the set of variables.  This index representation is a list of as many
    indices as the number of variables, defining which set of the partition
    each variables belongs to.
    On return is a tuple, with first the gamma used in the homotopy
    and then second, the string representations of the solutions
    computed at the end of the paths.
    Note: for mhom > 0 to work, the target, start system and solution
    must be provided in homogeneneous coordinates.
    """
    if vrblvl > 0:
        print('in double_double_track, with target :')
        for pol in target:
            print(pol)
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    set_double_double_target_system(target, vrblvl)
    set_double_double_start_system(start, vrblvl)
    set_double_double_start_solutions(len(target), startsols, vrblvl)
    usedgamma = get_gamma_constant(vrblvl)
    usedgamma = set_double_double_homotopy(usedgamma, pwt=1, vrblvl=vrblvl)
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 6)()
    apars[0] = c_int32(1)
    apars[1] = c_int32(len(filename))
    apars[2] = c_int32(vrblvl)
    apars[3] = c_int32(mhom)
    apars[4] = c_int32(0)
    if mhom < 2:
        apars[5] = c_int32(len(target))
    else:
        apars[5] = c_int32(len(partition))
    pars = pointer(apars)
    bname = str2int4a(filename, (vrblvl > 0))
    if mhom < 2:
        ccc = pointer(c_double(0.0))
    else:
        cpart = (c_double * len(partition))()
        for (idx, nbr) in enumerate(partition):
            cpart[idx] = c_double(nbr)
        ccc = pointer(cpart)
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_double_track calls phc', end='')
    retval = phc(739, pars, bname, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    sols = get_double_double_solutions(vrblvl)
    clear_double_double_solutions(vrblvl)
    return (usedgamma, sols)

def quad_double_track(target, start, startsols, filename="", \
    mhom=0, partition=None, vrblvl=0):
    r"""
    Does path tracking in quad double precision.
    On input are a target system, a start system with solutions.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in sols.
    The *startsols* is a list of strings representing start solutions.
    The first optional argument is the *filename*
    for writing extra diagnostics during the tracking.
    By default *mhom* is zero and all happens in the original coordinates,
    if *mhom* equals one, then 1-homogeneous coordinates are used, and
    if *mhom* is two or higher, then multi-homogenization applies and
    *partition* contains the index representation of the partition of
    the set of variables.  This index representation is a list of as many
    indices as the number of variables, defining which set of the partition
    each variables belongs to.
    On return is a tuple, with first the gamma used in the homotopy
    and then second, the string representations of the solutions
    computed at the end of the paths.
    Note: for mhom > 0 to work, the target, start system and solution
    must be provided in homogeneneous coordinates.
    """
    if vrblvl > 0:
        print('in quad_double_track, with target :')
        for pol in target:
            print(pol)
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    set_quad_double_target_system(target, vrblvl)
    set_quad_double_start_system(start, vrblvl)
    set_quad_double_start_solutions(len(target), startsols, vrblvl)
    usedgamma = get_gamma_constant(vrblvl)
    usedgamma = set_quad_double_homotopy(usedgamma, pwt=1, vrblvl=vrblvl)
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 6)()
    apars[0] = c_int32(2)
    apars[1] = c_int32(len(filename))
    apars[2] = c_int32(vrblvl)
    apars[3] = c_int32(mhom)
    apars[4] = c_int32(0)
    if mhom < 2:
        apars[5] = c_int32(len(target))
    else:
        apars[5] = c_int32(len(partition))
    pars = pointer(apars)
    bname = str2int4a(filename, (vrblvl > 0))
    if mhom < 2:
        ccc = pointer(c_double(0.0))
    else:
        cpart = (c_double * len(partition))()
        for (idx, nbr) in enumerate(partition):
            cpart[idx] = c_double(nbr)
        ccc = pointer(cpart)
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> quad_double_track calls phc', end='')
    retval = phc(739, pars, bname, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    sols = get_quad_double_solutions(vrblvl)
    clear_quad_double_solutions(vrblvl)
    return (usedgamma, sols)

def initialize_double_artificial_homotopy(target, start, \
    homogeneous=False, vrblvl=0):
    """
    Initializes the homotopy with the target and start system for a
    step-by-step run of the series-Pade tracker, in double precision.
    If homogeneous, then path tracking happens in projective space,
    otherwise the original affine coordinates are used.
    If vrblvl > 0, then extra output is written.
    Returns the failure code of the homotopy initializer.
    """
    if vrblvl > 0:
        print('in initialize_double_artificial_homotopy, homogeneous :', \
            homogeneous)
        print('the target system :')
        for pol in target:
            print(pol)
        print('the start system :')
        for pol in start:
            print(pol)
    set_double_target_system(target, vrblvl)
    set_double_start_system(start, vrblvl)
    phc = get_phcfun(vrblvl-1)
    aprc = pointer(c_int32(0))
    bpars = (c_int32 * 2)()
    bpars[0] = c_int32(vrblvl)
    bpars[1] = c_int32(int(homogeneous))
    bpar = pointer(bpars)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> initialize_double_artificial_homotopy calls phc', end='')
    retval = phc(860, aprc, bpar, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_double_double_artificial_homotopy(target, start, \
    homogeneous=False, vrblvl=0):
    """
    Initializes the homotopy with the target and start system for a
    step-by-step run of the series-Pade tracker, in double double precision.
    If homogeneous, then path tracking happens in projective space,
    otherwise the original affine coordinates are used.
    If vrblvl > 0, then extra output is written.
    Returns the failure code of the homotopy initializer.
    """
    if vrblvl > 0:
        print('in initialize_double_double_artificial_homotopy', end='')
        print(', homogeneous :', homogeneous)
        print('the target system :')
        for pol in target:
            print(pol)
        print('the start system :')
        for pol in start:
            print(pol)
    set_double_double_target_system(target, vrblvl)
    set_double_double_start_system(start, vrblvl)
    phc = get_phcfun(vrblvl-1)
    aprc = pointer(c_int32(1))
    bpars = (c_int32 * 2)()
    bpars[0] = c_int32(vrblvl)
    bpars[1] = c_int32(int(homogeneous))
    bpar = pointer(bpars)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> initialize_double_double_artificial_homotopy calls phc', \
            end='')
    retval = phc(860, aprc, bpar, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_quad_double_artificial_homotopy(target, start, \
    homogeneous=False, vrblvl=0):
    """
    Initializes the homotopy with the target and start system for a
    step-by-step run of the series-Pade tracker, in quad double precision.
    If homogeneous, then path tracking happens in projective space,
    otherwise the original affine coordinates are used.
    If vrblvl > 0, then extra output is written.
    Returns the failure code of the homotopy initializer.
    """
    if vrblvl > 0:
        print('in initialize_quad_double_artificial_homotopy', end='')
        print(', homogeneous :', homogeneous)
        print('the target system :')
        for pol in target:
            print(pol)
        print('the start system :')
        for pol in start:
            print(pol)
    set_quad_double_target_system(target, vrblvl)
    set_quad_double_start_system(start, vrblvl)
    phc = get_phcfun(vrblvl-1)
    aprc = pointer(c_int32(2))
    bpars = (c_int32 * 2)()
    bpars[0] = c_int32(vrblvl)
    bpars[1] = c_int32(int(homogeneous))
    bpar = pointer(bpars)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> initialize_quad_double_artificial_homotopy calls phc', \
            end='')
    retval = phc(860, aprc, bpar, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_double_parameter_homotopy(hom, idx, vrblvl=0):
    """
    Initializes the homotopy with the polynomials in hom for a
    step-by-step run of the series-Pade tracker, in double precision.
    The value idx gives the index of the continuation parameter
    and is the index of one of the variables in the homotopy hom.
    If vrblvl > 0, then extra output is written.
    Returns the failure code of the homotopy initializer.
    """
    if vrblvl > 0:
        print('in initialize_double_parameter_homotopy, idx :', idx)
        print('the natural parameter homotopy :')
        for pol in hom:
            print(pol)
    set_double_target_system(hom, vrblvl)
    phc = get_phcfun(vrblvl-1)
    aprc = pointer(c_int32(0))
    bpars = (c_int32 * 2)()
    bpars[0] = c_int32(idx)
    bpars[1] = c_int32(vrblvl)
    bpar = pointer(bpars)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> initialize_double_parameter_homotopy calls phc', end='')
    retval = phc(878, aprc, bpar, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_double_double_parameter_homotopy(hom, idx, vrblvl=0):
    """
    Initializes the homotopy with the polynomials in hom for a
    step-by-step run of the series-Pade tracker, in double double precision.
    The value idx gives the index of the continuation parameter
    and is the index of one of the variables in the homotopy hom.
    If vrblvl > 0, then extra output is written.
    Returns the failure code of the homotopy initializer.
    """
    if vrblvl > 0:
        print('in initialize_double_double_parameter_homotopy, idx :', idx)
        print('the natural parameter homotopy :')
        for pol in hom:
            print(pol)
    set_double_double_target_system(hom, vrblvl)
    phc = get_phcfun(vrblvl-1)
    aprc = pointer(c_int32(1))
    bpars = (c_int32 * 2)()
    bpars[0] = c_int32(idx)
    bpars[1] = c_int32(vrblvl)
    bpar = pointer(bpars)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> initialize_double_double_parameter_homotopy calls phc', \
            end='')
    retval = phc(878, aprc, bpar, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_quad_double_parameter_homotopy(hom, idx, vrblvl=0):
    """
    Initializes the homotopy with the polynomials in hom for a
    step-by-step run of the series-Pade tracker, in quad double precision.
    The value idx gives the index of the continuation parameter
    and is the index of one of the variables in the homotopy hom.
    If vrblvl > 0, then extra output is written.
    Returns the failure code of the homotopy initializer.
    """
    if vrblvl > 0:
        print('in initialize_quad_double_parameter_homotopy, idx :', idx)
        print('the natural parameter homotopy :')
        for pol in hom:
            print(pol)
    set_quad_double_target_system(hom, vrblvl)
    phc = get_phcfun(vrblvl-1)
    aprc = pointer(c_int32(2))
    bpars = (c_int32 * 2)()
    bpars[0] = c_int32(idx)
    bpars[1] = c_int32(vrblvl)
    bpar = pointer(bpars)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> initialize_quad_double_parameter_homotopy calls phc', end='')
    retval = phc(878, aprc, bpar, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_solution(nvr, sol, vrblvl=0):
    r"""
    Sets the start solution in *sol* for the step-by-step run of
    the series-Pade tracker, in double precision.
    The number of variables is in *nvr*.
    If vrblvl > 0, then extra output is written.
    """
    if vrblvl > 0:
        print('in initialize_double_solution, nvr :', nvr)
        print('the solution :')
        print(sol)
    set_double_solutions(nvr, [sol])
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 2)()
    apars[0] = c_int32(0)
    apars[1] = c_int32(1)
    apar = pointer(apars)
    bvrb = pointer(c_int32(vrblvl))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_solution calls phc', end='')
    retval = phc(861, apar, bvrb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_double_solution(nvr, sol, vrblvl=0):
    r"""
    Sets the start solution in *sol* for the step-by-step run of
    the series-Pade tracker, in double double precision.
    The number of variables is in *nvr*.
    If vrblvl > 0, then extra output is written.
    """
    if vrblvl > 0:
        print('in set_double_double_solution, nvr :', nvr)
        print('the solution :')
        print(sol)
    set_double_double_solutions(nvr, [sol])
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 2)()
    apars[0] = c_int32(1)
    apars[1] = c_int32(1)
    apar = pointer(apars)
    bvrb = pointer(c_int32(vrblvl))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_double_double_solution calls phc', end='')
    retval = phc(861, apar, bvrb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_quad_double_solution(nvr, sol, vrblvl=0):
    r"""
    Sets the start solution in *sol* for the step-by-step run of
    the series-Pade tracker, in quad double precision.
    The number of variables is in *nvr*.
    If vrblvl > 0, then extra output is written.
    """
    if vrblvl > 0:
        print('in set_quad_double_solution, nvr :', nvr)
        print('the solution :')
        print(sol)
    set_quad_double_solutions(nvr, [sol])
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 2)()
    apars[0] = c_int32(2)
    apars[1] = c_int32(1)
    apar = pointer(apars)
    bvrb = pointer(c_int32(vrblvl))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> set_quad_double_solution calls phc', end='')
    retval = phc(861, apar, bvrb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def get_double_solution(vrblvl=0):
    """
    Returns the current solution on the path, in double precision,
    which starts at the solution set with set_double_solution().
    If vrblvl > 0, then extra output is written.
    """
    if vrblvl > 0:
        print('in get_double_solution ...')
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 2)()
    apars[0] = c_int32(0)
    apars[1] = c_int32(1)
    apar = pointer(apars)
    bvrb = pointer(c_int32(vrblvl))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> get_double_solution calls phc', end='')
    retval = phc(863, apar, bvrb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    sol = get_next_double_solution(1, vrblvl)
    return sol

def get_double_double_solution(vrblvl=0):
    """
    Returns the current solution on the path, in double double precision,
    which starts at the solution set with set_double_double_solution().
    If vrblvl > 0, then extra output is written.
    """
    if vrblvl > 0:
        print('in get_double_double_solution ...')
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 2)()
    apars[0] = c_int32(1)
    apars[1] = c_int32(1)
    apar = pointer(apars)
    bvrb = pointer(c_int32(vrblvl))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> get_double_double_solution calls phc', end='')
    retval = phc(863, apar, bvrb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    sol = get_next_double_double_solution(1, vrblvl)
    return sol

def get_quad_double_solution(vrblvl=0):
    """
    Returns the current solution on the path, in quad double precision,
    which starts at the solution set with set_quad_double_solution().
    If vrblvl > 0, then extra output is written.
    """
    if vrblvl > 0:
        print('in get_quad_double_solution ...')
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 2)()
    apars[0] = c_int32(2)
    apars[1] = c_int32(1)
    apar = pointer(apars)
    bvrb = pointer(c_int32(vrblvl))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> get_quad_double_solution calls phc', end='')
    retval = phc(863, apar, bvrb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    sol = get_next_quad_double_solution(1, vrblvl)
    return sol

def get_double_predicted_solution(vrblvl=0):
    """
    Returns the predicted solution on the path, in double precision,
    which starts at the solution set with set_double_solution().
    If vrblvl > 0, then extra output is written.
    """
    if vrblvl > 0:
        print('in get_double_predicted_solution ...')
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 2)()
    apars[0] = c_int32(0)
    apars[1] = c_int32(1)
    apar = pointer(apars)
    bvrb = pointer(c_int32(vrblvl))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> get_double_predicted_solution calls phc', end='')
    retval = phc(919, apar, bvrb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    sol = get_next_double_solution(1, vrblvl)
    return sol

def get_double_double_predicted_solution(vrblvl=0):
    """
    Returns the predicted solution on the path, in double double precision,
    which starts at the solution set with set_double_double_solution().
    If vrblvl > 0, then extra output is written.
    """
    if vrblvl > 0:
        print('in get_double_double_predicted_solution ...')
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 2)()
    apars[0] = c_int32(1)
    apars[1] = c_int32(1)
    apar = pointer(apars)
    bvrb = pointer(c_int32(vrblvl))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> get_double_double_predicted_solution calls phc', end='')
    retval = phc(919, apar, bvrb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    sol = get_next_double_double_solution(1, vrblvl)
    return sol

def get_quad_double_predicted_solution(vrblvl=0):
    """
    Returns the predicted solution on the path, in quad double precision,
    which starts at the solution set with set_quad_double_solution().
    If vrblvl > 0, then extra output is written.
    """
    if vrblvl > 0:
        print('in get_quad_double_predicted_solution ...')
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 2)()
    apars[0] = c_int32(2)
    apars[1] = c_int32(1)
    apar = pointer(apars)
    bvrb = pointer(c_int32(vrblvl))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> get_quad_double_predicted_solution calls phc', end='')
    retval = phc(919, apar, bvrb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    sol = get_next_quad_double_solution(1, vrblvl)
    return sol

def double_predict_correct(vrblvl=0):
    """
    Performs one predictor and one corrector step on the set homotopy
    and the set solution, in double precision.
    If vrblvl > 0, then extra output is written.
    """
    if vrblvl > 0:
        print('in double_predict_correct ...')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(0))
    bvrb = pointer(c_int32(vrblvl))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_predict_correct calls phc', end='')
    retval = phc(862, apar, bvrb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_double_predict_correct(vrblvl=0):
    """
    Performs one predictor and one corrector step on the set homotopy
    and the set solution, in double double precision.
    If vrblvl > 0, then extra output is written.
    """
    if vrblvl > 0:
        print('in double_double_predict_correct ...')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(1))
    bvrb = pointer(c_int32(vrblvl))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_double_predict_correct calls phc', end='')
    retval = phc(862, apar, bvrb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def quad_double_predict_correct(vrblvl=0):
    """
    Performs one predictor and one corrector step on the set homotopy
    and the set solution, in quad double precision.
    If vrblvl > 0, then extra output is written.
    """
    if vrblvl > 0:
        print('in quad_double_predict_correct ...')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(2))
    bvrb = pointer(c_int32(vrblvl))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> quad_double_predict_correct calls phc', end='')
    retval = phc(862, apar, bvrb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_t_value(vrblvl=0):
    """
    Returns the current t value in the tracker in double precision.
    """
    if vrblvl > 0:
        print('in double_t_value ...')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(0))
    bvrb = pointer(c_int32(0))
    ctval = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_t_value calls phc', end='')
    retval = phc(867, apar, bvrb, ctval, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the t value :', ctval[0])
    return ctval[0]

def double_double_t_value(vrblvl=0):
    """
    Returns the current t value in the tracker in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_t_value ...')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(1))
    bvrb = pointer(c_int32(0))
    ctval = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_double_t_value calls phc', end='')
    retval = phc(867, apar, bvrb, ctval, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the t value :', ctval[0])
    return ctval[0]

def quad_double_t_value(vrblvl=0):
    """
    Returns the current t value in the tracker in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_t_value ...')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(2))
    bvrb = pointer(c_int32(0))
    ctval = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> quad_double_t_value calls phc', end='')
    retval = phc(867, apar, bvrb, ctval, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the t value :', ctval[0])
    return ctval[0]

def double_step_size(vrblvl=0):
    """
    Returns the step size in the tracker in double precision.
    """
    if vrblvl > 0:
        print('in double_step_size ...')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(0))
    bvrb = pointer(c_int32(0))
    cstep = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_step_size calls phc', end='')
    retval = phc(868, apar, bvrb, cstep, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the step size :', cstep[0])
    return cstep[0]

def double_double_step_size(vrblvl=0):
    """
    Returns the step size in the tracker in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_step_size ...')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(1))
    bvrb = pointer(c_int32(0))
    cstep = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_double_step_size calls phc', end='')
    retval = phc(868, apar, bvrb, cstep, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the step size :', cstep[0])
    return cstep[0]

def quad_double_step_size(vrblvl=0):
    """
    Returns the step size in the tracker in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_step_size ...')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(2))
    bvrb = pointer(c_int32(0))
    cstep = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> quad_double_step_size calls phc', end='')
    retval = phc(868, apar, bvrb, cstep, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the step size :', cstep[0])
    return cstep[0]

def double_pole_step(vrblvl=0):
    """
    Returns the pole step in the tracker in double precision.
    """
    if vrblvl > 0:
        print('in double_pole_step ...')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(0))
    bvrb = pointer(c_int32(0))
    cstep = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_pole_step calls phc', end='')
    retval = phc(886, apar, bvrb, cstep, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the step size :', cstep[0])
    return cstep[0]

def double_double_pole_step(vrblvl=0):
    """
    Returns the pole step in the tracker in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_pole_step ...')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(1))
    bvrb = pointer(c_int32(0))
    cstep = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_double_pole_step calls phc', end='')
    retval = phc(886, apar, bvrb, cstep, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the step size :', cstep[0])
    return cstep[0]

def quad_double_pole_step(vrblvl=0):
    """
    Returns the pole step in the tracker in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_pole_step ...')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(2))
    bvrb = pointer(c_int32(0))
    cstep = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> quad_double_pole_step calls phc', end='')
    retval = phc(886, apar, bvrb, cstep, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the step size :', cstep[0])
    return cstep[0]

def double_estimated_distance(vrblvl=0):
    """
    Returns the estimated distance to the closest solution computed
    by the tracker in double precision.
    """
    if vrblvl > 0:
        print('in double_estimated_distance ...')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(0))
    bvrb = pointer(c_int32(0))
    cdist = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_estimated_distance calls phc', end='')
    retval = phc(887, apar, bvrb, cdist, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the estimated distance :', cdist[0])
    return cdist[0]

def double_double_estimated_distance(vrblvl=0):
    """
    Returns the estimated distance to the closest solution computed
    by the tracker in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_estimated_distance ...')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(1))
    bvrb = pointer(c_int32(0))
    cdist = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_double_estimated_distance calls phc', end='')
    retval = phc(887, apar, bvrb, cdist, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the estimated distance :', cdist[0])
    return cdist[0]

def quad_double_estimated_distance(vrblvl=0):
    """
    Returns the estimated distance to the closest solution computed
    by the tracker in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_estimated_distance ...')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(2))
    bvrb = pointer(c_int32(0))
    cdist = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> quad_double_estimated_distance calls phc', end='')
    retval = phc(887, apar, bvrb, cdist, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the estimated distance :', cdist[0])
    return cdist[0]

def double_hessian_step(vrblvl=0):
    """
    Returns the Hessian step in the tracker in double precision.
    """
    if vrblvl > 0:
        print('in double_hessian_step ...')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(0))
    bvrb = pointer(c_int32(0))
    cstep = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_hessian_step calls phc', end='')
    retval = phc(888, apar, bvrb, cstep, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the step size :', cstep[0])
    return cstep[0]

def double_double_hessian_step(vrblvl=0):
    """
    Returns the Hessian step in the tracker in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_hessian_step ...')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(1))
    bvrb = pointer(c_int32(0))
    cstep = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_double_hessian_step calls phc', end='')
    retval = phc(888, apar, bvrb, cstep, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the step size :', cstep[0])
    return cstep[0]

def quad_double_hessian_step(vrblvl=0):
    """
    Returns the Hessian step in the tracker in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_hessian_step ...')
    phc = get_phcfun(vrblvl-1)
    apar = pointer(c_int32(2))
    bvrb = pointer(c_int32(0))
    cstep = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> quad_double_hessian_step calls phc', end='')
    retval = phc(888, apar, bvrb, cstep, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the step size :', cstep[0])
    return cstep[0]

def double_pole_radius(vrblvl=0):
    """
    Returns the smallest pole radius,
    used in the predictor in double precision.
    """
    if vrblvl > 0:
        print('in double_pole_radius ...')
    phc = get_phcfun(vrblvl-1)
    aprc = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    crad = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_pole_radius calls phc', end='')
    retval = phc(885, aprc, bbb, crad, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the pole radius :', crad[0])
    return crad[0]

def double_double_pole_radius(vrblvl=0):
    """
    Returns the smallest pole radius,
    used in the predictor in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_pole_radius ...')
    phc = get_phcfun(vrblvl-1)
    aprc = pointer(c_int32(1))
    bbb = pointer(c_int32(0))
    crad = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_double_pole_radius calls phc', end='')
    retval = phc(885, aprc, bbb, crad, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the pole radius :', crad[0])
    return crad[0]

def quad_double_pole_radius(vrblvl=0):
    """
    Returns the smallest pole radius,
    used in the predictor in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_pole_radius ...')
    phc = get_phcfun(vrblvl-1)
    aprc = pointer(c_int32(2))
    bbb = pointer(c_int32(0))
    crad = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> quad_double_pole_radius calls phc', end='')
    retval = phc(885, aprc, bbb, crad, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('the pole radius :', crad[0])
    return crad[0]

def double_closest_pole(vrblvl=0):
    """
    Returns a tuple with the real and imaginary part of the closest pole
    used in the predictor in double precision.
    """
    if vrblvl > 0:
        print('in double_closest_pole ...')
    phc = get_phcfun(vrblvl-1)
    aprc = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    pole = (c_double * 2)()
    cpole = pointer(pole)
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_closest_pole calls phc', end='')
    retval = phc(866, aprc, bbb, cpole, vrb)
    vals = cpole[:2]
    if vrblvl > 0:
        print(', return value :', retval)
        print('the closest pole :', (vals[0][0], vals[0][1]))
    return (vals[0][0], vals[0][1])

def double_double_closest_pole(vrblvl=0):
    """
    Returns a tuple with the real and imaginary part of the closest pole
    used in the predictor in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_closest_pole ...')
    phc = get_phcfun(vrblvl-1)
    aprc = pointer(c_int32(1))
    bbb = pointer(c_int32(0))
    pole = (c_double * 2)()
    cpole = pointer(pole)
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_double_closest_pole calls phc', end='')
    retval = phc(866, aprc, bbb, cpole, vrb)
    vals = cpole[:2]
    if vrblvl > 0:
        print(', return value :', retval)
        print('the closest pole :', (vals[0][0], vals[0][1]))
    return (vals[0][0], vals[0][1])

def quad_double_closest_pole(vrblvl=0):
    """
    Returns a tuple with the real and imaginary part of the closest pole
    used in the predictor in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_closest_pole ...')
    phc = get_phcfun(vrblvl-1)
    aprc = pointer(c_int32(2))
    bbb = pointer(c_int32(0))
    pole = (c_double * 2)()
    cpole = pointer(pole)
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> quad_double_closest_pole calls phc', end='')
    retval = phc(866, aprc, bbb, cpole, vrb)
    vals = cpole[:2]
    if vrblvl > 0:
        print(', return value :', retval)
        print('the closest pole :', (vals[0][0], vals[0][1]))
    return (vals[0][0], vals[0][1])

def double_series_coefficients(dim, vrblvl = 0):
    """
    Returns a list of lists with the coefficients of the series
    computed by the predictor in double precision.
    On entry in dim is the number of variables.
    """
    if vrblvl > 0:
        print('in double_series_coefficients, dim :', dim)
    deg = get_degree_of_numerator(vrblvl) + get_degree_of_denominator(vrblvl)
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 3)()
    apars[0] = 0
    apars[1] = 0
    apars[2] = 0
    apar = pointer(apars)
    bvrb = pointer(c_int32(vrblvl))
    coef = (c_double * 2)()
    ccff = pointer(coef)
    vrb = c_int32(vrblvl)
    result = []
    for row in range(1, dim+1):
        cfs = []
        for col in range(deg+1):
            apars[1] = row
            apars[2] = col
            if vrblvl > 0:
                print('-> double_series_coefficients calls phc', end='')
            retval = phc(869, apar, bvrb, ccff, vrb)
            if vrblvl > 0:
                print(', return value :', retval)
            vals = ccff[:2]
            rcf = vals[0][0]
            icf = vals[0][1]
            if vrblvl > 0:
                print('coefficients :', rcf, icf)
            cfs.append(complex(rcf, icf))
        result.append(cfs)
    return result

def double_double_series_coefficients(dim, vrblvl = 0):
    """
    Returns a list of lists with the coefficients of the series
    computed by the predictor in double double precision.
    The double coefficients are the highest parts of the double doubles.
    On entry in dim is the number of variables.
    """
    if vrblvl > 0:
        print('in double_double_series_coefficients, dim :', dim)
    deg = get_degree_of_numerator(vrblvl) + get_degree_of_denominator(vrblvl)
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 3)()
    apars[0] = 1
    apars[1] = 0
    apars[2] = 0
    apar = pointer(apars)
    bvrb = pointer(c_int32(vrblvl))
    coef = (c_double * 2)()
    ccff = pointer(coef)
    vrb = c_int32(vrblvl)
    result = []
    for row in range(1, dim+1):
        cfs = []
        for col in range(deg+1):
            apars[1] = row
            apars[2] = col
            if vrblvl > 0:
                print('-> double_double_series_coefficients calls phc', end='')
            retval = phc(869, apar, bvrb, ccff, vrb)
            if vrblvl > 0:
                print(', return value :', retval)
            vals = ccff[:2]
            rcf = vals[0][0]
            icf = vals[0][1]
            if vrblvl > 0:
                print('coefficients :', rcf, icf)
            cfs.append(complex(rcf, icf))
        result.append(cfs)
    return result

def quad_double_series_coefficients(dim, vrblvl = 0):
    """
    Returns a list of lists with the coefficients of the series
    computed by the predictor in quad double precision.
    The double coefficients are the highest parts of the quad doubles.
    On entry in dim is the number of variables.
    """
    if vrblvl > 0:
        print('in quad_double_series_coefficients, dim :', dim)
    deg = get_degree_of_numerator(vrblvl) + get_degree_of_denominator(vrblvl)
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 3)()
    apars[0] = 2
    apars[1] = 0
    apars[2] = 0
    apar = pointer(apars)
    bvrb = pointer(c_int32(vrblvl))
    coef = (c_double * 2)()
    ccff = pointer(coef)
    vrb = c_int32(vrblvl)
    result = []
    for row in range(1, dim+1):
        cfs = []
        for col in range(deg+1):
            apars[1] = row
            apars[2] = col
            if vrblvl > 0:
                print('-> quad_double_series_coefficients calls phc', end='')
            retval = phc(869, apar, bvrb, ccff, vrb)
            if vrblvl > 0:
                print(', return value :', retval)
            vals = ccff[:2]
            rcf = vals[0][0]
            icf = vals[0][1]
            if vrblvl > 0:
                print('coefficients :', rcf, icf)
            cfs.append(complex(rcf, icf))
        result.append(cfs)
    return result

def double_pade_coefficients(idx, vrblvl=0):
    """
    Returns a tuple of lists with the coefficients of the Pade approximants
    computed by the predictor in double precision.
    The first list in the tuple holds the coefficients of the numerator,
    the second list in the tuple holds the denominator coefficients.
    On entry in idx is the index of a variable.
    """
    if vrblvl > 0:
        print('in double_pade_coefficients, idx :', idx)
    numdeg = get_degree_of_numerator()
    dendeg = get_degree_of_denominator()
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 4)()
    apars[0] = 0 # precision
    apars[1] = 1
    apars[2] = 0
    apars[3] = 0
    apar = pointer(apars)
    bvrb = pointer(c_int32(vrblvl))
    coef = (c_double * 2)()
    ccff = pointer(coef)
    vrb = c_int32(vrblvl)
    numcfs = []
    for col in range(numdeg+1):
        apars[1] = 1
        apars[2] = idx
        apars[3] = col
        if vrblvl > 0:
            print('-> double_pade_coefficients calls phc', end='')
        retval = phc(870, apar, bvrb, ccff, vrb)
        if vrblvl > 0:
            print(', return value :', retval)
        vals = ccff[:2]
        rcf = vals[0][0]
        icf = vals[0][1]
        if vrblvl > 0:
            print('coefficients :', rcf, icf)
        numcfs.append(complex(rcf, icf))
    dencfs = []
    for col in range(dendeg+1):
        apars[1] = 0
        apars[2] = idx
        apars[3] = col
        if vrblvl > 0:
            print('-> double_pade_coefficients calls phc', end='')
        retval = phc(870, apar, bvrb, ccff, vrb)
        if vrblvl > 0:
            print(', return value :', retval)
        vals = ccff[:2]
        rcf = vals[0][0]
        icf = vals[0][1]
        if vrblvl > 0:
            print('coefficients :', rcf, icf)
        dencfs.append(complex(rcf, icf))
    return (numcfs, dencfs)

def double_double_pade_coefficients(idx, vrblvl=0):
    """
    Returns a tuple of lists with the coefficients of the Pade approximants
    computed by the predictor in double double precision.
    The first list in the tuple holds the coefficients of the numerator,
    the second list in the tuple holds the denominator coefficients.
    On entry in idx is the index of a variable.
    The double coefficients on return are the highest parts
    of the double doubles.
    """
    if vrblvl > 0:
        print('in double_double_pade_coefficients, idx :', idx)
    numdeg = get_degree_of_numerator()
    dendeg = get_degree_of_denominator()
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 4)()
    apars[0] = 1 # precision
    apars[1] = 1
    apars[2] = 0
    apars[3] = 0
    apar = pointer(apars)
    bvrb = pointer(c_int32(vrblvl))
    coef = (c_double * 2)()
    ccff = pointer(coef)
    vrb = c_int32(vrblvl)
    numcfs = []
    for col in range(numdeg+1):
        apars[1] = 1
        apars[2] = idx
        apars[3] = col
        if vrblvl > 0:
            print('-> double_double_pade_coefficients calls phc', end='')
        retval = phc(870, apar, bvrb, ccff, vrb)
        if vrblvl > 0:
            print(', return value :', retval)
        vals = ccff[:2]
        rcf = vals[0][0]
        icf = vals[0][1]
        if vrblvl > 0:
            print('coefficients :', rcf, icf)
        numcfs.append(complex(rcf, icf))
    dencfs = []
    for col in range(dendeg+1):
        apars[1] = 0
        apars[2] = idx
        apars[3] = col
        if vrblvl > 0:
            print('-> double_double_pade_coefficients calls phc', end='')
        retval = phc(870, apar, bvrb, ccff, vrb)
        if vrblvl > 0:
            print(', return value :', retval)
        vals = ccff[:2]
        rcf = vals[0][0]
        icf = vals[0][1]
        if vrblvl > 0:
            print('coefficients :', rcf, icf)
        dencfs.append(complex(rcf, icf))
    return (numcfs, dencfs)

def quad_double_pade_coefficients(idx, vrblvl=0):
    """
    Returns a tuple of lists with the coefficients of the Pade approximants
    computed by the predictor in quad double precision.
    The first list in the tuple holds the coefficients of the numerator,
    the second list in the tuple holds the denominator coefficients.
    On entry in idx is the index of a variable.
    The double coefficients on return are the highest parts
    of the quad doubles.
    """
    if vrblvl > 0:
        print('in quad_double_pade_coefficients, idx :', idx)
    numdeg = get_degree_of_numerator()
    dendeg = get_degree_of_denominator()
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 4)()
    apars[0] = 2 # precision
    apars[1] = 1
    apars[2] = 0
    apars[3] = 0
    apar = pointer(apars)
    bvrb = pointer(c_int32(vrblvl))
    coef = (c_double * 2)()
    ccff = pointer(coef)
    vrb = c_int32(vrblvl)
    numcfs = []
    for col in range(numdeg+1):
        apars[1] = 1
        apars[2] = idx
        apars[3] = col
        if vrblvl > 0:
            print('-> quad_double_pade_coefficients calls phc', end='')
        retval = phc(870, apar, bvrb, ccff, vrb)
        if vrblvl > 0:
            print(', return value :', retval)
        vals = ccff[:2]
        rcf = vals[0][0]
        icf = vals[0][1]
        if vrblvl > 0:
            print('coefficients :', rcf, icf)
        numcfs.append(complex(rcf, icf))
    dencfs = []
    for col in range(dendeg+1):
        apars[1] = 0
        apars[2] = idx
        apars[3] = col
        if vrblvl > 0:
            print('-> quad_double_pade_coefficients calls phc', end='')
        retval = phc(870, apar, bvrb, ccff, vrb)
        if vrblvl > 0:
            print(', return value :', retval)
        vals = ccff[:2]
        rcf = vals[0][0]
        icf = vals[0][1]
        if vrblvl > 0:
            print('coefficients :', rcf, icf)
        dencfs.append(complex(rcf, icf))
    return (numcfs, dencfs)

def double_pade_vector(dim, vrblvl=0):
    """
    Returns the list of all coefficients over all dim variables,
    computed by the predictor in double precision.
    """
    result = []
    for i in range(1, dim+1):
        result.append(double_pade_coefficients(i, vrblvl))
    return result

def double_double_pade_vector(dim, vrblvl=0):
    """
    Returns the list of all coefficients over all dim variables,
    computed by the predictor in double double precision.
    """
    result = []
    for i in range(1, dim+1):
        result.append(double_double_pade_coefficients(i, vrblvl))
    return result

def quad_double_pade_vector(dim, vrblvl=0):
    """
    Returns the list of all coefficients over all dim variables,
    computed by the predictor in quad double precision.
    """
    result = []
    for i in range(1, dim+1):
        result.append(quad_double_pade_coefficients(i, vrblvl))
    return result

def double_poles(dim, vrblvl=0):
    """
    Returns a list of lists of all poles of the vector of length dim,
    computed by the predictor in double precision.
    """
    if vrblvl > 0:
        print('in double_poles, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 3)()
    apars[0] = 0 # precision
    apars[1] = 0
    apars[2] = 0
    apar = pointer(apars)
    bvrb = pointer(c_int32(vrblvl))
    coef = (c_double * 2)()
    ccff = pointer(coef)
    vrb = c_int32(vrblvl)
    dendeg = get_degree_of_denominator()
    result = []
    for i in range(1, dim+1):
        poles = []
        for j in range(1, dendeg+1):
            apars[1] = i
            apars[2] = j
            if vrblvl > 0:
                print('-> double_poles calls phc', end='')
            retval = phc(871, apar, bvrb, ccff, vrb)
            if vrblvl > 0:
                print(', return value :', retval)
            vals = ccff[:2]
            repole = vals[0][0]
            impole = vals[0][1]
            if vrblvl > 0:
                print('pole coordinates :', repole, impole)
            poles.append(complex(repole, impole))
        result.append(poles)
    return result

def double_double_poles(dim, vrblvl=0):
    """
    Returns a list of lists of all poles of the vector of length dim,
    computed by the predictor in double double precision.
    The doubles on return are the highest parts of the double doubles.
    """
    if vrblvl > 0:
        print('in double_double_poles, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 3)()
    apars[0] = 1 # precision
    apars[1] = 0
    apars[2] = 0
    apar = pointer(apars)
    bvrb = pointer(c_int32(vrblvl))
    coef = (c_double * 2)()
    ccff = pointer(coef)
    vrb = c_int32(vrblvl)
    dendeg = get_degree_of_denominator()
    result = []
    for i in range(1, dim+1):
        poles = []
        for j in range(1, dendeg+1):
            apars[1] = i
            apars[2] = j
            if vrblvl > 0:
                print('-> double_double_poles calls phc', end='')
            retval = phc(871, apar, bvrb, ccff, vrb)
            if vrblvl > 0:
                print(', return value :', retval)
            vals = ccff[:2]
            repole = vals[0][0]
            impole = vals[0][1]
            if vrblvl > 0:
                print('pole coordinates :', repole, impole)
            poles.append(complex(repole, impole))
        result.append(poles)
    return result

def quad_double_poles(dim, vrblvl=0):
    """
    Returns a list of lists of all poles of the vector of length dim,
    computed by the predictor in quad double precision.
    The doubles on return are the highest parts of the quad doubles.
    """
    if vrblvl > 0:
        print('in quad_double_poles, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    apars = (c_int32 * 3)()
    apars[0] = 2 # precision
    apars[1] = 0
    apars[2] = 0
    apar = pointer(apars)
    bvrb = pointer(c_int32(vrblvl))
    coef = (c_double * 2)()
    ccff = pointer(coef)
    vrb = c_int32(vrblvl)
    dendeg = get_degree_of_denominator()
    result = []
    for i in range(1, dim+1):
        poles = []
        for j in range(1, dendeg+1):
            apars[1] = i
            apars[2] = j
            if vrblvl > 0:
                print('-> quad_double_poles calls phc', end='')
            retval = phc(871, apar, bvrb, ccff, vrb)
            if vrblvl > 0:
                print(', return value :', retval)
            vals = ccff[:2]
            repole = vals[0][0]
            impole = vals[0][1]
            if vrblvl > 0:
                print('pole coordinates :', repole, impole)
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
    (nuq, deq) = cff
    num = ['+' + str(nuq[k]) + (f"*t**{k}") for k in range(len(nuq))]
    numerator = ''.join(num)
    den = ['+' + str(deq[k]) + (f"*t**{k}") for k in range(len(deq))]
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

def clear_double_data(vrblvl=0):
    """
    Clears the data used by the tracker in double precision.
    """
    if vrblvl > 0:
        print('in clear_double_data ...')
    phc = get_phcfun(vrblvl-1)
    aprc = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> clear_double_data calls phc', end='')
    retval = phc(864, aprc, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_double_double_data(vrblvl=0):
    """
    Clears the data used by the tracker in double double precision.
    """
    if vrblvl > 0:
        print('in clear_double_double_data ...')
    phc = get_phcfun(vrblvl-1)
    aprc = pointer(c_int32(1))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> clear_double_double_data calls phc', end='')
    retval = phc(864, aprc, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def clear_quad_double_data(vrblvl=0):
    """
    Clears the data used by the tracker in quad double precision.
    """
    if vrblvl > 0:
        print('in clear_quad_double_data ...')
    phc = get_phcfun(vrblvl-1)
    aprc = pointer(c_int32(2))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> clear_quad_double_data calls phc', end='')
    retval = phc(864, aprc, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def next_double_track(target, start, sols, homogeneous=False, \
    interactive=False, vrblvl=0):
    r"""
    Runs the series-Pade tracker step by step in double precision,
    for an artificial-parameter homotopy.
    On input are a target system and a start system with solutions.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in *sols*.
    The *sols* is a list of strings representing start solutions.
    Prompts for each predictor-corrector step, if *interactive*.
    If *vrblvl* > 0, then extra output is written.
    If *homogeneous*, then path tracking happens in projective space,
    otherwise the original affine coordinates are used.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    if vrblvl > 0:
        print('in next_double_track, homogeneous :', homogeneous)
        print('interactive :', interactive)
        print('the target system :')
        for pol in target:
            print(pol)
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    result = []
    dim = number_of_symbols(start, vrblvl-1)
    initialize_double_artificial_homotopy(target, start, homogeneous, vrblvl-1)
    fmt = 'pole step : %.3e, estimated distance : %.3e, Hessian step : %.3e'
    for (idx, sol) in enumerate(sols):
        tval = 0.0
        if vrblvl > 0:
            print('tracking solution path', idx+1, '...')
        set_double_solution(dim, sol, vrblvl-1)
        for step in range(1, 101):
            if interactive:
                answer = input('next predictor-corrector step ? (y/n) ')
                if answer != 'y':
                    result.append(sol)
                    break
            else:
                if tval >= 1.0:
                    result.append(sol)
                    break
            double_predict_correct(vrblvl-1)
            polestep = double_pole_step(vrblvl-1)
            estidist = double_estimated_distance(vrblvl-1)
            curvstep = double_hessian_step(vrblvl-1)
            if vrblvl > 0:
                print('step', step, end=', ')
                print(fmt % (polestep, estidist, curvstep))
            previoustval = tval
            tval = double_t_value(vrblvl-1)
            step = double_step_size(vrblvl-1)
            frp = double_pole_radius(vrblvl-1)
            if vrblvl > 0:
                print(f't : {tval:.3e}, step : {step:.3e}, frp : {frp:.3e}')
            cfp = double_closest_pole(vrblvl-1)
            if vrblvl > 0:
                print('For the previous t value', previoustval, ':')
                print('1) closest pole : ', cfp)
                print('2) the series:', \
                    double_series_coefficients(dim, vrblvl-1))
                print('3) Pade vector:', double_pade_vector(dim, vrblvl-1))
                print('4) poles:', double_poles(dim, vrblvl-1))
            predsol = get_double_predicted_solution(vrblvl-1)
            sol = get_double_solution(vrblvl-1)
            if vrblvl > 0:
                print('the predicted solution :')
                print(predsol)
                print('the current solution :')
                print(sol)
    gamma = get_gamma_constant(vrblvl-1)
    clear_double_solutions(vrblvl-1)
    clear_double_data(vrblvl-1)
    return (gamma, result)

def next_double_double_track(target, start, sols, homogeneous=False, \
    interactive=False, vrblvl=0):
    r"""
    Runs the series-Pade tracker step by step in double double precision,
    for an artificial-parameter homotopy.
    On input are a target system and a start system with solutions.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in *sols*.
    The *sols* is a list of strings representing start solutions.
    Prompts for each predictor-corrector step, if *interactive*.
    If *vrblvl* > 0, then extra output is written.
    If *homogeneous*, then path tracking happens in projective space,
    otherwise the original affine coordinates are used.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    if vrblvl > 0:
        print('in next_double_double_track, homogeneous :', homogeneous)
        print('interactive :', interactive)
        print('the target system :')
        for pol in target:
            print(pol)
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    result = []
    dim = number_of_symbols(start, vrblvl-1)
    initialize_double_double_artificial_homotopy(target, start, \
        homogeneous, vrblvl-1)
    gamma = get_gamma_constant(vrblvl-1)
    fmt = 'pole step : %.3e, estimated distance : %.3e, Hessian step : %.3e'
    for (idx, sol) in enumerate(sols):
        tval = 0.0
        if vrblvl > 0:
            print('tracking solution path', idx+1, '...')
        set_double_double_solution(dim, sol, vrblvl-1)
        for step in range(1, 101):
            if interactive:
                answer = input('next predictor-corrector step ? (y/n) ')
                if answer != 'y':
                    result.append(sol)
                    break
            else:
                if tval >= 1.0:
                    result.append(sol)
                    break
            double_double_predict_correct(vrblvl-1)
            polestep = double_double_pole_step(vrblvl-1)
            estidist = double_double_estimated_distance(vrblvl-1)
            curvstep = double_double_hessian_step(vrblvl-1)
            if vrblvl > 0:
                print('step', step, end=', ')
                print(fmt % (polestep, estidist, curvstep))
            previoustval = tval
            tval = double_double_t_value(vrblvl-1)
            step = double_double_step_size(vrblvl-1)
            frp = double_double_pole_radius(vrblvl-1)
            if vrblvl > 0:
                print(f't : {tval:.3e}, step : {step:.3e}, frp : {frp:.3e}')
            cfp = double_double_closest_pole(vrblvl-1)
            if vrblvl > 0:
                print('For the previous t value', previoustval, ':')
                print('1) closest pole : ', cfp)
                print('2) the series:', \
                    double_double_series_coefficients(dim, vrblvl-1))
                print('3) Pade vector:', \
                    double_double_pade_vector(dim, vrblvl-1))
                print('4) poles:', double_double_poles(dim, vrblvl-1))
            predsol = get_double_double_predicted_solution(vrblvl-1)
            sol = get_double_double_solution(vrblvl-1)
            if vrblvl > 0:
                print('the predicted solution :')
                print(predsol)
                print('the current solution :')
                print(sol)
    clear_double_double_solutions(vrblvl-1)
    clear_double_double_data(vrblvl-1)
    return (gamma, result)

def next_quad_double_track(target, start, sols, homogeneous=False, \
    interactive=False, vrblvl=0):
    r"""
    Runs the series-Pade tracker step by step in quad double precision,
    for an artificial-parameter homotopy.
    On input are a target system and a start system with solutions.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in *sols*.
    The *sols* is a list of strings representing start solutions.
    Prompts for each predictor-corrector step, if *interactive*.
    If *vrblvl* > 0, then extra output is written.
    If *homogeneous*, then path tracking happens in projective space,
    otherwise the original affine coordinates are used.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    if vrblvl > 0:
        print('in next_quad_double_track, homogeneous :', homogeneous)
        print('interactive :', interactive)
        print('the target system :')
        for pol in target:
            print(pol)
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    result = []
    dim = number_of_symbols(start, vrblvl-1)
    initialize_quad_double_artificial_homotopy(target, start, \
        homogeneous, vrblvl-1)
    gamma = get_gamma_constant(vrblvl-1)
    fmt = 'pole step : %.3e, estimated distance : %.3e, Hessian step : %.3e'
    for (idx, sol) in enumerate(sols):
        tval = 0.0
        if vrblvl > 0:
            print('tracking solution path', idx+1, '...')
        set_quad_double_solution(dim, sol, vrblvl-1)
        for step in range(1, 101):
            if interactive:
                answer = input('next predictor-corrector step ? (y/n) ')
                if answer != 'y':
                    result.append(sol)
                    break
            else:
                if tval >= 1.0:
                    result.append(sol)
                    break
            quad_double_predict_correct(vrblvl-1)
            polestep = quad_double_pole_step(vrblvl-1)
            estidist = quad_double_estimated_distance(vrblvl-1)
            curvstep = quad_double_hessian_step(vrblvl-1)
            print(fmt % (polestep, estidist, curvstep))
            previoustval = tval
            tval = quad_double_t_value(vrblvl-1)
            step = quad_double_step_size(vrblvl-1)
            frp = quad_double_pole_radius(vrblvl-1)
            if vrblvl > 0:
                print('step', step, end=', ')
                print(f't : {tval:.3e}, step : {step:.3e}, frp : {frp:.3e}')
            cfp = quad_double_closest_pole(vrblvl-1)
            if vrblvl > 0:
                print('For the previous t value', previoustval, ':')
                print('1) closest pole : ', cfp)
                print('2) the series:', \
                    quad_double_series_coefficients(dim, vrblvl-1))
                print('3) Pade vector:', \
                    quad_double_pade_vector(dim, vrblvl-1))
                print('4) poles:', quad_double_poles(dim, vrblvl-1))
            predsol = get_quad_double_predicted_solution(vrblvl-1)
            sol = get_quad_double_solution(vrblvl-1)
            if vrblvl > 0:
                print('the predicted solution :')
                print(predsol)
                print('the current solution :')
                print(sol)
    clear_quad_double_solutions(vrblvl-1)
    clear_quad_double_data(vrblvl-1)
    return (gamma, result)

def next_double_loop(hom, idx, sols, interactive=False, vrblvl=0):
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
    Prompts before each predictor-corrector step, if *interactive*.
    If vrblvl > 0, then extra output is written.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    if vrblvl > 0:
        print('in next_double_loop, idx :', idx, end='')
        print(', interactive :', interactive)
        print('the homotopy :')
        for pol in hom:
            print(pol)
        print('the solutions :')
        for (solidx, sol) in enumerate(sols):
            print('Solution', solidx+1, ':')
            print(sol)
    result = []
    dim = number_of_symbols(hom) - 1
    initialize_double_parameter_homotopy(hom, idx, vrblvl-1)
    gamma = get_gamma_constant(vrblvl-1)
    fmt = 'pole step : %.3e, estimated distance : %.3e, Hessian step : %.3e'
    for (solidx, sol) in enumerate(sols):
        tval = 0.0
        if vrblvl > 0:
            print('tracking solution path', solidx+1, '...')
        set_double_solution(dim, sol, vrblvl-1)
        for step in range(1, 101):
            if interactive:
                answer = input('next predictor-corrector step ? (y/n) ')
                if answer != 'y':
                    result.append(sol)
                    break
            else:
                if tval >= 1.0:
                    result.append(sol)
                    break
            double_predict_correct(vrblvl-1)
            polestep = double_pole_step(vrblvl-1)
            estidist = double_estimated_distance(vrblvl-1)
            curvstep = double_hessian_step(vrblvl-1)
            if vrblvl > 0:
                print('step', step, end=', ')
                print(fmt % (polestep, estidist, curvstep))
            previoustval = tval
            tval = double_t_value(vrblvl-1)
            step = double_step_size(vrblvl-1)
            frp = double_pole_radius(vrblvl-1)
            if vrblvl > 0:
                print(f't : {tval:.3e}, step : {step:.3e}, frp : {frp:.3e}')
            cfp = double_closest_pole(vrblvl-1)
            if vrblvl > 0:
                print('For the previous t value', previoustval, ':')
                print('1) closest pole : ', cfp)
                print('2) the series:', \
                    double_series_coefficients(dim, vrblvl-1))
                print('3) Pade vector:', double_pade_vector(dim, vrblvl-1))
                print('4) poles:', double_poles(dim, vrblvl-1))
            predsol = get_double_predicted_solution(vrblvl-1)
            sol = get_double_solution(vrblvl-1)
            if vrblvl > 0:
                print('the predicted solution :')
                print(predsol)
                print('the current solution :')
                print(sol)
    clear_double_solutions(vrblvl-1)
    clear_double_data(vrblvl-1)
    return (gamma, result)

def next_double_double_loop(hom, idx, sols, interactive=False, vrblvl=0):
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
    Prompts before each predictor-corrector step, if *interactive*.
    If vrblvl > 0, then extra output is written.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    if vrblvl > 0:
        print('in next_double_double_loop, idx :', idx, end='')
        print(', interactive :', interactive)
        print('the homotopy :')
        for pol in hom:
            print(pol)
        print('the solutions :')
        for (solidx, sol) in enumerate(sols):
            print('Solution', solidx+1, ':')
            print(sol)
    result = []
    dim = number_of_symbols(hom) - 1
    initialize_double_double_parameter_homotopy(hom, idx, vrblvl-1)
    gamma = get_gamma_constant(vrblvl-1)
    fmt = 'pole step : %.3e, estimated distance : %.3e, Hessian step : %.3e'
    for (solidx, sol) in enumerate(sols):
        tval = 0.0
        if vrblvl > 0:
            print('tracking solution path', solidx+1, '...')
        set_double_double_solution(dim, sol, vrblvl-1)
        for step in range(1, 101):
            if interactive:
                answer = input('next predictor-corrector step ? (y/n) ')
                if answer != 'y':
                    result.append(sol)
                    break
            else:
                if tval >= 1.0:
                    result.append(sol)
                    break
            double_double_predict_correct(vrblvl-1)
            polestep = double_double_pole_step(vrblvl-1)
            estidist = double_double_estimated_distance(vrblvl-1)
            curvstep = double_double_hessian_step(vrblvl-1)
            if vrblvl > 0:
                print('step', step, end=', ')
                print(fmt % (polestep, estidist, curvstep))
            previoustval = tval
            tval = double_double_t_value(vrblvl-1)
            step = double_double_step_size(vrblvl-1)
            frp = double_double_pole_radius(vrblvl-1)
            if vrblvl > 0:
                print(f't : {tval:.3e}, step : {step:.3e}, frp : {frp:.3e}')
            cfp = double_double_closest_pole(vrblvl-1)
            if vrblvl > 0:
                print('For the previous t value', previoustval, ':')
                print('1) closest pole : ', cfp)
                print('2) the series:', \
                    double_double_series_coefficients(dim, vrblvl-1))
                print('3) Pade vector:', \
                    double_double_pade_vector(dim, vrblvl-1))
                print('4) poles:', double_double_poles(dim, vrblvl-1))
                predsol = get_double_double_predicted_solution(vrblvl-1)
                sol = get_double_double_solution(vrblvl-1)
                print('the predicted solution :')
                print(predsol)
                print('the current solution :')
                print(sol)
    clear_double_double_solutions(vrblvl-1)
    clear_double_double_data(vrblvl-1)
    return (gamma, result)

def next_quad_double_loop(hom, idx, sols, interactive=False, vrblvl=0):
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
    Prompts before each predictor-corrector step, if *interactive*.
    If vrblvl > 0, then extra output is written.
    On return are the string representations of the solutions
    computed at the end of the paths.
    """
    if vrblvl > 0:
        print('in next_quad_double_loop, idx :', idx, end='')
        print(', interactive :', interactive)
        print('the homotopy :')
        for pol in hom:
            print(pol)
        print('the solutions :')
        for (solidx, sol) in enumerate(sols):
            print('Solution', solidx+1, ':')
            print(sol)
    result = []
    dim = number_of_symbols(hom) - 1
    initialize_quad_double_parameter_homotopy(hom, idx, vrblvl-1)
    gamma = get_gamma_constant(vrblvl-1)
    fmt = 'pole step : %.3e, estimated distance : %.3e, Hessian step : %.3e'
    for (solidx, sol) in enumerate(sols):
        tval = 0.0
        if vrblvl > 0:
            print('tracking solution path', solidx+1, '...')
        set_quad_double_solution(dim, sol, vrblvl-1)
        for step in range(1, 101):
            if interactive:
                answer = input('next predictor-corrector step ? (y/n) ')
                if answer != 'y':
                    result.append(sol)
                    break
            else:
                if tval >= 1.0:
                    result.append(sol)
                    break
            quad_double_predict_correct(vrblvl-1)
            polestep = quad_double_pole_step(vrblvl-1)
            estidist = quad_double_estimated_distance(vrblvl-1)
            curvstep = quad_double_hessian_step(vrblvl-1)
            if vrblvl > 0:
                print('step', step, ':')
                print(fmt % (polestep, estidist, curvstep))
            previoustval = tval
            tval = quad_double_t_value(vrblvl-1)
            step = quad_double_step_size(vrblvl-1)
            frp = quad_double_pole_radius(vrblvl-1)
            if vrblvl > 0:
                print(f't : {tval:.3e}, step : {step:.3e}, frp : {frp:.3e}')
            cfp = quad_double_closest_pole(vrblvl-1)
            if vrblvl > 0:
                print('For the previous t value', previoustval, ':')
                print('1) closest pole : ', cfp)
                print('2) the series:', \
                    quad_double_series_coefficients(dim, vrblvl-1))
                print('3) Pade vector:', \
                    quad_double_pade_vector(dim, vrblvl-1))
                print('4) poles:', quad_double_poles(dim, vrblvl-1))
            predsol = get_quad_double_predicted_solution(vrblvl-1)
            sol = get_quad_double_solution(vrblvl-1)
            if vrblvl > 0:
                print('the predicted solution :')
                print(predsol)
                print('the current solution :')
                print(sol)
    clear_quad_double_solutions(vrblvl-1)
    clear_quad_double_data(vrblvl-1)
    return (gamma, result)

def test_tuning(vrblvl=0):
    """
    Tests the tuning of the parameters.
    """
    if vrblvl > 0:
        print('in test_tuning ...')
    set_gamma_constant(complex(1.2, 3.4), vrblvl)
    gamma = get_gamma_constant(vrblvl)
    print('the gamma constant :', gamma)
    set_degree_of_numerator(4, vrblvl)
    deg = get_degree_of_numerator(vrblvl)
    print('degree of numerator :', deg)
    set_degree_of_denominator(2, vrblvl)
    deg = get_degree_of_denominator(vrblvl)
    print('degree of denominator :', deg)
    set_maximum_step_size(0.2, vrblvl)
    maxstep = get_maximum_step_size(vrblvl)
    print('the maximum step size :', maxstep)
    set_minimum_step_size(1.0e-12, vrblvl)
    minstep = get_minimum_step_size(vrblvl)
    print('the minimum step size :', minstep)
    set_pole_radius_beta_factor(0.4, vrblvl)
    beta = get_pole_radius_beta_factor(vrblvl)
    print('pole radius beta factor :', beta)
    set_curvature_beta_factor(0.3, vrblvl)
    beta = get_curvature_beta_factor(vrblvl)
    print('curvature beta factor :', beta)
    set_predictor_residual_alpha(1.0e-8, vrblvl)
    alpha = get_predictor_residual_alpha(vrblvl)
    print('predictor residual alpha :', alpha)
    set_corrector_residual_tolerance(1.0e-12, vrblvl)
    tol = get_corrector_residual_tolerance(vrblvl)
    print('corrector residual tolerance :', tol)
    set_zero_series_coefficient_tolerance(1.0e-10, vrblvl)
    tol = get_zero_series_coefficient_tolerance(vrblvl)
    print('zero series coefficient tolerance :', tol)
    set_maximum_corrector_steps(2, vrblvl)
    maxsteps = get_maximum_corrector_steps(vrblvl)
    print('maximum corrector steps :', maxsteps)
    set_maximum_steps_on_path(1000, vrblvl)
    maxsteps = get_maximum_steps_on_path(vrblvl)
    print('maximum steps on path :', maxsteps)
    write_parameters(vrblvl)
    set_default_parameters(vrblvl)
    write_parameters(vrblvl)
    return 0

def test_double_track(vrblvl=0):
    """
    Runs on the mickey mouse example of two quadrics,
    in double precision.
    """
    if vrblvl > 0:
        print('in test_double_track ...')
    mickey = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    start, startsols = total_degree_start_system(mickey, vrblvl=vrblvl)
    if vrblvl > 0:
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    gamma, sols = double_track(mickey, start, startsols, vrblvl=vrblvl)
    if vrblvl > 0:
        print('gamma :', gamma)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(mickey, sols, vrblvl)
    if vrblvl > 0:
        print('the error sum :', err)
    if len(sols) == 4 and abs(err.real + err.imag) < 1.0e-10:
        if vrblvl > 0:
            print('Found 4 solutions and error is okay.')
        return 0
    if len(sols) != 4:
        if vrblvl > 0:
            print('Number of solutions is not 4 :', len(sols))
        return 1
    if abs(err.real + err.imag) >= 1.0e-10:
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def test_double_double_track(vrblvl=0):
    """
    Runs on the mickey mouse example of two quadrics,
    in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_track ...')
    mickey = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    start, startsols = total_degree_start_system(mickey, vrblvl=vrblvl)
    if vrblvl > 0:
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    gamma, sols = double_double_track(mickey, start, startsols, \
        vrblvl=vrblvl)
    if vrblvl > 0:
        print('gamma :', gamma)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(mickey, sols, vrblvl)
    if vrblvl > 0:
        print('the error sum :', err)
    if len(sols) == 4 and abs(err.real + err.imag) < 1.0e-10:
        if vrblvl > 0:
            print('Found 4 solutions and error is okay.')
        return 0
    if len(sols) != 4:
        if vrblvl > 0:
            print('Number of solutions is not 4 :', len(sols))
        return 1
    if abs(err.real + err.imag) >= 1.0e-10:
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def test_quad_double_track(vrblvl=0):
    """
    Runs on the mickey mouse example of two quadrics,
    in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_track ...')
    mickey = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    start, startsols = total_degree_start_system(mickey, vrblvl=vrblvl)
    if vrblvl > 0:
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    gamma, sols = quad_double_track(mickey, start, startsols, \
        vrblvl=vrblvl)
    if vrblvl > 0:
        print('gamma :', gamma)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(mickey, sols, vrblvl)
    if vrblvl > 0:
        print('the error sum :', err)
    if len(sols) == 4 and abs(err.real + err.imag) < 1.0e-10:
        if vrblvl > 0:
            print('Found 4 solutions and error is okay.')
        return 0
    if len(sols) != 4:
        if vrblvl > 0:
            print('Number of solutions is not 4 :', len(sols))
        return 1
    if abs(err.real + err.imag) >= 1.0e-10:
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def test_next_double_track(vrblvl=0):
    """
    Runs on the mickey mouse example of two quadrics,
    with a step-by-step tracker in double precision.
    """
    if vrblvl > 0:
        print('in test_next_double_track ...')
    mickey = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    start, startsols = total_degree_start_system(mickey, vrblvl=vrblvl-1)
    if vrblvl > 0:
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    gamma, sols = next_double_track(mickey, start, startsols, \
        interactive=False, vrblvl=vrblvl-1)
    if vrblvl > 0:
        print('gamma :', gamma)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(mickey, sols, vrblvl-1)
    if vrblvl > 0:
        print('the error sum :', err)
    if len(sols) == 4 and abs(err.real + err.imag) < 1.0e-10:
        if vrblvl > 0:
            print('Found 4 solutions and error is okay.')
        return 0
    if len(sols) != 4:
        if vrblvl > 0:
            print('Number of solutions is not 4 :', len(sols))
        return 1
    if abs(err.real + err.imag) >= 1.0e-10:
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def test_next_double_double_track(vrblvl=0):
    """
    Runs on the mickey mouse example of two quadrics,
    with a step-by-step tracker in double double precision.
    """
    if vrblvl > 0:
        print('in test_next_double_double_track ...')
    mickey = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    start, startsols = total_degree_start_system(mickey, vrblvl=vrblvl)
    if vrblvl > 0:
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    gamma, sols = next_double_double_track(mickey, start, startsols, \
        vrblvl=vrblvl)
    if vrblvl > 0:
        print('gamma :', gamma)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(mickey, sols, vrblvl)
    if vrblvl > 0:
        print('the error sum :', err)
    if len(sols) == 4 and abs(err.real + err.imag) < 1.0e-10:
        if vrblvl > 0:
            print('Found 4 solutions and error is okay.')
        return 0
    if len(sols) != 4:
        if vrblvl > 0:
            print('Number of solutions is not 4 :', len(sols))
        return 1
    if abs(err.real + err.imag) >= 1.0e-10:
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def test_next_quad_double_track(vrblvl=0):
    """
    Runs on the mickey mouse example of two quadrics,
    with a step-by-step tracker in quad double precision.
    """
    if vrblvl > 0:
        print('in test_next_quad_double_track ...')
    mickey = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    start, startsols = total_degree_start_system(mickey, vrblvl=vrblvl)
    if vrblvl > 0:
        print('the start system :')
        for pol in start:
            print(pol)
        print('the start solutions :')
        for (idx, sol) in enumerate(startsols):
            print('Solution', idx+1, ':')
            print(sol)
    gamma, sols = next_quad_double_track(mickey, start, startsols, \
        vrblvl=vrblvl)
    if vrblvl > 0:
        print('gamma :', gamma)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    err = verify(mickey, sols, vrblvl)
    if vrblvl > 0:
        print('the error sum :', err)
    if len(sols) == 4 and abs(err.real + err.imag) < 1.0e-10:
        if vrblvl > 0:
            print('Found 4 solutions and error is okay.')
        return 0
    if len(sols) != 4:
        if vrblvl > 0:
            print('Number of solutions is not 4 :', len(sols))
        return 1
    if abs(err.real + err.imag) >= 1.0e-10:
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def test_double_hyperbola(vrblvl=0):
    """
    Tests the step-by-step Pade tracker on a hyperbola,
    in double precision.
    """
    if vrblvl > 0:
        print('in test_double_hyperbola ...')
    par = 0.1
    xtp = ['x^2 - (t - 0.5)^2 - 0.01;']
    solx = sqrt(4*par**2+1)/2
    if vrblvl > 0:
        print('\nvalue of the first start solution :', solx)
    sol1 = make_solution(['x'], [solx])
    if vrblvl > 0:
        print('the first start solution :\n', sol1)
    sol2 = make_solution(['x'], [-solx])
    if vrblvl > 0:
        print('the second start solution :\n', sol2)
        print('tracking in double precision ...')
    gamma, sols = next_double_loop(xtp, 2, [sol1, sol2], \
        interactive=False, vrblvl=vrblvl)
    if vrblvl > 0:
        print('gamma :', gamma)
        print('the solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    endxtp = [xtp[0].replace('t', '1.0')]
    err = verify(endxtp, sols, vrblvl)
    if vrblvl > 0:
        print('the error sum :', err)
    if len(sols) == 2 and abs(err.real + err.imag) < 1.0e-10:
        if vrblvl > 0:
            print('Found 2 solutions and error is okay.')
        return 0
    if len(sols) != 2:
        if vrblvl > 0:
            print('Number of solutions is not 2 :', len(sols))
        return 1
    if abs(err.real + err.imag) >= 1.0e-10:
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def test_double_double_hyperbola(vrblvl=0):
    """
    Tests the step-by-step Pade tracker on a hyperbola,
    in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_hyperbola ...')
    par = 0.1
    xtp = ['x^2 - (t - 0.5)^2 - 0.01;']
    solx = sqrt(4*par**2+1)/2
    if vrblvl > 0:
        print('\nvalue of the first start solution :', solx)
    sol1 = make_solution(['x'], [solx])
    if vrblvl > 0:
        print('the first start solution :\n', sol1)
    sol2 = make_solution(['x'], [-solx])
    if vrblvl > 0:
        print('the second start solution :\n', sol2)
        print('tracking in double double precision ...')
    gamma, sols = next_double_double_loop(xtp, 2, [sol1, sol2], \
        interactive=False, vrblvl=vrblvl)
    if vrblvl > 0:
        print('gamma :', gamma)
    endxtp = [xtp[0].replace('t', '1.0')]
    err = verify(endxtp, sols, vrblvl)
    if vrblvl > 0:
        print('the error sum :', err)
    if len(sols) == 2 and abs(err.real + err.imag) < 1.0e-10:
        if vrblvl > 0:
            print('Found 2 solutions and error is okay.')
        return 0
    if len(sols) != 2:
        if vrblvl > 0:
            print('Number of solutions is not 2 :', len(sols))
        return 1
    if abs(err.real + err.imag) >= 1.0e-10:
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def test_quad_double_hyperbola(vrblvl=0):
    """
    Tests the step-by-step Pade tracker on a hyperbola,
    in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_hyperbola ...')
    par = 0.1
    xtp = ['x^2 - (t - 0.5)^2 - 0.01;']
    solx = sqrt(4*par**2+1)/2
    if vrblvl > 0:
        print('\nvalue of the first start solution :', solx)
    sol1 = make_solution(['x'], [solx])
    if vrblvl > 0:
        print('the first start solution :\n', sol1)
    sol2 = make_solution(['x'], [-solx])
    if vrblvl > 0:
        print('the second start solution :\n', sol2)
        print('tracking in quad double precision ...')
    gamma, sols = next_quad_double_loop(xtp, 2, [sol1, sol2], \
        interactive=False, vrblvl=vrblvl)
    if vrblvl > 0:
        print('gamma :', gamma)
    endxtp = [xtp[0].replace('t', '1.0')]
    err = verify(endxtp, sols, vrblvl)
    if vrblvl > 0:
        print('the error sum :', err)
    if len(sols) == 2 and abs(err.real + err.imag) < 1.0e-10:
        if vrblvl > 0:
            print('Found 2 solutions and error is okay.')
        return 0
    if len(sols) != 2:
        if vrblvl > 0:
            print('Number of solutions is not 2 :', len(sols))
        return 1
    if abs(err.real + err.imag) >= 1.0e-10:
        if vrblvl > 0:
            print('The error is too large.')
    return 1

def main():
    """
    Runs some tests on tuning and tracking.
    """
    lvl = 1
    fail = test_tuning(lvl)
    fail = fail + test_double_track(lvl)
    fail = fail + test_double_double_track(lvl)
    fail = fail + test_quad_double_track(lvl)
    fail = fail + test_next_double_track(lvl)
    fail = fail + test_next_double_double_track(lvl)
    fail = fail + test_next_quad_double_track(lvl)
    fail = fail + test_double_hyperbola(lvl)
    fail = fail + test_double_double_hyperbola(lvl)
    fail = fail + test_quad_double_hyperbola(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()
