"""
The module curves exports functions to approximate algebraic space curves
with rational expressions, in particular Pade approximants,
for use in a path tracker with apriori step size control.
"""
from ctypes import c_int32, c_double, pointer
from phcpy.version import get_phcfun
from phcpy.solutions import verify
from phcpy.homotopies import total_degree_start_system
from phcpy.trackers import set_double_target_system
from phcpy.trackers import set_double_start_system
from phcpy.trackers import set_double_start_solutions
from phcpy.trackers import get_double_solutions
from phcpy.trackers import set_double_homotopy

def set_default_parameters(vrblvl=0):
    """
    Sets the default values of the parameters.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in set_default_parameters ...')
    phc = get_phcfun()
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
    phc = get_phcfun()
    aidx = pointer(c_int32(idx))
    bval = pointer(c_int32(0))
    cpar = (c_double * 2)()
    if idx == 1:
        cpar[0] = c_double(value.real)
        cpar[1] = c_double(value.imag)
    else:
        if((idx == 2) or (idx == 3) or (idx == 11) or (idx == 12)):
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
    phc = get_phcfun()
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
    else:
        if((idx == 2) or (idx == 3) or (idx == 11) or (idx == 12)):
            return bval[0]
        else:
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

def double_track(target, start, startsols, vrblvl=0):
    r"""
    Does path tracking in double precision.
    On input are a target system, a start system with solutions.
    The *target* is a list of strings representing the polynomials
    of the target system (which has to be solved).
    The *start* is a list of strings representing the polynomials
    of the start system, with known solutions in sols.
    The *startsols* is a list of strings representing start solutions.
    On return is a tuple, with first the gamma used in the homotopy
    and then second, the string representations of the solutions
    computed at the end of the paths.
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
    phc = get_phcfun()
    apars = (c_int32 * 6)()
    apars[0] = c_int32(0)
    apars[1] = c_int32(0)
    apars[2] = c_int32(vrblvl)
    apars[3] = c_int32(0)
    apars[4] = c_int32(0)
    apars[5] = c_int32(len(target))
    pars = pointer(apars)
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl)
    if vrblvl > 0:
        print('-> double_track calls phc', end='')
    retval = phc(739, pars, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    sols = get_double_solutions(vrblvl)
    return (usedgamma, sols)

def test_tuning(vrblvl=0):
    """
    Tests the tuning of the parameters.
    """
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
    Runs on the mickey mouse example of two quadrics.
    """
    mickey = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    start, startsols = total_degree_start_system(mickey, vrblvl=vrblvl)
    print('the start system :')
    for pol in start:
        print(pol)
    print('the start solutions :')
    for (idx, sol) in enumerate(startsols):
        print('Solution', idx+1, ':')
        print(sol)
    gamma, sols = double_track(mickey, start, startsols, vrblvl=vrblvl)
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

def main():
    """
    Runs some tests on tuning and tracking.
    """
    lvl = 10
    fail = test_tuning(lvl)
    fail = fail + test_double_track(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()
