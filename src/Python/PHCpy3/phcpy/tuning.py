"""
The module tuning provides functions to tune the tolerances and
settings of the predictor and corrector parameters for the path trackers.
"""

def tune_track_parameters(difficulty=0, digits=16, \
                          interactive=False, silent=True):
    """
    Tunes the numerical continuation parameters for use
    during path tracking based on two parameters:
    the difficulty of the solution path (difficulty)
    and the number of decimal places (digits) for the
    accuracy of the approximations along a path.
    Increasing difficulty will decrease the step size
    and increase the number of steps along a path.
    Increasing digits will decrease the tolerances on
    the numerical approximations.
    If interactive is True, then the user can interactively
    adjust specific numerical continuation parameters.
    If silent is False, then the current values of the
    numerical continuation parameters are shown.
    """
    from phcpy.phcpy2c3 import py2c_autotune_continuation_parameters
    from phcpy.phcpy2c3 import py2c_tune_continuation_parameters
    from phcpy.phcpy2c3 import py2c_show_continuation_parameters
    py2c_autotune_continuation_parameters(difficulty, digits)
    if interactive:
        py2c_tune_continuation_parameters()
    elif not silent:
        py2c_show_continuation_parameters()

def condition_level_get():
    """
    Returns an integer that represents the difficulty level
    of the homotopy.  The default level equals zero,
    higher values lead to smaller tolerances.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return int(get(1))

def condition_level_set(lvl):
    """
    Sets the parameter that represents the difficulty level of the
    homotopy to the value of lvl.  The default level equals zero,
    higher values lead to smaller tolerances.
    On return is the failure code, which is zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(1, lvl)

def track_simultaneously_get():
    """
    Returns the number of paths that are tracked simultaneously, for the
    same discretization of the interval of the continuation parameter.
    The default value equals one.  Increasing this value avoids path crossing,
    also called path jumping.  This jumping happens when the approximated
    points on one path transition to approximated points on another path.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return int(get(2))

def track_simultaneously_set(nbr):
    """
    Sets the number of paths that are tracked simultaneously, for the
    same discretization of the interval of the continuation parameter,
    to the value of nbr.
    The default value equals one.  Increasing this value avoids path crossing,
    also called path jumping.  This jumping happens when the approximated
    points on one path transition to approximated points on another path.
    On return in the failure code, which is zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(2, nbr)

def max_steps_get():
    """
    Returns the maximum number of steps the path tracker will perform
    to reach the end of a solution path.  For paths that diverge to
    infinity are often truncated beore reaching extreme values.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return int(get(3))

def max_steps_set(mxs):
    """
    Sets the maximum number of steps the path tracker will perform
    to reach the end of a solution path to mxs.  For paths that diverge to
    infinity are often truncated beore reaching extreme values.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(3, mxs)

def distance_to_endgame_get():
    """
    Returns the distance to start the end game.  During the end game,
    the path tracker may apply tolerances that are more severe as the
    solution paths get more interesting near the end.
    The default value is 0.1.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return get(4)

def distance_to_endgame_set(dst):
    """
    Sets the distance to start the end game to the value of dst.
    During the end game, the path tracker may apply tolerances that are 
    more severe as the solution paths get more interesting near the end.
    The default value is 0.1.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(4, dst)

def order_endgame_extrapolator_get():
    """
    Returns the order of the extrapolator to estimate the winding number
    in a polyhedral end game. 
    If the order is zero, then no polyhedral end game will be applied.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return int(get(5))

def order_endgame_extrapolator_set(ord):
    """
    Sets the order of the extrapolator to estimate the winding number
    in a polyhedral end game to the value of ord.
    If the order is zero, then no polyhedral end game will be applied.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(5, ord)

def maxnum_reruns_get():
    """
    Returns the value of the maximum number of path reruns.
    If path jumping is detected, then the clustered paths are retracked
    with more severe values of the tolerances.
    The default value equals one.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return int(get(6))

def maxnum_reruns_set(mrr):
    """
    Sets the value of the maximum number of path reruns to the value of mrr.
    If path jumping is detected, then the clustered paths are retracked
    with more severe values of the tolerances.
    The default value equals one.
    On return is the failure code, which is zero when all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(6, mrr)

def predictor_typeonpath_get():
    """
    Returns an integer which represents the type of the predictor
    along a path, before the start of the end game.
    The integer on return takes values between 0 and 9, depending on
    the type for the solution x and for the continuation parameter t.
    The ten predictor types are
    0 : secant for x, real for t;
    1 : secant for x, complex for t;
    2 : secant for x, geometric for t;
    3 : tangent for x, real for t;
    4 : tangent for x, complex for t;
    5 : tangent for x, geometric for t;
    6 : Hermite for x, real for t;
    7 : quadratic for x, real for t;
    8 : cubic for x, real for t;
    9 : quartic for x, real for t.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return int(get(7))

def predictor_typeonpath_set(predtype):
    """
    Sets the type of the predictor along a path, before the end game,
    to what the value of predtype represents.
    A valid integer value for predtype lies between 0 and 9, setting
    the type for the solution x and for the continuation parameter t.
    The ten predictor types are
    0 : secant for x, real for t;
    1 : secant for x, complex for t;
    2 : secant for x, geometric for t;
    3 : tangent for x, real for t;
    4 : tangent for x, complex for t;
    5 : tangent for x, geometric for t;
    6 : Hermite for x, real for t;
    7 : quadratic for x, real for t;
    8 : cubic for x, real for t;
    9 : quartic for x, real for t.
    On return is the failure code, which is zero when all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(7, predtype)

def predictor_typeendgame_get():
    """
    Returns an integer which represents the type of the predictor
    along a path, before the start of the end game.
    The integer on return takes values between 0 and 9, depending on
    the type for the solution x and for the continuation parameter t.
    The ten predictor types are
    0 : secant for x, real for t;
    1 : secant for x, complex for t;
    2 : secant for x, geometric for t;
    3 : tangent for x, real for t;
    4 : tangent for x, complex for t;
    5 : tangent for x, geometric for t;
    6 : Hermite for x, real for t;
    7 : quadratic for x, real for t;
    8 : cubic for x, real for t;
    9 : quartic for x, real for t.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return int(get(8))

def predictor_typeendgame_set(predtype):
    """
    Sets the type of the predictor during the end game to what the
    value of predtype represents.
    A valid integer value for predtype lies between 0 and 9, setting
    the type for the solution x and for the continuation parameter t.
    The ten predictor types are
    0 : secant for x, real for t;
    1 : secant for x, complex for t;
    2 : secant for x, geometric for t;
    3 : tangent for x, real for t;
    4 : tangent for x, complex for t;
    5 : tangent for x, geometric for t;
    6 : Hermite for x, real for t;
    7 : quadratic for x, real for t;
    8 : cubic for x, real for t;
    9 : quartic for x, real for t.
    On return is the failure code, which is zero when all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(8, predtype)

def predictor_minsteponpath_get():
    """
    Returns the minimum value of the step size along a path, before the
    end game.  If the step size control cuts the step size to a value
    below this minimum, then the path tracking is aborted.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return get(9)

def predictor_minsteponpath_set(minstep):
    """
    Sets the minimum of the step size along a path before the end game
    to the value of minstep.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(9, minstep)

def predictor_minstependgame_get():
    """
    Returns the minimum value of the step size during the end game.
    If the step size control cuts the step size to a value
    below this minimum, then the path tracking is aborted.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return get(10)

def predictor_minstependgame_set(minstep):
    """
    Sets the minimum of the step size during the end game
    to the value of minstep.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(10, minstep)

def predictor_maxsteponpath_get():
    """
    Returns the maximum value of the step size along a path, before the
    end game.  The step size control will never increase the step size
    to a value above this maximum.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return get(11)

def predictor_maxsteponpath_set(maxstep):
    """
    Sets the maximum of the step size along a path before the end game
    to the value of maxstep.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(11, minstep)

def predictor_maxstependgame_get():
    """
    Returns the maximum value of the step size during the end game.
    The step size control will never increase the step size to a value
    above this maximum.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return get(12)

def predictor_maxtependgame_set(maxstep):
    """
    Sets the maximum of the step size during the end game
    to the value of maxstep.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(12, maxstep)

def predictor_redfaconpath_get():
    """
    Returns the value of the reduction factor to cut the step size in case
    of a failed corrector stage, along a path, before the end game.
    The reduction factor determines the speed at which the predictor reduces
    its step size when tracking a more difficult portion of the path.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return get(13)

def predictor_redfaconpath_set(redfac):
    """
    Sets the value of the reduction factor to cut the step size in case of
    a failed corrector step, along a path, before the end game, to the
    value of redfac.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(13, redfac)

def predictor_redfacendgame_get():
    """
    Returns the value of the reduction factor to cut the step size in case
    of a failed corrector stage, during the end game.
    The reduction factor determines the speed at which the predictor reduces
    its step size when tracking a more difficult portion of the path.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return get(14)

def predictor_redfacendgame_set(redfac):
    """
    Sets the value of the reduction factor to cut the step size in case of
    a failed corrector step during the end game, to the value of redfac.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(14, redfac)

def predictor_expfaconpath_get():
    """
    Returns the value of the expansion factor to increase the step size in
    case of a successful corrector stage, along a path, before the end game.
    The expansion factor determines the speed at which the predictor increases
    its step size when tracking an easier portion of the path.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return get(15)

def predictor_expfaconpath_set(expfac):
    """
    Sets the value of the expansion factor to increase the step size in
    case of a successful corrector stag, along a path, before the end game,
    to the value of expfac.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(15, expfac)

def predictor_expfacendgame_get():
    """
    Returns the value of the expansion factor to increase the step size in
    case of a successful corrector stage, during the end game.
    The expansion factor determines the speed at which the predictor increases
    its step size when tracking an easier portion of the path.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return get(16)

def predictor_expfacendgame_set(expfac):
    """
    Sets the value of the expansion factor to increase the step size in
    case of a successful corrector stag, during the end game,
    to the value of expfac.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(16, expfac)

def predictor_exptrsonpath_get():
    """
    Returns the value of the expansion threshold for the step size control,
    along a path, before the end game.
    The expansion threshold is the number of consecutive successful
    corrector stages that must be met before the step size is increased.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return int(get(17))

def predictor_exptrsonpath_set(exptrs):
    """
    Sets the value of the expansion threshold for the step size control,
    along a path, before the end game, to the value of exptrs.
    The expansion threshold is the number of consecutive successful
    corrector stages that must be met before the step size is increased.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(17, exptrs)

def predictor_exptrsendgame_get():
    """
    Returns the value of the expansion threshold for the step size control,
    during the end game.
    The expansion threshold is the number of consecutive successful
    corrector stages that must be met before the step size is increased.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return int(get(18))

def predictor_exptrsendgame_set(exptrs):
    """
    Sets the value of the expansion threshold for the step size control,
    during the end game, to the value of exptrs.
    The expansion threshold is the number of consecutive successful
    corrector stages that must be met before the step size is increased.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(18, exptrs)

def corrector_maxiteronpath_get():
    """
    Returns the maximum number of iterations the corrector does along
    a path, before the start of the end game.  The default equals 3.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return int(get(19))

def corrector_maxiteronpath_set(maxiter):
    """
    Sets the maximum number of iterations the corrector does along
    a path, before the start of the end game, to the value of maxiter.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(19, maxiiter)

def corrector_maxiterendgame_get():
    """
    Returns the maximum number of iterations the corrector does along
    a path, during the end game.  The default equals 3.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return int(get(20))

def corrector_maxiterendgame_set(maxiter):
    """
    Sets the maximum number of iterations the corrector does along
    a path, before the start of the end game, to the value of maxiter.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(20, maxiiter)

def corrector_relresonpath_get():
    """
    Returns the value of the tolerance on the relative residual for
    the corrector along a path, before the end game.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return int(get(21))

def corrector_relresonpath_set(tol):
    """
    Sets the tolerance on the relative residal for the corrector along
    a path, before the start of the end game, to the value of tol.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(21, tol)

def corrector_relresendgame_get():
    """
    Returns the value of the tolerance on the relative residual for
    the corrector during the end game.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return int(get(22))

def corrector_relresendgame_set(tol):
    """
    Sets the tolerance on the relative residal for the corrector during
    the end game, to the value of tol.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(22, tol)

def corrector_absresonpath_get():
    """
    Returns the value of the tolerance on the absolute residual for
    the corrector along a path, before the end game.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return int(get(23))

def corrector_absresonpath_set(tol):
    """
    Sets the tolerance on the absolute residal for the corrector along
    a path, before the start of the end game, to the value of tol.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(23, tol)

def corrector_absresendgame_get():
    """
    Returns the value of the tolerance on the absolute residual for
    the corrector during the end game.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return int(get(24))

def corrector_absresendgame_set(tol):
    """
    Sets the tolerance on the absolute residal for the corrector during
    the end game, to the value of tol.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(24, tol)

def corrector_relcoronpath_get():
    """
    Returns the value of the tolerance on the relative correction for
    the corrector along a path, before the end game.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return int(get(25))

def corrector_relcoronpath_set(tol):
    """
    Sets the tolerance on the relative correction for the corrector along
    a path, before the start of the end game, to the value of tol.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(25, tol)

def corrector_relcorendgame_get():
    """
    Returns the value of the tolerance on the relative correction for
    the corrector during the end game.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return int(get(26))

def corrector_relcorendgame_set(tol):
    """
    Sets the tolerance on the relative correction for the corrector during
    the end game, to the value of tol.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(26, tol)

def corrector_abscoronpath_get():
    """
    Returns the value of the tolerance on the absolute correction for
    the corrector along a path, before the end game.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return int(get(27))

def corrector_abscoronpath_set(tol):
    """
    Sets the tolerance on the absolute correction for the corrector along
    a path, before the start of the end game, to the value of tol.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(27, tol)

def corrector_abscorendgame_get():
    """
    Returns the value of the tolerance on the absolute correction for
    the corrector during the end game.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return int(get(28))

def corrector_absresendgame_set(tol):
    """
    Sets the tolerance on the absolute residal for the corrector during
    the end game, to the value of tol.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(28, tol)

def tolerance_rcondonpath_get():
    """
    Returns the tolerance on the inverse condition number of the Jacobian
    matrix of a solution along a path, before the end game, to decide
    whether a solution is singular.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return get(29)

def tolerance_rcondonpath_set(tol):
    """
    Sets the tolerance on the inverse condition number of the Jacobian
    matrix of a solution along a path, before the end game, to decide
    whether a solution is singular, to the value of tol.
    The failure code is returned, which is zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(29, tol)

def tolerance_rcondendgame_get():
    """
    Returns the tolerance on the inverse condition number of the Jacobian
    matrix of a solution, during the end game, to decide whether a solution
    is singular or not.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return get(30)

def tolerance_rcondendgame_set(tol):
    """
    Sets the tolerance on the inverse condition number of the Jacobian
    matrix of a solution along a path, before the end game, to decide
    whether a solution is singular, to the value of tol.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(30, tol)

def tolerance_clustsolonpath_get():
    """
    Returns the tolerance on the distance between two solutions to decide
    whether two solutions are clustered along a path, before the end game.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return get(31)

def tolerance_clustsolonpath_set(tol):
    """
    Sets the tolerance on the distance between two solutions to decide
    whether two solutions are clustered along a path, before the end game,
    to the value of tol.
    On return is the failure code which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(31, tol)

def tolerance_clustsolendgame_get():
    """
    Returns the tolerance on the distance between two solutions to decide
    whether two solutions are clustered during the end game.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return get(32)

def tolerance_clustsolendgame_set(tol):
    """
    Sets the tolerance on the distance between two solutions to decide
    whether two solutions are clustered during the end game,
    to the value of tol.
    On return is the failure code which equals zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(32, tol)

def tolerance_infsolonpath_get():
    """
    Returns the tolerance threshold to decide whether a solution path
    diverges to infinity, before the start of the end game.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return get(33)

def tolerance_infsolonpath_set(tol):
    """
    Sets the tolerance threshold to decice whether a solution path
    diverges to infinity, during the end game, to the value of tol.
    On return is the failure code which is zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(33, tol)

def tolerance_infsolendgame_get():
    """
    Returns the tolerance threshold to decide whether a solution path
    diverges to infinity, during the end game.
    """
    from phcpy.phcpy2c3 import py2c_get_value_of_continuation_parameter as get
    return get(34)

def tolerance_infsolendgame_set(tol):
    """
    Sets the tolerance threshold to decice whether a solution path
    diverges to infinity, during the end game, to the value of tol.
    On return is the failure code which is zero if all went well.
    """
    from phcpy.phcpy2c3 import py2c_set_value_of_continuation_parameter as set
    return set(34, tol)

def test():
    """
    Tests the tuning of the parameters.
    """
    print('\ntesting tuning of parameters...\n')
    tune_track_parameters(silent=False)

if __name__ == "__main__":
    test()
