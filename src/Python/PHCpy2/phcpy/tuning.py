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
    from phcpy.phcpy2c import py2c_autotune_continuation_parameters
    from phcpy.phcpy2c import py2c_tune_continuation_parameters
    from phcpy.phcpy2c import py2c_show_continuation_parameters
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
    from phcpy.phcpy2c import py2c_get_value_of_continuation_parameter as get
    return int(get(1))

def condition_level_set(lvl):
    """
    Sets the parameter that represents the difficulty level of the
    homotopy to the value of lvl.  The default level equals zero,
    higher values lead to smaller tolerances.
    On return is the failure code, which is zero if all went well.
    """
    from phcpy.phcpy2c import py2c_set_value_of_continuation_parameter as set
    return set(1, lvl)

def track_simultaneously_get():
    """
    Returns the number of paths that are tracked simultaneously, for the
    same discretization of the interval of the continuation parameter.
    The default value equals one.  Increasing this value avoids path crossing,
    also called path jumping.  This jumping happens when the approximated
    points on one path transition to approximated points on another path.
    """
    from phcpy.phcpy2c import py2c_get_value_of_continuation_parameter as get
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
    from phcpy.phcpy2c import py2c_set_value_of_continuation_parameter as set
    return set(2, nbr)

def max_steps_get():
    """
    Returns the maximum number of steps the path tracker will perform
    to reach the end of a solution path.  For paths that diverge to
    infinity are often truncated beore reaching extreme values.
    """
    from phcpy.phcpy2c import py2c_get_value_of_continuation_parameter as get
    return int(get(3))

def max_steps_set(mxs):
    """
    Sets the maximum number of steps the path tracker will perform
    to reach the end of a solution path to mxs.  For paths that diverge to
    infinity are often truncated beore reaching extreme values.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c import py2c_set_value_of_continuation_parameter as set
    return set(3, mxs)

def distance_to_endgame_get():
    """
    Returns the distance to start the end game.  During the end game,
    the path tracker may apply tolerances that are more severe as the
    solution paths get more interesting near the end.
    The default value is 0.1.
    """
    from phcpy.phcpy2c import py2c_get_value_of_continuation_parameter as get
    return get(4)

def distance_to_endgame_set(dst):
    """
    Sets the distance to start the end game to the value of dst.
    During the end game, the path tracker may apply tolerances that are 
    more severe as the solution paths get more interesting near the end.
    The default value is 0.1.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c import py2c_set_value_of_continuation_parameter as set
    return set(4, dst)

def order_endgame_extrapolator_get():
    """
    Returns the order of the extrapolator to estimate the winding number
    in a polyhedral end game. 
    If the order is zero, then no polyhedral end game will be applied.
    """
    from phcpy.phcpy2c import py2c_get_value_of_continuation_parameter as get
    return int(get(5))

def order_endgame_extrapolator_set(ord):
    """
    Sets the order of the extrapolator to estimate the winding number
    in a polyhedral end game to the value of ord.
    If the order is zero, then no polyhedral end game will be applied.
    On return is the failure code, which equals zero if all went well.
    """
    from phcpy.phcpy2c import py2c_set_value_of_continuation_parameter as set
    return set(5, ord)

def maxnum_reruns_get():
    """
    Returns the value of the maximum number of path reruns.
    If path jumping is detected, then the clustered paths are retracked
    with more severe values of the tolerances.
    The default value equals one.
    """
    from phcpy.phcpy2c import py2c_get_value_of_continuation_parameter as get
    return int(get(6))

def maxnum_reruns_set(mrr):
    """
    Sets the value of the maximum number of path reruns to the value of mrr.
    If path jumping is detected, then the clustered paths are retracked
    with more severe values of the tolerances.
    The default value equals one.
    On return is the failure code, which is zero when all went well.
    """
    from phcpy.phcpy2c import py2c_set_value_of_continuation_parameter as set
    return set(6, mrr)

def test():
    """
    Runs test_track(), test_next_track(), and test_monitored_track().
    """
    print('\ntesting tuning of parameters...\n')
    tune_track_parameters(silent=False)

if __name__ == "__main__":
    test()
