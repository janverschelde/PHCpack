"""
The module tuning provides functions to tune the tolerances and
settings of the predictor and corrector parameters for the path trackers
which apply an aposteriori step control algorithm.
"""
from ctypes import c_int, c_double, pointer
from phcpy.version import get_phcfun

def show_parameters(vrblvl=0):
    """
    Displays the current values of the continuation parameters.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in show_parameters ...')
    phc = get_phcfun()
    aaa = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    if vrblvl > 0:
        print('-> show_parameters calls phc', end='')
    retval = phc(194, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def autotune_parameters(difficulty_level, digits_of_precision, vrblvl=0):
    """
    Tunes the parameters given the difficulty level of the homotopy
    and the digits of working precision.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in autotune_parameters, difficulty :', difficulty_level, end='')
        print(', digits :', digits_of_precision)
    phc = get_phcfun()
    aaa = pointer(c_int(difficulty_level))
    bbb = pointer(c_int(digits_of_precision))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    if vrblvl > 0:
        print('-> set_parameter_value calls phc', end='')
    retval = phc(193, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def interactive_tune(vrblvl=0):
    """
    Interactive tuning of the parameters.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in interactive_tune ...')
    phc = get_phcfun()
    aaa = pointer(c_int(0))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    if vrblvl > 0:
        print('-> set_parameter_value calls phc', end='')
    retval = phc(70, aaa, bbb, ccc, vrb)
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
    aaa = pointer(c_int(idx))
    bbb = pointer(c_int(0))
    ccc = pointer(c_double(value))
    vrb = c_int(vrblvl)
    if vrblvl > 0:
        print('-> set_parameter_value calls phc', end='')
    retval = phc(190, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def get_parameter_value(idx, vrblvl=0):
    """
    Returns the value of the parameter with index idx.
    The verbose level is given by vrblvl.
    """
    if vrblvl > 0:
        print('in get_parameter_value, idx :', idx)
    phc = get_phcfun()
    aaa = pointer(c_int(idx))
    bbb = pointer(c_int(0))
    value = pointer(c_double(0.0))
    vrb = c_int(vrblvl)
    if vrblvl > 0:
        print('-> set_parameter_value calls phc', end='')
    retval = phc(189, aaa, bbb, value, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
        print('Value of parameter', idx, ':', value[0])
    return value[0]

def set_condition_level(level, vrblvl=0):
    r"""
    Sets the parameter that represents the difficulty level of the
    homotopy to the value of *level*.  The default level equals zero,
    higher values lead to smaller tolerances.
    On return is the failure code, which is zero if all went well.
    """
    return set_parameter_value(1, level, vrblvl)

def get_condition_level(vrblvl=0):
    """
    Returns the level of difficulty.
    The verbose level is given by vrblvl.
    """
    return int(get_parameter_value(1, vrblvl))

def main():
    """
    Runs some tests.
    """
    lvl = 10
    show_parameters(lvl)
    print('setting the condition level to 2 ...')
    set_condition_level(2, lvl)
    level = get_condition_level(lvl)
    print('the condition level :', level)
    autotune_parameters(level, 14, lvl)
    show_parameters(lvl)
    interactive_tune(lvl)
    show_parameters(lvl)

if __name__== '__main__':
    main()
