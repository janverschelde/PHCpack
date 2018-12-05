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

def test():
    """
    Returns the current values of the homotopy continuation parameters.
    """
    pars = [get_homotopy_continuation_parameter(k) for k in range(1, 13)]
    print(pars)
    tune_homotopy_continuation_parameters()

if __name__ == "__main__":
    test()
