"""
The module sweepers exports the definition of sweep homotopies and 
the tracking of solution paths defined by sweep homotopies.
A sweep homotopy is a polynomial system where some of the variables
are considered as parameters.  Given solutions for some parameters
and new values for the parameters, we can track the solution paths
starting at the given solutions and ending at the new solutions for
the new values of the parameters.
The sweep is controlled by a convex linear combination between the
list of start and target values for the parameters.
"""

def standard_sweep(pols, sols, pars, start, target):
    """
    For the polynomials in the list of strings pols
    and the solutions in sols for the values in the list start,
    a sweep through the parameter space will be peformed
    in standard double precision to the target values of
    the parameters in the list target.
    The list of symbols in pars contains the names of the variables
    in the polynomials pols that serve as parameters.
    The size of the lists pars, start, and target must be same.
    """
    return sols

def test():
    """
    Runs a sweep on two points on the unit circle.
    """
    from solutions import make_solution as makesol
    circle = ['x^2 + y^2 - 1;']
    first = makesol(['x', 'y'], [0, 1])
    second = makesol(['x', 'y'], [0, -1])
    startsols = [first, second]
    xpar = ['x']
    xstart = [0]
    xtarget = [2]
    newsols = standard_sweep(circle, startsols, xpar, xstart, xtarget)
    for sol in newsols:
        print sol

if __name__ == "__main__":
    test()
