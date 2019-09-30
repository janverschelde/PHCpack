"""
The module polynomials exports the definition of the class
Polynomials, as the object oriented interface to PHCpack.
"""

from phcpy import solver
from phcpy.solutions import Solution

class Polynomials(object):
    """
    An instance of this class is a list of polynomials,
    which represents a polynomial system in several variables,
    with coefficients viewed by default as complex numbers.
    """
    def __init__(self, pols):
        """
        A polynomial system is contructed as a list of strings.
        """
        self.pols = pols
        self.vars = solver.names_of_variables(pols)

    def __str__(self):
        """
        Returns the string representation of a polynomial system.
        """
        result = ''
        for pol in self.pols:
            result = result + pol + '\n'
        return result[:-1]

    def __repr__(self):
        """
        Defines the representation as the string representation.
        """
        return str(self)

    def variables(self):
        """
        Returns the list of the variables in the polynomials.
        """
        return self.vars

    def solve(self, nbtasks=0, verbose=True):
        """
        Applies the blackbox solver and returns a list of solutions.
        Multitasking is applied when the number of tasks in nbtasks
        is set to the number of tasks.
        By default, verbose is True, and root counts are written.
        """
        sols = solver.solve(self.pols, verbose=verbose, tasks=nbtasks)
        result = [Solution(sol) for sol in sols]
        return result

def test():
    """
    Tests the methods in the class Polynomials.
    """
    pols = solver.random_trinomials()
    p = Polynomials(pols)
    print(p)
    print('the variables :', p.variables())
    s = p.solve()
    for sol in s:
        print(sol)

if __name__ == "__main__":
    test()
