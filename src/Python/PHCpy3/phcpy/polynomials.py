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
        self.startpols = []
        self.startsols = []
        self.gamma = complex(0)
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

    def solve(self, verbose=True, nbtasks=0, mvfocus=0, \
        dictionary_output=False, verbose_level=0):
        r"""
        Applies the blackbox solver and returns a list of solutions.
        By default, *verbose* is True, and root counts are written.
        Multitasking is applied when the number of tasks in *nbtasks*
        is set to the number of tasks.
        To apply only mixed volumes and polyhedral homotopies,
        set the value for *mvfocus* to 1.
        If *dictionary_output*, then on return is a list of dictionaries,
        else the returned list is a list of strings.
        If *verbose_level* is larger than 0, then the names of the procedures
        called in the running of the blackbox solver will be listed.
        The solving happens in standard double precision arithmetic.
        """
        from phcpy.phcpy2c3 \
            import py2c_copy_standard_Laurent_start_system_to_container
        from phcpy.phcpy2c3 import py2c_copy_start_solutions_to_container
        from phcpy.phcpy2c3 import py2c_solcon_clear_standard_solutions
        from phcpy.phcpy2c3 import py2c_get_gamma_constant
        from phcpy.interface import load_standard_laurent_system as qload
        from phcpy.interface import load_standard_solutions as qsolsload
        sols = solver.solve(self.pols, verbose=verbose, tasks=nbtasks, \
            mvfocus=mvfocus, dictionary_output=dictionary_output, \
            verbose_level=verbose_level)
        py2c_copy_standard_Laurent_start_system_to_container()
        (regamma, imgamma) = py2c_get_gamma_constant(1, verbose_level)
        self.gamma = complex(regamma, imgamma)
        self.startpols = qload()
        py2c_solcon_clear_standard_solutions()
        py2c_copy_start_solutions_to_container()
        self.startsols = [Solution(sol) for sol in qsolsload()]
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
    q = p.startpols
    print('the start system :')
    for pol in q:
        print(pol)
    qsols = p.startsols
    print('the start solutions :')
    for sol in qsols:
        print(sol)

if __name__ == "__main__":
    test()
