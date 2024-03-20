"""
PHCpy --- a package for Polynomial Homotopy Continuation
========================================================

PHCpy is a collection of Python modules to compute solutions
of polynomial systems using PHCpack.

A homotopy defines the deformation of a start system
(system with known solutions) into the target system
(system that has to be solved).
Continuation or path tracking methods apply numerical
predictor-corrector techniques to track the solution paths
defined by the homotopy, starting at the known solutions of the
start system and ending at the solutions of the target system.

Available modules
-----------------
version
    retrieves the PHCpack version string and defines the basic
    interactions with ctypes.
dimension
    exports the dimension of the system of (Laurent) polynomials.
polynomials
    polynomials are represented as strings with coefficients
    in double, double double, or quad double precision;
    Laurent polynomials can have negative exponents.
solutions
    solutions are lists of PHCpack solution strings
    and this module exports operations to convert the solution
    strings into Python dictionaries, e.g. for evaluation.
solver
    exports the blackbox solver of PHCpack.
volumes
    exports functions to compute mixed volumes and stable mixed volumes,
    which provides generically sharp root counts.
examples
    contains some interesting examples from the research literature.
families
    polynomial system often occur in families and are defined for any
    number of equations and variables, e.g.: the cyclic n-roots system.
scaling
    equation and variable scaling of polynomials.
starters
    constructors of start systems for artificial parameter homotopies.
homotopies
    a homotopy connects the target system to a start system.
trackers
    path trackers with aposteriori step size control.
tropisms
    polyhedral end game with aposteriori step size control.
sweepers
    arc length parameter continuation in sweep homotopies.
series
    allows to compute series developments of solution curves
    defined by polynomial homotopies.
curves
    exports apriori step size control path trackers.
deflation
    restores quadratic convergence at an isolated singular solution.
schubert
    numerical Schubert calculus defines homotopies for enumerative geometry.
sets
    positive dimensional solution sets are represented by system with
    slack variables and a set of generic points.
cascades
    homotopies to compute candidate generic points on all components.
diagonal
    intersecting positive dimensional solution sets by diagonal homotopies.
factor
    factor positive dimensional solution sets into irreducible components.
decomposition
    solving is computing a numerical irreducible decomposition.
binomials
    a binomial system has exactly two terms in every equation
    and can be solved faster.

The main() of every module provides some basic examples and tests.
"""

def cite():
    """
    Displays the citation information for phcpy.
    """
    print("""
    To cite phcpy in publications use:

    Jan Verschelde: Modernizing PHCpack through phcpy.
    In Proceedings of the 6th European Conference on Python in Science
    (EuroSciPy 2013), edited by Pierre de Buyl and Nelle Varoquaux,
    pages 71-76, 2014, available at http://arxiv.org/abs/1310.0056.

    Jasmine Otto, Angus Forbes, and Jan Verschelde:
    Solving Polynomial Systems with phcpy.  In the Proceedings of the
    18th Python in Science Conference (SciPy 2019), edited by Chris Calloway,
    David Lippa, Dillon Niederhut, and David Shupe, pages 58-64, 2019.
    """)

from os import getcwd, chdir
from site import getsitepackages

def get_site_location(vrblvl=0):
    """
    Returns the location for site packages,
    handling the differences between linux and windows.
    """
    if vrblvl > 0:
        print('in get_site_location ...')
    sites = getsitepackages()
    if vrblvl > 0:
        print('sites :', sites)
    bools = ['site-packages' in x for x in sites]
    if len(bools) > 0:
        if True in bools:
            idx = bools.index(True)
            return sites[idx]
    if vrblvl > 0:
        print('no site-packages, look for dist-packages ...')
    bools = ['dist-packages' in x for x in sites]
    idx = bools.index(True)
    return sites[idx]

def set_phcfun(vrblvl=0):
    """
    Sets the variable phc to the function in libPHCpack
    and loads the modules.
    """
    if vrblvl > 0:
        print('in set_phcfun ...')
    cwd = getcwd()
    if vrblvl > 0:
        print('cwd =', cwd)
    location = get_site_location(vrblvl) + '/phcpy'
    if vrblvl > 0:
        print('location :')
        print(location)
    chdir(location)
    if vrblvl > 0:
        print('os.getcwd :')
        print(getcwd())
    try:
        from phcpy.version import version_string
        print(version_string(vrblvl=0) + ' works!')
        from phcpy.version import get_phcfun_fromlib
        result = get_phcfun_fromlib(vrblvl=0)
        chdir(cwd)
        return result
    except:
        print('Is the libPHCpack not suited for this platform?')
        return None

phc = set_phcfun()
from phcpy import version, dimension, polynomials, solutions
from phcpy import solver, volumes, examples, families
from phcpy import scaling, starters, homotopies, trackers
from phcpy import tropisms, sweepers
from phcpy import series, curves, deflation, schubert
from phcpy import sets, cascades, diagonal, factor
from phcpy import decomposition, binomials

# The version number is defined as a data attribute.
__version__ = '1.1.4'
