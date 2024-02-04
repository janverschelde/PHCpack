PHCpy is the latest development of phcpy.

Ctypes allows to import dynamically linked libraries into a Python session.

For the Python scripts to work, the file libPHCpack,
for the proper platform (extension .so, .dylib, or .dll) must be present
in the phcpy folder.

| module name    | description of each module                        |
|----------------|---------------------------------------------------|
| version        | gets the PHCpack version string                   |
| dimension      | dimension of the polynomial systems               |
| polynomials    | set and get of polynomials                        |
| solutions      | exports operations on solutions                   |
| volumes        | mixed volumes and stable mixed volumes            |
| solver         | exports the blackbox solver                       |
| examples       | runs some examples                                |
| families       | families of systems, for any dimension            |
| scaling        | equation and variable scaling of polynomials      |
| starters       | constructing start systems for homotopies         |
| homotopies     | homotopies are systems in one parameter           |
| trackers       | aposteriori step size control path trackers       |
| tropisms       | polyhedral end game with aposteriori step control |
| sweepers       | arc length parameter continuation for real        |
| series         | series expansions of solution curves              |
| curves         | apriori step size control path trackers           |
| deflation      | quadratic convergence at isolated singularity     |
| schubert       | homotopies for enumerative geometry               |
| sets           | representing positive dimensional solution sets   |
| cascades       | generic points on all solution components         |
| diagonal       | intersecting positive dimensional solutions sets  |
| factor         | factor sets in irreducible components             |
| decomposition  | a numerical irreducible decomposition             |
| binomials      | solving systems of two terms in every equation    |

The directed acyclic graph shows the dependencies of
the first seven modules:

    version
       |
       +----------------------------> solutions         
       |                                  |
       +--> dimension                     |
       |        |                         |
       +--------+--> polynomials          |
       |        |         |               |
       +------------------+--> volumes    |
       |        |         |       |       |
       +--------+---------+-------+-------+--> solver 
       |        |                                |
       +--------+--------------------------------+--> examples

At the root is the version module.  If version.py works,
then the interfacing with libPHCpack works.
Four other modules are needed for the blackbox solver.

To install, type
pip install .
