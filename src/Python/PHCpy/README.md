PHCpy is the latest development of phcpy.

Ctypes allows to import dynamically linked libraries into a Python session.

For the Python scripts to work, the file libPHCpack,
for the proper platform (extension .so, .dylib, or .dll) must be present
in the phcpy folder.

| file name              | description                                   |
|------------------------|-----------------------------------------------|
| version                | gets the PHCpack version string               |
| dimension              | dimension of the polynomial systems           |
| polynomials            | set and get of polynomials                    |
| solutions              | exports operations on solutions               |
| volumes                | mixed volumes and stable mixed volumes        |
| solver                 | exports the blackbox solver                   |
| examples               | runs some examples                            |
| series                 | series expansions of solution curves          |

The directed acyclic graph shows the module dependencies:

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
