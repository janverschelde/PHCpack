"""
To install phcpy, the file libPHCpack must be placed
to the same location where the modules of phcpy are installed.
"""

from distutils.core import setup

setup(
    name = 'PHCpy' ,
    author = 'Jan Verschelde' ,
    author_email = 'janv@uic.edu' ,
    description = 'a package for Polynomial Homotopy Continuation' ,
    url = 'https://github.com/janverschelde/PHCpack' ,
    version = '1.1.3' ,
    packages = ['phcpy'] ,
    py_modules = ['phcpy/version', 'phcpy/dimension', 'phcpy/polynomials', \
        'phcpy/solutions', 'phcpy/volumes', 'phcpy/solver', \
        'phcpy/examples', 'phcpy/series', 'phcpy/families', \
        'phcpy/homotopies', 'phcpy/trackers', 'phcpy/curves', \
        'phcpy/schubert' ],
    license = 'GNU GENERAL PUBLIC LICENSE version 3' ,
    package_data = {'phcpy':['libPHCpack.so', 'libPHCpack.dll', \
        'libPHCpack.dylib']} ,
    platforms = ['linux', 'windows', 'macosx'] ,
    long_description=open('README.md').read()
)
