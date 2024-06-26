"""
To install phcpy correctly, the shared object file must be copied
to the same location where the modules of phcpy are installed.
This is accomplished via the data_files entry in the setup below.
For the installation of the GPU accelerated path trackers in double,
double double and quad double precision, instead of phcpy2c3.so,
the shared objects that are needed are phcpy2c3path_d.so (double),
phcpy2c3path_dd.so (double double), and phcpy2c3path_qd.so (quad double).
"""

from distutils.core import setup

setup(
    name = 'PHCpy' ,
    author = 'Jan Verschelde' ,
    author_email = 'janv@uic.edu' ,
    description = 'a package for Polynomial Homotopy Continuation' ,
    url = 'https://github.com/janverschelde/PHCpack' ,
    version = '1.1.2' ,
    packages = ['phcpy'] ,
    py_modules = ['phcpy/interface', 'phcpy/solver', \
        'phcpy/solutions', 'phcpy/polynomials', \
        'phcpy/trackers', 'phcpy/sweepers', 'phcpy/tuning', \
        'phcpy/sets', 'phcpy/cascades', 'phcpy/factor', 'phcpy/diagonal', \
        'phcpy/polytopes', 'phcpy/tropisms', 'phcpy/maps', 'phcpy/series', \
        'phcpy/curves', 'phcpy/examples', 'phcpy/families', \
        'phcpy/schubert' , 'phcpy/dashboard', 'phcpy/server' ],
    license = 'GNU GENERAL PUBLIC LICENSE version 3' ,
    package_data = {'phcpy':['phcpy2c3.so']} ,
    platforms = ['linux2'] ,
    long_description=open('README.txt').read()
)
