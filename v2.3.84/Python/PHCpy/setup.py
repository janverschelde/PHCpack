"""
To install phcpy correctly, the shared object file must be copied
to the same location where the modules of phcpy are installed.
This is accomplished via the data_files entry in the setup below.
"""

from distutils.core import setup
from distutils.sysconfig import get_python_lib

setup(
   name = 'PHCpy' ,
   author = 'Jan Verschelde' ,
   author_email = 'jan@math.uic.edu' ,
   description = 'a package for Polynomial Homotopy Continuation' ,
   url = 'http://www.math.uic.edu/~jan/download.html' ,
   version = '0.1.4' ,
   packages = ['phcpy'] ,
   py_modules = ['phcpy/solver', 'phcpy/phcsols', 'phcpy/phcsets',
                 'phcpy/examples', 'phcpy/families', 'phcpy/schubert' ] ,
   license = 'GNU GENERAL PUBLIC LICENSE version 3' ,
   data_files = [(get_python_lib()+'/phcpy', ['phcpy/phcpy2c.so'])] ,
   platforms = ['linux2'] ,
   long_description=open('README.txt').read()
)
