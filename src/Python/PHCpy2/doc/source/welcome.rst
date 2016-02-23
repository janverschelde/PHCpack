welcome to phcpy's documentation!
=================================

This documentation describes a collection of Python modules
to compute solutions of polynomial systems using PHCpack.

The computation of the mixed volume in phcpy calls MixedVol
(ACM TOMS Algorithm 846 of T. Gao, T.Y. Li, M. Wu) 
as it is integrated in PHCpack.
For double double and quad double arithmetic, PHCpack incorporates
the QD library of Y. Hida, X.S. Li, and D.H. Bailey.
See the References section for pointers to the literature.

While PHCpack has been under development for over twenty years,
phcpy is still working through its proof-of-concept stage.
In its present state, working with phcpy will require persistence
and plenty of patience.

what is phcpy?
--------------

The main executable phc (polynomial homotopy continuation)
defined by the source code in PHCpack is a menu driven
and file oriented program.
The Python interface defined by phcpy replaces the files
with persistent objects allowing the user to work with
scripts or in interactive sessions.
The computationally intensive tasks such as path tracking
and mixed volume computations are executed as compiled code
so there will not be a loss of efficiency.

Both phcpy and PHCpack are free and open source software packages;
you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software Foundation;
either version 2 of the License, or (at your option) any later version.  

One application of phcpy is to run regression tests.
The Python interface phcpy to PHCpack is a programmer's interface.
The long-term goal is to make PHCpack more versatile,
at least for those programmers familiar 
with the Python scripting language.

installing phcpy
----------------

The source for PHCpack can be downloaded from
<http://www.math.uic.edu/~jan/download.html>
and is also available at
<https://github.com/janverschelde/PHCpack>.
For the installation from source, the gnu-ada compiler 
(available for free at <http://libre.adacore.com/>) is needed.  
Also your Python interpreter must most likely have been built with gcc.
In addition to the Python interpreter, the file Python.h of
of the developer's libraries will need to exist on your system.
Otherwise, PHCpack and phcpy are self contained
and do not require the installation of other software.

Up to version 0.3.7, phcpy was written with versions 2.6 and 2.7 of Python.
Version 0.3.8 of phcpy was ported to Python 3.5,
using a modified C interface phcpy2c3.c and the corresponding
shared object phcpy2c3.so.

The code runs on Red Hat Linux 64-bit, Ubuntu Linux, and on Mac OS X.
There is no support for Windows.
Below is a step-by-step installation procedure.

0. The gnu-ada compiler must be installed to compile the shared object file.
   Although several Linux distributions come with gcc that have Ada enabled,
   check whether ``gnatmake`` is in your execution path.
   In a terminal window, at the command prompt, type ``which gnatmake``.
   If the system answers with the location of ``gnatmake``,
   then the gnu-ada compiler is installed on your system.

   If you have multiple versions of gcc installed on your system,
   then the binaries of the gnu-ada compiler should appear first
   in your execution path.
   Typing ``gcc -v`` should show ``for GNAT GPL`` in the reply.

   If both ``which gnatmake`` and ``gcc -v`` gave satisfactory replies,
   then you can proceed to step 2 and skip the installation of the
   gnu-ada compiler. 

1. By default one needs to have superuser privileges to install
   the gnu-ada compiler at standard locations, but otherwise it is not hard
   to install the gnu-ada compiler in your own directory.

   Following the instructions of the gnu-ada compiler, the location
   with the binaries must be added in front of the execution path.
   You may have to edit ``.bashrc`` (for the Bourne shell)
   or ``.cshrc`` (for the C shell).

2. The source code directory of PHCpack contains the directory ``Objects``,
   a directory with file ``makefile`` in it.
   Depending on whether you are on a Unix-like system or on a mac,
   you will have to edit the ``makefile`` so the ``MAKEFILE`` variable
   either refers to ``makefile_unix`` or to ``makefile_mac``.
   Once the makefile is adjusted you could type, just as a test,
   ``make phc`` to compile the main executable program.
   Note that for phcpy, you will not need this executable.

3. To make the shared object file, your python system needs to have been
   installed with the development version, that is: the file ``Python.h``
   must be available on your disk.  Often, following binary installations
   of the Python interpreter, this ``Python.h`` might be absent.

   If packaged distributions for the development version of Python fail,
   then you may have to download the source code 
   from <http://www.python.org>,
   configure, compile, and install the Python system.
   An additional benefit of such a Python installation is that then the
   Python interpreter could be built with the gnu-ada compiler,
   so both the scripting environment as the compiled code are using
   the same version of gcc.

4. Once you have located the file ``Python.h`` on your system,
   you most likely will have to adjust the definitions in the files
   ``makefile_unix`` or ``makefile_mac``.  Assign to the variables
   ``PYTHON`` and ``PYTHON3`` the directories where ``Python.h`` is.

5. In the directory ``Objects`` of the PHCpack source distribution,
   type ``make phcpy2c2.so`` to make the shared object file for python2,
   or type ``make phcpy2c3.so`` for the python3 version of phcpy.
   If all goes well, the shared object ``phcpy2c2.so`` can then be
   found in ``Python/PHCpy2/phcpy`` or ``phcpy2c3.so`` is then in
   ``Python/PHCpy3/phcpy.``

   If you run the Python interpreter in a terminal window in the
   directory ``phcpy`` of the source code distribution, then the
   ``import phcpy`` should already work.

6. To extend your python2 installation with phcpy, 
   go to the ``Python/PHCpy2`` directory and run
   ``python2 setup.py install`` as superuser or as sudoer.
   For python3, go to ``PHCpy/phcpy3`` and run
   ``python3 setup.py install`` as superuser or as sudoer.

about this document
-------------------

This document arose as an exercise in exploring restructured text
and Sphinx.  Even with a wonderful tool like Sphinx, this documentation
is (just as phcpy) very much a work in progress...
