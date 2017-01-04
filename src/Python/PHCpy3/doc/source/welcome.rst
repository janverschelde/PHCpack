***************
Getting Started
***************

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
==============

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
================

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

extending Sage with phcpy
=========================

The current version of Sage uses python2.
So the instructions to extend Sage with phcpy work only with
the Python2 version of phcpy.

If you have installed Sage from source on your computer,
then this installation comes with its own python libraries
and interpreter.  Then it is not too much work any more
(in comparison to the steps in last section) to extend
the python interpreter of Sage with phcpy.

On Linux systems, locate the python interpreter of Sage.
Most likely this python is ``/local/bin`` of in the downloaded directory.
Use the absolute path name for the location of the Sage python
interpreter and navigate to the ``Python/PHCpy2`` directory which
contains the ``setup.py`` for phcpy.
Once in ``Python/PHCpy2``, type ``python setup.py install`` at
the command prompt.  This does not require superuser access,
but you must execute this setup with the same account you used
to install Sage with.

On Mac OS X, extending Sage with phcpy requires a bit more work as
the ``phcpy2c2.so`` must be compiled with the Python library that
comes with the Sage installation.  For this, the ``makefile_mac``
must be modified with the correct definition for the location of
the Python library of Sage, as defined by ``SAGEPYTHONLIB``.
With this definition, do ``make sage_phcpy2c2.so`` and then move
this file under the name ``phcpy2c2.so`` to the directory
``/Python/PHCpy2/phcpy``.  The installation is then similar as
for Linux, type ``python setup.py install`` at the command prompt
in the directory where ``setup.py`` exists and for ``python``
using the absolute file name of the executable, e.g., type
``/Users/jan/Downloads/sage-7.2/local/bin/python setup.py install``.

project history
===============

This section describes some milestones in the development history.

The Python interface to PHCpack got to a first start when
Kathy Piret met William Stein at the software for algebraic geometry
workshop at the IMA in the Fall of 2006.  
The first version of this interface is described
in the 2008 PhD Thesis of Kathy Piret.

The implementation of the Python bindings depend on the C interface
to PHCpack, developed for use with message passing on distributed
memory computers.

Version 0.0.1 originated at lecture 40 of MCS 507 in the Fall of 2012,
as an illustration of Sphinx.  In Spring of 2013, version 0.0.5 was
presented at a graduate computational algebraic geometry seminar.
Version 0.1.0 was prepared for presentation at EuroSciPy 2013 (August 2013).
Improvements using pylint led to version 0.1.1
and the module maps was added in version 0.1.2.
Version 0.1.4 added path trackers with a generator
so all solutions along a path are returned to the user.
Multicore path tracking was added in version 0.1.7.

The paper **Modernizing PHCpack through phcpy**
written for the EuroSciPy 2013 proceedings 
and available at <http://arxiv.org/abs/1310.0056>
describes the design of phcpy.

Version 0.2.9 coincides with version 2.4 of PHCpack and gives access
to the first version of the GPU accelerated path trackers.
Sweep homotopies to explore the parameter space with detection and
location of singularities along the solution paths were exported
in the module sweepers.py in version 0.3.3 of phcpy.
With the addition of a homotopy membership test in verion 0.3.7,
the sets.py module provides the key ingredients for a numerical
irreducible decomposition.  
Version 0.5.0 introduced Newton's method on power series.
Use cases were added to the documentation in versions 0.5.2, 0.5.3, and 0.5.4.

references
==========

1. T. Gao, T. Y. Li, M. Wu:
   **Algorithm 846: MixedVol: a software package for mixed-volume 
   computation.**
   *ACM Transactions on Mathematical Software*, 31(4):555-560, 2005.

2. Y. Hida, X.S. Li, and D.H. Bailey:
   **Algorithms for quad-double precision floating point arithmetic.**
   In *15th IEEE Symposium on Computer Arithmetic (Arith-15 2001)*,
   11-17 June 2001, Vail, CO, USA, pages 155-162.
   IEEE Computer Society, 2001.
   Shortened version of Technical Report LBNL-46996.

3. A. Leykin and J. Verschelde.
   **Interfacing with the numerical homotopy algorithms in PHCpack.**
   In N. Takayama and A. Iglesias, editors, *Proceedings of ICMS 2006*,
   volume 4151 of *Lecture Notes in Computer Science*,
   pages 354--360. Springer-Verlag, 2006.

4. K. Piret:
   **Computing Critical Points of Polynomial Systems 
   using PHCpack and Python.**
   PhD Thesis, University of Illinois at Chicago, 2008.

5. A. J. Sommese, J. Verschelde, and C. W. Wampler.
   **Numerical irreducible decomposition using PHCpack.**
   In *Algebra, Geometry, and Software Systems*, 
   edited by M. Joswig and N. Takayama,
   pages 109-130. Springer-Verlag, 2003.

6. J. Verschelde:
   **Algorithm 795: PHCpack: A general-purpose solver for polynomial
   systems by homotopy continuation.**
   *ACM Transactions on Mathematical Software*, 25(2):251--276, 1999.

7. J. Verschelde:
   **Modernizing PHCpack through phcpy.**
   In Proceedings of the 6th European Conference on Python in Science
   (EuroSciPy 2013), edited by Pierre de Buyl and Nelle Varoquaux,
   pages 71-76, 2014, available at
   <http://arxiv.org/abs/1310.0056>.

8. J. Verschelde and X. Yu:
   **Polynomial Homotopy Continuation on GPUs.**
   *ACM Communications in Computer Algebra*, 49(4):130-133, 2015.

acknowledgments
===============

This material is based upon work supported by the 
National Science Foundation under Grants 1115777 and 1440534.
Any opinions, findings, and conclusions or recommendations expressed 
in this material are those of the author(s) and do not necessarily 
reflect the views of the National Science Foundation. 

about this document
===================

This document arose as an exercise in exploring restructured text and Sphinx.
All good software document consists of four items: an installation guide,
a getting started, a tutorial, and a reference manual.
This document combines all four.
In its current state, phcpy is a collection of modules, with a focus
on exporting the functionality of PHCpack.  The design is functional.
The package does not define nor export an object oriented interface.
