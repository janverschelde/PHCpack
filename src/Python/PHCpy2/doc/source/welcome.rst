***************
Getting Started
***************

This documentation describes a collection of Python modules
to compute solutions of polynomial systems using PHCpack.

This work is licensed under 
a Creative Commons Attribution-Share Alike 3.0 License.

The computation of the mixed volume in phcpy calls MixedVol
(ACM TOMS Algorithm 846 of T. Gao, T.Y. Li, M. Wu) 
as it is integrated in PHCpack.
DEMiCs (Dynamic Enumeration of all Mixed Cells,
by T. Mizutani, A. Takeda, and M. Kojima) is faster than MixedVol
for larger systems with many different supports.
A function to compute mixed volumes with DEMiCs is available in phcpy.

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
   Typing ``gcc -v`` should show ``for GNAT GPL`` in the reply,
   or most recently (in 2018): ``GNAT Community 2018``.

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

2. Since version 0.6.4, the code depends on the quad double library QDlib,
   available at <http://crd-legacy.lbl.gov/~dhbailey/mpdist>.
   Version 2.3.17 is available at <https://github.com/scibuilder/QD>.

   On Linux systems, make sure to compile and install the library
   with the option -fPIC.  When configuring, run configure as
   ``./configure CXX=/usr/bin/g++ CXXFLAGS='-fPIC -O3'`` to set the flags
   of the c++ compiler.

   If you rather would not (or cannot) install QDlib, then it is possible
   to compile the library in your home directory.  All you need is to be
   aware of the location of the header files for the include statement
   and you need the ``qdlib.a`` file for linking.
   For an example, consider the makefile for Windows computers.
   The ``makefile_windows`` builds ``phc.exe`` with the QD library
   compiled in a home directory, not installed on the system.

3. The source code directory of PHCpack contains the directory ``Objects``,
   a directory with file ``makefile`` in it.
   Depending on whether you are on a Unix-like system or on a mac,
   you will have to edit the ``makefile`` so the ``MAKEFILE`` variable
   either refers to ``makefile_unix`` or to ``makefile_mac``.
   Once the makefile is adjusted you could type, just as a test,
   ``make phc`` to compile the main executable program.
   Note that for phcpy, you will not need this executable.

4. To make the shared object file, your python system needs to have been
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

5. Once you have located the file ``Python.h`` on your system,
   you most likely will have to adjust the definitions in the files
   ``makefile_unix`` or ``makefile_mac``.  Assign to the variables
   ``PYTHON`` and ``PYTHON3`` the directories where ``Python.h`` is.

6. In the directory ``Objects`` of the PHCpack source distribution,
   type ``make phcpy2c2.so`` to make the shared object file for python2,
   or type ``make phcpy2c3.so`` for the python3 version of phcpy.
   If all goes well, the shared object ``phcpy2c2.so`` can then be
   found in ``Python/PHCpy2/phcpy`` or ``phcpy2c3.so`` is then in
   ``Python/PHCpy3/phcpy.``

   If you run the Python interpreter in a terminal window in the
   directory ``phcpy`` of the source code distribution, then the
   ``import phcpy`` should already work.

7. To extend your python2 installation with phcpy, 
   go to the ``Python/PHCpy2`` directory and run
   ``python2 setup.py install`` as superuser or as sudoer.
   For python3, go to ``PHCpy/phcpy3`` and run
   ``python3 setup.py install`` as superuser or as sudoer.

The documentation is typeset with Sphinx.  Sphinx depends on the
default version of Python on your system.  If phcpy is installed
with a different version of Python than the version Sphinx depends on,
then this may cause problems to typeset the documentation.

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

We check the installation at the command prompt,
as shown in :numref:`figphcpyinsage`.

.. _figphcpyinsage:

.. figure:: ./figphcpyinsage.png
    :align: center

    Importing phcpy in a Sage terminal session.

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

Importing phcpy apparently changes the configuration of the signal 
handlers which may lead Sage to crash when exceptions occur.
Thanks to Marc Culler for reporting this problem
and for suggesting a work around:

::

   sage: import phcpy
   sage: from cysignals import init_cysignals
   sage: init_cysignals()
   sage: pari(1)/pari(0)

Without the ``init_cysignals()``,
the statement ``pari(1)/pari(0)`` crashes Sage.
With the ``init_cysignals()``, the ``PariError`` exception is handled
and the user can continue the Sage session.

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
With static linking, the dependencies on the gnat runtime libraries are removed
and the Sage python interpreter could be extended with version 0.6.2.
Better support of Laurent polynomial systems was added in version 0.6.8.
In version 0.6.9, the large module sets.py was divided up, leading to
the new modules cascades.py, factor.py, and diagonal.py.
Code snippets for jupyter notebook menu extensions were defined
in version 0.7.4.

Version 0.9.2 was presented in a talk at the Python devroom at FOSDEM 2019.
Another important milestone was the poster presentation
at the 18th Python in Science Conference, SciPy 2019,
with a paper which appeared in the proceedings (see the references below).

As the original goal of phcpy was on exporting the functionality of
PHCpack, its design is functional and phcpy is a collection of modules.
Since version 1.0.0, two class definitions were added,
one to represent systems of polynomials and another class to
represent solutions of polynomial systems.
The 1.x.y releases series will continue to develop phcpy
with an object-oriented interface to the functionality of PHCpack.

references
==========

1. T. Gao, T. Y. Li, M. Wu:
   **Algorithm 846: MixedVol: a software package for mixed-volume 
   computation.**
   *ACM Transactions on Mathematical Software*, 31(4):555-560, 2005.

#. Y. Hida, X.S. Li, and D.H. Bailey:
   **Algorithms for quad-double precision floating point arithmetic.**
   In *15th IEEE Symposium on Computer Arithmetic (Arith-15 2001)*,
   11-17 June 2001, Vail, CO, USA, pages 155-162.
   IEEE Computer Society, 2001.
   Shortened version of Technical Report LBNL-46996.

#. A. Leykin and J. Verschelde.
   **Interfacing with the numerical homotopy algorithms in PHCpack.**
   In N. Takayama and A. Iglesias, editors, *Proceedings of ICMS 2006*,
   volume 4151 of *Lecture Notes in Computer Science*,
   pages 354--360. Springer-Verlag, 2006.

#. T. Mizutani and A. Takeda.
   **DEMiCs: A software package for computing the mixed volume via
   dynamic enumeration of all mixed cells.**
   In M. E. Stillman, N. Takayama, and J. Verschelde, editors,
   *Software for Algebraic Geometry*, volume 148 of The IMA Volumes in
   Mathematics and its Applications, pages 59-79. Springer-Verlag, 2008.

#. T. Mizutani, A. Takeda, and M. Kojima.
   **Dynamic enumeration of all mixed cells.**
   *Discrete Comput. Geom.* 37(3):351-367, 2007.

#. J. Otto, A. Forbes, and J. Verschelde.
   **Solving Polynomial Systems with phcpy.**
   In the *Proceedings of the 18th Python in Science Conference (SciPy 2019)*,
   edited by Chris Calloway, David Lippa, Dillon Niederhut and David Shupe,
   pages 58-64, 2019. 

#. K. Piret:
   **Computing Critical Points of Polynomial Systems 
   using PHCpack and Python.**
   PhD Thesis, University of Illinois at Chicago, 2008.

#. A. J. Sommese, J. Verschelde, and C. W. Wampler.
   **Numerical irreducible decomposition using PHCpack.**
   In *Algebra, Geometry, and Software Systems*, 
   edited by M. Joswig and N. Takayama,
   pages 109-130. Springer-Verlag, 2003.

#. J. Verschelde:
   **Algorithm 795: PHCpack: A general-purpose solver for polynomial
   systems by homotopy continuation.**
   *ACM Transactions on Mathematical Software*, 25(2):251--276, 1999.

#. J. Verschelde:
   **Modernizing PHCpack through phcpy.**
   In Proceedings of the 6th European Conference on Python in Science
   (EuroSciPy 2013), edited by Pierre de Buyl and Nelle Varoquaux,
   pages 71-76, 2014, available at
   <http://arxiv.org/abs/1310.0056>.

#. J. Verschelde and X. Yu:
   **Polynomial Homotopy Continuation on GPUs.**
   *ACM Communications in Computer Algebra*, 49(4):130-133, 2015.

acknowledgments
===============

The PhD thesis of Kathy Piret (cited above) described the
development of a first Python interface to PHCpack.
The 2008 ``phcpy.py`` provided access to the blackbox solver,
the path trackers, and the mixed volume computation.

In the summer of 2017, Jasmine Otto helped with the setup of
jupyterhub and the definition of a SageMath kernel.
Code snippets with example uses of ``phcpy`` in a Jupyter notebook
were introduced during that summer.  The code snippets,
listed in a chapter of this document, provide another good way
to explore the capabilities of the software.

This material is based upon work supported by the 
National Science Foundation under Grants 1115777, 1440534, and 1854513.
Any opinions, findings, and conclusions or recommendations expressed 
in this material are those of the author(s) and do not necessarily 
reflect the views of the National Science Foundation. 

about this document
===================

This document arose as an exercise in exploring restructured text and Sphinx.
All good software documents contain the following four items: 
an installation guide, a getting started, a tutorial, and a reference manual.
This document combines all four.
