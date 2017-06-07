=================================
Welcome to phcpy's documentation!
=================================

This documentation describes a collection of Python modules
to compute solutions of polynomial systems using PHCpack.

Both phcpy and PHCpack are free and open source software packages;
you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software Foundation;
either version 2 of the License, or (at your option) any later version.

The main executable phc (polynomial homotopy continuation)
defined by the source code in PHCpack is a menu driven
and file oriented program.
The Python interface defined by phcpy replaces the files
with persistent objects allowing the user to work with
scripts or in interactive sessions.
The computationally intensive tasks such as path tracking
and mixed volume computations are executed as compiled code
so there will not be a loss of efficiency.
One application of phcpy is to run regression tests.

The source for PHCpack can be downloaded from
<http://www.math.uic.edu/~jan/download.html>

The Python interface phcpy to PHCpack is a programmer's interface.
You will need to compile the source code of PHCpack with
the GNU Ada compiler, available at <http://libre.adacore.com>.
In addition to the Python interpreter, the file Python.h of
of the developer's libraries will need to exist on your system.
Otherwise, PHCpack and phcpy are self contained
and do not require the installation of other software.

For Red Hat Linux 64-bit and some Mac OS X versions,
the distribution provides the shared object file phcpy2c3.so.
For phcpy to work, the file phcpy2c3.so must be present at the
same location as the Python modules.

See the Sphinx documentation in the directory doc
for examples on how to use phcpy and its modules.

Project History
===============

The Python interface to PHCpack got to a first start when
Kathy Piret met William Stein at the software for algebraic geometry
workshop at the IMA in the Fall of 2006.  
A preliminary version of this interface is described
in the 2008 PhD Thesis of Kathy Piret.

Version 0.0.1 originated at lecture 40 of MCS 507 in the Fall of 2012,
as an illustration of Sphinx.  In Spring of 2013, version 0.5.0 was
presented at a graduate computational algebraic geometry seminar.
Version 0.1.0 was prepared for presentation at EuroSciPy 2013 (August 2013).
Improvements using pylint led to version 0.1.1
and the module maps was added in version 0.1.2.
Version 0.1.4 added path trackers with a generator
so all solutions along a path are returned to the user.
Version 0.2.9 coincides with version 2.4 of PHCpack and gives access
to the first version of the GPU accelerated path trackers.

References
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

Acknowledgments
===============

This material is based upon work supported by the 
National Science Foundation under Grants 1115777 and 1440534.
Any opinions, findings, and conclusions or recommendations expressed 
in this material are those of the author(s) and do not necessarily 
reflect the views of the National Science Foundation. 
