PHCpack
=======

PHCpack is a software package to solve polynomial systems by homotopy continuation methods.

A polynomial system is given as a sequence of polynomials in several variables.
Homotopy continuation methods operate in two stages.  In the first stage, a family of polynomial systems
(the so-called homotopy) is constructed.  This homotopy contains a polynomial system with known solutions.
In the second stage, numerical continuation methods are applied to track the solution paths defined by
the homotopy, starting at the known solutions and leading to the solutions of the given polynomial system.

Version 1.0 of PHCpack has been archived by ACM Transactions of Mathematical Software (ACM TOMS) as Algorithm 795.
PHCpack incorporates MixedVol (Algorithm 846 of ACM TOMS by T. Gao, T.Y. Li, and M. Wu) to compute mixed volumes fast.
DEMiCs (Dynamic Enumeration of all Mixed Cells, by T. Mizutani, A. Takeda, and M. Kojima), computes mixed volumes 
at a faster pace than MixedVol for larger systems with many different supports.  DEMiCs is also integrated into PHCpack.
For its double double and quad double arithmetic, PHCpack contains QDlib
of Y. Hida, X.S. Li, and D.H. Bailey.
For triple double, and other multiple doubles (penta, octo, deca),
code generated from the CAMPARY software is used.
CAMPARY is the CudA Multiple Precision ARithmetic librarY,
by Mioara Joldes, Olivier Marty, Jean-Michel Muller,
Valentina Popescu and Warwick Tucker.

This material is based upon work supported by the National Science Foundation
under Grants No. 9804846, 0105739, 0134611, 0410036, 0713018, 1115777,
1440534, and 1854513.
Any opinions, findings, and conclusions or recommendations expressed in this material
are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.

Executable versions of the code for Linux, MacOS X, and Windows are
available via the assets of the releases (tags).

Links to the documentation:

* The documentation of PHCpack in html format is
  at <http://www.math.uic.edu/~jan/phcpack_doc_html>
  and its pdf version is
  at <http://www.math.uic.edu/~jan/PHCpack.pdf>.

* The documentation of phcpy in html format is
  at <http://www.math.uic.edu/~jan/phcpy_doc_html>
  and its pdf version is
  at <http://www.math.uic.edu/~jan/phcpydoc.pdf>.

The restructured text source for the documentation for PHCpack starts at
<https://github.com/janverschelde/PHCpack/tree/master/src/doc/source>
and for phcpy at
<https://github.com/janverschelde/PHCpack/tree/master/src/Python/PHCpy/doc/source>.

An alire (<https://alire.ada.dev>) workspace is provided.

The default configuration is for windows.
To change manually, replace the occurrences of the string "windows" 
by "linux" in the files in the config folder.

Type 'alr build' to build and 'alr clean' to clean.

To test the libPHCpack, first do

cp lib/libPHCpack.* src/lib

and then navigate to the src/Julia folder and do

julia showVersion.jl

or navigate to the src/Python/PHCpy/phcpy folder and do

python version.py
