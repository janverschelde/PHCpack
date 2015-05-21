Welcome to phcpy's documentation!
=================================

This documentation describes a collection of Python modules
to compute solutions of polynomial systems using PHCpack.
While PHCpack has been under development for over twenty years,
phcpy is still working through its proof-of-concept stage.
In its present state, working with phcpy will require persistence
and plenty of patience.

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

The computation of the mixed volume in phcpy calls MixedVol
(ACM TOMS Algorithm 846 of T. Gao, T.Y. Li, M. Wu) 
as it is integrated in PHCpack.
For double double and quad double arithmetic, PHCpack incorporates
the QD library of Y. Hida, X.S. Li, and D.H. Bailey.
See the References section for pointers to the literature.

For Red Hat Linux 64-bit and some Mac OS X versions,
the distribution provides the shared object file phcpy2c.so.
For phcpy to work, the file phcpy2c.so must be present at the
same location as the Python modules.

This document arose as an exercise in exploring restructured text
and Sphinx.  Even with a wonderful tool like Sphinx, this documentation
is (just as phcpy) very much a work in progress...
