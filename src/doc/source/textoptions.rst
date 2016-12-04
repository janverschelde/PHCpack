**********************************
Text Options of the Executable phc
**********************************

Text options of ``phc`` do not compute anything,
but provide text information about ``phc``.

While the regular options of ``phc`` start with a single dash ``-``,
the text options are double dashed, they start with ``--``.

phc --version : displays the current version string
===================================================

Typing at the command prompt ``phc --version`` displays
the :index:`version string` which includes the current :index:`version number`
and the :index:`release date`.

phc --license : writes the license to screen
============================================

The output of

::

   phc --license

is

::

   PHCpack is free and open source software.
   You can redistribute the code and/or modify it under
   the GNU General Pulic License as published by
   the Free Software Foundation.

phc --cite : how to cite PHCpack
================================

Typing ``phc --cite`` at the command prompt displays

::

   To cite PHCpack in publications use:

   Jan Verschelde:
   Algorithm 795: PHCpack: A general-purpose solver for polynomial
   systems by homotopy continuation.
   ACM Transactions on Mathematical Software, 25(2):251--276, 1999.

phc --help : writes helpful information to screen
=================================================

Typing at the command prompt ``phc --help`` provides some 
:index:`help` to get started with the quickest use, that is:
with the blackbox solver.

To obtain help about the blackbox solver, type ``phc -b --help``
or ``phc --help -b`` where the ``--help`` may be abbreviated by ``-h``.
