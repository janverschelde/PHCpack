the module phcpy.phcpy2c3
=========================

The Python scripts in the package phcpy call the wrappers
for the C interface to the Ada code in PHCpack.  Below is the list of all
functions exported by the shared object file phcpy2c3.so.
The source code provides more detailed documentation.

design of the Python to C interface
-----------------------------------

The design of phcpy depends on PHClib, a library of various
collections of C functions, originally developed for message passing
with the MPI library.  This design is sketched in :numref:`figdesign`.

.. _figdesign:

.. figure:: ./figdesign.png
    :align: center

    The design of phcpy through the C interface to Ada code.

PHClib interfaces to the Ada routines through one single
Ada procedure ``use_c2phc.adb``.
The collection of parallel distibuted memory programs (MPI2phc)
using message passing (MPI) depends on PHClib.
All C functions that are exported to the Python interface have
their prototypes in the header file ``phcpy2c.h``
while the definitions in ``phcpy2c3.c`` call the proper routines
in PHClib.

the interface to PHCpack
------------------------

The module interface collects the functions that parse the string
representations for polynomials and solutions to pass their data 
through the C interface to the Ada code of PHCpack.  
The reverse operations return the string representations for polynomials 
and solutions as stored internally in PHCpack.

The functions exported by **phcpy.interface** concern the movement
of data between Python and PHCpack.  The `store_` methods parse strings
representing polynomials and solutions into the corresponding internal
data structures.  The corresponding `load_` methods take the internal
data structures for polynomials and solutions, stored in containers,
and show their corresponding representations as Python strings.
For example, consider the session

::

   >>> from phcpy.interface import store_standard_system, load_standard_system
   >>> store_standard_system(['x^2 - 1/3;'])
   >>> load_standard_system()
   ['x^2 - 3.33333333333333E-01;']

The session above illustrates the parsing of a system one could use
to approximate the square root of 1/3.  With standard double precision,
the 1/3 is approximated to about 15 decimal places.

wrappers to the C interface to PHCpack
--------------------------------------

A basic application of the primitive operations in phcpy2c3
is an interactive reading of a polynomial system.
Assume the file example at /tmp/ contains a polynomial system,
then we can do the following:

::

   >>> from phcpy.phcpy2c3 import py2c_syscon_read_standard_system as readsys
   >>> from phcpy.phcpy2c3 import py2c_syscon_write_standard_system as writesys
   >>> readsys()

   Reading a polynomial system...
   Is the system on a file ? (y/n/i=info) y

   Reading the name of the input file.
   Give a string of characters : /tmp/example
   0
   >>> writesys()
    2
   x^2+4*y^2-4;
   2*y^2-x;
   0
   >>> from phcpy.phcpy2c3 import py2c_solve_system as solve
   >>> solve(0)

   ROOT COUNTS :

   total degree : 4
   general linear-product Bezout number : 4
     based on the set structure :
        { x y }{ x y }
        { x y }{ y }
   mixed volume : 4
   stable mixed volume : 4
   4
   >>> from phcpy.phcpy2c3 import py2c_solcon_write_standard_solutions as writesols
   >>> writesols()
   4 2
   ===========================================================================
   solution 1 :
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x :  1.23606797749979E+00  -9.91383530201425E-119
    y :  7.86151377757423E-01   4.95691765100713E-119
   == err :  1.567E-16 = rco :  3.067E-01 = res :  3.331E-16 ==
   solution 2 :
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x :  1.23606797749979E+00  -9.91383530201425E-119
    y : -7.86151377757423E-01  -4.95691765100713E-119
   == err :  1.567E-16 = rco :  3.067E-01 = res :  3.331E-16 ==
   solution 3 :
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -3.23606797749979E+00   3.17242729664456E-117
    y : -1.58621364832228E-117   1.27201964951407E+00
   == err :  3.703E-16 = rco :  1.515E-01 = res :  4.441E-16 ==
   solution 4 :
   t :  1.00000000000000E+00   0.00000000000000E+00
   m : 1
   the solution for t :
    x : -3.23606797749979E+00   3.17242729664456E-117
    y :  1.58621364832228E-117  -1.27201964951407E+00
   == err :  3.703E-16 = rco :  1.515E-01 = res :  4.441E-16 ==
   0
   >>> 

With these primitive operations in phcpy2c3 we can bypass the writing
and the parsing to strings.
