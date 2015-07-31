the module phcpy.phcpy2c
========================

The Python scripts in the package phcpy call the wrappers
for the C interface to PHCpack.  Below is the list of all
functions exported by the shared object file phcpy2c.so.
The source code provides more detailed documentation.

wrappers to the C interface to PHCpack
--------------------------------------

A basic application of the primitive operations in phcpy2c
is an interactive reading of a polynomial system.
Assume the file example at /tmp/ contains a polynomial system,
then we can do the following:

::

   >>> from phcpy.phcpy2c import py2c_syscon_read_standard_system as readsys
   >>> from phcpy.phcpy2c import py2c_syscon_write_standard_system as writesys
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
   >>> from phcpy.phcpy2c import py2c_solve_system as solve
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
   >>> from phcpy.phcpy2c import py2c_solcon_write_standard_solutions as writesols
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

With these primitive operations in phcpy2c we can bypass the writing
and the parsing to strings.

functions in the module
-----------------------

The module ``phcpy2c`` wraps the C functions in the C interface to PHCpack.
The C interface to PHCpack was developed in the application of message passing
(MPI) to run the path trackers on distributed memory multiprocessing computers.
All functions documented below have their counterpart in C 
that are therefore then also directly accessible from C programs.

.. automodule:: phcpy2c
   :members:
