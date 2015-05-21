interface to PHCpack
====================

The module interface collects the functions that parse the string
representations for polynomials and solutions to pass their data 
through the C interface of PHCpack.  The reverse operations return
the string representations for polynomials and solutions as stored
internally in PHCpack.

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

functions in the module
-----------------------

.. automodule:: interface
   :members:
