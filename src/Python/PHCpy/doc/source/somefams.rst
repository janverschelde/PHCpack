some families
=============

Polynomial systems often occur in families and are defined
for any number of equations and variables.

the cyclic n-roots problem
--------------------------

One such noteworthy family is the cyclic n-roots problem:

::

   >>> from phcpy.families import cyclic
   >>> c4 = cyclic(4)
   >>> for p in c4: print p
   ... 
   x0 + x1 + x2 + x3;
   x0*x1 + x1*x2 + x2*x3 + x3*x0;
   x0*x1*x2 + x1*x2*x3 + x2*x3*x0 + x3*x0*x1;
   x0*x1*x2*x3 - 1;
   >>> 

functions in the module
-----------------------

.. automodule:: families
   :members:
