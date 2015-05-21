monomial maps
=============

Systems that have exactly two monomials with nonzero coefficient
in every equation are called binomial systems.
Although such binomial systems are very particular,
because of their sparse structure, they can be solved much faster.

solving binomial systems
------------------------

The irreducible components of
positive dimensional solution sets of binomial systems
have coordinates that can be represented by maps of monomials 
in free independent variables.  In this representation, there
are as many free variables as the dimension of the solution set.
The module ``maps`` exports a solver for binomial systems.

In the example below, we consider a simple system
of two binomials in three variables:

::

   >>> f = [ 'x**2*y - z*x;', 'x**2*z - y**2*x;' ]
   >>> from phcpy.maps import binomial_solver
   >>> from phcpy.maps import solve_binomials
   >>> maps = solve_binomials(3,f)
   >>> for map in maps: print map
   ... 
   ['x - 0', 'y - (1+0j)*t1**1', 'z - (1+0j)*t2**1', 'dimension = 2', 'degree = 1']
   ['x - (1+0j)*t1**1', 'y - (1+0j)*t1**2', 'z - (1+0j)*t1**3', 'dimension = 1', 'degree = 3']
   ['x - (1+0j)*t1**1', 'y - 0', 'z - 0', 'dimension = 1', 'degree = 1']
   >>> 

In the output above we recognize the twisted cubic,
the x-axis, and the yz-plane as the three solution sets.

functions in the module
-----------------------

.. automodule:: maps
   :members:
