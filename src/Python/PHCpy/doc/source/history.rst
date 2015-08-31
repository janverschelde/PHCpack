Project History
===============

The Python interface to PHCpack got to a first start when
Kathy Piret met William Stein at the software for algebraic geometry
workshop at the IMA in the Fall of 2006.  
The first version of this interface is described
in the 2008 PhD Thesis of Kathy Piret.

The implementation of the Python bindings depend on the C interface
to PHCpack, developed for use with message passing on distributed
memory computers.

Version 0.0.1 originated at lecture 40 of MCS 507 in the Fall of 2012,
as an illustration of Sphinx.  In Spring of 2013, version 0.5.0 was
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
