about phcpy
===========

This section describes some milestones in the development history.

project history
---------------

The Python interface to PHCpack got to a first start when
Kathy Piret met William Stein at the software for algebraic geometry
workshop at the IMA in the Fall of 2006.  
The first version of this interface is described
in the 2008 PhD Thesis of Kathy Piret.

The implementation of the Python bindings depend on the C interface
to PHCpack, developed for use with message passing on distributed
memory computers.

Version 0.0.1 originated at lecture 40 of MCS 507 in the Fall of 2012,
as an illustration of Sphinx.  In Spring of 2013, version 0.0.5 was
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
Sweep homotopies to explore the parameter space with detection and
location of singularities along the solution paths were exported
in the module sweepers.py in version 0.3.3 of phcpy.
With the addition of a homotopy membership test in verion 0.3.7,
the sets.py module provides the key ingredients for a numerical
irreducible decomposition.

references
----------

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
   *ACM Communications in Computer Algebra*, to appear.

acknowledgments
---------------

This material is based upon work supported by the 
National Science Foundation under Grants 1115777 and 1440534.
Any opinions, findings, and conclusions or recommendations expressed 
in this material are those of the author(s) and do not necessarily 
reflect the views of the National Science Foundation. 
