***********
User Manual
***********

This chapter starts with a description of the blackbox solver,
provided by the solver module.  
If solutions of a start system for a polynomial homotopy
are available, then we may better call directly the path tracking
routines to solve the target system in the homotopy.
The path tracking functions are described in the second section.

The third section deals with the computation of positive dimensional
solution sets.  In a numerical irreducible decomposition of the
solution set of a polynomial system, generic points are computed
on all irreducible factors of the solution sets of all dimensions.
In a witness set representation, the number of generic points in
the witness set equals the degree of the pure dimensional solution set.
Solution sets can be computed in a top down fashion with cascade
homotopies or in a bottom up manner via diagonal homotopies.

The definitions of the polynomial systems which make interesting
examples and families of problems are illustrated in section four.

Pieri homotopies and Littlewood-Richardson homotopies 
solve problems in enumerative geometry.  The fifth section
of this chapter is concerned with numerical Schubert calculus.

Every polynomial in several variables has a Newton polytope,
spanned by the exponents of the monomials which occur with
a nonzero coefficient.  The mixed volume of a tuple of Newton
polytopes gives an upper bound for the number of isolated solutions.
Systems with exactly two monomials in every equation can be
solved fast, via unimodular coordinate transformations.
Section six of this chapter ends with an illustration of
the computation of power series solutions for algebraic curves.

Section seven and eight describe prototype modules for a
graphical user interface and a computational server.
The nineth and last section of this chapter sketches
the design of the C interface, the Python interface module,
and the wrappers to the C interface to PHCpack.

The last section collects the code snippets,
defined for the notebook extension of Jupyter,
as they pop up in the menus of the notebook.

.. toctree::
   :maxdepth: 2

   blackbox
   pathtrack
   posdimsols
   examfams
   numschub
   newtopes
   dashboard
   phcbwulf
   modphcpy2c3
   snippets
