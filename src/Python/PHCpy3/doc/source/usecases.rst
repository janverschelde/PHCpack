********
Tutorial
********

The purpose of this chapter is to introduce phcpy via some use cases.

In all cases, there are three stages:

1. Given a problem, formulate the polynomial equations which
   describe the problem.  This formulation results in a list
   of string representations for polynomials in several variables.
   This list is the input for the second stage.
   
2. If the problem leads to a polynomial system for which only
   the isolated solution are of interest, then the blackbox solver
   as in the solve method of the solver method will do.
   Otherwise, for positive dimensional solution sets,
   cascades of homotopies and diagonal homotopies are needed.

3. The solvers return their results as string representations for
   solutions.  To process the solutions we convert the string 
   representations into Python dictionaries.

In all use cases, plots of the solutions are made with matplotlib.
To formulate the polynomial equations we may use sympy,
as illustrated in the design of 4-bar mechanisms.

In all use cases, we distinguish between general instances of a problem
(where the numbers for the parameters are chosen at random),
and specific instances of the problem.
For these specific instances, singular solutions are likely to occur.

.. toctree::
   :maxdepth: 2

   appolonius
   fourbar
   fourlines
   touchcircle
