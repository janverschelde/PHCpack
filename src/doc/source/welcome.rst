.. PHCpack documentation master file, created by
   sphinx-quickstart on Sun Jan 27 13:05:16 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PHCpack's documentation!
===================================

PHCpack implements a collection of algorithms
to solve polynomial systems by homotopy continuation methods.

On input is a sequence of polynomials in several variables,
on output are the solutions to the polynomial system given on input.
The computational :index:`complexity`
of this problem is #P-hard because of
the exponential growth of the number of solutions as the number of
input polynomials increases.  For example, ten polynomials of degree two
may intersect in 1,024 isolated points (that is two to the power ten).
Twenty quadratic polynomials may lead to 1,048,576 solutions
(that is 1,024 times 1,024).  So it is not too difficult to write
down small input sequences that lead to a huge output.

Even as the computation of the total number of solutions may take
a long time, numerical homotopy continuation methods have the advantage
that they compute one solution after the other.  A homotopy is a family
of polynomial systems, connecting the system we want to solve with an
easier to solve system, which is called the start system.
Numerical continuation methods track the solution paths starting at
known solutions of an easier system to the system we want to solve.
We have an :index:`optimal homotopy`
if every path leads to a solution, that is: there are no divergent paths.

PHCpack offers optimal homotopies for systems that resemble 
linear-product structures, for geometric problems in enumerative geometry,
and for sparse polynomial systems with sufficiently generic choices 
of the coefficients.
While mathematically this sounds all good, most systems arising in
practical applications have their own peculiar structure
and so most homotopies will lead to diverging solution paths.
In general, a polynomial system may have solution sets of many
different dimensions, which renders the solving process challenging
but a the same time still very interesting.

Version 1.0 of PHCpack was archived as Algorithm 795
by ACM Transactions on Mathematical Software.  
PHCpack is open source and :index:`free software`
which gives any user the same rights as in free speech.
You can redistribute PHCpack and/or modify it under the terms of 
the GNU General Public :index:`License`
as published by the Free Software Foundation; 
either version 2 of the License, or (at your option) any later version.
