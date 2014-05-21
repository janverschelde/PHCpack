PHClab is a collection of scripts to call phc from within a MATLAB 
or Octave session.  It provides an interface to the blackbox solver
for finding isolated solutions.  PHClab also interfaces to the numerical
irreducible decomposition, giving access to the tools to represent, 
factor, and intersect positive dimensional solution sets.

PHClab is distributed under the GNU GENERAL PUBLIC LICENSE,
see the file COPYING.txt for version 3 of this license.

Release history:
Version 1.0.1 April 24, 2008
Version 1.0.2 March 30, 2013, num2str bug fix in form_poly.m
Version 1.0.3 May 21, 2014, cmplx2str of Bor Plestenjak in form_poly.m

The main functions of PHClab are divided into two categories: 

1. computing approximations to all isolated solutions

  solve_system   calls the blackbox solver
  refine_sols    refines the solutions by application of Newton's method
  mixed_volume   computes the mixed volumes and solve a random system
  track          given target and start system + solutions, tracks paths 
  deflation      applies deflation for isolated singularities

2. tools for a numerical irreducible decomposition

  embed          adds extra hyperplanes and slack variables
  cascade        performs a cascade of homotopies
  phc_filter     filters junk on higher dimensional components
  eqnbyeqn       an equation-by-equation solver
  decompose      partitions a witness set along irreducible factors
  intersection   intersection of two witness sets

PHClab was tested on Matlab 6.5 and Octave 2.1.64 on computers running
Windows and Linux. On an Apple laptop running Mac OS X version 10.3.7, 
we executed PHClab in Octave 2.1.57.

To install and use PHClab, execute the following steps:
1. download the appropriate executable of phc;
2. download the PHClab distribution;
3. add PHClab to the path of MATLAB (or Octave).

The first command of PHClab which one must execute is set_phcpath.
This command takes one argument: the full path name of where the
program phc can be found.  Type help "set_phcpath" for more info.

Note: please use "/" in the path even it is for Windows machine. 

Recommended Reading:
"PHClab: A MATLAB/Octave Interface to PHCpack"
by Yun Guan and Jan Verschelde.
In M.E. Stillman, N.~Takayama, and J.~Verschelde, editors,
Software for Algebraic Geometry, volume 148 of The IMA Volumes in
Mathematics and its Applications, pages 15-32. Springer-Verlag, 2008.
