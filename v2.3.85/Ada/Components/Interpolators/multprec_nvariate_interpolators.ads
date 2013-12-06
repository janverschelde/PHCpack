with text_io;                          use text_io;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Multprec_Floating_Numbers;        use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;         use Multprec_Complex_Numbers;
with Multprec_Complex_Vectors;         use Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;         use Multprec_Complex_VecVecs;
with Multprec_Complex_NesVecs;         use Multprec_Complex_NesVecs;
with Multprec_Complex_Polynomials;     use Multprec_Complex_Polynomials;

package Multprec_Nvariate_Interpolators is

-- DESCRIPTION :
--   This package provides an implementation of Newton interpolation with
--   divided differences in several variables over the complex numbers,
--   in multi-precision floating-point arithmetic.

-- CONSTRUCTORS :

  function Create ( n,d : natural32; x : VecVec; y : NesVec ) return NesVec;
  function Create_on_Square
                  ( n,d : natural32; x : VecVec; y : NesVec ) return NesVec;

  -- DESCRIPTION :
  --   Returns the n-dimensional matrix of all divided differences
  --   of order less than or equal than the degree d, over the grid
  --   defined by the vectors in x and the function values in y.
  --   The Create_on_Square assumes a square grid and computes also the
  --   higher order divided differences, which is good for testing.

  -- ON ENTRY :
  --   n          number of variables in the interpolating polynomial;
  --   d          degree of the interpolating polynomial;
  --   x          the grid consists of n vectors of range 0..d;
  --   y          function values at products of the grid vectors,
  --              as matrix of dimension (1+d)^n, n ranges of 0..d.

  -- ON RETURN :
  --   n-dimensional matrix f[x(1)(0)..x(1)(k_1);..;x(n)(0)..x(n)(k_n]
  --   of divided differences, where k_1 + .. + k_n <= d, representing
  --   the coefficients of (x(1) - x(1)(0))*..*(x(1) - x(1)(k_1))* ..
  --                   .. *(x(n) - x(n)(0))*..*(x(n) - x(n)(k_n)).

  function Expand ( f : NesVec; x : VecVec ) return Poly;
 
  -- DESCRIPTION :
  --   Writes the Newton form as expanded polynomial in n variables.

-- EVALUATORS :

  function Eval0 ( f : NesVec; x : VecVec; a : Vector ) return Complex_Number;
  function Eval0 ( file : file_type;
                   f : NesVec; x : VecVec; a : Vector ) return Complex_Number;

  function Eval ( f : NesVec; x : VecVec; a : Vector ) return Complex_Number;
  function Eval ( file : file_type;
                  f : NesVec; x : VecVec; a : Vector ) return Complex_Number;

  -- DESCRIPTION :
  --   Given the divided differences in f and the abscisses of the grid in x,
  --   Eval returns the function value of the interpolating polynomial at a.
  --   With the file as argument, extra information is printed on file.
  --   The f is created by the routines above using the same argument x.
  --   Eval0 is slower than Eval as it does not use Horner's scheme;
  --   the main use of Eval0 is to assert correctness, thus for testing.

-- DIAGNOSTICS :

  procedure Eval_Grid ( file : in file_type; x : in VecVec; y,f : in NesVec;
                        maxerr : out Floating_Number );

  -- DESCRIPTION :
  --   Writes the content of y and the evaluated grid points on file.

  -- ON ENTRY :
  --   file       must be opened for output;
  --   x          defines the grid of interpolation points;
  --   y          function values at products of the grid vectors;
  --   f          obtained as result of Create(n,d,x,y).

  -- ON RETURN :
  --   maxerr     maximal difference between sampled and evaluated point.

  function Maximal_Error ( x : VecVec; y,f : NesVec ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns the maximal difference between sampled and evaluated points.

end Multprec_Nvariate_Interpolators;
