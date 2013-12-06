with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;

package Standard_Numerical_Derivatives is

-- DESCRIPTION :
--   This package provides numerical differentiation routines for
--   functions defined over standard complex floating-point vectors.

-- FOR FUNCTIONS :

  type Complex_Multivariate_Function is 
    access function ( x : Vector ) return Complex_Number;

  function Diff1 ( f : Complex_Multivariate_Function; x : Vector;
                   i : integer32; h : double_float ) return Complex_Number;
  function Diff2 ( f : Complex_Multivariate_Function; x : Vector;
                   i : integer32; h : double_float ) return Complex_Number;
  function Diff3 ( f : Complex_Multivariate_Function; x : Vector;
                   i : integer32; h : double_float ) return Complex_Number;

  -- DESCRIPTION :
  --   Diff1 returns the central divided difference approximation to the
  --   partial derivative of f at the i-th component of x, with step h.
  --   Diff2 does one extrapolation and Diff3 extrapolates twice.

  -- REQUIRED : h > 0.

  -- ON RETURN :
  --   The error on return should be of order h^2, h^4, and h^6
  --   on the respective results of Diff1, Diff2, and Diff3.

-- FOR SYSTEMS :

  type Complex_Multivariate_System is
    access function ( x : Vector ) return Vector;

  function Diff1 ( f : Complex_Multivariate_System; n : integer32;
                   x : Vector; h : double_float ) return Matrix;
  function Diff2 ( f : Complex_Multivariate_System; n : integer32;
                   x : Vector; h : double_float ) return Matrix;
  function Diff3 ( f : Complex_Multivariate_System; n : integer32;
                   x : Vector; h : double_float ) return Matrix;

  -- DESCRIPTION :
  --   Returns an approximation for the Jacobian matrix of f at x,
  --   using the corresponding Diff's from above.

  -- REQUIRED : h > 0 and the vectors returned by f are of range 1..n.

end Standard_Numerical_Derivatives;
