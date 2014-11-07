with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Matrices;
with Double_Double_Matrices;
with Quad_Double_Matrices;
with Multprec_Floating_Matrices;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Multprec_Complex_Matrices;

package Random_Conditioned_Matrices is

-- DESCRIPTION :
--   A conditioned matrix is a matrix with a prescribed condition number.
--   This package offers functions to generated random matrices with
--   a given condition number.

  function Singular_Value_Matrix
             ( n : integer32; c : double_float )
             return Standard_Floating_Matrices.Matrix;
  function Singular_Value_Matrix
             ( n : integer32; c : double_float )
             return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns an n-dimensional diagonal matrix that has on its diagonal 
  --   the values starting at 1.0 and ending with c,
  --   with values in between on a logarithmic scale:
  --   for i in 2..n-1: 10**((i-1)*m/(n-1)), where m = log10(|c|).

  function Random_Conditioned_Matrix
             ( n : integer32; c : double_float )
             return Standard_Floating_Matrices.Matrix;
  function Random_Conditioned_Matrix
             ( n : integer32; c : double_float )
             return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a random n-dimensional matrix with condition number c,
  --   computed with standard floating-point double precision numbers.
  --   The matrix of return will have the condition number c,
  --   provided this condition number is less than 1.0E+16.

  function Random_Conditioned_Matrix
             ( n : integer32; c : double_float )
             return Double_Double_Matrices.Matrix;
  function Random_Conditioned_Matrix
             ( n : integer32; c : double_float )
             return DoblDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a random n-dimensional matrix with condition number c.
  --   computed with double double arithmetic.
  --   This will work for conditioned numbers no larger than 1.0E+32.
  --
  function Random_Conditioned_Matrix
             ( n : integer32; c : double_float )
             return Quad_Double_Matrices.Matrix;
  function Random_Conditioned_Matrix
             ( n : integer32; c : double_float )
             return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a random n-dimensional matrix with condition number c.
  --   computed with double double arithmetic.
  --   This will work for conditioned numbers no larger than 1.0E+64.
  --
  function Random_Conditioned_Matrix
             ( n : integer32; c : double_float )
             return Multprec_Floating_Matrices.Matrix;
  function Random_Conditioned_Matrix
             ( n : integer32; c : double_float )
             return Multprec_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a random n-dimensional matrix with condition number c.
  --   computed with arbitrary multiprecision arithmetic.
  --   The size of the multiprecision numbers is set as twice
  --   the magnitude of the condition number c.

end Random_Conditioned_Matrices;
