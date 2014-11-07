with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Integer_Matrices;
with Multprec_Integer64_Matrices;
with Multprec_Floating_Matrices;
with Multprec_Complex_Matrices;

package Multprec_Random_Matrices is

-- DESCRIPTION :
--   Offers routines to generate random matrices of multi-precision numbers.

  function Random_Matrix ( n,m : natural32; low,upp : integer32 )
                         return Multprec_Integer_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range 1..n,1..m with entries in [low,upp].

  function Random_Matrix ( n,m,sz : natural32 )
                         return Multprec_Integer_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range 1..n,1..m with entries of size sz.

  function Random_Matrix ( n,m,sz : natural32 )
                         return Multprec_Integer64_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range 1..n,1..m with entries of size sz.

  function Random_Matrix ( n,m,sz : natural32 )
                         return Multprec_Floating_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range 1..n,1..m with entries of size sz.

  function Orthogonalize ( mat : Multprec_Floating_Matrices.Matrix )
                         return Multprec_Floating_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the orthogonal matrix with the same column span as mat.

  function Random_Orthogonal_Matrix
             ( n,m,sz : natural32 ) return Multprec_Floating_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a random matrix where the columns form an orthonormal basis.

  function Random_Matrix ( n,m,sz : natural32 )
                         return Multprec_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range 1..n,1..m with entries of size sz.

  function Orthogonalize ( mat : Multprec_Complex_Matrices.Matrix )
                         return Multprec_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the orthogonal matrix with the same column span as mat.

  function Random_Orthogonal_Matrix
             ( n,m,sz : natural32 ) return Multprec_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a random matrix where the columns form an orthonormal basis.

end Multprec_Random_Matrices;
