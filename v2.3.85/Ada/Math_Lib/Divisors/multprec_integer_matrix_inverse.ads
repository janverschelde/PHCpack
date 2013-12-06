with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Integer_Matrices;

package Multprec_Integer_Matrix_Inverse is

-- DESCRIPTION :
--   This package offers routines to compute the inverse of a unimodular
--   integer matrix, using multiprecision arithmetic.

  function Is_Identity 
             ( A : Multprec_Integer_Matrices.Matrix ) return boolean;

  -- DESCRIPTION :
  --   Returns true if A is the identity matrix, false if otherwise.

  function Is_Inverse_Pair
             ( A,B : Multprec_Integer_Matrices.Matrix ) return boolean;

  -- DESCRIPTION :
  --   Returns true if A is the inverse of B and B the inverse of A,
  --   returns false if otherwise.
 
  function Inverse_of_Unimodular_Upper
             ( A : in Multprec_Integer_Matrices.Matrix )
             return Multprec_Integer_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the inverse of the unimodular upper triangular matrix in A.

  -- REQUIRED : A is upper triangular and unimodular.

  procedure Inverse ( A : in Multprec_Integer_Matrices.Matrix;
                      determinant : out Integer_Number;
                      B : out Multprec_Integer_Matrices.Matrix );

  -- DESCRIPTION :
  --   Computes the inverse of the square matrix A.
  --   If A is unimodular (determinant plus or minus 1),
  --   then B*A and A*B should be identity matrices.

  -- REQUIRED : A'range(1) = A'range(2).

end Multprec_Integer_Matrix_Inverse;
