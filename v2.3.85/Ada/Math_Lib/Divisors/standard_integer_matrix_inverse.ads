with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Matrices;
with Standard_Integer64_Matrices;

package Standard_Integer_Matrix_Inverse is

-- DESCRIPTION :
--   This package offers routines to compute the inverse of a unimodular
--   integer matrix, using standard 32-bit or 64-bit arithmetic.

  function Is_Identity 
              ( A : Standard_Integer_Matrices.Matrix ) return boolean;
  function Is_Identity 
              ( A : Standard_Integer64_Matrices.Matrix ) return boolean;

  -- DESCRIPTION :
  --   Returns true if A is the identity matrix, false if otherwise.

  function Is_Inverse_Pair
              ( A,B : Standard_Integer_Matrices.Matrix ) return boolean;
  function Is_Inverse_Pair
              ( A,B : Standard_Integer64_Matrices.Matrix ) return boolean;

  -- DESCRIPTION :
  --   Returns true if A is the inverse of B and B the inverse of A,
  --   returns false if otherwise.
 
  function Inverse_of_Unimodular_Upper
              ( A : in Standard_Integer_Matrices.Matrix )
              return Standard_Integer_Matrices.Matrix;
  function Inverse_of_Unimodular_Upper
              ( A : in Standard_Integer64_Matrices.Matrix )
              return Standard_Integer64_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the inverse of the unimodular upper triangular matrix in A.

  -- REQUIRED : A is upper triangular and unimodular.

  procedure Inverse ( A : in Standard_Integer_Matrices.Matrix;
                      determinant : out integer32;
                      B : out Standard_Integer_Matrices.Matrix );
  procedure Inverse ( A : in Standard_Integer64_Matrices.Matrix;
                      determinant : out integer64;
                      B : out Standard_Integer64_Matrices.Matrix );

  -- DESCRIPTION :
  --   Computes the inverse of the square matrix A.
  --   If A is unimodular (determinant plus or minus 1),
  --   then B*A and A*B should be identity matrices.

  -- REQUIRED : A'range(1) = A'range(2).

end Standard_Integer_Matrix_Inverse;
