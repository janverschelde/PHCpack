with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Complex_Vectors;          use Multprec_Complex_Vectors;
with Multprec_Complex_Matrices;         use Multprec_Complex_Matrices;

package Multprec_Numerical_Rank is

-- DESCRIPTION :
--   The singular value decomposition gives the numerical rank of a matrix,
--   using multiprecision complex floating-point arithmetic.

  function Numerical_Rank ( S : Vector; tol : double_float ) return natural32;

  -- DESCRIPTION :
  --   Returns the numerical rank of a matrix, based on its
  --   singular values in S, with respect to the tolerance tol.

  procedure Numerical_Rank
              ( A : in out Matrix; tol : in double_float;
                S : out Vector; U,V : out Matrix;
                rco : out Floating_Number; rank : out natural32 );

  -- DESCRIPTION :
  --   Returns numerical rank, estimate for the inverse condition number,
  --   and the singular value decomposition of A.

  -- ON ENTRY:
  --   A        matrix for which numerical rank is desired;
  --   tol      tolerance to decide the numerical rank.
  
  -- ON RETURN :
  --   A        modified version of the matrix;
  --   S        vector of range 1,..,min(A'last(1),A'last(2)+1)
  --            with the singular values of the matrix A.
  --   U        orthogonal square matrix of dimension A'range(1);
  --   V        orthogonal square matrix of dimension A'range(2);
  --   rco      estimate for inverse condition number;
  --   rank     numerical rank of the matrix.

  function Numerical_Rank ( A : Matrix; tol : double_float ) return natural32;

  -- DESCRIPTION :
  --   Computes the SVD of A and returns the numerical rank.

end Multprec_Numerical_Rank;
