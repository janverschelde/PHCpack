with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Integer_Matrices;          use Multprec_Integer_Matrices;

package Multprec_Smith_Normal_Form is

-- DESCRIPTION :
--   This package provides a diagonalization routine for a matrix of
--   multiprecision integer numbers, along with some utilities.

  function Identity ( n : natural32 ) return Matrix;

  -- DESCRIPTION :
  --   Returns the n-by-n identity matrix.

  function Diagonal ( mat : Matrix ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the matrix is diagonal, returns false otherwise.

  function Rank_of_Diagonal_Matrix ( d : matrix ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of nonzero consecutive elements on
  --   the diagonal matrix d.

  -- REQUIRED : d is a diagonal matrix and all its nonzero elements
  --   on the diagonal are consecutive.

  procedure Diagonalize ( u,a,v : in out Matrix );

  -- DESCRIPTION :
  --   Returns the Smith normal form of the matrix, a := u*a*v.

  -- ON ENTRY :
  --   u        identity matrix of range a'range(1);
  --   a        matrix to diagonalize;
  --   v        identity matrix of range a'range(2).

  -- ON RETURN :
  --   u        unimodular matrix of left multiplications;
  --   a        diagonal matrix, equals u*a*v;
  --   v        unimodular matrix of right multiplications.

end Multprec_Smith_Normal_Form;
